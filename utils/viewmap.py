"""
umap_annotator.py — interactive scatter annotation for Jupyter notebooks.

    annotateScatter(df, df_meta=None, image_paths={"sample1": "/path/a.tiff", ...})

  df        : DataFrame indexed by sample name, with columns "x","y" (any 2D
              embedding, e.g. UMAP). One point per row.
  df_meta   : optional DataFrame, same index, arbitrary columns -> shown in
              the hover tooltip.
  df_features : optional DataFrame, same index as df/df_meta, numeric
              feature columns used to train a one-vs-rest logistic
              regression classifier per label (only for labels with at
              least `min_label_count` positive AND negative examples).
              When provided, clicking a point fetches that sample's
              predicted probability for every label it doesn't have yet,
              and shows the ones that look likely as one-click "suggested"
              chips in the popup; confirming one adds it as a real label
              (same as typing it in) and retrains that label's classifier.
              Requires scikit-learn; the suggestion feature is silently
              disabled (with a printed note) if either is missing.
  augFunc   : optional (X, y) -> (X_aug, y_aug) hook, currently a passthrough
              stub -- reserved for future data augmentation before each
              label's classifier is (re)trained.
  clf_params : optional dict of LogisticRegression kwargs; defaults to
              {"penalty": "l2", "C": 10, "class_weight": "balanced",
               "solver": "liblinear", "max_iter": 1000}.
  min_label_count : minimum number of positive AND negative examples a
              label needs before a classifier is trained for it.
  image_paths : dict mapping sample name -> full path of that sample's
              thumbnail image (any format tifffile/PIL can read). Samples
              missing from this dict simply won't have an image in the
              popup. Each thumbnail's sibling "info.json" (same directory,
              i.e. path with "thumbnail.tiff" replaced by "info.json") is
              read on demand to locate that sample's whole-slide image
              (info["image"]) and its ROI/contour file (info["roifile"],
              same {"0": {"points": [x,...]}, "1": {"points": [y,...]}}
              format written by contour_annotator.py, x/y normalized 0-1)
              -- used for the second, larger whole-slide popup below.
  canvas_width, canvas_height : pixel size of the plot canvas. Can differ
              (a non-square canvas is fine) -- the x/y data scale is kept
              equal in both directions regardless, so the shape of the
              embedding itself is never distorted; the plot is just
              letterboxed/centered within whichever canvas you ask for.
  wsi_max_dim : max pixel dimension for the whole-slide overview image
              (read from the slide's lowest-resolution pyramid level).

Behavior:
  - Hover a point -> tooltip with df_meta fields.
  - Click a point -> popup docked to the right of the canvas (doesn't cover
    the whole plot) showing the sample's thumbnail — sized proportional to
    the image's own aspect ratio, with room below it for labels — and its
    labels, with +/- controls to add/remove single-word/short-phrase labels.
    Clicking the thumbnail itself opens a second, larger popup with the
    whole slide it was cropped from (its lowest-resolution pyramid level)
    with that sample's saved contour overlaid (black line, glowing yellow
    outline). Esc closes whichever popup is topmost -- the whole-slide
    popup first if it's open, then the small popup on a second Esc.
  - Labels are persisted to <savepath>/labels.json ({sample: [label, ...]}),
    colors to <savepath>/colors.json ({label: "#hex"}) -- colors are
    assigned once per label and stay consistent for the whole dataset /
    across reloads. Points are rendered as pie-sliced dots (grey/unlabeled
    behind, colored/labeled on top), one slice per label color, sized by
    the on-canvas size slider.
  - If df_features is given (and scikit-learn is installed), the popup also
    shows likely-label suggestions with a one-click "+" to confirm each;
    confirming adds it as a real label and retrains that label's model.
"""

import json, logging, os, socket, threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from io import BytesIO
from urllib.parse import urlparse, unquote

import numpy as np
import pandas as pd
import tifffile
from PIL import Image
from IPython.display import display, HTML

try:
    from sklearn.linear_model import LogisticRegression as LR
except ImportError:
    LR = None

# tifffile reports "shaped series shape does not match page shape" (for tiffs
# whose embedded shape metadata doesn't match the actual page -- e.g. written
# via a downsampled slice without updating that metadata) through its own
# logger, not warnings.warn, so it has to be silenced via logging, not
# warnings.filterwarnings. Harmless here since we always use the real decoded
# array shape, never the metadata.
try:
    tifffile.logger().setLevel(logging.ERROR)
except Exception:
    logging.getLogger("tifffile").setLevel(logging.ERROR)

_PALETTE = [
    "#e94560", "#0f9d58", "#4285f4", "#f4b400", "#9c27b0", "#00bcd4",
    "#ff7043", "#8bc34a", "#3f51b5", "#e91e63", "#009688", "#795548",
    "#607d8b", "#cddc39", "#ff5722", "#673ab7",
]

_DEFAULT_CLF_PARAMS = {
    "penalty": "l2", "C": 10, "class_weight": "balanced",
    "solver": "liblinear", "max_iter": 1000,
}

# Minimum predicted probability for a label to be assigned automatically
# by the "Classify" button (applied to every sample with a trained classifier).
AUTO_LABEL_THRESHOLD = 0.9


def annotateScatter(
    df: pd.DataFrame,
    df_meta: pd.DataFrame | None = None,
    df_features: pd.DataFrame | None = None,
    image_paths: dict | None = None,
    savepath: str = "./labels/",
    canvas_width: int = 700,
    canvas_height: int = 700,
    thumb_max_dim: int = 260,
    wsi_max_dim: int = 1000,
    augFunc = None,
    qs = None,
    alpha = 0.85,
    clf_params: dict | None = None,
    min_label_count: int = 3,
):
    assert {"x", "y"}.issubset(df.columns), "df must have 'x' and 'y' columns"
    image_paths = image_paths or {}
    samples = df.index.astype(str).tolist()
    pts = df[["x", "y"]].astype(float).values.tolist()
    missing = [s for s in samples if s not in image_paths]
    if missing:
        print(f"[annotateScatter] no image path for {len(missing)} sample(s) "
              f"(e.g. {missing[:3]}) -- popup will show no thumbnail for those.")

    meta = {}
    for s in samples:
        if df_meta is not None and s in df_meta.index:
            row = df_meta.loc[s]
            meta[s] = {k: (None if pd.isna(v) else (v.item() if hasattr(v, "item") else v))
                       for k, v in row.items()}
        else:
            meta[s] = {}

    savepath = savepath.rstrip("/")
    os.makedirs(savepath, exist_ok=True)
    labels_path = os.path.join(savepath, "labels.json")
    colors_path = os.path.join(savepath, "colors.json")
    labels = json.load(open(labels_path)) if os.path.isfile(labels_path) else {}
    colors = json.load(open(colors_path)) if os.path.isfile(colors_path) else {}
    auto_labels: dict = {}   # in-memory only; populated by /classify, cleared by /clear_auto

    lock = threading.Lock()

    def _save():
        with open(labels_path, "w") as f:
            json.dump(labels, f, indent=2)
        with open(colors_path, "w") as f:
            json.dump(colors, f, indent=2)

    def _color_for(label):
        if label not in colors:
            colors[label] = _PALETTE[len(colors) % len(_PALETTE)]
        return colors[label]

    # ── optional per-label suggestion classifiers (needs df_features + sklearn) ──
    _clf_params = clf_params or _DEFAULT_CLF_PARAMS
    sample_idx = {s: i for i, s in enumerate(samples)}
    have_features = df_features is not None
    X_all = None
    valid_row = None
    if have_features and LR is None:
        print("[annotateScatter] scikit-learn not installed -- label suggestions disabled.")
        have_features = False
    if have_features:
        try:
            # feat = df_features.reindex(samples)
            X_all = df_features.astype(float).values
            valid_row = ~np.isnan(X_all).any(axis=1)
        except Exception as e:
            print(f"[annotateScatter] could not use df_features ({e}) -- label suggestions disabled.")
            have_features = False

    def _augment(X, y):
        """Data augmentation ahead of each label's LR fit. Take (X, y) and return an augmented (X_aug, y_aug)."""
        if (augFunc is not None) and (qs is not None):
            curated_positive = np.where(y == 1)[0]
            curated_negative = np.where(y == 0)[0]
            names = df_features.columns
            dpos_v = []
            for s_pos in curated_positive:
                s_neg = np.random.choice(curated_negative)
                df_pos = pd.Series(X[s_pos], index=names).unstack()
                df_neg = pd.Series(X[s_neg], index=names).unstack()
                assert df_pos.index.equals(df_neg.index)
                acdf = augFunc(df_pos.index.values, df_pos.values,
                            df_neg.values, alpha=alpha, beta=1. - alpha)
                acdf = pd.DataFrame(index=df_pos.index, columns=df_pos.columns, data=acdf)
                dpos_v.append(acdf.T.sort_index().T.stack().rename(s_pos))
            dpos = pd.concat(dpos_v, axis=1).T
            X_ = np.concatenate([dpos.values, X[curated_negative]], axis=0)
            y_ = np.array([1] * len(dpos) + [0] * len(curated_negative))
            return X_, y_
        return X, y

    clfs: dict[str, "LR"] = {}   # label -> fitted classifier (only labels with enough data)

    def _train_one(label):
        """(Re)train label's one-vs-rest classifier, using only labelled samples, or drop it if there's
        not enough support in both classes. Assumes caller holds `lock`."""
        if not have_features:
            return
        y_pos = np.array([1 if label in labels.get(s, []) else 0 for s in samples])
        y_neg = np.array([1 if (len(labels.get(s, []))>0) and (label not in labels.get(s, [])) else 0 for s in samples])
        if np.sum(y_pos) < min_label_count or np.sum(y_neg) < min_label_count:
            clfs.pop(label, None)
            return
        mask = (y_pos > 0) | (y_neg > 0)
        X_train, y_train = X_all[mask], y_pos[mask]
        X_train, y_train = _augment(X_train, y_train)
        clf = LR(**_clf_params)
        clf.fit(X_train, y_train)
        clfs[label] = clf

    def _retrain_label(label):
        if not have_features:
            return
        with lock:
            _train_one(label)

    def _predict_candidates(sample):
        """Probability of each label the sample doesn't already have, for
        every label with a trained classifier. Returns [] if features are
        unavailable for this sample or no classifiers exist yet."""
        if not have_features or sample not in sample_idx or not clfs:
            return []
        i = sample_idx[sample]
        if not valid_row[i]:
            return []
        x = X_all[i:i + 1]
        current = set(labels.get(sample, []))
        out = []
        for label, clf in clfs.items():
            if label in current:
                continue
            prob = float(clf.predict_proba(x)[0, 1])
            out.append({"label": label, "prob": prob})
        out.sort(key=lambda d: -d["prob"])
        return [d for d in out if d["prob"] >= 0.15][:8]

    auto_labels_path = os.path.join(savepath, "auto-labels.json")

    def _classify_all():
        """Run every trained classifier on every sample; assign labels whose
        predicted probability >= AUTO_LABEL_THRESHOLD.  Returns a new dict
        mapping sample -> [label, ...].  Caller must hold `lock`."""
        result: dict = {}
        if not have_features or not clfs:
            return result
        for i, sample in enumerate(samples):
            if not valid_row[i]:
                continue
            x = X_all[i:i + 1]
            assigned = []
            for label, clf in clfs.items():
                prob = float(clf.predict_proba(x)[0, 1])
                if prob >= AUTO_LABEL_THRESHOLD:
                    assigned.append(label)
            if assigned:
                result[sample] = assigned
        with open(auto_labels_path, "w") as f:
            json.dump(result, f, indent=2)
        return result

    if have_features:
        for _label in {l for labs in labels.values() for l in labs}:
            try:
                _retrain_label(_label)
            except Exception as e:
                print(f"[annotateScatter] initial training failed for label '{_label}': {e}")

    _thumb_cache: dict[str, dict] = {}   # sample -> {'png': bytes, 'w': int, 'h': int}

    def _load_thumb(sample):
        """Decode + resize once per sample, cache both the PNG bytes and the
        resulting (post-resize) pixel dimensions together, so /thumb and
        /dims never disagree and the image is only ever decoded once."""
        if sample in _thumb_cache:
            return _thumb_cache[sample]
        with lock:
            if sample in _thumb_cache:      # re-check: another thread may have just filled it
                return _thumb_cache[sample]
            path = image_paths[sample]
            with tifffile.TiffFile(path) as tif:
                arr = tif.pages[0].asarray()
            if arr.ndim == 2:
                arr = np.stack([arr] * 3, axis=-1)
            if arr.dtype != np.uint8:
                mx = arr.max()
                arr = ((arr / mx) * 255).astype(np.uint8) if mx > 0 else arr.astype(np.uint8)
            im = Image.fromarray(arr[..., :3])
            r = thumb_max_dim / max(im.size)
            im = im.resize((max(1, round(im.width * r)), max(1, round(im.height * r))))
            buf = BytesIO()
            im.save(buf, format="PNG")
            _thumb_cache[sample] = {"png": buf.getvalue(), "w": im.width, "h": im.height}
            return _thumb_cache[sample]

    _wsi_cache: dict[str, dict] = {}   # sample -> {'png','w','h','points'} | None

    def _info_path(sample):
        return image_paths[sample].replace("thumbnail.tiff", "info.json")

    def _load_wsi(sample):
        """Lazily (per-click) load the whole slide this sample's thumbnail
        was cropped from -- lowest-res pyramid level only -- plus its saved
        contour, via the sibling info.json {"image":..., "roifile":...}."""
        if sample in _wsi_cache:
            return _wsi_cache[sample]
        with lock:
            if sample in _wsi_cache:
                return _wsi_cache[sample]
            with open(_info_path(sample)) as f:
                info = json.load(f)
            with tifffile.TiffFile(info["image"]) as tf:
                n_levels = len(tf.series[0].levels)
                level_idx = n_levels - 1 if n_levels <= 4 else 3
                level = tf.series[0].levels[level_idx]
                arr = level.asarray()
            if arr.ndim == 2:
                arr = np.stack([arr] * 3, axis=-1)
            if arr.dtype != np.uint8:
                mx = arr.max()
                arr = ((arr / mx) * 255).astype(np.uint8) if mx > 0 else arr.astype(np.uint8)
            im = Image.fromarray(arr[..., :3])
            r = wsi_max_dim / max(im.size)
            im = im.resize((max(1, round(im.width * r)), max(1, round(im.height * r))))
            buf = BytesIO()
            im.save(buf, format="PNG")

            points = None
            roifile = info.get("roifile").replace('/dev-komp/data-200/', '/roi/data_batch_0/')
            if roifile and os.path.isfile(roifile):
                try:
                    with open(roifile) as f:
                        d = json.load(f)
                    xs, ys = d["0"]["points"], d["1"]["points"]
                    points = [[float(x), float(y)] for x, y in zip(xs, ys)]
                except Exception:
                    points = None

            _wsi_cache[sample] = {"png": buf.getvalue(), "w": im.width, "h": im.height, "points": points}
            return _wsi_cache[sample]
    # ── tiny HTTP server: GET thumbnails/wsi, POST label add/remove ────────
    with socket.socket() as s:
        s.bind(("", 0))
        port = s.getsockname()[1]

    class Handler(BaseHTTPRequestHandler):
        def do_OPTIONS(self):
            self._ok(b"")

        def do_GET(self):
            p = urlparse(self.path)
            if p.path.startswith("/thumb/"):
                sample = unquote(p.path[len("/thumb/"):])
                try:
                    body = _load_thumb(sample)["png"]
                except Exception:
                    self.send_response(404)
                    self.end_headers()
                    return
                self.send_response(200)
                self.send_header("Content-Type", "image/png")
                self.send_header("Content-Length", str(len(body)))
                self.send_header("Access-Control-Allow-Origin", "*")
                self.end_headers()
                self.wfile.write(body)
            elif p.path.startswith("/dims/"):
                sample = unquote(p.path[len("/dims/"):])
                try:
                    t = _load_thumb(sample)
                    body = json.dumps({"ok": True, "w": t["w"], "h": t["h"]}).encode()
                except Exception:
                    body = json.dumps({"ok": False, "w": None, "h": None}).encode()
                self.send_response(200)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(body)))
                self.send_header("Access-Control-Allow-Origin", "*")
                self.end_headers()
                self.wfile.write(body)
            elif p.path.startswith("/wsi/"):
                sample = unquote(p.path[len("/wsi/"):])
                try:
                    body = _load_wsi(sample)["png"]
                except Exception:
                    self.send_response(404)
                    self.end_headers()
                    return
                self.send_response(200)
                self.send_header("Content-Type", "image/png")
                self.send_header("Content-Length", str(len(body)))
                self.send_header("Access-Control-Allow-Origin", "*")
                self.end_headers()
                self.wfile.write(body)
            elif p.path.startswith("/wsidims/"):
                sample = unquote(p.path[len("/wsidims/"):])
                try:
                    t = _load_wsi(sample)
                    body = json.dumps({"ok": True, "w": t["w"], "h": t["h"], "points": t["points"]}).encode()
                except Exception:
                    body = json.dumps({"ok": False, "w": None, "h": None, "points": None}).encode()
                self.send_response(200)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(body)))
                self.send_header("Access-Control-Allow-Origin", "*")
                self.end_headers()
                self.wfile.write(body)
            elif p.path.startswith("/predict/"):
                sample = unquote(p.path[len("/predict/"):])
                try:
                    with lock:
                        candidates = _predict_candidates(sample)
                    body = json.dumps({"ok": True, "candidates": candidates}).encode()
                except Exception:
                    body = json.dumps({"ok": False, "candidates": []}).encode()
                self.send_response(200)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(body)))
                self.send_header("Access-Control-Allow-Origin", "*")
                self.end_headers()
                self.wfile.write(body)
            elif p.path == "/classify":
                try:
                    with lock:
                        auto_labels.clear()
                        auto_labels.update(_classify_all())
                        body = json.dumps({"ok": True, "auto_labels": auto_labels}).encode()
                except Exception as exc:
                    body = json.dumps({"ok": False, "auto_labels": {}, "error": str(exc)}).encode()
                self.send_response(200)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(body)))
                self.send_header("Access-Control-Allow-Origin", "*")
                self.end_headers()
                self.wfile.write(body)
            elif p.path == "/clear_auto":
                with lock:
                    auto_labels.clear()
                body = json.dumps({"ok": True, "auto_labels": {}}).encode()
                self.send_response(200)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(body)))
                self.send_header("Access-Control-Allow-Origin", "*")
                self.end_headers()
                self.wfile.write(body)
            else:
                self.send_response(404)
                self.end_headers()

        def do_POST(self):
            n = int(self.headers.get("Content-Length", 0))
            payload = json.loads(self.rfile.read(n) or b"{}")
            action = payload.get("action")
            with lock:
                if action == "add":
                    s, l = payload["sample"], payload["label"].strip()
                    if l:
                        labels.setdefault(s, [])
                        if l not in labels[s]:
                            labels[s].append(l)
                        _color_for(l)
                        _save()
                        _train_one(l)
                elif action == "remove":
                    s, l = payload["sample"], payload["label"]
                    if s in labels and l in labels[s]:
                        labels[s].remove(l)
                        if not labels[s]:
                            del labels[s]
                        _save()
                        _train_one(l)
                self._ok(json.dumps({"ok": True, "labels": labels, "colors": colors, "auto_labels": auto_labels}).encode())

        def _ok(self, body):
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Access-Control-Allow-Origin", "*")
            self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
            self.send_header("Access-Control-Allow-Headers", "Content-Type")
            self.end_headers()
            if body:
                self.wfile.write(body)

        def log_message(self, *a):
            pass

    host = socket.gethostbyname(socket.gethostname())
    srv = ThreadingHTTPServer((host, port), Handler)
    threading.Thread(target=srv.serve_forever, daemon=True).start()

    display(HTML(_render(samples, pts, meta, labels, colors, host, port, canvas_width, canvas_height, have_features, auto_labels)))


def _render(samples, pts, meta, labels, colors, host, port, canvas_width, canvas_height, have_features, auto_labels):
    css = """
<style>
  #ua-root * { box-sizing:border-box; margin:0; padding:0; }
  #ua-root {
    font-family:'JetBrains Mono','Fira Code',Consolas,monospace;
    background:#1a1a2e; color:#eaeaea; padding:14px; border-radius:10px;
    display:inline-block; user-select:none;
  }
  #ua-wrap { position:relative; }
  #ua-canvas { display:block; border:1.5px solid #0f3460; border-radius:6px; cursor:crosshair; background:#0f0f18; touch-action:none; }
  #ua-reset {
    position:absolute; top:8px; left:8px; z-index:40;
    font-family:inherit; font-size:10px; color:#8892a4; background:rgba(15,15,24,0.75);
    border:1px solid #2a2a4a; border-radius:4px; padding:2px 7px; cursor:pointer;
  }
  #ua-reset:hover { color:#eaeaea; border-color:#4a4a6a; }
  #ua-zoomlabel {
    position:absolute; bottom:8px; left:8px; z-index:40; pointer-events:none;
    font-size:10px; color:#8892a4; background:rgba(15,15,24,0.75);
    border:1px solid #2a2a4a; border-radius:4px; padding:2px 7px;
  }
  #ua-sizectl {
    position:absolute; bottom:8px; right:8px; z-index:40;
    display:flex; align-items:center; gap:6px;
    font-size:10px; color:#8892a4; background:rgba(15,15,24,0.75);
    border:1px solid #2a2a4a; border-radius:4px; padding:2px 8px;
  }
  #ua-sizectl input[type=range] { width:70px; accent-color:#e94560; cursor:pointer; }
  #ua-tooltip {
    position:fixed; pointer-events:none; display:none; z-index:2147483647;
    background:#f3f1ea; color:#111; font-size:11px;
    border:1px solid #ccc; border-radius:6px; padding:6px 8px;
    box-shadow:0 4px 12px rgba(0,0,0,0.35); max-width:260px;
  }
  #ua-popup {
    position:absolute; top:38px; right:10px; display:none;
    flex-direction:column; overflow:hidden;
    width:max-content; min-width:170px; max-width:420px;
    /*POPUP_MAX_HEIGHT*/
    background:rgba(20,20,34,0.97); border:1px solid #0f3460; border-radius:8px;
    padding:10px; z-index:50; box-shadow:0 6px 20px rgba(0,0,0,0.6);
  }
  #ua-popup-close { position:absolute; top:8px; right:10px; cursor:pointer; color:#8892a4; font-size:13px; z-index:2; }
  #ua-popup-title { flex:0 0 auto; font-size:12px; color:#e94560; margin-bottom:6px; padding-right:16px; word-break:break-all; }
  #ua-popup img {
    flex:0 0 auto; display:block; margin:0 auto 8px; background:#000; border-radius:4px;
    max-width:100%; max-height:380px; object-fit:contain;
  }
  #ua-label-list, #ua-add-row, #ua-suggest-list { width:100%; min-width:150px; }
  #ua-suggest-list { flex:0 0 auto; display:flex; flex-wrap:wrap; gap:4px; margin-bottom:8px; }
  .ua-suggest-chip {
    display:inline-flex; align-items:center; gap:5px; font-size:11px;
    background:rgba(233,69,96,0.12); border:1px dashed #e94560; border-radius:9px;
    padding:2px 4px 2px 8px; color:#eaeaea;
  }
  .ua-suggest-prob { color:#8892a4; font-size:10px; }
  .ua-suggest-add {
    cursor:pointer; color:#e94560; font-weight:bold; font-size:13px;
    padding:0 4px; line-height:1;
  }
  .ua-suggest-add:hover { color:#fff; background:#e94560; border-radius:50%; }
  #ua-label-list { flex:0 0 auto; display:flex; flex-direction:column; gap:4px; margin-bottom:8px; max-height:120px; overflow-y:auto; }
  .ua-label-row { display:flex; align-items:center; gap:6px; font-size:11px; }
  .ua-dot { width:9px; height:9px; border-radius:50%; flex:0 0 auto; }
  .ua-label-text { flex:1; overflow:hidden; text-overflow:ellipsis; white-space:nowrap; }
  .ua-rm { cursor:pointer; color:#8892a4; font-size:12px; }
  .ua-rm:hover { color:#e94560; }
  #ua-add-row { flex:0 0 auto; display:flex; gap:4px; }
  #ua-add-input {
    flex:1; font-family:inherit; font-size:11px; background:#0f0f18; color:#eaeaea;
    border:1px solid #2a2a4a; border-radius:4px; padding:4px 6px;
  }
  #ua-add-btn {
    font-family:inherit; font-size:12px; background:#e94560; color:#fff;
    border:none; border-radius:4px; padding:4px 10px; cursor:pointer;
  }
  #ua-add-btn:hover { filter:brightness(1.15); }
  #ua-auto-list { flex:0 0 auto; display:flex; flex-direction:column; gap:3px; margin-bottom:6px; max-height:100px; overflow-y:auto; }
  .ua-auto-row { display:flex; align-items:center; gap:6px; font-size:11px; }
  #ua-classify-ctl {
    position:absolute; top:8px; right:8px; z-index:40;
    display:flex; align-items:center; gap:6px;
    font-size:10px; background:rgba(15,15,24,0.85);
    border:1px solid #2a2a4a; border-radius:4px; padding:3px 8px;
  }
  #ua-classify-btn, #ua-clear-auto-btn {
    font-family:inherit; font-size:10px; color:#eaeaea;
    border:1px solid #2a2a4a; border-radius:3px; padding:2px 8px; cursor:pointer;
  }
  #ua-classify-btn { background:#1b3a5c; }
  #ua-classify-btn:hover { background:#0f3460; border-color:#4a7ab5; }
  #ua-clear-auto-btn { background:#2a1a2e; }
  #ua-clear-auto-btn:hover { background:#3a1a3e; border-color:#8a4a9a; }
  #ua-auto-status { font-size:10px; color:#8892a4; white-space:nowrap; }
  #ua-popup-meta {
    flex:1 1 auto; min-height:0; overflow-y:auto; margin-top:8px;
    background:#f3f1ea; color:#000; border-radius:6px; padding:6px 8px; font-size:11px;
  }
  #ua-popup img { cursor: zoom-in; }
  #ua-wsi-backdrop {
    position:fixed; inset:0; background:rgba(0,0,0,0.6); z-index:2147483645; display:none;
  }
  #ua-wsi-popup {
    position:fixed; top:50%; left:50%; transform:translate(-50%,-50%);
    display:none; max-width:92vw; max-height:92vh; overflow:auto;
    background:rgba(10,10,18,0.98); border:1px solid #0f3460; border-radius:8px;
    padding:14px; z-index:2147483646; box-shadow:0 10px 40px rgba(0,0,0,0.85);
  }
  #ua-wsi-close { float:right; cursor:pointer; color:#8892a4; font-size:15px; }
  #ua-wsi-close:hover { color:#eaeaea; }
  #ua-wsi-title { font-size:12px; color:#e94560; margin-bottom:8px; padding-right:18px; }
  #ua-wsi-canvas { display:block; max-width:100%; height:auto; border-radius:4px; background:#000; }
  #ua-wsi-status { font-size:11px; color:#8892a4; margin-top:6px; }
</style>
<div id="ua-root">
  <div id="ua-wrap">
    <canvas id="ua-canvas"></canvas>
    <div id="ua-reset" title="reset zoom/pan">⤾ reset view</div>
    <div id="ua-zoomlabel">100%</div>
    <div id="ua-sizectl" title="dot/pie size">
      <span>size</span>
      <input type="range" id="ua-size-slider" min="2" max="14" step="1" value="5" />
    </div>
    <div id="ua-classify-ctl">
      <button id="ua-classify-btn">Classify</button>
      <button id="ua-clear-auto-btn">Clear</button>
      <span id="ua-auto-status"></span>
    </div>
    <div id="ua-popup">
      <span id="ua-popup-close">✕</span>
      <div id="ua-popup-title">—</div>
      <img id="ua-popup-img" src="" title="click to view whole slide" />
      <div id="ua-label-list"></div>
      <div id="ua-suggest-list"></div>
      <div id="ua-add-row">
        <input id="ua-add-input" placeholder="add label…" />
        <button id="ua-add-btn">+</button>
      </div>
      <div id="ua-auto-list"></div>
      <div id="ua-popup-meta"></div>
    </div>
  </div>
</div>
<div id="ua-tooltip"></div>
<div id="ua-wsi-backdrop"></div>
<div id="ua-wsi-popup">
  <span id="ua-wsi-close">✕</span>
  <div id="ua-wsi-title">—</div>
  <canvas id="ua-wsi-canvas"></canvas>
  <div id="ua-wsi-status"></div>
</div>
"""
    css = css.replace("/*POPUP_MAX_HEIGHT*/", f"max-height:{max(200, canvas_height - 20)}px;")
    js = (
        "<script>\n(function(){\n"
        f"  const SAMPLES = {json.dumps(samples)};\n"
        f"  const PTS     = {json.dumps(pts)};\n"
        f"  const META    = {json.dumps(meta)};\n"
        f"  let   LABELS     = {json.dumps(labels)};\n"
        f"  let   COLORS     = {json.dumps(colors)};\n"
        f"  let   AUTO_LABELS = {json.dumps(auto_labels)};\n"
        f"  const CW         = {canvas_width};\n"
        f"  const CH         = {canvas_height};\n"
        f"  const SRV        = 'http://{host}:{port}';\n"
        f"  const HAS_FEATURES = {json.dumps(have_features)};\n"
        + r"""
  const canvas  = document.getElementById('ua-canvas');
  const ctx     = canvas.getContext('2d');
  const tooltip = document.getElementById('ua-tooltip');
  const popup       = document.getElementById('ua-popup');
  const popupTitle  = document.getElementById('ua-popup-title');
  const popupImg    = document.getElementById('ua-popup-img');
  const labelList   = document.getElementById('ua-label-list');
  const suggestList  = document.getElementById('ua-suggest-list');
  const popupMeta    = document.getElementById('ua-popup-meta');
  const addInput     = document.getElementById('ua-add-input');
  const autoList     = document.getElementById('ua-auto-list');
  const classifyBtn  = document.getElementById('ua-classify-btn');
  const clearAutoBtn = document.getElementById('ua-clear-auto-btn');
  const autoStatus   = document.getElementById('ua-auto-status');
  const resetBtn     = document.getElementById('ua-reset');
  const zoomLabel   = document.getElementById('ua-zoomlabel');
  const sizeSlider  = document.getElementById('ua-size-slider');
  popupImg.addEventListener('error', () => { popupImg.style.display = 'none'; });
  const wsiBackdrop = document.getElementById('ua-wsi-backdrop');
  const wsiPopup    = document.getElementById('ua-wsi-popup');
  const wsiTitle    = document.getElementById('ua-wsi-title');
  const wsiCanvas   = document.getElementById('ua-wsi-canvas');
  const wsiStatus   = document.getElementById('ua-wsi-status');
  const wsiClose    = document.getElementById('ua-wsi-close');
  const GRAY = '#555a6e';
  let baseR = parseFloat(sizeSlider.value);
  const HIT = 8;

  // ── zoom/pan transform: screen = base*scale + offset ───────────────
  let vscale = 1, voffX = 0, voffY = 0;
  const MIN_SCALE = 0.2, MAX_SCALE = 40;
  const toScreen = (bx, by) => [bx*vscale + voffX, by*vscale + voffY];
  const toBase   = (sx, sy) => [(sx - voffX)/vscale, (sy - voffY)/vscale];

  canvas.width = CW; canvas.height = CH;

  // ── fit points to canvas: equal x/y scale (never distort the embedding),
  //    centered within whatever CW x CH the canvas turns out to be ─────────
  const PAD = 24;
  const xs = PTS.map(p => p[0]), ys = PTS.map(p => p[1]);
  const xmin = Math.min(...xs), xmax = Math.max(...xs);
  const ymin = Math.min(...ys), ymax = Math.max(...ys);
  const xrange = Math.max(xmax - xmin, 1e-9), yrange = Math.max(ymax - ymin, 1e-9);
  const scale = Math.min((CW - 2*PAD) / xrange, (CH - 2*PAD) / yrange);
  const offX = (CW - xrange * scale) / 2;
  const offY = (CH - yrange * scale) / 2;
  const px = x => offX + (x - xmin) * scale;
  const py = y => (CH - offY) - (y - ymin) * scale;   // flip so +y is up
  const PX = PTS.map(p => [px(p[0]), py(p[1])]);

  let hoverIdx = -1, selectedIdx = -1, popupSeq = 0;

  function draw() {
    ctx.clearRect(0, 0, CW, CH);
    // radius scales gently with zoom so points stay legible but don't
    // balloon at high zoom; user's slider sets the base size.
    const rBase = Math.min(Math.max(baseR * Math.sqrt(vscale), 2), 40);

    // When any auto-labels exist, dots are colored by auto-labels only;
    // samples without an auto-label are grey (even if manually labeled).
    // When no auto-labels exist at all, fall back to manual labels.
    const hasAnyAuto = Object.keys(AUTO_LABELS).length > 0;

    function drawPoint(i) {
      const [cx, cy] = toScreen(PX[i][0], PX[i][1]);
      if (cx < -20 || cx > CW+20 || cy < -20 || cy > CH+20) return;  // cull offscreen
      const labs = hasAnyAuto ? (AUTO_LABELS[SAMPLES[i]] || []) : (LABELS[SAMPLES[i]] || []);
      const rad = (i === hoverIdx || i === selectedIdx) ? rBase + 2 : rBase;
      if (labs.length === 0) {
        ctx.beginPath(); ctx.arc(cx, cy, rad, 0, 2*Math.PI);
        ctx.fillStyle = GRAY; ctx.fill();
      } else {
        const step = 2*Math.PI / labs.length;
        labs.forEach((l, j) => {
          ctx.beginPath(); ctx.moveTo(cx, cy);
          ctx.arc(cx, cy, rad, j*step, (j+1)*step);
          ctx.closePath(); ctx.fillStyle = COLORS[l] || GRAY; ctx.fill();
        });
      }
      if (i === selectedIdx) {
        ctx.beginPath(); ctx.arc(cx, cy, rad + 2.5, 0, 2*Math.PI);
        ctx.strokeStyle = '#ffff00'; ctx.lineWidth = 1.5; ctx.stroke();
      }
    }

    // unlabeled/grey points render first (behind), colored/labeled points
    // render second so they're always visually on top.
    const activeLabelMap = hasAnyAuto ? AUTO_LABELS : LABELS;
    for (let i = 0; i < SAMPLES.length; i++) {
      if ((activeLabelMap[SAMPLES[i]] || []).length === 0) drawPoint(i);
    }
    for (let i = 0; i < SAMPLES.length; i++) {
      if ((activeLabelMap[SAMPLES[i]] || []).length > 0) drawPoint(i);
    }
    zoomLabel.textContent = Math.round(vscale * 100) + '%';
  }

  function nearest(mx, my) {
    let best = -1, bestD = HIT;
    for (let i = 0; i < PX.length; i++) {
      const [cx, cy] = toScreen(PX[i][0], PX[i][1]);
      const d = Math.hypot(cx - mx, cy - my);
      if (d < bestD) { bestD = d; best = i; }
    }
    return best;
  }

  function canvasXY(e) {
    const r = canvas.getBoundingClientRect();
    return [e.clientX - r.left, e.clientY - r.top];
  }

  // ── zoom (wheel, centered on cursor) ─────────────────────────────────
  canvas.addEventListener('wheel', e => {
    e.preventDefault();
    const [mx, my] = canvasXY(e);
    const factor = Math.pow(1.0016, -e.deltaY);
    const newScale = Math.min(Math.max(vscale * factor, MIN_SCALE), MAX_SCALE);
    const [bx, by] = toBase(mx, my);
    voffX = mx - bx * newScale;
    voffY = my - by * newScale;
    vscale = newScale;
    tooltip.style.display = 'none';
    draw();
  }, { passive: false });

  // ── pan (drag) ────────────────────────────────────────────────────────
  let panStart = null, panOrigin = null, dragging = false, dragged = false;
  canvas.addEventListener('mousedown', e => {
    if (e.button !== 0) return;
    panStart = canvasXY(e); panOrigin = [voffX, voffY]; dragging = true; dragged = false;
  });
  window.addEventListener('mousemove', e => {
    if (!dragging) return;
    const [mx, my] = canvasXY(e);
    const dx = mx - panStart[0], dy = my - panStart[1];
    if (!dragged && Math.hypot(dx, dy) > 3) dragged = true;
    if (dragged) {
      voffX = panOrigin[0] + dx; voffY = panOrigin[1] + dy;
      tooltip.style.display = 'none';
      draw();
    }
  });
  window.addEventListener('mouseup', () => { dragging = false; });

  resetBtn.addEventListener('click', () => {
    vscale = 1; voffX = 0; voffY = 0; draw();
  });
  sizeSlider.addEventListener('input', () => {
    baseR = parseFloat(sizeSlider.value);
    draw();
  });
  window.addEventListener('keydown', e => {
    if (e.key === 'Escape') {
      if (wsiPopup.style.display === 'block') closeWsiPopup();
      else closePopup();
    }
  });

  // ── shared HTML builders (tooltip + popup both use these) ───────────
  function metaTableHtml(sample) {
    const m = META[sample] || {};
    const rows = Object.entries(m).map(([k, v]) =>
      `<tr><td style="color:#000;padding:1px 8px 1px 0;font-weight:bold;">${k}</td>
           <td style="color:#000;padding:1px 0;">${v}</td></tr>`).join('');
    return rows ? `<table style="border-collapse:collapse;">${rows}</table>` : '';
  }

  function labelChipsHtml(sample) {
    const labs = LABELS[sample] || [];
    if (!labs.length) return '';
    const chips = labs.map(l =>
      `<span style="display:inline-flex;align-items:center;gap:4px;background:rgba(0,0,0,0.06);
                    padding:1px 6px;border-radius:9px;margin:2px 4px 0 0;">
         <span style="width:7px;height:7px;border-radius:50%;background:${COLORS[l] || GRAY};display:inline-block;"></span>${l}
       </span>`).join('');
    return `<div style="margin-top:4px;">${chips}</div>`;
  }

  function tooltipLabelsHtml(sample) {
    const manual = LABELS[sample] || [];
    const auto   = AUTO_LABELS[sample] || [];
    if (!manual.length && !auto.length) return '';
    const chipStyle = 'display:inline-flex;align-items:center;gap:4px;background:rgba(0,0,0,0.06);padding:1px 6px;border-radius:9px;margin:2px 3px 0 0;';
    const dot = col => `<span style="width:7px;height:7px;border-radius:50%;background:${col};display:inline-block;"></span>`;
    let html = '';
    if (manual.length) {
      html += `<div style="margin-top:5px;font-size:10px;color:#777;letter-spacing:.04em;">MANUAL</div>`;
      html += `<div style="margin-top:2px;">${manual.map(l =>
        `<span style="${chipStyle}">${dot(COLORS[l] || GRAY)}${l}</span>`).join('')}</div>`;
    }
    if (auto.length) {
      html += `<div style="margin-top:5px;font-size:10px;color:#5a8aaa;letter-spacing:.04em;">AUTO</div>`;
      html += `<div style="margin-top:2px;">${auto.map(l =>
        `<span style="${chipStyle}">${dot(COLORS[l] || GRAY)}${l}</span>`).join('')}</div>`;
    }
    return html;
  }

  // ── hover tooltip ────────────────────────────────────────────────────
  canvas.addEventListener('mousemove', e => {
    if (dragging && dragged) return;   // panning takes priority over hover
    const [mx, my] = canvasXY(e);
    const idx = nearest(mx, my);
    if (idx !== hoverIdx) { hoverIdx = idx; draw(); }
    if (idx < 0) { tooltip.style.display = 'none'; return; }
    const sample = SAMPLES[idx];
    const meta = metaTableHtml(sample);
    tooltip.innerHTML = `<b style="color:#e94560;">${sample}</b>` +
      tooltipLabelsHtml(sample) +
      (meta ? `<div style="margin-top:4px;">${meta}</div>` : '');
    tooltip.style.display = 'block';
    tooltip.style.left = (e.clientX + 14) + 'px';
    tooltip.style.top  = (e.clientY + 14) + 'px';
  });
  canvas.addEventListener('mouseleave', () => { hoverIdx = -1; tooltip.style.display = 'none'; draw(); });

  // ── click -> popup ───────────────────────────────────────────────────
  function renderLabelList(sample) {
    const labs = LABELS[sample] || [];
    labelList.innerHTML = labs.map(l => `
      <div class="ua-label-row" data-label="${l}">
        <span class="ua-dot" style="background:${COLORS[l] || GRAY}"></span>
        <span class="ua-label-text" title="${l}">${l}</span>
        <span class="ua-rm" title="remove">✕</span>
      </div>`).join('') || '<div style="color:#666;font-size:11px;">no labels yet</div>';
    labelList.querySelectorAll('.ua-rm').forEach(el => {
      el.addEventListener('click', () => {
        const label = el.parentElement.dataset.label;
        postAction({ action: 'remove', sample, label });
      });
    });
  }

  const dimsCache = {};  // sample -> [w,h] | null, filled in lazily on demand

  function sizeImageTo(dims) {
    if (!dims) { popupImg.style.width = '200px'; popupImg.style.aspectRatio = '4 / 3'; return; }
    const [w, h] = dims;
    popupImg.style.width = Math.min(w, 380) + 'px';   // popup (max-content) wraps to this
    popupImg.style.aspectRatio = w + ' / ' + h;        // height follows automatically, true to source
  }

  function updateAutoStatus() {
    const total = Object.keys(AUTO_LABELS).length;
    autoStatus.textContent = total > 0
      ? `auto: ${Object.values(AUTO_LABELS).reduce((s,v)=>s+v.length,0)} in ${total} samples`
      : '';
  }

  function renderAutoList(sample) {
    const labs = AUTO_LABELS[sample] || [];
    if (!labs.length) { autoList.innerHTML = ''; return; }
    const header = '<div style="font-size:10px;color:#5a8aaa;letter-spacing:.04em;margin-bottom:3px;">AUTO</div>';
    autoList.innerHTML = header + labs.map(l =>
      `<div class="ua-auto-row">
         <span class="ua-dot" style="background:${COLORS[l] || GRAY}"></span>
         <span class="ua-label-text" title="${l}">${l}</span>
       </div>`).join('');
    updateAutoStatus();
  }

  function renderMeta(sample) {
    const html = metaTableHtml(sample);
    popupMeta.innerHTML = html || '<div style="color:#666;">no metadata</div>';
  }

  function renderSuggestions(sample) {
    suggestList.innerHTML = '';
    if (!HAS_FEATURES) return;
    fetch(SRV + '/predict/' + encodeURIComponent(sample))
      .then(r => r.json())
      .then(j => {
        if (selectedIdx < 0 || SAMPLES[selectedIdx] !== sample) return;   // popup moved on
        if (!j.ok || !j.candidates.length) return;
        suggestList.innerHTML = j.candidates.map(c => `
          <span class="ua-suggest-chip" data-label="${c.label}">
            ${c.label} <span class="ua-suggest-prob">${Math.round(c.prob * 100)}%</span>
            <span class="ua-suggest-add" title="confirm suggestion">+</span>
          </span>`).join('');
        suggestList.querySelectorAll('.ua-suggest-add').forEach(el => {
          el.addEventListener('click', () => {
            const chip = el.closest('.ua-suggest-chip');
            const label = chip.dataset.label;
            postAction({ action: 'add', sample, label });
            chip.remove();   // it's now a real label, shown in the label list instead
          });
        });
      })
      .catch(() => {});
  }

  function openPopup(idx) {
    selectedIdx = idx;
    const sample = SAMPLES[idx];
    const mySeq = ++popupSeq;
    popupTitle.textContent = sample;
    sizeImageTo(dimsCache[sample]);   // reasonable placeholder immediately; refined once real dims arrive
    popupImg.style.display = 'block';
    popupImg.src = SRV + '/thumb/' + encodeURIComponent(sample);  // loads in parallel
    renderLabelList(sample);
    renderSuggestions(sample);
    renderAutoList(sample);
    renderMeta(sample);
    addInput.value = '';
    popup.style.display = 'flex';
    draw();

    if (!(sample in dimsCache)) {
      fetch(SRV + '/dims/' + encodeURIComponent(sample))
        .then(r => r.json())
        .then(j => {
          dimsCache[sample] = j.ok ? [j.w, j.h] : null;
          if (mySeq !== popupSeq) return;              // popup moved on to another sample
          sizeImageTo(dimsCache[sample]);
        })
        .catch(() => { dimsCache[sample] = null; });
    }
  }

  function closePopup() {
    selectedIdx = -1; popup.style.display = 'none';
    closeWsiPopup();
    draw();
  }

  let wsiSeq = 0;

  function openWsiPopup() {
    if (selectedIdx < 0) return;
    const sample = SAMPLES[selectedIdx];
    const mySeq = ++wsiSeq;
    wsiTitle.textContent = sample + ' — whole slide';
    wsiStatus.textContent = 'Loading…';
    wsiCanvas.width = 400; wsiCanvas.height = 300;   // placeholder box while fetching
    wsiCanvas.getContext('2d').clearRect(0, 0, 400, 300);
    wsiBackdrop.style.display = 'block';
    wsiPopup.style.display = 'block';

    fetch(SRV + '/wsidims/' + encodeURIComponent(sample))
      .then(r => r.json())
      .then(j => {
        if (mySeq !== wsiSeq) return;
        if (!j.ok) { wsiStatus.textContent = '⚠ No whole-slide image available for this sample.'; return; }
        wsiCanvas.width = j.w; wsiCanvas.height = j.h;
        const img = new Image();
        img.onload = () => {
          if (mySeq !== wsiSeq) return;
          const c = wsiCanvas.getContext('2d');
          c.clearRect(0, 0, j.w, j.h);
          c.drawImage(img, 0, 0, j.w, j.h);
          if (j.points && j.points.length > 1) {
            const pts = j.points.map(([nx, ny]) => [nx * j.w, ny * j.h]);
            const tracePath = () => {
              c.beginPath();
              c.moveTo(pts[0][0], pts[0][1]);
              pts.slice(1).forEach(p => c.lineTo(p[0], p[1]));
              c.closePath();
            };
            // thick glowing yellow outline, then a crisp black line on top
            c.save();
            c.shadowColor = '#ffff00'; c.shadowBlur = 18;
            tracePath(); c.strokeStyle = '#ffff00'; c.lineWidth = 6; c.stroke();
            c.restore();
            tracePath(); c.strokeStyle = '#000'; c.lineWidth = 2.5; c.stroke();
            wsiStatus.textContent = '';
          } else {
            wsiStatus.textContent = '(no saved contour for this sample)';
          }
        };
        img.onerror = () => { wsiStatus.textContent = '⚠ Whole-slide image failed to load.'; };
        img.src = SRV + '/wsi/' + encodeURIComponent(sample);
      })
      .catch(() => { wsiStatus.textContent = '⚠ Request failed.'; });
  }

  function closeWsiPopup() {
    wsiSeq++;   // invalidate any in-flight load
    wsiBackdrop.style.display = 'none';
    wsiPopup.style.display = 'none';
  }

  popupImg.addEventListener('click', openWsiPopup);
  wsiClose.addEventListener('click', closeWsiPopup);
  wsiBackdrop.addEventListener('click', closeWsiPopup);

  async function postAction(body) {
    try {
      const resp = await fetch(SRV, {
        method: 'POST', headers: {'Content-Type': 'application/json'},
        body: JSON.stringify(body)
      });
      const json = await resp.json();
      if (json.ok) {
        LABELS = json.labels; COLORS = json.colors;
        if (json.auto_labels !== undefined) AUTO_LABELS = json.auto_labels;
        if (selectedIdx >= 0) {
          renderLabelList(SAMPLES[selectedIdx]);
          renderAutoList(SAMPLES[selectedIdx]);
        }
        draw();
      }
    } catch (e) { /* best effort */ }
  }

  canvas.addEventListener('click', e => {
    if (dragged) { dragged = false; return; }   // this click ended a pan, not a selection
    const [mx, my] = canvasXY(e);
    const idx = nearest(mx, my);
    if (idx >= 0) openPopup(idx); else closePopup();
  });
  document.getElementById('ua-popup-close').addEventListener('click', closePopup);
  document.getElementById('ua-add-btn').addEventListener('click', () => {
    if (selectedIdx < 0) return;
    const label = addInput.value.trim();
    if (!label) return;
    postAction({ action: 'add', sample: SAMPLES[selectedIdx], label });
    addInput.value = '';
  });
  addInput.addEventListener('keydown', e => {
    if (e.key === 'Enter') document.getElementById('ua-add-btn').click();
  });

  classifyBtn.addEventListener('click', async () => {
    classifyBtn.disabled = true;
    classifyBtn.textContent = 'Classifying…';
    autoStatus.textContent = 'running…';
    try {
      const resp = await fetch(SRV + '/classify');
      const json = await resp.json();
      if (json.ok) {
        AUTO_LABELS = json.auto_labels;
        updateAutoStatus();
        if (!Object.keys(AUTO_LABELS).length)
          autoStatus.textContent = 'no classifiers — add labels first';
        if (selectedIdx >= 0) renderAutoList(SAMPLES[selectedIdx]);
        draw();
      }
    } catch (e) { autoStatus.textContent = '⚠ failed'; }
    classifyBtn.disabled = false;
    classifyBtn.textContent = 'Classify';
  });

  clearAutoBtn.addEventListener('click', async () => {
    try {
      const resp = await fetch(SRV + '/clear_auto');
      const json = await resp.json();
      if (json.ok) {
        AUTO_LABELS = {};
        autoStatus.textContent = '';
        if (selectedIdx >= 0) renderAutoList(SAMPLES[selectedIdx]);
        else autoList.innerHTML = '';
        draw();
      }
    } catch (e) { /* best effort */ }
  });

  draw();
})();
</script>
"""
    )
    return css + js