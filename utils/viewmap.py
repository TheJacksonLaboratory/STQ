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
  - Hovering over the whole-slide popup shows a circular magnifying-glass
    (512px diameter) that fetches and displays a live 512x512 full-resolution
    tile from the *original* whole-slide image (not the downsampled overview),
    centered on the hovered location. Requires zarr (used via tifffile's
    aszarr() for lazy, region-only reads of the pyramidal TIFF); the
    magnifier is silently disabled (with a printed note) if zarr is missing.
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
from urllib.parse import urlparse, unquote, parse_qs

import numpy as np
import pandas as pd
import tifffile
from PIL import Image
from IPython.display import display, HTML
from matplotlib import pyplot as plt

try:
    from sklearn.linear_model import LogisticRegression as LR
except ImportError:
    LR = None

try:
    import zarr
except ImportError:
    zarr = None

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
AUTO_LABEL_THRESHOLD = 0.75

# Side length (px) of the live full-resolution tile fetched under the
# whole-slide magnifying glass. Kept square; magnifier circle diameter is 512, need to downsample.
TILE_SIZE = 1024


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

    have_magnifier = zarr is not None
    if not have_magnifier:
        print("[annotateScatter] zarr not installed -- whole-slide magnifier disabled.")


    def _augment(X, y):
        """Data augmentation ahead of each label's LR fit. Take (X, y) and return an augmented (X_aug, y_aug)."""
        if (augFunc is not None) and (qs is not None):
            curated_positive = np.where(y == 1)[0]
            curated_negative = np.where(y == 0)[0]
            X_ = X.copy()
            for s_pos in X.index[curated_positive]:
                s_neg = X.index[np.random.choice(curated_negative)]
                df_pos = X.loc[s_pos].unstack()
                df_neg = X.loc[s_neg].unstack()
                assert df_pos.index.equals(df_neg.index)
                acdf = augFunc(df_pos.index.values, df_pos.values,
                            df_neg.values, alpha=alpha, beta=1. - alpha)
                acdf = pd.DataFrame(index=df_pos.index, columns=df_pos.columns, data=acdf)
                acdf = acdf.T.sort_index().T.stack().rename(s_pos)
                X_.loc[s_pos] = acdf
            return X_, y
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
        X_train, y_train = df_features.loc[mask], y_pos[mask]
        X_train, y_train = _augment(X_train, y_train)
        clf = LR(**_clf_params)
        clf.fit(X_train.values, y_train)
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
        predicted probability >= AUTO_LABEL_THRESHOLD. Returns a new dict
        mapping sample -> [label, ...]. Caller must hold `lock`."""
        result: dict = {}
        if not have_features or not clfs:
            return result

        valid_idx = [i for i, ok in enumerate(valid_row) if ok]
        if not valid_idx:
            return result

        valid_samples = [samples[i] for i in valid_idx]
        X_valid = X_all[valid_idx]  # shape: (n_valid, n_features)

        # One batched predict_proba call per classifier instead of one per row
        label_probs = {}
        for label, clf in clfs.items():
            probs = clf.predict_proba(X_valid)[:, 1]  # shape: (n_valid,)
            label_probs[label] = probs

        # Assemble per-sample assigned labels
        assigned_lists = [[] for _ in valid_idx]
        for label, probs in label_probs.items():
            mask = probs >= AUTO_LABEL_THRESHOLD
            for row_pos in np.nonzero(mask)[0]:
                assigned_lists[row_pos].append(label)

        for sample, assigned in zip(valid_samples, assigned_lists):
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

    # ── live full-resolution tiles for the whole-slide magnifying glass ────
    _full_wsi_cache: dict[str, dict] = {}   # sample -> {'tf','z','h','w'} | None (open handles, kept alive)

    def _get_level0(sample):
        """Open (once) and cache a lazy zarr view onto the *full-resolution*
        (level 0) plane of the sample's whole-slide TIFF, so a single tile
        can be read without decoding the whole pyramid level. Returns None
        (and caches that) if unavailable, so repeated failed lookups are cheap."""
        if sample in _full_wsi_cache:
            return _full_wsi_cache[sample]
        with lock:
            if sample in _full_wsi_cache:
                return _full_wsi_cache[sample]
            try:
                with open(_info_path(sample)) as f:
                    info = json.load(f)

                store = tifffile.imread(info["image"], aszarr=True)
                z = zarr.open(store, mode='r')["0"]
                h, w = z.shape[0], z.shape[1]
                entry = {"z": z, "h": h, "w": w}
            except Exception as e:
                print(f"[annotateScatter] could not open full-res slide for '{sample}', {info['image']} ({e}) -- magnifier disabled for it.")
                entry = None
            _full_wsi_cache[sample] = entry
            return entry

    def _get_tile(sample, nx, ny):
        """Return a TILE_SIZE x TILE_SIZE PNG crop of the sample's *full
        resolution* whole-slide image, centered on the normalized (0-1, 0-1)
        location (nx, ny) within the slide -- same normalization the saved
        contour points use, so it lines up with the low-res overview shown
        in the whole-slide popup. Crops near an edge are clamped in-bounds
        (not padded), so the returned tile can be slightly off-center there."""
        entry = _get_level0(sample)
        if entry is None:
            return None
        z, h, w = entry["z"], entry["h"], entry["w"]
        cx, cy = nx * w, ny * h
        half = TILE_SIZE // 2
        x0 = int(np.clip(cx - half, 0, max(w - TILE_SIZE, 0)))
        y0 = int(np.clip(cy - half, 0, max(h - TILE_SIZE, 0)))
        x1, y1 = min(x0 + TILE_SIZE, w), min(y0 + TILE_SIZE, h)
        with lock:
            arr = np.asarray(z[y0:y1, x0:x1])
        if arr.ndim == 2:
            arr = np.stack([arr] * 3, axis=-1)
        arr = arr[..., :3]
        if arr.dtype != np.uint8:
            mx = arr.max()
            arr = ((arr / mx) * 255).astype(np.uint8) if mx > 0 else arr.astype(np.uint8)
        if arr.shape[0] != TILE_SIZE or arr.shape[1] != TILE_SIZE:
            canvas = np.zeros((TILE_SIZE, TILE_SIZE, 3), dtype=np.uint8)
            canvas[:arr.shape[0], :arr.shape[1]] = arr
            arr = canvas
        im = Image.fromarray(arr)
        buf = BytesIO()
        im.save(buf, format="PNG")
        return buf.getvalue()

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
            elif p.path.startswith("/tile/"):
                sample = unquote(p.path[len("/tile/"):])
                if not have_magnifier:
                    self.send_response(404)
                    self.end_headers()
                    return
                qsp = parse_qs(p.query)
                try:
                    nx = float(qsp.get("nx", ["0.5"])[0])
                    ny = float(qsp.get("ny", ["0.5"])[0])
                    body = _get_tile(sample, nx, ny)
                    if body is None:
                        raise ValueError("tile unavailable")
                except Exception:
                    self.send_response(404)
                    self.end_headers()
                    return
                self.send_response(200)
                self.send_header("Content-Type", "image/png")
                self.send_header("Content-Length", str(len(body)))
                self.send_header("Access-Control-Allow-Origin", "*")
                self.send_header("Cache-Control", "no-store")
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

    display(HTML(_render(samples, pts, meta, labels, colors, host, port, canvas_width, canvas_height, have_features, auto_labels, have_magnifier)))


def _render(samples, pts, meta, labels, colors, host, port, canvas_width, canvas_height, have_features, auto_labels, have_magnifier):
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
    padding:2px 8px; color:#eaeaea; cursor:pointer;
  }
  .ua-suggest-chip:hover { background:#e94560; color:#fff; border-style:solid; }
  .ua-suggest-prob { color:#8892a4; font-size:10px; }
  .ua-suggest-chip:hover .ua-suggest-prob { color:#fff; }
  #ua-label-list { flex:0 0 auto; display:flex; flex-direction:column; gap:4px; margin-bottom:8px; max-height:120px; overflow-y:auto; }
  .ua-label-row { display:flex; align-items:center; gap:6px; font-size:11px; }
  .ua-dot { width:9px; height:9px; border-radius:50%; flex:0 0 auto; }
  .ua-label-text { flex:1; overflow:hidden; text-overflow:ellipsis; white-space:nowrap; }
  .ua-rm { cursor:pointer; color:#8892a4; font-size:12px; }
  .ua-rm:hover { color:#e94560; }
  #ua-add-row { flex:0 0 auto; display:flex; gap:4px; position:relative; }
  #ua-add-input {
    flex:1; font-family:inherit; font-size:11px; background:#0f0f18; color:#eaeaea;
    border:1px solid #2a2a4a; border-radius:4px; padding:4px 6px;
  }
  #ua-add-btn {
    font-family:inherit; font-size:12px; background:#e94560; color:#fff;
    border:none; border-radius:4px; padding:4px 10px; cursor:pointer;
  }
  #ua-add-btn:hover { filter:brightness(1.15); }
  #ua-add-dropdown {
    position:absolute; top:100%; left:0; right:26px; margin-top:2px; display:none;
    max-height:120px; overflow-y:auto; z-index:60;
    background:#0f0f18; border:1px solid #2a2a4a; border-radius:4px;
    box-shadow:0 4px 12px rgba(0,0,0,0.5);
  }
  .ua-add-dropdown-item {
    display:flex; align-items:center; gap:6px; font-size:11px; padding:4px 7px; cursor:pointer;
  }
  .ua-add-dropdown-item:hover, .ua-add-dropdown-item.active { background:#1b3a5c; }
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
  #ua-wsi-canvas { display:block; max-width:100%; height:auto; border-radius:4px; background:#000; cursor:none; }
  #ua-wsi-status { font-size:11px; color:#8892a4; margin-top:6px; }
  #ua-magnifier {
    position:fixed; width:512px; height:512px; border-radius:50%;
    border:0px solid #000000; box-shadow:0 8px 28px rgba(0,0,0,0.7), 0 0 0 1px rgba(255,255,255,0.08) inset;
    pointer-events:none; display:none; z-index:2147483647; overflow:hidden; background:#000;
  }
  #ua-magnifier canvas { display:block; width:512px; height:512px; image-rendering:pixelated; }
  #ua-magnifier-status {
    position:absolute; inset:0; display:flex; align-items:center; justify-content:center;
    color:#8892a4; font-size:11px; text-align:center; padding:0 20px;
  }
  #ua-wsi-crosshair {
    position:fixed; width:26px; height:26px; left:0; top:0;
    pointer-events:none; display:none; z-index:2147483647;
  }
  #ua-wsi-crosshair::before, #ua-wsi-crosshair::after {
    content:''; position:absolute; background:#ffff00; box-shadow:0 0 3px 0 rgba(0,0,0,0.9);
  }
  #ua-wsi-crosshair::before { left:0; right:0; top:50%; height:2px; transform:translateY(-60px); }
  #ua-wsi-crosshair::after { top:0; bottom:0; left:50%; width:2px; transform:translateX(-1px); transform:translateY(-60px); }
</style>
<div class="ua-instance">
<div id="ua-root">
  <div id="ua-wrap">
    <canvas id="ua-canvas"></canvas>
    <div id="ua-reset" title="reset zoom/pan">⤾ reset view</div>
    <div id="ua-zoomlabel">100%</div>
    <div id="ua-sizectl" title="dot/pie size">
      <span>size</span>
      <input type="range" id="ua-size-slider" min="2" max="14" step="1" value="3" />
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
        <input id="ua-add-input" autocomplete="off" placeholder="add label…" />
        <div id="ua-add-dropdown"></div>
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
<div id="ua-magnifier">
  <canvas id="ua-magnifier-canvas" width="512" height="512"></canvas>
  <div id="ua-magnifier-status"></div>
</div>
<div id="ua-wsi-crosshair"></div>
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
        f"  const HAS_MAGNIFIER = {json.dumps(have_magnifier)};\n"
        + r"""
  // Scope all element lookups to THIS widget instance's own container, rather
  // than the global document -- Jupyter keeps previous cell outputs in the
  // DOM, so with repeated re-runs the same ids exist multiple times on the
  // page and a bare document.getElementById() would silently grab a stale
  // element from an earlier run instead of the one just rendered.
  const container = document.currentScript.previousElementSibling;
  const getEl = id => container.querySelector('#' + id);
  const canvas  = getEl('ua-canvas');
  const ctx     = canvas.getContext('2d');
  const tooltip = getEl('ua-tooltip');
  const popup       = getEl('ua-popup');
  const popupTitle  = getEl('ua-popup-title');
  const popupImg    = getEl('ua-popup-img');
  const labelList   = getEl('ua-label-list');
  const suggestList  = getEl('ua-suggest-list');
  const popupMeta    = getEl('ua-popup-meta');
  const addInput     = getEl('ua-add-input');
  const addDropdown  = getEl('ua-add-dropdown');
  const autoList     = getEl('ua-auto-list');
  const classifyBtn  = getEl('ua-classify-btn');
  const clearAutoBtn = getEl('ua-clear-auto-btn');
  const autoStatus   = getEl('ua-auto-status');
  const resetBtn     = getEl('ua-reset');
  const zoomLabel   = getEl('ua-zoomlabel');
  const sizeSlider  = getEl('ua-size-slider');
  popupImg.addEventListener('error', () => { popupImg.style.display = 'none'; });
  const wsiBackdrop = getEl('ua-wsi-backdrop');
  const wsiPopup    = getEl('ua-wsi-popup');
  const wsiTitle    = getEl('ua-wsi-title');
  const wsiCanvas   = getEl('ua-wsi-canvas');
  const wsiStatus   = getEl('ua-wsi-status');
  const wsiClose    = getEl('ua-wsi-close');
  const magnifier       = getEl('ua-magnifier');
  const magCanvas       = getEl('ua-magnifier-canvas');
  const magCtx          = magCanvas.getContext('2d');
  const magStatus       = getEl('ua-magnifier-status');
  const wsiCrosshair    = getEl('ua-wsi-crosshair');
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
      else if (addDropdown.style.display === 'block') hideDropdown();
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
          <span class="ua-suggest-chip" data-label="${c.label}" title="click to add this label">
            ${c.label} <span class="ua-suggest-prob">${Math.round(c.prob * 100)}%</span>
          </span>`).join('');
        suggestList.querySelectorAll('.ua-suggest-chip').forEach(chip => {
          chip.addEventListener('click', () => {
            const label = chip.dataset.label;
            postAction({ action: 'add', sample, label });
            chip.remove();   // it's now a real label, shown in the label list instead
          });
        });
      })
      .catch(() => {});
  }

  // ── label autocomplete dropdown ──────────────────────────────────────
  let ddItems = [], ddActive = -1;

  function allKnownLabels() {
    return Object.keys(COLORS).sort((a, b) => a.localeCompare(b));
  }

  function hideDropdown() {
    addDropdown.style.display = 'none';
    addDropdown.innerHTML = '';
    ddItems = []; ddActive = -1;
  }

  function renderDropdown(filter) {
    const q = filter.trim().toLowerCase();
    const all = allKnownLabels();
    ddItems = q ? all.filter(l => l.toLowerCase().includes(q)) : all;
    if (!ddItems.length) { hideDropdown(); return; }
    ddActive = -1;
    addDropdown.innerHTML = ddItems.map((l, i) => `
      <div class="ua-add-dropdown-item" data-idx="${i}">
        <span class="ua-dot" style="background:${COLORS[l] || GRAY}"></span>
        <span class="ua-label-text" title="${l}">${l}</span>
      </div>`).join('');
    addDropdown.querySelectorAll('.ua-add-dropdown-item').forEach(el => {
      el.addEventListener('mousedown', ev => {
        // mousedown (not click) so it fires before the input's blur hides the dropdown
        ev.preventDefault();
        addInput.value = ddItems[parseInt(el.dataset.idx, 10)];
        hideDropdown();
        addInput.focus();
      });
    });
    addDropdown.style.display = 'block';
  }

  function setActiveItem(i) {
    const rows = addDropdown.querySelectorAll('.ua-add-dropdown-item');
    rows.forEach(r => r.classList.remove('active'));
    if (i >= 0 && i < rows.length) {
      rows[i].classList.add('active');
      rows[i].scrollIntoView({ block: 'nearest' });
    }
    ddActive = i;
  }

  addInput.addEventListener('input', () => renderDropdown(addInput.value));
  addInput.addEventListener('focus', () => renderDropdown(addInput.value));
  addInput.addEventListener('blur', () => hideDropdown());

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
    hideDropdown();
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
    popupImg.src = '';               // clear stale thumbnail so it doesn't flash for the next sample
    popupImg.style.display = 'none';
    hideDropdown();
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
    // clear the canvas so a stale slide/contour doesn't flash when the next sample opens
    wsiCanvas.getContext('2d').clearRect(0, 0, wsiCanvas.width, wsiCanvas.height);
    hideMagnifier();
  }

  popupImg.addEventListener('click', openWsiPopup);
  wsiClose.addEventListener('click', closeWsiPopup);
  wsiBackdrop.addEventListener('click', closeWsiPopup);

  // ── magnifying glass: live full-resolution tile under the cursor ────────
  // Fetches are throttled to "one in flight at a time, always chasing the
  // latest cursor position" rather than firing on every mousemove, so fast
  // motion doesn't queue up a burst of stale requests.
  let magSeq = 0, magPending = null, magBusy = false, magLastUrl = null;

  function wsiCanvasXY(e) {
    const r = wsiCanvas.getBoundingClientRect();
    // canvas has a real pixel buffer (wsiCanvas.width/height, set to the
    // whole-slide overview's own dimensions) that may be CSS-scaled to fit
    // the popup (max-width:100%), so convert client coords through both.
    const scaleX = wsiCanvas.width / r.width;
    const scaleY = wsiCanvas.height / r.height;
    return [(e.clientX - r.left) * scaleX, (e.clientY - r.top) * scaleY];
  }

  function positionMagnifier(clientX, clientY) {
    const size = 512;
    const gap = 28;   // keep clear of the sampling point itself
    // prefer up-and-to-the-right of the cursor so the magnifier never covers
    // the region it's currently sampling
    let left = clientX + gap;
    let top  = clientY - size - gap;
    if (top < 4) top = clientY + gap;                 // flip below if no room above
    if (left + size + 4 > window.innerWidth) left = clientX - size - gap;  // flip left if no room right
    left = Math.max(4, Math.min(left, window.innerWidth - size - 4));
    top  = Math.max(4, Math.min(top, window.innerHeight - size - 4));
    magnifier.style.left = left + 'px';
    magnifier.style.top  = top + 'px';
  }

  function runMagnifierFetch() {
    if (!magPending) { magBusy = false; return; }
    const { sample, nx, ny } = magPending;
    magPending = null;
    magBusy = true;
    const mySeq = ++magSeq;
    fetch(SRV + '/tile/' + encodeURIComponent(sample) + `?nx=${nx.toFixed(5)}&ny=${ny.toFixed(5)}`)
      .then(r => { if (!r.ok) throw new Error('tile unavailable'); return r.blob(); })
      .then(blob => {
        if (mySeq !== magSeq) return;
        const url = URL.createObjectURL(blob);
        const img = new Image();
        img.onload = () => {
          if (mySeq !== magSeq) { URL.revokeObjectURL(url); return; }
          magCtx.clearRect(0, 0, 512, 512);
          magCtx.drawImage(img, 0, 0, 512, 512);
          magStatus.textContent = '';
          if (magLastUrl) URL.revokeObjectURL(magLastUrl);
          magLastUrl = url;
        };
        img.src = url;
      })
      .catch(() => { if (mySeq === magSeq) magStatus.textContent = '⚠ tile unavailable'; })
      .finally(runMagnifierFetch);
  }

  function scheduleMagnifierFetch(sample, nx, ny) {
    magPending = { sample, nx, ny };
    if (!magBusy) runMagnifierFetch();
  }

  function hideMagnifier() {
    magnifier.style.display = 'none';
    wsiCrosshair.style.display = 'none';
    // Clear image of the magnifier so it doesn't flash a stale tile when the next sample opens.
    magCtx.clearRect(0, 0, 512, 512);
    magSeq++;   // invalidate any in-flight/queued fetch
    magPending = null;
  }

  wsiCanvas.addEventListener('mousemove', e => {
    if (!HAS_MAGNIFIER || selectedIdx < 0) return;
    const sample = SAMPLES[selectedIdx];
    const [px, py] = wsiCanvasXY(e);
    if (px < 0 || py < 0 || px > wsiCanvas.width || py > wsiCanvas.height) { hideMagnifier(); return; }
    const nx = px / wsiCanvas.width, ny = py / wsiCanvas.height;
    positionMagnifier(e.clientX, e.clientY);
    wsiCrosshair.style.left = (e.clientX - 13) + 'px';
    wsiCrosshair.style.top  = (e.clientY - 13) + 'px';
    wsiCrosshair.style.display = 'block';
    magnifier.style.display = 'block';
    scheduleMagnifierFetch(sample, nx, ny);
  });
  wsiCanvas.addEventListener('mouseleave', hideMagnifier);

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
  getEl('ua-popup-close').addEventListener('click', closePopup);
  getEl('ua-add-btn').addEventListener('click', () => {
    if (selectedIdx < 0) return;
    const label = addInput.value.trim();
    if (!label) return;
    postAction({ action: 'add', sample: SAMPLES[selectedIdx], label });
    addInput.value = '';
    hideDropdown();
  });
  addInput.addEventListener('keydown', e => {
    const rows = addDropdown.style.display === 'block' ? addDropdown.querySelectorAll('.ua-add-dropdown-item') : [];
    if (e.key === 'ArrowDown' && rows.length) {
      e.preventDefault();
      setActiveItem(ddActive < rows.length - 1 ? ddActive + 1 : 0);
    } else if (e.key === 'ArrowUp' && rows.length) {
      e.preventDefault();
      setActiveItem(ddActive > 0 ? ddActive - 1 : rows.length - 1);
    } else if (e.key === 'Enter') {
      if (rows.length && ddActive >= 0) {
        e.preventDefault();
        addInput.value = ddItems[ddActive];
        hideDropdown();
      } else {
        getEl('ua-add-btn').click();
      }
    } else if (e.key === 'Escape') {
      // handled globally, but stop it from also bubbling to close the popup
      // when the dropdown is open (global handler already closes just the dropdown)
      if (rows.length) e.stopPropagation();
    }
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