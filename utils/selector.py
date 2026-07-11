"""
contour_annotator.py — freehand + magic-wand contour annotation for Jupyter notebooks.

Save format (one file per contour):
  <savepath>/<sample>.oid<N>.json
  {
    "0": { "points": [0.104, 0.244, ...] },   <- x, normalised 0-1
    "1": { "points": [0.160, 0.328, ...] }    <- y, normalised 0-1
  }

Magic wand:
  Draw a small rough scribble anywhere inside a tissue region. The scribble is
  sent to the local Python server as a *seed*; the server grows it into the
  full tissue boundary (tissue = any pixel not background-colored, restricted
  to the connected component touching the seed, with automatic crop-window
  growth if the region hits the crop border), fills in any enclosed holes,
  and stops exactly at background. Background is assumed to be a uniform
  value (default pure white, 255) -- see `bg_value`/`bg_tol` on
  `_grow_tissue_region` if yours differs or has slight compression/
  anti-aliasing noise. The refined contour replaces the scribble before it's
  drawn/saved. Toggle the "magic wand" checkbox off to fall back to pure
  free-hand drawing, or hold Cmd (Mac) / Ctrl (Win/Linux) while drawing a
  single stroke to disable it just for that stroke.

Navigation:
  Use "Next ->" / "<- Previous" to move between samples. Revisiting a sample
  reloads any contours already saved to disk for it, so you can review or
  keep editing earlier work.

Clearing contours:
  "Clear last" removes the most recently drawn contour *and* deletes its
  saved JSON file from disk, if one was written for it. "Clear all" removes
  every contour for the *current sample only* and deletes every JSON file on
  disk that belongs to that sample (it never touches other samples' files).

Images are loaded lazily: each sample's image is only encoded and sent to the
browser the first time it's actually displayed (via a small HTTP GET to the
local server), instead of encoding every sample up front and embedding all of
them as base64 in the notebook cell's output. This makes the initial cell
render effectively O(1) instead of O(number of samples).
"""

import base64, glob, json, os, re, socket, threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from io import BytesIO
from urllib.parse import urlparse

import cv2
import numpy as np
from PIL import Image
from IPython.display import display, HTML
import tifffile
from tqdm import tqdm

# ── internal tuning constants (not exposed in the UI) ───────────────────────
# Dilate the final grown tissue contour outward by this many pixels
# (full-resolution image space) before it's returned. Set to 0 to disable.
CONTOUR_GROW_PX = 10
# Thickness (full-res px) of separator-line barriers burned into the tissue mask
# before connected-component analysis, so the wand stops at drawn separators.
SEP_BARRIER_PX = 5

def loadAndFilter(paths, spath=None, downsample=1, maskfunc=None):
    if maskfunc is None:
        print('Provide a mask function, such as wsiMask.getInTissuePixelMask from STQ')
        raise
    imgs = []
    samples = []
    for i, p in tqdm(enumerate(paths), total=len(paths)):
        sample = p.split('/')[-1]
        sname = f"{spath}/{sample}.tiff"
        if not os.path.isfile(sname):
            img = tifffile.imread(p, level=2)
            whp = (img>250).all(axis=2)
            img[whp] = np.quantile(img[~whp], 0.99)
            m = maskfunc(img, None, kernel_size=7)
            wh = m.T==0
            img[wh, ...] = 255
            if not spath is None:
                tifffile.imwrite(sname, img, compression='zlib')
        else:
            img = tifffile.imread(sname)
        imgs.append(img[::downsample, ::downsample])
        samples.append(sample)
    return samples, imgs

def setNotebookWidth(widthPercent=100):
    """Set the notebook container width in a Jupyter environment."""
    display(HTML(f"""<style>:root {{ --jp-notebook-max-width: {widthPercent}%; }}
body .container, div.container {{ width: {widthPercent}% !important; }}</style>"""))
    return

def _fill_holes(mask_u8: np.ndarray) -> np.ndarray:
    """
    Fill enclosed background pockets inside a binary foreground mask.
    mask_u8: single-channel 0/255 uint8 array.
    Classic trick: flood-fill the *background* starting from a corner, then
    anything NOT reached by that flood-fill was an enclosed hole. Pads by one
    pixel first so the corner is guaranteed to be background even if the
    foreground mask itself touches the array edge.
    """
    padded = np.pad(mask_u8, 1, mode='constant', constant_values=0)
    h, w = padded.shape
    flood = padded.copy()
    ff_mask = np.zeros((h + 2, w + 2), np.uint8)
    cv2.floodFill(flood, ff_mask, (0, 0), 255)
    holes_filled = padded | cv2.bitwise_not(flood)
    return holes_filled[1:-1, 1:-1]


# ────────────────────────────────────────────────────────────────────────────
# Region-growing ("magic wand") core
# ────────────────────────────────────────────────────────────────────────────

def _grow_tissue_region(
    img_rgb: np.ndarray,
    seed_pts_px: np.ndarray,
    separator_lines: list | None = None,
    bg_value: int = 255,
    bg_tol: int = 0,
    margin_factor: float = 4.0,
    min_margin: float = 40.0,
    max_crop: float = 6000.0,
    max_expansions: int = 3,
    close_iters: int = 1,
    open_iters: int = 1,
):
    """
    Grow a rough seed scribble into a full tissue-region contour, stopping
    exactly at background pixels (all channels >= bg_value - bg_tol).

    img_rgb       : full-resolution HxWx3 uint8 image the seed was drawn on
    seed_pts_px   : Nx2 array of (x, y) pixel coords of the scribble, in img_rgb space
    bg_value      : background pixel value (default 255, i.e. pure white)
    bg_tol        : tolerance below bg_value still counted as background
                    (raise this a little, e.g. 2-5, if edges are anti-aliased/compressed)

    Returns an Mx2 float array of (x, y) pixel coords (in img_rgb space)
    describing the refined boundary, or None if no region could be found.
    """
    seed_pts_px = np.asarray(seed_pts_px, dtype=np.float64)
    if len(seed_pts_px) == 0:
        return None

    cx, cy = seed_pts_px.mean(axis=0)
    x0s, y0s = seed_pts_px.min(axis=0)
    x1s, y1s = seed_pts_px.max(axis=0)
    seed_w, seed_h = max(x1s - x0s, 1.0), max(y1s - y0s, 1.0)

    H, W = img_rgb.shape[:2]
    component = None
    x0 = y0 = 0

    for attempt in range(max_expansions + 1):
        grow = 1.6 ** attempt
        half_w = min(max(seed_w * margin_factor, min_margin) * grow, max_crop / 2)
        half_h = min(max(seed_h * margin_factor, min_margin) * grow, max_crop / 2)

        x0 = int(max(cx - half_w, 0)); x1 = int(min(cx + half_w, W))
        y0 = int(max(cy - half_h, 0)); y1 = int(min(cy + half_h, H))
        crop = img_rgb[y0:y1, x0:x1]
        if crop.size == 0:
            return None

        # tissue = any channel not (near-)background; background = all
        # channels >= bg_value - bg_tol (e.g. pure white == 255)
        is_bg = np.all(crop >= (bg_value - bg_tol), axis=-1)
        mask = (~is_bg).astype(np.uint8) * 255

        if close_iters or open_iters:
            kernel = np.ones((5, 5), np.uint8)
            if close_iters:  # fill tiny background specks inside tissue
                mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE, kernel, iterations=close_iters)
            if open_iters:   # drop tiny tissue specks/dust in background
                mask = cv2.morphologyEx(mask, cv2.MORPH_OPEN, kernel, iterations=open_iters)

        seed_x = int(np.clip(cx - x0, 0, mask.shape[1] - 1))
        seed_y = int(np.clip(cy - y0, 0, mask.shape[0] - 1))

        # Burn separator lines into the mask as background barriers so the
        # connected-component flood cannot cross them.
        if separator_lines:
            for _sep in separator_lines:
                if len(_sep) < 2:
                    continue
                _pts = np.round(
                    np.asarray(_sep, dtype=np.float64) - np.array([[x0, y0]])
                ).astype(np.int32)
                for _j in range(len(_pts) - 1):
                    _p1 = (int(np.clip(_pts[_j,     0], 0, mask.shape[1] - 1)),
                           int(np.clip(_pts[_j,     1], 0, mask.shape[0] - 1)))
                    _p2 = (int(np.clip(_pts[_j + 1, 0], 0, mask.shape[1] - 1)),
                           int(np.clip(_pts[_j + 1, 1], 0, mask.shape[0] - 1)))
                    cv2.line(mask, _p1, _p2, 0, thickness=SEP_BARRIER_PX)

        n_labels, labels = cv2.connectedComponents(mask)
        seed_label = labels[seed_y, seed_x]

        if seed_label == 0:
            # seed itself landed on a background-classified pixel (e.g. right
            # on the boundary) -- nudge to the nearest tissue pixel nearby.
            ys, xs = np.where(mask > 0)
            if len(xs) == 0:
                return None
            d2 = (xs - seed_x) ** 2 + (ys - seed_y) ** 2
            nearest = np.argmin(d2)
            seed_label = labels[ys[nearest], xs[nearest]]

        # Collect ALL component labels touched by any point in the scribble
        # path, so tissues that the stroke crosses are all included.
        touched_labels = set()
        for _sp in seed_pts_px:
            _sx = int(np.clip(_sp[0] - x0, 0, mask.shape[1] - 1))
            _sy = int(np.clip(_sp[1] - y0, 0, mask.shape[0] - 1))
            _lbl = labels[_sy, _sx]
            if _lbl != 0:
                touched_labels.add(int(_lbl))
        # Always include the centroid-seed label as fallback
        if seed_label != 0:
            touched_labels.add(int(seed_label))
        if not touched_labels:
            return None

        component = np.isin(labels, list(touched_labels))

        touches_border = (
            component[0, :].any() or component[-1, :].any()
            or component[:, 0].any() or component[:, -1].any()
        )
        if touches_border and attempt < max_expansions:
            continue  # crop was too small to contain the whole tissue piece
        break

    if component is None or not component.any():
        return None

    comp_u8 = (component.astype(np.uint8)) * 255
    comp_u8 = _fill_holes(comp_u8)  # no holes allowed inside a tissue contour

    if CONTOUR_GROW_PX > 0:
        # pad so the dilation has room to expand without being clipped by
        # the crop's own array edges
        pad = CONTOUR_GROW_PX + 5
        comp_u8 = np.pad(comp_u8, pad, mode='constant', constant_values=0)
        x0 -= pad
        y0 -= pad
        kernel = cv2.getStructuringElement(
            cv2.MORPH_ELLIPSE, (2 * CONTOUR_GROW_PX + 1, 2 * CONTOUR_GROW_PX + 1)
        )
        comp_u8 = cv2.dilate(comp_u8, kernel)

    contours, _ = cv2.findContours(comp_u8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if not contours:
        return None
    # When multiple tissues were selected, take the convex-hull-merged outer
    # boundary: find the single largest contour by area, but if there are
    # several, also draw the connecting stroke pixels onto comp_u8 first so
    # the regions are bridged, then re-extract.
    if len(contours) > 1:
        # Draw the original seed path into comp_u8 as a thick bridge
        bridge_pts = np.round(
            seed_pts_px - np.array([[x0, y0]])
        ).astype(np.int32)
        for _j in range(len(bridge_pts) - 1):
            _p1 = (int(np.clip(bridge_pts[_j,     0], 0, comp_u8.shape[1] - 1)),
                   int(np.clip(bridge_pts[_j,     1], 0, comp_u8.shape[0] - 1)))
            _p2 = (int(np.clip(bridge_pts[_j + 1, 0], 0, comp_u8.shape[1] - 1)),
                   int(np.clip(bridge_pts[_j + 1, 1], 0, comp_u8.shape[0] - 1)))
            cv2.line(comp_u8, _p1, _p2, 255, thickness=max(SEP_BARRIER_PX, 4))
        comp_u8 = _fill_holes(comp_u8)
        contours, _ = cv2.findContours(comp_u8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        if not contours:
            return None
    contour = max(contours, key=cv2.contourArea)
    if cv2.contourArea(contour) < 25:
        return None

    peri = cv2.arcLength(contour, True)
    # Use a small absolute epsilon so all tissues get equally dense point
    # sampling — large tissues were coarser when epsilon scaled with perimeter.
    approx = cv2.approxPolyDP(contour, epsilon=2.0, closed=True)
    pts_full = approx.reshape(-1, 2).astype(np.float64)
    pts_full[:, 0] += x0
    pts_full[:, 1] += y0
    pts_full[:, 0] = np.clip(pts_full[:, 0], 0, W - 1)
    pts_full[:, 1] = np.clip(pts_full[:, 1], 0, H - 1)
    return pts_full


def annotateContours(
    samples: list[str],
    images: list[np.ndarray],
    canvas_max_dim: int = 800,
    savepath: str = "./contours/",
    mpp: float | None = None,
):
    assert len(samples) == len(images), "samples and images must have the same length"
    savepath = savepath.rstrip('/')
    os.makedirs(savepath, exist_ok=True)

    sizes = [(int(img.shape[1]), int(img.shape[0])) for img in images]

    # Encoded-image cache: images are only PNG-encoded the first time they're
    # actually requested by the browser (lazy), not all up front. This is
    # what keeps initial cell display fast regardless of how many samples
    # there are.
    _img_cache: dict[int, bytes] = {}

    def _encode_image(idx: int) -> bytes:
        if idx not in _img_cache:
            arr = images[idx]
            if arr.dtype != np.uint8:
                arr = np.clip(arr, 0, 255).astype(np.uint8)
            buf = BytesIO()
            Image.fromarray(arr).save(buf, format="PNG")
            _img_cache[idx] = buf.getvalue()
        return _img_cache[idx]

    _OID_RE = re.compile(r'\.oid(\d+)\.json$')

    def _sample_files(sample: str) -> list[str]:
        return glob.glob(os.path.join(savepath, sample + '.oid*.json'))

    def _oid_of(fname: str) -> int:
        m = _OID_RE.search(fname)
        return int(m.group(1)) if m else -1

    # ── tiny HTTP server ─────────────────────────────────────────────────────
    # JS talks to this over GET (images) and POST (save/delete/list/wand).
    # Runs in a daemon thread — dies automatically when the notebook kernel stops.

    with socket.socket() as s:
        s.bind(('', 0))
        port = s.getsockname()[1]

    class Handler(BaseHTTPRequestHandler):
        def do_OPTIONS(self):          # preflight CORS
            self._ok(b'')

        def do_GET(self):
            parsed = urlparse(self.path)
            if parsed.path.startswith('/image/'):
                try:
                    idx = int(parsed.path.rsplit('/', 1)[-1])
                    body = _encode_image(idx)
                except Exception:
                    self.send_response(404)
                    self.end_headers()
                    return
                self.send_response(200)
                self.send_header('Content-Type', 'image/png')
                self.send_header('Content-Length', str(len(body)))
                self.send_header('Cache-Control', 'public, max-age=604800, immutable')
                self.send_header('Access-Control-Allow-Origin', '*')
                self.end_headers()
                self.wfile.write(body)
            else:
                self.send_response(404)
                self.end_headers()

        def do_POST(self):
            n    = int(self.headers.get('Content-Length', 0))
            body = self.rfile.read(n)
            try:
                payload = json.loads(body)
                action  = payload.get('action', 'save')
                handler = {
                    'save':       self._handle_save,
                    'wand':       self._handle_wand,
                    'delete':     self._handle_delete,
                    'delete_all': self._handle_delete_all,
                    'list':       self._handle_list,
                    'save_sep':   self._handle_save_sep,
                    'list_sep':   self._handle_list_sep,
                }.get(action, self._handle_save)
                handler(payload)
            except Exception as e:
                self._ok(json.dumps({'ok': False, 'error': str(e)}).encode())

        def _handle_save(self, payload):
            sample = payload['sample']
            data   = payload['data']
            used = {_oid_of(f) for f in _sample_files(sample)}
            oid = 0
            while oid in used:
                oid += 1
            fname = os.path.join(savepath, f'{sample}.oid{oid}.json')
            with open(fname, 'w') as fh:
                json.dump(data, fh, indent=2)
            self._ok(json.dumps({'ok': True, 'file': fname}).encode())

        def _handle_delete(self, payload):
            fname = payload.get('file')
            # safety check: only ever delete files that live directly inside
            # savepath and look like our own <sample>.oid<N>.json files.
            if fname and os.path.dirname(os.path.abspath(fname)) == os.path.abspath(savepath) \
                    and _OID_RE.search(fname):
                try:
                    os.remove(fname)
                except FileNotFoundError:
                    pass
            self._ok(b'{"ok":true}')

        def _handle_delete_all(self, payload):
            sample = payload['sample']
            for f in _sample_files(sample):
                try:
                    os.remove(f)
                except FileNotFoundError:
                    pass
            self._ok(b'{"ok":true}')

        def _handle_save_sep(self, payload):
            sample = payload['sample']
            lines  = payload.get('lines', [])
            fname  = os.path.join(savepath, f'{sample}.sep.json')
            if lines:
                with open(fname, 'w') as fh:
                    json.dump({'lines': lines}, fh, indent=2)
            else:
                try:
                    os.remove(fname)
                except FileNotFoundError:
                    pass
            self._ok(b'{"ok":true}')

        def _handle_list_sep(self, payload):
            sample = payload['sample']
            fname  = os.path.join(savepath, f'{sample}.sep.json')
            try:
                with open(fname) as fh:
                    d = json.load(fh)
                self._ok(json.dumps({'ok': True, 'lines': d.get('lines', [])}).encode())
            except FileNotFoundError:
                self._ok(b'{"ok":true,"lines":[]}')
            except Exception as e:
                self._ok(json.dumps({'ok': False, 'error': str(e)}).encode())

        def _handle_list(self, payload):
            sample = payload['sample']
            files = sorted(_sample_files(sample), key=_oid_of)
            out = []
            for f in files:
                try:
                    with open(f) as fh:
                        d = json.load(fh)
                    xs = d["0"]["points"]
                    ys = d["1"]["points"]
                    out.append({'file': f, 'points': [[x, y] for x, y in zip(xs, ys)]})
                except Exception:
                    continue
            self._ok(json.dumps({'ok': True, 'contours': out}).encode())

        def _handle_wand(self, payload):
            idx              = int(payload['idx'])
            canvas_w         = float(payload['canvas_w'])
            canvas_h         = float(payload['canvas_h'])
            pts              = payload['points']
            sep_lines_canvas = payload.get('sep_lines', [])

            full_w, full_h = sizes[idx]
            sx, sy = full_w / canvas_w, full_h / canvas_h
            seed_full = np.array([[x * sx, y * sy] for x, y in pts])
            sep_full  = [
                np.array([[p[0] * sx, p[1] * sy] for p in line])
                for line in sep_lines_canvas
            ]

            result = _grow_tissue_region(
                images[idx], seed_full,
                separator_lines=sep_full if sep_full else None,
            )
            if result is None:
                self._ok(json.dumps({'ok': False, 'error': 'no region found'}).encode())
                return

            canvas_pts = [[x / sx, y / sy] for x, y in result]
            self._ok(json.dumps({'ok': True, 'points': canvas_pts}).encode())

        def _ok(self, body):
            self.send_response(200)
            self.send_header('Content-Type',                 'application/json')
            self.send_header('Content-Length',               str(len(body)))
            self.send_header('Access-Control-Allow-Origin',  '*')
            self.send_header('Access-Control-Allow-Methods', 'GET, POST, OPTIONS')
            self.send_header('Access-Control-Allow-Headers', 'Content-Type')
            self.end_headers()
            if body:
                self.wfile.write(body)

        def log_message(self, *a):
            pass

    host = socket.gethostbyname(socket.gethostname())
    srv = ThreadingHTTPServer((host, port), Handler)
    t   = threading.Thread(target=srv.serve_forever, daemon=True)
    t.start()

    # ── HTML + JS ────────────────────────────────────────────────────────────
    css = """
<style>
  :root {
    --bg:#1a1a2e; --accent:#e94560; --accent2:#0f3460;
    --text:#eaeaea; --muted:#8892a4; --radius:6px;
    --font:'JetBrains Mono','Fira Code','Consolas',monospace;
  }
  #ca-root * { box-sizing:border-box; margin:0; padding:0; }
  #ca-root {
    font-family:var(--font); background:var(--bg); color:var(--text);
    padding:16px; border-radius:10px; display:inline-flex;
    flex-direction:column; gap:12px; user-select:none;
  }
  #ca-header { display:flex; align-items:center; gap:12px; }
  #ca-title { font-size:11px; letter-spacing:.12em; text-transform:uppercase; color:var(--accent); }
  #ca-sample-label {
    font-size:13px; color:var(--text); background:var(--accent2);
    padding:3px 10px; border-radius:var(--radius);
  }
  #ca-counter { margin-left:auto; font-size:11px; color:var(--muted); }
  #ca-canvas-wrap {
    position:relative; border:1.5px solid var(--accent2);
    border-radius:var(--radius); overflow:hidden; cursor:crosshair;
  }
  #ca-img-canvas, #ca-draw-canvas { display:block; }
  #ca-draw-canvas { position:absolute; top:0; left:0; outline:none; }
  #ca-toolbar { display:flex; gap:8px; align-items:center; flex-wrap:wrap; }
  .ca-btn {
    font-family:var(--font); font-size:12px; letter-spacing:.06em;
    padding:6px 14px; border:none; border-radius:var(--radius);
    cursor:pointer; transition:filter .15s;
  }
  .ca-btn:hover  { filter:brightness(1.15); }
  .ca-btn:active { filter:brightness(.9); }
  .ca-btn:disabled { opacity:.35; cursor:default; filter:none; }
  #btn-prev, #btn-next { background:var(--accent2); color:var(--text); }
  #btn-save       { background:var(--accent);  color:#fff; }
  #btn-clear-last { background:#2a2a4a; color:var(--text); }
  #btn-clear-all  { background:#2a2a4a; color:var(--muted); }
  #ca-autosave-wrap, #ca-wand-wrap, #ca-ruler-wrap {
    display:flex; align-items:center; gap:6px;
    font-size:11px; color:var(--muted); margin-left:4px;
  }
  #ca-autosave-wrap input, #ca-wand-wrap input, #ca-ruler-wrap input { accent-color:var(--accent); width:14px; height:14px; cursor:pointer; }
  #ca-ruler-wrap.disabled { opacity:.4; }
  #ca-ruler-wrap.disabled input { cursor:not-allowed; pointer-events:none; }
  #ca-status { font-size:11px; color:var(--muted); min-height:16px; }
  #ca-contour-list { font-size:10px; color:var(--muted); max-height:52px; overflow-y:auto; line-height:1.6; }
</style>
<div id="ca-root">
  <div id="ca-header">
    <span id="ca-title">Contour Annotator</span>
    <span id="ca-sample-label">—</span>
    <span id="ca-counter">—</span>
  </div>
  <div id="ca-canvas-wrap">
    <canvas id="ca-img-canvas"></canvas>
    <canvas id="ca-draw-canvas" tabindex="0"></canvas>
  </div>
  <div id="ca-toolbar">
    <button class="ca-btn" id="btn-prev">← Previous</button>
    <button class="ca-btn" id="btn-save">Save contour</button>
    <button class="ca-btn" id="btn-clear-last">Clear last</button>
    <button class="ca-btn" id="btn-clear-all">Clear all</button>
    <button class="ca-btn" id="btn-next">Next →</button>
    <label id="ca-wand-wrap">
      <input type="checkbox" id="ca-wand" checked> magic wand
    </label>
    <label id="ca-autosave-wrap">
      <input type="checkbox" id="ca-autosave" checked> autosave
    </label>
    <label id="ca-ruler-wrap">
      <input type="checkbox" id="ca-ruler"> ruler
    </label>
  </div>
  <div id="ca-status">Loading…</div>
  <div id="ca-contour-list"></div>
</div>
"""

    js = (
        "<script>\n(function(){\n"
        f"  const SIZES    = {json.dumps(sizes)};\n"
        f"  const SAMPLES  = {json.dumps(samples)};\n"
        f"  const MAX_DIM  = {canvas_max_dim};\n"
        f"  const SRV_PORT = {port};\n"
        f"  const SRV_HOST = {json.dumps(host)};\n"
        f"  const MPP      = {json.dumps(mpp)};\n"
        + r"""
  let imgIdx      = 0;
  let contours    = [];   // [{points: [[x,y],...], file: string|null}]
  let separators  = [];   // [{points: [[x,y],...]}]  — open barrier lines
  let selectedSepIdx = -1;
  let drawMode    = 'contour';  // 'contour' | 'separator'
  let currentPath = [];
  let drawing     = false;
  let strokeWandDisabled = false;   // Cmd/Ctrl held at mousedown -> skip wand for this stroke
  const prefetched = new Set();     // sample indices whose image bytes we've already requested

  const imgCanvas   = document.getElementById('ca-img-canvas');
  const drawCanvas  = document.getElementById('ca-draw-canvas');
  const ic = imgCanvas.getContext('2d', {willReadFrequently: true});
  const dc = drawCanvas.getContext('2d');
  const sampleLabel = document.getElementById('ca-sample-label');
  const counterEl   = document.getElementById('ca-counter');
  const statusEl    = document.getElementById('ca-status');
  const contourList = document.getElementById('ca-contour-list');
  const autosaveCb  = document.getElementById('ca-autosave');
  const wandCb      = document.getElementById('ca-wand');
  const rulerCb     = document.getElementById('ca-ruler');
  const rulerWrap   = document.getElementById('ca-ruler-wrap');
  const prevBtn     = document.getElementById('btn-prev');
  const nextBtn     = document.getElementById('btn-next');

  // Ruler: enabled only when MPP was provided by Python.
  if (MPP === null) {
    rulerCb.disabled = true;
    rulerWrap.classList.add('disabled');
    rulerWrap.title = 'Provide mpp to enable ruler';
  } else {
    rulerCb.checked = true;
  }

  const srvUrl = () => 'http://' + SRV_HOST + ':' + SRV_PORT;

  // ── save / delete / list via local HTTP server ──────────────────────
  async function doSave(idx) {
    const sample = SAMPLES[imgIdx];
    const pts    = contours[idx].points;
    const [cw, ch] = [drawCanvas.width, drawCanvas.height];
    const data  = {
      "0": { "points": pts.map(p => parseFloat((p[0] / cw).toFixed(4))) },
      "1": { "points": pts.map(p => parseFloat((p[1] / ch).toFixed(4))) }
    };
    statusEl.textContent = '💾 Saving contour #' + idx + '…';
    try {
      const resp = await fetch(srvUrl(), {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({action: 'save', sample: sample, data: data})
      });
      const json = await resp.json();
      if (json.ok) {
        contours[idx].file = json.file;
        statusEl.textContent = '✅ Saved → ' + json.file + '  (' + pts.length + ' pts)';
      } else {
        statusEl.textContent = '⚠ Server error: ' + json.error;
      }
    } catch(e) {
      statusEl.textContent = '⚠ Fetch failed: ' + e;
    }
  }

  async function deleteFile(fname) {
    try {
      await fetch(srvUrl(), {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({action: 'delete', file: fname})
      });
    } catch (e) { /* best effort */ }
  }

  async function deleteAllForSample(sample) {
    try {
      await fetch(srvUrl(), {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({action: 'delete_all', sample: sample})
      });
    } catch (e) { /* best effort */ }
  }

  async function listSavedContours(sample) {
    try {
      const resp = await fetch(srvUrl(), {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({action: 'list', sample: sample})
      });
      const json = await resp.json();
      return json.ok ? json.contours : [];
    } catch (e) {
      return [];
    }
  }

  // ── separator lines: save / load ─────────────────────────────────
  async function listSavedSeparators(sample) {
    try {
      const resp = await fetch(srvUrl(), {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({action: 'list_sep', sample: sample})
      });
      const json = await resp.json();
      return json.ok ? json.lines : [];
    } catch (e) {
      return [];
    }
  }

  async function doSaveSeparators() {
    const sample   = SAMPLES[imgIdx];
    const [cw, ch] = [drawCanvas.width, drawCanvas.height];
    const lines = separators.map(s =>
      s.points.map(([x, y]) => [parseFloat((x / cw).toFixed(4)), parseFloat((y / ch).toFixed(4))])
    );
    try {
      await fetch(srvUrl(), {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({action: 'save_sep', sample: sample, lines: lines})
      });
      statusEl.textContent = '✅ Separator saved (' + separators.length + ' line' + (separators.length !== 1 ? 's' : '') + ' total).';
    } catch(e) {
      statusEl.textContent = '⚠ Could not save separator: ' + e;
    }
  }

  function commitSeparator(pts) {
    separators.push({points: pts});
    redraw(); updateList();
    doSaveSeparators();
  }

  // ── magic wand: send scribble as seed, get refined contour back ────
  async function growRegion(pts) {
    const resp = await fetch(srvUrl(), {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({
        action: 'wand',
        idx: imgIdx,
        canvas_w: drawCanvas.width,
        canvas_h: drawCanvas.height,
        points: pts,
        sep_lines: separators.map(s => s.points)
      })
    });
    return await resp.json();
  }

  // ── prefetch: warm the server's encode cache + the browser's HTTP cache
  // for a neighboring sample while the user is still looking at this one.
  function prefetchImage(idx) {
    if (idx < 0 || idx >= SIZES.length || prefetched.has(idx)) return;
    prefetched.add(idx);
    const im = new Image();
    im.crossOrigin = 'anonymous';
    im.src = srvUrl() + '/image/' + idx;
  }

  // ── geometry ──────────────────────────────────────────────────────
  function scaledSize(w, h) {
    const r = Math.min(MAX_DIM / w, MAX_DIM / h, 1);
    return [Math.round(w * r), Math.round(h * r)];
  }

  async function loadImage(idx) {
    if (!SIZES.length || idx < 0 || idx >= SIZES.length || !SIZES[idx]) {
      statusEl.textContent = '⚠ No image data — re-run the cell to reload.';
      return;
    }
    const [w, h]   = SIZES[idx];
    const [cw, ch] = scaledSize(w, h);
    imgCanvas.width = drawCanvas.width = cw;
    imgCanvas.height = drawCanvas.height = ch;
    const img = new Image();
    img.crossOrigin = 'anonymous';
    img.onload = () => ic.drawImage(img, 0, 0, cw, ch);
    // Lazy: only fetched from the server the moment it's actually shown.
    img.src = srvUrl() + '/image/' + idx;
    prefetchImage(idx + 1);
    prefetchImage(idx - 1);
    sampleLabel.textContent = SAMPLES[idx];
    counterEl.textContent   = (idx + 1) + ' / ' + SIZES.length;
    prevBtn.disabled = idx === 0;
    currentPath = [];
    contours = []; separators = []; selectedSepIdx = -1;
    redraw(); updateList();
    statusEl.textContent = 'Loading saved contours…';

    // Restore contours and separator lines already saved to disk.
    const results = await Promise.all([
      listSavedContours(SAMPLES[idx]),
      listSavedSeparators(SAMPLES[idx])
    ]);
    const saved    = results[0] || [];
    const savedSep = results[1] || [];
    if (idx === imgIdx) {   // guard against a fast double-click racing us
      contours = saved.map(c => ({
        points: c.points.map(([nx, ny]) => [nx * cw, ny * ch]),
        file: c.file
      }));
      separators = savedSep.map(line => ({
        points: line.map(([nx, ny]) => [nx * cw, ny * ch])
      }));
      redraw(); updateList();
    }

    statusEl.textContent = wandCb.checked
      ? 'Scribble inside a tissue — it will snap to the tissue boundary. (Hold ⌘/Ctrl to draw free-hand.)'
      : (autosaveCb.checked
          ? 'Draw a contour — saved automatically on mouse-up.'
          : 'Draw a contour, then click Save contour.');
  }

  function canvasXY(e) {
    const r = drawCanvas.getBoundingClientRect();
    const x = Math.min(Math.max(e.clientX - r.left, 0), drawCanvas.width);
    const y = Math.min(Math.max(e.clientY - r.top, 0), drawCanvas.height);
    return [x, y];
  }

  // ── drawing ───────────────────────────────────────────────────────
  function redraw() {
    dc.clearRect(0, 0, drawCanvas.width, drawCanvas.height);
    // Draw separator barrier lines first (below contours)
    separators.forEach((sep, i) => {
      const pts = sep.points;
      if (pts.length < 2) return;
      const sel = i === selectedSepIdx;
      dc.save();
      dc.beginPath();
      dc.moveTo(pts[0][0], pts[0][1]);
      pts.slice(1).forEach(p => dc.lineTo(p[0], p[1]));
      if (sel) { dc.shadowColor = '#ffff00'; dc.shadowBlur = 16; }
      dc.strokeStyle = sel ? '#ffff00' : 'rgba(255,210,0,0.9)';
      dc.lineWidth   = sel ? 3 : 2;
      dc.setLineDash([8, 4]);
      dc.stroke();
      dc.restore();
      // endpoint dots for clarity
      dc.save();
      dc.fillStyle = sel ? '#ffff00' : 'rgba(255,210,0,0.9)';
      [pts[0], pts[pts.length - 1]].forEach(([px, py]) => {
        dc.beginPath(); dc.arc(px, py, sel ? 4 : 3, 0, 2 * Math.PI); dc.fill();
      });
      dc.restore();
    });
    contours.forEach((c, i) => {
      const pts = c.points;
      if (pts.length < 2) return;
      [['rgba(255,255,255,0.65)', 4], ['#000', 2]].forEach(([col, lw]) => {
        dc.beginPath(); dc.moveTo(pts[0][0], pts[0][1]);
        pts.slice(1).forEach(p => dc.lineTo(p[0], p[1]));
        dc.closePath(); dc.strokeStyle = col; dc.lineWidth = lw;
        dc.setLineDash([]); dc.stroke();
      });
      dc.font = 'bold 16px monospace'; dc.setLineDash([]);
      dc.lineWidth = 3; dc.strokeStyle = 'white';
      dc.strokeText('#' + i, pts[0][0] + 4, pts[0][1] - 4);
      dc.fillStyle = '#000'; dc.fillText('oid' + i, pts[0][0] + 4, pts[0][1] - 4);
    });
    if (currentPath.length > 1) {
      dc.beginPath(); dc.moveTo(currentPath[0][0], currentPath[0][1]);
      currentPath.slice(1).forEach(p => dc.lineTo(p[0], p[1]));
      dc.strokeStyle = '#000'; dc.lineWidth = 1.5;
      dc.setLineDash([5, 4]); dc.stroke(); dc.setLineDash([]);
    }
    if (MPP !== null && rulerCb.checked) _drawRuler();
  }

  // ── ruler ─────────────────────────────────────────────────────────
  function _drawRuler() {
    const cw = drawCanvas.width, ch = drawCanvas.height;
    const [fw, fh] = (SIZES[imgIdx] || [cw, ch]);
    // canvas pixels per 1 mm  (MPP is µm/px at full resolution)
    const pxPerMm_x = (1000 / MPP) * (cw / fw);
    const pxPerMm_y = (1000 / MPP) * (ch / fh);
    if (pxPerMm_x < 1.5 || pxPerMm_y < 1.5) return;  // too zoomed out

    const MINOR = 5;   // mm tick height px
    const MAJOR = 11;  // cm tick height px
    const BAND  = 5;  // background strip thickness
    const COL   = 'rgba(0,0,255,0.93)';
    const SHAD  = 'rgba(0,0,0,0.6)';

    dc.save();
    dc.setLineDash([]);

    // Draw one ruler strip along one edge.
    // axis 'x': ticks at regular x positions, edge at y=edgeY, inward=±1
    // axis 'y': ticks at regular y positions, edge at x=edgeX, inward=±1
    function _edge(axis, edgePx, inward) {
      const pxPerMm = axis === 'x' ? pxPerMm_x : pxPerMm_y;
      const totalPx = axis === 'x' ? cw : ch;
      const nMm     = Math.ceil(totalPx / pxPerMm);

      // background strip
      dc.fillStyle = 'rgba(15,15,30,0.55)';
      if (axis === 'x') {
        const y0 = inward > 0 ? edgePx : edgePx - BAND;
        dc.fillRect(0, y0, cw, BAND);
      } else {
        const x0 = inward > 0 ? edgePx : edgePx - BAND;
        dc.fillRect(x0, 0, BAND, ch);
      }

      for (let mm = 0; mm <= nMm; mm++) {
        const isCm    = mm % 10 === 0;
        const tickLen = isCm ? MAJOR : MINOR;
        const pos     = mm * pxPerMm;
        if (pos > totalPx + 1) break;

        // tick line
        dc.beginPath();
        if (axis === 'x') {
          dc.moveTo(pos, edgePx);
          dc.lineTo(pos, edgePx + inward * tickLen);
        } else {
          dc.moveTo(edgePx, pos);
          dc.lineTo(edgePx + inward * tickLen, pos);
        }
        dc.strokeStyle = SHAD; dc.lineWidth = 2; dc.stroke();
        dc.strokeStyle = COL;  dc.lineWidth = 1; dc.stroke();

        // cm label
        if (isCm && mm > 0) {
          const label = (mm / 10) + ' cm';
          dc.font = 'normal 14px halfspace';
          if (axis === 'x') {
            dc.textAlign = 'center';
            dc.textBaseline = inward > 0 ? 'top' : 'bottom';
            const lx = pos, ly = edgePx + inward * (MAJOR + 1);
            dc.fillStyle = SHAD; dc.fillText(label, lx + 0.5, ly + 0.5);
            dc.fillStyle = COL;  dc.fillText(label, lx, ly);
          } else {
            dc.textAlign = inward > 0 ? 'left' : 'right';
            dc.textBaseline = 'middle';
            const lx = edgePx + inward * (MAJOR + 1), ly = pos;
            dc.fillStyle = SHAD; dc.fillText(label, lx + 0.5, ly + 0.5);
            dc.fillStyle = COL;  dc.fillText(label, lx, ly);
          }
        }
      }
    }

    _edge('x', 0,  +1);   // top
    //_edge('x', ch, -1);   // bottom
    _edge('y', 0,  +1);   // left
    //_edge('y', cw, -1);   // right
    dc.restore();
  }

  function updateList() {
    const cParts = contours.map((c, i)   => '#' + i + ': ' + c.points.length + ' pts' + (c.file ? '' : ' (unsaved)'));
    const sParts = separators.map((s, i) => 'sep#' + i + ': ' + s.points.length + ' pts');
    contourList.textContent = [...cParts, ...sParts].join('  ·  ');
  }

  function commitContour(pts) {
    contours.push({points: pts, file: null});
    const idx = contours.length - 1;
    redraw(); updateList();
    if (autosaveCb.checked) doSave(idx);
    else statusEl.textContent = 'Contour drawn. Click Save contour to write file.';
  }

  // ── mouse ─────────────────────────────────────────────────────────
  drawCanvas.addEventListener('mousedown', e => {
    const [x, y] = canvasXY(e);
    // Detect background: if the pixel under the cursor is near-white, start a
    // separator line rather than a tissue contour.
    const px = ic.getImageData(Math.max(0, Math.floor(x)), Math.max(0, Math.floor(y)), 1, 1).data;
    drawMode = (px[0] >= 230 && px[1] >= 230 && px[2] >= 230) ? 'separator' : 'contour';
    drawing = true;
    currentPath = [[x, y]];
    strokeWandDisabled = e.metaKey || e.ctrlKey;   // hold Cmd/Ctrl to draw free-hand just this once
  });
  // Listen on window (not just the canvas) so the stroke keeps going when the
  // cursor briefly crosses the canvas edge, instead of cutting the contour
  // short there; coordinates are clamped to the canvas in canvasXY().
  window.addEventListener('mousemove', e => { if (drawing) { currentPath.push(canvasXY(e)); redraw(); } });
  window.addEventListener('mouseup',   finishDraw);
  window.addEventListener('blur',      finishDraw);   // e.g. alt-tab while dragging

  // ── separator selection (double-click) ────────────────────────────
  function _distToSeg(px, py, x1, y1, x2, y2) {
    const dx = x2 - x1, dy = y2 - y1, len2 = dx * dx + dy * dy;
    if (len2 < 1e-10) return Math.hypot(px - x1, py - y1);
    const t = Math.max(0, Math.min(1, ((px - x1) * dx + (py - y1) * dy) / len2));
    return Math.hypot(px - (x1 + t * dx), py - (y1 + t * dy));
  }

  drawCanvas.addEventListener('dblclick', e => {
    const [x, y] = canvasXY(e);
    let best = -1, bestDist = 12;
    separators.forEach((sep, i) => {
      const pts = sep.points;
      for (let j = 0; j < pts.length - 1; j++) {
        const d = _distToSeg(x, y, pts[j][0], pts[j][1], pts[j+1][0], pts[j+1][1]);
        if (d < bestDist) { bestDist = d; best = i; }
      }
    });
    selectedSepIdx = best;
    redraw();
    statusEl.textContent = best >= 0
      ? '✏ Separator #' + best + ' selected — Del/Backspace to delete, Esc to deselect.'
      : '';
  });

  async function finishDraw() {
    if (!drawing) return; drawing = false;
    if (currentPath.length <= 1) { currentPath = []; return; }

    const scribble = [...currentPath];
    const useWand  = wandCb.checked && !strokeWandDisabled;
    currentPath = [];

    if (drawMode === 'separator') {
      commitSeparator(scribble);
      return;
    }

    if (useWand) {
      statusEl.textContent = '🪄 Growing region…';
      try {
        const res = await growRegion(scribble);
        if (res.ok) {
          commitContour(res.points);
          return;
        }
        statusEl.textContent = '⚠ Wand: ' + res.error + ' — used raw scribble instead.';
      } catch (e) {
        statusEl.textContent = '⚠ Wand request failed: ' + e + ' — used raw scribble instead.';
      }
    } else if (strokeWandDisabled && wandCb.checked) {
      statusEl.textContent = '✋ Wand skipped (⌘/Ctrl held) — used free-hand contour.';
    }
    commitContour(scribble);
  }

  // ── buttons ───────────────────────────────────────────────────────
  document.getElementById('btn-save').addEventListener('click', () => {
    if (!contours.length) { statusEl.textContent = '⚠ No contour — draw one first.'; return; }
    doSave(contours.length - 1);
  });
  document.getElementById('btn-clear-last').addEventListener('click', async () => {
    if (currentPath.length) {
      currentPath = [];
      redraw();
      statusEl.textContent = 'In-progress stroke cleared.';
      return;
    }
    if (!contours.length) return;
    const last = contours.pop();
    redraw(); updateList();
    if (last.file) {
      statusEl.textContent = '🗑 Deleting ' + last.file + '…';
      await deleteFile(last.file);
      statusEl.textContent = '✅ Removed contour and deleted ' + last.file;
    } else {
      statusEl.textContent = 'Last (unsaved) contour removed.';
    }
  });
  document.getElementById('btn-clear-all').addEventListener('click', async () => {
    const sample = SAMPLES[imgIdx];
    contours = []; currentPath = [];
    redraw(); updateList();
    statusEl.textContent = '🗑 Deleting all saved contours for ' + sample + '…';
    await deleteAllForSample(sample);
    statusEl.textContent = '✅ All contours cleared for ' + sample + '.';
  });
  nextBtn.addEventListener('click', () => {
    if (imgIdx < SIZES.length - 1) { imgIdx++; loadImage(imgIdx); }
    else statusEl.textContent = '✓ All images annotated.';
  });
  prevBtn.addEventListener('click', () => {
    if (imgIdx > 0) { imgIdx--; loadImage(imgIdx); }
  });
  autosaveCb.addEventListener('change', () => {
    statusEl.textContent = autosaveCb.checked
      ? 'Autosave ON — saved automatically on mouse-up.'
      : 'Autosave OFF — click Save contour after drawing.';
  });
  wandCb.addEventListener('change', () => {
    statusEl.textContent = wandCb.checked
      ? 'Magic wand ON — scribbles snap to the tissue boundary.'
      : 'Magic wand OFF — free-hand drawing.';
  });
  rulerCb.addEventListener('change', () => { redraw(); });
  // ── keyboard ──────────────────────────────────────────────────────
  // JupyterLab intercepts keys in the capture phase on the window, before
  // any bubble-phase listener on the canvas can see them.  The only reliable
  // fix is to also register in the CAPTURE phase on window, check our own
  // hover flag, and stop propagation there so the notebook never sees the event.
  let canvasHovered = false;
  drawCanvas.addEventListener('mouseenter', () => {
    canvasHovered = true;
    try { if (Jupyter && Jupyter.keyboard_manager) Jupyter.keyboard_manager.disable(); } catch(_) {}
    drawCanvas.focus();
  });
  drawCanvas.addEventListener('mouseleave', () => {
    canvasHovered = false;
    try { if (Jupyter && Jupyter.keyboard_manager) Jupyter.keyboard_manager.enable(); } catch(_) {}
    drawCanvas.blur();
  });

  function _handleKey(e) {
    if (!canvasHovered) return;
    // Swallow the event entirely so the notebook/CodeMirror never sees it.
    e.stopPropagation();
    e.stopImmediatePropagation();
    if (e.key === 'ArrowRight') { e.preventDefault(); nextBtn.click(); }
    else if (e.key === 'ArrowLeft') { e.preventDefault(); prevBtn.click(); }
    else if (e.key === 'z') { e.preventDefault(); document.getElementById('btn-clear-last').click(); }
    else if (e.key === 's' && (e.ctrlKey || e.metaKey)) { e.preventDefault(); document.getElementById('btn-save').click(); }
    else if (e.key === 'w') { wandCb.checked = !wandCb.checked; wandCb.dispatchEvent(new Event('change')); }
    else if ((e.key === 'Delete' || e.key === 'Backspace') && selectedSepIdx >= 0) {
      e.preventDefault();
      separators.splice(selectedSepIdx, 1);
      selectedSepIdx = -1;
      redraw(); updateList();
      doSaveSeparators();
      statusEl.textContent = '🗑 Separator deleted.';
    } else if (e.key === 'Escape' && selectedSepIdx >= 0) {
      selectedSepIdx = -1;
      redraw();
      statusEl.textContent = 'Selection cleared.';
    }
  }
  // Register in CAPTURE phase (true) so we beat JupyterLab's own capture listeners.
  window.addEventListener('keydown', _handleKey, true);

  loadImage(0);
})();
</script>
"""
    )

    display(HTML(css + js))