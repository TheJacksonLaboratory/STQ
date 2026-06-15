"""
contour_annotator.py — freehand contour annotation for Jupyter notebooks.

Save format (one file per contour):
  <savepath>/<sample>.oid<N>.json
  {
    "0": { "points": [0.104, 0.244, ...] },   <- x, normalised 0-1
    "1": { "points": [0.160, 0.328, ...] }    <- y, normalised 0-1
  }
"""

import base64, json, os, socket, threading
from http.server import BaseHTTPRequestHandler, HTTPServer
from io import BytesIO

import numpy as np
from PIL import Image
from IPython.display import display, HTML


def annotate_contours(
    samples: list[str],
    images: list[np.ndarray],
    canvas_max_dim: int = 800,
    savepath: str = "./contours/",
):
    assert len(samples) == len(images), "samples and images must have the same length"
    os.makedirs(savepath, exist_ok=True)

    # ── encode images ────────────────────────────────────────────────────────
    def _b64(arr):
        if arr.dtype != np.uint8:
            arr = np.clip(arr, 0, 255).astype(np.uint8)
        buf = BytesIO()
        Image.fromarray(arr).save(buf, format="PNG")
        return base64.b64encode(buf.getvalue()).decode()

    b64s  = [_b64(img) for img in images]
    sizes = [(int(img.shape[1]), int(img.shape[0])) for img in images]

    # ── tiny HTTP server ─────────────────────────────────────────────────────
    # JS POSTs {file, data} JSON here; Python writes the file and replies 200.
    # Runs in a daemon thread — dies automatically when the notebook kernel stops.

    with socket.socket() as s:
        s.bind(('', 0))
        port = s.getsockname()[1]

    class Handler(BaseHTTPRequestHandler):
        def do_OPTIONS(self):          # preflight CORS
            self._ok(b'')
        def do_POST(self):
            n    = int(self.headers.get('Content-Length', 0))
            body = self.rfile.read(n)
            try:
                payload = json.loads(body)
                fname   = payload['file']
                data    = payload['data']
                dname   = os.path.dirname(fname)
                if dname:
                    os.makedirs(dname, exist_ok=True)
                with open(fname, 'w') as fh:
                    json.dump(data, fh, indent=2)
                self._ok(b'{"ok":true}')
            except Exception as e:
                self._ok(json.dumps({'ok': False, 'error': str(e)}).encode())
        def _ok(self, body):
            self.send_response(200)
            self.send_header('Content-Type',                 'application/json')
            self.send_header('Content-Length',               str(len(body)))
            self.send_header('Access-Control-Allow-Origin',  '*')
            self.send_header('Access-Control-Allow-Methods', 'POST, OPTIONS')
            self.send_header('Access-Control-Allow-Headers', 'Content-Type')
            self.end_headers()
            if body:
                self.wfile.write(body)
        def log_message(self, *a):
            pass

    host = socket.gethostbyname(socket.gethostname())
    srv = HTTPServer((host, port), Handler)
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
  #ca-draw-canvas { position:absolute; top:0; left:0; }
  #ca-toolbar { display:flex; gap:8px; align-items:center; flex-wrap:wrap; }
  .ca-btn {
    font-family:var(--font); font-size:12px; letter-spacing:.06em;
    padding:6px 14px; border:none; border-radius:var(--radius);
    cursor:pointer; transition:filter .15s;
  }
  .ca-btn:hover  { filter:brightness(1.15); }
  .ca-btn:active { filter:brightness(.9); }
  #btn-save       { background:var(--accent);  color:#fff; }
  #btn-next       { background:var(--accent2); color:var(--text); }
  #btn-clear-last { background:#2a2a4a; color:var(--text); }
  #btn-clear-all  { background:#2a2a4a; color:var(--muted); }
  #ca-autosave-wrap {
    display:flex; align-items:center; gap:6px;
    font-size:11px; color:var(--muted); margin-left:4px;
  }
  #ca-autosave-wrap input { accent-color:var(--accent); width:14px; height:14px; cursor:pointer; }
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
    <canvas id="ca-draw-canvas"></canvas>
  </div>
  <div id="ca-toolbar">
    <button class="ca-btn" id="btn-save">Save contour</button>
    <button class="ca-btn" id="btn-clear-last">Clear last</button>
    <button class="ca-btn" id="btn-clear-all">Clear all</button>
    <button class="ca-btn" id="btn-next">Next →</button>
    <label id="ca-autosave-wrap">
      <input type="checkbox" id="ca-autosave" checked> autosave
    </label>
  </div>
  <div id="ca-status">Loading…</div>
  <div id="ca-contour-list"></div>
</div>
"""

    js = (
        "<script>\n(function(){\n"
        f"  const IMAGES   = {json.dumps(b64s)};\n"
        f"  const SIZES    = {json.dumps(sizes)};\n"
        f"  const SAMPLES  = {json.dumps(samples)};\n"
        f"  const SAVEPATH = {json.dumps(savepath)};\n"
        f"  const MAX_DIM  = {canvas_max_dim};\n"
        f"  const SRV_PORT = {port};\n"
        f"  const SRV_HOST = {json.dumps(host)};\n"
        + r"""
  let imgIdx      = 0;
  let contours    = [];
  let currentPath = [];
  let drawing     = false;
  let oidCounters = {};

  const imgCanvas   = document.getElementById('ca-img-canvas');
  const drawCanvas  = document.getElementById('ca-draw-canvas');
  const ic = imgCanvas.getContext('2d');
  const dc = drawCanvas.getContext('2d');
  const sampleLabel = document.getElementById('ca-sample-label');
  const counterEl   = document.getElementById('ca-counter');
  const statusEl    = document.getElementById('ca-status');
  const contourList = document.getElementById('ca-contour-list');
  const autosaveCb  = document.getElementById('ca-autosave');

  // ── save via local HTTP server ────────────────────────────────────
  async function doSave(pts) {
    const sample = SAMPLES[imgIdx];
    if (!(sample in oidCounters)) oidCounters[sample] = 0;
    const oid   = oidCounters[sample]++;
    const fname = SAVEPATH + sample + '.oid' + oid + '.json';
    const [cw, ch] = [drawCanvas.width, drawCanvas.height];
    const data  = {
      "0": { "points": pts.map(p => parseFloat((p[0] / cw).toFixed(4))) },
      "1": { "points": pts.map(p => parseFloat((p[1] / ch).toFixed(4))) }
    };
    statusEl.textContent = '💾 Saving → ' + fname + '…';
    try {
      const resp = await fetch('http://' + SRV_HOST + ':' + SRV_PORT, {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({file: fname, data: data})
      });
      const json = await resp.json();
      if (json.ok) {
        statusEl.textContent = '✅ Saved → ' + fname + '  (' + pts.length + ' pts)';
      } else {
        statusEl.textContent = '⚠ Server error: ' + json.error;
      }
    } catch(e) {
      statusEl.textContent = '⚠ Fetch failed: ' + e;
    }
  }

  // ── geometry ──────────────────────────────────────────────────────
  function scaledSize(w, h) {
    const r = Math.min(MAX_DIM / w, MAX_DIM / h, 1);
    return [Math.round(w * r), Math.round(h * r)];
  }

  function loadImage(idx) {
    const [w, h]   = SIZES[idx];
    const [cw, ch] = scaledSize(w, h);
    imgCanvas.width = drawCanvas.width = cw;
    imgCanvas.height = drawCanvas.height = ch;
    const img = new Image();
    img.onload = () => ic.drawImage(img, 0, 0, cw, ch);
    img.src = 'data:image/png;base64,' + IMAGES[idx];
    sampleLabel.textContent = SAMPLES[idx];
    counterEl.textContent   = (idx + 1) + ' / ' + IMAGES.length;
    contours = []; currentPath = [];
    redraw(); updateList();
    statusEl.textContent = autosaveCb.checked
      ? 'Draw a contour — saved automatically on mouse-up.'
      : 'Draw a contour, then click Save contour.';
  }

  function canvasXY(e) {
    const r = drawCanvas.getBoundingClientRect();
    return [e.clientX - r.left, e.clientY - r.top];
  }

  // ── drawing ───────────────────────────────────────────────────────
  function redraw() {
    dc.clearRect(0, 0, drawCanvas.width, drawCanvas.height);
    contours.forEach((pts, i) => {
      if (pts.length < 2) return;
      [['rgba(255,255,255,0.65)', 4], ['#000', 2]].forEach(([col, lw]) => {
        dc.beginPath(); dc.moveTo(pts[0][0], pts[0][1]);
        pts.slice(1).forEach(p => dc.lineTo(p[0], p[1]));
        dc.closePath(); dc.strokeStyle = col; dc.lineWidth = lw;
        dc.setLineDash([]); dc.stroke();
      });
      dc.font = 'bold 11px monospace'; dc.setLineDash([]);
      dc.lineWidth = 3; dc.strokeStyle = 'white';
      dc.strokeText('#' + i, pts[0][0] + 4, pts[0][1] - 4);
      dc.fillStyle = '#000'; dc.fillText('#' + i, pts[0][0] + 4, pts[0][1] - 4);
    });
    if (currentPath.length > 1) {
      dc.beginPath(); dc.moveTo(currentPath[0][0], currentPath[0][1]);
      currentPath.slice(1).forEach(p => dc.lineTo(p[0], p[1]));
      dc.strokeStyle = '#000'; dc.lineWidth = 1.5;
      dc.setLineDash([5, 4]); dc.stroke(); dc.setLineDash([]);
    }
  }

  function updateList() {
    contourList.textContent = !contours.length ? '' :
      contours.map((c, i) => '#' + i + ': ' + c.length + ' pts').join('  ·  ');
  }

  // ── mouse ─────────────────────────────────────────────────────────
  drawCanvas.addEventListener('mousedown', e => { drawing = true; currentPath = [canvasXY(e)]; });
  drawCanvas.addEventListener('mousemove', e => { if (drawing) { currentPath.push(canvasXY(e)); redraw(); } });
  function finishDraw() {
    if (!drawing) return; drawing = false;
    if (currentPath.length > 1) {
      contours.push([...currentPath]); currentPath = []; redraw(); updateList();
      if (autosaveCb.checked) doSave(contours[contours.length - 1]);
      else statusEl.textContent = 'Contour drawn. Click Save contour to write file.';
    } else { currentPath = []; }
  }
  drawCanvas.addEventListener('mouseup',    finishDraw);
  drawCanvas.addEventListener('mouseleave', finishDraw);

  // ── buttons ───────────────────────────────────────────────────────
  document.getElementById('btn-save').addEventListener('click', () => {
    if (!contours.length) { statusEl.textContent = '⚠ No contour — draw one first.'; return; }
    doSave(contours[contours.length - 1]);
  });
  document.getElementById('btn-clear-last').addEventListener('click', () => {
    if (currentPath.length) {
      currentPath = [];
    } else if (contours.length) {
      contours.pop();
      const sample = SAMPLES[imgIdx];
      oidCounters[sample] = Math.max((oidCounters[sample] || 1) - 1, 0);
    }
    redraw(); updateList(); statusEl.textContent = 'Last contour removed.';
  });
  document.getElementById('btn-clear-all').addEventListener('click', () => {
    contours = []; currentPath = []; redraw(); updateList();
    statusEl.textContent = 'All contours cleared.';
    oidCounters[SAMPLES[imgIdx]] = 0;
  });
  document.getElementById('btn-next').addEventListener('click', () => {
    if (imgIdx < IMAGES.length - 1) { imgIdx++; loadImage(imgIdx); }
    else statusEl.textContent = '✓ All images annotated.';
  });
  autosaveCb.addEventListener('change', () => {
    statusEl.textContent = autosaveCb.checked
      ? 'Autosave ON — saved automatically on mouse-up.'
      : 'Autosave OFF — click Save contour after drawing.';
  });
  document.addEventListener('keydown', e => {
    if (e.key === 'ArrowRight') document.getElementById('btn-next').click();
    if (e.key === 's' && (e.ctrlKey || e.metaKey)) { e.preventDefault(); document.getElementById('btn-save').click(); }
    if (e.key === 'z' && (e.ctrlKey || e.metaKey)) { e.preventDefault(); document.getElementById('btn-clear-last').click(); }
  });

  loadImage(0);
})();
</script>
"""
    )

    display(HTML(css + js))