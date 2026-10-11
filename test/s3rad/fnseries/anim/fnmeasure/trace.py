"""Trace the hydrides on a micrograph and measure F_n the way Plyasov et al.
(Phys. At. Nucl. 86 (2023) 2604) define it: total length of hydrides within
45 deg of the radial direction over the total length of all hydrides.

    python trace.py panel.png out.json
"""
import sys, json
import numpy as np
from PIL import Image
from scipy import ndimage as ndi
from skimage.morphology import skeletonize, remove_small_objects
from skimage.measure import approximate_polygon

src, out = sys.argv[1], sys.argv[2]
rgb = np.asarray(Image.open(src).convert('RGB')).astype(float)
g = rgb.mean(2)
H, W = g.shape

# ---- tube region and its centre (outer surface on top, inner at the bottom)
bright = ndi.binary_opening(g > 150, iterations=2)
top, bot, xs = [], [], []
for x in range(W):
    col = np.where(bright[:, x])[0]
    if len(col) > H * 0.4:
        xs.append(x); top.append(col.min()); bot.append(col.max())
xs, top, bot = map(np.array, (xs, top, bot))


def fit_circle(x, y):
    A = np.c_[2 * x, 2 * y, np.ones_like(x)]
    b = x ** 2 + y ** 2
    cx, cy, c = np.linalg.lstsq(A, b, rcond=None)[0]
    return cx, cy, np.sqrt(c + cx ** 2 + cy ** 2)
co = fit_circle(xs.astype(float), top.astype(float))
ci = fit_circle(xs.astype(float), bot.astype(float))
cx, cy = (co[0] + ci[0]) / 2, (co[1] + ci[1]) / 2
Ro = np.hypot(xs - cx, top - cy).mean()
Ri = np.hypot(xs - cx, bot - cy).mean()
print('outer circle', co, 'inner', ci, '-> centre (%.0f, %.0f), Ro %.1f Ri %.1f, wall %.1f px'
      % (cx, cy, Ro, Ri, Ro - Ri))

yy, xx = np.mgrid[0:H, 0:W]
rr = np.hypot(xx - cx, yy - cy)
inside = (rr < Ro - 4) & (rr > Ri + 4)

# ---- hydrides: darker than the local background
bg = ndi.median_filter(g, size=17)
THR = float(sys.argv[3]) if len(sys.argv) > 3 else 28
hyd = (bg - g > THR) & inside
hyd = remove_small_objects(hyd, max_size=6)
sk = skeletonize(hyd)
sk = remove_small_objects(sk, max_size=5, connectivity=2)

# ---- split the skeleton into branches at junctions, then into straight pieces
K = np.ones((3, 3)); K[1, 1] = 0
nb = ndi.convolve(sk.astype(int), K, mode='constant') * sk
junction = sk & (nb >= 3)
branches = sk & ~ndi.binary_dilation(junction, structure=np.ones((3, 3)))
lab, n = ndi.label(branches, structure=np.ones((3, 3)))
segs = []
for i in range(1, n + 1):
    pts = np.argwhere(lab == i)
    if len(pts) < 3:
        continue
    # order the pixels along the branch: start at an end, walk neighbours
    S = {tuple(p) for p in pts}
    deg = {p: sum((p[0] + dy, p[1] + dx) in S for dy in (-1, 0, 1) for dx in (-1, 0, 1)
                  if dy or dx) for p in S}
    start = min(S, key=lambda p: deg[p])
    path, seen = [start], {start}
    while True:
        y, x = path[-1]
        nxt = [(y + dy, x + dx) for dy in (-1, 0, 1) for dx in (-1, 0, 1)
               if (dy or dx) and (y + dy, x + dx) in S and (y + dy, x + dx) not in seen]
        if not nxt:
            break
        nxt.sort(key=lambda q: abs(q[0] - y) + abs(q[1] - x))
        path.append(nxt[0]); seen.add(nxt[0])
    # the junction pixel the branch was cut from: reconnect both ends to it
    P = np.array(path, float)
    for end in (0, -1):
        y, x = P[end].astype(int)
        jn = [(y + dy, x + dx) for dy in (-2, -1, 0, 1, 2) for dx in (-2, -1, 0, 1, 2)
              if 0 <= y + dy < H and 0 <= x + dx < W and junction[y + dy, x + dx]]
        if jn:
            q = np.array(min(jn, key=lambda q: abs(q[0] - y) + abs(q[1] - x)), float)
            P = np.vstack([q, P]) if end == 0 else np.vstack([P, q])
    poly = approximate_polygon(P, tolerance=1.3)
    if np.hypot(*np.diff(poly, axis=0).T).sum() < 4:      # the whole branch is too short
        continue
    for a, b in zip(poly[:-1], poly[1:]):
        L = np.hypot(*(b - a))
        if L < 0.5:
            continue
        m = (a + b) / 2
        r = np.array([m[1] - cx, m[0] - cy]); r /= np.linalg.norm(r)
        d = np.array([b[1] - a[1], b[0] - a[0]]) / L
        ang = np.degrees(np.arccos(min(1., abs(r @ d))))      # 0 = along R
        segs.append(dict(x0=a[1], y0=a[0], x1=b[1], y1=b[0], L=L, ang=ang))

L = np.array([s['L'] for s in segs]); A = np.array([s['ang'] for s in segs])
Fn = L[A <= 45].sum() / L.sum()
hist = [round(100 * L[(A >= a) & (A < a + 5)].sum() / L.sum(), 1) for a in range(0, 90, 5)]
print('segments %d, total length %.0f px, F_n = %.3f' % (len(segs), L.sum(), Fn))
print('L(alpha) %, 5 deg bins from R:', hist)
json.dump(dict(W=W, H=H, cx=cx, cy=cy, Ro=Ro, Ri=Ri, Fn=Fn, segs=segs), open(out, 'w'))
Image.fromarray((hyd * 255).astype(np.uint8)).save(out.replace('.json', '_mask.png'))
