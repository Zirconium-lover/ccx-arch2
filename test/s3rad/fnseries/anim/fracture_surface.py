"""Pictures of the fracture surface of a broken run (slide "Поверхности
разрушения"): the block at the end of the run is split into its two halves
(tets connected through faces; a fully debonded interface does not hold),
and the surface of the left half that faces the crack is drawn in 3D,
coloured by how it was made:

    ductile   wall of a deleted matrix element      (metal)
    cleavage  wall of a deleted hydride element     (orange)
    debond    a face of the debonded interface      (green, as in the GIFs)

The outer faces of the half are dark.  Lambert shading, painter's order,
back faces culled.

    python fracture_surface.py out_dir [run ...]
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import to_rgb
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from fndata import Run, FACES

COL = {'ductile': '#A9B4C2', 'cleavage': '#FFA040', 'debond': '#3CC9A0', 'outer': '#2C333C'}
VIEW = np.array([1.0, -0.55, 0.62])      # towards the camera: from +x (the crack side), front, above


D_CUT, DEB_CUT = 0.5, 0.2      # at the end the load is ~0: what has lost half its strength is crack


def halves(r, t_end):
    """Connected parts of the tets that hold at the end of the run (not
    deleted, D < D_CUT); an interface face links its two tets only while the
    facet has lost less than DEB_CUT of its strength.  Returns labels per tet
    (-1: part of the crack) and the labels of the parts at x = 0 and x = L."""
    alive = (r.tdel > t_end) & (np.nan_to_num(r.D_end, nan=1.) < D_CUT)
    a, b, _ = r._pairs()
    # interface pairs: the two tets carry different phases (matrix / hydride)
    s = r.at(r.eps[-1])
    keyf = {}
    geo = r.geo
    for f, nodes in enumerate(r.coh[:, :3]):
        keyf[tuple(sorted(geo[nodes]))] = f
    F = r.T[:, FACES]
    link = alive[a] & alive[b]
    ifc = (r.phase[a] == 0) != (r.phase[b] == 0)
    if np.any(ifc):
        # find the face shared by a and b, look its facet up
        idx = np.flatnonzero(link & ifc)
        for j in idx:
            fa = {tuple(sorted(geo[q])) for q in F[a[j]]}
            fb = {tuple(sorted(geo[q])) for q in F[b[j]]}
            k = (fa & fb).pop()
            f = keyf.get(k)
            if f is not None and s.debond[f] >= DEB_CUT:
                link[j] = False
    n = len(r.T)
    g = coo_matrix((np.ones(link.sum()), (a[link], b[link])), shape=(n, n))
    _, lab = connected_components(g, directed=False)
    lab = np.where(alive, lab, -1)
    cen = r.X0[r.T].mean(1)
    left = np.bincount(lab[(cen[:, 0] < r.lo[0] + 0.3) & alive]).argmax()
    right = np.bincount(lab[(cen[:, 0] > r.hi[0] - 0.3) & alive]).argmax()
    return lab, left, right, s


def surface_of(r, part, lab, s):
    """Faces of one part with their kind (outer / ductile / cleavage / debond)."""
    mask = lab == part
    faces, owner, kind = r.surface(mask, s.X)
    out = np.full(len(faces), 'outer', object)
    inner = kind == 1
    if inner.any():
        # what lies across each inner face: a deleted tet (which phase) or another part
        geo = r.geo
        key = {}
        F = r.T[:, FACES].reshape(-1, 3)
        own = np.repeat(np.arange(len(r.T)), 4)
        ks = np.sort(geo[F], axis=1)
        for i, k in enumerate(map(tuple, ks)):
            key.setdefault(k, []).append(own[i])
        for i in np.flatnonzero(inner):
            k = tuple(sorted(geo[faces[i]]))
            nb = [e for e in key.get(k, []) if e != owner[i]]
            if not nb:
                out[i] = 'outer'
                continue
            e = nb[0]
            if lab[e] < 0:                                   # deleted neighbour
                out[i] = 'ductile' if r.phase[e] == 0 else 'cleavage'
            else:                                            # across a debonded interface
                out[i] = 'debond'
    return faces, out


def render(r, out_png, px=(1200, 900)):
    t_end = r.t[-1]
    lab, left, right, s = halves(r, t_end)
    faces, kind = surface_of(r, left, lab, s)
    X = s.X - (r.lo + r.hi) / 2
    w = VIEW / np.linalg.norm(VIEW)
    up = np.array([0., 0., 1.])
    u = np.cross(up, w); u /= np.linalg.norm(u)
    v = np.cross(w, u)
    P = X[faces]
    nrm = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])
    nrm /= np.linalg.norm(nrm, axis=1, keepdims=True) + 1e-30
    front = nrm @ w > 0                                       # back faces culled
    P, nrm, kind = P[front], nrm[front], kind[front]
    depth = (P.mean(1) @ w)
    o = np.argsort(depth)
    P, nrm, kind = P[o], nrm[o], kind[o]
    light = np.array([0.35, -0.5, 0.8]); light /= np.linalg.norm(light)
    lam = 0.38 + 0.62 * np.clip(nrm @ light, 0, 1)
    base = np.array([to_rgb(COL[k]) for k in kind])
    cols = np.clip(base * lam[:, None], 0, 1)
    xy = np.stack([P @ u, P @ v], -1)
    fig = plt.figure(figsize=(px[0] / 100, px[1] / 100), dpi=100)
    fig.patch.set_alpha(0)
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_axis_off(); ax.patch.set_alpha(0)
    ax.add_collection(PolyCollection(xy, facecolors=cols, edgecolors=cols, linewidths=0.3))
    lo, hi = xy.reshape(-1, 2).min(0), xy.reshape(-1, 2).max(0)
    c, half = (lo + hi) / 2, (hi - lo) / 2 * 1.04
    k = px[0] / px[1]
    if half[0] / half[1] < k:
        half[0] = half[1] * k
    else:
        half[1] = half[0] / k
    ax.set_xlim(c[0] - half[0], c[0] + half[0]); ax.set_ylim(c[1] - half[1], c[1] + half[1])
    fig.savefig(out_png, transparent=True)
    plt.close(fig)
    areas = {}
    for kk in ('ductile', 'cleavage', 'debond'):
        m = kind == kk
        areas[kk] = float(0.5 * np.linalg.norm(np.cross(P[m, 1] - P[m, 0], P[m, 2] - P[m, 0]), axis=1).sum())
    return areas


if __name__ == '__main__':
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    for name in sys.argv[2:] or ['fn0', 'fn25v1e-3', 'fn60', 'fn100']:
        r = Run(name)
        print(name, render(r, os.path.join(out, 'fs_%s.png' % name)))


def projected_brittle(r, px_per_mm=200):
    """Fractography view: the fracture surface of each half seen along T (the
    load axis), i.e. projected on the R-L plane with hidden parts removed.
    Brittle = debonding + hydride cleavage, ductile = matrix.  Returns the
    brittle share of the projected area (mean of the two halves) and the
    per-half shares."""
    t_end = r.t[-1]
    lab, left, right, s = halves(r, t_end)
    X = r.X0
    shares = []
    for part, sgn in ((left, 1.), (right, -1.)):          # the left half is seen from +x
        faces, kind = surface_of(r, part, lab, s)
        P = X[faces]
        nrm = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])
        front = sgn * nrm[:, 0] > 0
        keep = front & (kind != 'outer')
        P, kind = P[keep], kind[keep]
        o = np.argsort(sgn * P[:, :, 0].mean(1))           # far first, near last
        P, kind = P[o], kind[o]
        cid = {'ductile': (1., 0., 0.), 'debond': (0., 1., 0.), 'cleavage': (0., 0., 1.)}
        cols = np.array([cid[k] for k in kind])
        wpx, hpx = int((r.hi[1] - r.lo[1]) * px_per_mm), int((r.hi[2] - r.lo[2]) * px_per_mm)
        fig = plt.figure(figsize=(wpx / 100, hpx / 100), dpi=100)
        ax = fig.add_axes([0, 0, 1, 1]); ax.set_axis_off()
        ax.add_collection(PolyCollection(P[:, :, 1:], facecolors=cols, edgecolors='none', antialiased=False))
        ax.set_xlim(r.lo[1], r.hi[1]); ax.set_ylim(r.lo[2], r.hi[2])
        fig.canvas.draw()
        img = np.asarray(fig.canvas.buffer_rgba())[..., :3].astype(int)
        plt.close(fig)
        duct = np.sum((img[..., 0] == 255) & (img[..., 1] == 0) & (img[..., 2] == 0))
        brit = np.sum((img[..., 0] == 0) & ((img[..., 1] == 255) ^ (img[..., 2] == 255)))
        shares.append(brit / max(duct + brit, 1))
    return float(np.mean(shares)), shares
