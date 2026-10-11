"""Contact sheet of the packed runs: a check that the data draw right, and a
menu of what the animations can use.

    python data_preview.py [out.png]
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from fndata import Run, SERIES, project

RADC, TANC, ACCC, MATC = '#E8604C', '#2CB1D9', '#FFD24A', '#5A6B80'
PHI, TH = np.radians(62), np.radians(-62)


def shade(faces, X, base, light=(0.3, -0.5, 0.8)):
    a, b, c = X[faces[:, 0]], X[faces[:, 1]], X[faces[:, 2]]
    n = np.cross(b - a, c - a)
    n /= np.linalg.norm(n, axis=1, keepdims=True) + 1e-30
    k = 0.35 + 0.65 * np.abs(n @ np.array(light) / np.linalg.norm(light))
    return np.clip(np.asarray(base)[:, :3] * k[:, None], 0, 1)


def view3d(ax, r, eps, title):
    s = r.at(eps)
    X = s.X - (r.lo + r.hi) / 2
    groups = []
    hyd = s.alive & (r.phase > 0)
    f, o, _ = r.surface(hyd, s.X)
    col = np.array([matplotlib.colors.to_rgb(RADC if r.plate_axis[r.plate[i]] == 0 else TANC) for i in o])
    groups.append((f, col))
    f, o, k = r.surface(s.alive, s.X)
    crack = f[k == 1]
    groups.append((crack, np.tile(matplotlib.colors.to_rgb(ACCC), (len(crack), 1))))
    F = np.vstack([g[0] for g in groups]); C = np.vstack([g[1] for g in groups])
    P2, z = project(X, PHI, TH)
    depth = z[F].mean(1)
    order = np.argsort(depth)
    polys = P2[F[order]]
    ax.add_collection(PolyCollection(polys, facecolors=shade(F[order], X, C[order]), edgecolors='none'))
    # box edges
    lo, hi = r.lo - (r.lo + r.hi) / 2, r.hi - (r.lo + r.hi) / 2
    cs = np.array([[x, y, zz] for x in (lo[0], hi[0]) for y in (lo[1], hi[1]) for zz in (lo[2], hi[2])])
    for i in range(8):
        for j in range(i + 1, 8):
            if np.sum(cs[i] != cs[j]) == 1:
                p, _ = project(cs[[i, j]], PHI, TH)
                ax.plot(p[:, 0], p[:, 1], color='#8899AA', lw=0.8)
    ax.set_aspect('equal'); ax.autoscale(); ax.axis('off')
    ax.set_title(title, color='w', fontsize=10)


def section(ax, r, eps, field, title, vmax=1.0, cmap='inferno'):
    s = r.at(eps, deform=1.0)
    polys, own, pv = r.slice('z', (r.lo[2] + r.hi[2]) / 2, s, field)
    v = np.array([np.mean(p) for p in pv])
    pc = PolyCollection([p[:, :2] for p in polys], array=v, cmap=cmap, edgecolors='none')
    pc.set_clim(0, vmax)
    ax.add_collection(pc)
    hp, ho, _ = r.slice('z', (r.lo[2] + r.hi[2]) / 2, s, None, s.alive & (r.phase > 0))
    ax.add_collection(PolyCollection([p[:, :2] for p in hp], facecolors='none',
                                     edgecolors=[RADC if r.plate_axis[r.plate[i]] == 0 else TANC for i in ho], lw=0.3))
    ax.set_aspect('equal'); ax.autoscale(); ax.axis('off')
    ax.set_title(title, color='w', fontsize=10)
    return pc


def main(out):
    plt.style.use('dark_background')
    fig = plt.figure(figsize=(16, 13), dpi=110)
    gs = fig.add_gridspec(3, 4, hspace=0.25, wspace=0.12)
    runs = [Run(fn) for fn, _ in SERIES]
    cm = plt.get_cmap('coolwarm')
    ax = fig.add_subplot(gs[0, :2])
    for r in runs:
        ax.plot(100 * r.curve[:, 0], r.curve[:, 1], color=cm(r.k_radial / 12), lw=2,
                label='k = %d / 12' % r.k_radial)
    ax.set_xlabel('ε, %'); ax.set_ylabel('σ, МПа'); ax.legend(fontsize=8, ncol=2)
    ax.set_title('кривые σ–ε (curve)', fontsize=10)
    ax = fig.add_subplot(gs[0, 2:])
    k = np.array([r.k_radial / 12 for r in runs])
    M = [r.mechanisms(r.eps[-1]) for r in runs]
    for key, c, lab in (('matrix', MATC, 'удалено матрицы'), ('hydride', ACCC, 'удалено гидрида'),
                        ('debond_rad', RADC, 'отслоение рад.'), ('debond_tan', TANC, 'отслоение танг.')):
        ax.plot(k, [m[key] for m in M], 'o-', color=c, label=lab)
    ax.set_xlabel('Fn = k/12'); ax.set_ylabel('доля'); ax.legend(fontsize=8)
    ax.set_title('механизмы в конце расчёта (mechanisms)', fontsize=10)
    for j, (r, frac) in enumerate(((runs[0], 1.0), (runs[0], 0.6), (runs[-1], 1.0), (runs[-1], 0.6))):
        e = r.eps_peak + frac * (r.eps[-1] - r.eps_peak)
        view3d(fig.add_subplot(gs[1, j]), r, e, '%s, ε = %.1f %%: трещина + пластины' % (r.name, 100 * e))
    r = runs[4]
    for j, (e, fld, vm, cmap) in enumerate(((r.eps_peak, 'omega', 1.0, 'inferno'), (r.eps[-1], 'omega', 1.0, 'inferno'),
                                            (r.eps_peak, 'peeq', 0.3, 'viridis'), (r.eps_peak, 'sxx', 900, 'magma'))):
        a = fig.add_subplot(gs[2, j])
        pc = section(a, r, e, fld, '%s, сечение z = середина, %s, ε = %.1f %%' % (r.name, fld, 100 * e), vm, cmap)
        fig.colorbar(pc, ax=a, fraction=0.04, pad=0.01)
    fig.savefig(out, facecolor='black')


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'preview.png')
