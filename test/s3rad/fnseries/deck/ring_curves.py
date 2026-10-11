"""Tensile diagrams of the ring specimens (Saburov et al., YaFI 2025, 16 (4)
442-454, fig. 8: load F against crosshead travel x, after hydride
reorientation under 0...140 MPa) redrawn in the deck's style.

The curves are read off the figure's embedded raster (page 8): the grid
gives the scales (x: 1...4 mm, F: 500...2000 N), each curve is found by its
legend colour.  Before fracture the curves lie on one another and only the
topmost is seen, so a curve hidden in a column takes the common course
there; its fracture is the long vertical run of its colour.

    python ring_curves.py saburov2025.pdf out.png [out.csv]
"""
import os, subprocess, sys, tempfile
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

STRESS = [0, 50, 70, 90, 110, 140]                  # MPa, legend order
INK = np.array([(2, 3, 3), (236, 4, 137), (238, 38, 56), (239, 96, 33), (216, 150, 30), (206, 219, 32)], float)
GRID_X = (192.5, 175.5)                             # px of x = 1 mm, px per mm
GRID_F = (587.0, 141.5)                             # px of F = 500 N, px per 500 N
LEGEND = (330, 290)                                 # legend box: x < 330, y < 290 px

BG, AX, TXT = '#000000', '#7C848E', '#B8BEC6'
TAN, MID, RAD = np.array([0x2C, 0xB1, 0xD9]), np.array([0xE6, 0xE6, 0xE6]), np.array([0xE8, 0x60, 0x4C])


def raster(pdf):
    with tempfile.TemporaryDirectory() as d:
        subprocess.run(['pdfimages', '-png', '-f', '8', '-l', '8', pdf, os.path.join(d, 'p')], check=True)
        return np.asarray(Image.open(os.path.join(d, 'p-000.png')).convert('RGB')).astype(float)


def to_xF(px, py):
    return (px - GRID_X[0]) / GRID_X[1] + 1, (GRID_F[0] - py) / GRID_F[1] * 500 + 500


def smooth(y, n):
    return np.convolve(np.r_[[y[0]] * (n // 2), y, [y[-1]] * (n // 2)], np.ones(n) / n, 'valid')


def digitize(a):
    H, W, _ = a.shape
    mx, mn = a.max(2), a.min(2)
    ink = ((mx - mn > 70) | (mx < 80)) & (mn < 200)                 # coloured or black, not grid / paper
    ink[:LEGEND[1], :LEGEND[0]] = False
    d = np.linalg.norm(a[..., None, :] - INK[None, None], axis=-1)
    lab = np.where(ink & (d.min(-1) < 90), d.argmin(-1), -1)

    # the common course (all curves lie on one another up to x ~ 2.5 mm): centres of the ink band,
    # scanned by columns where it is flat and by rows where it is steep, drops left out
    x0 = int(GRID_X[0] - GRID_X[1]) + 2
    band = lab >= 0
    drops = np.zeros_like(band)
    for x in range(W):                                               # long vertical runs: the drops
        r = np.flatnonzero(band[:, x])
        if len(r):
            for seg in np.split(r, np.flatnonzero(np.diff(r) > 1) + 1):
                if len(seg) > 60:                                    # the band itself runs < 30 px
                    drops[seg, x] = True
    band &= ~drops
    near_drop = np.convolve(drops.any(0), np.ones(5), 'same') > 0     # columns beside a drop
    band[:, near_drop] = False
    # up to the peak the course rises in x and F alike, so every anti-diagonal x - y = u crosses the band
    # in one short piece, steep part or flat: its middle is a point of the centre line
    xe = int(GRID_X[0] + 1.9 * GRID_X[1])                            # x = 2.9 mm
    ys_, xs_ = np.nonzero(band[:, :xe])
    u = xs_ - ys_
    pts = []
    for uu in np.unique(u):
        sel = u == uu
        xx, yy = xs_[sel], ys_[sel]
        o = np.argsort(xx)
        xx, yy = xx[o], yy[o]
        piece = np.split(np.arange(len(xx)), np.flatnonzero(np.diff(xx) > 2) + 1)[0]
        if len(piece) <= 16:
            pts.append((xx[piece].mean(), yy[piece].mean()))
    pts = np.array(sorted(pts))
    grid = np.arange(x0, W)
    cy = np.interp(grid, pts[:, 0], smooth(pts[:, 1], 5))
    split = GRID_X[0] + 1.5 * GRID_X[1]                             # x = 2.5 mm: the curves part
    out = []
    for k in range(len(STRESS)):
        m = lab == k
        runs = []
        for x in range(W):
            r = np.flatnonzero(m[:, x])
            seg = np.split(r, np.flatnonzero(np.diff(r) > 2) + 1)[0] if len(r) else r
            runs.append((seg.mean() if len(seg) else None, len(seg)))
        xf = int(np.argmax([n for _, n in runs]))                   # fracture: the long vertical run
        tops = []                                                    # where the drop starts: the load at fracture
        for x in range(xf - 3, xf + 4):
            r = np.flatnonzero(m[:, x])
            for seg in np.split(r, np.flatnonzero(np.diff(r) > 1) + 1) if len(r) else []:
                if len(seg) > 60:
                    tops.append(seg[0])
        y_f = min(tops) + 1.5
        own = np.array([x for x, (y, n) in enumerate(runs)              # own colour past 2.5 mm
                        if y is not None and n < 12 and split <= x < xf - 6])
        xs = np.arange(x0, xf + 1)
        ys = np.interp(xs, grid, cy)
        if len(own) > 5:
            oy = smooth(np.array([runs[x][0] for x in own]), 5)
            sel = xs >= own[0]
            w = np.clip((xs[sel] - own[0]) / 15., 0, 1)              # blend into the own course
            ys[sel] = (1 - w) * ys[sel] + w * np.interp(xs[sel], own, oy)
            tail = own[-8:]                                          # up to the drop: the last slope
            c = np.polyfit(tail, oy[-8:], 1)
            ys[xs > own[-1]] = np.polyval(c, xs[xs > own[-1]])
        w = np.clip(1 - (xs[-1] - xs) / 10., 0, 1)                    # the last 10 px lead into the drop
        ys = (1 - w) * ys + w * (ys + (y_f - ys[-1]))
        x_mm, F = to_xF(xs, ys)
        out.append((STRESS[k], x_mm, F))
    return out


def colour(t):
    c = TAN + (MID - TAN) * t / 0.5 if t <= 0.5 else MID + (RAD - MID) * (t - 0.5) / 0.5
    return c / 255


def plot(curves, out_png):
    # a serif with Cyrillic close to the deck's Cambria (Caladea has no Cyrillic)
    for f in ('Regular', 'Italic', 'Bold'):
        p = '/usr/share/fonts/truetype/noto/NotoSerif-%s.ttf' % f
        if os.path.exists(p):
            font_manager.fontManager.addfont(p)
    plt.rcParams.update({'font.family': 'Noto Serif', 'mathtext.fontset': 'custom',
                         'mathtext.rm': 'Noto Serif', 'mathtext.it': 'Noto Serif:italic',
                         'mathtext.bf': 'Noto Serif:bold', 'mathtext.cal': 'Noto Serif',
                         'mathtext.sf': 'Noto Serif', 'font.size': 17})
    fig = plt.figure(figsize=(7.2, 4.6), dpi=300)
    fig.patch.set_facecolor(BG)
    ax = fig.add_axes([0.13, 0.15, 0.6, 0.78])
    ax.set_facecolor(BG)
    for i, (s, x, F) in enumerate(curves):
        c = colour(s / STRESS[-1])
        # as in the paper, the curve that breaks first lies on top of the common course
        ax.plot(x, F, color=c, lw=2.6, solid_capstyle='round', zorder=2 + i, label='%d МПа' % s)
        ax.plot([x[-1], x[-1]], [F[-1], 0], color=c, lw=1.3, alpha=0.7, zorder=1)
        ax.plot(x[-1], F[-1], 'o', color=c, ms=6, zorder=20)
    ax.set_xlim(0, 5); ax.set_ylim(0, 2500)
    ax.set_xticks(range(6)); ax.set_yticks(range(0, 2501, 500))
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    for sp in ('left', 'bottom'):
        ax.spines[sp].set_color(AX); ax.spines[sp].set_linewidth(1.4)
    ax.tick_params(colors=TXT, labelsize=15, length=5, width=1.2, color=AX)
    ax.set_xlabel(r'$x$, мм', color=TXT, fontsize=19, labelpad=6)
    ax.set_ylabel(r'$F$, Н', color=TXT, fontsize=19, labelpad=8)
    leg = ax.legend(loc='upper left', bbox_to_anchor=(1.03, 1.0), frameon=False, fontsize=15,
                    title='напряжение\nпереориентации', title_fontsize=15, labelcolor=TXT,
                    handlelength=1.6, borderaxespad=0)
    leg.get_title().set_color(TXT); leg._legend_box.align = 'left'
    fig.savefig(out_png, facecolor=BG, bbox_inches='tight', pad_inches=0.12)
    plt.close(fig)


if __name__ == '__main__':
    a = raster(sys.argv[1])
    cv = digitize(a)
    plot(cv, sys.argv[2])
    if len(sys.argv) > 3:
        with open(sys.argv[3], 'w') as f:
            f.write('stress_MPa,x_mm,F_N\n')
            for s, x, F in cv:
                for p, q in zip(x, F):
                    f.write('%d,%.4f,%.1f\n' % (s, p, q))
    for s, x, F in cv:
        print('%3d MPa: x_f = %.2f mm, F_max = %.0f N, F at fracture = %.0f N' % (s, x[-1], F.max(), F[-1]))
