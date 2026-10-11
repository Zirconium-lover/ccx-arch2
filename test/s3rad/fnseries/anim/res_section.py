"""Looping GIF from the computed model: a section of the block (load along
x, the T-R plane, at the height that cuts the most plates) while the grip
strain grows.  The matrix is coloured by omega (PAL picks the scale), the
hydrides by their orientation with a dark outline, deleted elements leave
the crack black.  Right: the stress-strain curve of the run with the current
point.  Same axes and scale for every run, so the GIFs compare side by side.

    RUN=fn40 PAL=white manim -r 1600,900 --fps 25 --disable_caching res_section.py ResSection
"""
import os
from mech_common import *
from fndata import Run, SERIES

RUN = os.environ.get('RUN', 'fn40')
T_LOOP = 12.0
T_RUN = (0.6, 8.8)                    # the grip strain grows
T_DIP = (9.8, 10.5, 10.8, 11.5)
NB = 24                               # colour levels of omega
EPS_MAX, SIG_MAX = 0.135, 650.        # common axes of the series
U_MAX = 0.55                          # largest grip displacement of the series, mm
STEEL = ManimColor('#34404F')         # the metal of the mechanism GIFs
C = ManimColor
# How omega is drawn on the matrix.  kind: smooth (a colour scale), bands
# (steps), lines (isolines over a flat matrix), bands_lines, dots (white
# specks whose density grows with omega, like the pore cloud of the ductile
# mechanism GIF).  'f' maps omega to what is drawn (dots: probability).
GREY = [(0.0, C('#2E2E2E')), (0.5, C('#5C5C5C')), (1.0, C('#F2F2F2'))]
STYLES = {
    # the first rounds
    'white': dict(kind='smooth', stops=[(0.0, STEEL), (0.45, C('#3E4A5A')), (0.8, C('#8D99A8')), (1.0, C('#FFFFFF'))]),
    'violet': dict(kind='smooth', stops=[(0.0, STEEL), (0.45, C('#403F63')), (0.8, C('#9156C4')), (1.0, C('#F2D2FF'))]),
    'rainbow': dict(kind='smooth', stops=[(0.0, C('#30123B')), (0.13, C('#4662D7')), (0.27, C('#36AAF9')),
                                          (0.4, C('#1AE4B6')), (0.53, C('#72FE5E')), (0.66, C('#C8EF34')),
                                          (0.78, C('#FABA39')), (0.89, C('#F66B19')), (1.0, C('#C42503'))]),
    'coolwarm': dict(kind='smooth', stops=[(0.0, C('#3B4CC0')), (0.5, C('#DDDDDD')), (1.0, C('#B40426'))]),
    'teal_pink': dict(kind='smooth', stops=[(0.0, C('#0B3C49')), (1.0, C('#F0609E'))]),
    'viridis': dict(kind='smooth', stops=[(0.0, C('#440154')), (0.25, C('#3B528B')), (0.5, C('#21918C')),
                                          (0.75, C('#5EC962')), (1.0, C('#FDE725'))]),
    'magma': dict(kind='smooth', stops=[(0.0, C('#1C1044')), (0.35, C('#51127C')), (0.6, C('#B73779')),
                                        (0.82, C('#FC8961')), (1.0, C('#FCFDBF'))]),
    # neutral: the hydrides keep the only colours
    'grey': dict(kind='smooth', stops=GREY),
    'grey_bands': dict(kind='bands', stops=GREY, levels=[0, 0.2, 0.4, 0.6, 0.8, 0.95, 1.0001]),
    'isolines': dict(kind='lines', base=C('#2B2F35'), levels=[0.2, 0.4, 0.6, 0.8], top=0.97),
    'bands_lines': dict(kind='bands_lines', stops=[(0.0, C('#26292E')), (1.0, C('#9A9DA3'))],
                        levels=[0, 0.2, 0.4, 0.6, 0.8, 1.0001]),
    'dots': dict(kind='dots', base=C('#2B2F35'), f=lambda w: smoothstep((w - 0.35) / 0.6)),
    'haze': dict(kind='smooth', stops=[(0.0, C('#2B2F35')), (0.6, C('#2F3339')), (0.85, C('#8E9196')),
                                       (1.0, C('#FFFFFF'))]),
    'sand': dict(kind='smooth', stops=[(0.0, C('#2A2622')), (0.55, C('#5E554B')), (0.85, C('#B9AC98')),
                                       (1.0, C('#F6F0E6'))]),
    'mint': dict(kind='smooth', stops=[(0.0, C('#1D2624')), (0.55, C('#2C5A50')), (0.85, C('#58B79C')),
                                       (1.0, C('#D9F7EC'))]),
}
# steps + isolines in several colour schemes: (low, high) or full stops, line colour
BL = [0, 0.2, 0.4, 0.6, 0.8, 1.0001]
for name, stops, line in (
        ('bl_grey', [(0.0, C('#26292E')), (1.0, C('#9A9DA3'))], '#E6E6E6'),
        ('bl_contrast', [(0.0, C('#1C1C1C')), (1.0, C('#E4E4E4'))], '#000000'),
        ('bl_sand', [(0.0, C('#2A2622')), (1.0, C('#CDBFA8'))], '#F6F0E6'),
        ('bl_mint', [(0.0, C('#1D2624')), (1.0, C('#7FC7B0'))], '#E3FFF5'),
        ('bl_violet', [(0.0, C('#24203A')), (1.0, C('#A58BD6'))], '#EFE7FF'),
        ('bl_steel', [(0.0, C('#1F2833')), (1.0, C('#8FA3BA'))], '#EAF2FA'),
        ('bl_olive', [(0.0, C('#262A1F')), (1.0, C('#A9B07A'))], '#F3F6E2'),
        ('bl_rainbow', STYLES['rainbow']['stops'], '#111111'),
        ('bl_viridis', STYLES['viridis']['stops'], '#111111'),
        ('bl_magma', STYLES['magma']['stops'], '#F2F2F2'),
        ('bl_coolwarm', STYLES['coolwarm']['stops'], '#2A2A2A'),
        ('bl_teal_pink', STYLES['teal_pink']['stops'], '#F2F2F2')):
    STYLES[name] = dict(kind='bands_lines', stops=stops, levels=BL, line=line)
# rainbow laid semi-transparent over a dark grey matrix: muted colours that
# leave the red / blue hydrides readable (steps or smooth, two opacities)
# (mixed into the colours beforehand: semi-transparent triangles would show
# their seams)
for name, kind, al in (('rb_bands_50', 'bands', 0.5), ('rb_bands_70', 'bands', 0.7),
                       ('rb_smooth_50', 'smooth', 0.5), ('rb_smooth_70', 'smooth', 0.7)):
    STYLES[name] = dict(kind=kind, levels=BL,
                        stops=[(a, interpolate_color(C('#2B2F35'), c, al)) for a, c in STYLES['rainbow']['stops']])
# glass: the rainbow mixed into a light grey (tinted frosted glass instead of a
# darkened field); sheen adds a soft light from the upper left and a faint
# diagonal streak, as on a glass pane
for name, base, frac, sheen in (('gl_light', '#C9CED6', 0.55, 0.), ('gl_light_sheen', '#C9CED6', 0.55, 0.22),
                                ('gl_mid_sheen', '#7D848E', 0.6, 0.22), ('gl_vivid_sheen', '#FFFFFF', 0.72, 0.18)):
    STYLES[name] = dict(kind='smooth', sheen=sheen,
                        stops=[(a, interpolate_color(C(base), c, frac)) for a, c in STYLES['rainbow']['stops']])
STYLE = os.environ.get('PAL', 'white')
ST = STYLES[STYLE]
DOT_DENS = 1300.                      # specks per mm^2 of section at full omega


def fn_tex(k):
    v = '0' if k == 0 else '1' if k == 12 else ('%.2f' % (k / 12)).replace('.', '{,}')
    return r'F_n = ' + v


PX = 1600 / config.frame_width        # pixels per scene unit at 1600 x 900


class FieldRaster:
    """Picture of a nodal field on the section polygons in one of the STYLES,
    rendered by matplotlib into an RGBA array of the given box (scene
    units), transparent where there is no material."""
    def __init__(self, box, style):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap, ListedColormap
        self.box, self.st = box, style
        self.w, self.h = int(round(box.w * PX)), int(round(box.h * PX))
        self.fig = plt.figure(figsize=(self.w / 100, self.h / 100), dpi=100)
        self.fig.patch.set_alpha(0)
        self.ax = self.fig.add_axes([0, 0, 1, 1])
        if 'stops' in style:
            self.cmap = LinearSegmentedColormap.from_list('w', [(a, c.to_hex()) for a, c in style['stops']])
        if 'base' in style:
            self.flat = ListedColormap([style['base'].to_hex()])

    @staticmethod
    def mesh(polys, vals):
        """One connected triangulation (shared corners merged), so that
        isolines and bands run on without seams."""
        import matplotlib.tri as mtri
        pts, tris, vv, n = [], [], [], 0
        for q, v in zip(polys, vals):
            pts.append(q); vv.append(v)
            tris.append([n, n + 1, n + 2])
            if len(q) == 4:
                tris.append([n, n + 2, n + 3])
            n += len(q)
        pts = np.vstack(pts); vv = np.concatenate(vv)
        _, idx, inv = np.unique(np.round(pts / 1e-4).astype(np.int64), axis=0, return_index=True,
                                return_inverse=True)
        T = inv.ravel()[np.array(tris)]
        P = pts[idx]
        # merging can collapse a sliver; contouring loops forever on those
        A = (P[T[:, 1], 0] - P[T[:, 0], 0]) * (P[T[:, 2], 1] - P[T[:, 0], 1]) - \
            (P[T[:, 2], 0] - P[T[:, 0], 0]) * (P[T[:, 1], 1] - P[T[:, 0], 1])
        ok = (T[:, 0] != T[:, 1]) & (T[:, 1] != T[:, 2]) & (T[:, 0] != T[:, 2]) & (np.abs(A) > 1e-12)
        return mtri.Triangulation(P[:, 0], P[:, 1], T[ok]), np.clip(vv[idx], 0, 1)

    def render(self, polys, vals, dots=None):
        """dots: (xy (m,2), omega (m,), u (m,)) - a speck shows when f(omega) > u."""
        st, ax = self.st, self.ax
        ax.clear(); ax.set_axis_off(); ax.patch.set_alpha(0)
        ax.set_xlim(self.box.x0, self.box.x1); ax.set_ylim(self.box.y0, self.box.y1)
        if polys:
            tr, v = self.mesh(polys, vals)
            kind = st['kind']
            al = st.get('alpha', 1.)                       # < 1: over the flat base colour
            if al < 1:
                ax.tripcolor(tr, v * 0, shading='gouraud', cmap=self.flat, vmin=0, vmax=1)
            if kind == 'smooth':
                ax.tripcolor(tr, v, shading='gouraud', cmap=self.cmap, vmin=0, vmax=1, alpha=al)
            elif kind in ('bands', 'bands_lines'):
                lv = st['levels']
                mids = [(a + b) / 2 for a, b in zip(lv[:-1], lv[1:])]
                cols = [self.cmap(m) for m in mids]
                if al == 1:
                    ax.tripcolor(tr, v * 0, shading='gouraud', cmap=self.cmap, vmin=0, vmax=1)  # under the seams
                ax.tricontourf(tr, v, levels=lv, colors=cols, antialiased=True, alpha=al)
                if kind == 'bands_lines':
                    ax.tricontour(tr, v, levels=lv[1:-1], colors=st.get('line', '#E6E6E6'), linewidths=0.9,
                                  alpha=0.8)
            elif kind in ('lines', 'dots'):
                ax.tripcolor(tr, v * 0, shading='gouraud', cmap=self.flat, vmin=0, vmax=1)
                if kind == 'lines':
                    ax.tricontour(tr, v, levels=st['levels'], colors='#D8D8D8', linewidths=1.1, alpha=0.75)
                    ax.tricontour(tr, v, levels=[st['top']], colors='#FFFFFF', linewidths=2.6)
        if dots is not None and len(dots[0]):
            xy, w, u = dots
            on = st['f'](w) > u
            ax.scatter(xy[on, 0], xy[on, 1], s=3.2, c='#F4F4F4', linewidths=0)
        self.fig.canvas.draw()
        img = np.asarray(self.fig.canvas.buffer_rgba()).copy()
        sh = st.get('sheen', 0.)
        if sh > 0:                                     # screen a light gradient over the material
            hh, ww = img.shape[:2]
            u = np.linspace(0, 1, ww)[None, :]; v = np.linspace(0, 1, hh)[:, None]
            d = (u + v) / 2                            # 0 at the upper left
            g = sh * (1 - d) ** 1.5 + 0.5 * sh * np.exp(-((d - 0.32) / 0.07) ** 2)
            rgb = img[..., :3].astype(float)
            img[..., :3] = (rgb + (255 - rgb) * np.clip(g, 0, 1)[..., None]).astype(np.uint8)
        return img


def sample_dots(r, z, rng):
    """Fixed specks in the material of the section: owner tet and its
    barycentric weights, a random threshold each."""
    polys, own, _ = r.slice('z', z, None, None, r.phase == 0)
    tets, W = [], []
    for q, e in zip(polys, own):
        tri = [(q[0], q[i], q[i + 1]) for i in range(1, len(q) - 1)]
        for a, b, c in tri:
            A = 0.5 * abs((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]))
            n = rng.poisson(A * DOT_DENS)
            if not n:
                continue
            s_, t_ = rng.random(n), rng.random(n)
            fl = s_ + t_ > 1
            s_[fl], t_[fl] = 1 - s_[fl], 1 - t_[fl]
            p = a + s_[:, None] * (b - a) + t_[:, None] * (c - a)
            X = r.X0[r.T[e]]                                   # barycentric in the tet
            M = np.vstack([X.T, np.ones(4)])
            lam = np.linalg.solve(M, np.vstack([p.T, np.ones(n)])).T
            tets.append(np.full(n, e)); W.append(lam)
    tets, W = np.concatenate(tets), np.vstack(W)
    return tets, W, rng.random(len(tets))


class ResSection(Scene):
    def construct(self):
        self.camera.background_color = BG
        r = Run(RUN)
        k = r.k_radial
        Z_CUT = float(os.environ.get('ZCUT', 0)) or r.best_section('z')
        print('section z = %.2f' % Z_CUT)
        e_end = float(r.curve[-1, 0])
        e1 = 0.9 * r.eps_peak

        cellL = Box(CONTENT.x0, CONTENT.x0 + 7.7, CONTENT.y0, CONTENT.y1)
        cellR = Box(cellL.x1 + 0.6, CONTENT.x1, CONTENT.y0, CONTENT.y1)

        # ---------------- left: Fn on top, the section below ----------------
        head = tex(fn_tex(k), L1)
        AL, AG = 0.7, 0.12                               # arrow length, gap (scene)
        LX, LY = r.hi[0] - r.lo[0], r.hi[1] - r.lo[1]
        secbox = Box(cellL.x0, cellL.x1, cellL.y0, cellL.y1 - head.height - 0.3)
        sc = min((secbox.w - 2 * (AL + AG)) / (LX + U_MAX), secbox.h / LY)
        x_left = secbox.x0 + AL + AG
        y_mid = secbox.c[1]
        place(head, Box(x_left, cellL.x1, secbox.y1 + 0.3, cellL.y1), UL)
        field_box = Box(x_left - 0.1, x_left + (LX + U_MAX) * sc + 0.1, y_mid - LY * sc / 2 - 0.1,
                        y_mid + LY * sc / 2 + 0.1)
        raster = FieldRaster(field_box, ST)
        dots = sample_dots(r, Z_CUT, np.random.default_rng(7)) if ST['kind'] == 'dots' else None

        def P(xy):
            xy = np.asarray(xy)
            return np.c_[x_left + (xy[..., 0] - r.lo[0]) * sc, y_mid + (xy[..., 1] - (r.lo[1] + r.hi[1]) / 2) * sc,
                         np.zeros(xy.shape[:-1])].reshape(*xy.shape[:-1], 3)

        def eps_of(t):
            if t <= T_RUN[0]:
                return 0.
            tau = smoothstep((min(t, T_RUN[1]) - T_RUN[0]) / (T_RUN[1] - T_RUN[0]))
            return float(np.interp(tau, [0, 0.3, 1], [0, e1, e_end]))

        clk = ValueTracker(0.)

        def dip(t):
            if t < T_DIP[0]:
                return 1.
            return 1 - ramp(t, T_DIP[0], T_DIP[1]) + ramp(t, T_DIP[2], T_DIP[3])

        def te_of(t):
            return t if t < T_DIP[1] + 0.15 else 0.

        def polys_vm(polys, color, sw=1.2):
            """Many polygons as one path; a stroke of the fill colour hides
            the seams between neighbours."""
            pts = []
            for q in polys:
                b = np.roll(q, -1, axis=0)
                pts.append(np.stack([q, q + (b - q) / 3, q + 2 * (b - q) / 3, b], 1).reshape(-1, 3))
            vm = VMobject(stroke_color=color, stroke_width=sw, fill_color=color, fill_opacity=1)
            if pts:
                vm.set_points(np.vstack(pts))
            return vm

        def picture():
            t = clk.get_value()
            g = dip(t)
            eps = eps_of(te_of(t))
            s = r.at(eps)
            polys, own, pv = r.slice('z', Z_CUT, s, 'omega')
            out = Group()
            ph = r.phase[own] if len(own) else np.zeros(0, int)
            mat = np.flatnonzero(ph == 0)
            dd = None
            if dots is not None:
                te_, W, u = dots
                al = s.alive[te_]
                nodes = r.T[te_[al]]
                xy = np.einsum('mk,mkj->mj', W[al], s.X[nodes])[:, :2]
                w = np.einsum('mk,mk->m', W[al], s.omega[nodes])
                dd = (P(xy)[:, :2], w, u[al])
            img = raster.render([P(polys[i][:, :2])[:, :2] for i in mat], [pv[i] for i in mat], dd)
            if g < 1:
                img[..., 3] = (img[..., 3] * g).astype(np.uint8)
            im = ImageMobject(img)
            im.stretch_to_fit_width(field_box.w).stretch_to_fit_height(field_box.h).move_to(field_box.c)
            out.add(im)
            vec = VGroup()
            for code, col in ((1, RAD), (2, TAN)):
                sel = np.flatnonzero(ph == code)
                if len(sel):
                    vec.add(polys_vm([P(polys[i][:, :2]) for i in sel], col))
            # one dark outline around each plate's cut (its outer edges only)
            hp = [P(polys[i][:, :2]) for i in np.flatnonzero(ph > 0)]
            cnt = {}
            for q in hp:
                for a, b in zip(q, np.roll(q, -1, axis=0)):
                    ka, kb = tuple(np.round(a[:2], 4)), tuple(np.round(b[:2], 4))
                    key = (ka, kb) if ka < kb else (kb, ka)
                    cnt[key] = (cnt[key][0] + 1, a, b) if key in cnt else (1, a, b)
            seg = [np.array([a, a + (b - a) / 3, a + 2 * (b - a) / 3, b]) for n_, a, b in cnt.values() if n_ == 1]
            if seg:
                vec.add(VMobject(stroke_color=BLACK, stroke_width=2.5).set_points(np.vstack(seg)))
            # load arrows on the current faces
            xr = r.hi[0] + r.grip * r.t_of(eps)
            for yy in np.linspace(0.2, 0.8, 3) * LY + r.lo[1]:
                for x0, sg in ((r.lo[0], -1), (xr, 1)):
                    p0 = P(np.array([[x0, yy]]))[0] + sg * RIGHT * AG
                    vec.add(Arrow(p0, p0 + sg * RIGHT * AL, buff=0, color=WHITE, stroke_width=6,
                                  max_tip_length_to_length_ratio=0.3))
            out.add(fade_all(vec, g) if g < 1 else vec)
            return out

        pic = Group()
        pic.add_updater(lambda m: m.set_submobjects(list(picture())))

        # ---------------- right: curve, omega bar ----------------
        lab_s, lab_e = tex(r'\sigma', L1), tex(r'\varepsilon', L1)
        BW_, BH_ = 3.0, 0.34
        lab_w = tex(r'\omega', L1)
        bar = Rectangle(width=BW_, height=BH_, stroke_color=GREY_B, stroke_width=2)
        z0, z1 = tex('0', L3M, GREY_A), tex('1', L3M, GREY_A)
        wrow = VGroup(lab_w, bar).arrange(RIGHT, buff=0.3)
        z0.next_to(bar, DOWN, buff=0.12).align_to(bar, LEFT)
        z1.next_to(bar, DOWN, buff=0.12).align_to(bar, RIGHT)
        wgrp = Group(wrow, z0, z1)
        hg = cellR.h - wgrp.height - 0.6
        rows = cellR.stack(hg, wgrp.height, gap=0.6)
        gbox = rows[0]
        ax = Axes(x_range=[0, EPS_MAX, EPS_MAX], y_range=[0, SIG_MAX, SIG_MAX],
                  x_length=gbox.w - lab_s.width - 0.3 - lab_e.width - 0.3,
                  y_length=gbox.h - 0.12, tips=True,
                  axis_config=dict(color=GREY_B, stroke_width=3, include_ticks=False, tip_width=0.2, tip_height=0.2))
        ax.move_to(np.array([gbox.x0 + lab_s.width + 0.3 + ax.width / 2, gbox.y1 - 0.06 - ax.height / 2, 0.]))
        lab_s.next_to(ax.y_axis, LEFT, buff=0.2).align_to(ax.y_axis, UP)
        lab_e.next_to(ax.x_axis, RIGHT, buff=0.15)
        place(wgrp, rows[1], LEFT)
        bb = Box(bar.get_left()[0], bar.get_right()[0], bar.get_bottom()[1], bar.get_top()[1])
        braster = FieldRaster(bb, ST)
        nq = 60
        xs = np.linspace(bb.x0, bb.x1, nq + 1)
        qs = [np.array([[xs[i], bb.y0], [xs[i + 1], bb.y0], [xs[i + 1], bb.y1], [xs[i], bb.y1]]) for i in range(nq)]
        qv = [np.array([i, i + 1, i + 1, i]) / nq for i in range(nq)]
        bd = None
        if ST['kind'] == 'dots':
            rng = np.random.default_rng(3)
            m = rng.poisson(bb.w * bb.h / sc ** 2 * DOT_DENS)
            bxy = np.c_[bb.x0 + bb.w * rng.random(m), bb.y0 + bb.h * rng.random(m)]
            bd = (bxy, (bxy[:, 0] - bb.x0) / bb.w, rng.random(m))
        bimg = ImageMobject(braster.render(qs, qv, bd))
        bimg.stretch_to_fit_width(bb.w).stretch_to_fit_height(bb.h).move_to(bb.c)
        wgrp.add_to_back(bimg)
        cpts = [ax.c2p(0, 0)] + [ax.c2p(e, s_) for e, s_ in r.curve]
        ghost = VMobject(stroke_color=GREY_D, stroke_width=3).set_points_as_corners(cpts)
        cx = np.r_[0, r.curve[:, 0]]; cy = np.r_[0, r.curve[:, 1]]

        def traced():
            t = clk.get_value()
            eps = eps_of(te_of(t))
            s_now = float(np.interp(eps, cx, cy))
            g = VGroup()
            if eps > 1e-4:
                m = cx < eps
                pts = [ax.c2p(e, s_) for e, s_ in zip(cx[m], cy[m])] + [ax.c2p(eps, s_now)]
                g.add(VMobject(stroke_color=WHITE, stroke_width=6).set_points_as_corners(pts))
            g.add(Dot(ax.c2p(eps, s_now), radius=0.11, color=WHITE))
            return fade_all(g, dip(t)) if dip(t) < 1 else g
        tr = VGroup()
        tr.add_updater(lambda m: m.set_submobjects(list(traced())))
        sec_rect = Rectangle(width=(LX + U_MAX) * sc + 2 * (AL + AG), height=LY * sc).move_to(
            np.array([secbox.x0 + ((LX + U_MAX) * sc + 2 * (AL + AG)) / 2, y_mid, 0.]))
        report({'section': sec_rect, 'head': head, 'graph': VGroup(ax, lab_s, lab_e), 'omega': wgrp})

        self.add(pic, ghost, ax, lab_s, lab_e, tr, head, wgrp)
        if os.environ.get('STILL'):
            clk.set_value(float(os.environ['STILL']))
            self.update_mobjects(0)
            self.wait(0.04)
            return
        self.play(clk.animate.set_value(T_LOOP), run_time=T_LOOP, rate_func=linear)
