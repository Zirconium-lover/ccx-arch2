"""Looping GIF from the computed model: a section through the middle of the
block (load along x, the T-R plane) while the grip strain grows.  The matrix
is coloured by omega (the colours of the mechanism GIFs), the hydrides by
their orientation, debonded interfaces glow green, deleted elements leave
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
Z_CUT = 0.8                           # the middle of the block (z from 0 to 1.6)
STEEL = ManimColor('#34404F')         # the metal of the mechanism GIFs
# omega: the matrix stays steel while omega is low, only the band near
# fracture lights up (no yellow fields)
PALETTES = {
    'white': ([(0.0, STEEL), (0.45, ManimColor('#3E4A5A')), (0.8, ManimColor('#8D99A8')),
               (1.0, ManimColor('#FFFFFF'))], WHITE),
    'violet': ([(0.0, STEEL), (0.45, ManimColor('#403F63')), (0.8, ManimColor('#9156C4')),
                (1.0, ManimColor('#F2D2FF'))], ManimColor('#C38BEA')),
}
PAL = os.environ.get('PAL', 'white')
STOPS, OM_LABEL = PALETTES[PAL]


def field_color(w):
    w = float(np.clip(w, 0, 1))
    for (a, ca), (b, cb) in zip(STOPS[:-1], STOPS[1:]):
        if w <= b:
            return interpolate_color(ca, cb, (w - a) / (b - a))
    return STOPS[-1][1]


def field_bar(width, height=0.3):
    n = 48
    g = VGroup(*[Rectangle(width=width / n + 0.005, height=height, stroke_width=0,
                           fill_color=field_color((i + 0.5) / n), fill_opacity=1)
                 .move_to(np.array([(i + 0.5) * width / n - width / 2, 0., 0.])) for i in range(n)])
    return VGroup(g, Rectangle(width=width, height=height, stroke_color=GREY_B, stroke_width=2))


def fn_tex(k):
    v = '0' if k == 0 else '1' if k == 12 else r'{%d \over 12}' % k
    return r'F_n = ' + v


PX = 1600 / config.frame_width        # pixels per scene unit at 1600 x 900


class FieldRaster:
    """Smooth (Gouraud) picture of a nodal field on the section polygons,
    rendered by matplotlib into an RGBA array of the given box (scene units),
    transparent where there is no material."""
    def __init__(self, box):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap
        self.plt, self.box = plt, box
        self.w, self.h = int(round(box.w * PX)), int(round(box.h * PX))
        self.fig = plt.figure(figsize=(self.w / 100, self.h / 100), dpi=100)
        self.fig.patch.set_alpha(0)
        self.ax = self.fig.add_axes([0, 0, 1, 1])
        self.ax.set_axis_off(); self.ax.patch.set_alpha(0)
        self.cmap = LinearSegmentedColormap.from_list(
            'omega', [(a, c.to_hex()) for a, c in STOPS])

    def render(self, polys, vals):
        import matplotlib.tri as mtri
        ax = self.ax
        ax.clear(); ax.set_axis_off(); ax.patch.set_alpha(0)
        ax.set_xlim(self.box.x0, self.box.x1); ax.set_ylim(self.box.y0, self.box.y1)
        if polys:
            pts, tris, vv = [], [], []
            n = 0
            for q, v in zip(polys, vals):
                pts.append(q); vv.append(v)
                tris.append([n, n + 1, n + 2])
                if len(q) == 4:
                    tris.append([n, n + 2, n + 3])
                n += len(q)
            pts = np.vstack(pts); vv = np.concatenate(vv)
            tr = mtri.Triangulation(pts[:, 0], pts[:, 1], np.array(tris))
            ax.tripcolor(tr, np.clip(vv, 0, 1), shading='gouraud', cmap=self.cmap, vmin=0, vmax=1)
        self.fig.canvas.draw()
        return np.asarray(self.fig.canvas.buffer_rgba()).copy()


class ResSection(Scene):
    def construct(self):
        self.camera.background_color = BG
        r = Run(RUN)
        k = r.k_radial
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
        raster = FieldRaster(field_box)

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
            img = raster.render([P(polys[i][:, :2])[:, :2] for i in mat], [pv[i] for i in mat])
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
            # debonded interface: green segments, brighter as the traction drops
            segs, fi = r.slice_facets('z', Z_CUT, s)
            lam = s.debond[fi]
            for lo_, hi_, op in ((0.5, 0.95, 0.45), (0.95, 1.01, 1.0)):
                sel = np.flatnonzero((lam >= lo_) & (lam < hi_))
                if len(sel):
                    pts = []
                    for i in sel:
                        a, b = P(segs[i][:, :2])
                        pts.append(np.array([a, a + (b - a) / 3, a + 2 * (b - a) / 3, b]))
                    vec.add(VMobject(stroke_color=IFC, stroke_width=5, stroke_opacity=op).set_points(np.vstack(pts)))
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

        # ---------------- right: curve, omega bar, debond key ----------------
        lab_s, lab_e = tex(r'\sigma', L1), tex(r'\varepsilon', L1)
        bar = field_bar(3.0, 0.3)
        lab_w = tex(r'\omega', L1, OM_LABEL)
        z0, z1 = tex('0', L3M, GREY_A), tex('1', L3M, GREY_A)
        wrow = VGroup(lab_w, bar).arrange(RIGHT, buff=0.3)
        z0.next_to(bar, DOWN, buff=0.12).align_to(bar, LEFT)
        z1.next_to(bar, DOWN, buff=0.12).align_to(bar, RIGHT)
        wgrp = VGroup(wrow, z0, z1)
        sw = Line(ORIGIN, RIGHT * lab_w.width, stroke_width=6, color=IFC)
        deb = VGroup(sw, txt('отслоение', L2, IFC)).arrange(RIGHT, buff=0.3)
        hg = cellR.h - wgrp.height - deb.height - 2 * 0.5
        rows = cellR.stack(hg, wgrp.height, deb.height, gap=0.5)
        gbox = rows[0]
        ax = Axes(x_range=[0, EPS_MAX, EPS_MAX], y_range=[0, SIG_MAX, SIG_MAX],
                  x_length=gbox.w - lab_s.width - 0.3 - lab_e.width - 0.3,
                  y_length=gbox.h - 0.12, tips=True,
                  axis_config=dict(color=GREY_B, stroke_width=3, include_ticks=False, tip_width=0.2, tip_height=0.2))
        ax.move_to(np.array([gbox.x0 + lab_s.width + 0.3 + ax.width / 2, gbox.y1 - 0.06 - ax.height / 2, 0.]))
        lab_s.next_to(ax.y_axis, LEFT, buff=0.2).align_to(ax.y_axis, UP)
        lab_e.next_to(ax.x_axis, RIGHT, buff=0.15)
        place(wgrp, rows[1], LEFT)
        place(deb, rows[2], LEFT)
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
        report({'section': sec_rect, 'head': head, 'graph': VGroup(ax, lab_s, lab_e), 'omega': wgrp, 'debond': deb})

        self.add(pic, ghost, ax, lab_s, lab_e, tr, head, wgrp, deb)
        if os.environ.get('STILL'):
            clk.set_value(float(os.environ['STILL']))
            self.update_mobjects(0)
            self.wait(0.04)
            return
        self.play(clk.animate.set_value(T_LOOP), run_time=T_LOOP, rate_func=linear)
