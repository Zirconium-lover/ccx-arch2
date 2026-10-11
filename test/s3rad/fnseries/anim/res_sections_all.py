"""Looping GIF from the computed model: the sections of all eight runs side
by side (one plane, z = 0.61 mm, the one that cuts the most plates; one
scale), stretched together by a common grip strain; right: their
stress-strain curves running with it.  omega on the matrix as grey steps
with isolines (style bl_grey of res_section.py); hydrides in their colours;
deleted elements leave the cracks black.  Each section's Fn is written in
the colour of its curve.

    manim -r 1600,900 --fps 25 --disable_caching res_sections_all.py ResSectionsAll
"""
import os
from mech_common import *
from fndata import Run, SERIES
from res_section import FieldRaster, STYLES, U_MAX
from res_curves import fn_color, SIG_MAX

T_LOOP = 13.5
T_RUN = (0.6, 10.2)                   # the common grip strain grows
T_DIP = (11.3, 12.0, 12.3, 13.0)
EPS_END = 0.131                       # past the last fracture (Fn = 0, eps_f10 12.9 %)
EPS_AX = 0.135
ST = STYLES[os.environ.get('PAL', 'bl_grey')]


def frac_tex(k):
    v = '0' if k == 0 else '1' if k == 12 else '%d/12' % k
    return r'F_n = ' + v


class ResSectionsAll(Scene):
    def construct(self):
        self.camera.background_color = BG
        runs = [Run(fn) for fn, _ in SERIES]
        F = [r.k_radial / 12 for r in runs]
        z_cut = float(os.environ.get('ZCUT', 0)) or runs[0].best_section('z')
        r0 = runs[0]
        LX, LY = r0.hi[0] - r0.lo[0], r0.hi[1] - r0.lo[1]

        # ---------------- layout: a 3 x 3 grid of sections | the curves ----------------
        labs = [tex(frac_tex(r.k_radial), L3M, fn_color(f)) for r, f in zip(runs, F)]
        hl = max(m.height for m in labs)
        GX, GY, GL = 0.3, 0.28, 0.1                     # gaps: columns, rows, label - section
        ch = (CONTENT.h - 2 * GY) / 3 - hl - GL           # section height
        sc = ch / LY
        cw = (LX + U_MAX) * sc                            # cell width (room for the stretch)
        grid_w = 3 * cw + 2 * GX
        cellR = Box(CONTENT.x0 + grid_w + 0.7, CONTENT.x1, CONTENT.y0, CONTENT.y1)
        origins = []                                      # lower-left of each section (scene)
        for i in range(9):
            c, rw = i % 3, i // 3
            x0 = CONTENT.x0 + c * (cw + GX)
            ytop = CONTENT.y1 - rw * (hl + GL + ch + GY)
            origins.append(np.array([x0, ytop - hl - GL - ch, 0.]))
        for i, lab in enumerate(labs):
            o = origins[i]
            lab.move_to(np.array([o[0], o[1] + ch + GL, 0.]), aligned_edge=DL)

        def Pk(i, xy):
            xy = np.asarray(xy)
            o = origins[i]
            return np.c_[o[0] + (xy[..., 0] - r0.lo[0]) * sc, o[1] + (xy[..., 1] - r0.lo[1]) * sc,
                         np.zeros(xy.shape[:-1])].reshape(*xy.shape[:-1], 3)

        grid_box = Box(CONTENT.x0 - 0.05, CONTENT.x0 + grid_w + 0.05, CONTENT.y0 - 0.05, CONTENT.y1 + 0.05)
        raster = FieldRaster(grid_box, ST)

        # the omega key in the ninth cell
        o9 = origins[8]
        key_box = Box(o9[0] + 0.15, o9[0] + cw - 0.1, o9[1], o9[1] + ch + hl + GL)
        lab_w = tex(r'\omega', L2M)
        bw = key_box.w - lab_w.width - 0.25
        bar = Rectangle(width=bw, height=0.3, stroke_color=GREY_B, stroke_width=2)
        z0, z1 = tex('0', L3M, GREY_A), tex('1', L3M, GREY_A)
        wrow = VGroup(lab_w, bar).arrange(RIGHT, buff=0.25)
        z0.next_to(bar, DOWN, buff=0.1).align_to(bar, LEFT)
        z1.next_to(bar, DOWN, buff=0.1).align_to(bar, RIGHT)
        key = Group(wrow, z0, z1)
        place(key, key_box)
        bb = Box(bar.get_left()[0], bar.get_right()[0], bar.get_bottom()[1], bar.get_top()[1])
        nq = 60
        xs = np.linspace(bb.x0, bb.x1, nq + 1)
        qs = [np.array([[xs[i], bb.y0], [xs[i + 1], bb.y0], [xs[i + 1], bb.y1], [xs[i], bb.y1]]) for i in range(nq)]
        qv = [np.array([i, i + 1, i + 1, i]) / nq for i in range(nq)]
        bimg = ImageMobject(FieldRaster(bb, ST).render(qs, qv))
        bimg.stretch_to_fit_width(bb.w).stretch_to_fit_height(bb.h).move_to(bb.c)
        key.add_to_back(bimg)

        # ---------------- the curves ----------------
        ls, le = tex(r'\sigma', L1), tex(r'\varepsilon', L1)
        ax = Axes(x_range=[0, EPS_AX, EPS_AX], y_range=[0, SIG_MAX, SIG_MAX],
                  x_length=cellR.w - ls.width - 0.3 - le.width - 0.3, y_length=cellR.h - 0.2, tips=True,
                  axis_config=dict(color=GREY_B, stroke_width=3, include_ticks=False, tip_width=0.2, tip_height=0.2))
        ax.move_to(np.array([cellR.x0 + ls.width + 0.3 + ax.width / 2, cellR.y1 - 0.06 - ax.height / 2, 0.]))
        ls.next_to(ax.y_axis, LEFT, buff=0.2).align_to(ax.y_axis, UP)
        le.next_to(ax.x_axis, RIGHT, buff=0.15)
        cx = [np.r_[0, r.curve[:, 0]] for r in runs]
        cy = [np.r_[0, r.curve[:, 1]] for r in runs]

        clk = ValueTracker(0.)

        def dip(t):
            if t < T_DIP[0]:
                return 1.
            return 1 - ramp(t, T_DIP[0], T_DIP[1]) + ramp(t, T_DIP[2], T_DIP[3])

        def te_of(t):
            return t if t < T_DIP[1] + 0.15 else 0.

        def eps_of(te):
            if te <= T_RUN[0]:
                return 0.
            tau = smoothstep((min(te, T_RUN[1]) - T_RUN[0]) / (T_RUN[1] - T_RUN[0]))
            return float(np.interp(tau, [0, 0.15, 1], [0, 0.03, EPS_END]))  # little time on the uniform part

        def polys_vm(polys, color, sw=0.5):
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
            out = Group()
            allp, allv, vec = [], [], VGroup()
            for i, r in enumerate(runs):
                s = r.at(eps)
                polys, own, pv = r.slice('z', z_cut, s, 'omega')
                ph = r.phase[own]
                for j in np.flatnonzero(ph == 0):
                    allp.append(Pk(i, polys[j][:, :2])[:, :2]); allv.append(pv[j])
                for code, col in ((1, RAD), (2, TAN)):
                    sel = np.flatnonzero(ph == code)
                    if len(sel):
                        vec.add(polys_vm([Pk(i, polys[j][:, :2]) for j in sel], col))
            img = raster.render(allp, allv)
            if g < 1:
                img[..., 3] = (img[..., 3] * g).astype(np.uint8)
            im = ImageMobject(img)
            im.stretch_to_fit_width(grid_box.w).stretch_to_fit_height(grid_box.h).move_to(grid_box.c)
            out.add(im)
            # the curves, Fn = 1 on top; a dot on each while it runs
            for x, y, f in zip(cx, cy, F):
                e = min(eps, x[-1])
                if e > 1e-4:
                    m = x < e
                    pts = [ax.c2p(a, b) for a, b in zip(x[m], y[m])] + [ax.c2p(e, np.interp(e, x, y))]
                    vec.add(VMobject(stroke_color=fn_color(f), stroke_width=5).set_points_as_corners(pts))
                if eps < x[-1]:
                    vec.add(Dot(ax.c2p(eps, np.interp(eps, x, y)), radius=0.08, color=fn_color(f)))
            out.add(fade_all(vec, g) if g < 1 else vec)
            return out

        pic = Group()
        pic.add_updater(lambda m: m.set_submobjects(list(picture())))
        grid_rect = Rectangle(width=grid_w, height=CONTENT.h).move_to(
            np.array([CONTENT.x0 + grid_w / 2, CONTENT.c[1], 0.]))
        report({'grid': grid_rect, 'curves': VGroup(ax, ls, le)})

        self.add(pic, ax, ls, le, *labs, key)
        if os.environ.get('STILL'):
            clk.set_value(float(os.environ['STILL']))
            self.update_mobjects(0)
            self.wait(0.04)
            return
        self.play(clk.animate.set_value(T_LOOP), run_time=T_LOOP, rate_func=linear)
