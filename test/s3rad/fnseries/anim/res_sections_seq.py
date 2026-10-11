"""Looping GIF from the computed model: the sections of seven runs (Fn = 0.33
left out: its section cuts the same plates as 0.42)
(z = 0.61 mm, one plane and one scale), loaded one after another, Fn growing:
the active section is framed in the colour of its curve and its curve runs
on the right; the finished ones stay broken with their curves, the waiting
ones stand intact and dimmed.  omega as grey steps with isolines (style
bl_grey of res_section.py).

    manim -r 1600,900 --fps 25 --disable_caching res_sections_seq.py ResSectionsSeq
"""
import os
from mech_common import *
from fndata import Run, SERIES
from res_section import FieldRaster, STYLES, U_MAX
from res_curves import fn_color, SIG_MAX
from res_sections_all import frac_tex

SKIP = (0.33,)                        # Fn = 4/12: its section shows the same plates as 5/12
N = sum(1 for fn, _ in SERIES if fn not in SKIP)
T0 = 0.6                              # the first run starts
D, G = 2.0, 0.3                       # one run, the pause after it
T_END = T0 + N * (D + G)
T_DIP = (T_END + 1.4, T_END + 2.1, T_END + 2.4, T_END + 3.1)
T_LOOP = T_DIP[3] + 0.4
EPS_AX = 0.135
DIM = 0.4                             # waiting sections
ST = STYLES[os.environ.get('PAL', 'rb_smooth_75')]


class ResSectionsSeq(Scene):
    def construct(self):
        self.camera.background_color = BG
        runs = [Run(fn) for fn, _ in SERIES if fn not in SKIP]
        F = [r.k_radial / 12 for r in runs]
        z_cut = float(os.environ.get('ZCUT', 0)) or runs[0].best_section('z')
        r0 = runs[0]
        LX, LY = r0.hi[0] - r0.lo[0], r0.hi[1] - r0.lo[1]

        # ---------------- layout: a 3 x 3 grid of sections | the curves ----------------
        labs = [tex(frac_tex(r.k_radial), L3M, fn_color(f)) for r, f in zip(runs, F)]
        hl = max(m.height for m in labs)
        GX, GY, GL = 0.3, 0.28, 0.1                     # gaps: columns, rows, label - section
        CI = CONTENT.inset(0.08)                          # room for the frame of the active section
        ch = (CI.h - 2 * GY) / 3 - hl - GL                # section height
        sc = ch / LY
        cw = (LX + U_MAX) * sc                            # cell width (room for the stretch)
        grid_w = 3 * cw + 2 * GX
        cellR = Box(CI.x0 + grid_w + 0.7, CONTENT.x1, CONTENT.y0, CONTENT.y1)
        origins = []                                      # lower-left of each section (scene)
        for i in range(9):
            c, rw = i % 3, i // 3
            x0 = CI.x0 + c * (cw + GX)
            ytop = CI.y1 - rw * (hl + GL + ch + GY)
            origins.append(np.array([x0, ytop - hl - GL - ch, 0.]))
        for i, lab in enumerate(labs):
            o = origins[i]
            lab.move_to(np.array([o[0], o[1] + ch + GL, 0.]), aligned_edge=DL)

        def Pk(i, xy):
            xy = np.asarray(xy)
            o = origins[i]
            return np.c_[o[0] + (xy[..., 0] - r0.lo[0]) * sc, o[1] + (xy[..., 1] - r0.lo[1]) * sc,
                         np.zeros(xy.shape[:-1])].reshape(*xy.shape[:-1], 3)

        boxes = [Box(o[0] - 0.04, o[0] + cw + 0.04, o[1] - 0.04, o[1] + ch + 0.04) for o in origins[:N]]
        rasters = [FieldRaster(b, ST) for b in boxes]

        # the omega key in the free cells of the last row
        key_box = Box(origins[N][0] + 0.15, origins[8][0] + cw - 0.1, origins[8][1], origins[8][1] + ch + hl + GL)
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

        # ---------------- the curves: axes with titles, units and numbers ----------------
        # Fn of each finished curve, in a row above the graph (values in the curve's colour)
        vals = [tex('0' if f == 0 else '1' if f == 1 else ('%.2f' % f).replace('.', '{,}'), L3M, fn_color(f))
                for f in F]
        row_h = max(v.height for v in vals)
        ls = tex(r'\sigma,\ \text{МПа}', L2M).rotate(PI / 2)     # upright along the axis: the top stays free
        le = tex(r'\varepsilon,\ \%', L2M)
        YT, XT = [0, 100, 200, 300, 400, 500], [2, 4, 6, 8, 10, 12]       # one 0, at the corner
        ynum = [tex(str(v), L3M, GREY_A) for v in YT]
        xnum = [tex(str(v), L3M, GREY_A) for v in XT]
        yw = max(m.width for m in ynum)
        xh = max(m.height for m in xnum)
        x0 = cellR.x0 + ls.width + 0.2 + yw + 0.15        # the sigma axis
        y0 = CONTENT.y0 + xh + 0.12 + le.height + 0.12    # the eps axis: numbers, then its title
        y_top = CONTENT.y1 - row_h - 0.55                 # room for the row and the leaders
        ax = Axes(x_range=[0, EPS_AX, 0.02], y_range=[0, SIG_MAX, 100],
                  x_length=CONTENT.x1 - 0.25 - x0, y_length=y_top - 0.2 - y0, tips=True,
                  axis_config=dict(color=GREY_B, stroke_width=3, include_ticks=True, tick_size=0.06,
                                   tip_width=0.2, tip_height=0.2))
        ax.shift(np.array([x0, y0, 0.]) - ax.c2p(0, 0))
        for v, m in zip(YT, ynum):
            m.next_to(ax.c2p(0, v), LEFT, buff=0.15)
        for v, m in zip(XT, xnum):
            m.next_to(ax.c2p(v / 100, 0), DOWN, buff=0.12)
        ls.move_to(np.array([cellR.x0 + ls.width / 2, ax.c2p(0, SIG_MAX / 2)[1], 0.]))
        le.move_to(np.array([CONTENT.x1, CONTENT.y0, 0.]), aligned_edge=DR)   # under the numbers, at the end
        nums = VGroup(*ynum, *xnum)
        cx = [np.r_[0, r.curve[:, 0]] for r in runs]
        cy = [np.r_[0, r.curve[:, 1]] for r in runs]
        e_end = [x[-1] for x in cx]
        # one row, Fn = 1 ... 0 (the order of the knees of the curves along eps), packed with
        # an even gap between the labels and centred over the knees; the dashed leaders
        # (knee -> its value) hardly cross
        order = sorted(range(N), key=lambda i: -F[i])
        knees = [ax.c2p(runs[i].eps_peak, 0)[0] for i in order]
        gap = 0.24
        ws = [vals[i].width for i in order]
        lefts = np.cumsum([0.] + [w + gap for w in ws[:-1]])
        span = lefts[-1] + ws[-1]
        x_l = (min(knees) + max(knees)) / 2 - span / 2
        x_l = min(max(x_l, cellR.x0 + 0.1), CONTENT.x1 - 0.1 - span)
        for j, i in enumerate(order):
            vals[i].move_to(np.array([x_l + lefts[j] + ws[j] / 2, CONTENT.y1 - row_h / 2, 0.]))
            vals[i].set_background_stroke(color=BLACK, width=8)   # cuts the leaders passing close by
        leaders = []
        for i, r in enumerate(runs):
            p0 = ax.c2p(r.eps_peak, r.peak) + UP * 0.08
            leaders.append(DashedLine(p0, vals[i].get_bottom() + DOWN * 0.1, dash_length=0.07, stroke_width=2.5,
                                      color=fn_color(F[i])))

        clk = ValueTracker(0.)

        def dip(t):
            if t < T_DIP[0]:
                return 1.
            return 1 - ramp(t, T_DIP[0], T_DIP[1]) + ramp(t, T_DIP[2], T_DIP[3])

        def te_of(t):
            return t if t < T_DIP[1] + 0.15 else 0.

        def eps_run(i, te):
            """Grip strain of run i at time te: 0 before its turn, its last
            strain after it; quick through the uniform start."""
            a = T0 + i * (D + G)
            if te <= a:
                return 0.
            tau = smoothstep((min(te, a + D) - a) / D)
            return float(np.interp(tau, [0, 0.2, 1], [0, 0.03, e_end[i]]))

        def polys_vm(polys, color, sw=0.5):
            pts = []
            for q in polys:
                b = np.roll(q, -1, axis=0)
                pts.append(np.stack([q, q + (b - q) / 3, q + 2 * (b - q) / 3, b], 1).reshape(-1, 3))
            vm = VMobject(stroke_color=color, stroke_width=sw, fill_color=color, fill_opacity=1)
            if pts:
                vm.set_points(np.vstack(pts))
            return vm

        cache = {}

        def section(i, eps):
            """Raster of the matrix and the hydride paths of section i at
            eps; the intact and the broken states are kept."""
            key_ = (i, round(eps, 7))
            if key_ in cache:
                return cache[key_]
            r = runs[i]
            s = r.at(eps)
            polys, own, pv = r.slice('z', z_cut, s, 'omega')
            ph = r.phase[own]
            mat = np.flatnonzero(ph == 0)
            img = rasters[i].render([Pk(i, polys[j][:, :2])[:, :2] for j in mat], [pv[j] for j in mat])
            hyd = [(polys_vm([Pk(i, polys[j][:, :2]) for j in np.flatnonzero(ph == code)], col))
                   for code, col in ((1, RAD), (2, TAN)) if np.any(ph == code)]
            # a thin dark outline around each plate's cut (its outer edges only)
            cnt = {}
            for j in np.flatnonzero(ph > 0):
                q = Pk(i, polys[j][:, :2])
                for a, b in zip(q, np.roll(q, -1, axis=0)):
                    ka, kb = tuple(np.round(a[:2], 4)), tuple(np.round(b[:2], 4))
                    kk = (ka, kb) if ka < kb else (kb, ka)
                    cnt[kk] = (cnt[kk][0] + 1, a, b) if kk in cnt else (1, a, b)
            seg = [np.array([a, a + (b - a) / 3, a + 2 * (b - a) / 3, b]) for n_, a, b in cnt.values() if n_ == 1]
            if seg:
                hyd.append(VMobject(stroke_color=BLACK, stroke_width=1.2).set_points(np.vstack(seg)))
            out = (img, hyd)
            if eps == 0. or eps == e_end[i]:
                cache[key_] = out
            return out

        def picture():
            t = clk.get_value()
            g = dip(t)
            te = te_of(t)
            out = Group()
            vec = VGroup()
            for i in range(N):
                eps = eps_run(i, te)
                a = T0 + i * (D + G)
                k = g * (DIM + (1 - DIM) * ramp(te, a - 0.3, a))       # lights up just before its turn
                img, hyd = section(i, eps)
                if k < 1:
                    img = img.copy(); img[..., 3] = (img[..., 3] * k).astype(np.uint8)
                im = ImageMobject(img)
                im.stretch_to_fit_width(boxes[i].w).stretch_to_fit_height(boxes[i].h).move_to(boxes[i].c)
                out.add(im)
                for h in hyd:
                    vec.add(fade_all(h.copy(), k) if k < 1 else h.copy())
                # the frame of the active section
                on = ramp(te, a - 0.3, a) * (1 - ramp(te, a + D, a + D + G + 0.2))
                if on > 0.01:
                    fr = Rectangle(width=cw + 0.12, height=ch + 0.12, stroke_color=fn_color(F[i]),
                                   stroke_width=4).move_to(np.array([origins[i][0] + cw / 2, origins[i][1] + ch / 2, 0.]))
                    vec.add(fade_all(fr, on * g))
                # its curve
                if eps > 1e-4:
                    m = cx[i] < eps
                    pts = [ax.c2p(x_, y_) for x_, y_ in zip(cx[i][m], cy[i][m])] + \
                          [ax.c2p(eps, np.interp(eps, cx[i], cy[i]))]
                    active = te < a + D
                    ln = VMobject(stroke_color=fn_color(F[i]), stroke_width=6 if active else 4).set_points_as_corners(pts)
                    vec.add(fade_all(ln, g) if g < 1 else ln)
                    if active:
                        d = Dot(ax.c2p(eps, np.interp(eps, cx[i], cy[i])), radius=0.1, color=fn_color(F[i]))
                        vec.add(fade_all(d, g) if g < 1 else d)
            out.add(vec)
            return out

        pic = Group()
        pic.add_updater(lambda m: m.set_submobjects(list(picture())))

        def lab_up(i):
            def f(m):
                te = te_of(clk.get_value())
                a = T0 + i * (D + G)
                m.set_opacity((DIM + (1 - DIM) * ramp(te, a - 0.3, a)) * dip(clk.get_value()))
            return f
        for i, lab in enumerate(labs):
            lab.add_updater(lab_up(i))

        def done_up(i, stroke):
            def f(m):
                t = clk.get_value()
                a = T0 + i * (D + G)
                k = ramp(te_of(t), a + D, a + D + 0.4) * dip(t)
                if stroke:
                    m.set_stroke(opacity=k)
                else:
                    m.set_opacity(k)
            return f
        for i in range(N):
            vals[i].add_updater(done_up(i, False))
            leaders[i].add_updater(done_up(i, True))
        grid_rect = Rectangle(width=grid_w + 0.12, height=CI.h + 0.12).move_to(
            np.array([CI.x0 + grid_w / 2, CI.c[1], 0.]))
        report({'grid': grid_rect, 'curves': VGroup(ax, ls, le, nums), 'row': VGroup(*vals)})

        self.add(pic, ax, ls, le, nums, *labs, key, *leaders, *vals)
        if os.environ.get('STILL'):
            clk.set_value(float(os.environ['STILL']))
            self.update_mobjects(0)
            self.wait(0.04)
            return
        self.play(clk.animate.set_value(T_LOOP), run_time=T_LOOP, rate_func=linear)
