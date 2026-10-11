"""Looping GIF from the computed model: the stress-strain curves of the whole
series drawn together (colour: blue Fn = 0 -> red Fn = 1).  Each curve
leaves a dot at its fracture strain (eps_f10: the force has dropped to 10 %
of its peak) on the right panel, eps_f against Fn.

    manim -r 1600,900 --fps 25 --disable_caching res_curves.py ResCurves
"""
import os
from mech_common import *
from fndata import Run, SERIES

T_LOOP = 12.0
T_RUN = (0.6, 7.6)                    # the strain sweep
T_LINK = (8.0, 8.8)                   # the dots are joined
T_DIP = (9.8, 10.5, 10.8, 11.5)
EPS_MAX, SIG_MAX = 0.135, 650.
EF_MAX = 0.155
MID = ManimColor('#E6E6E6')


def fn_color(f):
    """Blue (tangential) -> light -> red (radial), like the hydride colours."""
    return interpolate_color(TAN, MID, f / 0.5) if f <= 0.5 else interpolate_color(MID, RAD, (f - 0.5) / 0.5)


class ResCurves(Scene):
    def construct(self):
        self.camera.background_color = BG
        runs = [Run(fn) for fn, _ in SERIES]
        F = [r.k_radial / 12 for r in runs]

        cellL = Box(CONTENT.x0, CONTENT.x0 + 7.6, CONTENT.y0, CONTENT.y1)
        cellR = Box(cellL.x1 + 0.7, CONTENT.x1, CONTENT.y0, CONTENT.y1)
        axcfg = dict(color=GREY_B, stroke_width=3, include_ticks=False, tip_width=0.2, tip_height=0.2)

        # both graphs stand on one baseline (room under it for the 0 / 1 of Fn)
        t0, t1 = tex('0', L3M, GREY_A), tex('1', L3M, GREY_A)
        y_base = CONTENT.y0 + t0.height + 0.3

        # ---------------- left: sigma - eps ----------------
        ls, le = tex(r'\sigma', L1), tex(r'\varepsilon', L1)
        axL = Axes(x_range=[0, EPS_MAX, EPS_MAX], y_range=[0, SIG_MAX, SIG_MAX],
                   x_length=cellL.w - ls.width - 0.3 - le.width - 0.3, y_length=cellL.y1 - 0.06 - y_base,
                   tips=True, axis_config=axcfg)
        axL.move_to(np.array([cellL.x0 + ls.width + 0.3 + axL.width / 2, y_base + axL.height / 2, 0.]))
        axL.shift(UP * (y_base - axL.c2p(0, 0)[1]))
        ls.next_to(axL.y_axis, LEFT, buff=0.2).align_to(axL.y_axis, UP)
        le.next_to(axL.x_axis, RIGHT, buff=0.15)

        # ---------------- right: eps_f - Fn ----------------
        lf, lF = tex(r'\varepsilon_f', L1), tex(r'F_n', L1)
        gw = cellR.w - lf.width - 0.3 - lF.width - 0.3
        axR = Axes(x_range=[0, 1.12, 1.12], y_range=[0, EF_MAX, EF_MAX], x_length=gw,
                   y_length=min(cellL.y1 - 0.06 - y_base, 1.45 * gw), tips=True, axis_config=axcfg)
        axR.move_to(np.array([cellR.x0 + lf.width + 0.3 + axR.width / 2, y_base + axR.height / 2, 0.]))
        axR.shift(UP * (y_base - axR.c2p(0, 0)[1]))
        lf.next_to(axR.y_axis, LEFT, buff=0.2).align_to(axR.y_axis, UP)
        lF.next_to(axR.x_axis, RIGHT, buff=0.15)
        t0.next_to(axR.c2p(0, 0), DOWN, buff=0.15)
        t1.next_to(axR.c2p(1, 0), DOWN, buff=0.15)
        tick1 = Line(axR.c2p(1, 0) + DOWN * 0.08, axR.c2p(1, 0) + UP * 0.08, color=GREY_B, stroke_width=3)

        clk = ValueTracker(0.)

        def dip(t):
            if t < T_DIP[0]:
                return 1.
            return 1 - ramp(t, T_DIP[0], T_DIP[1]) + ramp(t, T_DIP[2], T_DIP[3])

        def te_of(t):
            return t if t < T_DIP[1] + 0.15 else 0.

        def eps_of(te):
            return EPS_MAX * smoothstep((te - T_RUN[0]) / (T_RUN[1] - T_RUN[0]))

        cx = [np.r_[0, r.curve[:, 0]] for r in runs]
        cy = [np.r_[0, r.curve[:, 1]] for r in runs]

        def picture():
            t = clk.get_value()
            te = te_of(t)
            eps = eps_of(te)
            g = VGroup()
            # curves, Fn = 1 drawn last (on top)
            for r, x, y, f in zip(runs, cx, cy, F):
                e = min(eps, x[-1])
                if e <= 1e-4:
                    continue
                m = x < e
                pts = [axL.c2p(a, b) for a, b in zip(x[m], y[m])] + [axL.c2p(e, np.interp(e, x, y))]
                g.add(VMobject(stroke_color=fn_color(f), stroke_width=5).set_points_as_corners(pts))
            # the dots: born when a curve passes its fracture strain
            link = ramp(te, *T_LINK)
            done = [(f, r.eps_f10) for r, f in zip(runs, F) if eps >= r.eps_f10]
            if link > 0 and len(done) == len(runs):
                pts = [axR.c2p(f, e) for f, e in sorted(done)]
                ln = VMobject(stroke_color=GREY_B, stroke_width=4).set_points_smoothly(pts)
                g.add(ln.pointwise_become_partial(ln.copy(), 0, link) if link < 1 else ln)
            for r, f in zip(runs, F):
                if eps >= r.eps_f10:
                    a = ramp(eps, r.eps_f10, r.eps_f10 + 0.006)
                    pop = 1 + 0.6 * np.sin(np.pi * a)
                    g.add(Dot(axR.c2p(f, r.eps_f10), radius=0.13 * pop, color=fn_color(f)))
                    g.add(Dot(axL.c2p(r.eps_f10, 0.1 * r.peak), radius=0.09, color=fn_color(f)))
            k = dip(t)
            return fade_all(g, k) if k < 1 else g

        pic = VGroup()
        pic.add_updater(lambda m: m.set_submobjects(list(picture())))
        report({'left': VGroup(axL, ls, le), 'right': VGroup(axR, lf, lF, t0, t1)})

        self.add(axL, ls, le, axR, lf, lF, t0, t1, tick1, pic)
        if os.environ.get('STILL'):
            clk.set_value(float(os.environ['STILL']))
            self.update_mobjects(0)
            self.wait(0.04)
            return
        self.play(clk.animate.set_value(T_LOOP), run_time=T_LOOP, rate_func=linear)
