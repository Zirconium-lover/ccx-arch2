"""Looping GIF: stress reorientation of hydrides in a tube wall (the procedure
of Plyasov et al. 2023 / Egorochkin 2023): dissolve at 400 C, apply hoop
stress above the threshold (sigma_th = 45 +/- 9 MPa for E635), cool slowly -
the hydrides re-precipitate radially and stay radial when the load is off.
No claim is made about later stress-free cycles (memory effect); the loop
returns to the start by an explicit cross-fade labelled as the initial state.

Fn ~ 0.1 initially and ~0.45 after one reorientation at 140 MPa are the
values of the E635 rings (Egorochkin 2023, tab. 3.1).

    manim -r 1600,900 --fps 25 reorient.py Reorient
"""
import numpy as np
from style import *

# (time s, temperature C, hoop stress MPa)
SCHED = [(0, 20, 0), (1.0, 20, 0), (3.2, 400, 0), (3.4, 400, 0), (3.8, 400, 140),
         (4.0, 400, 140), (7.4, 20, 140), (7.9, 20, 140), (8.4, 20, 0), (12.0, 20, 0)]
ts = np.array([s[0] for s in SCHED])
T_END = SCHED[-1][0]
RESET = 1.2                     # cross-fade back to the initial state


def temp(t):
    return float(np.interp(t, ts, [s[1] for s in SCHED]))


def stress(t):
    return float(np.interp(t, ts, [s[2] for s in SCHED]))


def status(t):
    if t < 1.0:
        return 'исходная структура', TAN
    if t < 3.4:
        return 'нагрев до 400 °C: растворение', HOT
    if t < 4.0:
        return 'нагрузка выше порога', ACC
    if t < 7.9:
        return 'охлаждение под нагрузкой', ACC
    if t < 8.4:
        return 'нагрузка снята', WHITE
    return 'радиальные остаются', RAD


class Reorient(Scene):
    def construct(self):
        self.camera.background_color = BG
        rng = np.random.default_rng(3)
        tt = ValueTracker(0.)
        wA = ValueTracker(0.)          # weight of the initial set during the reset
        cur_op = ValueTracker(1.)

        # ---------------- tube wall ----------------
        C = np.array([-3.55, -0.45, 0.])
        r1, r2 = 1.5, 3.05
        ring = Annulus(inner_radius=r1, outer_radius=r2, fill_color=WALL, fill_opacity=1,
                       stroke_color=GREY_B, stroke_width=3).move_to(C)
        N = 78
        pos, tsol = [], []
        while len(pos) < N:
            r = rng.uniform(r1 + 0.22, r2 - 0.22)
            a = rng.uniform(0, 2 * PI)
            p = C + r * np.array([np.cos(a), np.sin(a), 0.])
            if all(np.linalg.norm(p - q) > 0.4 for q in pos):
                pos.append(p)
                tsol.append(rng.uniform(220, 385))
        tsol = np.array(tsol)

        def make_set(p_radial):
            radial = np.zeros(N, bool)
            radial[rng.permutation(N)[:int(round(p_radial * N))]] = True
            return np.where(radial, 90 + rng.normal(0, 10, N), rng.normal(0, 10, N)), radial
        setA, radA = make_set(0.10)
        setB, radB = make_set(0.45)

        def bars(ang, T, op):
            g = VGroup()
            for i, p in enumerate(pos):
                if T >= tsol[i] or op <= 0.01:
                    continue
                rh = (p - C) / np.linalg.norm(p - C)
                nrm = np.arctan2(rh[1], rh[0]) + np.radians(ang[i])
                d = np.array([-np.sin(nrm), np.cos(nrm), 0.])
                w = np.clip(np.sin(np.radians(ang[i])) ** 2, 0, 1)
                f = np.clip((tsol[i] - T) / 25., 0, 1)
                g.add(Line(p - d * 0.2 * f, p + d * 0.2 * f, stroke_width=8,
                           color=interpolate_color(TAN, RAD, w)).set_opacity(op))
            return g

        def hydrides():
            t, w = tt.get_value(), wA.get_value()
            T = temp(t)
            if w > 0:                                  # reset cross-fade
                return VGroup(bars(setB, 20, 1 - w), bars(setA, 20, w))
            return bars(setB if t >= 3.4 else setA, T, 1.)
        hyd = always_redraw(hydrides)

        def load():
            s = stress(tt.get_value()) / 140.
            g = VGroup()
            if s < 1e-3:
                return g
            for k in range(8):
                a = k * PI / 4 + PI / 8
                u = np.array([np.cos(a), np.sin(a), 0.])
                tng = np.array([-u[1], u[0], 0.])
                base = C + u * (r2 + 0.3)
                for sg in (1, -1):
                    g.add(Arrow(base, base + sg * tng * 0.6, buff=0, color=ACC, stroke_width=6,
                                max_tip_length_to_length_ratio=0.38).set_opacity(s))
                g.add(Arrow(C + u * 0.45, C + u * (r1 - 0.1), buff=0, color=GREY_A,
                            stroke_width=5, max_tip_length_to_length_ratio=0.3).set_opacity(0.85 * s))
            return g
        arrows = always_redraw(load)

        # ---------------- T(t), sigma(t) ----------------
        X0 = 1.85
        axT = Axes(x_range=[0, T_END, 2], y_range=[0, 450, 100], x_length=4.4, y_length=2.35,
                   tips=False, axis_config=dict(color=GREY_B, stroke_width=3, include_ticks=False)
                   ).move_to(np.array([X0 + 1.9, 1.35, 0.]))
        axS = Axes(x_range=[0, T_END, 2], y_range=[0, 160, 50], x_length=4.4, y_length=1.35,
                   tips=False, axis_config=dict(color=GREY_B, stroke_width=3, include_ticks=False)
                   ).move_to(np.array([X0 + 1.9, -1.0, 0.]))
        tl = tex(r'T,\ ^\circ\mathrm{C}', FS_SM).next_to(axT.y_axis, LEFT, buff=0.2).shift(UP * 0.75)
        sl = tex(r'\sigma_\theta,\ \text{МПа}', FS_SM).next_to(axS.y_axis, LEFT, buff=0.2).shift(UP * 0.35)
        t400 = tex('400', FS_SM - 6, GREY_A).next_to(axT.c2p(T_END, 400), RIGHT, buff=0.15)
        s140 = tex(r'140', FS_SM - 6, GREY_A).next_to(axS.c2p(T_END, 140), RIGHT,
                                                                buff=0.15).shift(DOWN * 0.05)
        Tcurve = VMobject(stroke_color=HOT, stroke_width=5).set_points_as_corners(
            [axT.c2p(s[0], s[1]) for s in SCHED])
        Scurve = VMobject(stroke_color=ACC, stroke_width=5).set_points_as_corners(
            [axS.c2p(s[0], s[2]) for s in SCHED])
        th = DashedLine(axS.c2p(0, 45), axS.c2p(T_END, 45), stroke_width=3, color=GREY_A)
        thl = tex(r'\sigma_{th}\approx 45', FS_SM - 8, GREY_A).next_to(axS.c2p(8.7, 45),
                                                                      UP, buff=0.08, aligned_edge=LEFT)
        cursor = always_redraw(lambda: VGroup(
            Dot(axT.c2p(tt.get_value(), temp(tt.get_value())), radius=0.12, color=HOT),
            Dot(axS.c2p(tt.get_value(), stress(tt.get_value())), radius=0.12, color=ACC),
            DashedLine(axT.c2p(tt.get_value(), 0), axT.c2p(tt.get_value(), 450),
                       stroke_width=2, color=GREY_C),
            DashedLine(axS.c2p(tt.get_value(), 0), axS.c2p(tt.get_value(), 160),
                       stroke_width=2, color=GREY_C)).set_opacity(cur_op.get_value()))

        def stat():
            if wA.get_value() > 0:
                return txt('↺ исходная структура', TX + 3, TAN).move_to(np.array([X0 + 1.9, -2.55, 0.]))
            s, c = status(tt.get_value())
            return txt(s, TX + 3, c).move_to(np.array([X0 + 1.9, -2.55, 0.]))
        st = always_redraw(stat)

        fn_lab = tex(r'F_n =', FS_BIG).move_to(np.array([X0 + 1.2, -3.45, 0.]))
        fn_num = DecimalRU(0, num_decimal_places=2, font_size=FS_BIG, color=WHITE)

        def fn_upd(m):
            t, w = tt.get_value(), wA.get_value()
            if w > 0:
                v = (1 - w) * radB.mean() + w * radA.mean()
                m.set_value(v).set_opacity(1)
            else:
                rad = radB if t >= 3.4 else radA
                vis = temp(t) < tsol
                n = vis.sum()
                m.set_value(rad[vis].sum() / n if n else 0.)
                m.set_opacity(1. if n > 0.2 * N else 0.)
            m.next_to(fn_lab, RIGHT, buff=0.2)
        fn_num.add_updater(fn_upd)
        fn_num.update()

        self.add(ring, hyd, arrows, axT, axS, tl, sl, t400, s140, Tcurve, Scurve, th, thl,
                 cursor, st, fn_lab, fn_num)
        self.play(tt.animate.set_value(T_END), run_time=T_END, rate_func=linear)
        # back to the start: fade the cursor, cross-fade the structure
        self.play(cur_op.animate.set_value(0.), wA.animate.set_value(1.), run_time=RESET * 0.6)
        tt.set_value(0.)
        wA.set_value(0.)
        self.play(cur_op.animate.set_value(1.), run_time=RESET * 0.4 - 0.04)
