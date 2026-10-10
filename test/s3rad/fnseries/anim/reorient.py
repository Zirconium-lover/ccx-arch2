"""Looping GIF: stress reorientation of hydrides in a tube (the procedure of
Plyasov et al. 2023 / Egorochkin 2023: dissolve at 400 C, cool slowly).

Cycle 1: heat to 400 C - the hydrides dissolve; cool under hoop stress above
the threshold (sigma_th = 45 +/- 9 MPa for E635 guide tubes) - they
re-precipitate radially, and stay radial when the load is removed.
Cycle 2: heat again, cool without load - tangential again.  End = start.

Each hydride has its own solvus temperature, so on heating they vanish one by
one and on cooling they appear one by one; the orientation set is switched at
400 C, when none is visible.  Fn counts the visible ones by the 45 deg rule.
Fn ~ 0.1 without load and ~0.45 at 140 MPa are the values of the E635 rings
(Egorochkin 2023, tab. 3.1).

    manim -r 1280,720 --fps 25 reorient.py Reorient
"""
import numpy as np
from manim import *

config.tex_template = TexTemplate(preamble=r'''\usepackage[T2A]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage[english,russian]{babel}
\usepackage{amsmath}
\usepackage{amssymb}''')

BG = ManimColor('#000000')
RAD = ManimColor('#E8604C')
TAN = ManimColor('#2CB1D9')
HOT = ManimColor('#FF8A3D')
STRESS = ManimColor('#FFD24A')
WALL = ManimColor('#1A2430')
FONT = 'Noto Sans'

# schedule: (time s, temperature C, hoop stress MPa)
SCHED = [(0, 20, 0), (1.0, 20, 0), (3.5, 400, 0), (4.0, 400, 140), (4.5, 400, 140),
         (8.0, 20, 140), (8.5, 20, 140), (9.0, 20, 0), (11.5, 400, 0), (12.0, 400, 0),
         (15.5, 20, 0), (16.0, 20, 0)]
T_LOOP = 16.0
ts = np.array([s[0] for s in SCHED])


def temp(t):
    return float(np.interp(t, ts, [s[1] for s in SCHED]))


def stress(t):
    return float(np.interp(t, ts, [s[2] for s in SCHED]))


class DecimalRU(DecimalNumber):
    def _get_num_string(self, number):
        return super()._get_num_string(number).replace('.', ',')


class Reorient(Scene):
    def construct(self):
        self.camera.background_color = BG
        rng = np.random.default_rng(3)
        tt = ValueTracker(0.)

        # ---------------- tube ----------------
        C = np.array([-3.35, -0.2, 0.])
        r1, r2 = 1.55, 2.75
        ring = Annulus(inner_radius=r1, outer_radius=r2, fill_color=WALL, fill_opacity=1,
                       stroke_color=GREY_B, stroke_width=2).move_to(C)
        N = 64
        pos, tsol = [], []
        while len(pos) < N:
            r = rng.uniform(r1 + 0.18, r2 - 0.18)
            a = rng.uniform(0, 2 * PI)
            p = C + r * np.array([np.cos(a), np.sin(a), 0.])
            if all(np.linalg.norm(p - q) > 0.34 for q in pos):
                pos.append(p)
                tsol.append(rng.uniform(220, 385))
        tsol = np.array(tsol)
        # orientation sets: angle of the plate normal from the local radius
        def make_set(p_radial):
            radial = rng.random(N) < p_radial
            # radial plate = normal tangential (90 deg from R)
            return np.where(radial, 90 + rng.normal(0, 12, N), rng.normal(0, 12, N)), radial
        setA, radA = make_set(0.10)
        setB, radB = make_set(0.45)

        def hydrides():
            t = tt.get_value()
            T = temp(t)
            ang = setB if 3.5 <= t < 11.5 else setA
            g = VGroup()
            for i, p in enumerate(pos):
                if T >= tsol[i]:
                    continue
                rh = (p - C) / np.linalg.norm(p - C)
                a_r = np.arctan2(rh[1], rh[0])
                nrm = a_r + np.radians(ang[i])                 # plate normal
                d = np.array([-np.sin(nrm), np.cos(nrm), 0.])    # plate direction
                w = np.clip(np.sin(np.radians(ang[i])) ** 2, 0, 1)
                fade = np.clip((tsol[i] - T) / 25., 0, 1)      # grows in over 25 C
                g.add(Line(p - d * 0.16 * fade, p + d * 0.16 * fade, stroke_width=5,
                           color=interpolate_color(TAN, RAD, w)))
            return g
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
                base = C + u * (r2 + 0.32)
                g.add(Arrow(base, base + tng * 0.55, buff=0, color=STRESS, stroke_width=4,
                            max_tip_length_to_length_ratio=0.35).set_opacity(s),
                      Arrow(base, base - tng * 0.55, buff=0, color=STRESS, stroke_width=4,
                            max_tip_length_to_length_ratio=0.35).set_opacity(s))
                g.add(Arrow(C + u * 0.55, C + u * (r1 - 0.08), buff=0, color=GREY_A,
                            stroke_width=3, max_tip_length_to_length_ratio=0.3).set_opacity(0.8 * s))
            return g
        arrows = always_redraw(load)
        p_lab = always_redraw(lambda: MathTex(r'p', font_size=30, color=GREY_A).move_to(C).set_opacity(
            stress(tt.get_value()) / 140.))
        s_lab = always_redraw(lambda: MathTex(r'\sigma_\theta = 140\ \text{МПа}', font_size=34,
                                              color=STRESS).move_to(C + UP * 3.35).set_opacity(
            stress(tt.get_value()) / 140.))

        fn_lab = MathTex(r'F_n =', font_size=36).move_to(C + DOWN * 3.3 + LEFT * 0.4)
        fn_num = DecimalRU(0, num_decimal_places=2, font_size=36, color=WHITE)

        def fn_upd(m):
            t = tt.get_value()
            T = temp(t)
            rad = radB if 3.5 <= t < 11.5 else radA
            vis = T < tsol
            n = vis.sum()
            m.set_value(rad[vis].sum() / n if n else 0.)
            m.next_to(fn_lab, RIGHT, buff=0.15)
            m.set_opacity(1. if n > 0.15 * N else 0.)
        fn_num.add_updater(fn_upd)
        fn_num.update()

        # ---------------- T(t) and sigma(t) ----------------
        axT = Axes(x_range=[0, T_LOOP, 2], y_range=[0, 450, 100], x_length=5.6, y_length=2.6,
                   tips=False, axis_config=dict(color=GREY_B, stroke_width=2, include_ticks=False),
                   ).move_to(np.array([3.35, 1.35, 0.]))
        axS = Axes(x_range=[0, T_LOOP, 2], y_range=[0, 160, 50], x_length=5.6, y_length=1.3,
                   tips=False, axis_config=dict(color=GREY_B, stroke_width=2, include_ticks=False),
                   ).move_to(np.array([3.35, -1.75, 0.]))
        tl = MathTex(r'T,\ ^\circ\mathrm{C}', font_size=28).next_to(axT.y_axis, UP, buff=0.12)
        sl = MathTex(r'\sigma_\theta,\ \text{МПа}', font_size=28).next_to(axS.y_axis, UP, buff=0.12)
        t400 = MathTex('400', font_size=24, color=GREY_A).next_to(axT.c2p(0, 400), LEFT, buff=0.1)
        t20 = MathTex('20', font_size=24, color=GREY_A).next_to(axT.c2p(0, 20), LEFT, buff=0.1)
        Tcurve = VMobject(stroke_color=HOT, stroke_width=3).set_points_as_corners(
            [axT.c2p(s[0], s[1]) for s in SCHED])
        Scurve = VMobject(stroke_color=STRESS, stroke_width=3).set_points_as_corners(
            [axS.c2p(s[0], s[2]) for s in SCHED])
        th = DashedLine(axS.c2p(0, 45), axS.c2p(T_LOOP, 45), stroke_width=2, color=GREY_A)
        thl = MathTex(r'\sigma_{th} \approx 45\ \text{МПа}', font_size=24, color=GREY_A).next_to(
            th, DOWN, buff=0.06).align_to(th, RIGHT)
        s140 = MathTex('140', font_size=24, color=GREY_A).next_to(axS.c2p(0, 140), LEFT, buff=0.1)
        cursor = always_redraw(lambda: VGroup(
            DashedLine(axT.c2p(tt.get_value(), 0), axT.c2p(tt.get_value(), 450),
                       stroke_width=1.5, color=GREY_C),
            DashedLine(axS.c2p(tt.get_value(), 0), axS.c2p(tt.get_value(), 160),
                       stroke_width=1.5, color=GREY_C),
            Dot(axT.c2p(tt.get_value(), temp(tt.get_value())), radius=0.08, color=HOT),
            Dot(axS.c2p(tt.get_value(), stress(tt.get_value())), radius=0.08, color=STRESS)))
        note1 = Text('растворение', font=FONT, color=HOT).scale(0.32).next_to(
            axT.c2p(4.0, 400), UP, buff=0.08)
        note2 = Text('охлаждение\nпод нагрузкой', font=FONT, color=STRESS,
                     line_spacing=0.8).scale(0.3).move_to(axT.c2p(8.1, 300))
        note3 = Text('охлаждение\nбез нагрузки', font=FONT, color=GREY_A,
                     line_spacing=0.8).scale(0.3).move_to(axT.c2p(14.6, 330))
        r1lab = Text('радиальные остаются', font=FONT, color=RAD).scale(0.3).next_to(
            axS.c2p(8.7, 0), DOWN, buff=0.12)

        self.add(ring, hyd, arrows, p_lab, s_lab, fn_lab, fn_num,
                 axT, axS, tl, sl, t400, t20, Tcurve, Scurve, th, thl, s140, cursor,
                 note1, note2, note3, r1lab)
        self.play(tt.animate.set_value(T_LOOP - 0.04), run_time=T_LOOP - 0.04, rate_func=linear)
