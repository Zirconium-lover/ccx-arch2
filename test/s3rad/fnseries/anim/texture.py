"""Looping GIF: how the tube texture sets the hydride orientation.

Each grain's basal pole (c axis) sits at phi + spread from the radial
direction R; the hydride lies in the basal plane, i.e. normal to c.  phi
sweeps 0 -> 90 deg -> 0: basal poles along R give tangential hydrides,
basal poles along T give radial ones.

Right: (0001) pole figure with R at the centre (stereographic, r = tan(psi/2)),
and Fn(phi) = share of hydrides whose normal is > 45 deg from R (the
experimental definition), computed for the same spread - Zaimovsky et al.'s
step.  f_R = <cos^2 psi_R> is the Kearns parameter of the same grains.

    manim -r 1280,720 --fps 25 texture.py Texture
"""
import numpy as np
from math import erf
from manim import *

config.tex_template = TexTemplate(preamble=r'''\usepackage[T2A]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage[english,russian]{babel}
\usepackage{amsmath}
\usepackage{amssymb}''')

BG = ManimColor('#000000')
RAD = ManimColor('#E8604C')
TAN = ManimColor('#2CB1D9')
POLE = ManimColor('#FFD24A')
WALL = ManimColor('#1A2430')
FONT = 'Noto Sans'
SPREAD = 15.0            # deg, texture spread of the basal poles
T_LOOP = 14.0


class DecimalRU(DecimalNumber):
    """DecimalNumber with a decimal comma."""
    def _get_num_string(self, number):
        return super()._get_num_string(number).replace('.', ',')


def fn_of(phi_deg):
    """P(|phi + N(0, SPREAD)| > 45) - the expected Fn for this texture."""
    s = SPREAD * np.sqrt(2)
    return 0.5 * (1 - erf((45 - phi_deg) / s)) + 0.5 * (1 - erf((45 + phi_deg) / s))


class Texture(Scene):
    def construct(self):
        self.camera.background_color = BG
        rng = np.random.default_rng(7)

        # ---------------- tube wall with grains ----------------
        O = np.array([-3.3, -6.3, 0.])
        R1, R2 = 5.6, 7.4
        A0, A1 = 66, 114
        wall = AnnularSector(inner_radius=R1, outer_radius=R2, angle=np.radians(A1 - A0),
                             start_angle=np.radians(A0), arc_center=O,
                             fill_color=WALL, fill_opacity=1, stroke_color=GREY_B,
                             stroke_width=2)
        pts = []
        tries = 0
        while len(pts) < 18 and tries < 20000:
            tries += 1
            r = rng.uniform(R1 + 0.3, R2 - 0.3)
            t = np.radians(rng.uniform(A0 + 3, A1 - 3))
            p = O + r * np.array([np.cos(t), np.sin(t), 0.])
            if all(np.linalg.norm(p - q) > 0.58 for q in pts):
                pts.append(p)
        dev = rng.normal(0, SPREAD, len(pts))
        phi = ValueTracker(0.)

        def grain(i):
            p = pts[i]
            rhat = (p - O) / np.linalg.norm(p - O)
            ang_r = np.arctan2(rhat[1], rhat[0])

            def build():
                psi = np.radians(phi.get_value() + dev[i])      # c axis from R
                cdir = np.array([np.cos(ang_r + psi), np.sin(ang_r + psi), 0.])
                pdir = np.array([-cdir[1], cdir[0], 0.])          # basal trace
                # colour: how far the plate normal (= c) is from R
                w = np.clip(abs(np.sin(psi)) ** 2, 0, 1)
                plate = Line(p - pdir * 0.24, p + pdir * 0.24, stroke_width=7,
                             color=interpolate_color(TAN, RAD, w))
                ax = Line(p - cdir * 0.2, p + cdir * 0.2, stroke_width=1.6,
                          color=POLE, stroke_opacity=0.85)
                return VGroup(ax, plate)
            return always_redraw(build)

        grains = VGroup(*[grain(i) for i in range(len(pts))])
        mid = O + (R2 + 0.15) * np.array([np.cos(np.radians(90)), np.sin(np.radians(90)), 0.])
        r_arrow = Arrow(mid, mid + UP * 0.8, buff=0, color=WHITE, stroke_width=4)
        r_lab = MathTex('R', font_size=34).next_to(r_arrow, LEFT, buff=0.1)
        t_arrow = Arrow(mid + RIGHT * 0.12, mid + RIGHT * 0.95, buff=0, color=WHITE,
                        stroke_width=4)
        t_lab = MathTex('T', font_size=34).next_to(t_arrow, RIGHT, buff=0.1)
        legend = VGroup(
            VGroup(Line(LEFT * 0.18, RIGHT * 0.18, color=POLE, stroke_width=2),
                   Text('базисный полюс (ось c)', font=FONT).scale(0.32)).arrange(RIGHT, buff=0.12),
            VGroup(Line(LEFT * 0.18, RIGHT * 0.18, color=GREY_A, stroke_width=7),
                   Text('гидрид в базисной плоскости', font=FONT).scale(0.32)).arrange(RIGHT, buff=0.12),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.12).to_corner(DL, buff=0.35)

        # ---------------- pole figure ----------------
        pc = np.array([3.55, 1.75, 0.])
        Rp = 1.25
        pf = Circle(radius=Rp, stroke_color=GREY_B, stroke_width=2).move_to(pc)
        cross = VGroup(Line(pc + LEFT * Rp, pc + RIGHT * Rp, stroke_width=1, color=GREY_D),
                       Line(pc + DOWN * Rp, pc + UP * Rp, stroke_width=1, color=GREY_D))
        pf_lab = MathTex('(0001)', font_size=30).next_to(pf, UP, buff=0.1)
        pR = MathTex('R', font_size=26, color=GREY_A).move_to(pc + DR * 0.17)
        pT = MathTex('T', font_size=26, color=GREY_A).next_to(pf, RIGHT, buff=0.42)
        pL = MathTex('L', font_size=26, color=GREY_A).next_to(pf, DOWN, buff=0.05)

        def blobs():
            x = np.tan(np.radians(phi.get_value()) / 2) * Rp
            g = VGroup()
            for sgn in (-1, 1):
                for k, (rr, op) in enumerate(((0.36, 0.18), (0.24, 0.35), (0.13, 0.9))):
                    g.add(Ellipse(width=2 * rr, height=2 * rr * 1.25, stroke_width=0,
                                  fill_color=POLE, fill_opacity=op).move_to(pc + RIGHT * sgn * x))
            return g
        pole = always_redraw(blobs)

        # ---------------- Fn(phi) curve ----------------
        ax = Axes(x_range=[0, 90, 15], y_range=[0, 1, 0.5], x_length=4.2, y_length=2.1,
                  tips=False, axis_config=dict(color=GREY_B, stroke_width=2, include_ticks=True,
                                               font_size=22),
                  ).move_to(np.array([3.75, -1.75, 0.]))
        ax.add_coordinates([0, 45, 90], [0, 0.5, 1])
        xl = MathTex(r'\varphi\ (\text{полюс}\angle R),\ ^\circ', font_size=26).next_to(
            ax.x_axis, DOWN, buff=0.38)
        yl = MathTex('F_n', font_size=32).next_to(ax.y_axis, LEFT, buff=0.15).shift(UP * 0.7)
        curve = ax.plot(lambda x: fn_of(x), x_range=[0, 90], color=WHITE, stroke_width=3)
        dot = always_redraw(lambda: Dot(ax.c2p(phi.get_value(), fn_of(phi.get_value())),
                                        radius=0.09, color=POLE))
        fn_lab = MathTex(r'F_n =', font_size=36, color=POLE).move_to(np.array([2.7, -0.1, 0.]))
        fn_num = DecimalRU(0, num_decimal_places=2, font_size=36, color=POLE)
        fn_num.add_updater(lambda m: m.set_value(fn_of(phi.get_value())).next_to(fn_lab, RIGHT, buff=0.15))
        fr_lab = MathTex(r'f_R =', font_size=30, color=GREY_A).move_to(np.array([4.75, -0.1, 0.]))
        fr_num = DecimalRU(0, num_decimal_places=2, font_size=30, color=GREY_A)
        fr_num.add_updater(lambda m: m.set_value(
            np.mean(np.cos(np.radians(phi.get_value() + dev)) ** 2)).next_to(fr_lab, RIGHT, buff=0.12))
        fn_num.update()
        fr_num.update()
        fn_val = VGroup(fn_lab, fn_num)
        fr_val = VGroup(fr_lab, fr_num)

        self.add(wall, grains, r_arrow, r_lab, t_arrow, t_lab, legend,
                 pf, cross, pole, pf_lab, pR, pT, pL,
                 ax, xl, yl, curve, dot, fn_val, fr_val)
        self.wait(1.5)
        self.play(phi.animate.set_value(90.), run_time=5.0, rate_func=smooth)
        self.wait(1.5)
        self.play(phi.animate.set_value(0.), run_time=5.0, rate_func=smooth)
        self.wait(T_LOOP - 13.0 - 0.04)
