"""Looping GIF: the tube texture sets the hydride orientation.

Each grain's basal pole (c axis) sits at phi + spread from the radial
direction R; the hydride lies in the basal plane, i.e. normal to c.  phi
sweeps 0 -> 90 deg -> 0: basal poles along R give tangential hydrides,
along T radial ones.  Right: (0001) pole figure with R at the centre
(stereographic, r = tan(psi/2)) and Fn(phi) = share of hydrides whose normal
is > 45 deg from R (the experimental rule) for the same spread - the step of
Zaimovsky et al.

    manim -r 1600,900 --fps 25 texture.py Texture
"""
import numpy as np
from math import erf
from style import *

SPREAD = 15.0
POLE = ACC


def fn_of(phi_deg):
    s = SPREAD * np.sqrt(2)
    return 0.5 * (1 - erf((45 - phi_deg) / s)) + 0.5 * (1 - erf((45 + phi_deg) / s))


class Texture(Scene):
    def construct(self):
        self.camera.background_color = BG
        rng = np.random.default_rng(7)
        phi = ValueTracker(0.)

        # ---------------- tube wall ----------------
        O = np.array([-3.45, -9.0, 0.])
        R1, R2 = 7.2, 10.25
        A0, A1 = 70.5, 109.5
        wall = AnnularSector(inner_radius=R1, outer_radius=R2, angle=np.radians(A1 - A0),
                             start_angle=np.radians(A0), arc_center=O, fill_color=WALL,
                             fill_opacity=1, stroke_color=GREY_B, stroke_width=3)
        pts = []
        tries = 0
        while len(pts) < 22 and tries < 50000:
            tries += 1
            r = rng.uniform(R1 + 0.4, R2 - 0.4)
            t = np.radians(rng.uniform(A0 + 2.5, A1 - 2.5))
            p = O + r * np.array([np.cos(t), np.sin(t), 0.])
            if all(np.linalg.norm(p - q) > 0.85 for q in pts):
                pts.append(p)
        dev = rng.normal(0, SPREAD, len(pts))
        # the highlighted grain: the one nearest the wall's middle, deviation 0
        star = int(np.argmin([np.linalg.norm(p - (O + (R1 + R2) / 2 * UP)) for p in pts]))
        dev[star] = 0.

        def grain_geom(i):
            p = pts[i]
            rhat = (p - O) / np.linalg.norm(p - O)
            psi = np.arctan2(rhat[1], rhat[0]) + np.radians(phi.get_value() + dev[i])
            c = np.array([np.cos(psi), np.sin(psi), 0.])
            return p, c, np.array([-c[1], c[0], 0.]), np.sin(np.radians(phi.get_value() + dev[i])) ** 2

        def grains():
            g = VGroup()
            for i in range(len(pts)):
                p, c, d, w = grain_geom(i)
                L = 0.34 if i != star else 0.46
                g.add(Line(p - c * 0.27, p + c * 0.27, stroke_width=3, color=POLE),
                      Line(p - d * L, p + d * L, stroke_width=10 if i != star else 13,
                           color=interpolate_color(TAN, RAD, w)))
            p, c, d, w = grain_geom(star)
            g.add(Circle(radius=0.68, stroke_color=WHITE, stroke_width=2.5).move_to(p))
            return g
        gr = always_redraw(grains)

        ps = pts[star]
        lab_c = txt('ось c', TX, POLE).move_to(np.array([ps[0] - 1.6, 2.2, 0.]))
        lab_h = txt('гидрид', TX, WHITE).move_to(np.array([ps[0] + 1.7, 2.2, 0.]))
        lead = always_redraw(lambda: VGroup(
            Line(lab_c.get_bottom() + DOWN * 0.05, grain_geom(star)[0] + grain_geom(star)[1] * 0.3,
                 stroke_width=2, color=POLE),
            Line(lab_h.get_bottom() + DOWN * 0.05, grain_geom(star)[0] + grain_geom(star)[2] * 0.42,
                 stroke_width=2, color=WHITE)))

        top = O + R2 * UP
        ax0 = np.array([-6.55, 1.35, 0.])
        r_ar = Arrow(ax0, ax0 + UP * 0.95, buff=0, color=WHITE, stroke_width=6)
        r_l = tex('R', FS_SM).next_to(r_ar, RIGHT, buff=0.12).shift(UP * 0.25)
        t_ar = Arrow(ax0, ax0 + RIGHT * 0.95, buff=0, color=WHITE, stroke_width=6)
        t_l = tex('T', FS_SM).next_to(t_ar, RIGHT, buff=0.12)

        def caption():
            ph = phi.get_value()
            if ph < 45:
                s1, s2, col, op = 'полюсы вдоль R', '→ тангенциальные гидриды', TAN, np.clip((25 - ph) / 15, 0, 1)
            else:
                s1, s2, col, op = 'полюсы вдоль T', '→ радиальные гидриды', RAD, np.clip((ph - 65) / 15, 0, 1)
            return VGroup(txt(s1, TX + 3, WHITE), txt(s2, TX + 3, col)).arrange(
                DOWN, buff=0.18).move_to(np.array([-3.45, -3.2, 0.])).set_opacity(op)
        cap = always_redraw(caption)

        # ---------------- pole figure + numbers ----------------
        pc = np.array([1.85, 1.3, 0.])
        Rp = 1.3
        pf = Circle(radius=Rp, stroke_color=GREY_B, stroke_width=3).move_to(pc)
        cross = VGroup(Line(pc + LEFT * Rp, pc + RIGHT * Rp, stroke_width=1.5, color=GREY_D),
                       Line(pc + DOWN * Rp, pc + UP * Rp, stroke_width=1.5, color=GREY_D))
        pf_lab = tex('(0001)', FS_SM).move_to(pc + np.array([-Rp - 0.25, Rp - 0.05, 0.]))
        pT = tex('T', FS_SM - 6, GREY_A).move_to(pc + RIGHT * (Rp + 0.25) + DOWN * 0.38)
        pL = tex('L', FS_SM - 6, GREY_A).next_to(pf, UP, buff=0.08)

        def blobs():
            x = np.tan(np.radians(phi.get_value()) / 2) * Rp
            g = VGroup()
            for sgn in (-1, 1):
                for rr, op in ((0.42, 0.18), (0.28, 0.38), (0.15, 0.95)):
                    g.add(Ellipse(width=2 * rr, height=2 * rr * 1.25, stroke_width=0,
                                  fill_color=POLE, fill_opacity=op).move_to(pc + RIGHT * sgn * x))
            return g
        pole = always_redraw(blobs)
        pR = tex('R', FS_SM - 6, GREY_A).move_to(pc + DOWN * 0.32 + RIGHT * 0.3)

        fn_lab = tex(r'F_n =', FS_BIG).move_to(np.array([4.85, 1.25, 0.]))
        fn_num = DecimalRU(0, num_decimal_places=2, font_size=FS_BIG, color=POLE)
        fn_num.add_updater(lambda m: m.set_value(fn_of(phi.get_value())).next_to(
            fn_lab, RIGHT, buff=0.22))
        fn_num.update()

        # ---------------- Fn(phi) ----------------
        ax = Axes(x_range=[0, 90, 15], y_range=[0, 1, 0.5], x_length=5.3, y_length=2.15,
                  tips=False, axis_config=dict(color=GREY_B, stroke_width=3, include_ticks=True,
                                               font_size=34)).move_to(np.array([3.75, -1.55, 0.]))
        ax.add_coordinates([0, 45, 90], [])
        ylab = VGroup(tex('0', 34).next_to(ax.c2p(0, 0), LEFT, buff=0.15),
                      tex('1', 34).next_to(ax.c2p(0, 1), LEFT, buff=0.15))
        xl = txt('угол полюса от R, °', TX_SM, GREY_A).next_to(ax.x_axis, DOWN, buff=0.62)
        yl = tex('F_n', FS_SM).next_to(ax.y_axis, LEFT, buff=0.3).shift(UP * 0.4)
        curve = ax.plot(lambda x: fn_of(x), x_range=[0, 90], color=WHITE, stroke_width=5)
        dot = always_redraw(lambda: Dot(ax.c2p(phi.get_value(), fn_of(phi.get_value())),
                                        radius=0.13, color=POLE))

        self.add(wall, gr, lead, lab_c, lab_h, r_ar, r_l, t_ar, t_l, cap,
                 pf, cross, pole, pf_lab, pT, pL, pR, fn_lab, fn_num,
                 ax, ylab, xl, yl, curve, dot)
        self.wait(2.0)
        self.play(phi.animate.set_value(90.), run_time=5.0, rate_func=smooth)
        self.wait(2.0)
        self.play(phi.animate.set_value(0.), run_time=5.0, rate_func=smooth)
        self.wait(0.04)
