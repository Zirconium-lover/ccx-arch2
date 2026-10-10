"""Looping GIF: nanoscale delta-hydride platelets on the alpha-Zr basal plane
stack like a deck of cards into a mesoscale hydride ~14.7 deg off the basal
plane, i.e. along (10-17) (Motta, Capolungo, Chen, JNM 518 (2019)).

Left, 3D: hcp alpha-Zr cell (c/a = 1.593), basal plane, a delta platelet with
(0001)_alpha || {111}_delta.  The cell turns 120 deg per loop - the
3-fold symmetry of the hcp cell - so the loop has no seam.
Right, 2D: section through a grain, basal traces, platelets appearing one by
one, the mesoscale envelope and its angle.

    manim -r 1280,720 --fps 25 deck_of_cards.py DeckOfCards
"""
import numpy as np
from manim import *

config.tex_template = TexTemplate(preamble=r'''\usepackage[T2A]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage[english,russian]{babel}
\usepackage{amsmath}
\usepackage{amssymb}''')

BG = ManimColor('#000000')
ZR = ManimColor('#7FB2D9')        # alpha-Zr lattice
BASAL = ManimColor('#2CB1D9')     # basal plane
HYD = ManimColor('#D8D8D8')       # delta hydride
ACC = YELLOW
FONT = 'Noto Sans'
COA = 1.593
ANG = np.degrees(np.arctan(2 / np.sqrt(3) * COA / 7))   # 14.7 deg, (10-17)
T_LOOP = 14.0


class DeckOfCards(ThreeDScene):
    def construct(self):
        self.camera.background_color = BG
        TH = -30 * DEGREES
        self.set_camera_orientation(phi=68 * DEGREES, theta=TH, zoom=1.0)
        # screen-right direction in world coordinates for this camera; the
        # cell is moved along it so that it sits on the left of the frame
        r_scr = np.array([np.cos(TH + PI / 2), np.sin(TH + PI / 2), 0.])
        shift3d = -3.7 * r_scr + IN * 0.85

        # ---------------- left: the hcp cell ----------------
        a = 1.45
        c = COA * a
        hexa = [np.array([a * np.cos(k * PI / 3), a * np.sin(k * PI / 3), 0.]) for k in range(6)]
        cell = VGroup()
        for z in (0., c):
            cell.add(Polygon(*[p + OUT * z for p in hexa], stroke_color=ZR, stroke_width=2.5,
                             fill_opacity=0))
        for p in hexa:
            cell.add(Line(p, p + OUT * c, stroke_color=ZR, stroke_width=2.5))
        atoms = VGroup()
        for z in (0., c):
            for p in hexa + [ORIGIN]:
                atoms.add(Dot3D(p + OUT * z, radius=0.07, color=ZR))
        for k in (0, 2, 4):                          # the B layer at c/2
            q = (hexa[k] + hexa[k + 1] + ORIGIN) / 3
            atoms.add(Dot3D(q + OUT * c / 2, radius=0.07, color=ZR))
        basal = Polygon(*[p + OUT * c / 2 for p in hexa], stroke_width=0,
                        fill_color=BASAL, fill_opacity=0.28)
        caxis = Arrow3D(OUT * (c + 0.05), OUT * (c + 0.9), color=GREY_A, thickness=0.012)
        disc = Circle(radius=0.75 * a, stroke_color=WHITE, stroke_width=1.5,
                      fill_color=HYD, fill_opacity=0.95).move_to(OUT * (c / 2 + 0.01))
        left = VGroup(cell, atoms, basal, caxis)
        left.shift(shift3d)
        disc.shift(shift3d)

        lab_c = MathTex(r'[0001]', font_size=30, color=GREY_A)
        lab_rel = MathTex(r'(0001)_{\alpha}\parallel\{111\}_{\delta}', font_size=38)
        lab_d = MathTex(r'\delta\text{-ZrH}_{1{,}66}', font_size=38, color=HYD)
        lab_a = MathTex(r'\alpha\text{-Zr}', font_size=38, color=ZR)
        for m in (lab_c, lab_rel, lab_d, lab_a):
            self.add_fixed_in_frame_mobjects(m)
        lab_c.move_to(np.array([-3.05, 2.75, 0]))
        lab_a.to_corner(UL, buff=0.45)
        lab_d.next_to(lab_a, DOWN, buff=0.22, aligned_edge=LEFT)
        lab_rel.move_to(np.array([-3.7, -2.75, 0]))
        self.remove(lab_c, lab_rel, lab_d)

        # ---------------- right: the stack in a section ----------------
        o = np.array([0.55, -1.8, 0.])
        traces = VGroup(*[Line(o + UP * (0.42 * j), o + UP * (0.42 * j) + RIGHT * 5.4,
                               stroke_color=BASAL, stroke_width=1.2, stroke_opacity=0.45)
                          for j in range(10)])
        n = 9
        d = np.array([np.cos(np.radians(ANG)), np.sin(np.radians(ANG)), 0.])
        L = 4.9 / np.cos(np.radians(ANG))
        start = o + RIGHT * 0.25 + UP * 0.35
        plates = VGroup()
        for j in range(n):
            # platelets lie in the basal plane (horizontal); the stack climbs
            # one basal spacing per platelet, so its envelope is inclined
            centre = start + d * (L * (j + 0.5) / n)
            plates.add(RoundedRectangle(width=0.95, height=0.11, corner_radius=0.05,
                                        stroke_width=0, fill_color=HYD, fill_opacity=1)
                       .move_to(centre))
        p0, p1 = start, start + d * L
        nrm = np.array([-d[1], d[0], 0.])
        env = Polygon(p0 - nrm * 0.32, p1 - nrm * 0.32, p1 + nrm * 0.32, p0 + nrm * 0.32,
                      stroke_color=ACC, stroke_width=2.5, fill_opacity=0)
        env = DashedVMobject(env, num_dashes=60)
        base = DashedLine(p0, p0 + RIGHT * 2.6, stroke_color=WHITE, stroke_width=2)
        arc = Arc(radius=2.2, start_angle=0, angle=np.radians(ANG), arc_center=p0,
                  stroke_color=ACC, stroke_width=3)
        ang = MathTex(r'14{,}7^\circ', font_size=36, color=ACC).move_to(
            p0 + RIGHT * 2.0 + DOWN * 0.38)
        t_basal = MathTex(r'\alpha(0001)', font_size=30, color=BASAL).next_to(
            traces[0].get_end(), DOWN, buff=0.12).shift(LEFT * 0.7)
        t_hab = MathTex(r'\alpha(10\bar{1}7)', font_size=32, color=ACC).next_to(
            p1 + nrm * 0.32, UP, buff=0.1)
        t_nano = Text('наногидриды δ', font=FONT, color=HYD).scale(0.4).move_to(
            p0 + RIGHT * 1.1 + UP * 1.05)
        t_meso = Text('мезомасштабный гидрид', font=FONT, color=ACC, weight=BOLD).scale(0.45)
        t_meso.move_to(o + RIGHT * 2.6 + UP * 4.2)
        right = [traces, plates, env, base, arc, ang, t_basal, t_hab, t_nano, t_meso]
        for m in right:
            self.add_fixed_in_frame_mobjects(m)
        self.remove(*right)

        # ---------------- timeline ----------------
        self.add(left, lab_a)
        self.add_fixed_in_frame_mobjects(traces)
        # the cell turns 120 deg per loop: the hcp cell with its B layer has
        # 3-fold symmetry about c, so the last frame matches the first
        rate = (2 * PI / 3) / T_LOOP
        axis_pt = left[2].get_center()
        spin = VGroup(left, disc)
        spin.add_updater(lambda m, dt: m.rotate(rate * dt, axis=OUT, about_point=axis_pt))
        self.add(spin)

        self.play(GrowFromCenter(disc), FadeIn(lab_d), FadeIn(lab_c), run_time=1.2)
        self.play(Write(lab_rel), FadeIn(t_basal), run_time=1.0)
        self.play(LaggedStart(*[FadeIn(p, scale=0.4) for p in plates], lag_ratio=0.55),
                  FadeIn(t_nano), run_time=3.6)
        self.play(Create(env), Create(base), run_time=1.2)
        self.play(Create(arc), Write(ang), FadeIn(t_hab), run_time=1.0)
        self.play(FadeIn(t_meso, shift=DOWN * 0.15), plates.animate.set_fill(ACC, 0.95),
                  run_time=0.8)
        self.wait(2.4)
        self.play(*[FadeOut(m) for m in (plates, env, base, arc, ang, t_hab, t_nano,
                                          t_meso, lab_rel, t_basal, lab_d, lab_c)],
                  disc.animate.scale(0.01), run_time=1.2)
        self.remove(disc)
        used = 1.2 + 1.0 + 3.6 + 1.2 + 1.0 + 0.8 + 2.4 + 1.2
        self.wait(max(0.04, T_LOOP - used))
