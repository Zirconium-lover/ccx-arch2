"""Looping GIF: nanoscale delta-hydride platelets on the alpha-Zr basal plane
stack like a deck of cards into a mesoscale hydride ~14.7 deg off the basal
plane, i.e. along (10-17) (Motta, Capolungo, Chen, JNM 518 (2019)).

Left (3D): hcp alpha-Zr cell (c/a = 1.593) with a delta platelet on the
basal plane, (0001)_alpha || {111}_delta; the cell turns 120 deg per loop,
its own 3-fold symmetry, so the loop has no seam.
Right: a grain with a mesoscale hydride lying 14.7 deg off the basal traces;
a zoom shows it is a stack of nanoscale platelets, each in the basal plane.

    manim -r 1600,900 --fps 25 deck_of_cards.py DeckOfCards
"""
import numpy as np
from style import *

COA = 1.593
ANG = np.degrees(np.arctan(2 / np.sqrt(3) * COA / 7))   # 14.7 deg, (10-17)
BASAL = TAN
T_LOOP = 14.0


class DeckOfCards(ThreeDScene):
    def fixed(self, *ms):
        self.add_fixed_in_frame_mobjects(*ms)
        self.remove(*ms)
        return ms

    def construct(self):
        self.camera.background_color = BG
        TH = -30 * DEGREES
        PHI = 68 * DEGREES
        self.set_camera_orientation(phi=PHI, theta=TH, zoom=1.0)
        r_scr = np.array([np.cos(TH + PI / 2), np.sin(TH + PI / 2), 0.])

        # ---------------- left: the hcp cell ----------------
        a = 1.8
        c = COA * a
        shift3d = -4.05 * r_scr + IN * (c / 2 + 0.25 / np.sin(PHI))
        hexa = [np.array([a * np.cos(k * PI / 3), a * np.sin(k * PI / 3), 0.]) for k in range(6)]
        cell = VGroup()
        for z in (0., c):
            cell.add(Polygon(*[p + OUT * z for p in hexa], stroke_color=ZR, stroke_width=4,
                             fill_opacity=0))
        for p in hexa:
            cell.add(Line(p, p + OUT * c, stroke_color=ZR, stroke_width=4))
        atoms = VGroup()
        for z in (0., c):
            for p in hexa + [ORIGIN]:
                atoms.add(Dot3D(p + OUT * z, radius=0.1, color=ZR))
        for k in (0, 2, 4):
            q = (hexa[k] + hexa[k + 1] + ORIGIN) / 3
            atoms.add(Dot3D(q + OUT * c / 2, radius=0.1, color=ZR))
        basal = Polygon(*[p + OUT * c / 2 for p in hexa], stroke_width=0,
                        fill_color=BASAL, fill_opacity=0.28)
        caxis = Arrow3D(OUT * (c + 0.05), OUT * (c + 1.0), color=GREY_A, thickness=0.016)
        disc = Circle(radius=0.78 * a, stroke_color=WHITE, stroke_width=2,
                      fill_color=HYD, fill_opacity=0.95).move_to(OUT * (c / 2 + 0.01))
        left = VGroup(cell, atoms, basal, caxis)
        left.shift(shift3d)
        disc.shift(shift3d)

        lab_a = self.fixed(tex(r'\alpha\text{-Zr}', FS, ZR))[0]
        lab_a.move_to(np.array([-6.2, 2.45, 0.]))
        lab_c = self.fixed(tex('c', FS_SM, GREY_A))[0]
        lab_c.move_to(np.array([-3.75, 2.2, 0.]))
        lab_d = self.fixed(tex(r'\delta\text{-ZrH}_{1{,}66}', FS_SM, HYD))[0]
        lab_d.move_to(np.array([-0.85, 0.8, 0.]))
        lead_d = self.fixed(Line(np.array([-1.55, 0.5, 0.]), np.array([-2.55, -0.2, 0.]),
                                 stroke_width=2.5, color=HYD))[0]
        lab_rel = self.fixed(tex(r'(0001)_{\alpha}\parallel\{111\}_{\delta}', FS))[0]
        lab_rel.move_to(np.array([-3.95, -3.35, 0.]))

        # ---------------- right top: grain with a mesoscale hydride ----------------
        g0, g1 = np.array([0.45, 0.45, 0.]), np.array([6.85, 2.85, 0.])
        grain = RoundedRectangle(width=g1[0] - g0[0], height=g1[1] - g0[1], corner_radius=0.2,
                                 stroke_color=GREY_B, stroke_width=2.5, fill_color=WALL,
                                 fill_opacity=1).move_to((g0 + g1) / 2)
        gtr = VGroup(*[Line(np.array([g0[0] + 0.15, y, 0.]), np.array([g1[0] - 0.15, y, 0.]),
                            stroke_width=1.5, color=BASAL, stroke_opacity=0.4)
                       for y in np.arange(g0[1] + 0.3, g1[1] - 0.1, 0.32)])
        d = np.array([np.cos(np.radians(ANG)), np.sin(np.radians(ANG)), 0.])
        m0 = np.array([1.05, 0.95, 0.])
        meso = Line(m0, m0 + d * 5.55, stroke_width=11, color=ACC)
        zc = m0 + d * 2.3
        zoom = Square(side_length=0.75, stroke_color=WHITE, stroke_width=3).move_to(zc)
        t_meso = txt('мезомасштабный гидрид', TX, ACC).move_to(np.array([3.05, 2.45, 0.]))

        # ---------------- right bottom: the zoom - a deck of cards ----------------
        z0, z1 = np.array([0.45, -3.85, 0.]), np.array([6.85, 0.05, 0.])
        zbox = Rectangle(width=z1[0] - z0[0], height=z1[1] - z0[1], stroke_color=WHITE,
                         stroke_width=3).move_to((z0 + z1) / 2)
        zl1 = DashedLine(zoom.get_corner(DL), zbox.get_corner(UL), stroke_width=2, color=GREY_B)
        zl2 = DashedLine(zoom.get_corner(DR), zbox.get_corner(UR), stroke_width=2, color=GREY_B)
        ztr = VGroup(*[Line(np.array([z0[0] + 0.15, y, 0.]), np.array([z1[0] - 0.15, y, 0.]),
                            stroke_width=1.5, color=BASAL, stroke_opacity=0.45)
                       for y in np.arange(z0[1] + 0.35, z1[1] - 0.1, 0.42)])
        n = 8
        L = 5.4 / np.cos(np.radians(ANG))
        p0 = np.array([0.95, -3.05, 0.])
        plates = VGroup(*[
            RoundedRectangle(width=1.3, height=0.17, corner_radius=0.085, stroke_width=0,
                             fill_color=HYD, fill_opacity=1).move_to(p0 + d * (L * (j + 0.5) / n))
            for j in range(n)])
        nrm = np.array([-d[1], d[0], 0.])
        p1 = p0 + d * L
        env = DashedVMobject(Polygon(p0 - nrm * 0.42, p1 - nrm * 0.42, p1 + nrm * 0.42,
                                     p0 + nrm * 0.42, stroke_color=ACC, stroke_width=3.5),
                             num_dashes=70)
        base = DashedLine(p0, p0 + RIGHT * 3.2, stroke_color=WHITE, stroke_width=2.5)
        arc = Arc(radius=2.9, start_angle=0, angle=np.radians(ANG), arc_center=p0,
                  stroke_color=ACC, stroke_width=4)
        t_ang = tex(r'14{,}7^\circ', FS, ACC).move_to(p0 + RIGHT * 2.35 + DOWN * 0.45)
        t_hab = tex(r'(10\bar{1}7)', FS_SM, ACC).next_to(p1 + nrm * 0.42, UP, buff=0.12).shift(LEFT * 0.4)
        t_bas = tex(r'(0001)', FS_SM - 4, BASAL).move_to(np.array([5.95, -3.5, 0.]))
        t_nano = txt('наногидриды δ', TX, HYD).move_to(np.array([2.1, -1.0, 0.]))

        static_r = [grain, gtr, zbox, ztr, zl1, zl2, t_bas]
        dyn_r = [meso, zoom, t_meso, plates, env, base, arc, t_ang, t_hab, t_nano]
        self.fixed(*static_r, *dyn_r)

        # ---------------- timeline ----------------
        rate = (2 * PI / 3) / T_LOOP
        axis_pt = left[2].get_center()
        spin = VGroup(left, disc)
        spin.add_updater(lambda m, dt: m.rotate(rate * dt, axis=OUT, about_point=axis_pt))
        self.add(spin)
        self.remove(disc)
        self.add(lab_a, lab_c, *static_r)

        self.play(GrowFromCenter(disc), FadeIn(lab_d), Create(lead_d), run_time=1.2)
        self.play(Write(lab_rel), run_time=0.9)
        self.play(LaggedStart(*[FadeIn(p, scale=0.4) for p in plates], lag_ratio=0.55),
                  FadeIn(t_nano), run_time=3.4)
        self.play(Create(env), Create(base), run_time=1.0)
        self.play(Create(arc), Write(t_ang), FadeIn(t_hab), run_time=1.0)
        self.play(Create(meso), FadeIn(zoom), FadeIn(t_meso, shift=DOWN * 0.1),
                  plates.animate.set_fill(ACC, 1), run_time=1.1)
        self.wait(2.6)
        self.play(*[FadeOut(m) for m in (meso, zoom, t_meso, plates, env, base, arc, t_ang,
                                          t_hab, t_nano, lab_rel, lab_d, lead_d)],
                  disc.animate.scale(0.01), run_time=1.2)
        self.remove(disc)
        used = 1.2 + 0.9 + 3.4 + 1.0 + 1.0 + 1.1 + 2.6 + 1.2
        self.wait(max(0.04, T_LOOP - used))
