"""Looping GIF: nanoscale delta-hydride platelets on the alpha-Zr basal plane
stack like a deck of cards into a mesoscale hydride ~14.7 deg off the basal
plane, i.e. along (10-17) (Motta, Capolungo, Chen, JNM 518 (2019)).

Left (3D): hcp alpha-Zr cell (c/a = 1.593) with a delta platelet on the
basal plane, (0001)_alpha || {111}_delta; the cell turns 120 deg per loop,
its own 3-fold symmetry, so the loop has no seam.
Right: a grain with a mesoscale hydride lying 14.7 deg off the basal traces;
a zoom shows it is a stack of nanoscale platelets, each in the basal plane.

Layout from measurement: the 3D cell's size and position come from the
camera projection of its vertices over the whole 120 deg turn; the 2D labels
are checked against that projected outline and against the hydride band.

    manim -r 1600,900 --fps 25 deck_of_cards.py DeckOfCards
"""
from style import *

COA = 1.593
ANG = np.degrees(np.arctan(2 / np.sqrt(3) * COA / 7))   # 14.7 deg, (10-17)
BASAL = TAN
T_LOOP = 14.0
TH, PHI = -30 * DEGREES, 68 * DEGREES


def clear_of(m, pts, pad=0.08):
    """True if no point of pts lies in m's bounding box grown by pad."""
    b = bbox(m)
    return not np.any((pts[:, 0] > b[0] - pad) & (pts[:, 0] < b[1] + pad) &
                      (pts[:, 1] > b[2] - pad) & (pts[:, 1] < b[3] + pad))


def band_pts(p0, p1, half, n=200):
    """Sample points filling the band p0-p1 of half-width half."""
    d = (p1 - p0) / np.linalg.norm(p1 - p0)
    nrm = np.array([-d[1], d[0], 0.])
    s, t = np.meshgrid(np.linspace(0, 1, n), np.linspace(-1, 1, 9))
    return (p0[None, :] + s.ravel()[:, None] * (p1 - p0)[None, :]
            + t.ravel()[:, None] * half * nrm[None, :])


class DeckOfCards(ThreeDScene):
    def fixed(self, *ms):
        self.add_fixed_in_frame_mobjects(*ms)
        self.remove(*ms)
        return ms

    def construct(self):
        self.camera.background_color = BG
        self.set_camera_orientation(phi=PHI, theta=TH, zoom=1.0)
        cam = self.camera
        cam.reset_rotation_matrix()          # else it still holds the default view
        Rm = cam.get_rotation_matrix()

        # ---------------- fixed-size labels ----------------
        lab_a = tex(r'\alpha\text{-Zr}', L2M, ZR)
        lab_d = tex(r'\delta\text{-ZrH}_{1{,}66}', L2M, HYD)
        lab_c = tex('c', L2M, GREY_A)
        lab_rel = tex(r'(0001)_{\alpha}\parallel\{111\}_{\delta}', L2M)
        t_meso = VGroup(txt('мезомасштабный', L2, ACC), txt('гидрид', L2, ACC)).arrange(
            DOWN, aligned_edge=LEFT, buff=0.14)
        t_nano = txt('наногидриды', L2, HYD)
        t_ang = tex(r'14{,}7^\circ', L2M, ACC)
        t_hab = tex(r'(10\bar{1}7)', L2M, ACC)
        t_bas = tex(r'(0001)', L2M, BASAL)

        # ---------------- grid ----------------
        side_w = max(lab_a.width, lab_d.width)
        wL = 6.9
        cellL = Box(CONTENT.x0, CONTENT.x0 + wL, CONTENT.y0, CONTENT.y1)
        cellR = Box(cellL.x1 + 0.5, CONTENT.x1, CONTENT.y0, CONTENT.y1)
        rCell = Box(cellL.x0, cellL.x1 - side_w - 0.35, cellL.y0 + lab_rel.height + 0.35, cellL.y1)
        rG, rZ = cellR.stack(2.7, cellR.h - 2.7 - 0.5, gap=0.5)

        # ---------------- the hcp cell: measure its projection ----------------
        hexa0 = [np.array([np.cos(k * PI / 3), np.sin(k * PI / 3), 0.]) for k in range(6)]

        def outline(a, shift, ang):
            c = COA * a
            rot = np.array([[np.cos(ang), -np.sin(ang), 0], [np.sin(ang), np.cos(ang), 0], [0, 0, 1]])
            v = [rot @ (a * p) + OUT * z for p in hexa0 for z in (-c / 2, c / 2)]
            v.append(OUT * (c / 2 + 0.6 * a))                 # c-axis tip
            return cam.project_points(np.array(v) + shift)[:, :2]

        def extent(a, shift):
            pts = np.vstack([outline(a, shift, g) for g in np.radians(np.arange(0, 121, 4))])
            return pts, pts.min(0), pts.max(0)
        a, shift = 1., ORIGIN
        for _ in range(6):                  # perspective is not linear in a: iterate
            _, lo, hi = extent(a, shift)
            a *= 0.98 * min(rCell.w / (hi - lo)[0], rCell.h / (hi - lo)[1])
            _, lo, hi = extent(a, shift)
            dsc = rCell.c[:2] - (lo + hi) / 2
            shift = shift + Rm.T @ np.array([dsc[0], dsc[1], 0.])
        cell_pts, lo, hi = extent(a, shift)
        print('LAYOUT cell a=%.3f  projected x %.2f..%.2f y %.2f..%.2f  (box %s)'
              % (a, lo[0], hi[0], lo[1], hi[1], rCell))

        c = COA * a
        hexa = [a * p for p in hexa0]
        cell = VGroup()
        for z in (-c / 2, c / 2):
            cell.add(Polygon(*[p + OUT * z for p in hexa], stroke_color=ZR, stroke_width=4,
                             fill_opacity=0))
        for p in hexa:
            cell.add(Line(p + IN * c / 2, p + OUT * c / 2, stroke_color=ZR, stroke_width=4))
        atoms = VGroup()
        for z in (-c / 2, c / 2):
            for p in hexa + [ORIGIN]:
                atoms.add(Dot3D(p + OUT * z, radius=0.055 * a, color=ZR))
        for k in (0, 2, 4):
            atoms.add(Dot3D((hexa[k] + hexa[k + 1]) / 3, radius=0.055 * a, color=ZR))
        basal = Polygon(*hexa, stroke_width=0, fill_color=BASAL, fill_opacity=0.28)
        caxis = Arrow3D(OUT * (c / 2 + 0.03 * a), OUT * (c / 2 + 0.6 * a), color=GREY_A,
                        thickness=0.009 * a)
        disc = Circle(radius=0.78 * a, stroke_color=WHITE, stroke_width=2,
                      fill_color=HYD, fill_opacity=0.95).move_to(OUT * 0.01)
        left = VGroup(cell, atoms, basal, caxis).shift(shift)
        disc.shift(shift)

        # labels around the cell, from the projected points
        tip = cam.project_points(np.array([shift + OUT * (c / 2 + 0.6 * a)]))[0]
        rim = cam.project_points(np.array([shift + 0.78 * a * np.array([np.cos(t), np.sin(t), 0.])
                                           for t in np.linspace(0, 2 * PI, 72)]))
        rim_r = rim[np.argmax(rim[:, 0])]
        dc = cam.project_points(np.array([shift]))[0]
        lab_c.next_to(np.array([tip[0], tip[1], 0.]), RIGHT, buff=0.15).align_to(
            np.array([0., tip[1], 0.]), UP)
        lab_d.move_to(np.array([cellL.x1, dc[1], 0.]), aligned_edge=RIGHT)
        lab_a.move_to(np.array([cellL.x1, rCell.y1, 0.]), aligned_edge=UR)
        lead_d = Line(lab_d.get_left() + LEFT * 0.12, np.array([rim_r[0] - 0.25, rim_r[1], 0.]),
                      stroke_width=2.5, color=HYD)
        lab_rel.move_to(np.array([rCell.c[0], cellL.y0, 0.]), aligned_edge=DOWN)

        # ---------------- right top: grain with a mesoscale hydride ----------------
        d = np.array([np.cos(np.radians(ANG)), np.sin(np.radians(ANG)), 0.])
        grain = RoundedRectangle(width=rG.w, height=rG.h, corner_radius=0.2, stroke_color=GREY_B,
                                 stroke_width=2.5, fill_color=WALL, fill_opacity=1).move_to(rG.c)
        gtr = VGroup(*[Line(np.array([rG.x0 + 0.15, y, 0.]), np.array([rG.x1 - 0.15, y, 0.]),
                            stroke_width=1.5, color=BASAL, stroke_opacity=0.4)
                       for y in np.arange(rG.y0 + 0.3, rG.y1 - 0.1, 0.32)])
        m0 = np.array([rG.x0 + 0.5, rG.y0 + 0.35, 0.])
        m1 = m0 + d * (rG.w - 1.0) / np.cos(np.radians(ANG))
        meso = Line(m0, m1, stroke_width=11, color=ACC)
        t_meso.move_to(np.array([rG.x0 + 0.3, rG.y1 - 0.2, 0.]), aligned_edge=UL)
        zc = m0 + (m1 - m0) * 0.38
        zoom = Square(side_length=0.75, stroke_color=WHITE, stroke_width=3).move_to(zc)
        zl1 = DashedLine(zoom.get_corner(DL), np.array([rZ.x0, rZ.y1, 0.]), stroke_width=2,
                         color=GREY_B)
        zl2 = DashedLine(zoom.get_corner(DR), np.array([rZ.x1, rZ.y1, 0.]), stroke_width=2,
                         color=GREY_B)

        # ---------------- right bottom: the zoom - a deck of cards ----------------
        zbox = Rectangle(width=rZ.w, height=rZ.h, stroke_color=WHITE, stroke_width=3).move_to(rZ.c)
        ztr = VGroup(*[Line(np.array([rZ.x0 + 0.15, y, 0.]), np.array([rZ.x1 - 0.15, y, 0.]),
                            stroke_width=1.5, color=BASAL, stroke_opacity=0.45)
                       for y in np.arange(rZ.y0 + 0.35, rZ.y1 - 0.1, 0.42)])
        n, pw, half = 8, 1.3, 0.42
        p0 = np.array([rZ.x0 + 0.5, rZ.y0 + 0.85, 0.])
        L = (rZ.w - 1.0 - pw / 2) / np.cos(np.radians(ANG))
        p1 = p0 + d * L
        plates = VGroup(*[
            RoundedRectangle(width=pw, height=0.17, corner_radius=0.085, stroke_width=0,
                             fill_color=HYD, fill_opacity=1).move_to(p0 + d * (L * (j + 0.5) / n))
            for j in range(n)])
        nrm = np.array([-d[1], d[0], 0.])
        env = DashedVMobject(Polygon(p0 - nrm * half, p1 - nrm * half, p1 + nrm * half,
                                     p0 + nrm * half, stroke_color=ACC, stroke_width=3.5),
                             num_dashes=70)
        base = DashedLine(p0, p0 + RIGHT * 3.2, stroke_color=WHITE, stroke_width=2.5)
        arc = Arc(radius=2.9, start_angle=0, angle=np.radians(ANG), arc_center=p0,
                  stroke_color=ACC, stroke_width=4)
        t_ang.next_to(p0 + RIGHT * 2.4, DOWN, buff=0.18)
        t_hab.next_to(p1 + nrm * half, UP, buff=0.15).align_to(np.array([rZ.x1 - 0.3, 0., 0.]), RIGHT)
        t_nano.move_to(np.array([rZ.x0 + 0.3, rZ.y1 - 0.25, 0.]), aligned_edge=UL)
        t_bas.move_to(np.array([rZ.x1 - 0.3, rZ.y0 + 0.2, 0.]), aligned_edge=DR)

        # ---------------- checks ----------------
        report({'cell': VGroup(*[Dot(np.array([x, y, 0.]), radius=0.01)
                                 for x, y in ((lo[0], lo[1]), (hi[0], hi[1]))]),
                'lab_a': lab_a, 'lab_d': lab_d, 'lab_rel': lab_rel,
                'grain': grain, 'zbox': zbox})
        cp = np.hstack([cell_pts, np.zeros((len(cell_pts), 1))])
        hull_band = []
        for i in range(0, len(cp), 13):                     # each angle: 12 vertices + tip
            v = cp[i:i + 12]
            for j in range(6):                              # vertical edges and both hexagons
                for q0, q1 in ((v[2 * j], v[2 * j + 1]), (v[2 * j], v[(2 * j + 2) % 12]),
                               (v[2 * j + 1], v[(2 * j + 3) % 12])):
                    hull_band.append(np.linspace(q0, q1, 12))
        edge_pts = np.vstack(hull_band)
        deck_pts = np.vstack([band_pts(p0, p1, half), band_pts(p0 - RIGHT * pw / 2, p1 + RIGHT * pw / 2, 0.1)])
        meso_pts = band_pts(m0, m1, 0.06)
        for name, m, pts in (('lab_a', lab_a, edge_pts), ('lab_d', lab_d, edge_pts),
                             ('lab_c', lab_c, edge_pts), ('lab_rel', lab_rel, edge_pts),
                             ('t_meso', t_meso, meso_pts), ('t_nano', t_nano, deck_pts),
                             ('t_hab', t_hab, deck_pts), ('t_bas', t_bas, deck_pts),
                             ('t_ang', t_ang, deck_pts)):
            print('LAYOUT %-8s %s' % (name, 'clear' if clear_of(m, pts) else 'HITS the drawing'))
        for name, m, box in (('t_meso', t_meso, rG), ('t_nano', t_nano, rZ), ('t_hab', t_hab, rZ),
                             ('t_bas', t_bas, rZ), ('t_ang', t_ang, rZ)):
            b = bbox(m)
            ok = b[0] >= box.x0 and b[1] <= box.x1 and b[2] >= box.y0 and b[3] <= box.y1
            print('LAYOUT %-8s %s' % (name, 'inside its panel' if ok else 'OUTSIDE its panel'))

        static_r = [grain, gtr, zbox, ztr, t_bas]
        dyn_r = [meso, zoom, zl1, zl2, t_meso, plates, env, base, arc, t_ang, t_hab, t_nano]
        self.fixed(lab_a, lab_c, lab_d, lead_d, lab_rel, *static_r, *dyn_r)

        # ---------------- timeline ----------------
        rate = (2 * PI / 3) / T_LOOP
        spin = VGroup(left, disc)
        spin.add_updater(lambda m, dt: m.rotate(rate * dt, axis=OUT, about_point=shift))
        self.add(spin)
        self.remove(disc)
        self.add(lab_a, lab_c, *static_r)

        self.play(GrowFromCenter(disc), FadeIn(lab_d), Create(lead_d), run_time=1.2)
        self.play(Write(lab_rel), run_time=0.9)
        self.play(LaggedStart(*[FadeIn(p, scale=0.4) for p in plates], lag_ratio=0.55),
                  FadeIn(t_nano), run_time=3.4)
        self.play(Create(env), Create(base), run_time=1.0)
        self.play(Create(arc), Write(t_ang), FadeIn(t_hab), run_time=1.0)
        self.play(Create(meso), FadeIn(zoom), Create(zl1), Create(zl2),
                  FadeIn(t_meso, shift=DOWN * 0.1), plates.animate.set_fill(ACC, 1), run_time=1.1)
        self.wait(2.6)
        self.play(*[FadeOut(m) for m in (meso, zoom, zl1, zl2, t_meso, plates, env, base, arc,
                                          t_ang, t_hab, t_nano, lab_rel, lab_d, lead_d)],
                  disc.animate.scale(0.01), run_time=1.2)
        self.remove(disc)
        used = 1.2 + 0.9 + 3.4 + 1.0 + 1.0 + 1.1 + 2.6 + 1.2
        self.wait(max(0.04, T_LOOP - used))
