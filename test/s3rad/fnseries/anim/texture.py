"""Looping GIF: the tube texture sets the hydride orientation.

Each grain's basal pole (c axis) sits at phi + spread from the radial
direction R; the hydride lies in the basal plane, i.e. normal to c.  phi
sweeps 0 -> 90 deg -> 0: the hydrides go from tangential to radial and back.
Middle: a zoom on the circled grain - its hcp cell (as on the deck-of-cards
slide) turns with the grain: the c axis is the crystal's, the hydride lies in
the basal plane across it.
Right: the (0001) pole figure, R at the centre (stereographic, r = tan(psi/2)).
Qualitative: no numbers.

Layout: CONTENT = [wall cell: axes + labels band / wall / caption | zoom
column | pole figure column]; the wall is built in nominal units and mapped
into its cell by the measured scale k.

    manim -r 1600,900 --fps 25 texture.py Texture
"""
import os
from style import *

SPREAD = 15.0
POLE = ACC
ZR = ManimColor('#7FB2D9')     # alpha-Zr lattice, as on the deck-of-cards slide
COA = 1.593                    # c / a of alpha-Zr
TILT = np.radians(24)          # the zoomed cell is seen a little from above its basal plane
CAP = 32                       # caption Text size (the wall cell is narrower with the zoom)


class Texture(Scene):
    def construct(self):
        self.camera.background_color = BG
        rng = np.random.default_rng(7)
        phi = ValueTracker(0.)

        # ---------------- measure fixed-size labels ----------------
        caps = [txt('тангенциальные гидриды', CAP, TAN), txt('радиальные гидриды', CAP, RAD)]
        cap_h = max(m.height for m in caps)
        lab_c = txt('ось c', L2, POLE)
        lab_h = txt('гидрид', L2, WHITE)
        r_l, t_l = tex('R', L2M), tex('T', L2M)
        AR = 0.8                                   # axis arrow length
        band_h = max(AR + 0.1, lab_c.height)
        pT, pL, pR = tex('T', L2M, GREY_A), tex('L', L2M, GREY_A), tex('R', L2M, GREY_A)
        pf_lab = tex('(0001)', L2M)

        lab_z = txt('ячейка α-Zr', L2, ZR)

        # ---------------- grid ----------------
        pole_w, zoom_w, gp = 3.7, 2.75, 0.35
        cellW = Box(CONTENT.x0 + 0.03, CONTENT.x1 - pole_w - zoom_w - 2 * gp, CONTENT.y0, CONTENT.y1)
        cellZ = Box(cellW.x1 + gp, cellW.x1 + gp + zoom_w, CONTENT.y0, CONTENT.y1)
        cellP = Box(CONTENT.x1 - pole_w, CONTENT.x1, CONTENT.y0, CONTENT.y1)

        # nominal wall: centre O at the origin
        R1, R2, A0, A1 = 7.2, 10.25, 70.5, 109.5
        half = np.radians((A1 - A0) / 2)
        nb = Box(-R2 * np.sin(half), R2 * np.sin(half), R1 * np.cos(half), R2)
        gap1, gap2 = 0.3, 0.35
        avail_h = cellW.h - band_h - cap_h - gap1 - gap2
        k = min(cellW.w / nb.w, avail_h / nb.h)
        used = band_h + gap1 + nb.h * k + gap2 + cap_h
        top = cellW.y1 - (cellW.h - used) / 2
        rBand = Box(cellW.x0, cellW.x1, top - band_h, top)
        rWall = Box(cellW.c[0] - nb.w * k / 2, cellW.c[0] + nb.w * k / 2,
                    rBand.y0 - gap1 - nb.h * k, rBand.y0 - gap1)
        rCap = Box(cellW.x0, cellW.x1, rWall.y0 - gap2 - cap_h, rWall.y0 - gap2)
        O = rWall.c + k * (np.zeros(3) - nb.c)

        # ---------------- tube wall ----------------
        wall = AnnularSector(inner_radius=R1 * k, outer_radius=R2 * k, angle=np.radians(A1 - A0),
                             start_angle=np.radians(A0), arc_center=O, fill_color=WALL,
                             fill_opacity=1, stroke_color=GREY_B, stroke_width=3)
        pts, tries = [], 0
        while len(pts) < 22 and tries < 50000:
            tries += 1
            r = rng.uniform(R1 + 0.4, R2 - 0.4) * k
            t = np.radians(rng.uniform(A0 + 2.5, A1 - 2.5))
            p = O + r * np.array([np.cos(t), np.sin(t), 0.])
            if all(np.linalg.norm(p - q) > 0.85 * k for q in pts):
                pts.append(p)
        dev = rng.normal(0, SPREAD, len(pts))
        star = int(np.argmin([np.linalg.norm(p - (O + (R1 + R2) / 2 * k * UP)) for p in pts]))
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
                L = (0.34 if i != star else 0.46) * k
                g.add(Line(p - c * 0.27 * k, p + c * 0.27 * k, stroke_width=3, color=POLE),
                      Line(p - d * L, p + d * L, stroke_width=10 if i != star else 13,
                           color=interpolate_color(TAN, RAD, w)))
            g.add(Circle(radius=0.68 * k, stroke_color=WHITE, stroke_width=2.5).move_to(pts[star]))
            return g
        gr = always_redraw(grains)

        # band: axes on the left, the two labels over the highlighted grain
        ax0 = np.array([rBand.x0 + 0.1, rBand.y0 + 0.05, 0.])
        r_ar = Arrow(ax0, ax0 + UP * AR, buff=0, color=WHITE, stroke_width=6,
                     max_tip_length_to_length_ratio=0.3)
        t_ar = Arrow(ax0, ax0 + RIGHT * AR, buff=0, color=WHITE, stroke_width=6,
                     max_tip_length_to_length_ratio=0.3)
        r_l.next_to(r_ar.get_end(), RIGHT, buff=0.12).align_to(r_ar, UP)
        t_l.next_to(t_ar.get_end(), RIGHT, buff=0.12).align_to(t_ar, DOWN)
        ps = pts[star]
        lab_c.move_to(np.array([ps[0] - 1.25, rBand.c[1], 0.]))
        lab_h.move_to(np.array([ps[0] + 1.25, rBand.c[1], 0.]))
        lead = always_redraw(lambda: VGroup(
            Line(lab_c.get_bottom() + DOWN * 0.06, grain_geom(star)[0] + grain_geom(star)[1] * 0.3 * k,
                 stroke_width=2.5, color=POLE),
            Line(lab_h.get_bottom() + DOWN * 0.06, grain_geom(star)[0] - grain_geom(star)[2] * 0.46 * k,
                 stroke_width=2.5, color=WHITE)))

        def caption():
            ph = phi.get_value()
            m, op = (caps[0], np.clip((25 - ph) / 15, 0, 1)) if ph < 45 else \
                    (caps[1], np.clip((ph - 65) / 15, 0, 1))
            return place(m.copy(), rCap).set_opacity(op)
        cap = always_redraw(caption)

        # ---------------- zoom: the circled grain's hcp cell ----------------
        # a magnifier level with the grain, dashed lines from the small circle; inside,
        # the cell in parallel projection: c in the T-R plane along the grain's c axis,
        # tilted by TILT about the in-plane normal so that the basal plane shows
        Rz = min(cellZ.w / 2 - 0.05, 1.45)
        rs = 0.68 * k
        Cz = np.array([cellZ.c[0], ps[1], 0.])
        mag = Circle(radius=Rz, stroke_color=WHITE, stroke_width=2.5).move_to(Cz)
        links = VGroup(*[DashedLine(ps + sg * UP * rs, Cz + sg * UP * Rz, stroke_width=2, color=GREY_B,
                                    dash_length=0.1) for sg in (1, -1)])
        lab_z.next_to(mag, DOWN, buff=0.18)
        lab_cz = tex('c', L2M, POLE)
        A_ = 0.5 * Rz                                  # hexagon radius in the zoom
        H_ = COA * A_
        ARR = 0.36 * Rz

        def cell_zoom():
            _, c, d, w = grain_geom(star)
            ct, st = np.cos(TILT), np.sin(TILT)

            def P(sc, sd, sz):                         # (along c, along d, out of screen) -> screen
                return Cz + (sc * ct - sz * st) * c + sd * d

            def depth(sc, sz):
                return sc * st + sz * ct
            ang = [np.radians(15 + 60 * j) for j in range(6)]
            hexa = [(A_ * np.cos(a), A_ * np.sin(a)) for a in ang]          # (sd, sz)
            side_vis = [np.sin(np.radians(45 + 60 * j)) > 0 for j in range(6)]  # face j: between j, j+1
            back, front = VGroup(), VGroup()
            for j in range(6):
                j1 = (j + 1) % 6
                (d0, z0), (d1, z1) = hexa[j], hexa[j1]
                # top hexagon: always seen; bottom: seen with its side face
                front.add(Line(P(H_ / 2, d0, z0), P(H_ / 2, d1, z1), stroke_color=ZR, stroke_width=3.5))
                (front if side_vis[j] else back).add(
                    Line(P(-H_ / 2, d0, z0), P(-H_ / 2, d1, z1), stroke_color=ZR, stroke_width=3.5))
                vis = side_vis[j] or side_vis[j - 1]
                (front if vis else back).add(
                    Line(P(-H_ / 2, d0, z0), P(H_ / 2, d0, z0), stroke_color=ZR, stroke_width=3.5))
            back.set_stroke(opacity=0.35)
            atoms = VGroup()
            for sc in (-H_ / 2, H_ / 2):
                for (dd, zz) in hexa + [(0., 0.)]:
                    seen = sc > 0 or depth(sc, zz) > depth(-H_ / 2, 0.)
                    atoms.add(Dot(P(sc, dd, zz), radius=0.05, color=ZR).set_opacity(1 if seen else 0.4))
            # the hydride: a plate in the basal plane at mid-height, coloured as the grain's
            plate = Polygon(*[P(0., 0.92 * dd, 0.92 * zz) for dd, zz in hexa], stroke_width=0,
                            fill_color=interpolate_color(TAN, RAD, w), fill_opacity=0.92)
            tip = P(H_ / 2 + ARR + 0.06, 0., 0.)
            arrow = Arrow(P(H_ / 2 + 0.04, 0., 0.), tip, buff=0, color=POLE, stroke_width=6,
                          max_tip_length_to_length_ratio=0.35)
            # the c axis through the cell, thin
            axis = DashedLine(P(-H_ / 2 - 0.12, 0., 0.), P(H_ / 2, 0., 0.), stroke_width=2,
                              color=POLE, dash_length=0.06)
            lab = lab_cz.copy().move_to(tip + c * 0.24 + d * 0.0)
            return VGroup(back, atoms[:len(hexa) + 1], axis, plate, front, atoms[len(hexa) + 1:], arrow, lab)
        zcell = always_redraw(cell_zoom)

        # ---------------- pole figure ----------------
        side_w = pT.width + 0.15
        cb = Box(cellP.x0, cellP.x1 - side_w, cellP.y0 + pf_lab.height + 0.3,
                 cellP.y1 - pL.height - 0.15)
        Rp = min(cb.w, cb.h) / 2
        pc = cb.c
        pf = Circle(radius=Rp, stroke_color=GREY_B, stroke_width=3).move_to(pc)
        cross = VGroup(Line(pc + LEFT * Rp, pc + RIGHT * Rp, stroke_width=1.5, color=GREY_D),
                       Line(pc + DOWN * Rp, pc + UP * Rp, stroke_width=1.5, color=GREY_D))
        pT.next_to(pf, RIGHT, buff=0.15)
        pL.next_to(pf, UP, buff=0.15)
        pf_lab.next_to(pf, DOWN, buff=0.3)
        pR.move_to(pc + (DOWN + RIGHT) * 0.42 * Rp / 1.3)
        pR.set_stroke(BG, width=8, background=True)

        disk = Circle(radius=Rp).move_to(pc)

        def blobs():
            # a pole at the rim spreads past it: clip to the projection circle
            x = np.tan(np.radians(phi.get_value()) / 2) * Rp
            g = VGroup()
            for sgn in (-1, 1):
                for rr, op in ((0.42, 0.18), (0.28, 0.38), (0.15, 0.95)):
                    rr *= Rp / 1.3
                    e = Ellipse(width=2 * rr, height=2 * rr * 1.25).move_to(pc + RIGHT * sgn * x)
                    g.add(Intersection(e, disk, stroke_width=0, fill_color=POLE, fill_opacity=op))
            return g
        pole = always_redraw(blobs)

        report({'axes': VGroup(r_ar, t_ar, r_l, t_l), 'lab_c': lab_c, 'lab_h': lab_h,
                'wall': wall, 'caption': place(caps[0].copy(), rCap),
                'zoom': VGroup(mag, lab_z), 'pf': VGroup(pf, pT, pL, pf_lab)})

        self.add(wall, gr, links, lead, lab_c, lab_h, r_ar, r_l, t_ar, t_l, cap,
                 mag, zcell, lab_z, pf, cross, pole, pT, pL, pR, pf_lab)
        if os.environ.get('STILL'):
            phi.set_value(float(os.environ['STILL']))
            self.update_mobjects(0)
            self.wait(0.04)
            return
        self.wait(2.0)
        self.play(phi.animate.set_value(90.), run_time=5.0, rate_func=smooth)
        self.wait(2.0)
        self.play(phi.animate.set_value(0.), run_time=5.0, rate_func=smooth)
        self.wait(0.04)
