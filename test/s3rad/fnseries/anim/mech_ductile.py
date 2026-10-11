"""Looping GIF, fracture mechanism 1: ductile fracture of the metal
(illustrative, no FE mesh).

As in the solver (calcdamage.f, DE1): the initiation variable
omega = int d(eps_p) / eps_f(eta) accumulates; at omega = 1 damage starts and
D grows from 0 with the further plastic strain (sigma = (1 - D) sigma~); at
D = 1 the element is deleted.  Illustration: tension across a radial plate,
the metal glows (omega) at the tips; at a marked point ahead of the tip the
two bars on the right fill in turn: omega -> 1 (a void appears), D -> 1 (the
void links with the crack, which runs on to the free edge).

Layout: CONTENT = [square cell: block, plate, load arrows | column: formula,
omega bar, D bar, their captions].  One clock t drives everything; the loop
dips the picture to black and returns to the start.

    manim -r 1600,900 --fps 25 mech_ductile.py MechDuctile
"""
import os
from mech_common import *

T_LOOP = 12.0
L = 3.0                              # block half-size
PW, PH = 0.30, 2.2                   # radial plate: width (along the load) x height
T_LOAD = (0.6, 3.0)                  # the glow builds up at the tips
T_CRACK = (3.0, 7.6)                 # voids link up, the crack runs to the edges
T_DIP = (9.8, 10.5, 10.8, 11.5)      # fade out, black, fade in
MARK = 2                             # the tracked material point: void no. 2 above the top tip
DCOL = ManimColor('#C9CED6')         # D bar


def void_sites(sign, rng):
    """Void centres ahead of one tip, from the tip to the edge."""
    ys = PH / 2 + 0.28 + np.arange(6) * 0.33
    xs = np.cumsum(rng.normal(0, 0.07, len(ys)))
    return [np.array([x, sign * y, 0.]) for x, y in zip(xs, ys) if y < L - 0.12]


class MechDuctile(Scene):
    def construct(self):
        self.camera.background_color = BG
        rng = np.random.default_rng(4)

        # ---------------- layout ----------------
        side = CONTENT.h
        cellL = Box(CONTENT.x0, CONTENT.x0 + side, CONTENT.y0, CONTENT.y1)
        cellR = Box(cellL.x1 + 0.7, CONTENT.x1, CONTENT.y0, CONTENT.y1)

        # ---------------- the picture: domain units, then one scale + shift ----------------
        blk = block(L)
        plate = RoundedRectangle(width=PW, height=PH, corner_radius=PW / 2, stroke_width=0,
                                 fill_color=RAD, fill_opacity=1)
        arrows = load_arrows(L)
        whole = VGroup(blk, plate, arrows)
        w0 = whole.width
        fit(whole, cellL)
        sc = whole.width / w0
        sh = blk.get_center()

        def P(p):                                       # domain -> scene
            return sh + np.asarray(p) * sc
        tips = [np.array([0., PH / 2, 0.]), np.array([0., -PH / 2, 0.])]
        sites = [void_sites(1, rng), void_sites(-1, rng)]
        n_s = len(sites[0])
        t_link = np.linspace(T_CRACK[0] + 0.5, T_CRACK[1], n_s)   # D = 1: linked into the crack
        t_nuc = t_link - 1.6                                       # omega = 1: a void appears
        clk = ValueTracker(0.)

        def dip(t):
            if t < T_DIP[0]:
                return 1.
            return 1 - ramp(t, T_DIP[0], T_DIP[1]) + ramp(t, T_DIP[2], T_DIP[3])

        def te_of(t):
            return t if t < T_DIP[1] + 0.15 else 0.

        def omega_D(te):
            """omega and D at the tracked point."""
            w = ramp(te, T_LOAD[0], t_nuc[MARK])
            D = ramp(te, t_nuc[MARK], t_link[MARK])
            return w, D

        def picture():
            t = clk.get_value()
            g, te = dip(t), te_of(t)
            out = VGroup()
            I_tip = ramp(te, *T_LOAD)
            for side_i, (tip, ss) in enumerate(zip(tips, sites)):
                n_link = int(np.sum(te >= t_link))
                front = tip if n_link == 0 else ss[n_link - 1]
                sgn = 1 if side_i == 0 else -1
                if n_link < n_s:                        # plastic zone ahead of the crack front
                    out.add(glow(P(front + np.array([0, sgn * 0.3, 0])), I_tip, 0.62 * sc, 0.9 * sc, clip=blk))
                pts = [tip] + ss[:n_link]
                if n_link == n_s:
                    pts.append(np.array([ss[-1][0], sgn * L, 0.]))
                if len(pts) > 1:
                    out.add(VMobject(stroke_color=BLACK, stroke_width=9 + 5 * n_link / n_s)
                            .set_points_as_corners([P(p) for p in pts]))
                for j, c in enumerate(ss):              # voids: omega = 1 .. D = 1
                    a = ramp(te, t_nuc[j], t_link[j])
                    if a <= 0 or te >= t_link[j] + 0.3:
                        continue
                    out.add(Ellipse(width=(0.05 + 0.2 * a) * sc, height=(0.04 + 0.11 * a) * sc,
                                    stroke_width=0, fill_color=BLACK, fill_opacity=1).move_to(P(c)))
            out.add(plate.copy())
            w, D = omega_D(te)                          # the tracked point
            ring_col = ACC if D <= 0 else interpolate_color(ACC, WHITE, D)
            out.add(Circle(radius=0.22 * sc, stroke_color=ring_col, stroke_width=3).move_to(
                P(sites[0][MARK])))
            return fade_all(out, g) if g < 1 else out
        pic = always_redraw(picture)
        blk_d = always_redraw(lambda: fade_all(blk.copy(), dip(clk.get_value())))

        # ---------------- right column: formula, omega bar, D bar ----------------
        frm = tex(r'\omega = \int \frac{d\bar{\varepsilon}^{\,p}}{\varepsilon_f(\eta)}', L1, ACC)
        bw = min(cellR.w - 1.2, 4.2)
        lab_w, lab_D = tex(r'\omega', L1, ACC), tex('D', L1, DCOL)
        cap_w = txt('ω = 1: начало повреждения', L2, GREY_B)
        cap_D = txt('D = 1: элемент удаляется', L2, GREY_B)
        bh = 0.36
        frame_w = Rectangle(width=bw, height=bh, stroke_color=GREY_B, stroke_width=2)
        frame_D = frame_w.copy()
        row_w = VGroup(lab_w, frame_w).arrange(RIGHT, buff=0.3)
        row_D = VGroup(lab_D, frame_D).arrange(RIGHT, buff=0.3)
        frame_D.align_to(frame_w, LEFT)
        lab_D.next_to(frame_D, LEFT, buff=0.3)
        grp_w = VGroup(row_w, cap_w).arrange(DOWN, aligned_edge=LEFT, buff=0.22)
        grp_D = VGroup(row_D, cap_D).arrange(DOWN, aligned_edge=LEFT, buff=0.22)
        col = VGroup(frm, grp_w, grp_D).arrange(DOWN, aligned_edge=LEFT, buff=0.6)
        place(col, cellR, LEFT)
        frame_D.align_to(frame_w, LEFT)
        lab_D.next_to(frame_D, LEFT, buff=0.3)
        cap_D.align_to(cap_w, LEFT)

        def bars():
            te = te_of(clk.get_value())
            w, D = omega_D(te)
            g = VGroup()
            if w > 1e-3:
                g.add(Rectangle(width=bw * w, height=bh, stroke_width=0, fill_color=omega_color(w),
                                fill_opacity=1).move_to(frame_w.get_left() + RIGHT * bw * w / 2))
            if D > 1e-3:
                g.add(Rectangle(width=bw * D, height=bh, stroke_width=0, fill_color=DCOL,
                                fill_opacity=1).move_to(frame_D.get_left() + RIGHT * bw * D / 2))
            return g
        fills = always_redraw(bars)

        def light(m, which):
            te = te_of(clk.get_value())
            w, D = omega_D(te)
            on = (w >= 0.999) if which == 'w' else (D >= 0.999)
            m.set_color(WHITE if on else GREY_B)
        cap_w.add_updater(lambda m: light(m, 'w'))
        cap_D.add_updater(lambda m: light(m, 'D'))
        report({'picture': VGroup(blk, arrows), 'formula': frm, 'omega': grp_w, 'D': grp_D})

        self.add(blk_d, pic, arrows, col, fills)
        if os.environ.get('STILL'):                  # layout check at a given time
            clk.set_value(float(os.environ['STILL']))
            self.update_mobjects(0)                  # -s skips the updaters otherwise
            self.wait(0.04)
            return
        self.play(clk.animate.set_value(T_LOOP), run_time=T_LOOP, rate_func=linear)
