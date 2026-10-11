"""Looping GIF, fracture mechanism 3: debonding of the hydride / metal
interface (illustrative).

Tension across a radial plate: the bonds on its faces (cohesive "springs")
stretch; the bond force T follows the bilinear law - it rises to the
strength T0, then falls to zero; the area under it is the debonding energy
G_c.  The faces open from the middle out, the springs snap, and the plate
works as a crack.  Right: the T - delta law with the state of the middle
spring running along it.

    manim -r 1600,900 --fps 25 mech_debond.py MechDebond
"""
import os
from mech_common import *

T_LOOP = 12.0
L = 2.2                               # a smaller block: the interface is the subject
PW, PH = 0.42, 3.0                    # a bigger plate: the faces are the subject here
A = PH / 2                            # half-length of the bonded faces
D0, DC = 0.05, 0.32                   # opening at the strength peak, at full separation
T_OPEN = (0.6, 7.0)
T_DIP = (9.8, 10.5, 10.8, 11.5)
NS = 7                                # springs per face


def opening(t):
    """Opening of the middle of the faces."""
    return 0.85 * ramp(t, *T_OPEN) ** 1.6


def law(d):
    """Bilinear traction-separation, T / T0."""
    if d <= D0:
        return d / D0
    if d < DC:
        return (DC - d) / (DC - D0)
    return 0.


class MechDebond(Scene):
    def construct(self):
        self.camera.background_color = BG

        side = CONTENT.h
        cellL = Box(CONTENT.x0, CONTENT.x0 + side, CONTENT.y0, CONTENT.y1)
        cellR = Box(cellL.x1 + 0.7, CONTENT.x1, CONTENT.y0, CONTENT.y1)

        blk = block(L)
        plate = RoundedRectangle(width=PW, height=PH, corner_radius=PW / 2, stroke_width=0,
                                 fill_color=RAD, fill_opacity=1)
        arrows = load_arrows(L)
        whole = VGroup(blk, plate, arrows)
        w0 = whole.width
        fit(whole, cellL)
        sc = whole.width / w0
        sh = blk.get_center()

        def P(x, y):
            return sh + np.array([x, y, 0.]) * sc
        ys = np.linspace(-0.8 * A, 0.8 * A, NS)
        clk = ValueTracker(0.)

        def dip(t):
            if t < T_DIP[0]:
                return 1.
            return 1 - ramp(t, T_DIP[0], T_DIP[1]) + ramp(t, T_DIP[2], T_DIP[3])

        def te_of(t):
            return t if t < T_DIP[1] + 0.15 else 0.

        def picture():
            t = clk.get_value()
            g, te = dip(t), te_of(t)
            dm = opening(te)
            out = VGroup()
            yy = np.linspace(-A, A, 60)
            dd = dm * np.sqrt(np.clip(1 - (yy / A) ** 2, 0, 1))
            for sg in (-1, 1):                          # the gaps on both faces
                if dm > 1e-3:
                    pts = [P(sg * PW / 2, y) for y in yy] + \
                          [P(sg * (PW / 2 + d), y) for y, d in zip(yy[::-1], dd[::-1])]
                    out.add(Polygon(*pts, stroke_width=0, fill_color=BLACK, fill_opacity=1))
                for y in ys:                            # the springs
                    d = dm * np.sqrt(max(0., 1 - (y / A) ** 2))
                    if d >= DC:
                        continue
                    dmg = 0. if d <= D0 else (d - D0) / (DC - D0)
                    col = interpolate_color(IFC, GREY_D, dmg)
                    x0, x1 = sg * PW / 2, sg * (PW / 2 + max(d, 0.1))
                    zz = [P(x0 + (x1 - x0) * q, y + (0.07 if i % 2 else -0.07) * (0 < i < 6))
                          for i, q in enumerate(np.linspace(0, 1, 7))]
                    out.add(VMobject(stroke_color=col, stroke_width=5 - 2.5 * dmg)
                            .set_points_as_corners(zz))
            out.add(plate.copy())
            return fade_all(out, g) if g < 1 else out
        pic = always_redraw(picture)
        blk_d = always_redraw(lambda: fade_all(blk.copy(), dip(clk.get_value())))

        # ---------------- right column: the T - delta law ----------------
        cap = VGroup(txt('граница рвётся —', L2), txt('пластина работает', L2),
                     txt('как трещина', L2)).arrange(DOWN, aligned_edge=LEFT, buff=0.16)
        lab_T, lab_d = tex('T', L1), tex(r'\delta', L1)
        gw = cellR.w - lab_T.width - lab_d.width - 0.6
        gh = cellR.h - cap.height - 0.9 - 0.45
        ax = Axes(x_range=[0, 1.25 * DC, DC], y_range=[0, 1.5, 1], x_length=gw, y_length=gh,
                  tips=True, axis_config=dict(color=GREY_B, stroke_width=3, include_ticks=False,
                                               tip_width=0.2, tip_height=0.2))
        ax.move_to(np.array([cellR.x0 + lab_T.width + 0.3 + gw / 2, cellR.y1 - 0.15 - gh / 2, 0.]))
        lab_T.next_to(ax.y_axis, LEFT, buff=0.25).align_to(ax.y_axis, UP)
        lab_d.next_to(ax.x_axis, RIGHT, buff=0.15)
        tri = VMobject(stroke_color=IFC, stroke_width=5).set_points_as_corners(
            [ax.c2p(0, 0), ax.c2p(D0, 1), ax.c2p(DC, 0)])
        t0 = tex('T^0', L2M, GREY_A).next_to(ax.c2p(0, 1), LEFT, buff=0.2)
        t0_line = DashedLine(ax.c2p(0, 1), ax.c2p(D0, 1), stroke_width=2, color=GREY_B)
        gc = tex('G_c', L1, IFC).move_to(ax.c2p(0.42 * DC, 0.32))
        place(cap, Box(cellR.x0, cellR.x1, cellR.y0, cellR.y0 + cap.height), LEFT)

        def area():
            d = min(opening(te_of(clk.get_value())), DC)
            g = VGroup()
            if d > 1e-4:
                xs = np.linspace(0, d, 40)
                pts = [ax.c2p(x, law(x)) for x in xs] + [ax.c2p(d, 0), ax.c2p(0, 0)]
                g.add(Polygon(*pts, stroke_width=0, fill_color=IFC, fill_opacity=0.35))
                g.add(Dot(ax.c2p(d, law(d)), radius=0.12, color=WHITE))
            return g
        fill = always_redraw(area)
        gc.add_updater(lambda m: m.set_opacity(ramp(opening(te_of(clk.get_value())), 0.7 * DC, DC)))
        report({'picture': VGroup(blk, arrows), 'graph': VGroup(ax, lab_T, lab_d, t0), 'caption': cap})

        self.add(blk_d, pic, arrows, ax, lab_T, lab_d, fill, tri, t0_line, t0, gc, cap)
        if os.environ.get('STILL'):
            clk.set_value(float(os.environ['STILL']))
            self.update_mobjects(0)
            self.wait(0.04)
            return
        self.play(clk.animate.set_value(T_LOOP), run_time=T_LOOP, rate_func=linear)
