"""Looping GIF, fracture mechanism 2: cleavage of the hydride (illustrative).

Tension along a tangential plate: the plate stretches together with the
metal.  The hydride follows the same damage law as the metal but with a
much smaller ductility, so its omega reaches 1 first: the plate cracks
across in a few places, the pieces move apart with the metal, and the metal
starts to glow (omega) at the crack mouths.  Right: omega of the hydride and
of the metal side by side.

    manim -r 1600,900 --fps 25 mech_cleave.py MechCleave
"""
import os
from mech_common import *

T_LOOP = 12.0
L = 2.4                               # a smaller block: the plate is the subject
PL, PT = 3.4, 0.36                    # tangential plate: length (along the load) x thickness
T_LOAD = (0.6, 6.8)                   # the stretch keeps growing
T_H1 = 3.6                            # omega of the hydride reaches 1
XC = [0.08, -0.92, 0.95]              # cleavage sites, in the order they crack
DT_C = 0.3                            # time between the cracks
GAP = 0.18                            # final opening of each crack
T_DIP = (9.8, 10.5, 10.8, 11.5)
EPS = 0.06                            # stretch of the picture at full load (exaggerated)


def omega_h(t):
    return ramp(t, T_LOAD[0], T_H1)


def omega_m(t):
    """The metal at the crack mouths: slow, then faster once the plate is cracked."""
    return 0.45 * ramp(t, T_LOAD[0], T_H1) + 0.45 * ramp(t, T_H1, 8.0)


class MechCleave(Scene):
    def construct(self):
        self.camera.background_color = BG

        side = CONTENT.h
        cellL = Box(CONTENT.x0, CONTENT.x0 + side, CONTENT.y0, CONTENT.y1)
        cellR = Box(cellL.x1 + 0.7, CONTENT.x1, CONTENT.y0, CONTENT.y1)

        blk = block(L)
        arrows = load_arrows(L)
        whole = VGroup(blk, arrows)
        w0 = whole.width
        fit(whole, cellL)
        sc = whole.width / w0
        sh = blk.get_center()

        def P(x, y):
            return sh + np.array([x, y, 0.]) * sc
        t_c = [T_H1 + i * DT_C for i in range(len(XC))]
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
            e = EPS * ramp(te, *T_LOAD)
            out = VGroup()
            cracked = [x for x, tc in zip(XC, t_c) if te >= tc]
            opens = {x: GAP * ramp(te, tc, tc + 2.5) for x, tc in zip(XC, t_c) if te >= tc}

            def X(x0):                                   # material x -> current x
                return x0 * (1 + e) + sum(o / 2 * np.sign(x0 - xc) for xc, o in opens.items())
            # the metal glows at the mouths of the cracks
            for xc in cracked:
                for sg in (1, -1):
                    out.add(glow(P(X(xc), sg * (PT / 2 + 0.25)), omega_m(te),
                                 0.55 * sc, 0.7 * sc, n=8, clip=blk))
            # plate pieces between the cracks
            edges = [-PL / 2] + sorted(cracked) + [PL / 2]
            col = interpolate_color(TAN, WHITE, 0.35 * omega_h(te) * (len(cracked) == 0))
            for a, b in zip(edges[:-1], edges[1:]):
                xa, xb = X(a) + (0.0 if a == -PL / 2 else 0.01), X(b) - (0.0 if b == PL / 2 else 0.01)
                if a in opens:
                    xa = X(a + 1e-6)
                if b in opens:
                    xb = X(b - 1e-6)
                r = PT / 2 if (a == -PL / 2 and b == PL / 2) else 0.05
                out.add(RoundedRectangle(width=(xb - xa) * sc, height=PT * sc, corner_radius=r * sc,
                                         stroke_width=0, fill_color=col, fill_opacity=1)
                        .move_to(P((xa + xb) / 2, 0)))
            # a flash where a crack has just appeared
            for xc, tc in zip(XC, t_c):
                f = 1 - ramp(te, tc, tc + 0.45) if te >= tc else 0.
                if f > 0.01:
                    out.add(Line(P(X(xc), -PT * 1.4), P(X(xc), PT * 1.4), stroke_color=WHITE,
                                 stroke_width=6).set_stroke(opacity=f))
            return fade_all(out, g) if g < 1 else out
        pic = always_redraw(picture)
        blk_d = always_redraw(lambda: fade_all(blk.copy(), dip(clk.get_value())))

        # ---------------- right column: omega of the hydride and of the metal ----------------
        head = VGroup(txt('исчерпание пластичности', L3, GREY_A), tex(r'\omega', L1, ACC)).arrange(RIGHT, buff=0.25)
        lab_h, lab_m = txt('гидрид', L2, TAN), txt('металл', L2, GREY_A)
        lw = max(lab_h.width, lab_m.width) + 0.3
        bw = min(cellR.w - lw - 0.1, 4.0)
        bh = 0.36
        fr_h = Rectangle(width=bw, height=bh, stroke_color=GREY_B, stroke_width=2)
        fr_m = fr_h.copy()
        cap = VGroup(txt('пластичность гидрида', L2), txt('ниже, чем у металла —', L2),
                     txt('он раскалывается первым', L2)).arrange(DOWN, aligned_edge=LEFT, buff=0.16)
        rows = VGroup(VGroup(lab_h, fr_h), VGroup(lab_m, fr_m))
        for lab, fr in rows:
            fr.next_to(lab, RIGHT, buff=0.3)
        rows.arrange(DOWN, aligned_edge=LEFT, buff=0.45)
        fr_m.align_to(fr_h, LEFT)
        col = VGroup(head, rows, cap).arrange(DOWN, aligned_edge=LEFT, buff=0.6)
        place(col, cellR, LEFT)
        fr_m.align_to(fr_h, LEFT)

        def bars():
            te = te_of(clk.get_value())
            g = VGroup()
            for w, fr in ((omega_h(te), fr_h), (omega_m(te), fr_m)):
                if w > 1e-3:
                    g.add(Rectangle(width=bw * w, height=bh, stroke_width=0, fill_color=omega_color(w),
                                    fill_opacity=1).move_to(fr.get_left() + RIGHT * bw * w / 2))
            return g
        fills = always_redraw(bars)
        one_h = tex('1', L2M, GREY_A).next_to(fr_h, RIGHT, buff=0.12)
        report({'picture': whole, 'column': col, 'one': one_h})

        self.add(blk_d, pic, arrows, col, fills, one_h)
        if os.environ.get('STILL'):
            clk.set_value(float(os.environ['STILL']))
            self.update_mobjects(0)
            self.wait(0.04)
            return
        self.play(clk.animate.set_value(T_LOOP), run_time=T_LOOP, rate_func=linear)
