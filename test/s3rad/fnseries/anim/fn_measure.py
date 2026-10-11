"""Looping GIF: how F_n is measured on a micrograph.

The micrograph is Fig. 4 (110 MPa) of Plyasov et al., Phys. At. Nucl. 86
(2023) 2604 - unirradiated E635 tube after a reorientation test, R up, T to
the right.  The hydrides were traced with fnmeasure/trace.py (370 straight
pieces); a piece is radial if it lies within 45 deg of the local radial
direction, and F_n is the radial share of the total LENGTH, as in the paper
(F_n = 0.42 here; the paper gives 0.28-0.44 at 110 MPa).

Story: the rule (a needle turning from R to T changes colour at 45 deg), then
a scan colours every hydride and fills a length bar, then the formula.

Layout: CONTENT = [micrograph | column: rule, legend, bar, formula]; the
micrograph is pre-scaled with Lanczos to its exact on-screen pixel size.

    manim -r 1600,900 --fps 25 fn_measure.py FnMeasure
"""
import json
from PIL import Image
from style import *

HERE = '/home/user/fn-series/anim/fnmeasure/'
IMG = HERE + 'plyasov2023_fig4_110MPa.png'
CROP = (0, 12, 389, 282)                      # x0, y0, x1, y1 in panel pixels
PX = 1600 / config.frame_width                # output pixels per scene unit


class FnMeasure(Scene):
    def construct(self):
        self.camera.background_color = BG
        D = json.load(open(HERE + 'p110.json'))
        segs = D['segs']
        cw, ch = CROP[2] - CROP[0], CROP[3] - CROP[1]

        # ---------------- measure the right column ----------------
        leg = legend(L2)
        frm = tex(r'F_n = \frac{\sum L_{\text{рад}}}{\sum L}', L2M)
        val = tex(r'\approx 0{,}4', L1, ACC)
        col_w = max(leg.width, frm.width + 0.25 + val.width) + 0.05
        gap = 0.5
        img_w = CONTENT.w - gap - col_w - 0.03
        img_h = img_w * ch / cw
        if img_h > CONTENT.h:
            img_h = CONTENT.h
            img_w = img_h * cw / ch
        cellI = Box(CONTENT.x0 + 0.03, CONTENT.x0 + 0.03 + img_w, CONTENT.c[1] - img_h / 2, CONTENT.c[1] + img_h / 2)
        cellR = Box(CONTENT.x1 - col_w, CONTENT.x1, CONTENT.y0, CONTENT.y1)

        # ---------------- micrograph at its exact pixel size ----------------
        wpx, hpx = int(round(img_w * PX)), int(round(img_h * PX))
        pil = Image.open(IMG).convert('RGB').crop(CROP).resize((wpx, hpx), Image.LANCZOS)
        img = ImageMobject(np.asarray(pil))
        img.stretch_to_fit_width(img_w).stretch_to_fit_height(img_h).move_to(cellI.c)
        dim = Rectangle(width=img_w, height=img_h, stroke_width=0, fill_color=BLACK,
                        fill_opacity=1).move_to(cellI.c).set_opacity(0)
        frame = Rectangle(width=img_w, height=img_h, stroke_color=GREY_B, stroke_width=2,
                          fill_opacity=0).move_to(cellI.c)
        # black strips just outside the micrograph: the glow stops at its edge
        clip = VGroup(*[Rectangle(width=0.3, height=img_h + 0.4, stroke_width=0, fill_color=BG,
                                  fill_opacity=1).move_to(np.array([x, cellI.c[1], 0.]))
                        for x in (cellI.x0 - 0.15, cellI.x1 + 0.15)])

        def to_scene(x, y):
            return np.array([cellI.x0 + (x - CROP[0]) / cw * img_w,
                             cellI.y1 - (y - CROP[1]) / ch * img_h, 0.])
        lines, xs_mid, Ls, rad = VGroup(), [], [], []
        for s in segs:
            r = s['ang'] <= 45
            ln = neon(to_scene(s['x0'], s['y0']), to_scene(s['x1'], s['y1']), RAD if r else TAN, 3.5)
            lines.add(ln)
            xs_mid.append(ln.get_center()[0]); Ls.append(s['L']); rad.append(r)
        xs_mid, Ls, rad = map(np.array, (xs_mid, Ls, rad))
        Ltot = Ls.sum()
        print('LAYOUT image %.2f x %.2f (%d x %d px), column %.2f; F_n = %.3f'
              % (img_w, img_h, wpx, hpx, col_w, Ls[rad].sum() / Ltot))

        # ---------------- right column: rule, legend, bar, formula ----------------
        lab_R, lab_T = tex('R', L2M, GREY_A), tex('T', L2M, GREY_A)
        t45 = tex(r'45^\circ', L3M, GREY_A)
        bar_h = 0.34
        rg = 0.28
        key_h = cellR.h - leg.height - bar_h - frm.height - 3 * rg
        rK, rLg, rB, rF = cellR.stack(key_h, leg.height, bar_h, frm.height, gap=rg)
        rad_k = min((rK.h - lab_R.height - 0.12) / 2, (rK.w - 2 * (lab_T.width + 0.15)) / 2)
        kc = np.array([rK.c[0], rK.y0 + rad_k, 0.])
        sectors = VGroup()
        for a0, col in ((45, RAD), (225, RAD), (-45, TAN), (135, TAN)):
            sectors.add(Sector(radius=rad_k, start_angle=np.radians(a0), angle=PI / 2, arc_center=kc,
                               fill_color=col, fill_opacity=0.22, stroke_width=0))
        circ = Circle(radius=rad_k, stroke_color=GREY_B, stroke_width=2).move_to(kc)
        bounds = VGroup(*[DashedLine(kc - rad_k * np.array([np.cos(a), np.sin(a), 0.]),
                                     kc + rad_k * np.array([np.cos(a), np.sin(a), 0.]),
                                     stroke_width=2, color=GREY_B, dash_length=0.08)
                          for a in (np.radians(45), np.radians(135))])
        lab_R.next_to(circ, UP, buff=0.12)
        lab_T.next_to(circ, RIGHT, buff=0.15)
        ang45 = Arc(radius=0.55 * rad_k, start_angle=PI / 2, angle=-PI / 4, arc_center=kc,
                    stroke_color=GREY_A, stroke_width=3)
        u45 = np.array([np.cos(PI / 4), np.sin(PI / 4), 0.])
        t45.next_to(kc + rad_k * u45, UR, buff=0.04)
        print('LAYOUT key radius %.2f' % rad_k)
        phi = ValueTracker(0.)                     # needle angle from R, degrees

        def needle():
            a = np.radians(90 - phi.get_value())
            u = np.array([np.cos(a), np.sin(a), 0.])
            return neon(kc - u * rad_k * 0.88, kc + u * rad_k * 0.88,
                        RAD if phi.get_value() <= 45 else TAN, 6.)
        ndl = always_redraw(needle)

        place(leg, rLg, LEFT)
        scan = ValueTracker(cellI.x0)
        bar_op = ValueTracker(1.)
        bar_w = rB.w
        bar_frame = Rectangle(width=bar_w, height=bar_h, stroke_color=GREY_B, stroke_width=2).move_to(rB.c)

        def bar():
            sel = xs_mid <= scan.get_value()
            lr, lt = Ls[sel & rad].sum() / Ltot * bar_w, Ls[sel & ~rad].sum() / Ltot * bar_w
            g = VGroup()
            if lr > 1e-3:
                g.add(Rectangle(width=lr, height=bar_h, stroke_width=0, fill_color=RAD,
                                fill_opacity=1).move_to(np.array([rB.x0 + lr / 2, rB.c[1], 0.])))
            if lt > 1e-3:
                g.add(Rectangle(width=lt, height=bar_h, stroke_width=0, fill_color=TAN,
                                fill_opacity=1).move_to(np.array([rB.x0 + lr + lt / 2, rB.c[1], 0.])))
            return g.set_opacity(bar_op.get_value())
        bars = always_redraw(bar)
        frm.move_to(np.array([rF.x0, rF.c[1], 0.]), aligned_edge=LEFT)
        val.next_to(frm, RIGHT, buff=0.25)
        form = VGroup(frm, val)

        lop = ValueTracker(1.)                     # overall opacity of the traced sticks

        def show_lines(m):
            x, k = scan.get_value(), lop.get_value()
            for ln, xm in zip(m, xs_mid):
                set_neon(ln, k if xm <= x else 0.)
        lines.add_updater(show_lines)
        scan_line = always_redraw(lambda: Line(
            np.array([scan.get_value(), cellI.y0, 0.]), np.array([scan.get_value(), cellI.y1, 0.]),
            stroke_width=3, color=WHITE).set_opacity(
            0.8 if cellI.x0 + 0.01 < scan.get_value() < cellI.x1 - 0.01 else 0.))

        report({'image': frame, 'key': VGroup(circ, lab_T), 'R': lab_R, '45': t45, 'legend': leg,
                'bar': bar_frame, 'formula': form})

        # ---------------- timeline (first frame = last frame) ----------------
        self.add(img, dim, lines, clip, frame, sectors, bounds, circ, ang45, t45, lab_R, lab_T, ndl,
                 leg, bar_frame, bars, scan_line)
        self.wait(1.0)
        self.play(phi.animate.set_value(90.), run_time=1.8, rate_func=smooth)
        self.play(phi.animate.set_value(0.), run_time=1.8, rate_func=smooth)
        self.play(dim.animate.set_opacity(0.6), run_time=0.6)
        self.play(scan.animate.set_value(cellI.x1), run_time=6.0, rate_func=linear)
        self.play(Write(frm), FadeIn(val, shift=LEFT * 0.2), run_time=1.0)
        self.wait(2.6)
        self.play(FadeOut(form), lop.animate.set_value(0.), dim.animate.set_opacity(0),
                  bar_op.animate.set_value(0.), run_time=1.0)
        scan.set_value(cellI.x0)
        bar_op.set_value(1.)
        lop.set_value(1.)
        self.wait(0.04)
