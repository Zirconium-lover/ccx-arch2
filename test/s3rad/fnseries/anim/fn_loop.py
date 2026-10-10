"""Looping GIF: the 12 plates turn radial one by one (Fn 0 -> 100 %), hold,
and turn back in reverse order (100 -> 0 %).  First and last frames are the
same, so the GIF loops without a seam for as long as the slide is shown.

Plate centres and the order in which they turn are those of the series
decks (mkseeddeck.py --seed 5).

    manim -r 1280,720 --fps 25 fn_loop.py FnLoop
"""
import sys
import numpy as np
from manim import *

sys.path.insert(0, '/home/user/ccx-fn/test/s3rad')
import mkseeddeck as gen  # noqa: E402

config.tex_template = TexTemplate(preamble=r'''\usepackage[T2A]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage[english,russian]{babel}
\usepackage{amsmath}
\usepackage{amssymb}''')

RAD = ManimColor('#E8604C')
TAN = ManimColor('#2CB1D9')
BG = ManimColor('#000000')
LX, LY = 4.0, 3.6
S = 1.45                 # scene units per mm
C = gen.place_any(5)
ORDER = list(np.random.default_rng(10000 + 5).permutation(12))
STEP, HOLD0, HOLD1 = 0.42, 1.2, 1.6


class FnLoop(Scene):
    def construct(self):
        self.camera.background_color = BG
        centre = np.array([-1.6, 0.25, 0.])

        def at(x, y):
            return centre + np.array([(x - LX / 2) * S, (y - LY / 2) * S, 0.])

        box = Rectangle(width=LX * S, height=LY * S, stroke_color=GREY_B,
                        stroke_width=2.5).move_to(centre)
        arrows = VGroup()
        for dy in (-1.65, 0., 1.65):
            arrows.add(Arrow(box.get_left() + UP * dy, box.get_left() + UP * dy + LEFT * 0.85,
                             buff=0, color=WHITE, stroke_width=6),
                       Arrow(box.get_right() + UP * dy, box.get_right() + UP * dy + RIGHT * 0.85,
                             buff=0, color=WHITE, stroke_width=6))
        sig = MathTex(r'\sigma_{\theta}', font_size=44).next_to(arrows[5], UP, buff=0.12)
        r_ax = VGroup(Arrow(ORIGIN, UP * 0.85, buff=0, color=GREY_A, stroke_width=5),
                      MathTex('r', font_size=36, color=GREY_A))
        r_ax[1].next_to(r_ax[0], RIGHT, buff=0.08)
        r_ax.next_to(box.get_corner(UL), DR, buff=0.16)

        plates = VGroup(*[
            RoundedRectangle(width=0.6 * S, height=0.09 * S, corner_radius=0.045 * S,
                             stroke_width=0, fill_color=TAN, fill_opacity=1).move_to(at(c[0], c[1]))
            for c in C])

        k = ValueTracker(0)
        right = np.array([3.0, 0., 0.]) + RIGHT * 1.15
        lab = MathTex(r'F_n', '=', font_size=60).move_to(right + UP * 0.9 + LEFT * 0.9)
        frac = always_redraw(lambda: MathTex(r'{%d \over 12}' % round(k.get_value()),
                                             font_size=60).next_to(lab, RIGHT, buff=0.2))
        pct = always_redraw(lambda: DecimalNumber(
            100 * round(k.get_value()) / 12, num_decimal_places=0, unit=r'\,\%',
            font_size=64, color=YELLOW).next_to(lab, DOWN, buff=0.75).align_to(lab, LEFT))
        bar_bg = Rectangle(width=3.0, height=0.14, stroke_color=GREY_B, stroke_width=1.5)
        bar_bg.next_to(pct, DOWN, buff=0.45).align_to(lab, LEFT)
        bar = always_redraw(lambda: Rectangle(
            width=max(1e-3, 3.0 * k.get_value() / 12), height=0.14, stroke_width=0,
            fill_color=YELLOW, fill_opacity=1).align_to(bar_bg, LEFT).align_to(bar_bg, DOWN))

        self.add(box, arrows, sig, r_ax, plates, lab, frac, pct, bar_bg, bar)
        self.wait(HOLD0)

        def turn(i, j, to_radial):
            p = plates[i]
            start = p.copy()
            c0, c1 = (TAN, RAD) if to_radial else (RAD, TAN)
            sgn = 1 if to_radial else -1

            def upd(m, a, start=start):
                m.become(start.copy().rotate(sgn * a * PI / 2)
                         .set_fill(interpolate_color(c0, c1, a)))
            return AnimationGroup(UpdateFromAlphaFunc(p, upd), k.animate.set_value(j),
                                  run_time=STEP, rate_func=smooth)

        for j, i in enumerate(ORDER, start=1):
            self.play(turn(i, j, True))
        self.wait(HOLD1)
        for j, i in zip(range(11, -1, -1), reversed(ORDER)):
            self.play(turn(i, j, False))
        # the hold at the start of the next loop completes the pause
        self.wait(0.04)
