"""Looping GIF: the 12 plates of the series turn radial one by one
(Fn 0 -> 100 %), hold, and turn back (100 -> 0 %).  First frame = last frame.

    manim -r 1600,900 --fps 25 fn_loop.py FnLoop
"""
import sys
import numpy as np
from style import *

sys.path.insert(0, '/home/user/ccx-fn/test/s3rad')
import mkseeddeck as gen  # noqa: E402

LX, LY = 4.0, 3.6
S = 1.55                       # scene units per mm
C = gen.place_any(5)
ORDER = list(np.random.default_rng(10000 + 5).permutation(12))
STEP, HOLD0, HOLD1 = 0.42, 1.2, 1.6


class FnLoop(Scene):
    def construct(self):
        self.camera.background_color = BG
        centre = np.array([-3.0, -0.45, 0.])

        def at(x, y):
            return centre + np.array([(x - LX / 2) * S, (y - LY / 2) * S, 0.])

        box = Rectangle(width=LX * S, height=LY * S, stroke_color=GREY_B,
                        stroke_width=3).move_to(centre)
        arrows = VGroup()
        for dy in (-1.9, 0., 1.9):
            arrows.add(Arrow(box.get_left() + UP * dy, box.get_left() + UP * dy + LEFT * 0.75,
                             buff=0, color=WHITE, stroke_width=8),
                       Arrow(box.get_right() + UP * dy, box.get_right() + UP * dy + RIGHT * 0.75,
                             buff=0, color=WHITE, stroke_width=8))
        sig = tex(r'\sigma_{\theta}', FS).next_to(arrows[5], UP, buff=0.15)
        r_ax = VGroup(Arrow(ORIGIN, UP * 1.1, buff=0, color=GREY_A, stroke_width=7),
                      tex('r', FS_SM, GREY_A))
        r_ax[1].next_to(r_ax[0], RIGHT, buff=0.1)
        r_ax.next_to(box.get_corner(UL), DR, buff=0.2)

        plates = VGroup(*[
            RoundedRectangle(width=0.6 * S, height=0.12 * S, corner_radius=0.06 * S,
                             stroke_width=0, fill_color=TAN, fill_opacity=1).move_to(at(c[0], c[1]))
            for c in C])

        k = ValueTracker(0)
        x0 = 2.0
        lab = tex(r'F_n =', FS_BIG).move_to(np.array([x0 + 0.9, 1.75, 0.]))
        frac = always_redraw(lambda: tex(r'{%d \over 12}' % round(k.get_value()), FS_BIG)
                             .next_to(lab, RIGHT, buff=0.25))
        pct = always_redraw(lambda: Integer(
            round(100 * round(k.get_value()) / 12), unit=r'\,\%', font_size=FS_HUGE + 8,
            color=ACC).move_to(np.array([x0 + 2.45, -0.2, 0.])))
        bar_bg = Rectangle(width=4.6, height=0.28, stroke_color=GREY_B, stroke_width=2)
        bar_bg.move_to(np.array([x0 + 2.45, -1.55, 0.]))
        bar = always_redraw(lambda: Rectangle(
            width=max(1e-3, 4.6 * k.get_value() / 12), height=0.28, stroke_width=0,
            fill_color=ACC, fill_opacity=1).align_to(bar_bg, LEFT).align_to(bar_bg, DOWN))
        legend = VGroup(
            VGroup(RoundedRectangle(width=0.16, height=0.75, corner_radius=0.06, stroke_width=0,
                                    fill_color=RAD, fill_opacity=1),
                   txt('радиальные', TX, RAD)).arrange(RIGHT, buff=0.3),
            VGroup(RoundedRectangle(width=0.75, height=0.16, corner_radius=0.06, stroke_width=0,
                                    fill_color=TAN, fill_opacity=1),
                   txt('тангенциальные', TX, TAN)).arrange(RIGHT, buff=0.3),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.35).move_to(np.array([x0 + 2.45, -2.9, 0.]))

        self.add(box, arrows, sig, r_ax, plates, lab, frac, pct, bar_bg, bar, legend)
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
        self.wait(0.04)
