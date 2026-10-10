"""Looping GIF: the 12 plates of the series turn radial one by one
(Fn 0 -> 1), hold, and turn back.  First frame = last frame.

Layout from measurement: the right column is exactly as wide as the widest
state of the formula (L1) or the legend (L2); the specimen gets the rest and
is fitted into its cell.

    manim -r 1600,900 --fps 25 fn_loop.py FnLoop
"""
import sys
from style import *

sys.path.insert(0, '/home/user/ccx-fn/test/s3rad')
import mkseeddeck as gen  # noqa: E402

LX, LY = 4.0, 3.6
C = gen.place_any(5)
ORDER = list(np.random.default_rng(10000 + 5).permutation(12))
STEP, HOLD0, HOLD1 = 0.42, 1.2, 1.6


def fn_tex(k):
    v = '0' if k == 0 else '1' if k == 12 else ('%.2f' % (k / 12)).replace('.', '{,}')
    return r'F_n = {%d \over 12} = %s' % (k, v)


class FnLoop(Scene):
    def construct(self):
        self.camera.background_color = BG

        # ---- measure the right column ----
        widest = max((tex(fn_tex(k), L1) for k in range(13)), key=lambda m: m.width)
        leg = legend(L2)
        col_w = max(widest.width, leg.width) + 0.1
        gap = 0.6
        left = Box(CONTENT.x0, CONTENT.x1 - col_w - gap, CONTENT.y0, CONTENT.y1)
        right = Box(CONTENT.x1 - col_w, CONTENT.x1, CONTENT.y0, CONTENT.y1)

        # ---- specimen, built in mm and fitted into the left cell ----
        box = Rectangle(width=LX, height=LY, stroke_color=GREY_B, stroke_width=3)
        arrows = VGroup()
        for dy in (-1.15, 0., 1.15):
            for side, d in ((box.get_left(), LEFT), (box.get_right(), RIGHT)):
                arrows.add(Arrow(side + UP * dy, side + UP * dy + d * 0.5, buff=0, color=WHITE,
                                 stroke_width=7, max_tip_length_to_length_ratio=0.4))
        plates = VGroup(*[
            RoundedRectangle(width=0.6, height=0.12, corner_radius=0.06, stroke_width=0,
                             fill_color=TAN, fill_opacity=1).move_to(
                np.array([c[0] - LX / 2, c[1] - LY / 2, 0.])) for c in C])
        r_ax = VGroup(Arrow(ORIGIN, UP * 0.7, buff=0, color=GREY_A, stroke_width=6,
                            max_tip_length_to_length_ratio=0.3))
        r_ax.next_to(box.get_corner(UL), DR, buff=0.12)
        spec = VGroup(box, arrows, plates, r_ax)
        fit(spec, left)
        # labels at their final size (not scaled with the specimen)
        r_l = tex('r', L2M, GREY_A).next_to(r_ax, RIGHT, buff=0.12).align_to(r_ax, UP)
        sig = tex(r"\sigma_\theta", L2M).next_to(arrows[5], UP, buff=0.12)

        # ---- right column: formula over legend, centred as a block ----
        k = ValueTracker(0)
        block = VGroup(widest.copy(), leg).arrange(DOWN, buff=1.1, aligned_edge=LEFT)
        place(block, right, ORIGIN)
        anchor = block[0].get_left()
        formula = always_redraw(lambda: tex(fn_tex(round(k.get_value())), L1)
                                .move_to(anchor, aligned_edge=LEFT))
        report({'spec': spec, 'formula': block[0], 'legend': leg, 'sig': sig, 'r': r_l})

        self.add(spec, r_l, sig, formula, leg)
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
