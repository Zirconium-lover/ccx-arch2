"""Demo: the 12-plate model, Fn swept from 0 to 100 % (top view, x-y plane).

Plate centres and the nested order in which plates turn radial are the
exact ones of the series decks (mkseeddeck.py --seed 5).

    manim-slides render -qm demo_fn.py FnSweep
    manim-slides convert FnSweep fn_sweep.pptx --to=pptx
"""
import sys
import numpy as np
from manim import *
from manim_slides import Slide

sys.path.insert(0, '/home/user/ccx-fn/test/s3rad')
import mkseeddeck as gen  # noqa: E402

RU = TexTemplate(preamble=r'''\usepackage[T2A]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage[english,russian]{babel}
\usepackage{amsmath}
\usepackage{amssymb}''')
config.tex_template = RU

RAD = '#E8604C'      # radial plates (as in the presentation)
TAN = '#2CB1D9'      # tangential plates
BG = '#0E0E10'
LX, LY = 4.0, 3.6
S = 1.55             # scene units per mm
FONT = 'Noto Sans'

C = gen.place_any(5)
ORDER = list(np.random.default_rng(10000 + 5).permutation(12))  # = axes_for()


def to_scene(x, y):
    return np.array([(x - LX / 2) * S - 2.2, (y - LY / 2) * S - 0.25, 0.])


class FnSweep(Slide):
    def construct(self):
        self.camera.background_color = ManimColor(BG)
        title = Text('Доля радиальных гидридов', font=FONT, weight=BOLD).scale(0.62)
        title.to_corner(UL, buff=0.45)

        box = Rectangle(width=LX * S, height=LY * S, stroke_color=GREY_B,
                        stroke_width=2).move_to(to_scene(LX / 2, LY / 2))
        arrows = VGroup(*[
            VGroup(Arrow(box.get_left() + UP * dy, box.get_left() + UP * dy + LEFT * 0.9,
                         buff=0, color=WHITE, stroke_width=5),
                   Arrow(box.get_right() + UP * dy, box.get_right() + UP * dy + RIGHT * 0.9,
                         buff=0, color=WHITE, stroke_width=5))
            for dy in (-1.6, 0., 1.6)])
        sigma = MathTex(r'\sigma_{\theta}', font_size=40).next_to(arrows[2][1], UP, buff=0.15)
        r_ax = VGroup(Arrow(ORIGIN, UP * 0.9, buff=0, color=GREY_A, stroke_width=4),
                      MathTex('r', font_size=34, color=GREY_A))
        r_ax[1].next_to(r_ax[0], RIGHT, buff=0.08)
        r_ax.next_to(box.get_corner(UL), DR, buff=0.18)

        plates = VGroup()
        for c in C:
            p = RoundedRectangle(width=0.6 * S, height=0.09 * S, corner_radius=0.045 * S,
                                 stroke_width=0, fill_color=TAN, fill_opacity=1)
            p.move_to(to_scene(c[0], c[1]))
            plates.add(p)

        self.play(FadeIn(title, shift=DOWN * 0.2), Create(box), run_time=1.2)
        self.play(LaggedStart(*[GrowFromCenter(p) for p in plates], lag_ratio=0.08),
                  run_time=1.8)
        self.play(LaggedStart(*[GrowArrow(a) for pair in arrows for a in pair], lag_ratio=0.05),
                  Write(sigma), FadeIn(r_ax), run_time=1.2)
        self.next_slide()

        # the definition, on the right
        k = ValueTracker(0)
        frac = MathTex(r'F_n', '=', r'{N_\text{рад}', r'\over', r'N}', font_size=46)
        frac.to_edge(RIGHT, buff=1.0).shift(UP * 1.6)
        num = always_redraw(lambda: MathTex(
            r'{%d \over 12}' % round(k.get_value()), r'=', font_size=46,
        ).next_to(frac, DOWN, buff=0.45).align_to(frac, LEFT))
        pct = always_redraw(lambda: DecimalNumber(
            100 * k.get_value() / 12, num_decimal_places=0, unit=r'\,\%',
            font_size=46, color=YELLOW).next_to(num, RIGHT, buff=0.2))
        legend = VGroup(
            VGroup(Square(0.22, fill_color=RAD, fill_opacity=1, stroke_width=0),
                   Text('радиальные', font=FONT).scale(0.38)).arrange(RIGHT, buff=0.15),
            VGroup(Square(0.22, fill_color=TAN, fill_opacity=1, stroke_width=0),
                   Text('тангенциальные', font=FONT).scale(0.38)).arrange(RIGHT, buff=0.15),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.18).next_to(num, DOWN, buff=0.7).align_to(frac, LEFT)
        bar_bg = Rectangle(width=3.2, height=0.16, stroke_color=GREY_B, stroke_width=1.5)
        bar_bg.next_to(legend, DOWN, buff=0.6).align_to(frac, LEFT)
        bar = always_redraw(lambda: Rectangle(
            width=max(1e-3, 3.2 * k.get_value() / 12), height=0.16, stroke_width=0,
            fill_color=YELLOW, fill_opacity=1).align_to(bar_bg, LEFT).align_to(bar_bg, DOWN))
        self.play(Write(frac), FadeIn(num), FadeIn(pct), FadeIn(legend), Create(bar_bg),
                  FadeIn(bar), run_time=1.5)
        self.next_slide()

        # plates turn radial one by one, in the series' nested order;
        # a pause at the points of the computed series (0/20/40/60/80/100 %)
        stops = {2, 5, 7, 10, 12}
        for j, i in enumerate(ORDER, start=1):
            p = plates[i]
            start = p.copy()

            def turn(m, a, start=start):
                m.become(start.copy().rotate(a * PI / 2)
                         .set_fill(interpolate_color(ManimColor(TAN), ManimColor(RAD), a)))

            self.play(UpdateFromAlphaFunc(p, turn), k.animate.set_value(j),
                      Flash(p.get_center(), color=RAD, line_length=0.15, flash_radius=0.35,
                            run_time=0.6),
                      run_time=0.7 if j not in stops else 1.0, rate_func=smooth)
            if j in stops:
                self.next_slide()

        note = Text('то же расположение — меняется только ориентация', font=FONT,
                    slant=ITALIC, color=GREY_A).scale(0.36).next_to(box, DOWN, buff=0.3)
        self.play(FadeIn(note, shift=UP * 0.1))
        self.wait(1)
