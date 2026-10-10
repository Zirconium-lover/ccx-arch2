"""Shared look of the presentation GIFs: one GIF per slide, 16:9, black, the
top band left empty for the slide's own title, text readable from the back
of a hall (>= ~20 pt on the slide)."""
from manim import *

config.tex_template = TexTemplate(preamble=r'''\usepackage[T2A]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage[english,russian]{babel}
\usepackage{amsmath}
\usepackage{amssymb}''')

FONT = 'Noto Serif'            # the deck's font
BG = ManimColor('#000000')
RAD = ManimColor('#E8604C')    # radial hydrides
TAN = ManimColor('#2CB1D9')    # tangential hydrides
HYD = ManimColor('#D8D8D8')    # delta hydride (neutral)
ZR = ManimColor('#7FB2D9')     # alpha-Zr lattice
ACC = ManimColor('#FFD24A')    # accent / numbers
HOT = ManimColor('#FF8A3D')    # temperature
WALL = ManimColor('#1A2430')   # tube wall fill

TOP = 2.95                     # content stays below this (title band above)
BOTTOM = -3.85
FS_HUGE, FS_BIG, FS, FS_SM = 96, 66, 50, 42      # MathTex sizes
TX_BIG, TX, TX_SM = 32, 27, 23                    # Text sizes (Pango runs larger)


def txt(s, size=TX, color=WHITE, **kw):
    return Text(s, font=FONT, font_size=size, color=color, **kw)


def tex(s, size=FS, color=WHITE, **kw):
    return MathTex(s, font_size=size, color=color, **kw)


class DecimalRU(DecimalNumber):
    """DecimalNumber with a decimal comma."""
    def _get_num_string(self, number):
        return super()._get_num_string(number).replace('.', ',')
