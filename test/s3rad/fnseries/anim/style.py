"""Shared look and layout of the presentation GIFs.

One GIF per slide, 16:9 (1600x900), black.  The top band is left empty for
the slide's own title.  Everything else is placed from measured bounding
boxes: a scene splits CONTENT into cells and fits each element into its cell
with fit(); report() lists overlaps and anything outside CONTENT.

Type scale, measured at 1600x900 (900 px = 7.5 in = 540 pt, 1 px = 0.6 pt):
  L1  MathTex 64          cap 51 px   key formula / number
  L2  Text 30, MathTex 44 cap ~34 px  labels (~20 pt on the slide)
  L3  Text 26, MathTex 38 cap ~29 px  secondary - nothing smaller
"""
import numpy as np
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
ACC = ManimColor('#FFD24A')    # accent
HOT = ManimColor('#FF8A3D')    # temperature
WALL = ManimColor('#1A2430')   # tube wall fill

L1, L2M, L3M = 64, 44, 38      # MathTex sizes
L2, L3 = 30, 26                # Text sizes
STROKE = 4                     # default line width of drawings


class Box:
    """An axis-aligned cell of the frame, in scene units."""
    def __init__(self, x0, x1, y0, y1):
        self.x0, self.x1, self.y0, self.y1 = x0, x1, y0, y1

    w = property(lambda s: s.x1 - s.x0)
    h = property(lambda s: s.y1 - s.y0)
    c = property(lambda s: np.array([(s.x0 + s.x1) / 2, (s.y0 + s.y1) / 2, 0.]))

    def cols(self, *fr, gap=0.5):
        """Split left to right; fr are fractions of the width left after the gaps."""
        free = self.w - gap * (len(fr) - 1)
        out, x = [], self.x0
        for f in fr:
            out.append(Box(x, x + free * f, self.y0, self.y1))
            x += free * f + gap
        return out

    def rows(self, *fr, gap=0.3):
        """Split top to bottom."""
        free = self.h - gap * (len(fr) - 1)
        out, y = [], self.y1
        for f in fr:
            out.append(Box(self.x0, self.x1, y - free * f, y))
            y -= free * f + gap
        return out

    def stack(self, *hs, gap=0.3):
        """Split top to bottom into rows of the given absolute heights."""
        out, y = [], self.y1
        for h in hs:
            out.append(Box(self.x0, self.x1, y - h, y))
            y -= h + gap
        return out

    def inset(self, d):
        return Box(self.x0 + d, self.x1 - d, self.y0 + d, self.y1 - d)

    def __repr__(self):
        return 'Box(x %.2f..%.2f, y %.2f..%.2f)' % (self.x0, self.x1, self.y0, self.y1)


FRAME = Box(-config.frame_width / 2, config.frame_width / 2, -4., 4.)
CONTENT = Box(-6.71, 6.71, -3.65, 2.85)


def fit(m, box, align=ORIGIN, grow=True):
    """Scale m uniformly to the largest size that fits box (only down unless
    grow), then place it at box's centre or edge given by align."""
    k = min(box.w / max(m.width, 1e-6), box.h / max(m.height, 1e-6))
    if grow or k < 1:
        m.scale(k)
    target = box.c + np.array([align[0] * box.w / 2, align[1] * box.h / 2, 0.])
    m.move_to(target, aligned_edge=align)
    return m


def place(m, box, align=ORIGIN):
    """Place m without scaling."""
    target = box.c + np.array([align[0] * box.w / 2, align[1] * box.h / 2, 0.])
    return m.move_to(target, aligned_edge=align)


def bbox(m):
    return np.array([m.get_left()[0], m.get_right()[0], m.get_bottom()[1], m.get_top()[1]])


def report(named, box=CONTENT, tol=0.02):
    """Print overlaps between named mobjects and anything outside box."""
    items = list(named.items())
    bad = 0
    for n, m in items:
        b = bbox(m)
        if b[0] < box.x0 - tol or b[1] > box.x1 + tol or b[2] < box.y0 - tol or b[3] > box.y1 + tol:
            print('LAYOUT outside: %-12s x %.2f..%.2f y %.2f..%.2f' % (n, *b))
            bad += 1
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            a, b = bbox(items[i][1]), bbox(items[j][1])
            ox = min(a[1], b[1]) - max(a[0], b[0])
            oy = min(a[3], b[3]) - max(a[2], b[2])
            if ox > tol and oy > tol:
                print('LAYOUT overlap: %s / %s (%.2f x %.2f)' % (items[i][0], items[j][0], ox, oy))
                bad += 1
    print('LAYOUT %d problem(s)' % bad)


def txt(s, size=L2, color=WHITE, **kw):
    return Text(s, font=FONT, font_size=size, color=color, **kw)


def tex(s, size=L2M, color=WHITE, **kw):
    return MathTex(s, font_size=size, color=color, **kw)


class DecimalRU(DecimalNumber):
    """DecimalNumber with a decimal comma."""
    def _get_num_string(self, number):
        return super()._get_num_string(number).replace('.', ',')


def legend(size=L2):
    """Radial / tangential swatches with their names (swatches share one
    0.62 x 0.62 slot so the names line up)."""
    def row(w, h, col, name):
        slot = Square(0.62, stroke_width=0, fill_opacity=0)
        sw = RoundedRectangle(width=w, height=h, corner_radius=0.06, stroke_width=0,
                              fill_color=col, fill_opacity=1).move_to(slot)
        return VGroup(VGroup(slot, sw), txt(name, size, col)).arrange(RIGHT, buff=0.3)
    return VGroup(row(0.14, 0.62, RAD, 'радиальные'),
                  row(0.62, 0.14, TAN, 'тангенциальные')).arrange(DOWN, aligned_edge=LEFT, buff=0.3)
