"""Looping GIF from the computed model: what the new fracture surface is made
of, against Fn.  One stacked bar per run at its Fn = k/12: walls left by the
deleted matrix (ductile fracture), by the deleted hydride (cleavage), and the
faces of the debonded interface (Run.fracture_surface at eps_f10).  The bars
rise one after another, Fn growing.

    manim -r 1600,900 --fps 25 --disable_caching res_mech.py ResMech
"""
import os
from mech_common import *
from fndata import Run, SERIES

T_LOOP = 10.0
T_RISE = 0.6                          # the first bar starts
DT, T_BAR = 0.45, 0.9                 # stagger, time a bar takes to rise
T_DIP = (7.8, 8.5, 8.8, 9.5)
DUCT = ManimColor('#8D99A8')          # matrix (steel, as the omega glow of the section GIFs)
KEYS = (('ductile', DUCT, 'вязкий разрыв\nматрицы'), ('cleavage', TAN, 'скол гидрида'),
        ('debond', IFC, 'отслоение'))


class ResMech(Scene):
    def construct(self):
        self.camera.background_color = BG
        runs = [Run(fn) for fn, _ in SERIES]
        F = np.array([r.k_radial / 12 for r in runs])
        share = []
        for r in runs:
            f = r.fracture_surface(r.eps_f10)
            tot = sum(f.values())
            share.append([f[k] / tot for k, _, _ in KEYS])
        share = np.array(share)

        # ---------------- legend (right), measured first ----------------
        rows = []
        for _, col, name in KEYS:
            sw = Square(0.5, stroke_width=0, fill_color=col, fill_opacity=1)
            lab = VGroup(*[txt(w, L2, col if col != DUCT else GREY_A) for w in name.split('\n')]).arrange(
                DOWN, aligned_edge=LEFT, buff=0.15)
            rows.append(VGroup(sw, lab).arrange(RIGHT, buff=0.3))
        leg = VGroup(*rows[::-1]).arrange(DOWN, aligned_edge=LEFT, buff=0.45)   # top row = top segment
        cellR = Box(CONTENT.x1 - leg.width, CONTENT.x1, CONTENT.y0, CONTENT.y1)
        cellL = Box(CONTENT.x0, cellR.x0 - 0.8, CONTENT.y0, CONTENT.y1)

        # ---------------- bars over an Fn axis ----------------
        lF = tex(r'F_n', L1)
        t0, t1 = tex('0', L3M, GREY_A), tex('1', L3M, GREY_A)
        y_base = cellL.y0 + t0.height + 0.3
        H = cellL.y1 - 0.1 - y_base
        x_len = cellL.w - lF.width - 0.3 - 0.6
        BW = 0.62 * x_len / 12
        x_of = lambda f: cellL.x0 + 0.4 + f * (x_len - 0.8)
        axis = Arrow(np.array([cellL.x0, y_base, 0.]), np.array([cellL.x0 + x_len, y_base, 0.]), buff=0,
                     color=GREY_B, stroke_width=3, max_tip_length_to_length_ratio=0.03)
        lF.next_to(axis, RIGHT, buff=0.15)
        t0.next_to(np.array([x_of(0), y_base, 0.]), DOWN, buff=0.15)
        t1.next_to(np.array([x_of(1), y_base, 0.]), DOWN, buff=0.15)
        place(leg, Box(cellR.x0, cellR.x1, y_base, y_base + H), LEFT)

        clk = ValueTracker(0.)

        def dip(t):
            if t < T_DIP[0]:
                return 1.
            return 1 - ramp(t, T_DIP[0], T_DIP[1]) + ramp(t, T_DIP[2], T_DIP[3])

        def te_of(t):
            return t if t < T_DIP[1] + 0.15 else 0.

        def picture():
            t = clk.get_value()
            te = te_of(t)
            g = VGroup()
            for i, (f, sh) in enumerate(zip(F, share)):
                h = ramp(te, T_RISE + i * DT, T_RISE + i * DT + T_BAR)   # the bar rises bottom up
                if h <= 0:
                    continue
                y0, top = 0., h
                for s_, (_, col, _) in zip(sh, KEYS):
                    a, b = y0, min(y0 + s_, top)
                    if b > a + 1e-4:
                        g.add(Rectangle(width=BW, height=(b - a) * H, stroke_width=0, fill_color=col,
                                        fill_opacity=1).move_to(np.array([x_of(f), y_base + (a + b) / 2 * H, 0.])))
                    y0 += s_
            k = dip(t)
            return fade_all(g, k) if k < 1 else g

        pic = VGroup()
        pic.add_updater(lambda m: m.set_submobjects(list(picture())))
        report({'bars': VGroup(axis, lF, t0, t1, Rectangle(width=x_len, height=H).move_to(
            np.array([cellL.x0 + x_len / 2, y_base + H / 2, 0.]))), 'legend': leg})

        self.add(pic, axis, lF, t0, t1, leg)
        if os.environ.get('STILL'):
            clk.set_value(float(os.environ['STILL']))
            self.update_mobjects(0)
            self.wait(0.04)
            return
        self.play(clk.animate.set_value(T_LOOP), run_time=T_LOOP, rate_func=linear)
