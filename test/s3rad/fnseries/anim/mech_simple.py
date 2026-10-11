"""Three simple square looping GIFs of the fracture mechanisms, in the spirit
of the old slide 7: the picture only, no bars or laws (the slide text
explains).  Render square:

    manim -r 900,900 --fps 25 mech_simple.py Ductile   # voids coalesce, crack runs
    manim -r 900,900 --fps 25 mech_simple.py Cleave    # the hydride plate cleaves
    manim -r 900,900 --fps 25 mech_simple.py Debond    # the interface opens

Each scene is a pure function of one clock t in [0, T_LOOP); the loop dips to
black and comes back to the untouched state.
"""
import os
from mech_common import *

config.frame_height = 8.0                # square frame for the -r 900,900 render
config.frame_width = 8.0

T_LOOP = 9.0
T_DIP = (7.3, 7.9, 8.1, 8.7)
L = 2.6                                  # block half-size
SQ = Box(-3.7, 3.7, -3.7, 3.7)           # square frame minus margins
STEEL = ManimColor('#34404F')            # lighter metal: black cracks and gaps must read


def dip(t):
    if t < T_DIP[0]:
        return 1.
    return 1 - ramp(t, T_DIP[0], T_DIP[1]) + ramp(t, T_DIP[2], T_DIP[3])


def te_of(t):
    return t if t < T_DIP[1] + 0.1 else 0.


class _Base(Scene):
    def setup_picture(self, *extra):
        self.camera.background_color = BG
        self.blk = block(L, STEEL)
        self.arrows = load_arrows(L, length=0.6)
        whole = VGroup(self.blk, self.arrows, *extra)
        w0 = whole.width
        fit(whole, SQ)
        self.sc = whole.width / w0
        self.sh = self.blk.get_center()
        self.clk = ValueTracker(0.)

    def P(self, x, y):
        return self.sh + np.array([x, y, 0.]) * self.sc

    def run(self, picture):
        pic = always_redraw(picture)
        blk_d = always_redraw(lambda: fade_all(self.blk.copy(), dip(self.clk.get_value())))
        self.add(blk_d, pic, self.arrows)
        if os.environ.get('STILL'):
            self.clk.set_value(float(os.environ['STILL']))
            self.update_mobjects(0)
            self.wait(0.04)
            return
        self.play(self.clk.animate.set_value(T_LOOP), run_time=T_LOOP, rate_func=linear)


class Ductile(_Base):
    """Radial plate across the load: the metal glows at the tips, voids
    nucleate ahead of them, grow and coalesce, the crack runs to the edges."""
    def construct(self):
        PW, PH = 0.3, 2.0
        plate = RoundedRectangle(width=PW, height=PH, corner_radius=PW / 2, stroke_width=0,
                                 fill_color=RAD, fill_opacity=1)
        self.setup_picture()
        plate.scale(self.sc).move_to(self.sh)
        rng = np.random.default_rng(4)
        sites = []
        for sg in (1, -1):
            ys = PH / 2 + 0.25 + np.arange(5) * 0.3
            xs = np.cumsum(rng.normal(0, 0.06, len(ys)))
            sites.append([(x, sg * y) for x, y in zip(xs, ys) if y < L - 0.1])
        n = len(sites[0])
        t_join = np.linspace(2.2, 5.6, n)            # void j joins the crack
        t_nuc = t_join - 1.4                          # void j appears

        def picture():
            t = self.clk.get_value()
            g, te = dip(t), te_of(t)
            out = VGroup()
            I = ramp(te, 0.4, 2.0)
            for sg, ss in zip((1, -1), sites):
                k = int(np.sum(te >= t_join))
                front = (0., sg * PH / 2) if k == 0 else ss[k - 1]
                if k < n:                              # plastic zone ahead of the front
                    out.add(glow(self.P(front[0], front[1] + sg * 0.3), I, 0.5 * self.sc,
                                 0.7 * self.sc, clip=self.blk))
                pts = [(0., sg * PH / 2)] + ss[:k]
                if k == n:
                    pts.append((ss[-1][0], sg * L))
                if len(pts) > 1:
                    out.add(VMobject(stroke_color=BLACK, stroke_width=8 + 4 * k / n)
                            .set_points_as_corners([self.P(*p) for p in pts]))
                for j, (x, y) in enumerate(ss):        # voids grow and stretch toward the crack
                    a = ramp(te, t_nuc[j], t_join[j])
                    if a <= 0 or te >= t_join[j]:
                        continue
                    out.add(Ellipse(width=(0.06 + 0.16 * a) * self.sc,
                                    height=(0.05 + 0.24 * a ** 2) * self.sc, stroke_width=0,
                                    fill_color=BLACK, fill_opacity=1).move_to(self.P(x, y)))
            out.add(plate.copy())
            return fade_all(out, g) if g < 1 else out
        self.run(picture)


class Cleave(_Base):
    """Tangential plate along the load: it stretches with the metal and
    cleaves across in a few places; the pieces move apart."""
    def construct(self):
        PL, PT = 3.0, 0.32
        self.setup_picture()
        XC = [0.06, -0.82, 0.86]
        t_c = [3.0, 3.25, 3.5]
        GAP = 0.12

        def picture():
            t = self.clk.get_value()
            g, te = dip(t), te_of(t)
            e = 0.05 * ramp(te, 0.4, 6.0)
            out = VGroup()
            opens = {x: GAP * ramp(te, tc, tc + 2.0) for x, tc in zip(XC, t_c) if te >= tc}

            def X(x0):
                return x0 * (1 + e) + sum(o / 2 * np.sign(x0 - xc) for xc, o in opens.items())
            edges = [-PL / 2] + sorted(opens) + [PL / 2]
            whole = len(opens) == 0
            for a, b in zip(edges[:-1], edges[1:]):
                xa, xb = X(a + 1e-6), X(b - 1e-6)
                r = PT / 2 if whole else 0.04
                out.add(RoundedRectangle(width=(xb - xa) * self.sc, height=PT * self.sc,
                                         corner_radius=r * self.sc, stroke_width=0, fill_color=TAN,
                                         fill_opacity=1).move_to(self.P((xa + xb) / 2, 0)))
            for xc, tc in zip(XC, t_c):                 # a short flash as each crack appears
                f = 1 - ramp(te, tc, tc + 0.35) if te >= tc else 0.
                if f > 0.01:
                    out.add(Line(self.P(X(xc), -PT * 1.3), self.P(X(xc), PT * 1.3),
                                 stroke_color=WHITE, stroke_width=5).set_stroke(opacity=f))
            return fade_all(out, g) if g < 1 else out
        self.run(picture)


class Debond(_Base):
    """Radial plate across the load: the interface comes apart along the
    faces from the middle outward; a thin gap opens and the plate works as
    a crack."""
    def construct(self):
        PW, PH = 0.3, 2.4
        A = PH / 2
        self.setup_picture()
        plate = RoundedRectangle(width=PW, height=PH, corner_radius=PW / 2, stroke_width=0,
                                 fill_color=RAD, fill_opacity=1).scale(self.sc).move_to(self.sh)
        DMAX = 0.15                                     # final opening of each face - modest

        def picture():
            t = self.clk.get_value()
            g, te = dip(t), te_of(t)
            front = A * ramp(te, 1.0, 5.0)              # debonded half-length
            dm = DMAX * ramp(te, 1.0, 5.6)
            out, bonds = VGroup(), VGroup()
            # the bonded interface: a thin green line on both faces where still intact
            for sg in (1, -1):
                xf = sg * PW / 2
                for y0, y1 in ((front, A - PW / 2), (-(A - PW / 2), -front)):
                    if y1 - y0 > 0.02:
                        bonds.add(Line(self.P(xf, y0), self.P(xf, y1), stroke_color=IFC, stroke_width=4))
                if front > 0.02 and dm > 1e-3:          # the gap on the debonded part
                    yy = np.linspace(-front, front, 40)
                    dd = dm * np.sqrt(np.clip(1 - (yy / front) ** 2, 0, 1))
                    pts = [self.P(xf, y) for y in yy] + \
                          [self.P(xf + sg * d, y) for y, d in zip(yy[::-1], dd[::-1])]
                    out.add(Polygon(*pts, stroke_width=0, fill_color=BLACK, fill_opacity=1))
            out.add(plate.copy(), bonds)
            return fade_all(out, g) if g < 1 else out
        self.run(picture)
