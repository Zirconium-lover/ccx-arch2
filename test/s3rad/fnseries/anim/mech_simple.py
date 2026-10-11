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


ZONE_COLS = [ManimColor('#5E5332'), ManimColor('#8A7232'), ManimColor('#B8902E'), ManimColor('#E0B23A')]


def butterfly(n=90):
    """Plane-strain mode-I plastic zone, r(theta) ~ cos^2(theta/2) ((1-2 nu)^2
    + 3 sin^2(theta/2)), theta from the crack growth direction; max r = 1."""
    th = np.linspace(-PI, PI, n)
    r = np.cos(th / 2) ** 2 * ((1 - 2 * 0.3) ** 2 + 3 * np.sin(th / 2) ** 2)
    return th, r / r.max()


class Ductile(_Base):
    """Radial plate across the load: a plastic zone (two-lobed, banded like a
    FE contour plot) forms at each tip; a cloud of pores nucleates in it and
    grows; the pores on the crack line coalesce and the crack runs to the
    edges, the zone moving with its tip."""
    def construct(self):
        PW, PH = 0.3, 2.0
        A = PH / 2
        self.setup_picture()
        plate = RoundedRectangle(width=PW, height=PH, corner_radius=PW / 2, stroke_width=0,
                                 fill_color=RAD, fill_opacity=1).scale(self.sc).move_to(self.sh)
        rng = np.random.default_rng(11)
        T_RUN = (2.4, 6.2)                          # the front goes from the tip to the edge
        th, rb = butterfly()

        def t_front(y):                             # when the front reaches |y|
            return T_RUN[0] + (abs(y) - A) / (L - A) * (T_RUN[1] - T_RUN[0])
        sides = []
        for sg in (1, -1):
            ys = A + 0.2 + np.arange(12) * 0.22
            ys = ys[ys < L - 0.08]
            xs = np.cumsum(rng.normal(0, 0.045, len(ys)))
            chain = [(x, sg * y) for x, y in zip(xs, ys)]
            cloud = []
            tries = 0
            while len(cloud) < 16 and tries < 4000:
                tries += 1
                y = rng.uniform(A + 0.08, L - 0.08)
                x = rng.normal(np.interp(y, ys, xs), 0.2)
                if abs(x) > L - 0.1:
                    continue
                p = np.array([x, sg * y])
                if all(np.hypot(*(p - np.array(q[:2]))) > 0.13 for q in cloud + chain):
                    cloud.append((x, sg * y, rng.uniform(0.022, 0.055),
                                  t_front(y) - rng.uniform(0.5, 0.9)))
            sides.append((sg, chain, cloud))

        def picture():
            t = self.clk.get_value()
            g, te = dip(t), te_of(t)
            load = ramp(te, 0.4, 2.0)
            out = VGroup()
            for sg, chain, cloud in sides:
                k = sum(1 for (x, y) in chain if te >= t_front(y))
                # current front: the last merged pore, or the plate tip
                fx, fy = (0., sg * A) if k == 0 else chain[k - 1]
                done = k == len(chain)
                # plastic zone around the front, banded
                r0 = 0.85 * (0.45 + 0.55 * load) * (1 - ramp(te, T_RUN[1] - 0.4, T_RUN[1] + 0.2))
                if r0 > 0.02:
                    for kk, col in zip((1.0, 0.74, 0.5, 0.3), ZONE_COLS):
                        pts = [self.P(fx + r0 * kk * r * np.sin(a), fy + sg * r0 * kk * r * np.cos(a))
                               for a, r in zip(th, rb)]
                        out.add(Intersection(Polygon(*pts), self.blk, stroke_width=0, fill_color=col,
                                             fill_opacity=0.9 * min(1., 3 * load)))
                # cloud pores: nucleate, grow (a little wider along the load)
                for x, y, rm, tn in cloud:
                    a = ramp(te, tn, tn + 1.5)
                    if a > 0:
                        r = rm * (0.3 + 0.7 * a)
                        out.add(Ellipse(width=2.6 * r * self.sc, height=2 * r * self.sc, stroke_width=0,
                                        fill_color=BLACK, fill_opacity=1).move_to(self.P(x, y)))
                # pores on the crack line grow, then coalesce into the crack
                for j, (x, y) in enumerate(chain):
                    tf = t_front(y)
                    a = ramp(te, tf - 1.3, tf)
                    if a > 0:
                        r = 0.03 + 0.05 * a
                        out.add(Ellipse(width=2.4 * r * self.sc, height=2 * r * self.sc, stroke_width=0,
                                        fill_color=BLACK, fill_opacity=1).move_to(self.P(x, y)))
                pts = [(0., sg * A)] + chain[:k]
                if done:
                    pts.append((chain[-1][0], sg * L))
                if len(pts) > 1:
                    out.add(VMobject(stroke_color=BLACK, stroke_width=7).set_points_as_corners(
                        [self.P(*p) for p in pts]))
            out.add(plate.copy())
            return fade_all(out, g) if g < 1 else out
        self.run(picture)


class Cleave(_Base):
    """Tangential plate along the load: it stretches with the metal; a thin
    zigzag crack runs across it, then a second one, and the pieces move
    slightly apart."""
    def construct(self):
        PL, PT = 3.0, 0.32
        self.setup_picture()
        rng = np.random.default_rng(5)
        CR = [(0.12, 3.0), (-0.86, 4.1)]            # (x, time) of the two cracks
        GAP, NZ = 0.06, 6
        zig = {xc: np.r_[0., rng.uniform(-0.045, 0.045, NZ - 1), 0.] for xc, _ in CR}
        yz = np.linspace(PT / 2, -PT / 2, NZ + 1)

        def picture():
            t = self.clk.get_value()
            g, te = dip(t), te_of(t)
            e = 0.03 * ramp(te, 0.4, 6.0)
            run = {xc: ramp(te, tc, tc + 0.25) for xc, tc in CR}      # crack across the thickness
            gap = {xc: GAP * ramp(te, tc + 0.25, tc + 1.8) for xc, tc in CR}

            def X(x0):
                return x0 * (1 + e) + sum(gp / 2 * np.sign(x0 - xc) for xc, gp in gap.items())
            out = VGroup()
            cut = sorted(xc for xc in run if run[xc] >= 1)
            # pieces between the cracks that have gone through; outline = right
            # boundary top -> bottom, then left boundary bottom -> top
            edges = [-PL / 2] + cut + [PL / 2]
            r = PT / 2
            for a, b in zip(edges[:-1], edges[1:]):
                if b == PL / 2:
                    cx = X(b) - r
                    rb = [self.P(cx + r * np.cos(q), r * np.sin(q)) for q in np.linspace(PI / 2, -PI / 2, 14)]
                else:
                    xb = X(b - 1e-6)
                    rb = [self.P(xb + dz, y) for dz, y in zip(zig[b], yz)]
                if a == -PL / 2:
                    cx = X(a) + r
                    lb = [self.P(cx + r * np.cos(q), r * np.sin(q)) for q in np.linspace(1.5 * PI, PI / 2, 14)]
                else:
                    xa = X(a + 1e-6)
                    lb = [self.P(xa + dz, y) for dz, y in zip(zig[a][::-1], yz[::-1])]
                out.add(Polygon(*(rb + lb), stroke_width=0, fill_color=TAN, fill_opacity=1))
            # a crack still running across: a thin dark zigzag over the plate
            for xc, _ in CR:
                f = run[xc]
                if 0 < f < 1:
                    m = max(2, int(np.ceil(f * NZ)) + 1)
                    pts = [self.P(X(xc) + dz, y) for dz, y in zip(zig[xc][:m], yz[:m])]
                    out.add(VMobject(stroke_color=BLACK, stroke_width=4).set_points_as_corners(pts))
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
