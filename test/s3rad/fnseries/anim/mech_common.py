"""Shared pieces of the three fracture-mechanism GIFs (illustrative, no FE
mesh): the metal block, the omega colour scale and glow, load arrows, and the
loop clock.  Every scene is a pure function of one clock t in [0, T_LOOP),
so the first and the last frame are the same state.
"""
from style import *

MESH = ManimColor('#1A2430')        # undamaged matrix
EDGE = ManimColor('#33435A')        # mesh lines
IFC = ManimColor('#3CC9A0')         # hydride / matrix interface (old deck: green)
OM_STOPS = [(0.0, MESH), (0.35, ManimColor('#5A4A1E')), (0.7, ManimColor('#E0A526')),
            (0.9, ACC), (1.0, ManimColor('#FFF6D0'))]


def omega_color(w):
    """Dark (omega = 0) -> amber -> yellow -> white-hot (omega = 1)."""
    w = float(np.clip(w, 0, 1))
    for (a, ca), (b, cb) in zip(OM_STOPS[:-1], OM_STOPS[1:]):
        if w <= b:
            return interpolate_color(ca, cb, (w - a) / (b - a))
    return OM_STOPS[-1][1]


def omega_bar(width, height=0.32):
    """Horizontal colour bar for omega, 0 at the left, 1 at the right."""
    n = 48
    g = VGroup(*[Rectangle(width=width / n + 0.005, height=height, stroke_width=0,
                           fill_color=omega_color((i + 0.5) / n), fill_opacity=1)
                 .move_to(np.array([(i + 0.5) * width / n - width / 2, 0., 0.]))
                 for i in range(n)])
    frame = Rectangle(width=width, height=height, stroke_color=GREY_B, stroke_width=2)
    return VGroup(g, frame)


def tri_mesh(L, h):
    """Triangles of a structured union-jack mesh on [-L, L]^2: list of 3x3
    vertex arrays (scene units, z = 0) and their centroids."""
    n = int(round(2 * L / h))
    xs = np.linspace(-L, L, n + 1)
    tris = []
    for i in range(n):
        for j in range(n):
            p00 = np.array([xs[i], xs[j], 0.]); p10 = np.array([xs[i + 1], xs[j], 0.])
            p01 = np.array([xs[i], xs[j + 1], 0.]); p11 = np.array([xs[i + 1], xs[j + 1], 0.])
            if (i + j) % 2 == 0:
                tris += [np.array([p00, p10, p11]), np.array([p00, p11, p01])]
            else:
                tris += [np.array([p00, p10, p01]), np.array([p10, p11, p01])]
    cen = np.array([t.mean(0) for t in tris])
    return tris, cen


def load_arrows(L, n=3, length=0.75, gap=0.12):
    """Tension along x: arrows pointing out of the left and right faces."""
    g = VGroup()
    for y in np.linspace(-L * 0.6, L * 0.6, n):
        for sg in (-1, 1):
            p0 = np.array([sg * (L + gap), y, 0.])
            g.add(Arrow(p0, p0 + sg * RIGHT * length, buff=0, color=WHITE, stroke_width=6,
                        max_tip_length_to_length_ratio=0.3))
    return g


def smoothstep(x):
    x = np.clip(x, 0, 1)
    return x * x * (3 - 2 * x)


def ramp(t, t0, t1):
    """0 before t0, 1 after t1, smooth in between."""
    return smoothstep((t - t0) / (t1 - t0))


METAL = ManimColor('#18212C')       # the metal block


def block(L, color=METAL):
    return RoundedRectangle(width=2 * L, height=2 * L, corner_radius=0.18, stroke_color=GREY_B,
                            stroke_width=3, fill_color=color, fill_opacity=1)


def glow(center, intensity, rx, ry, n=9, clip=None):
    """Plastic zone: nested ellipses, dim amber outside, white-hot core.
    intensity in [0, 1] scales both colour and size; clip (a shape) keeps the
    glow inside the metal."""
    g = VGroup()
    if intensity <= 1e-3:
        return g
    for k in range(n):
        f = 1 - k / n                                  # 1 (outer) .. 1/n (inner)
        e = Ellipse(width=2 * rx * f * (0.4 + 0.6 * intensity),
                    height=2 * ry * f * (0.4 + 0.6 * intensity)).move_to(center)
        col = omega_color(0.35 + 0.65 * intensity * (1 - 0.75 * f))   # amber at the faint end, never dark
        op = (0.18 + 0.5 * (1 - f)) * min(1., 3 * intensity)
        if clip is not None:
            e = Intersection(e, clip)
        g.add(e.set_fill(col, op).set_stroke(width=0))
    return g


def fade_all(g, k):
    """Scale every part's own fill and stroke opacity by k (set_opacity would
    give outline-only shapes a fill)."""
    for m in g.family_members_with_points():          # each part once: an Arrow's own
        m.set_fill(opacity=m.get_fill_opacity() * k, family=False)     # (empty) fill must
        m.set_stroke(opacity=m.get_stroke_opacity() * k, family=False)  # not reach its tip
    return g
