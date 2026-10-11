"""Looping GIF for the problem statement: the computational block in 3D with
its 12 hydride discs (seed 5, the deck geometry: block 4 x 3.6 x 1.6 mm,
discs of radius 0.3 and thickness 0.09).  The discs turn from tangential
(normal along R) to radial (normal along T, the load) one by one in the
series' nested order, F_n = k/12 follows; then they turn back.

Model axes: x = T (hoop, the load), y = R, z = L.  Scene: X = T, Y = R (into
the picture), Z (up) = L - the tube axis stands roughly upright and the discs
turn about it.  The camera sways a little; the axis and load labels follow
the camera projection every frame.  Size and position come from the
projection of the block, arrows and triad over the whole sway.

    manim -r 1600,900 --fps 25 fn_block.py FnBlock
"""
import sys
from style import *

sys.path.insert(0, '/home/user/ccx-fn/test/s3rad')
import mkseeddeck as gen  # noqa: E402

LX, LY, LZ = gen.LX, gen.LY, gen.LZ
CEN = gen.place_any(5)
ORDER = list(np.random.default_rng(10000 + 5).permutation(12))
STEP, HOLD0, HOLD1 = 0.42, 1.2, 1.6
PHI, TH0, SWAY = 56 * DEGREES, -62 * DEGREES, 0 * DEGREES   # no sway: steady labels, small GIF
T_LOOP = HOLD0 + 12 * STEP + HOLD1 + 12 * STEP + 0.04


def fn_tex(k):
    v = '0' if k == 0 else '1' if k == 12 else ('%.2f' % (k / 12)).replace('.', '{,}')
    return r'F_n = {%d \over 12} = %s' % (k, v)


def m2s(p):
    """Model point (x=T, y=R, z=L), mm -> scene (X=T, Y=R, Z=L up), centred."""
    x, y, z = p
    return np.array([x - LX / 2, y - LY / 2, z - LZ / 2])


class FnBlock(ThreeDScene):
    def fixed(self, *ms):
        self.add_fixed_in_frame_mobjects(*ms)
        return ms

    def construct(self):
        self.camera.background_color = BG
        self.set_camera_orientation(phi=PHI, theta=TH0)
        cam = self.camera

        # ---------------- right column, measured ----------------
        widest = max((tex(fn_tex(k), L1) for k in range(13)), key=lambda m: m.width)
        leg = legend(L2)
        col_w = max(widest.width, leg.width) + 0.1
        left = Box(CONTENT.x0, CONTENT.x1 - col_w - 0.6, CONTENT.y0, CONTENT.y1)
        right = Box(CONTENT.x1 - col_w, CONTENT.x1, CONTENT.y0, CONTENT.y1)

        # ---------------- geometry in scene units (scale s, shift sh) ----------------
        corners = [m2s((x, y, z)) for x in (0, LX) for y in (0, LY) for z in (0, LZ)]
        AR = 0.75                                 # load-arrow length, mm
        extra = [m2s((-AR - 0.05, LY / 2, LZ / 2)), m2s((LX + AR + 0.05, LY / 2, LZ / 2))]
        pts0 = np.array(corners + extra)
        def proj(P, th):
            self.set_camera_orientation(phi=PHI, theta=th)
            cam.reset_rotation_matrix()
            return cam.project_points(np.atleast_2d(P))[:, :2]

        def extent(s, sh):
            q = np.vstack([proj(pts0 * s + sh, TH0 + SWAY * np.sin(a))
                           for a in np.linspace(0, 2 * PI, 13)])
            return q.min(0), q.max(0)
        box = left.inset(0.35)                    # room for the labels around
        s, sh = 1., ORIGIN
        for _ in range(6):
            lo, hi = extent(s, sh)
            s *= 0.98 * min(box.w / (hi - lo)[0], box.h / (hi - lo)[1])
            lo, hi = extent(s, sh)
            d = box.c[:2] - (lo + hi) / 2
            self.set_camera_orientation(phi=PHI, theta=TH0)
            cam.reset_rotation_matrix()
            sh = sh + cam.get_rotation_matrix().T @ np.array([d[0], d[1], 0.])
        lo, hi = extent(s, sh)
        self.set_camera_orientation(phi=PHI, theta=TH0)
        cam.reset_rotation_matrix()
        print('LAYOUT scale %.3f, projected x %.2f..%.2f y %.2f..%.2f (cell %s)'
              % (s, lo[0], hi[0], lo[1], hi[1], left))

        def S(p):
            return m2s(p) * s + sh

        # block edges
        edges = VGroup()
        for a in range(8):
            for b in range(a + 1, 8):
                pa, pb = corners[a], corners[b]
                if np.sum(np.abs(pa - pb) > 1e-9) == 1:
                    edges.add(Line(pa * s + sh, pb * s + sh, stroke_color=GREY_B, stroke_width=3))
        floor = Polygon(*[S(p) for p in ((0, 0, 0), (LX, 0, 0), (LX, LY, 0), (0, LY, 0))],
                        stroke_width=0, fill_color=WALL, fill_opacity=0.6)

        # discs: tangential = normal along R (scene Y), radial = normal along T (scene X)
        discs = []
        for c in CEN:
            d = Circle(radius=gen.RAD * s, stroke_color=interpolate_color(TAN, WHITE, 0.4),
                       stroke_width=2, fill_color=TAN, fill_opacity=0.95)
            d.rotate(PI / 2, axis=RIGHT).move_to(S(c))
            discs.append(d)
        dg = VGroup(*discs)

        # load arrows on the two T faces
        arrows = VGroup()
        for dy in (-1.15, 0., 1.15):
            for x0, sg in ((0., -1), (LX, 1)):
                p0 = S((x0, LY / 2 + dy, LZ / 2))
                arrows.add(Arrow3D(p0, p0 + sg * np.array([AR * s, 0, 0]), color=WHITE,
                                   thickness=0.022, height=0.26, base_radius=0.09))

        sig_at = S((LX + AR + 0.05, LY / 2 + 1.15, LZ / 2))
        # a 2D axis triad in the free corner, drawn along the projected axes
        o2 = np.array([left.x0 + 0.45, left.y0 + 0.8, 0.])        # free corner, bottom left
        TRI = 0.85
        def projected(P):
            cam.reset_rotation_matrix()
            q = cam.project_points(np.atleast_2d(P))[0]
            return np.array([q[0], q[1], 0.])

        labs = VGroup(tex('T', L2M, GREY_A), tex('R', L2M, GREY_A), tex('L', L2M, GREY_A))
        sig = tex(r'\sigma_\theta', L2M)

        def triad():
            cam.reset_rotation_matrix()
            q0 = cam.project_points(np.atleast_2d(sh))[0]
            g = VGroup()
            for e, lab in zip((RIGHT, UP, OUT), labs):     # T = X, R = Y, L = Z
                q = cam.project_points(np.atleast_2d(sh + e))[0] - q0
                u = np.array([q[0], q[1], 0.]); u /= np.linalg.norm(u)
                g.add(Arrow(o2, o2 + u * TRI, buff=0, color=GREY_A, stroke_width=5,
                            max_tip_length_to_length_ratio=0.25))
                lab.move_to(o2 + u * (TRI + 0.28))
            return g
        tri2 = triad()

        def follow(m):
            tri2.become(triad())
            cam.fixed_in_frame_mobjects.update(tri2.get_family())   # new parts stay 2D
            sig.next_to(projected(sig_at), UP, buff=0.12)
        follow(None)

        # ---------------- right column ----------------
        k = ValueTracker(0)
        blk = VGroup(widest.copy(), leg).arrange(DOWN, buff=1.1, aligned_edge=LEFT)
        place(blk, right)
        anchor = blk[0].get_left()
        states = VGroup(*[tex(fn_tex(j), L1).move_to(anchor, aligned_edge=LEFT) for j in range(13)])

        def show_state(m):
            j = int(round(k.get_value()))
            for i, st in enumerate(m):
                st.set_opacity(1. if i == j else 0.)
        show_state(states)
        states.add_updater(show_state)
        report({'block': VGroup(*[Dot(np.array([x, y, 0.]), radius=0.01)
                                  for x, y in ((lo[0], lo[1]), (hi[0], hi[1]))]),
                'triad': VGroup(tri2, labs),
                'formula': blk[0], 'legend': leg})

        self.add(floor, edges, dg, arrows)
        self.fixed(labs, sig, leg, states, tri2)
        labs.add_updater(lambda m: follow(m))

        # camera sway over the whole loop, back where it started
        t_now = ValueTracker(0.)
        self.camera.theta_tracker.add_updater(
            lambda m: m.set_value(TH0 + SWAY * np.sin(2 * PI * t_now.get_value() / T_LOOP)))
        t_now.add_updater(lambda m, dt: m.increment_value(dt))
        self.add(t_now, self.camera.theta_tracker)

        def turn(i, j, to_radial):
            dsk = discs[i]
            start = dsk.copy()
            c0, c1 = (TAN, RAD) if to_radial else (RAD, TAN)
            sgn = 1 if to_radial else -1
            ctr = start.get_center()

            def upd(m, a, start=start):
                m.become(start.copy().rotate(-sgn * a * PI / 2, axis=OUT, about_point=ctr)
                         .set_fill(interpolate_color(c0, c1, a))
                         .set_stroke(interpolate_color(interpolate_color(c0, c1, a), WHITE, 0.4)))
            return AnimationGroup(UpdateFromAlphaFunc(dsk, upd), k.animate.set_value(j),
                                  run_time=STEP, rate_func=smooth)

        self.wait(HOLD0)
        for j, i in enumerate(ORDER, start=1):
            self.play(turn(i, j, True))
        self.wait(HOLD1)
        for j, i in zip(range(11, -1, -1), reversed(ORDER)):
            self.play(turn(i, j, False))
        self.wait(0.04)
