"""Looping GIF: stress reorientation of hydrides in a tube wall.

Heat - the hydrides dissolve; apply hoop stress and cool - they precipitate
radially and stay radial when the load is off.  The loop returns to the start
by an explicit cross-fade labelled as the initial structure (no claim about
stress-free re-cooling).  Qualitative: no temperatures or stresses printed.

Layout: CONTENT = [square cell for the ring | column of: T graph, sigma graph,
phase caption, legend]; graph width = column width - axis-label width.

    manim -r 1600,900 --fps 25 reorient.py Reorient
"""
from style import *

# (time s, temperature 0..1, stress 0..1)
SCHED = [(0, 0, 0), (1.2, 0, 0), (3.2, 1, 0), (3.6, 1, 0), (4.2, 1, 1), (4.8, 1, 1),
         (8.0, 0, 1), (8.5, 0, 1), (9.0, 0, 0), (12.5, 0, 0)]
ts = np.array([s[0] for s in SCHED])
T_END = SCHED[-1][0]
RESET = 1.2


def temp(t):
    return float(np.interp(t, ts, [s[1] for s in SCHED]))


def stress(t):
    return float(np.interp(t, ts, [s[2] for s in SCHED]))


# each caption stays >= 1 s
PHASES = [(1.2, 'исходная структура', TAN), (3.6, 'гидриды растворяются', HOT),
          (4.8, 'нагрузка', ACC), (8.5, 'охлаждение под нагрузкой', ACC),
          (9.5, 'нагрузка снята', WHITE), (99, 'радиальные гидриды остаются', RAD)]


def phase(t):
    for t1, s, c in PHASES:
        if t < t1:
            return s, c


class Reorient(Scene):
    def construct(self):
        self.camera.background_color = BG
        rng = np.random.default_rng(3)
        tt, wA, cur_op = ValueTracker(0.), ValueTracker(0.), ValueTracker(1.)

        # ---------------- grid ----------------
        captions = [txt(s, L2, c) for _, s, c in PHASES] + [txt('исходная структура', L2, TAN)]
        cap_w = max(m.width for m in captions)
        cap_h = max(m.height for m in captions)
        leg = legend(L2)
        gap = 0.6
        side = min(CONTENT.h, CONTENT.w - gap - max(cap_w, leg.width))
        cellL = Box(CONTENT.x0, CONTENT.x0 + side, CONTENT.y0, CONTENT.y1)
        cellR = Box(cellL.x1 + gap, CONTENT.x1, CONTENT.y0, CONTENT.y1)
        rg = 0.42
        free = cellR.h - cap_h - leg.height - 3 * rg
        rT, rS, rC, rL = cellR.stack(0.6 * free, 0.4 * free, cap_h, leg.height, gap=rg)
        print('LAYOUT ring side %.2f, column %.2f wide, graphs %.2f / %.2f high' % (
            side, cellR.w, rT.h, rS.h))

        # ---------------- ring ----------------
        C = cellL.c
        r2 = side / 2 - 0.42
        r1 = 0.5 * r2
        ring = Annulus(inner_radius=r1, outer_radius=r2, fill_color=WALL, fill_opacity=1,
                       stroke_color=GREY_B, stroke_width=3).move_to(C)
        pos, tsol, tries = [], [], 0
        while len(pos) < 48 and tries < 200000:
            tries += 1
            r = rng.uniform(r1 + 0.22, r2 - 0.22)
            a = rng.uniform(0, 2 * PI)
            p = C + r * np.array([np.cos(a), np.sin(a), 0.])
            if all(np.linalg.norm(p - q) > 0.40 for q in pos):
                pos.append(p)
                tsol.append(rng.uniform(0.55, 0.95))     # solvus, in units of the peak T
        tsol = np.array(tsol)
        N = len(pos)

        def make_set(p_radial):
            radial = np.zeros(N, bool)
            radial[rng.permutation(N)[:int(round(p_radial * N))]] = True
            return np.where(radial, 90 + rng.normal(0, 10, N), rng.normal(0, 10, N)), radial
        setA, radA = make_set(0.10)
        setB, radB = make_set(0.45)

        def bars(ang, T, op):
            g = VGroup()
            for i, p in enumerate(pos):
                if T >= tsol[i] or op <= 0.01:
                    continue
                rh = (p - C) / np.linalg.norm(p - C)
                nrm = np.arctan2(rh[1], rh[0]) + np.radians(ang[i])
                d = np.array([-np.sin(nrm), np.cos(nrm), 0.])
                w = np.clip(np.sin(np.radians(ang[i])) ** 2, 0, 1)
                f = np.clip((tsol[i] - T) / 0.07, 0, 1)
                g.add(Line(p - d * 0.19 * f, p + d * 0.19 * f, stroke_width=8,
                           color=interpolate_color(TAN, RAD, w)).set_opacity(op))
            return g

        def hydrides():
            t, w = tt.get_value(), wA.get_value()
            if w > 0:
                return VGroup(bars(setB, 0, 1 - w), bars(setA, 0, w))
            return bars(setB if t >= 3.6 else setA, temp(t), 1.)
        hyd = always_redraw(hydrides)

        def load():
            s = stress(tt.get_value())
            g = VGroup()
            if s < 1e-3:
                return g
            for k in range(8):
                a = k * PI / 4 + PI / 8
                u = np.array([np.cos(a), np.sin(a), 0.])
                tng = np.array([-u[1], u[0], 0.])
                base = C + u * (r2 + 0.2)
                for sg in (1, -1):
                    g.add(Arrow(base, base + sg * tng * 0.55, buff=0, color=ACC, stroke_width=6,
                                max_tip_length_to_length_ratio=0.38).set_opacity(s))
                g.add(Arrow(C + u * 0.4, C + u * (r1 - 0.08), buff=0, color=GREY_A,
                            stroke_width=5, max_tip_length_to_length_ratio=0.3).set_opacity(0.85 * s))
            return g
        arrows = always_redraw(load)

        # ---------------- graphs ----------------
        labT = tex('T', L1, HOT)
        labS = tex(r'\sigma_\theta', L1, ACC)
        lab_w = max(labT.width, labS.width) + 0.3
        gx0 = cellR.x0 + lab_w
        gw = cellR.x1 - gx0 - 0.45                      # room for the time-axis tip + t

        def axes(box):
            ax = Axes(x_range=[0, T_END, 1], y_range=[0, 1.15, 1], x_length=gw,
                      y_length=box.h - 0.12, tips=True,
                      axis_config=dict(color=GREY_B, stroke_width=3, include_ticks=False,
                                       tip_width=0.2, tip_height=0.2))
            ax.move_to(np.array([gx0 + gw / 2, box.c[1] - 0.06, 0.]))
            return ax
        axT, axS = axes(rT), axes(rS)
        labT.next_to(axT.y_axis, LEFT, buff=0.3).align_to(axT.y_axis, UP)
        labS.next_to(axS.y_axis, LEFT, buff=0.3).align_to(axS.y_axis, UP)
        tT = tex('t', L3M, GREY_A).next_to(axT.x_axis.get_end(), DOWN, buff=0.12)
        tS = tex('t', L3M, GREY_A).next_to(axS.x_axis.get_end(), DOWN, buff=0.12)
        Tcurve = VMobject(stroke_color=HOT, stroke_width=6).set_points_as_corners(
            [axT.c2p(s[0], s[1]) for s in SCHED])
        Scurve = VMobject(stroke_color=ACC, stroke_width=6).set_points_as_corners(
            [axS.c2p(s[0], s[2]) for s in SCHED])
        cursor = always_redraw(lambda: VGroup(
            DashedLine(axT.c2p(tt.get_value(), 0), axT.c2p(tt.get_value(), 1.12),
                       stroke_width=2, color=GREY_C),
            DashedLine(axS.c2p(tt.get_value(), 0), axS.c2p(tt.get_value(), 1.12),
                       stroke_width=2, color=GREY_C),
            Dot(axT.c2p(tt.get_value(), temp(tt.get_value())), radius=0.13, color=HOT),
            Dot(axS.c2p(tt.get_value(), stress(tt.get_value())), radius=0.13, color=ACC),
        ).set_opacity(cur_op.get_value()))

        def cap():
            if wA.get_value() > 0:
                m = captions[-1].copy()
            else:
                m = captions[[s for _, s, _ in PHASES].index(phase(tt.get_value())[0])].copy()
            return place(m, rC, LEFT)
        caption = always_redraw(cap)
        place(leg, rL, LEFT)

        report({'ring': VGroup(ring, *[Dot(C + r2 * np.array([np.cos(a), np.sin(a), 0.]))
                                        for a in np.linspace(0, 2 * PI, 9)]),
                'axT': VGroup(axT, labT, tT), 'axS': VGroup(axS, labS, tS),
                'caption': place(max(captions, key=lambda m: m.width).copy(), rC, LEFT),
                'legend': leg})

        self.add(ring, hyd, arrows, axT, axS, labT, labS, tT, tS, Tcurve, Scurve, cursor,
                 caption, leg)
        self.play(tt.animate.set_value(T_END), run_time=T_END, rate_func=linear)
        self.play(cur_op.animate.set_value(0.), wA.animate.set_value(1.), run_time=RESET * 0.6)
        tt.set_value(0.)
        wA.set_value(0.)
        self.play(cur_op.animate.set_value(1.), run_time=RESET * 0.4 - 0.04)
