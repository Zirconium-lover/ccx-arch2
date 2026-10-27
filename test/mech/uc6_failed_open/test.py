"""ACCEPTANCE: a failed UC6 facet carries no tension and no shear (CCX_UC6_FAILED_OPEN=1).

The default law keeps a broken facet pulling with g_min Kn delta in tension
and shear.  That floor is a numerical device, and test/mech/size_effect
showed what it costs when g_min Kn is not far below E: at g_min = 1e-3 the
broken ligament of a notched plate carried 1.44 Tn0 and the plate never
fractured.  The switch (developer, ffec13c; interface agreed 27.10):
    floor      removed: g = max(0, 1 - D_vis), so traction and tangent go
               to 0 CONTINUOUSLY, reaching it exactly at delta_f (D = 1);
               with mu > 0 the facet is failed when D_vis = 1
    tension    T_n = 0 and its tangent 0 once failed
    shear      T_s = 0 and its tangent 0, at ANY sign of delta_n (no friction)
    contact    T_n = Kn delta_n, as before
    off        the law is bit-identical to the default
References (closed form, exact on these paths): T = 0 to the printed digit
for delta >= delta_f; the work to delta_f is Gc exactly (the default law adds
the floor's sliver, 1.00059 Gc); compression Kn delta_n; the notched plate of
size_effect at D/l_ch = 1 breaks at the same sigma_N with g_min = 1e-3 as the
default law does with g_min = 1e-6 (0.581 Tn0).

Facet: uc6_law/uc6deck.py, Kn = 1e6, Tn0 = 400, Ts0 = 300, Gc = 2, g_min =
1e-3 (large on purpose, so the floor is visible): delta_0 = 4e-4 and delta_f =
1e-2 lie on the increments, so the trapezoidal work is exact (1e-5).

HISTORY.  The developer's first version opened the facet with a JUMP when
1 - D reached g_min (delta_c = 0.009766 here): the traction fell by
g_min Kn delta (~1 % of Tn0) at once on many facets, and the size_effect
plate failed before its peak (theta 0.302, sigma_N 0.306 Tn0).  This test
then expected W = 0.99943 Gc; the continuous form replaced both.

Measured on e2d474f, a binary WITHOUT the switch: seven checks red, as an
acceptance test must be before its code - residual tractions 0.03 Tn0, W
1.00061 Gc, the plate at 0.976 Tn0 with g_min = 1e-3 - and both control
checks green.  Not in the gate until the developer's commit is merged.

RED PATH: the control arms run the same paths with the switch OFF and must
show the floor (T = g_min Kn delta after failure) - on a binary without the
switch every "with the switch" check fails the same way.
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "uc6_law"))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "size_effect"))
from mechlib import Check, ccx, totals, trapz, write  # noqa: E402
import uc6deck  # noqa: E402

WHAT = "ACCEPTANCE CCX_UC6_FAILED_OPEN: failed facet T = 0 in tension/shear, Kn in contact, W = Gc; plate breaks at gmin=1e-3"
SOURCE = "closed-form bilinear law without residual traction; size_effect plate (sigma_N 0.581 Tn0 at D/l_ch = 1)"

KN, TN0, TS0, GC, GMIN = 1e6, 400.0, 300.0, 2.0, 1e-3
BETA = (TS0 / TN0) ** 2
D0, DF = TN0 / KN, 2 * GC / TN0
AREA = 0.5
SW = {"CCX_UC6_FAILED_OPEN": "1"}


def facet(work, name, path, env, dt=0.01):
    d = os.path.join(work, name)
    write(os.path.join(d, "u.inp"), uc6deck.deck(path, (KN, TN0, TS0, GC, GMIN, 0.0), dt))
    rc, _ = ccx(d, "u", EXE, env)
    dat = os.path.join(d, "u.dat")
    f = [totals(dat, "NPLUS", c) for c in range(3)] if rc == 0 else [[], [], []]
    return rc, [(t, (fx / AREA, fy / AREA, fz / AREA)) for (t, fx), (_, fy), (_, fz) in zip(*f)]


def lerp(path, t):
    for (t0, j0), (t1, j1) in zip(path, path[1:]):
        if t0 - 1e-12 <= t <= t1 + 1e-12:
            s = (t - t0) / (t1 - t0)
            return tuple(j0[c] + s * (j1[c] - j0[c]) for c in range(3))
    return path[-1][1]


def deff(j):
    return math.sqrt(max(j[0], 0) ** 2 + BETA * (j[1] ** 2 + j[2] ** 2))


def work(tr, path, upto_df=True):
    w, pj, pt = 0.0, (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)
    for t, T in tr:
        j = lerp(path, t)
        w += sum(0.5 * (T[c] + pt[c]) * (j[c] - pj[c]) for c in range(3))
        pj, pt = j, T
        if upto_df and deff(j) >= DF * (1 - 1e-12):
            break
    return w


EXE = None


def run(wk, exe):
    global EXE
    EXE = exe
    ck = []
    s = 1 / math.sqrt(BETA)
    p1 = [(0, (0, 0, 0)), (1.2, (0.012, 0, 0))]
    p2 = [(0, (0, 0, 0)), (1.2, (0, 0.012 * s, 0))]
    p3 = [(0, (0, 0, 0)), (1.2, (0.012, 0, 0)), (1.4, (-0.001, 0, 0)), (1.6, (-0.001, 0.006, 0))]
    for tag, env in (("switch", SW), ("default", {})):
        r1, t1 = facet(wk, "modeI_" + tag, p1, env)
        r2, t2 = facet(wk, "modeII_" + tag, p2, env)
        r3, t3 = facet(wk, "contact_" + tag, p3, env)
        ck.append(Check("%s: runs (rc)" % tag, max(abs(r1), abs(r2), abs(r3)), 0, tol=0))
        if r1 or r2 or r3:
            continue
        post1 = [abs(T[0]) for t, T in t1 if deff(lerp(p1, t)) >= DF * (1 - 1e-12)]
        post2 = [abs(T[1]) for t, T in t2 if deff(lerp(p2, t)) >= DF * (1 - 1e-12)]
        shear3 = [abs(T[1]) for t, T in t3 if t > 1.4 + 1e-9]
        tn3 = [T[0] / (KN * lerp(p3, t)[0]) for t, T in t3 if t > 1.3 and lerp(p3, t)[0] < -1e-5]
        if tag == "switch":
            ck.append(Check("mode I after failure: max |T_n| / Tn0", max(post1) / TN0, 0.0, tol=1e-9))
            ck.append(Check("mode II after failure: max |T_s| / Ts0", max(post2) / TS0, 0.0, tol=1e-9))
            ck.append(Check("failed, compressed, sheared: max |T_s| / Ts0 (no friction)", max(shear3) / TS0,
                            0.0, tol=1e-9))
            ck.append(Check("failed, compressed: T_n / (Kn delta_n)", min(tn3), 1.0, tol=1e-5))
            ck.append(Check("failed, compressed: T_n / (Kn delta_n) (max)", max(tn3), 1.0, tol=1e-5))
            ck.append(Check("mode I: W to delta_f / Gc (no g_min sliver)", work(t1, p1) / GC, 1.0, tol=1e-5))
            ck.append(Check("mode II: W to delta_f / Gc", work(t2, p2) / GC, 1.0, tol=1e-5))
            ck.append(Check("mode I: work after failure (to 1.2 delta_f)", work(t1, p1, False) - work(t1, p1), 0.0,
                            tol=1e-9, why="the facet does no work once failed"))
        else:
            ck.append(Check("CONTROL switch off: residual T_n / (g_min Kn delta) at 1.2 delta_f",
                            t1[-1][1][0] / (GMIN * KN * 0.012), 1.0, tol=1e-4,
                            why="the floor the switch removes; also the red path"))
            ck.append(Check("CONTROL switch off: failed+compressed shear traction > 0", max(shear3) / TS0,
                            1e-4, cmp="ge"))
    # the switch together with the mid-surface frame: mode I, one arm
    r4, t4 = facet(wk, "modeI_midsurface", p1, dict(SW, CCX_UC6_MIDSURFACE="1"))
    if r4 == 0:
        ck.append(Check("FAILED_OPEN + MIDSURFACE, mode I: W to delta_f / Gc", work(t4, p1) / GC, 1.0, tol=1e-5))
    else:
        ck.append(Check("FAILED_OPEN + MIDSURFACE: run rc", r4, 0, tol=0))
    # the notched plate of size_effect at D/l_ch = 1 with g_min = 1e-3
    import plate
    d = os.path.join(wk, "plate_gmin1e-3")
    w = 10.0
    ue = 1.4 * 70.0 / 1e5 * 3 * w * 1.15
    write(os.path.join(d, "p.inp"), plate.deck(w, a_over_w=0.3, nx=40, ny=16, uend=ue, dt=1 / 50.0, gmin=1e-3))
    rc, _ = ccx(d, "p", exe, SW, timeout=300)
    c = totals(os.path.join(d, "p.dat"), "NTOP", 1)
    sn = max(f for _, f in c) / (w * w / 40) / 100.0 if c else float("nan")
    ck.append(Check("plate D/l_ch=1, g_min=1e-3, switch: sigma_N/Tn0 vs g_min=1e-6 default (0.5809)", sn, 0.5809,
                    rtol=0.03, why="without the switch the broken ligament carried 1.44 Tn0"))
    return ck
