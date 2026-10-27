"""UC6 cohesive law through the solver, one facet, against the closed form.

The facet is driven by its own nodes (uc6deck.py), so what comes back is the
law as the solver assembles, integrates and reports it: three Gauss points,
area weights, local frame, state update - not the routine called directly,
which an earlier driver already did (W/Gc = 1.001 there).

Closed form (src/cohesive_uc6.f header; Camanho and Davila 2002 bilinear law
with a single Gc and a quadratic effective separation):
    beta = (Ts0/Tn0)^2, d_eff = sqrt(<dn>^2 + beta (ds1^2 + ds2^2))
    d0 = Tn0/Kn, df = 2 Gc/Tn0, D = df (dmax - d0) / (dmax (df - d0))
    Tn = g Kn dn (dn >= 0), Kn dn (dn < 0); Ts = g beta Kn ds; g = max(gmin, 1 - Dv)
    Dv: Dv_n = Dv_{n-1} + dt/(mu+dt) (D_n - Dv_{n-1})   (backward Euler of
        Duvaut-Lions dDv/dt = (D - Dv)/mu)
so on any monotonic path, the work per unit area to d_eff = df is Gc - plus
the sliver the gmin floor adds just before df (0.061 % here), which w_exact()
integrates in closed form.

Parameters Kn=1e6, Tn0=400, Ts0=300 (beta=0.5625, so mode II is NOT a
relabelled mode I), Gc=2, gmin=1e-3: d0=4e-4, df=0.01.  Paths are sampled so
that d0 and df fall on increments, which makes the trapezoidal work exact for
the piecewise-linear law; the tolerances are then the 7 digits the .dat file
prints (1e-5 relative), not a discretisation allowance - except for the work,
whose floor crossing dc = 0.009766 is NOT on the grid: the trapezoid across
that kink is off by at most (slope jump)*dd^2/8 = 42667*1e-8/8 = 5e-5 of
W, measured 2e-5, so the work gets 5e-5.

RED PATH, run: MECH_RED_ENV=CCX_UC6_CONTACT_SMOOTH=1 (the compression slope
blended) turns the four compression checks red (T(-0.001) = 0.8125 Kn d);
MECH_RED_ENV="CCX_UC6_RESCUE_ELEMS=1 CCX_UC6_RESCUE_G=0.01" (the floor raised
to 1e-2 on this facet) turns eleven red: every T(d) history, W +4.78 % on
all four energy paths, the residual traction x10.

Viscosity: mu = 0.05 against a ramp crossing d0..df in 0.96 time units.  The
deck's Dv follows the backward-Euler recursion to the printed digits; against
the continuous Duvaut-Lions ODE (RK4, h=1e-4) the error is 0.031 at dt=0.02
and 0.010 at dt=0.005 - ratio 3.0 for a 4x smaller step, first order with
the kinks of D(d) at d0 and df taking the rest (35 % allowed on the ratio).

Known limitation, kept as XFAIL (acceptance for the co-rotated frame the
developer is writing; set ROT_ENV / ROT_FIXED below when it lands): the frame
and area are the REFERENCE ones, so a rigidly rotated open facet is read in
the old axes.  Checked at 30 / 60 / 90 deg: the force along the turned normal
is 0.891 / 0.672 / 0.5625 of Kn d A, the force across it -0.189 / -0.189 / 0,
and with Ts0 > Tn0 the elastic facet is damaged from 60 deg on (d_eff =
d sqrt(cos^2 + beta sin^2) = 4.17e-4 and 4.38e-4 against d0 = 4e-4).  Turned 90 deg about
z with an elastic opening of 2e-4, the reaction must be Kn*d*A = 100 along
the new normal; the code reports beta*Kn*d*A = 56.25 (a shear reading).  With
Ts0 > Tn0 the same rotation DAMAGES an elastic facet (d_eff = sqrt(beta)*d >
d0), which an objective law cannot do.
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mechlib import Check, ccx, totals, blocks, write  # noqa: E402
import uc6deck  # noqa: E402

WHAT = "UC6 cohesive law in a deck: modes I/II/45deg, Gc, unloading, compression, viscosity, gmin, rotation"
SOURCE = "closed-form bilinear law (cohesive_uc6.f header; Camanho & Davila 2002); Duvaut-Lions ODE"

KN, TN0, TS0, GC, GMIN = 1e6, 400.0, 300.0, 2.0, 1e-3
BETA = (TS0 / TN0) ** 2
D0, DF = TN0 / KN, 2 * GC / TN0
AREA = 0.5


def dback(dmax):
    if dmax <= D0:
        return 0.0
    if dmax >= DF:
        return 1.0
    return DF * (dmax - D0) / (dmax * (DF - D0))


def w_exact():
    """Work per unit area to d_eff = df of the law WITH its gmin floor:
    near df the softening line falls below gmin*Kn*d, and from the crossing
    dc on the facet pulls with the floor.  Exact for the piecewise law."""
    s = TN0 / (DF - D0)                       # softening slope magnitude
    dc = s * DF / (GMIN * KN + s)             # gmin*Kn*d = s*(df - d)
    return (KN * D0 ** 2 / 2 + (TN0 + s * (DF - dc)) / 2 * (dc - D0)
            + GMIN * KN * (DF ** 2 - dc ** 2) / 2)


def law(path_j, mu=0.0, dts=None):
    """Exact traction history for a list of local jumps (one per increment)."""
    out, dmax, dv = [], 0.0, 0.0
    for k, j in enumerate(path_j):
        de = math.sqrt(max(j[0], 0) ** 2 + BETA * (j[1] ** 2 + j[2] ** 2))
        dmax = max(dmax, de)
        db = dback(dmax)
        a = 1.0 if mu == 0 else dts[k] / (mu + dts[k])
        dv = max(dv, min(1.0, dv + a * (db - dv)))
        g = max(GMIN, 1 - dv)
        tn = KN * j[0] if j[0] < 0 else g * KN * j[0]
        out.append((tn, g * BETA * KN * j[1], g * BETA * KN * j[2], dv))
    return out


def lerp(path, t):
    for (t0, j0), (t1, j1) in zip(path, path[1:]):
        if t0 - 1e-12 <= t <= t1 + 1e-12:
            s = (t - t0) / (t1 - t0)
            return tuple(j0[c] + s * (j1[c] - j0[c]) for c in range(3))
    return path[-1][1]


def drive(work, name, path, dt, exe, mu=0.0, rot=None):
    d = os.path.join(work, name)
    write(os.path.join(d, "u.inp"), uc6deck.deck(path, (KN, TN0, TS0, GC, GMIN, mu), dt, rot=rot))
    rc, _ = ccx(d, "u", exe)
    dat = os.path.join(d, "u.dat")
    rf = [(t, (f0,)) for t, f0 in totals(dat, "NPLUS", 0)] if rc == 0 else []
    rfy = dict(totals(dat, "NPLUS", 1)) if rc == 0 else {}
    rfz = dict(totals(dat, "NPLUS", 2)) if rc == 0 else {}
    st = blocks(dat, "stresses") if rc == 0 else {}
    f = [(t, (fx[0], rfy[t], rfz[t])) for t, fx in rf]
    return rc, f, st


def compare(name, path, dt, exe, work, mu=0.0, energy=True):
    rc, f, st = drive(work, name, path, dt, exe, mu)
    ck = [Check("%s: run rc" % name, rc, 0, tol=0)]
    if rc:
        return ck, None
    ts = [t for t, _ in f]
    js = [lerp(path, t) for t in ts]
    dts = [ts[0]] + [b - a for a, b in zip(ts, ts[1:])]
    ref = law(js, mu, dts)
    err = max(abs(F[c] / AREA - R[c]) for (_, F), R in zip(f, ref) for c in range(3))
    ck.append(Check("%s: max |T - T_exact| / Tn0" % name, err / TN0, 0.0, tol=1e-5))
    # the three Gauss points must agree: the jump is uniform over the facet
    gp = max(abs(r[c] - rows[0][c]) for rows in st.values() for r in rows for c in range(2, 8))
    ck.append(Check("%s: Gauss points agree (max diff)" % name, gp, 0.0, tol=1e-6 * TN0))
    w = None
    if energy:
        # work per unit area up to the first instant d_eff reaches df
        w, prev_j, prev_t = 0.0, (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)
        for (t, F), j in zip(f, js):
            T = tuple(x / AREA for x in F)
            w += sum(0.5 * (T[c] + prev_t[c]) * (j[c] - prev_j[c]) for c in range(3))
            prev_j, prev_t = j, T
            if math.sqrt(max(j[0], 0) ** 2 + BETA * (j[1] ** 2 + j[2] ** 2)) >= DF * (1 - 1e-12):
                break
    return ck, (w, f, js, st)


def run(work, exe):
    ck = []
    s = 1 / math.sqrt(BETA)
    # mode I: d0 at t=0.04, df at t=1.0, then on to 1.2 (residual gmin*Kn*d)
    c, r = compare("mode I", [(0, (0, 0, 0)), (1.2, (0.012, 0, 0))], 0.01, exe, work)
    ck += c
    if r:
        ck.append(Check("mode I: W/W_exact", r[0] / w_exact(), 1.0, tol=5e-5))
        ck.append(Check("mode I: peak traction / Tn0", max(F[0] for _, F in r[1]) / AREA / TN0,
                        1.0, tol=1e-5))
        ck.append(Check("mode I: residual T at 1.2 df / (gmin Kn d)",
                        r[1][-1][1][0] / AREA / (GMIN * KN * 0.012), 1.0, tol=1e-4,
                        why="a broken facet still pulls with gmin*Kn*d in tension"))
    # mode II along s1: d_eff = sqrt(beta) ds, so ds = d_eff/sqrt(beta)
    c, r = compare("mode II", [(0, (0, 0, 0)), (1.2, (0, 0.012 * s, 0))], 0.01, exe, work)
    ck += c
    if r:
        ck.append(Check("mode II: W/W_exact", r[0] / w_exact(), 1.0, tol=5e-5))
        ck.append(Check("mode II: peak shear / Ts0", max(F[1] for _, F in r[1]) / AREA / TS0,
                        1.0, tol=1e-5))
    # 45 deg in the (n, s2) plane: dn = ds2 = d_eff/sqrt(1+beta)
    q = 0.012 / math.sqrt(1 + BETA)
    c, r = compare("mixed 45deg", [(0, (0, 0, 0)), (1.2, (q, 0, q))], 0.01, exe, work)
    ck += c
    if r:
        ck.append(Check("mixed 45deg: W/W_exact", r[0] / w_exact(), 1.0, tol=5e-5,
                        why="one Gc for every mode mix"))
    # unloading to zero along the secant, reloading on it, then softening on
    c, r = compare("unload-reload", [(0, (0, 0, 0)), (0.4, (0.004, 0, 0)), (0.8, (0, 0, 0)),
                                     (1.2, (0.008, 0, 0))], 0.01, exe, work, energy=False)
    ck += c
    if r:
        f = dict(r[1])
        ck.append(Check("unload: T at d=0 after damage", f[0.8][0], 0.0, tol=1e-6,
                        why="secant unloading returns to the origin"))
        ck.append(Check("unload: T(d=0.002)/((1-D(0.004)) Kn d)", f[0.6][0] / AREA /
                        ((1 - dback(0.004)) * KN * 0.002), 1.0, tol=1e-5))
    # compression after damage: full Kn
    c, r = compare("compression", [(0, (0, 0, 0)), (0.4, (0.004, 0, 0)), (0.6, (-0.001, 0, 0))],
                   0.01, exe, work, energy=False)
    ck += c
    if r:
        ck.append(Check("compression: T(-0.001)/(Kn d)", r[1][-1][1][0] / AREA / (KN * -0.001),
                        1.0, tol=1e-5, why="no damage in compression"))
    # shear under compression: no friction, so the shear law and Gc are untouched
    c, r = compare("shear under compression", [(0, (0, 0, 0)), (0.1, (-0.001, 0, 0)),
                                                (1.3, (-0.001, 0.012 * s, 0))], 0.01, exe, work,
                   energy=False)
    ck += c
    if r:
        tn = [F[0] / AREA for _, F in r[1][10:]]
        ck.append(Check("shear+compr.: Tn stays Kn*dn", max(abs(x / (KN * -0.001) - 1) for x in tn),
                        0.0, tol=1e-5))
        w, pj, pt = 0.0, 0.0, 0.0
        for (t, F), j in zip(r[1][10:], r[2][10:]):
            w += 0.5 * (F[1] / AREA + pt) * (j[1] - pj)
            pj, pt = j[1], F[1] / AREA
            if math.sqrt(BETA) * j[1] >= DF * (1 - 1e-12):
                break
        ck.append(Check("shear+compr.: W_shear/W_exact", w / w_exact(), 1.0, tol=5e-5,
                        why="compression does not raise the shear work: no friction in UC6"))
    # viscosity: exact backward-Euler recursion and convergence to the ODE
    mu = 0.05
    path = [(0, (0, 0, 0)), (1.2, (0.012, 0, 0))]
    errs = []
    for dt in (0.02, 0.005):
        c, r = compare("viscous dt=%g" % dt, path, dt, exe, work, mu=mu, energy=False)
        ck += c
        if r:
            # Duvaut-Lions ODE, integrated finely (RK4, 1e-4 steps), at the output times
            dv, t, k, e = 0.0, 0.0, 0, 0.0
            h = 1e-4
            fd = lambda tt, x: (dback(max(D0 * 0, 0.01 * tt)) - x) / mu
            for (tt, F) in r[1]:
                while t < tt - 1e-12:
                    hh = min(h, tt - t)
                    k1 = fd(t, dv); k2 = fd(t + hh / 2, dv + hh / 2 * k1)
                    k3 = fd(t + hh / 2, dv + hh / 2 * k2); k4 = fd(t + hh, dv + hh * k3)
                    dv += hh / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
                    t += hh
                dvs = [row[5] for row in r[3][tt]]
                e = max(e, abs(dvs[0] - min(dv, 1.0)))
            errs.append(e)
    if len(errs) == 2:
        ck.append(Check("viscous: Dv error vs ODE, dt 0.02 -> 0.005 (ratio)", errs[0] / errs[1],
                        4.0, rtol=0.35,
                        why="backward Euler is first order: a 4x smaller step, ~4x smaller error"))
        ck.append(Check("viscous: Dv error vs ODE at dt=0.005", errs[1], 0.0, tol=0.02))
    ck.append(Check("W_exact/Gc (the gmin sliver)", w_exact() / GC, None, cmp="info"))
    # rigid rotation of an open, elastic facet: the traction must turn with it
    ck += rotation(work, exe)
    return ck


# Switch of the co-rotated frame, once the developer names it; empty = the
# reference-frame law, where the rotation checks are XFAIL.
ROT_ENV = {}
ROT_FIXED = False


def rotation(work, exe):
    """Open a facet elastically by jn along its normal, then turn the whole
    facet (both faces) rigidly by th about z, an axis in its plane.  An
    objective law gives the same traction in the turned axes: |RF| = Kn jn A
    along the turned normal (cos th, sin th, 0), nothing across it, and no
    damage.  The reference-frame law reads the turned opening as
    (jn cos th, jn sin th) in the old axes: |RF| = Kn jn A sqrt(cos^2 +
    beta^2 sin^2), off the normal, and with beta > 1 an effective separation
    sqrt(beta) jn sin th that can pass d0 - an elastic facet damaged by a
    rotation."""
    ck = []
    for ts0, jn, tag in ((TS0, 2e-4, "beta<1"), (500.0, 3.5e-4, "beta>1")):
        b = (ts0 / TN0) ** 2
        for deg in (30, 60, 90):
            th = math.radians(deg)
            d = os.path.join(work, "rot_%s_%d" % (tag, deg))
            write(os.path.join(d, "u.inp"), uc6deck.deck(None, (KN, TN0, ts0, GC, GMIN, 0.0), 0.1,
                                                         rot=(th, (jn, 0, 0))))
            rc, _ = ccx(d, "u", exe, ROT_ENV)
            if rc:
                ck.append(Check("rotation %s %d deg: run rc" % (tag, deg), rc, 0, tol=0))
                continue
            fx = totals(os.path.join(d, "u.dat"), "NPLUS", 0)[-1][1]
            fy = totals(os.path.join(d, "u.dat"), "NPLUS", 1)[-1][1]
            fn = fx * math.cos(th) + fy * math.sin(th)          # along the turned normal
            ft = -fx * math.sin(th) + fy * math.cos(th)         # across it, in the plane of rotation
            ref = KN * jn * AREA
            if tag == "beta<1":
                ck.append(Check("rotation %d deg: RF along turned normal / (Kn d A)" % deg, fn / ref, 1.0, tol=1e-4,
                                expect_fail=not ROT_FIXED,
                                why="" if ROT_FIXED else "reference axes: the turned opening is read partly as shear"))
                # old axes: across = (beta - 1) sin cos - zero at 90 deg, where
                # the whole opening is read as shear along the turned normal
                ck.append(Check("rotation %d deg: RF across turned normal / (Kn d A)" % deg, ft / ref, 0.0, tol=1e-4,
                                expect_fail=(not ROT_FIXED) and deg != 90))
            else:
                st = blocks(os.path.join(d, "u.dat"), "stresses")
                dmax = max(row[7] for rows in st.values() for row in rows)
                old = jn * math.sqrt(math.cos(th) ** 2 + b * math.sin(th) ** 2)
                ck.append(Check("rotation %d deg, Ts0>Tn0: no damage (d_max <= d0)" % deg, dmax, D0, cmp="le",
                                expect_fail=(not ROT_FIXED) and old > D0,
                                why="old axes: d_eff = d sqrt(cos^2 + beta sin^2) = %.3g against d0 = %.3g"
                                % (old, D0)))
    return ck
