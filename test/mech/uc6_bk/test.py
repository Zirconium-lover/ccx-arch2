"""ACCEPTANCE: UC6 mixed mode by Benzeggagh-Kenane (Turon et al. 2006/2010).

Interface agreed with the developer (27.10): *USER SECTION with CONSTANTS=8,
the six old constants in place (Kn, Tn0, Ts0, G_Ic in slot 4, g_min, mu) and
7 = G_IIc, 8 = eta; a 6-constant section runs the old law bit for bit.
Law (library notes/cohesive-uc6.md sec. 4, cohesive_ref.py class BK):
    Ks = Kn (G_Ic/G_IIc)(Ts0/Tn0)^2   (Turon 2010), kappa = Ks/Kn
    lambda = sqrt(<dn>^2 + kappa s^2),  B = kappa s^2 / lambda^2
    lambda0(B) = (Tn0/Kn) sqrt(1 + (G_IIc/G_Ic - 1) B^eta)
    m = 2 G_Ic Kn / Tn0^2,  h = max lambda/lambda0(B),  D = m(h-1)/(h(m-1))
References for proportional loading at fixed B (closed form):
    work to failure  W = G_c(B) = G_Ic + (G_IIc - G_Ic) B^eta   (BK 1996)
    peak traction    T_n = Tn0 (mode I), T_s = Ts0 (mode II)
    reduction        G_IIc = G_Ic gives Ks = beta Kn, the 6-constant law,
                     to round-off (~1e-13 by the developer's analysis; here
                     the 7 printed digits)
Parameters of the library's check: Kn = 1e5, Tn0 = 50, Ts0 = 80, G_Ic = 0.3,
G_IIc = 1.0, eta = 2, g_min = 1e-6 -> G_c(B) = 0.300 / 0.344 / 0.475 /
0.694 / 1.000 at B = 0 / 0.25 / 0.5 / 0.75 / 1.

The path is lambda = t along a fixed direction, delta = t (cos th, sin th /
sqrt(kappa), 0), B = sin^2 th, to 1.2 lambda_f(B) in 1000 increments.  The
kinks at lambda0 and lambda_f are not on the grid: trapezoidal error <= 1e-4
of W; g_min = 1e-6 adds < 3e-5.  Tolerance 1e-3.

Acceptance test: red until the developer's B-K code is in.  Measured on
e2d474f (no B-K): the 8-constant section is read as the old law with Gc =
G_Ic, so W = 0.300 at every B - 0.873 / 0.632 / 0.433 / 0.300 of G_c(B) at
B = 0.25 / 0.5 / 0.75 / 1 - and the reduction checks pass (that is what they
must keep doing).  On ab4720b: 15/15, W/G_c(B) = 0.999993 at every B (the
path is scaled to lambda_f(B), so the trapezoidal error is the same at every
mix), peaks 0.99965, reduction identical to the printed digit.
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "uc6_law"))
from mechlib import Check, ccx, totals, write  # noqa: E402
import uc6deck  # noqa: E402

WHAT = "ACCEPTANCE UC6 B-K mixed mode: W = G_c(B) at B 0..1, peaks Tn0/Ts0, reduction to the old law at G_IIc = G_Ic"
SOURCE = "Benzeggagh & Kenane 1996; Turon et al. 2006 eqs. 33-53, 2010 (Ks); library cohesive_ref.py BK"

KN, TN0, TS0, GIC, GIIC, ETA, GMIN = 1e5, 50.0, 80.0, 0.3, 1.0, 2.0, 1e-6
AREA, N = 0.5, 1000


def gc(b, giic=GIIC):
    return GIC + (giic - GIC) * b ** ETA


def facet(work, name, path, props, exe):
    d = os.path.join(work, name)
    write(os.path.join(d, "u.inp"), uc6deck.deck(path, props, path[-1][0] / N))
    rc, _ = ccx(d, "u", exe)
    dat = os.path.join(d, "u.dat")
    f = [totals(dat, "NPLUS", c) for c in range(3)] if rc == 0 else [[], [], []]
    return rc, [(t, (fx / AREA, fy / AREA, fz / AREA)) for (t, fx), (_, fy), (_, fz) in zip(*f)]


def run(work, exe):
    ck = []
    kappa = (GIC / GIIC) * (TS0 / TN0) ** 2
    l0i = TN0 / KN
    m = 2 * GIC * KN / TN0 ** 2
    props = (KN, TN0, TS0, GIC, GMIN, 0.0, GIIC, ETA)
    for b in (0.0, 0.25, 0.5, 0.75, 1.0):
        th = math.asin(math.sqrt(b))
        lf = m * l0i * math.sqrt(1 + (GIIC / GIC - 1) * b ** ETA)
        u = (math.cos(th), math.sin(th) / math.sqrt(kappa), 0.0)
        tend = 1.2 * lf
        path = [(0.0, (0, 0, 0)), (tend, tuple(tend * x for x in u))]
        rc, tr = facet(work, "B%g" % b, path, props, exe)
        ck.append(Check("B=%g: run rc" % b, rc, 0, tol=0))
        if rc or not tr:
            continue
        w, pt, pl = 0.0, (0.0, 0.0, 0.0), 0.0
        for t, T in tr:                      # lambda = t on this path
            w += sum(0.5 * (T[c] + pt[c]) * (t - pl) * u[c] for c in range(3))
            pt, pl = T, t
            if t >= lf * (1 - 1e-12):
                break
        ck.append(Check("B=%g: W / G_c(B) = %.4f" % (b, gc(b)), w / gc(b), 1.0, tol=1e-3))
        if b == 0.0:
            ck.append(Check("mode I: peak T_n / Tn0", max(T[0] for _, T in tr) / TN0, 1.0, tol=5e-3))
        if b == 1.0:
            ck.append(Check("mode II: peak T_s / Ts0", max(T[1] for _, T in tr) / TS0, 1.0, tol=5e-3,
                            why="the peak falls between increments, lambda0 ~35 increments in"))
    # reduction: G_IIc = G_Ic with 8 constants against the 6-constant law, B = 0.5 path
    th = math.asin(math.sqrt(0.5))
    beta = (TS0 / TN0) ** 2
    u = (math.cos(th), math.sin(th) / math.sqrt(beta), 0.0)
    tend = 1.2 * m * l0i
    path = [(0.0, (0, 0, 0)), (tend, tuple(tend * x for x in u))]
    r6, t6 = facet(work, "red_6const", path, (KN, TN0, TS0, GIC, GMIN, 0.0), exe)
    r8, t8 = facet(work, "red_8const", path, (KN, TN0, TS0, GIC, GMIN, 0.0, GIC, ETA), exe)
    ck.append(Check("reduction G_IIc = G_Ic: runs (rc)", max(abs(r6), abs(r8)), 0, tol=0))
    if not (r6 or r8) and t6 and t8:
        dmax = max(abs(a[1][c] - b_[1][c]) for a, b_ in zip(t6, t8) for c in range(3)) / TN0
        ck.append(Check("reduction G_IIc = G_Ic: max |T_8 - T_6| / Tn0", dmax, 0.0, tol=1e-6,
                        why="Ks = beta Kn: the old law, to the printed digits"))
        ck.append(Check("reduction: same number of increments", len(t6) == len(t8), True, cmp="true"))
    return ck
