"""Fracture work of the C3D4 bar against the 1-D reference: the neck, not G_f.

The bar of nl_benchmark (library B1: 10 x 1 x 1, one weak slice, H = 0,
eps_f -> 0, u_f = 0.1, G_f = 5) cut into Kuhn tetrahedra - the family every
Fn run used.  The 1-D reference W/G_f = 1.008 assumes the plastic flow stays
in one layer.  In 3-D it cannot: an isochoric stretch needs lateral
contraction the elastic neighbours forbid, so the neck is ~0.6 of the
section wide (test/crackband README, "the floor is the specimen") and
    LOCAL crack band:  every element layer of the neck dissipates its own G_f,
                       W ~ G_f * w_neck / h  - grows as the mesh is refined;
    NONLOCAL:          the neck's plastic work is bounded by l_d and converges
                       in h, to a multiple of G_f set by w_neck / l_d.
Measured on e04c124 (viscosity 1e-3 and the separation end, so every arm
ruptures; the displacement per increment is the same 7.5e-4 on every arm,
because the lag error of the nonlocal scheme scales with it):
    LOCAL CHARLEN=0 (6V)^(1/3):  h=0.25 1.941, h=0.125 2.111
    LOCAL CHARLEN=1 projection:  h=0.25 5.509, h=0.125 8.464
    INTEGRAL ell=0.5:            h=ell/2 1.586 (rc 201, 7 % of peak left),
                                 ell/4 1.389, ell/8 1.335
THE LEVER, dimensionless.  With w_D = sum(D V)/A the damage-weighted band
width (bandwidth.py; deleted elements at D = 1) and L the length in the law
(h for the local projection, (6V)^(1/3) for CHARLEN=0, l_d = sqrt(pi) ell for
INTEGRAL), every arm obeys
    W / G_f = c * w_D / L,   c = 1.37 ... 1.57
    arm                 W/G_f   w_D     w_D/L   c
    LOCAL CL=0 h=0.25   1.941   0.829   1.315   1.48
    LOCAL CL=0 h=0.125  2.111   0.693   1.385   1.52
    LOCAL CL=1 h=0.25   5.509   0.908   3.633   1.52
    LOCAL CL=1 h=0.125  8.464   0.672   5.377   1.57
    INTEGRAL ell/4      1.389   0.898   1.013   1.37
    INTEGRAL ell/8      1.335   0.838   0.945   1.41
So the fracture work is one G_f per length L of damaged band, times c ~ 1.5
for the plastic work of the neck around it.  w_D is set by the SPECIMEN (0.7
to 0.9 of the unit section here): for the local band the lever is
section / h, which grows without bound as the mesh is refined; for the
nonlocal model it is section / l_d, fixed by the material length - and that
is the ratio a change of specimen scale moves.  c is checked to 15 % around
1.5; w_D/L is recorded.

The growth of the local arms and the nonlocal multiplier are KNOWN PROPERTIES
of the model in 3-D: pinned here (growth, convergence) and XFAIL against the
1-D reference.  Numbers measured without the developer's 390c3bf (output
pass read-only); rebaselined with the scheme the team chooses.

Tolerances: the INTEGRAL multiplier at most 1.6 (measured 1.335); growth of the local arms >= 5 % per halving of h (measured 8.8
and 54 %); INTEGRAL ell/4 against ell/8 within 6 % (measured 4.0 %).

RED PATH, run: MECH_RED_ENV=CCX_DAMAGE_NLWIDTH=0 (the element width in the
nonlocal law instead of l_d) gives the INTEGRAL arm W/G_f 5.79 at ell/8 and
the ceiling check goes red; the h-convergence check alone does NOT (5.55
against 5.79 at ell/4 and ell/8, 4.2 % apart) - which is why the ceiling is
there.
"""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "nl_benchmark"))
from mechlib import Check, ROOT, ccx, totals, trapz, write  # noqa: E402
import bar1d  # noqa: E402

WHAT = "C3D4 bar vs 1-D B1: local crack band W grows with 1/h, nonlocal converges to ~1.35 G_f"
SOURCE = "library notes/nonlocal-benchmarks.md sec. 1 (1-D reference 1.008); measured anatomy in nl_benchmark"

GF, REF = 5.0, 1.008


def arm(work, mode, charlen, nx, ell, uend, exe):
    d = os.path.join(work, "%s_cl%s_%d" % (mode, charlen, nx))
    write(os.path.join(d, "t.inp"), bar1d.deck(nx, ell=ell, eltype="C3D4", uend=uend, dt=0.00075 / uend))
    env = {"CCX_DAMAGE_CHARLEN": charlen, "CCX_FRACTURE_TERMINATION": "NLEFT:NRIGHT",
           "CCX_DAMAGE_VISCOSITY": "1.e-3"}
    if mode != "LOCAL":
        env["CCX_DAMAGE_NONLOCAL_MODE"] = mode
    rc, _ = ccx(d, "t", exe, env)
    c = [(0.0, 0.0)] + [(uend * t, f) for t, f in totals(os.path.join(d, "t.dat"))]
    sys.path.insert(0, os.path.join(ROOT, "test", "crackband"))
    import bandwidth
    try:
        wd = bandwidth.width(d, area=1.0)[0]       # sum(D V)/A, deleted elements at D = 1
    except Exception:
        wd = float("nan")
    h = 10.0 / nx
    lw = {("LOCAL", "1"): h, ("LOCAL", "0"): h ** (1.0 / 3)}.get((mode, charlen), math.sqrt(math.pi) * (ell or 0))
    return rc, trapz(c) / GF, c[-1][1], wd, lw


def run(work, exe):
    ck = []
    lever = {}
    for cl, ue in (("0", 0.3), ("1", 0.8)):
        a1, a2 = (arm(work, "LOCAL", cl, nx, None, ue, exe) for nx in (40, 80))
        (rc1, w1), (rc2, w2) = a1[:2], a2[:2]
        tag = "LOCAL CHARLEN=%s" % cl
        lever.update({tag + " h=0.25": a1, tag + " h=0.125": a2})
        ck.append(Check("%s: runs rupture (rc)" % tag, max(abs(rc1), abs(rc2)), 0, tol=0))
        ck.append(Check("%s: W/G_f grows when h halves (h=0.25 -> 0.125)" % tag, w2 / w1, 1.05, cmp="ge",
                        why="W/G_f %.3f -> %.3f: each layer of the neck dissipates its own G_f" % (w1, w2)))
        ck.append(Check("%s h=0.125: W/G_f = 1-D reference" % tag, w2, REF, rtol=0.1, expect_fail=True,
                        why="the neck is ~0.6 of the section, not one element layer"))
    r = {nx: arm(work, "INTEGRAL", "1", nx, 0.5, 0.3, exe) for nx in (40, 80, 160)}
    ck.append(Check("INTEGRAL ell=0.5 h=ell/2: W/G_f (coarse, may stop early)", r[40][1], None, cmp="info"))
    ck.append(Check("INTEGRAL ell/4, ell/8: runs rupture (rc)", max(abs(r[80][0]), abs(r[160][0])), 0, tol=0))
    ck.append(Check("INTEGRAL: W/G_f at ell/4 vs ell/8 (converges in h)", r[80][1], r[160][1], rtol=0.06))
    ck.append(Check("INTEGRAL h=ell/8: W/G_f multiplier ceiling (measured 1.335)", r[160][1], 1.6, cmp="le",
                    why="pins the nonlocal multiplier so a growth of it is seen"))
    lever.update({"INTEGRAL h=ell/4": r[80], "INTEGRAL h=ell/8": r[160]})
    for tag, a in lever.items():
        ck.append(Check("%s: w_D/L (the lever)" % tag, a[3] / a[4], None, cmp="info"))
        ck.append(Check("%s: (W/G_f) / (w_D/L)" % tag, a[1] / (a[3] / a[4]), 1.5, rtol=0.15,
                        why="one G_f per length L of damaged band, times the neck's plastic work"))
    ck.append(Check("INTEGRAL h=ell/8: W/G_f = 1-D reference", r[160][1], REF, rtol=0.1, expect_fail=True,
                    why="neck plastic work bounded by l_d, not zero: converges to ~1.35 G_f"))
    return ck
