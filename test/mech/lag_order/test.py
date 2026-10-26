"""Accuracy of the LAGGED (staggered) nonlocal scheme against the implicit one.

The nonlocal driver is built from the last converged increment and read,
rescaled by dt_cur/dt_ref, on every Newton iteration of the next one.
FROZEN runs exactly that scheme with no averaging, so FROZEN against LOCAL
(implicit, same deck, same steps) isolates the error of the lag alone.  The
question it answers is whether FROZEN may stand in for LOCAL in a series:
it may if its distance from LOCAL vanishes as dt -> 0, at a known order, and
is small against the time-discretisation error LOCAL already carries (the
viscous regularisation, beta = dt/(eta+dt), is itself first order in dt).

Notched bar of test/crackband at h = 0.25, dtmax = 0.01 / 2^k, dt0 = dtmax/2.

Measured on 71ed70d, W_total:
  dtmax     LOCAL     FROZEN    FROZEN-LOCAL   GRADIENT(ell=0.5)
  0.01      15.2925   15.8666*  +0.574         13.9918
  0.005     15.8483   16.0706   +0.2223        14.0248
  0.0025    15.9902   16.0980   +0.1078        14.0731
  0.00125   16.0808   16.1134   +0.0326        14.0940
  0.000625  16.1286   16.1445   +0.0159        14.1118   (not run here: 11 s per arm)
  * rc=201 at u=0.077 with 2.6 % of the peak force left: FROZEN does not
    finish at the coarsest step where LOCAL does.  Kept as XFAIL.
So the lag error halves per halving (observed order 1.0-1.7, 1.4 over the
three middle steps) and is below LOCAL's own distance from its dt->0
limit at the same step (0.60 at dtmax 0.0025 against the Richardson limit
16.171): the lag is not the dominant time error.

CCX_DAMAGE_NONLOCAL_INNER=1 (refresh inside Newton) brings FROZEN to within
0.12 % of LOCAL at dtmax 0.0025 (15.9902 vs 16.0087), 0.5 % allowed.

RED PATH, run: a scratch build whose lagged increment is multiplied by 1.1
after the dt rescaling (damnonlocal.f, the line after val=val*dtcur/dtref) -
a lag that is no longer consistent - gives FROZEN-LOCAL of -8.5/-9.0/-9.5 %,
order -0.09, ratio 7.9, INNER 9.5 % off: three checks red.  The opposite
switch, CCX_DAMAGE_NONLOCAL_DTSCALE=0, does NOT turn it red - with a fixed
dtmax there are few cutbacks for the rescaling to matter on - and it is the
arm on which the coarse FROZEN run FINISHES (XPASS), see below.

Tolerances.  Order >= 0.8: first order is what a one-step lag can give and
the measured 1.4 has margin either side; a scheme that did not converge
(order ~0) or was not consistent fails it.  GRADIENT self-convergence 2 %:
the measured spread between dtmax 0.01 and 0.00125 is 0.7 %.
"""
import math
import os
from mechlib import Check, cross_bar

WHAT = "lagged nonlocal scheme: FROZEN -> implicit LOCAL as dt->0 (order), GRADIENT converges in dt"
SOURCE = "self-convergence against the implicit local solution of the same deck (no closed form)"

DTS = (0.005, 0.0025, 0.00125)


def arm(work, mode, dtm, exe):
    ell = {"LOCAL": None, "FROZEN": 0.01, "GRADIENT": 0.5}[mode]
    return cross_bar(os.path.join(work, "%s_%g" % (mode, dtm)), 24, ell, mode, exe=exe,
                     gen=("--dt0", repr(dtm / 2), "--dtmax", repr(dtm)))


def run(work, exe):
    ck = []
    loc = [arm(work, "LOCAL", d, exe) for d in DTS]
    frz = [arm(work, "FROZEN", d, exe) for d in DTS]
    diff = [f["wtot"] - l["wtot"] for f, l in zip(frz, loc)]
    for d, f, l, x in zip(DTS, frz, loc, diff):
        ck.append(Check("dtmax=%g FROZEN-LOCAL W rel." % d, x / l["wtot"], None, cmp="info"))
    ck.append(Check("all arms finish (rc)", max(abs(r["rc"]) for r in loc + frz), 0, tol=0))
    p = math.log(abs(diff[0]) / max(abs(diff[-1]), 1e-30), 2) / (len(DTS) - 1)
    ck.append(Check("observed order of the lag error", p, 0.8, cmp="ge",
                    why="FROZEN-LOCAL in W_total over dtmax 0.005 -> 0.00125"))
    # how large the lag error is against LOCAL's own time error at the same
    # step; LOCAL's dt->0 limit by Richardson with order 1 from the two finest
    linf = 2 * loc[-1]["wtot"] - loc[-2]["wtot"]
    ck.append(Check("lag error / LOCAL time error at dtmax=0.0025",
                    abs(diff[1]) / abs(linf - loc[1]["wtot"]), 1.0, cmp="le",
                    why="measured 0.60: the lag adds less than the step already costs LOCAL"))
    inner = cross_bar(os.path.join(work, "FROZEN_INNER"), 24, 0.01, "FROZEN",
                      env={"CCX_DAMAGE_NONLOCAL_INNER": "1"}, exe=exe,
                      gen=("--dt0", "0.00125", "--dtmax", "0.0025"))
    ck.append(Check("FROZEN+INNER=1 equals LOCAL at dtmax=0.0025", inner["wtot"], loc[1]["wtot"],
                    rtol=0.005, why="refresh once per Newton iteration removes the lag: "
                    "measured +0.12 % against +0.67 % lagged"))
    grd = [arm(work, "GRADIENT", d, exe) for d in (0.01,) + DTS]
    ck.append(Check("GRADIENT W spread dtmax 0.01..0.00125",
                    (max(r["wtot"] for r in grd) - min(r["wtot"] for r in grd)) / grd[-1]["wtot"],
                    0.02, cmp="le"))
    coarse = arm(work, "FROZEN", 0.01, exe)
    ck.append(Check("FROZEN finishes at dtmax=0.01 (LOCAL does)", coarse["rc"], 0, tol=0,
                    expect_fail=True,
                    why="measured rc=201 at u=0.077, 2.6 % of peak force left; reported to the developer"))
    return ck
