"""Nodal viscous damper (CCX_DAMAGE_VISCOUS_DAMPING = c): force law, audit, c -> 0.

The dashpot adds -c m (u - u_ini)/dt to the static residual at every ACTIVE
degree of freedom, m = the node's quarter share of the reference C3D4
volume (unit density) - nonlingeo.c, calcresidual.c.

1. Force law and audit, exactly.  An elastic, un-notched bar (mkcross.py,
   nu = 0, uend = 0.006 so nothing yields) stretches as u_i = x_i U(t)/L, so
   over the step the dashpot dissipates
       E = c * U_end^2 / L^2 * sum_i m_i x_i^2 * sum_inc dt/T  = c U_end^2 S / L^2
   independent of the increments, the sum running over the dofs the force
   acts on.  The external work must exceed the undamped run by exactly E:
   W(c) - W(0) is a measurement of the force law that does not use the
   audit at all.  And the printed audit must report E.
2. Convergence to the static solution on a full fracture run (the notched
   bar of test/crackband, local model): W_total(c) -> W_total(0) as c -> 0,
   and the printed dissipated fraction falls in proportion to c.

The first increment is 1e-4 of the step: the dashpot force jumps from 0 at
t = 0, so the trapezoidal work under-counts half of increment 1's viscous
share; a first increment of 1e-4 keeps that at 1e-4.

Measured on 71ed70d.  Part 1, c = 0.2, h = 0.25: W(c)-W(0) is the active-dof
E to 0.01 %, so the FORCE is right.  The audit is not, twice, both XFAIL and
reported: it sums m du^2 over ALL dofs including the prescribed grip (x-dof
of NRIGHT), where no dashpot force acts, +6.7 % here; and it never audits
the LAST increment of the step (the audit of increment n runs at the top of
increment n+1), so it misses dt_last/T of a uniform-rate step: -5.00 % with
steps of 0.05, -0.415 % here.  The "as coded" check pins both to 1e-3.  Part 2: W - W(0) = +0.0663 / +0.0086 / +0.0020 / +0.0001 at
c = 0.2 / 0.02 / 0.002 / 0.0002 (W(0) = 15.2925); fraction 4.53e-3 / 4.53e-4
/ 4.53e-5 / 4.53e-6 - exactly linear.

RED PATH, run: MECH_RED_ENV=CCX_DAMAGE_VISCOUS_DAMPING=0.02 (one c forced
on every arm, the "undamped" ones included) turns four checks red: W(c) -
W(0) = 0 against 1.35e-5, the audit -90 %, both fraction ratios 1 against 0.1.

Tolerances.  Part 1: the bar is elastic and the field linear to the
solver's precision, so 1e-3 relative on the audit; the W difference is a
difference of two works of 0.18 that agree to 1e-4, so it gets 1 %.  Part 2: W within
0.1 % of the undamped run at c <= 2e-3 (measured 0.013 %); fraction linear
in c to 1 %.
"""
import os
import re
import subprocess

from mechlib import Check, ROOT, ccx, cross_bar, totals, trapz

WHAT = "viscous damper: force law and energy audit on an elastic bar (exact); W -> static as c -> 0"
SOURCE = "closed form for a linear displacement field; self-convergence c -> 0"

UEND, LBAR = 0.006, 6.0


def deck_masses(inp):
    co, el, sets, cur = {}, [], {}, None
    for line in open(inp):
        t = line.strip()
        if not t or t.startswith("**"):
            continue
        if t.startswith("*"):
            u, cur = t.upper(), None
            if u.startswith("*NODE") and "PRINT" not in u:
                cur = "N"
            elif u.startswith("*ELEMENT"):
                cur = "E"
            elif u.startswith("*NSET"):
                cur = re.search(r"NSET=(\w+)", u).group(1)
                sets[cur] = []
            continue
        f = [v for v in t.split(",") if v.strip()]
        if cur == "N":
            co[int(f[0])] = tuple(map(float, f[1:4]))
        elif cur == "E":
            el.append([int(v) for v in f[1:5]])
        elif cur:
            sets[cur] += [int(v) for v in f]
    m = {}
    for e in el:
        p = [co[n] for n in e]
        a, b, c = ([p[k][i] - p[0][i] for i in range(3)] for k in (1, 2, 3))
        v = abs(a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6
        for n in e:
            m[n] = m.get(n, 0.0) + v / 4
    return co, m, sets


def audit(rundir):
    tot = [l for l in open(os.path.join(rundir, "run.log")) if "step total: dissipated" in l]
    if not tot:
        return None, None
    f = re.findall(r"[-+]?\d\.\d+e[-+]\d+", tot[-1])
    return float(f[0]), float(f[1])


def last_dt(d):
    rows = [l.split() for l in open(os.path.join(d, "t.sta")) if l.strip()[:1].isdigit()]
    return float(rows[-1][6])


def elastic(work, c, exe, dt0=1e-4, dtmax=0.01):
    d = os.path.join(work, "elastic_c%g_dt%g" % (c, dt0))
    os.makedirs(d, exist_ok=True)
    subprocess.check_call(["python3", os.path.join(ROOT, "test/crackband/mkcross.py"),
                           "--nslice", "24", "--ltot", str(LBAR), "--ncross", "1", "--alldamage",
                           "--trigger", "1.0", "--notch", "0", "--uend", str(UEND),
                           "--dt0", repr(dt0), "--dtmax", repr(dtmax),
                           "-o", os.path.join(d, "t.inp")], stdout=subprocess.DEVNULL)
    env = {"CCX_DAMAGE_CHARLEN": "1"}
    if c > 0:
        env["CCX_DAMAGE_VISCOUS_DAMPING"] = repr(c)
    rc, _ = ccx(d, "t", exe, env)
    w = trapz([(0.0, 0.0)] + [(UEND * t, f) for t, f in totals(os.path.join(d, "t.dat"))])
    return d, rc, w


def run(work, exe):
    ck = []
    c = 0.2
    d, rc, wc = elastic(work, c, exe)
    _, rc0, w0 = elastic(work, 0.0, exe)
    ck.append(Check("elastic bar: runs (rc)", max(abs(rc), abs(rc0)), 0, tol=0))
    co, m, sets = deck_masses(os.path.join(d, "t.inp"))
    held = set(sets["NLEFT"]) | set(sets["NRIGHT"])
    s_act = sum(m[n] * co[n][0] ** 2 for n in m if n not in held)
    s_all = sum(m[n] * co[n][0] ** 2 for n in m)
    e_act = c * UEND ** 2 * s_act / LBAR ** 2
    e_all = c * UEND ** 2 * s_all / LBAR ** 2
    ck.append(Check("W(c) - W(0) = c U^2 sum_active(m x^2)/L^2", wc - w0, e_act, rtol=0.01,
                    why="the force law, measured through the external work alone"))
    ed, _ = audit(d)
    ck.append(Check("audit: printed dissipation = active-dof E", ed, e_act, rtol=1e-3,
                    expect_fail=True,
                    why="the audit sums m du^2 over ALL dofs, the grip's prescribed one too (reported)"))
    ck.append(Check("audit = all-dof sum without the last increment (as coded)", ed,
                    e_all * (1 - last_dt(d)), rtol=1e-3))
    # equal steps of 0.05: the last increment is 5 % of E
    d5, rc5, _ = elastic(work, c, exe, dt0=0.05, dtmax=0.05)
    ed5, _ = audit(d5)
    ck.append(Check("audit counts the last increment (dt = 0.05 steps)", ed5, e_all, rtol=1e-3,
                    expect_fail=True, why="the audit of an increment runs at the top of the next "
                    "one, which the last never gets: -5.00 % = dt_last/T (reported)"))
    # part 2: c -> 0 on a fracture run
    runs = {}
    for cc in (0.0, 0.2, 0.02, 0.002):
        runs[cc] = cross_bar(os.path.join(work, "frac_c%g" % cc), 24, None, "LOCAL", exe=exe,
                             env={"CCX_DAMAGE_VISCOUS_DAMPING": repr(cc)} if cc else None)
        runs[cc]["audit"] = audit(os.path.join(work, "frac_c%g" % cc))
    ck.append(Check("fracture runs (rc)", max(abs(r["rc"]) for r in runs.values()), 0, tol=0))
    w0 = runs[0.0]["wtot"]
    for cc in (0.2, 0.02, 0.002):
        ck.append(Check("c=%g: W_total - W(0)" % cc, runs[cc]["wtot"] - w0, None, cmp="info"))
    ck.append(Check("c=0.002: W_total within 0.1 % of static", runs[0.002]["wtot"], w0, rtol=1e-3))
    fr = [runs[cc]["audit"][0] / runs[cc]["audit"][1] for cc in (0.2, 0.02, 0.002)]
    ck.append(Check("dissipated fraction ratio c=0.02/0.2", fr[1] / fr[0], 0.1, rtol=0.01))
    ck.append(Check("dissipated fraction ratio c=0.002/0.02", fr[2] / fr[1], 0.1, rtol=0.01))
    return ck
