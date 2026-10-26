"""Shared helpers for the mechanics tests in test/mech/.

Each test is a directory test/mech/<name>/ holding test.py, which defines

    WHAT    one line: the mechanism and the reference it is checked against
    SOURCE  where the reference comes from (paper and equation, or derivation)
    def run(work, exe) -> list of Check

run() generates its decks under `work`, runs the solver and returns checks.
A Check carries the measured value, the reference and the tolerance, so the
table printed by run.py is the record; a test never prints "ok" on its own.
"""
import math
import os
import re
import subprocess
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class Check:
    """One assertion: |got - ref| <= tol (abs) or <= rtol*|ref| (rel).

    cmp="info" records a measured number with no pass/fail criterion yet
    (a reference still being established); it never fails the test.

    expect_fail marks a KNOWN LIMITATION of the code: the check is still
    evaluated and printed, and the test is counted as passing only while it
    keeps failing (XFAIL).  If it starts passing, that is reported (XPASS)
    so the limitation is re-examined instead of being silently forgotten.
    """

    def __init__(self, name, got, ref, tol=None, rtol=None, why="",
                 expect_fail=False, cmp="abs"):
        self.name, self.got, self.ref = name, got, ref
        self.tol, self.rtol, self.why = tol, rtol, why
        self.expect_fail = expect_fail
        self.cmp = cmp            # "abs", "le" (got<=ref), "ge" (got>=ref)

    @property
    def ok(self):
        g, r = self.got, self.ref
        if g is None or (isinstance(g, float) and math.isnan(g)):
            return False
        if self.cmp == "le":
            return g <= r
        if self.cmp == "ge":
            return g >= r
        if self.cmp == "true":
            return bool(g)
        if self.cmp == "info":
            return True
        lim = 0.0
        if self.tol is not None:
            lim = max(lim, self.tol)
        if self.rtol is not None:
            lim = max(lim, self.rtol * abs(r))
        return abs(g - r) <= lim

    @property
    def status(self):
        if self.cmp == "info":
            return "INFO"
        if self.expect_fail:
            return "XFAIL" if not self.ok else "XPASS"
        return "PASS" if self.ok else "FAIL"

    @property
    def counts_as_pass(self):
        return self.status in ("PASS", "XFAIL", "INFO")

    def line(self):
        def f(v):
            if isinstance(v, bool) or v is None:
                return str(v)
            if isinstance(v, (int, float)):
                return "%.6g" % v
            return str(v)
        if self.cmp == "abs":
            tol = []
            if self.tol is not None:
                tol.append("+-%.3g" % self.tol)
            if self.rtol is not None:
                tol.append("+-%.3g%%" % (100 * self.rtol))
            rel = ""
            if isinstance(self.ref, (int, float)) and self.ref and self.got is not None:
                rel = " (%+.2f%%)" % (100.0 * (self.got - self.ref) / abs(self.ref))
            crit = "ref %s %s%s" % (f(self.ref), " or ".join(tol), rel)
        elif self.cmp == "true":
            crit = "must hold"
        elif self.cmp == "info":
            crit = "(measured, no criterion)" if self.ref is None else "(compare %s)" % f(self.ref)
        else:
            crit = "%s %s" % ({"le": "<=", "ge": ">="}[self.cmp], f(self.ref))
        s = "  %-5s %-46s got %-12s %s" % (self.status, self.name, f(self.got), crit)
        if self.why:
            s += "\n        " + self.why
        return s


def ccx(rundir, job, exe, env=None, timeout=120):
    """Run the solver on <rundir>/<job>.inp with 1 thread; return (rc, seconds).

    Thread counts are fixed because the damage trajectory is not reproducible
    across them (test/regress/cases.json says why)."""
    e = dict(os.environ)
    for k in list(e):
        if k.startswith("CCX_"):
            del e[k]               # a stray switch in the caller's shell must not leak in
    e.update(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    e.update(env or {})
    # MECH_RED_ENV="NAME=VALUE ..." is applied LAST, over the test's own
    # switches: it is how a test is shown red on a wrong variant (a switch
    # turned off, a parameter changed) without editing the test.
    for kv in os.environ.get("MECH_RED_ENV", "").split():
        k, v = kv.split("=", 1)
        e[k] = v
    t0 = time.time()
    with open(os.path.join(rundir, "run.log"), "w") as log:
        try:
            rc = subprocess.call([exe, job], cwd=rundir, env=e, stdout=log,
                                 stderr=subprocess.STDOUT, timeout=timeout)
        except subprocess.TimeoutExpired:
            rc = -999
    return rc, time.time() - t0


def totals(datpath, nset=None, comp=0):
    """[(time, total force component)] from '*NODE PRINT ..., TOTALS=ONLY' blocks."""
    out, pend = [], None
    for line in open(datpath):
        m = re.search(r"total force \(fx,fy,fz\) for set (\S+) and time\s+([-\dEe.+]+)", line)
        if m:
            pend = float(m.group(2)) if (nset is None or m.group(1).upper() == nset.upper()) else None
            continue
        if pend is not None:
            f = line.split()
            if len(f) == 3:
                try:
                    out.append((pend, float(f[comp])))
                except ValueError:
                    pass
                pend = None
    return out


def blocks(datpath, key):
    """{time: [[floats], ...]} for every .dat block whose header contains key."""
    out, cur = {}, None
    for line in open(datpath):
        if key in line and "time" in line:
            t = float(line.split()[-1])
            cur = out.setdefault(t, [])
            continue
        s = line.split()
        if not s:
            continue
        try:
            vals = [float(v) for v in s]
        except ValueError:
            cur = None
            continue
        if cur is not None:
            cur.append(vals)
    return out


def trapz(xy):
    w = 0.0
    for (x0, y0), (x1, y1) in zip(xy, xy[1:]):
        w += 0.5 * (y0 + y1) * (x1 - x0)
    return w


def interp(xy, x):
    if x <= xy[0][0]:
        return xy[0][1]
    for (x0, y0), (x1, y1) in zip(xy, xy[1:]):
        if x1 >= x:
            return y0 if x1 == x0 else y0 + (y1 - y0) * (x - x0) / (x1 - x0)
    return xy[-1][1]


def deleted_count(rundir, job):
    n = 0
    p = os.path.join(rundir, "run.log")
    if os.path.exists(p):
        for l in open(p, errors="replace"):
            m = re.search(r"\[DAMAGE COMMIT\] batch=\d+ .* deleted=(\d+)", l)
            if m:
                n += int(m.group(1))
    return n


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(text)


# --------------------------------------------------------------------------
# The notched bar of test/crackband (mkcross.py), shared by the nonlocal
# tests.  6 x 1 x 1, C3D4, refined ALONG the axis at fixed length, uniform
# material, a 10 % linear-tent notch over the middle third (kinks at x=2,3,4,
# which are nodes for every nslice divisible by 6, so every mesh is the same
# solid - run_nlobjectivity.sh preflights that).  sigma_0 = 408, u_f = 0.02,
# G_f = sigma_0 u_f / 2 = 4.08 per unit area; notch area 0.81.
BAR_LEN, BAR_UEND, BAR_AREA, BAR_GF = 6.0, 0.30, 0.81, 4.08


def cross_bar(d, nslice, ell=None, mode=None, env=None, exe=None, gen=()):
    """Generate and run one arm; return a dict of what the tests compare.

    ell=None is the LOCAL model.  The arm always runs with the directional
    width (CHARLEN=1), the viscosity 1e-3 of run_nlobjectivity.sh, and the
    separation end of gate case cross-separation-end, so an arm that has
    actually severed the bar ends rc=0 instead of walling on the last
    ligament (measured: without it INTEGRAL at h=ell/4 stops rc=201 at 0.04 %
    of the peak force, after the work is done)."""
    os.makedirs(d, exist_ok=True)
    subprocess.check_call(
        ["python3", os.path.join(ROOT, "test/crackband/mkcross.py"),
         "--nslice", str(nslice), "--ltot", str(BAR_LEN), "--ncross", "1",
         "--alldamage", "--trigger", "1.0", "--notch", "0.10",
         "--notchlen", "0.33333333333333333", "--uend", str(BAR_UEND),
         "-o", os.path.join(d, "t.inp")] + list(gen), stdout=subprocess.DEVNULL)
    e = {"CCX_DAMAGE_CHARLEN": "1", "CCX_DAMAGE_VISCOSITY": "1.e-3",
         "CCX_FRACTURE_TERMINATION": "NLEFT:NRIGHT"}
    if ell is not None:
        p = os.path.join(d, "t.inp")
        s = open(p).read()
        n = s.count("EVOLUTION=DISPLACEMENT\n")
        assert n == 2, "mkcross.py card changed: %d damage cards" % n
        s = s.replace("EVOLUTION=DISPLACEMENT\n",
                      "EVOLUTION=DISPLACEMENT, NONLOCAL=%r\n" % ell)
        open(p, "w").write(s)
        e["CCX_DAMAGE_NONLOCAL_MODE"] = mode
    e.update(env or {})
    rc, sec = ccx(d, "t", exe, e, timeout=300)
    c = [(0.0, 0.0)] + [(BAR_UEND * t, f) for t, f in totals(os.path.join(d, "t.dat"))]
    ip = max(range(len(c)), key=lambda i: c[i][1])
    r = dict(rc=rc, sec=sec, curve=c, fmax=c[ip][1], upk=c[ip][0],
             wtot=trapz(c), wpost=trapz(c[ip:]), fend=c[-1][1],
             ndel=deleted_count(d, "t"))
    try:
        import sys
        sys.path.insert(0, os.path.join(ROOT, "test", "crackband"))
        import bandwidth
        r["width"] = bandwidth.width(d)[0]
    except Exception:
        r["width"] = float("nan")
    return r


def spread(vals):
    """(max - min) / |last|: the last value is the finest mesh / step."""
    return (max(vals) - min(vals)) / abs(vals[-1])
