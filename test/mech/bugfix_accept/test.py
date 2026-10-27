"""ACCEPTANCE for three developer bug fixes, each a minimal reproducer.

(1) Nonlocal neighbour radius with a RELATIVE tolerance.  The B1 bar of
    nl_benchmark in C3D4 at ell = 0.5, h = 0.125: the averaging radius 2 ell
    is exactly 8 h, so element centroids sit ON it.  The same mesh at x1e-3,
    x1 and x1e3 (ell scaled with it) must report the same mean neighbour
    count - the geometry is the same in other units.  Found in scale_bar:
    82.7 against 81.8 at x1/25, and rc=201 against 0.
(2) Deleting the last element of the model.  One C3D4 under uniaxial stress
    (bulk_law/onelem.py, no second element) driven to D = 1: the run must end
    cleanly, rc = 0, instead of SIGSEGV in results() (reported on 71ed70d).
(3) The dashpot audit line [DAMAGE VISCOUS DAMPING] on a run with deletions:
    every number it prints must be finite (it once read freed memory).
Each check is XFAIL while its fix is pending (PENDING below); a fix turns it
XPASS, which fails the gate until PENDING is flipped - the gate notices.
Measured on 61d1c52: (1) 82.7250 at x1 and x1e3, 82.4625 at x1e-3 (spread
0.2625): XFAIL.  (2) rc = -11, SIGSEGV: XFAIL.  (3) 8 numbers in 3 audit
lines, all finite, on a run with 8 deletions: PASS already - the
use-after-free was fixed in 71ed70d; the audit's two accounting flaws are
pinned in test damper.
"""
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "nl_benchmark"))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "bulk_law"))
from mechlib import Check, ccx, cross_bar, write  # noqa: E402
import bar1d  # noqa: E402
import onelem  # noqa: E402

WHAT = "ACCEPTANCE bug fixes: neighbour count scale-free at 2ell = 8h; lone element deleted cleanly; finite dashpot audit"
SOURCE = "reproducers from scale_bar, bulk_law and the damper audit"

PENDING = {"radius": True, "lone": True, "audit": False}


def neighbours(work, k, exe):
    d = os.path.join(work, "nb_%g" % k)
    write(os.path.join(d, "t.inp"), bar1d.deck(80, length=10 * k, width=k, ell=0.5 * k, uf=0.1 * k,
                                               eltype="C3D4", uend=1e-4 * k, dt=1.0))
    rc, _ = ccx(d, "t", exe, {"CCX_DAMAGE_CHARLEN": "1", "CCX_DAMAGE_NONLOCAL_MODE": "INTEGRAL"})
    m = re.search(r"mean neighbours=\s*(\S+)", open(os.path.join(d, "run.log"), errors="replace").read())
    return rc, float(m.group(1)) if m else float("nan")


def run(work, exe):
    ck = []
    nb = {k: neighbours(work, k, exe) for k in (1e-3, 1.0, 1e3)}
    ck.append(Check("radius: runs (rc)", max(abs(r) for r, _ in nb.values()), 0, tol=0))
    spread = max(v for _, v in nb.values()) - min(v for _, v in nb.values())
    ck.append(Check("radius: mean neighbours at x1e-3 / x1 / x1e3 identical", spread, 0.0, tol=0.0,
                    expect_fail=PENDING["radius"],
                    why="counts %s" % " / ".join("%.4f" % nb[k][1] for k in (1e-3, 1.0, 1e3))))
    # (2) one C3D4, nothing else, driven to deletion
    mat = ["*MATERIAL, NAME=M", "*ELASTIC", "200000., 0.3", "*PLASTIC", "400., 0.",
           "*DAMAGE INITIATION, CRITERION=RICETRACEY, EVOLUTION=DISPLACEMENT", "0.005, 1.0, 0.01, 0."]
    d = os.path.join(work, "lone")
    write(os.path.join(d, "o.inp"), onelem.deck("C3D4", (1, None, None), 0.06, mat, 0.004, a=0.5, dummy=False))
    rc, _ = ccx(d, "o", exe)
    ck.append(Check("lone C3D4 deleted: rc = 0 (no SIGSEGV)", rc, 0, tol=0, expect_fail=PENDING["lone"],
                    why="rc=%d (-11 is SIGSEGV)" % rc))
    # (3) dashpot audit on a run with deletions
    r = cross_bar(os.path.join(work, "audit"), 24, None, "LOCAL", exe=exe,
                  env={"CCX_DAMAGE_VISCOUS_DAMPING": "2e-4"})
    lines = [l for l in open(os.path.join(work, "audit", "run.log"), errors="replace")
             if "[DAMAGE VISCOUS DAMPING]" in l]
    nums = [float(x) for l in lines for x in re.findall(r"[-+]?\d+\.\d+e[-+]\d+", l)]
    ok = bool(nums) and r["ndel"] > 0 and all(math.isfinite(x) and abs(x) < 1e30 for x in nums)
    ck.append(Check("audit: %d numbers in %d lines finite, %d deletions" % (len(nums), len(lines), r["ndel"]),
                    ok, True, cmp="true", expect_fail=PENDING["audit"]))
    return ck
