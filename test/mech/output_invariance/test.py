"""The solution must not depend on how often results are written.

Writing output calls results() at the end of an increment.  If that call
re-integrates the damage law and leaves its answer in the committed state,
the next increment starts from a state that exists only because output was
asked for, and two runs of one deck that differ only in FREQUENCY follow
different trajectories.  Found by the developer on 71ed70d: the output
results() (nonlingeo.c ~16827, iout=2) rewrites the committed dam, in FROZEN
with a fresh dpsave; the notched bar's LOCAL arm leaves the frequency-1 run
at increment 23, the FROZEN arm ends rc=201 at frequency 1 and rc=0 at 1000.

Check: each deck is run with every output card at FREQUENCY=1 (as generated)
and at FREQUENCY=1000; the .sta (increment, attempt, iterations, step time)
and the .damage deletion history must be BYTE-identical, and the exit codes
equal.  Byte identity is the right bar: the two runs do the same arithmetic
on the same inputs, and nothing in a correct solver reads the output
frequency.  Arms: the notched bar of test/crackband (h = 0.25) LOCAL, LOCAL
without viscosity, FROZEN, GRADIENT ell=0.5; and the UC6 bar of
test/pathfollow/mkcohesive.py (short enough not to snap back), inviscid and
with mu > 0.

This is the developer's acceptance test: it is expected to FAIL until the
output pass stops writing the committed state.

Measured on e656a02 (71ed70d + e1be59f): every bulk arm differs - LOCAL from
.sta line 27 (.damage line 6), LOCAL without viscosity from line 26, FROZEN
from line 20 and rc 201 at frequency 1 against rc 0 at 1000, GRADIENT from
line 18.  Both UC6 arms are byte-identical with the facet driven to D = 1,
which is the green path: the comparison can pass, and the defect is in the
bulk damage state only.
"""
import filecmp
import os
import re
import subprocess

from mechlib import Check, ROOT, ccx

WHAT = "output invariance: .sta and .damage byte-identical at output FREQUENCY 1 and 1000"
SOURCE = "a correct solver does not read the output frequency (no reference value needed)"

OUTCARDS = re.compile(r"^\*(NODE PRINT|EL PRINT|NODE FILE|EL FILE|NODE OUTPUT|ELEMENT OUTPUT|CONTACT PRINT|CONTACT FILE)\b",
                      re.I)


def sparse(text, freq):
    out = []
    for line in text.splitlines():
        if OUTCARDS.match(line.strip()) and "FREQUENCY" not in line.upper():
            line = line.rstrip() + ", FREQUENCY=%d" % freq
        out.append(line)
    return "\n".join(out) + "\n"


def bar(path, ell=None):
    subprocess.check_call(["python3", os.path.join(ROOT, "test/crackband/mkcross.py"),
                           "--nslice", "24", "--ltot", "6.0", "--ncross", "1", "--alldamage",
                           "--trigger", "1.0", "--notch", "0.10", "--notchlen",
                           "0.33333333333333333", "--uend", "0.30", "-o", path],
                          stdout=subprocess.DEVNULL)
    if ell is not None:
        s = open(path).read().replace("EVOLUTION=DISPLACEMENT\n",
                                      "EVOLUTION=DISPLACEMENT, NONLOCAL=%r\n" % ell)
        open(path, "w").write(s)


def cohesive(path, mu):
    # L = 4: r = (Tn0/(df-d0))/(E/L) = 0.83 < 1, so displacement control
    # follows the whole softening branch and the facet reaches gmin
    subprocess.check_call(["python3", os.path.join(ROOT, "test/pathfollow/mkcohesive.py"),
                           "--nslice", "4", "--isplit", "2", "--uend", "0.02",
                           "--mu", repr(mu), "-o", path], stdout=subprocess.DEVNULL)


BAR_ENV = {"CCX_DAMAGE_CHARLEN": "1", "CCX_DAMAGE_VISCOSITY": "1.e-3",
           "CCX_FRACTURE_TERMINATION": "NLEFT:NRIGHT"}
ARMS = (
    ("bar LOCAL", lambda p: bar(p), BAR_ENV),
    ("bar LOCAL no viscosity", lambda p: bar(p), {k: v for k, v in BAR_ENV.items()
                                                  if k != "CCX_DAMAGE_VISCOSITY"}),
    ("bar FROZEN", lambda p: bar(p, 0.01), dict(BAR_ENV, CCX_DAMAGE_NONLOCAL_MODE="FROZEN")),
    ("bar GRADIENT l=0.5", lambda p: bar(p, 0.5), dict(BAR_ENV, CCX_DAMAGE_NONLOCAL_MODE="GRADIENT")),
    ("UC6 bar", lambda p: cohesive(p, 0.0), {}),
    ("UC6 bar mu=1e-3", lambda p: cohesive(p, 1e-3), {}),
)


def run(work, exe):
    ck = []
    for name, gen, env in ARMS:
        tag = re.sub(r"[^A-Za-z0-9]+", "_", name)
        res = {}
        for freq in (1, 1000):
            d = os.path.join(work, "%s_f%d" % (tag, freq))
            os.makedirs(d, exist_ok=True)
            p = os.path.join(d, "t.inp")
            gen(p)
            if freq != 1:
                text = sparse(open(p).read(), freq)
                open(p, "w").write(text)
            rc, _ = ccx(d, "t", exe, env)
            res[freq] = (d, rc)
        (d1, rc1), (d2, rc2) = res[1], res[1000]
        n1 = sum(1 for l in open(os.path.join(d1, "t.sta")) if l.strip()[:1].isdigit())
        ck.append(Check("%s: rc equal (%d vs %d), %d inc" % (name, rc1, rc2, n1), rc1 == rc2, True,
                        cmp="true"))
        for ext in ("sta", "damage"):
            a, b = os.path.join(d1, "t." + ext), os.path.join(d2, "t." + ext)
            if not os.path.exists(a) and not os.path.exists(b):
                continue
            same = os.path.exists(a) and os.path.exists(b) and filecmp.cmp(a, b, shallow=False)
            first = None
            if not same and os.path.exists(a) and os.path.exists(b):
                la, lb = open(a).read().splitlines(), open(b).read().splitlines()
                first = next((i + 1 for i, (x, y) in enumerate(zip(la, lb)) if x != y),
                             min(len(la), len(lb)) + 1)
            ck.append(Check("%s: t.%s identical" % (name, ext), same, True, cmp="true",
                            why="" if same else "first differing line %s" % first))
    return ck
