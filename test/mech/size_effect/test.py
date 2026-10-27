"""Size effect of a notched cohesive plate: strength limit, LEFM limit, Bazant's law.

M(T) plate (plate.py): half-width W, centre crack 2a with a/W = 0.3, height
2H = 3W, one C3D8 layer (plane stress), the ligament a < x < W a UC6 plane;
elastic bulk E = 1e5, nu = 0.3; Kn = 1e8, Tn0 = 100, Gc = 1, gmin = 1e-6, so
l_ch = E Gc / Tn0^2 = 10.  Four sizes W = 1, 10, 30, 100: D/l_ch = 0.1 ... 10.
Every size is the same mesh (40 x 2*16), so only the physics changes.

References.
    small D:  sigma_N -> Tn0 (1 - a/W)           net-section strength
    large D:  sigma_N -> K_IC / (sqrt(pi a) F),  K_IC = sqrt(E Gc) (plane
              stress), F = sqrt(sec(pi a / 2W))  (Feddersen; Tada, Paris and
              Irwin, The Stress Analysis of Cracks Handbook, M(T) specimen)
    between:  Bazant sigma_N = B Tn0 / sqrt(1 + D/d0)  (Bazant 1984; fitted
              as 1/sigma^2 linear in D - INFO, d0 compared with l_ch)

Measured on d324ea1 (sigma_N/Tn0, peak of the converged increments):
    W/l_ch  0.1      1        3        10
            0.6849   0.5807   0.4468   0.2842
  0.1: 0.978 of the net-section strength 0.70.
  10:  0.924 of LEFM 0.3075.  Below, not above: the cohesive zone acts as a
       crack extension c_f, sigma_N = K_IC/sqrt(pi (a + c_f)) roughly, and
       (0.3075/0.2842)^2 - 1 = 0.17 of a gives c_f = 0.52 l_ch - the order of
       the process zone of a linear-softening law (~0.4-0.6 l_ch).  So the
       large-size end approaches LEFM from below, as it must for a notched
       body, and the 7.6 % is size effect not yet exhausted, not error.
  Mesh: nx = 80 gives 0.2841 at W/l_ch = 10 and 0.4467 at 3 against 0.2842 and
  0.4468 - converged.  Load: displacement control; the larger plates snap
  back at the peak and the run stops rc=201 after cutting its step to ~2e-6
  of the step at the limit point - so the last converged force IS the peak,
  which the "peak captured" check verifies from the .sta.

Earlier attempt, kept as a lesson: gmin = 1e-3 made gmin*Kn = 1e5 = E, and
the "broken" ligament carried load like the bulk - sigma_N 1.44 Tn0, no
fracture at all.  That is the library's group gmin Kn delta/Tn0 (5 there,
0.02 at delta_f now), and it is checked below (<= 0.05).

RED PATH, run: the same deck with gmin = 1e-3 (plate.deck(..., gmin=1e-3)),
measured at W/l_ch = 1: sigma_N 1.44 Tn0 - the ligament never breaks, and
the strength, monotonicity and LEFM checks cannot hold.

Tolerances: strength limit 5 %; the LEFM end between 0.85 and 1.0 of LEFM
(from below, see above); monotonic decrease; the Bazant fit within 5 % at
every point.  LONG: four runs, ~90 s - not in the gate; run.py --long.
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mechlib import Check, ccx, totals, write  # noqa: E402
import plate  # noqa: E402

WHAT = "size effect, notched UC6 plate D/l_ch 0.1..10: strength and LEFM limits, Bazant fit (LONG)"
SOURCE = "net-section strength; LEFM M(T) (Tada, Feddersen secant); Bazant 1984 size-effect law"
LONG = True

E, TN0, GC, AW, LCH = 1.0e5, 100.0, 1.0, 0.3, 10.0


def lefm(w):
    return math.sqrt(E * GC) / (math.sqrt(math.pi * AW * w) * math.sqrt(1 / math.cos(math.pi * AW / 2)))


def run(work, exe):
    ck = []
    sig = {}
    for w in (1.0, 10.0, 30.0, 100.0):
        d = os.path.join(work, "W%g" % w)
        sp = min(TN0 * (1 - AW), 1.2 * lefm(w))
        ue = 1.4 * sp / E * 3 * w * 1.15          # ~1.4x the predicted peak over 2H = 3W
        write(os.path.join(d, "p.inp"), plate.deck(w, a_over_w=AW, nx=40, ny=16, uend=ue, dt=1 / 50.0))
        rc, _ = ccx(d, "p", exe, timeout=600)
        c = totals(os.path.join(d, "p.dat"), "NTOP", 1)
        fm = max(f for _, f in c) if c else float("nan")
        ipk = [f for _, f in c].index(fm) if c else 0
        # the peak is captured if the force falls after it, or if the run cut
        # its step down to the limit point (snap-back past the peak): the last
        # converged increment is then < 1e-4 of the step from the instability
        rows = [l.split() for l in open(os.path.join(d, "p.sta")) if l.strip()[:1].isdigit()]
        conv = [r_ for r_ in rows if not r_[2].endswith("U")]
        last_dt = float(conv[-1][6]) if conv else 1.0
        ck.append(Check("W/l_ch=%g: peak captured (falls after it, or limit point within 1e-4)" % (w / LCH),
                        ipk < len(c) - 1 or last_dt < 1e-4, True, cmp="true",
                        why="rc=%d, last converged increment %.1e of the step" % (rc, last_dt)))
        sig[w] = fm / (w * (w / 40)) / TN0
        ck.append(Check("W/l_ch=%g: sigma_N/Tn0" % (w / LCH), sig[w], None, cmp="info"))
    ck.append(Check("D/l_ch=0.1: sigma_N / net-section strength Tn0(1-a/W)", sig[1.0], 1 - AW, rtol=0.05))
    r = sig[100.0] / (lefm(100.0) / TN0)
    ck.append(Check("D/l_ch=10: sigma_N / sigma_LEFM (from below)", r, 0.85, cmp="ge"))
    ck.append(Check("D/l_ch=10: sigma_N / sigma_LEFM <= 1", r, 1.0, cmp="le"))
    ck.append(Check("D/l_ch=10: equivalent crack extension c_f / l_ch", (1 / r ** 2 - 1) * AW * 100.0 / LCH, None,
                    cmp="info", why="sigma_N = K_IC/sqrt(pi(a+c_f)) with F held"))
    ws = sorted(sig)
    ck.append(Check("sigma_N decreases with size", all(sig[a] > sig[b] for a, b in zip(ws, ws[1:])), True,
                    cmp="true"))
    # Bazant: 1/s^2 = (1 + D/d0)/B^2, least squares in D
    xs, ys = ws, [1 / sig[w] ** 2 for w in ws]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
    a = my - b * mx
    bb, d0 = 1 / math.sqrt(a), a / b
    ck.append(Check("Bazant fit: B", bb, None, cmp="info"))
    ck.append(Check("Bazant fit: d0 / l_ch", d0 / LCH, None, cmp="info"))
    worst = max(abs(bb / math.sqrt(1 + w / d0) / sig[w] - 1) for w in ws)
    ck.append(Check("Bazant fit: worst relative misfit over the 4 sizes", worst, 0.05, cmp="le"))
    ck.append(Check("gmin Kn delta_f / Tn0 (residual traction group) << 1",
                    1e-6 * 1e8 * (2 * GC / TN0) / TN0, 0.05, cmp="le",
                    why="at 1e-3 instead of 1e-6 it is 5 and the ligament never breaks"))
    return ck
