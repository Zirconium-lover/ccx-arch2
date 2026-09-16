#!/usr/bin/env python3
"""Measure crack-band objectivity as agreement BETWEEN runs of a sweep.

Usage:
  check_objectivity.py UEND LABEL:L:dat [LABEL:L:dat ...]

Every run in a sweep shares the specimen, the material and the loading, and
differs only in how the mesh is refined in directions the crack band does
not run in.  A correct model must give the same answer on all of them, so
no analytical reference is needed and no normalisation either: the
force-displacement curves are directly comparable.

Two numbers are reported.

  W = integral F du, the external work.  This is the one the crack-band
  approach makes a claim about - it exists precisely to keep the dissipated
  energy independent of element size - so the spread of W across the sweep
  is the headline result.  Each run's W is also shown against 1/L, because
  if the characteristic length is the cause then W must scale as 1/L.

  The pointwise spread of F(u), split at the peak.  The pre-peak spread
  should be ~0 on any mesh that represents the uniform field exactly; it is
  reported so that a post-peak spread cannot be dismissed as the two curves
  merely being sampled differently.
"""
import re, sys


def curve(path, uend, area=1.0):
    out, pend = [], None
    for line in open(path):
        m = re.search(r"total force .* and time\s+([-\dEe.+]+)", line)
        if m:
            pend = float(m.group(1))
            continue
        if pend is not None:
            f = line.split()
            if len(f) == 3:
                try:
                    out.append((uend * pend, float(f[0]) / area))
                except ValueError:
                    pass
                pend = None
    return out


def interp(c, u):
    if u <= c[0][0]:
        return c[0][1]
    if u >= c[-1][0]:
        return c[-1][1]
    for i in range(1, len(c)):
        if c[i][0] >= u:
            u0, s0 = c[i - 1]
            u1, s1 = c[i]
            return s1 if u1 == u0 else s0 + (s1 - s0) * (u - u0) / (u1 - u0)
    return c[-1][1]


def work(c):
    return sum(0.5 * (c[i][1] + c[i - 1][1]) * (c[i][0] - c[i - 1][0])
               for i in range(1, len(c)))


def main():
    uend = float(sys.argv[1])
    runs = []
    for spec in sys.argv[2:]:
        label, ell, path = spec.split(":", 2)
        runs.append((label, float(ell), curve(path, uend)))
    if len(runs) < 2:
        print("need at least two runs")
        return 2

    peak = max(max(f for _, f in c) for _, _, c in runs)
    upk = max(max(c, key=lambda p: p[1])[0] for _, _, c in runs)
    umax = min(c[-1][0] for _, _, c in runs)

    ws = [work(c) for _, _, c in runs]
    base, lbase = ws[0], runs[0][1]
    print("  W = integral F du   (the quantity crack-band scaling is for)")
    for (label, ell, _), w in zip(runs, ws):
        print("    %-10s L=%.6f  W=%10.5f  W/W0=%6.3f   1/L ratio=%6.3f"
              % (label, ell, w, w / base, lbase / ell))
    print("    SPREAD of W across the sweep : %6.2f %% of W0"
          % (100.0 * (max(ws) - min(ws)) / base))

    pre = post = 0.0
    for i in range(1, 801):
        u = umax * i / 800
        v = [interp(c, u) for _, _, c in runs]
        s = max(v) - min(v)
        if u <= upk:
            pre = max(pre, s)
        else:
            post = max(post, s)
    print("  pointwise spread of F(u), peak F=%.4f at u=%.6f, over u<=%.6f"
          % (peak, upk, umax))
    print("    before peak : %9.5f = %6.2f %% of peak" % (pre, 100 * pre / peak))
    print("    AFTER  peak : %9.5f = %6.2f %% of peak" % (post, 100 * post / peak))
    return 0


if __name__ == "__main__":
    sys.exit(main())
