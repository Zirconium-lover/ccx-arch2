#!/usr/bin/env python3
"""Measure crack-band objectivity as agreement BETWEEN runs of a sweep.

Usage:
  check_objectivity.py UEND LABEL:L:dat [LABEL:L:dat ...]

Every run in a sweep shares the specimen, the material and the loading, and
differs only in how the mesh is refined in directions the crack band does
not run in.  A correct model must give the same answer on all of them, so
no analytical reference is needed and no normalisation either: the
force-displacement curves are directly comparable.

Three numbers are reported.

  W = integral F du, the external work.  This is the one the crack-band
  approach makes a claim about - it exists precisely to keep the dissipated
  energy independent of element size - so the spread of W across the sweep
  is the headline result.  Each run's W is also shown against 1/L, because
  if the characteristic length is the cause then W must scale as 1/L.

  The fit of W = A + B/L across the sweep, which turns that "if" into a
  test.  A spread caused by the characteristic length lies ON this law; a
  spread that does not lie on it has some other cause, however big it is.
  Reporting the spread alone invites the spread to be blamed on the width
  by default, which is wrong in the one sweep this file exists for: refine
  a cross-section and the band stays one slice thick, so under a width
  that is measured across the band, L does not move at all and the law
  predicts no spread whatever.  The fit says so in that case instead of
  being computed.

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


def work(c, umax=None):
    """External work, integrated to umax rather than to wherever the run died.

    WHY THE LIMIT.  Without it this returns less work for an arm that
    stopped early simply BECAUSE it stopped early, and the caller then
    compares runs at different stages of failure while believing it is
    comparing the same quantity.  Every headline number in this file so
    far came from arms that all reached the end, so the omission never
    showed - which is exactly how a latent trap survives.

    It stops being latent the moment a nonlocal arm fails to converge
    where its local reference does, which is the normal case for the
    experiment this tool is now being asked to support.
    """
    tot = 0.0
    for i in range(1, len(c)):
        u0, s0 = c[i - 1]
        u1, s1 = c[i]
        if umax is not None:
            if u0 >= umax:
                break
            if u1 > umax:                     # partial trapezium to umax
                s1 = s0 + (s1 - s0) * (umax - u0) / (u1 - u0)
                u1 = umax
        tot += 0.5 * (s0 + s1) * (u1 - u0)
    return tot


def lawfit(runs, ws):
    """Fit W = A + B/L, the crack-band scaling, across the sweep.

    Crack-band theory fixes the dissipation per unit volume at G_f/L, so
    the part of the external work the band dissipates must vary as 1/L
    while the rest of the response does not (Bazant and Oh 1983).  That
    makes the fit a test of ATTRIBUTION and not just of size: a spread
    that is the characteristic length lies on this law, and a spread that
    does not lie on it is something else, however large it is.

    A sweep whose L never changes cannot be fitted and does not need to
    be.  There the law predicts no spread at all, so whatever spread is
    measured is by construction not the width - which is the whole point
    of refining a cross-section, where the band stays one slice thick.
    """
    x = [1.0 / ell for _, ell, _ in runs]
    print("  W = A + B/L, the crack-band scaling (Bazant and Oh 1983)")
    if max(x) - min(x) <= 1.0e-12 * max(1.0, max(x)):
        print("    L is the same on every mesh (1/L = %.4f), so this law"
              % max(x))
        print("    predicts a spread of exactly ZERO.  Whatever spread is")
        print("    reported above is therefore NOT the band width.")
        return
    n = len(x)
    sx, sy = sum(x), sum(ws)
    sxx = sum(v * v for v in x)
    sxy = sum(a * b for a, b in zip(x, ws))
    den = n * sxx - sx * sx
    if abs(den) <= 1.0e-30:
        print("    the sweep does not span enough of 1/L to fit")
        return
    b = (n * sxy - sx * sy) / den
    a = (sy - b * sx) / n
    worst = 0.0
    for (label, ell, _), w in zip(runs, ws):
        f = a + b / ell
        r = abs(w - f) / w * 100.0
        worst = max(worst, r)
        print("    %-10s 1/L=%6.3f  W=%10.5f  fit=%10.5f  off by %5.2f %%"
              % (label, 1.0 / ell, w, f, r))
    print("    A=%10.4f  the mesh independent part" % a)
    print("    B=%10.4f  the part dissipated in the band" % b)
    if n < 3:
        print("    two points fit a two-parameter law exactly, so the")
        print("    departure below is not evidence; sweep at least three.")
    print("    WORST departure from the 1/L law : %6.2f %%" % worst)


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

    ws = [work(c, umax) for _, _, c in runs]
    base, lbase = ws[0], runs[0][1]
    ends = [c[-1][0] for _, _, c in runs]
    print("  W = integral F du to u=%.6f, the FURTHEST ALL ARMS REACHED" % umax)
    if max(ends) - min(ends) > 1e-12:
        print("    arms ended at %s, so the comparison is made on the"
              % ", ".join("%.4f" % e for e in ends))
        print("    common window; an arm that stopped early is not credited")
        print("    with less work for having stopped.")
    for (label, ell, _), w in zip(runs, ws):
        print("    %-10s L=%.6f  W=%10.5f  W/W0=%6.3f   1/L ratio=%6.3f"
              % (label, ell, w, w / base, lbase / ell))
    print("    SPREAD of W across the sweep : %6.2f %% of W0"
          % (100.0 * (max(ws) - min(ws)) / base))
    lawfit(runs, ws)

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
