#!/usr/bin/env python3
"""Was CCX_DAMAGE_NLWIDTH inside its domain on this run?

    check_nlwidth_domain.py RUNDIR [RUNDIR ...]

WHY THIS EXISTS.  The substitution makes the softening law charge 2*ell as the
band's width, and that accounting is right only while the band really is at
least that wide.  Measured on this tree (test/crackband/README.md, forum
2026-09-17), at h=0.25 with the three arms that reached rupture:

    w/(2 ell) = 1.07   the dissipation is right to  +1 %
                0.88                                +42 %
                0.77                                +95 %

The error is monotone in w/(2*ell) and crosses zero where the band stops being
as wide as the length the law charges.  So the condition is w >= 2*ell.

AND IT CANNOT BE A GUARD.  w is the width of the band that formed, which is
known only after the run; the criterion the solver CAN check - whether the mesh
resolves ell - is satisfied on all four of those arms, including the one wrong
by 95 per cent.  A switch whose validity condition is unknowable at the moment
it acts needs the condition checked afterwards instead, which is this, rather
than a silent hope that the user is inside it.

It reads the run's own files and takes nothing on faith from the caller: ell
from the solver's banner, the switch's state from the switch log it prints, the
width from the damage field and the deletion history.  A run with the switch
OFF is reported as not applicable rather than as passing, because it was never
in the domain's scope.
"""
import os, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bandwidth import width, stage_ok


def banner_ell(rundir):
    """The per-material internal lengths, as the solver itself printed them."""
    ells = []
    p = os.path.join(rundir, 'run.log')
    if not os.path.exists(p):
        return ells
    for line in open(p):
        m = re.search(r"material\s+\d+\s+ell=\s*([-\d.eED+]+)", line)
        if m:
            try:
                ells.append(float(m.group(1).replace('D', 'E')))
            except ValueError:
                pass
    return ells


def switch_on(rundir):
    p = os.path.join(rundir, 'run.log')
    if not os.path.exists(p):
        return None
    for line in open(p):
        m = re.search(r"CCX_DAMAGE_NLWIDTH\s*=\s*(\S+)", line)
        if m:
            return m.group(1).strip() == '1'
    return False


def report(rundir):
    on = switch_on(rundir)
    if on is None:
        print("  %-34s no run.log - cannot say" % rundir)
        return 2
    if not on:
        print("  %-34s NLWIDTH was off - not applicable" % rundir)
        return 0
    ok, why = stage_ok(rundir)
    if not ok:
        print("  %-34s not comparable: %s" % (rundir, why))
        return 2
    ells = set(banner_ell(rundir))
    if not ells:
        print("  %-34s NLWIDTH on but no ell in the banner" % rundir)
        return 2
    if len(ells) > 1:
        print("  %-34s several lengths %s - the domain is per material and"
              % (rundir, sorted(ells)))
        print("  %-34s this check does not split by material yet" % "")
        return 2
    ell = ells.pop()
    w = width(rundir)[0]
    r = w / (2.0 * ell)
    verdict = "INSIDE" if r >= 1.0 else "OUTSIDE"
    print("  %-34s ell=%.4f  2ell=%.4f  w=%.4f  w/2ell=%.3f  %s"
          % (rundir, ell, 2.0 * ell, w, r, verdict))
    return 0 if r >= 1.0 else 1


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    # PRECEDENCE, AND NOT max().  A definite OUTSIDE (1) must outrank a
    # cannot-judge (2), because one says the answer is wrong by a known amount
    # and the other only says to look closer.  Taking the maximum of the two
    # codes reported the weaker finding and hid the stronger one, which is what
    # the first version of this function did.
    outside = unjudged = 0
    for d in sys.argv[1:]:
        r = report(d)
        if r == 1:
            outside += 1
        elif r:
            unjudged += 1
    print()
    if outside:
        print("  %d run(s) OUTSIDE: the law charged a width the specimen did"
              % outside)
        print("  not sustain, so the dissipation is too high by roughly the")
        print("  overshoot - 0.88 measured +42 %, 0.77 measured +95 %.")
        if unjudged:
            print("  %d more could not be judged." % unjudged)
        return 1
    if unjudged:
        print("  %d run(s) could not be judged; see the lines above" % unjudged)
        return 2
    print("  every applicable run was inside the substitution's domain")
    return 0


if __name__ == '__main__':
    sys.exit(main())
