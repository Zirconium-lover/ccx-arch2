#!/usr/bin/env python3
"""Band width as damage-weighted volume, with no threshold anywhere.

    bandwidth.py RUNDIR [--area A]

WHY NOT A COUNT OVER A THRESHOLD.  The obvious measure is "how many
integration points have D above something", and it was the measure used
here first.  It does not survive its own sensitivity test.  On one sweep,
the same runs give three different verdicts:

    threshold   local                    NONLOCAL=0.5
    D > 0.1     1.000 0.667 0.917        1.667 1.333 1.667
    D > 0.5     1.000 0.667 0.417        1.000 0.833 1.083
    D > 0.9     1.000 0.250 0.083        1.000 0.500 0.333

At 0.1 the local band holds too; at 0.9 the nonlocal band collapses too.
Only the middle row says the two behave differently.  A verdict that
depends on a number nobody derived is not a measurement of the model, and
the three rows are not even measuring one quantity: D > 0.1 includes
material that has barely begun to damage, so it tracks the pre-localisation
spread the notch imposes, while elements past D > 0.9 are deleted within a
few increments, so that count tracks the deletion rate and is a snapshot of
the crack tip rather than a width.

WHAT THIS MEASURES INSTEAD.  Crack-band theory makes its claim about the
volume over which the fracture energy is dissipated, so that is what is
measured:

    w = ( sum_e  D_e * V_e )  /  A

with V_e the element's REFERENCE volume, A the reference cross-section
where the band forms, and the sum over every element - including the ones
deleted on the way, which dissipated their full G_f and are counted at
D = 1.  No threshold enters, and nothing is maximised over history: a
deleted element stays counted, so the measure only grows and there is no
"largest simultaneous" to choose.

WHAT THE NORMALISATION ASSUMES, AND WHERE IT IS WRONG.  A is the reference
cross-section at the notch, its minimum, because that is where the band forms
and what the fracture energy is dissipated through.  The sum, however, runs
over every damaged element, including the diffuse low-level damage outside the
notch, whose own cross-section is larger than A.  Those contributions are
therefore divided by a smaller area than their own and overstate their share
of the width.  The error is bounded by the notch depth - here max/min = 1/0.81,
so at most 23 per cent, and only on the part of the sum that lies outside the
narrowed segment - and it is IDENTICAL for every run of a sweep, since the
specimen is.  So it cancels in ratios and in any comparison across a sweep, and
it does not cancel in an absolute width.  Read the absolute number as "the
length of fully damaged notch-section material that would dissipate the same
energy", which is what it is, rather than as the distance between two points.

CROSS-CHECKED BY A SECOND IMPLEMENTATION.  Agent 1 wrote this formula from
scratch without reading this file - deliberately, so that a shared mistake
could not survive - and ran it on the same sweep.  All three nonlocal/local
ratios agreed exactly (1.34, 1.72, 2.13) and every absolute value differed by
exactly 0.810, which is min/max of this deck's cross-section: they had
normalised by the gross area.  A constant factor is the signature of a
normalisation difference rather than a measurement one, and it cancels in
precisely the ratios the conclusion rests on.

The inputs are files the solver already writes.  The VTK snapshot carries
per-cell DE1_D with an ELEMENT_ID array, so a cell needs no reconstruction
to be identified; the .damage history lists what was deleted and when; the
deck supplies the reference coordinates.  Deformed coordinates are in the
VTK too and are deliberately NOT used - a band measured in the deformed
configuration grows with the stretch inside it, and arms that rupture at
different displacements would then be compared at different stretches.
"""
import os, re, sys


def deck_geometry(path):
    """Reference coordinates and C3D4 connectivity, from the deck itself."""
    co, el, mode = {}, {}, None
    for line in open(path):
        t = line.strip()
        if t.startswith('**') or not t:
            continue
        if t.startswith('*'):
            u = t.upper()
            if u.startswith('*NODE'):
                mode = 'n'
            elif u.startswith('*ELEMENT'):
                mode = 'e' if 'C3D4' in u else None
            else:
                mode = None
            continue
        f = [v.strip() for v in t.split(',') if v.strip() != '']
        if mode == 'n' and len(f) >= 4:
            co[int(f[0])] = (float(f[1]), float(f[2]), float(f[3]))
        elif mode == 'e' and len(f) >= 5:
            el[int(f[0])] = [int(v) for v in f[1:5]]
    return co, el


def vol(co, nodes):
    p, q, r, s = (co[n] for n in nodes)
    a = [q[i] - p[i] for i in range(3)]
    b = [r[i] - p[i] for i in range(3)]
    c = [s[i] - p[i] for i in range(3)]
    det = (a[0] * (b[1] * c[2] - b[2] * c[1])
           - a[1] * (b[0] * c[2] - b[2] * c[0])
           + a[2] * (b[0] * c[1] - b[1] * c[0]))
    return abs(det) / 6.0


def vtk_damage(path):
    """{element id: D} from the snapshot, by its own ELEMENT_ID array."""
    lines = open(path).read().split('\n')
    n, dam, eid = None, [], []
    i = 0
    while i < len(lines):
        s = lines[i]
        if s.startswith('CELL_DATA'):
            n = int(s.split()[1])
        elif s.startswith('SCALARS DE1_D') and n:
            dam = [float(x) for x in lines[i + 2:i + 2 + n]]
            i += 1 + n
        elif s.startswith('SCALARS ELEMENT_ID') and n:
            eid = [int(x) for x in lines[i + 2:i + 2 + n]]
            i += 1 + n
        i += 1
    if not eid or len(eid) != len(dam):
        raise SystemExit("%s: no usable ELEMENT_ID/DE1_D pair" % path)
    return dict(zip(eid, dam))


def deleted(path):
    out = set()
    if not os.path.exists(path):
        return out
    for line in open(path):
        if line.startswith('#'):
            continue
        f = line.split()
        if f:
            out.add(int(f[0]))
    return out


def notch_area(deckpath):
    """The reference cross-section where the band forms.

    Read from the deck's own header rather than assumed: the tent notch
    makes the area a function of x, and the band forms at its minimum.
    """
    red, a = 0.0, 1.0
    for line in open(deckpath):
        if not line.startswith('**'):
            break
        m = re.search(r"specimen: *[\d.eE+-]+ *x *([\d.eE+-]+) *x *([\d.eE+-]+)",
                      line)
        if m:
            a = float(m.group(1)) * float(m.group(2))
        m = re.search(r"reduced by up to ([\d.eE+-]+)", line)
        if m:
            red = float(m.group(1))
    return a * (1.0 - red) ** 2


def width(rundir, area=None):
    deck = os.path.join(rundir, 't.inp')
    co, el = deck_geometry(deck)
    dam = vtk_damage(os.path.join(rundir, 't.de1.vtk'))
    gone = deleted(os.path.join(rundir, 't.damage'))
    if area is None:
        area = notch_area(deck)
    tot = 0.0
    for e, nodes in el.items():
        d = 1.0 if e in gone else dam.get(e, 0.0)
        tot += d * vol(co, nodes)
    missing = [e for e in el if e not in gone and e not in dam]
    return tot / area, area, len(el), len(gone), len(missing)


def stage_ok(rundir):
    """Did this run reach the same stage as the others - full rupture?

    WHY A WIDTH NEEDS THIS AND A WORK INTEGRAL DOES NOT.  An external-work
    integral can be taken to a window every arm reached, so an arm that
    stopped early is compared on equal terms.  A width has no such window:
    it is read off whatever state the run got to, and an arm that stopped
    before rupture contributes the width of a half-formed band with nothing
    in the number to say so.  That produced a non-monotone scaling law once
    (0.882, 0.741, 1.582 in ell, with the middle arm stopped) and, separately,
    a published agreement of -4 per cent from an arm that had reached 15 per
    cent of its step.

    THREE SIGNS, and they were found by two people independently, which is
    why all three are kept rather than the shortest set.  rc!=0 is the
    solver's own verdict, recorded by the sweep as an UNCONVERGED marker
    because the exit status is not in any output file.  A step time short of
    1.0, and the U suffix the solver puts on an unconverged attempt, are two
    ways a run fails while still exiting 0.  And nothing deleted means the
    bar was loaded but never cut, which all of the others miss.

    Reads t.sta as STEP INC ATT ITRS TOT_TIME STEP_TIME INC_TIME - column 6
    is the progress and column 3 is the attempt counter.  Testing column 3
    against 1.0 and calling it theta is what this function did first; it
    passes a finished run only because the last increment usually converges
    in one attempt.

    This lives in the module both sweeps import, and not in either of them,
    because the first version existed twice and the two copies had already
    started to differ.

    WHAT THIS IS NOT.  It asks whether a run reached FULL RUPTURE, which is
    the right comparability criterion only where the comparison needs it -
    a band width, which has no common window.  It is NOT a general validity
    test, and applying it where the comparison is made at a common
    displacement will reject sound runs: the energy-equivalence decks are
    driven to a tenth of working strain and are not meant to break at all,
    and the cross-section sweep compares force-displacement curves on the
    window every arm reached.  Those are handled by work(c, umax) instead.
    Read a False from this as "not comparable BY RUPTURE", nothing wider.
    """
    if os.path.exists(os.path.join(rundir, 'UNCONVERGED')):
        return False, 'solver exited non-zero'
    att = stime = None
    try:
        for line in open(os.path.join(rundir, 't.sta')):
            f = line.split()
            if len(f) > 5 and f[0].isdigit():
                att, stime = f[2], f[5]
    except Exception:
        return False, 'no t.sta'
    if stime is None:
        return False, 'no step record'
    if 'U' in att:
        return False, 'last attempt did not converge (att=%s)' % att
    try:
        if abs(float(stime) - 1.0) > 1.0e-9:
            return False, 'step reached %.4f of 1.0' % float(stime)
    except ValueError:
        return False, 'unreadable step time %r' % stime
    ndel = 0
    try:
        for line in open(os.path.join(rundir, 't.damage')):
            if not line.startswith('#'):
                ndel += 1
    except Exception:
        pass
    if ndel == 0:
        return False, 'step completed but nothing deleted, bar never cut'
    return True, 'step complete, %d deleted' % ndel


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    area = None
    for a in sys.argv[1:]:
        if a.startswith('--area'):
            area = float(a.split('=', 1)[1])
    if not args:
        raise SystemExit(__doc__)
    for d in args:
        w, a, ne, ng, nm = width(d, area)
        print("%-40s w=%.4f  (A=%.4f, %d elements, %d deleted%s)"
              % (d, w, a, ne, ng,
                 ", %d UNACCOUNTED" % nm if nm else ""))
    return 0


if __name__ == '__main__':
    sys.exit(main())
