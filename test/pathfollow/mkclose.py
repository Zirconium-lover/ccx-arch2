#!/usr/bin/env python3
"""Close a DAMAGED cohesive facet, so the NORMAL law can be checked directly.

The mixed-mode and snap-back benchmarks in this directory only ever OPEN
their facet, so the compression branch of `cohesive_uc6.f` has no analytical
test at all - and that branch is where the law's tangent jumps by `1/gmin`,
which was measured to be a wall.  A fix to it was checked only by its effect
on a structure, which is evidence about the fix and not about the law.

This deck fixes that.  It reuses the bar from `mkcohesive.py` verbatim by
importing it, so the two specimens cannot drift apart, and replaces only the
step cards with two steps:

  1. pull the single facet PAST failure, so `g` reaches `gmin` and the
     tangent jump is at its full five orders;
  2. push the same facet back through zero separation and deep into contact,
     far enough that the far-field compression slope is `Kn` and not still
     inside any blend band.

Every increment prints what `check_close.py` needs and nothing else has to be
inferred: the separation from `U` at the two coincident interface node sets,
the damage from `SDV 2`, the traction from `S` on the facet.  The identity
checked is on the constitutive routine alone - no bulk model, no structural
branch, no finite-strain correction enters it.

    ./mkclose.py -o close.inp
    ccx_2.23_pardiso -i close
    ./check_close.py <rundir> --zeta 0
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mkcohesive


DOC = """** This deck exists to verify the NORMAL law of cohesive_uc6.f directly,
** not through a structural response.  Step 1 pulls the single facet past
** failure; step 2 pushes it back through zero separation into contact.
**
**   sharp     T = kn*d                        d < 0
**             T = g*kn*d                      d >= 0
**   blended   T = g*kn*d - (1-g)*kn*psi(-d)
**             psi(u) = 0, u<=0 ; u^2/(2w), 0<u<w ; u-w/2, u>=w ; w = zeta*d0
**
** with g = max(gmin, 1-dvisc) from SDV 2.  See check_close.py."""


def build_close(a, iface_minus, iface_plus):
    """The bar from mkcohesive, with the step cards replaced."""
    out = []
    for ln in mkcohesive.build(a).split("\n"):
        if ln.startswith("*STEP"):
            break
        out.append(ln)
    out.append("** the coincident interface nodes: the minus side, then the")
    out.append("** duplicates the right block uses.")
    out.append("*NSET, NSET=NIFM")
    out.append(", ".join(str(x) for x in iface_minus))
    out.append("*NSET, NSET=NIFP")
    out.append(", ".join(str(x) for x in iface_plus))
    out.append(DOC)
    for uend, dt in ((a.upull, a.dtpull), (-a.uclose, a.dtclose)):
        out.append("*STEP, INC=%d, NLGEOM" % a.inc)
        out.append("*STATIC")
        out.append("%g, 1.0, %g, %g" % (dt, a.dtmin, dt))
        out.append("*BOUNDARY")
        out.append("NRIGHT, 1, 1, %g" % uend)
        for nset in ("NIFM", "NIFP"):
            out.append("*NODE PRINT, NSET=%s, FREQUENCY=1" % nset)
            out.append("U")
        out.append("*EL PRINT, ELSET=ECOH, FREQUENCY=1")
        out.append("S, SDV")
        out.append("*END STEP")
    return "\n".join(out) + "\n"


def interface_nodes(nslice, isplit):
    """The node numbers mkcohesive.build assigns to the interface, both sides.

    Recomputed here the way it numbers them rather than parsed back out of
    the deck, so a change to its numbering breaks this loudly instead of
    silently printing the wrong nodes."""
    n = 0
    nid = {}
    for i in range(nslice + 1):
        for j in (0, 1):
            for k in (0, 1):
                n += 1
                nid[(i, j, k)] = n
    minus = [nid[(isplit, j, k)] for j in (0, 1) for k in (0, 1)]
    plus = [n + 1 + t for t in range(4)]
    return minus, plus


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-o", "--out", default="close.inp")
    for k in ("e", "nu", "kn", "tn0", "ts0", "gc", "gmin", "mu", "h"):
        p.add_argument("--" + k, type=float, default=getattr(mkcohesive, k.upper()))
    p.add_argument("--nslice", type=int, default=4)
    p.add_argument("--isplit", type=int, default=2)
    p.add_argument("--upull", type=float, default=1.2e-2,
                   help="end displacement of step 1; past df so g reaches gmin")
    p.add_argument("--uclose", type=float, default=1.2e-1,
                   help="how far into compression step 2 pushes.  It has to be "
                        "large: the closed facet is Kn-stiff, so nearly all of "
                        "it goes into the bar and only u/(1+Kn*L/E) reaches "
                        "the facet")
    p.add_argument("--dtpull", type=float, default=2.0e-2)
    p.add_argument("--dtclose", type=float, default=5.0e-3)
    p.add_argument("--inc", type=int, default=4000)
    p.add_argument("--dtmin", type=float, default=1e-9)
    # mkcohesive.build reads these; they only affect the step cards it writes,
    # which this deck replaces, but they must exist.
    p.add_argument("--uend", type=float, default=mkcohesive.UEND)
    p.add_argument("--dt0", type=float, default=0.01)
    p.add_argument("--dtmax", type=float, default=0.02)
    a = p.parse_args()
    minus, plus = interface_nodes(a.nslice, a.isplit)
    open(a.out, "w").write(build_close(a, minus, plus))
    d0 = a.tn0 / a.kn
    df = 2.0 * a.gc / a.tn0
    print("wrote %s: d0=%.6g df=%.6g; step 1 to u=%g (past df, so g -> gmin), "
          "step 2 to u=%g" % (a.out, d0, df, a.upull, -a.uclose))


if __name__ == "__main__":
    main()
