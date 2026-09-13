#!/usr/bin/env python3
"""Mixed-mode snap-back benchmark: elastic bar cut by ONE INCLINED cohesive
interface.

Why a second benchmark
----------------------
`mkcohesive.py` loads its interface in pure Mode I, so it cannot tell a
normal-opening control coordinate from an effective-opening one: both
reduce to the same functional there.  The real target has substantial
shear on the active front, so the mixed-mode direction has to be
exercised by a problem whose answer is known.

Construction
------------
The specimen is a rectangular box 0<=x<=L, and the interface is an
INCLINED internal plane

    x = x0 + t*(y - W/2)        normal  n = (1,-t,0)/sqrt(1+t^2)

The mesh is a structured grid whose layers are sheared by a hat function
that is 1 at the interface layer and 0 at both ends, so the interface is
one exactly planar quad while the two loaded faces stay perpendicular to
x.  Grading does not matter: the exact solution is uniform strain inside
each block and linear tetrahedra represent that exactly.

Every node is held in y and z, and only u_x is prescribed, so the model is
still an exact 1-D chain - but the interface now sees the axial jump
resolved onto its own inclined frame:

    jump = (D,0,0)  ->  dn = D*c ,  |ds| = D*sqrt(1-c^2) ,  c = n_x

    deff = D*kappa ,  kappa = sqrt(c^2 + beta*(1-c^2)) ,  beta=(Ts0/Tn0)^2

so a FIXED, non-trivial mode mix is loaded through the whole branch.  The
fraction of deff^2 carried by shear is beta*(1-c^2)/kappa^2, printed in
the header of the generated deck.

Closed form
-----------
Axial force balance across the inclined facet (area A/c against a
cross-section A) gives

    sigma = g*Kn*kappa^2*D/c                      (bulk axial stress)
    u_end = D + sigma*L/E = D*(1 + (L/E)*(kappa^2/c)*g*Kn)

and with the bilinear law  g*Kn*deff = Tn0*(df-deff)/(df-d0),

    sigma(deff) = (kappa/c)*Tn0*(df-deff)/(df-d0)
    u_end(deff) = deff/kappa + (L/E)*sigma(deff)

    peak     deff = d0 ,  sigma = (kappa/c)*Tn0
    failure  deff = df ,  sigma = 0

    snap-back indicator  r = L*kappa^2*Tn0 / (E*c*(df-d0))   (>1 required)
"""

import argparse
import math

E      = 200000.0
NU     = 0.0
KN     = 1.0e6
TN0    = 400.0
TS0    = 300.0     # Ts0 != Tn0, so beta != 1 and shear is not a rescaling
GC     = 2.0
GMIN   = 1.0e-6
MU     = 0.0
H      = 1.0
NSLICE = 20
ISPLIT = 10
TILT   = 1.0       # tan of the interface tilt; 1.0 = 45 degrees
UEND   = 0.060


def kuhn(c):
    v0, v1, v2, v3, v4, v5, v6, v7 = c
    return [(v0, v1, v2, v6), (v0, v2, v3, v6), (v0, v3, v7, v6),
            (v0, v7, v4, v6), (v0, v4, v5, v6), (v0, v5, v1, v6)]


def vol6(p, q, r, s):
    ax, ay, az = q[0] - p[0], q[1] - p[1], q[2] - p[2]
    bx, by, bz = r[0] - p[0], r[1] - p[1], r[2] - p[2]
    cx, cy, cz = s[0] - p[0], s[1] - p[1], s[2] - p[2]
    return (ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx)
            + az * (bx * cy - by * cx))


def build(a):
    # shear of layer i: a hat that is 1 at the interface layer and 0 at
    # both loaded faces, so those stay perpendicular to x.
    def hat(i):
        if i <= a.isplit:
            return float(i) / float(a.isplit)
        return float(a.nslice - i) / float(a.nslice - a.isplit)

    def xco(i, j):
        return i * a.h + a.tilt * (j * a.h - 0.5 * a.h) * hat(i)

    # ---- nodes -----------------------------------------------------
    nid, coord = {}, {}
    n = 0
    for i in range(a.nslice + 1):
        for j in (0, 1):
            for k in (0, 1):
                n += 1
                nid[(i, j, k)] = n
                coord[n] = (xco(i, j), j * a.h, k * a.h)

    dup = {}
    for j in (0, 1):
        for k in (0, 1):
            n += 1
            dup[(j, k)] = n
            coord[n] = (xco(a.isplit, j), j * a.h, k * a.h)

    def node(i, j, k, right):
        if (i == a.isplit) and right:
            return dup[(j, k)]
        return nid[(i, j, k)]

    # ---- bulk tetrahedra -------------------------------------------
    bulk = []
    eid = 0
    for i in range(a.nslice):
        right = (i >= a.isplit)
        c = [node(i, 0, 0, right), node(i + 1, 0, 0, right),
             node(i + 1, 1, 0, right), node(i, 1, 0, right),
             node(i, 0, 1, right), node(i + 1, 0, 1, right),
             node(i + 1, 1, 1, right), node(i, 1, 1, right)]
        for t in kuhn(c):
            eid += 1
            t = list(t)
            if vol6(*[coord[x] for x in t]) < 0.0:
                t[1], t[2] = t[2], t[1]
            bulk.append((eid, t))

    # ---- cohesive elements -----------------------------------------
    tris = [((0, 0), (1, 0), (1, 1)), ((0, 0), (1, 1), (0, 1))]
    coh = []
    for tri in tris:
        m = [nid[(a.isplit, j, k)] for (j, k) in tri]
        p = [dup[(j, k)] for (j, k) in tri]
        pa, pb, pc = (coord[m[0]], coord[m[1]], coord[m[2]])
        e1 = [pb[q] - pa[q] for q in range(3)]
        e2 = [pc[q] - pa[q] for q in range(3)]
        nx = e1[1] * e2[2] - e1[2] * e2[1]
        if nx < 0.0:
            m[1], m[2] = m[2], m[1]
            p[1], p[2] = p[2], p[1]
        eid += 1
        coh.append((eid, m + p))

    # ---- closed form ------------------------------------------------
    d0 = a.tn0 / a.kn
    df = 2.0 * a.gc / a.tn0
    beta = (a.ts0 / a.tn0) ** 2
    cnx = 1.0 / math.sqrt(1.0 + a.tilt ** 2)
    kap2 = cnx ** 2 + beta * (1.0 - cnx ** 2)
    kap = math.sqrt(kap2)
    shearfrac = beta * (1.0 - cnx ** 2) / kap2
    Ltot = a.nslice * a.h
    r = Ltot * kap2 * a.tn0 / (a.e * cnx * (df - d0))

    def sig(deff):
        return (kap / cnx) * a.tn0 * (df - deff) / (df - d0)

    def uend(deff):
        return deff / kap + Ltot * sig(deff) / a.e

    upk, ufl = uend(d0), uend(df)

    L = []
    w = L.append
    w("** Mixed-mode snap-back benchmark: elastic bar + ONE INCLINED")
    w("** cohesive interface.  Generated by mkmixed.py - do not edit.")
    w("**")
    w("**   E=%g nu=%g  Kn=%g Tn0=%g Ts0=%g Gc=%g gmin=%g mu=%g"
      % (a.e, a.nu, a.kn, a.tn0, a.ts0, a.gc, a.gmin, a.mu))
    w("**   d0=Tn0/Kn=%.6g   df=2Gc/Tn0=%.6g   L=%g" % (d0, df, Ltot))
    w("**   interface tilt tan=%g (%.3f deg), n_x=%.8f"
      % (a.tilt, math.degrees(math.atan(a.tilt)), cnx))
    w("**   beta=(Ts0/Tn0)^2=%.8f  kappa^2=%.8f" % (beta, kap2))
    w("**   MIXED MODE: shear carries %.4f of deff^2 at every UC6 point"
      % shearfrac)
    w("**   snap-back indicator r = L*kappa^2*Tn0/(E*c*(df-d0)) = %.4f (>1)"
      % r)
    w("**   exact peak     deff=%.6g  u=%.8f (lambda=%.8f) sigma=%.6f"
      % (d0, upk, upk / a.uend, sig(d0)))
    w("**   exact failure  deff=%.6g  u=%.8f (lambda=%.8f) sigma=0"
      % (df, ufl, ufl / a.uend))
    w("**   u must DECREASE by %.8f while the interface fails" % (upk - ufl))
    w("**")
    w("** The bulk is linear elastic and every node is held in y and z, so")
    w("** the exact solution is uniform strain in each block with a")
    w("** constant jump across the inclined facet.  Linear tetrahedra")
    w("** represent that exactly, whatever the layer grading.")

    w("*NODE, NSET=NALL")
    for i in range(1, n + 1):
        x, y, z = coord[i]
        w("%d, %.10g, %.10g, %.10g" % (i, x, y, z))

    w("*ELEMENT, TYPE=C3D4, ELSET=EBULK")
    for e, t in bulk:
        w("%d, %d, %d, %d, %d" % (e, t[0], t[1], t[2], t[3]))
    w("*USER ELEMENT, TYPE=UC6, NODES=6, INTEGRATION POINTS=3, MAXDOF=3")
    w("*ELEMENT, TYPE=UC6, ELSET=ECOH")
    for e, t in coh:
        w("%d, %d, %d, %d, %d, %d, %d"
          % (e, t[0], t[1], t[2], t[3], t[4], t[5]))

    left = [nid[(0, j, k)] for j in (0, 1) for k in (0, 1)]
    right = [nid[(a.nslice, j, k)] for j in (0, 1) for k in (0, 1)]
    w("*NSET, NSET=NLEFT")
    w(", ".join(str(x) for x in left))
    w("*NSET, NSET=NRIGHT")
    w(", ".join(str(x) for x in right))

    w("*MATERIAL, NAME=ELAST")
    w("*ELASTIC")
    w("%g, %g" % (a.e, a.nu))
    w("*MATERIAL, NAME=COH")
    w("*ELASTIC")
    w("%g, %g" % (a.e, a.nu))
    w("*DEPVAR")
    w("4")
    w("*SOLID SECTION, ELSET=EBULK, MATERIAL=ELAST")
    w("*USER SECTION, ELSET=ECOH, MATERIAL=COH, CONSTANTS=6")
    w("%g, %g, %g, %g, %g, %g" % (a.kn, a.tn0, a.ts0, a.gc, a.gmin, a.mu))

    w("*BOUNDARY")
    w("NALL, 2, 3, 0.")
    w("NLEFT, 1, 1, 0.")

    w("*STEP, INC=%d, NLGEOM" % a.inc)
    w("*STATIC")
    w("%g, 1.0, %g, %g" % (a.dt0, a.dtmin, a.dtmax))
    w("*BOUNDARY")
    w("NRIGHT, 1, 1, %g" % a.uend)
    w("*NODE PRINT, NSET=NRIGHT, TOTALS=ONLY")
    w("RF")
    w("*NODE PRINT, NSET=NRIGHT")
    w("U")
    w("** E holds the three LOCAL separations (dn,ds1,ds2), dvisc, deff,")
    w("** dmax for a UC6 point - see resultsmech_uc6.f.")
    w("*EL PRINT, ELSET=ECOH")
    w("SDV, E")
    w("*NODE FILE")
    w("U, RF")
    w("*END STEP")
    return "\n".join(L) + "\n"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-o", "--out", default="mixed.inp")
    for k, v in [("e", E), ("nu", NU), ("kn", KN), ("tn0", TN0),
                 ("ts0", TS0), ("gc", GC), ("gmin", GMIN), ("mu", MU),
                 ("h", H), ("uend", UEND), ("tilt", TILT)]:
        p.add_argument("--" + k, type=float, default=v)
    p.add_argument("--nslice", type=int, default=NSLICE)
    p.add_argument("--isplit", type=int, default=ISPLIT)
    p.add_argument("--inc", type=int, default=5000)
    p.add_argument("--dt0", type=float, default=0.01)
    p.add_argument("--dtmin", type=float, default=1e-9)
    p.add_argument("--dtmax", type=float, default=0.02)
    a = p.parse_args()
    open(a.out, "w").write(build(a))
    print("wrote", a.out)


if __name__ == "__main__":
    main()
