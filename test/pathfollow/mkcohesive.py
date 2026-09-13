#!/usr/bin/env python3
"""Clean snap-back benchmark: elastic bar cut by ONE cohesive interface.

Why this replaces the plastic-bar benchmark
-------------------------------------------
In the earlier bar the whole specimen was elastic-plastic, so the bar
dissipated everywhere before the damaging slice localised.  Any
dissipation-based control then arms on distributed plastic work rather
than on the failure that is of interest, and the measured dG mixes two
sources.  Here the ONLY nonlinearity and the ONLY dissipation in the model
is a single zero-thickness cohesive interface (UC6).  Everything else is
linear elastic, so:

  * the limit point and the whole descending branch are known in closed
    form and are exact for this mesh, not approximate;
  * dG measured by any path-following scheme is, by construction, the
    energy released by the crack and nothing else;
  * there is no plastic history to interact with the constraint.

Model
-----
A bar along x, cross-section h x h, N unit cubes each split into the
standard 6-tetrahedron (Kuhn) decomposition.  The mesh is cut at one
interior plane: the nodes there are duplicated, the right block uses the
duplicates, and the two triangles covering the square face carry one UC6
cohesive element each.

nu = 0 and every node is restrained in y and z.  With pure axial loading
the exact solution has zero lateral displacement, so that restraint is
consistent (it introduces no spurious stress) and it makes the model an
exact 1-D chain: linear elastic bar in series with one cohesive spring.
Linear tetrahedra represent the resulting uniform strain field exactly.

Cohesive law (verified against cohesive_uc6.f, not assumed)
-----------------------------------------------------------
    d0 = Tn0/Kn ,  df = 2*Gc/Tn0
    D(d) = df*(d-d0)/(d*(df-d0))            for d0 < d < df
    T    = (1-D)*Kn*d = Tn0*(df-d)/(df-d0)

i.e. exactly linear softening from (d0,Tn0) to (df,0).

Closed form
-----------
    u(sigma) = sigma*L/E + d(sigma) ,  d = df - (df-d0)*sigma/Tn0
    du/dsigma = L/E - (df-d0)/Tn0

so the branch snaps back - u must DECREASE while the bar breaks - iff

    r = (Tn0/(df-d0)) / (E/L) > 1

(the interface softening slope is steeper than the elastic release the
bar can supply).  Key points:

    peak     sigma = Tn0 ,  u_peak = Tn0*L/E + d0
    failure  sigma = 0   ,  u_fail = df
"""

import argparse

E     = 200000.0
NU    = 0.0
KN    = 1.0e6      # cohesive penalty stiffness
TN0   = 400.0      # normal strength
TS0   = 400.0      # shear strength (not exercised in pure mode I)
GC    = 2.0        # fracture energy
GMIN  = 1.0e-6     # residual stiffness fraction
MU    = 0.0        # viscous relaxation time (0 = inviscid)
H     = 1.0
NSLICE = 20
ISPLIT = 10        # interface sits at x = ISPLIT*h
UEND  = 0.050


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
    # ---- nodes -----------------------------------------------------
    nid, coord = {}, {}
    n = 0
    for i in range(a.nslice + 1):
        for j in (0, 1):
            for k in (0, 1):
                n += 1
                nid[(i, j, k)] = n
                coord[n] = (i * a.h, j * a.h, k * a.h)

    # duplicated interface nodes, used by the RIGHT block
    dup = {}
    for j in (0, 1):
        for k in (0, 1):
            n += 1
            dup[(j, k)] = n
            coord[n] = (a.isplit * a.h, j * a.h, k * a.h)

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
    # The x=const face of the Kuhn cube is covered by the triangles
    # (v1,v2,v6) and (v1,v6,v5) in local numbering, i.e. in face
    # coordinates (0,0),(1,0),(1,1) and (0,0),(1,1),(0,1).
    tris = [((0, 0), (1, 0), (1, 1)), ((0, 0), (1, 1), (0, 1))]
    coh = []
    for tri in tris:
        m = [nid[(a.isplit, j, k)] for (j, k) in tri]
        p = [dup[(j, k)] for (j, k) in tri]
        # cohesive_uc6 builds its normal as (x2-x1) x (x3-x1); it must
        # point from the minus side to the plus side, i.e. +x here.
        pa, pb, pc = (coord[m[0]], coord[m[1]], coord[m[2]])
        e1 = [pb[i] - pa[i] for i in range(3)]
        e2 = [pc[i] - pa[i] for i in range(3)]
        nx = e1[1] * e2[2] - e1[2] * e2[1]
        if nx < 0.0:
            m[1], m[2] = m[2], m[1]
            p[1], p[2] = p[2], p[1]
        eid += 1
        coh.append((eid, m + p))

    # ---- closed form ------------------------------------------------
    d0 = a.tn0 / a.kn
    df = 2.0 * a.gc / a.tn0
    Ltot = a.nslice * a.h
    r = (a.tn0 / (df - d0)) / (a.e / Ltot)
    upk = a.tn0 * Ltot / a.e + d0
    ufl = df

    L = []
    w = L.append
    w("** Clean snap-back benchmark: elastic bar + ONE cohesive interface.")
    w("** Generated by mkcohesive.py - do not edit by hand.")
    w("**")
    w("**   E=%g nu=%g  Kn=%g Tn0=%g Ts0=%g Gc=%g gmin=%g mu=%g"
      % (a.e, a.nu, a.kn, a.tn0, a.ts0, a.gc, a.gmin, a.mu))
    w("**   d0=Tn0/Kn=%.6g   df=2Gc/Tn0=%.6g   L=%g" % (d0, df, Ltot))
    w("**   snap-back indicator r = (Tn0/(df-d0)) / (E/L) = %.4f  (>1)" % r)
    w("**   exact peak     u = %.6f  (lambda = %.6f) at sigma = %g"
      % (upk, upk / a.uend, a.tn0))
    w("**   exact failure  u = %.6f  (lambda = %.6f) at sigma = 0"
      % (ufl, ufl / a.uend))
    w("**   u must DECREASE by %.6f while the interface fails" % (upk - ufl))
    w("**")
    w("** The ONLY nonlinearity and the ONLY dissipation in this model is")
    w("** the cohesive interface: the bulk is linear elastic, so any dG a")
    w("** path-following scheme measures is crack energy and nothing else.")

    w("*NODE, NSET=NALL")
    for i in range(1, n + 1):
        x, y, z = coord[i]
        w("%d, %.10g, %.10g, %.10g" % (i, x, y, z))

    w("*ELEMENT, TYPE=C3D4, ELSET=EBULK")
    for e, t in bulk:
        w("%d, %d, %d, %d, %d" % (e, t[0], t[1], t[2], t[3]))
    w("** UC6 is a user element and must be declared before it is used;")
    w("** e_c3d_u.f dispatches on lakon(2:3)=='C6' to e_c3d_uc6.")
    w("*USER ELEMENT, TYPE=UC6, NODES=6, INTEGRATION POINTS=3, MAXDOF=3")
    w("*ELEMENT, TYPE=UC6, ELSET=ECOH")
    for e, t in coh:
        w("%d, %d, %d, %d, %d, %d, %d" % (e, t[0], t[1], t[2], t[3], t[4],
                                          t[5]))

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

    w("** y and z are held everywhere: with nu=0 and pure axial loading the")
    w("** exact solution has no lateral displacement, so this adds no")
    w("** stress and makes the model an exact 1-D chain.")
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
    w("*EL PRINT, ELSET=ECOH")
    w("SDV")
    w("*NODE FILE")
    w("U, RF")
    w("*END STEP")
    return "\n".join(L) + "\n"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-o", "--out", default="cohesive.inp")
    for k, v in [("e", E), ("nu", NU), ("kn", KN), ("tn0", TN0),
                 ("ts0", TS0), ("gc", GC), ("gmin", GMIN), ("mu", MU),
                 ("h", H), ("uend", UEND)]:
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
