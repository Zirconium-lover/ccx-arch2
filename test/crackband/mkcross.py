#!/usr/bin/env python3
"""Generate the crack-band CROSS-SECTION REFINEMENT objectivity deck.

The question
------------
The DE1/DM2 softening law spreads the fracture energy over a band of width
L, D = D_base + L*eps_p/u_f, so L must be the width of the band that
actually forms.  Does the global response depend on how the mesh is
refined in directions the band does not run in?  It must not.

The specimen is FIXED
---------------------
A prismatic bar along x, length Ltot = nslice*t, cross-section h x h.
Geometry, material, loading and the number of slices along x are the same
in every run of the sweep, so the band that forms - one slice thick, t -
is the same physical band.  ONLY the cross-section discretisation changes:
each slice is cut into ncross x ncross cells, each cell into the standard
6-tetrahedron (Kuhn) decomposition.

    element size  = t  x  (h/ncross)  x  (h/ncross)
    true band width                  = t          (independent of ncross)
    legacy   L = (6V)^(1/3) = (t*(h/ncross)^2)^(1/3)   SHRINKS with ncross
    projected L                      = t          for every Kuhn tet,
        because all six share the c0-c6 diagonal and so span the full
        slice thickness along x.

This is why the sweep is decisive.  Refining the cross-section cannot
change the answer of a correct model - it does not even change the band.
The legacy length nevertheless falls as ncross^(-2/3), so the damage
D = L*eps_p/u_f accumulates more slowly on the finer mesh and the bar
appears TOUGHER the more it is refined.  With the projection L stays t and
the curves must lie on top of one another.

Note that at ncross = 1 the cell is a cube and the two lengths AGREE
exactly, which is the case a naive refinement study would report as
"objective".  The defect only appears once the element stops being a cube
- which is the point Jirasek and Bauer 2012 make in section 5.3.1, and the
error they size as "comparable to a misprediction of the fracture energy by
50 % or even more" (section 7).

No closed form is claimed.  Objectivity here is agreement BETWEEN runs of
the sweep, which needs no analytical reference: the runs share a specimen,
so their force-displacement curves are directly comparable with no
normalisation at all.  r = L_e*sigma_0/(E*u_f) is kept below 1 so plain
displacement control traces the branch and no path-following method enters
the measurement.
"""
import argparse

E, NU = 200000.0, 0.0        # nu=0: no elastic lateral mismatch
SIGY, HMOD, EPSH = 400.0, 800.0, 0.05    # ONE law for every slice
EPS0 = 0.01                  # Rice-Tracey reference strain -> initiation
UF = 0.02                    # DE1 failure plastic displacement
T, H = 1.0, 1.0              # slice thickness along x, cross-section edge
NSLICE, IWEAK, NCROSS = 6, 3, 1
UEND = 0.30


def kuhn(c):
    v0, v1, v2, v3, v4, v5, v6, v7 = c
    return [(v0, v1, v2, v6), (v0, v2, v3, v6), (v0, v3, v7, v6),
            (v0, v7, v4, v6), (v0, v4, v5, v6), (v0, v5, v1, v6)]


def signed_vol6(p, q, r, s):
    ax, ay, az = q[0] - p[0], q[1] - p[1], q[2] - p[2]
    bx, by, bz = r[0] - p[0], r[1] - p[1], r[2] - p[2]
    cx, cy, cz = s[0] - p[0], s[1] - p[1], s[2] - p[2]
    return (ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx)
            + az * (bx * cy - by * cx))


def predict(a):
    cell = a.h / a.ncross
    sig0 = a.sigy + a.hmod * a.eps0
    ltot = a.nslice * a.t
    return dict(cell=cell, sig0=sig0, ltot=ltot,
                l_legacy=(a.t * cell * cell) ** (1.0 / 3.0), l_proj=a.t,
                u_peak=sig0 * ltot / a.e + ltot * a.eps0,
                r=(a.nslice - 1) * a.t * sig0 / (a.e * a.uf))


def build(a):
    p = predict(a)
    nc, cell = a.ncross, p["cell"]
    nid, coord, n = {}, {}, 0
    for i in range(a.nslice + 1):
        for j in range(nc + 1):
            for k in range(nc + 1):
                n += 1
                nid[(i, j, k)] = n
                coord[n] = (i * a.t, j * cell, k * cell)

    bulk, weak, eid = [], [], 0
    for i in range(a.nslice):
        for j in range(nc):
            for k in range(nc):
                c = [nid[(i, j, k)], nid[(i + 1, j, k)],
                     nid[(i + 1, j + 1, k)], nid[(i, j + 1, k)],
                     nid[(i, j, k + 1)], nid[(i + 1, j, k + 1)],
                     nid[(i + 1, j + 1, k + 1)], nid[(i, j + 1, k + 1)]]
                if a.eltype == "C3D8":
                    # The cell IS the element.  Its extent along x is the
                    # slice thickness, so the projected width is t here too,
                    # while the legacy length is not even defined for this
                    # family - which is the point of running the sweep on it.
                    eid += 1
                    (weak if i == a.iweak else bulk).append((eid, c))
                else:
                    for tet in kuhn(c):
                        eid += 1
                        tet = list(tet)
                        if signed_vol6(*[coord[x] for x in tet]) < 0.0:
                            tet[1], tet[2] = tet[2], tet[1]
                        (weak if i == a.iweak else bulk).append((eid, tet))

    left = [nid[(0, j, k)] for j in range(nc + 1) for k in range(nc + 1)]
    right = [nid[(a.nslice, j, k)] for j in range(nc + 1)
             for k in range(nc + 1)]

    L = []
    add = L.append
    add("** Crack-band CROSS-SECTION REFINEMENT objectivity deck.")
    add("** Generated by test/crackband/mkcross.py - do not edit by hand.")
    add("**")
    add("**   specimen: %g x %g x %g   slices=%d  weak slice=%d"
        % (p["ltot"], a.h, a.h, a.nslice, a.iweak))
    add("**   eltype=%s  ncross=%d -> cell %g x %g x %g, %d elements,"
        " %d nodes" % (a.eltype, nc, a.t, cell, cell, eid, n))
    add("**   true band width            = %.10g  (same on every mesh)" % a.t)
    if a.eltype == "C3D4":
        add("**   legacy    L=(6V)^(1/3)     = %.10g" % p["l_legacy"])
    else:
        add("**   legacy    L                = NOT DEFINED for %s -"
            " CCX_DAMAGE_CHARLEN=0 refuses this deck" % a.eltype)
    add("**   projected L                = %.10g" % p["l_proj"])
    add("**   snap-back indicator r      = %.4f  (<1, displacement control)"
        % p["r"])
    add("**   exact peak end displacement u = %.8f  (precedes damage)"
        % p["u_peak"])
    add("**   evolution=%s  sigma_0=%.6f  u_f=%.6f  G_f=%.6f"
        % (a.evolution, p["sig0"], a.uf, p["sig0"] * a.uf / 2.0))
    add("**   The SPECIMEN is identical across the sweep, so the")
    add("**   force-displacement curves are directly comparable with no")
    add("**   normalisation; any spread is mesh dependence and nothing else.")
    add("*NODE, NSET=NALL")
    for i in range(1, n + 1):
        add("%d, %.10g, %.10g, %.10g" % ((i,) + coord[i]))
    fmt = "%d," + ",".join(["%d"] * (8 if a.eltype == "C3D8" else 4))
    for setname, elems in (("EBULK", bulk), ("EWEAK", weak)):
        add("*ELEMENT, TYPE=%s, ELSET=%s" % (a.eltype, setname))
        for e, nodes in elems:
            add(fmt % ((e,) + tuple(nodes)))

    def nset(name, ids):
        add("*NSET, NSET=%s" % name)
        for i in range(0, len(ids), 8):
            add(", ".join(str(x) for x in ids[i:i + 8]))

    nset("NLEFT", left)
    nset("NRIGHT", right)
    add("*NSET, NSET=NPIN")
    add("%d" % nid[(0, 0, 0)])
    add("*NSET, NSET=NROLL")
    add("%d" % nid[(0, nc, 0)])
    sigh = a.sigy + a.hmod * a.epsh
    add("** EVERY slice carries the same elastic-plastic law, so the whole")
    add("** bar flows together and the field stays uniformly uniaxial.  A")
    add("** single yielding slice between elastic blocks would be laterally")
    add("** constrained and carry an axial stress well above sigma_y - and")
    add("** nu=0 does NOT avoid that, because J2 plastic flow is")
    add("** incompressible, so a yielding slice contracts laterally while")
    add("** elastic neighbours do not.  The slices differ ONLY in that the")
    add("** weak one may damage, which is what localises failure.")
    for name, dam in (("BULK", False), ("WEAK", True)):
        add("*MATERIAL, NAME=%s" % name)
        add("*ELASTIC")
        add("%g, %g" % (a.e, a.nu))
        add("*PLASTIC")
        add("%g, 0." % a.sigy)
        add("%g, %g" % (sigh, a.epsh))
        if dam:
            if a.evolution == "ENERGY":
                # G_f = sigma_0*u_f/2 for the linear law, with sigma_0 the
                # flow stress at initiation.  Passing the EQUIVALENT G_f
                # makes the two cards describe the same material, which is
                # what the equivalence check tests.
                add("*DAMAGE INITIATION, CRITERION=RICETRACEY,"
                    " EVOLUTION=ENERGY")
                gf = a.gf if a.gf > 0.0 else p["sig0"] * a.uf / 2.0
                add("%g, 1.0, %g, 0." % (a.eps0, gf))
            else:
                add("*DAMAGE INITIATION, CRITERION=RICETRACEY,"
                    " EVOLUTION=DISPLACEMENT")
                add("%g, 1.0, %g, 0." % (a.eps0, a.uf))
    add("*SOLID SECTION, ELSET=EBULK, MATERIAL=BULK")
    add("*SOLID SECTION, ELSET=EWEAK, MATERIAL=WEAK")
    add("*BOUNDARY")
    add("NLEFT, 1, 1, 0.")
    add("NPIN, 2, 3, 0.")
    add("NROLL, 3, 3, 0.")
    add("*STEP, INC=%d, NLGEOM" % a.inc)
    add("*STATIC")
    add("%g, 1.0, %g, %g" % (a.dt0, a.dtmin, a.dtmax))
    add("*BOUNDARY")
    add("NRIGHT, 1, 1, %g" % a.uend)
    add("*NODE PRINT, NSET=NRIGHT, TOTALS=ONLY")
    add("RF")
    add("*END STEP")
    return "\n".join(L) + "\n"


def main():
    q = argparse.ArgumentParser()
    q.add_argument("-o", "--out", default="cross.inp")
    q.add_argument("--evolution", choices=("DISPLACEMENT", "ENERGY"),
                   default="DISPLACEMENT")
    q.add_argument("--eltype", choices=("C3D4", "C3D8"), default="C3D4")
    q.add_argument("--gf", type=float, default=-1.0,
                   help="fracture energy for EVOLUTION=ENERGY; "
                        "default derives the equivalent of --uf")
    for nm, tp, df in (("e", float, E), ("nu", float, NU),
                       ("sigy", float, SIGY), ("hmod", float, HMOD),
                       ("epsh", float, EPSH), ("eps0", float, EPS0),
                       ("uf", float, UF), ("t", float, T), ("h", float, H),
                       ("nslice", int, NSLICE), ("iweak", int, IWEAK),
                       ("ncross", int, NCROSS), ("uend", float, UEND),
                       ("inc", int, 20000), ("dt0", float, 0.005),
                       ("dtmin", float, 1e-10), ("dtmax", float, 0.01)):
        q.add_argument("--" + nm, type=tp, default=df)
    a = q.parse_args()
    open(a.out, "w").write(build(a))
    p = predict(a)
    print("wrote %s  ncross=%d cell=%.5f  L_legacy=%.6f L_proj=%.6f  r=%.3f"
          % (a.out, a.ncross, p["cell"], p["l_legacy"], p["l_proj"], p["r"]))


if __name__ == "__main__":
    main()
