#!/usr/bin/env python3
"""Generate the analytical snap-back benchmark deck.

Geometry
--------
A prismatic bar along x, cross-section h x h, built from N unit cubes each
split into the standard 6-tetrahedron (Kuhn) decomposition.  DE1 damage
evolution is implemented for C3D4 only, hence tetrahedra.

EVERY slice carries the same elastic-plastic law, so the whole bar yields
together and contracts laterally without mismatch.  That matters: a single
plastic slice sandwiched between elastic blocks is laterally constrained,
develops triaxiality and carries an axial stress well above sigma_y, which
destroys the 1-D reference.  With the bar uniformly plastic and nu = 0 the
stress field is uniformly uniaxial, and a uniform strain field is
represented EXACTLY by linear tetrahedra.  The slices differ only in that
one of them is allowed to damage, which is what localises the failure.

Closed-form 1-D reference
-------------------------
The damaging slice has edge h, so its DE1 characteristic length is
L_char = (6*V)^(1/3) = h.  After initiation the degradation grows as
D = delta_p / u_f with delta_p = h*eps_p the plastic displacement of the
slice, so the slice carries

    sigma = (1 - D) * sigma_0 ,   sigma_0 = sigma_y + Hmod*eps0

Up to initiation the bar flows uniformly, so with total length Ltot

    u_peak = sigma_0*Ltot/E + Ltot*eps0

Past initiation the damaging slice softens while the remaining L_e =
(N-1)*h unloads ELASTICALLY, hence

    u(sigma) = sigma*Ltot/E + Ltot*eps0 + u_f*(1 - sigma/sigma_0)
    du/dsigma = Ltot/E - u_f/sigma_0

so the response SNAPS BACK - the controlled end displacement must
DECREASE while the bar fails - once

    r = L_e*sigma_0 / (E*u_f) > 1 .

    failure   sigma = 0 ,  u_fail = Ltot*eps0 + u_f

With the defaults below r ~ 3.8, u_peak ~ 0.0500, u_fail ~ 0.0200: the end
displacement has to travel BACKWARDS by ~0.030 while the bar breaks.  No
monotone displacement- or step-time-controlled solver can follow that; a
path-following method can.

Note that the ELEMENT-level crack-band indicator printed by CalculiX,
L_char*sigma_y/(u_f*E) = 0.2 here, is below one: the element's own
softening is not steeper than its own elastic release.  The snap-back is
STRUCTURAL - it comes from the elastic release of the other 19 slices - so
the local indicator cannot see it.  That is deliberate.
"""

import argparse

# ---------------------------------------------------------------- defaults

E      = 200000.0   # Young's modulus
NU     = 0.0        # keeps the elastic field a uniform uniaxial strain
SIGY   = 400.0      # yield stress (elastic - perfectly plastic)
EPS0   = 2.0e-5     # Rice-Tracey reference strain -> initiation at PEEQ~eps0.
                    # Kept small on purpose: the Gutierrez constraint measures
                    # TOTAL dissipation, so a long uniform plastic plateau ahead
                    # of localisation would spend the increment budget on bulk
                    # plasticity instead of on the branch of interest.
UF     = 0.010      # DE1 failure plastic displacement
HMOD   = 800.0      # plastic hardening modulus (keeps the pre-peak branch
                    # well defined; a perfectly plastic plateau is a neutral
                    # branch and tells nothing about path following)
EPSH   = 0.05       # end of the tabulated hardening range
H      = 1.0        # cube edge
NSLICE = 20         # number of slices along the bar
IWEAK  = 10         # 0-based index of the damaging slice
UEND   = 0.055      # prescribed end displacement at lambda = 1


def kuhn(c):
    """Standard 6-tetrahedron decomposition of a hexahedron.

    c is the list of the 8 corner node ids in the CalculiX C3D8 order
    (bottom face 0..3 counter-clockwise, then the matching top face 4..7).
    All six tetrahedra share the c0-c6 diagonal and each has volume V/6.
    """
    v0, v1, v2, v3, v4, v5, v6, v7 = c
    return [(v0, v1, v2, v6), (v0, v2, v3, v6), (v0, v3, v7, v6),
            (v0, v7, v4, v6), (v0, v4, v5, v6), (v0, v5, v1, v6)]


def signed_vol6(p, q, r, s):
    ax, ay, az = q[0] - p[0], q[1] - p[1], q[2] - p[2]
    bx, by, bz = r[0] - p[0], r[1] - p[1], r[2] - p[2]
    cx, cy, cz = s[0] - p[0], s[1] - p[1], s[2] - p[2]
    return (ax * (by * cz - bz * cy)
            - ay * (bx * cz - bz * cx)
            + az * (bx * cy - by * cx))


def build(args):
    nid = {}
    coord = {}
    n = 0
    for i in range(args.nslice + 1):
        for j in (0, 1):
            for k in (0, 1):
                n += 1
                nid[(i, j, k)] = n
                coord[n] = (i * args.h, j * args.h, k * args.h)

    bulk, weak = [], []
    eid = 0
    for i in range(args.nslice):
        c = [nid[(i, 0, 0)], nid[(i + 1, 0, 0)],
             nid[(i + 1, 1, 0)], nid[(i, 1, 0)],
             nid[(i, 0, 1)], nid[(i + 1, 0, 1)],
             nid[(i + 1, 1, 1)], nid[(i, 1, 1)]]
        for t in kuhn(c):
            eid += 1
            t = list(t)
            if signed_vol6(*[coord[x] for x in t]) < 0.0:
                t[1], t[2] = t[2], t[1]
            (weak if i == args.iweak else bulk).append((eid, t))

    left = [nid[(0, j, k)] for j in (0, 1) for k in (0, 1)]
    right = [nid[(args.nslice, j, k)] for j in (0, 1) for k in (0, 1)]

    L = []
    a = L.append
    a("** Analytical snap-back benchmark for dissipation path following.")
    a("** Generated by mkbar.py - do not edit by hand.")
    a("**")
    a("**   E=%g nu=%g sigma_y=%g eps0=%g u_f=%g h=%g slices=%d weak=%d"
      % (args.e, args.nu, args.sigy, args.eps0, args.uf, args.h,
         args.nslice, args.iweak))
    Le = (args.nslice - 1) * args.h
    sig0 = args.sigy + args.hmod * args.eps0     # stress at damage initiation
    Ltot = args.nslice * args.h
    upk = sig0 * Ltot / args.e + Ltot * args.eps0  # uniform flow up to eps0
    ufl = upk - sig0 * Ltot / args.e + args.uf     # elastic release of the bar
    r = Le * sig0 / (args.e * args.uf)
    a("**   snap-back indicator r = L_e*sigma_y/(E*u_f) = %.4f  (>1 required)"
      % r)
    a("**   exact peak      end displacement u = %.6f at sigma = %g" % (upk, args.sigy))
    a("**   exact failure   end displacement u = %.6f at sigma = 0" % ufl)
    a("**   the end displacement must DECREASE by %.6f while the bar fails"
      % (upk - ufl))
    a("*NODE, NSET=NALL")
    for i in range(1, n + 1):
        x, y, z = coord[i]
        a("%d, %.10g, %.10g, %.10g" % (i, x, y, z))

    a("*ELEMENT, TYPE=C3D4, ELSET=EBULK")
    for e, t in bulk:
        a("%d, %d, %d, %d, %d" % (e, t[0], t[1], t[2], t[3]))
    a("*ELEMENT, TYPE=C3D4, ELSET=EWEAK")
    for e, t in weak:
        a("%d, %d, %d, %d, %d" % (e, t[0], t[1], t[2], t[3]))

    a("*NSET, NSET=NLEFT")
    a(", ".join(str(x) for x in left))
    a("*NSET, NSET=NRIGHT")
    a(", ".join(str(x) for x in right))
    a("*NSET, NSET=NPIN")
    a("%d" % nid[(0, 0, 0)])
    a("*NSET, NSET=NROLL")
    a("%d" % nid[(0, 1, 0)])

    sigh = args.sigy + args.hmod * args.epsh
    a("** Both materials share one elastic-plastic law, so the WHOLE bar")
    a("** yields together and contracts laterally without mismatch.  That")
    a("** keeps the stress field uniformly uniaxial, which linear tetrahedra")
    a("** represent exactly, so the 1-D reference above is the FE answer and")
    a("** not just an estimate.  A single plastic slice between elastic")
    a("** blocks would instead be laterally constrained and would carry an")
    a("** axial stress well above sigma_y.")
    a("** The slices differ ONLY in that the weak one can damage.")
    a("*MATERIAL, NAME=BULK")
    a("*ELASTIC")
    a("%g, %g" % (args.e, args.nu))
    a("*PLASTIC")
    a("%g, 0." % args.sigy)
    a("%g, %g" % (sigh, args.epsh))
    a("*MATERIAL, NAME=WEAK")
    a("*ELASTIC")
    a("%g, %g" % (args.e, args.nu))
    a("*PLASTIC")
    a("%g, 0." % args.sigy)
    a("%g, %g" % (sigh, args.epsh))
    a("*DAMAGE INITIATION, CRITERION=RICETRACEY, EVOLUTION=DISPLACEMENT")
    a("%g, 1.0, %g, 0." % (args.eps0, args.uf))
    a("*SOLID SECTION, ELSET=EBULK, MATERIAL=BULK")
    a("*SOLID SECTION, ELSET=EWEAK, MATERIAL=WEAK")

    a("*BOUNDARY")
    a("NLEFT, 1, 1, 0.")
    a("NPIN, 2, 3, 0.")
    a("NROLL, 3, 3, 0.")

    a("*STEP, INC=%d, NLGEOM" % args.inc)
    a("*STATIC")
    a("%g, 1.0, %g, %g" % (args.dt0, args.dtmin, args.dtmax))
    a("*BOUNDARY")
    a("NRIGHT, 1, 1, %g" % args.uend)
    a("*NODE PRINT, NSET=NRIGHT, TOTALS=ONLY")
    a("RF")
    a("*NODE PRINT, NSET=NRIGHT")
    a("U")
    a("*EL PRINT, ELSET=EWEAK")
    a("SDV")
    a("*NODE FILE")
    a("U, RF")
    a("*EL FILE")
    a("S, PEEQ")
    a("*END STEP")
    return "\n".join(L) + "\n"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-o", "--out", default="snapback.inp")
    p.add_argument("--e", type=float, default=E)
    p.add_argument("--nu", type=float, default=NU)
    p.add_argument("--sigy", type=float, default=SIGY)
    p.add_argument("--eps0", type=float, default=EPS0)
    p.add_argument("--uf", type=float, default=UF)
    p.add_argument("--hmod", type=float, default=HMOD)
    p.add_argument("--epsh", type=float, default=EPSH)
    p.add_argument("--h", type=float, default=H)
    p.add_argument("--nslice", type=int, default=NSLICE)
    p.add_argument("--iweak", type=int, default=IWEAK)
    p.add_argument("--uend", type=float, default=UEND)
    p.add_argument("--inc", type=int, default=2000)
    p.add_argument("--dt0", type=float, default=0.01)
    p.add_argument("--dtmin", type=float, default=1e-9)
    p.add_argument("--dtmax", type=float, default=0.02)
    a = p.parse_args()
    open(a.out, "w").write(build(a))
    print("wrote", a.out)


if __name__ == "__main__":
    main()
