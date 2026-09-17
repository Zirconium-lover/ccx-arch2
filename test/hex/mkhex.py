#!/usr/bin/env python3
"""A C3D8 bar with DE1 ductile damage.

Exists to answer one question that no deck in this tree can ask: what
happens to CB1's projected crack-band width on an element that is not a
linear tetrahedron.  Before W1 every DE1 path refused anything but C3D4,
so there was nothing to test; CCX_DAMAGE_CHARLEN=1 lifts that, and the
volume estimate that used to be the fallback is filled for C3D4 only.

--mode uniaxial    pull along x.  Principal stresses are well separated,
                   so damcbevec always finds a direction.
--mode hydro       pull all three faces equally.  Principal stresses are
                   equal by symmetry, which is exactly the state
                   damcbwidth refuses with iok=0.
"""
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("-o", required=True)
ap.add_argument("--mode", choices=("uniaxial", "hydro"), default="uniaxial")
ap.add_argument("--nx", type=int, default=4)
ap.add_argument("--ny", type=int, default=2)
ap.add_argument("--nz", type=int, default=2)
# NONLOCAL= on the material card rather than CCX_DAMAGE_NONLOCAL in the
# environment.  The two routes went to different places until they were
# joined: every guard asked ellsave, which only the environment writes, so
# a card-only deck never entered the nonlocal block and ran unregularised
# while looking healthy.  A deck that can ask for the length on the CARD is
# what lets a test hold that seam shut.
ap.add_argument("--nonlocal-ell", type=float, default=0.0,
                help="put NONLOCAL=<ell> on the *DAMAGE INITIATION cards")
a = ap.parse_args()

nx, ny, nz = a.nx, a.ny, a.nz
hx, hy, hz = 1.0, 1.0, 1.0            # equal edges: a cube element


def nid(i, j, k):
    return 1 + i + (nx + 1) * (j + (ny + 1) * k)


L = []
L.append("*Node")
for k in range(nz + 1):
    for j in range(ny + 1):
        for i in range(nx + 1):
            L.append("%d, %.6f, %.6f, %.6f" % (nid(i, j, k), i*hx, j*hy, k*hz))

L.append("*Element, Type=C3D8, Elset=BAR")
e = 0
for k in range(nz):
    for j in range(ny):
        for i in range(nx):
            e += 1
            n = [nid(i, j, k), nid(i+1, j, k), nid(i+1, j+1, k), nid(i, j+1, k),
                 nid(i, j, k+1), nid(i+1, j, k+1), nid(i+1, j+1, k+1),
                 nid(i, j+1, k+1)]
            L.append("%d, %s" % (e, ", ".join(str(x) for x in n)))

x0 = [nid(0, j, k) for k in range(nz+1) for j in range(ny+1)]
xL = [nid(nx, j, k) for k in range(nz+1) for j in range(ny+1)]
y0 = [nid(i, 0, k) for k in range(nz+1) for i in range(nx+1)]
yL = [nid(i, ny, k) for k in range(nz+1) for i in range(nx+1)]
z0 = [nid(i, j, 0) for j in range(ny+1) for i in range(nx+1)]
zL = [nid(i, j, nz) for j in range(ny+1) for i in range(nx+1)]

for name, s in (("X0", x0), ("XL", xL), ("Y0", y0), ("YL", yL),
                ("Z0", z0), ("ZL", zL)):
    L.append("*Nset, Nset=%s" % name)
    for i in range(0, len(s), 8):
        L.append(", ".join(str(x) for x in s[i:i+8]))

# One weak element in the middle, so damage has somewhere to start and the
# run does not depend on round-off to pick a site.
weak = 1 + (nx // 2) + nx * ((ny // 2) + ny * (nz // 2))
L.append("*Elset, Elset=WEAK")
L.append("%d" % weak)
L.append("*Elset, Elset=STRONG, Generate")
L.append("1, %d, 1" % e)

DMG = "*Damage Initiation, Criterion=Ductile, Evolution=Displacement, Npoints=6"
if a.nonlocal_ell > 0.0:
    DMG += ", NONLOCAL=%g" % a.nonlocal_ell

L += [
    "*Material, Name=STEEL",
    "*Elastic",
    "210000., 0.3",
    "*Plastic",
    "300., 0.", "330., 0.010", "360., 0.050", "380., 0.150",
    DMG,
    "1.0, 0.0150, 0.00, 0.3000, 0.33, 0.1200, 0.50, 0.0700",
    "0.67, 0.0425, 1.00, 0.0250, 1.50, 0.0150",
    "*Material, Name=WEAKMAT",
    "*Elastic",
    "210000., 0.3",
    "*Plastic",
    "240., 0.", "260., 0.010", "280., 0.050", "300., 0.150",
    DMG,
    "1.0, 0.0150, 0.00, 0.1500, 0.33, 0.0600, 0.50, 0.0350",
    "0.67, 0.0210, 1.00, 0.0120, 1.50, 0.0070",
    "*Solid Section, Elset=STRONG, Material=STEEL",
    "*Solid Section, Elset=WEAK, Material=WEAKMAT",
    "*Step, Nlgeom, Inc=100000",
    "*Static",
    "0.01, 1.0, 1.e-9, 0.02",
    "*Boundary",
    "X0, 1, 1", "Y0, 2, 2", "Z0, 3, 3",
]
if a.mode == "uniaxial":
    L += ["*Boundary", "XL, 1, 1, %.4f" % (0.30 * nx * hx)]
else:
    L += ["*Boundary",
          "XL, 1, 1, %.4f" % (0.30 * nx * hx),
          "YL, 2, 2, %.4f" % (0.30 * ny * hy),
          "ZL, 3, 3, %.4f" % (0.30 * nz * hz)]
L += [
    "*Node File", "U",
    "*El File", "S, E, PEEQ",
    "*End Step",
]
open(a.o, "w").write("\n".join(L) + "\n")
print("wrote %s: %d nodes, %d C3D8, mode=%s, weak element %d"
      % (a.o, (nx+1)*(ny+1)*(nz+1), e, a.mode, weak))
