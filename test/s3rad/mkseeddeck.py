#!/usr/bin/env python3
"""Generate an s<N>rad deck: the s3rad specimen with the hydride plates
placed by random seed N.

What is kept from the committed deck, and how
---------------------------------------------
Everything except the mesh is copied VERBATIM from
m12_s3rad_gc24_w.inp (hash checked): the materials, the two cohesive
user sections, the step, the boundary conditions and the output cards.
The mesh is rebuilt to the geometry and the grading measured on that
deck (tools used: deckgeom.py / deckiface.py in the session log, forum
entry of 2026-09-25):

  block            4 x 3.6 x 1.6, load along x (FACE_X0 / FACE_XL)
  plates           12 "pills": radius 0.3, thickness 0.09, the rim a
                   half circle of radius 0.045; axis along x, i.e. the
                   plate normal is the load direction
  placement        centres in x 1.55..2.45, y 0.4..3.2, z 0.4..1.2;
                   clearance >= 0.1 between plates (s3rad: min 0.104)
  interface        every plate surface carries UC6 facets; nodes 1-3 on
                   the matrix side, 4-6 duplicated nodes on the hydride
                   side, (x2-x1)x(x3-x1) pointing into the hydride
                   (cohesive_uc6.f: positive normal separation = opening)
  seed facet       INTERFACE_SEED = the facet nearest the centre of the
                   -x face of the plate with the smallest x - this is
                   where s3rad has it (plate 8, r = 0.017)
  grading          size field on the distance to the plates, calibrated
                   against s3rad's element counts and h(distance)

The original generator and its random number stream are not in any of
the repositories, and the s3rad centres do not follow numpy's or
Python's stream for seed 3, so s<N> here is NOT the s<N> the original
generator would have produced.  --centres-from DECK rebuilds the mesh
around another deck's plate centres; run on the committed deck, that
is how this generator is checked against s3rad.

    ./mkseeddeck.py --seed 4 -o s4rad.inp
    ./mkseeddeck.py --centres-from m12_s3rad_gc24_w.inp -o s3regen.inp

Mixed orientation (--fn)
------------------------
--fn P makes P percent of the plates RADIAL (axis along x, the load, as in
every deck above) and the rest TANGENTIAL (axis along y: the plate lies in
the x-z plane, parallel to the load).  round(P/100*12) plates are radial.
The centres are placed so that EVERY pair keeps the 0.1 clearance in
EVERY combination of the two orientations, so one seed gives one set of
centres for the whole series and only the orientation changes between
its members.  Which plates are radial is a fixed random order per seed,
taken from the front, so the radial sets are nested: the plates radial at
P=40 are radial at P=60 too.  Without --fn nothing changes - the decks are
byte-identical to what this script wrote before the option existed.

    for p in 0 20 40 60 80 100; do ./mkseeddeck.py --seed 5 --fn $p -o fn$p.inp; done
"""

import argparse
import hashlib
import os
import sys

import numpy as np

SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   'm12_s3rad_gc24_w.inp')
SRC_SHA = ('2fb0cf4e3554282e1f85cf641c788939bc7a2fa5918dd842f54d38844df'
           'a5391')

LX, LY, LZ = 4.0, 3.6, 1.6
NPLATE = 12
RAD, THK = 0.3, 0.09
RIM = THK / 2.
BOX_X = (1.55, 2.45)
BOX_Y = (0.4, LY - 0.4)
BOX_Z = (0.4, LZ - 0.4)
CLEAR = 0.1


def pill_gap(a, b):
    """Clearance between two coaxial (x) pills with centres a, b."""
    dx = abs(a[0] - b[0])
    rho = np.hypot(a[1] - b[1], a[2] - b[2])
    core = RAD - RIM
    return np.hypot(dx, max(rho - 2. * core, 0.)) - THK


AX_X = np.array([1., 0., 0.])
AX_Y = np.array([0., 1., 0.])


def _basis(a):
    t = np.array([0., 0., 1.]) if abs(a[2]) < 0.9 else np.array([1., 0., 0.])
    u = np.cross(a, t)
    u /= np.linalg.norm(u)
    return u, np.cross(a, u)


def _disk_pts(c, a, n_r=8, n_t=48):
    """Points on the core disc of a pill (centre c, axis a)."""
    u, v = _basis(a)
    core = RAD - RIM
    P = [c]
    th = np.linspace(0., 2. * np.pi, n_t, endpoint=False)
    for r in np.linspace(core / n_r, core, n_r):
        P += list(c + r * (np.outer(np.cos(th), u) + np.outer(np.sin(th), v)))
    return np.array(P)


def _pt_disk(P, c, a):
    """Exact distance from points P to the core disc (c, a)."""
    d = P - c
    h = d @ a
    rho = np.linalg.norm(d - np.outer(h, a), axis=1)
    return np.hypot(h, np.maximum(rho - (RAD - RIM), 0.))


def pill_gap_axes(c1, a1, c2, a2):
    """Clearance between two pills of any axes: the pill is its core disc
    swollen by RIM, so the gap is the disc-to-disc distance less THK.  The
    distance is sampled on one disc (resolution ~0.005) and taken exactly
    to the other, both ways; for coaxial x pills it reproduces pill_gap()
    to three figures on the five seed decks."""
    return min(_pt_disk(_disk_pts(c2, a2), c1, a1).min(),
               _pt_disk(_disk_pts(c1, a1), c2, a2).min()) - THK


def pill_gap_any(a, b):
    """Smallest clearance over all four radial/tangential combinations."""
    if np.linalg.norm(a - b) > 2. * RAD + CLEAR + 1e-9:
        return 1.
    return min(pill_gap_axes(a, p, b, q) for p in (AX_X, AX_Y)
               for q in (AX_X, AX_Y))


def place_any(seed):
    """As place(), but the clearance holds for any orientation of any
    plate (--fn).  A different stream from place(): the two are not
    meant to give the same centres."""
    rng = np.random.default_rng(seed)
    lo = np.array([BOX_X[0], BOX_Y[0], BOX_Z[0]])
    hi = np.array([BOX_X[1], BOX_Y[1], BOX_Z[1]])
    for attempt in range(1000):
        c = []
        tries = 0
        while len(c) < NPLATE and tries < 100000:
            tries += 1
            p = lo + (hi - lo) * rng.random(3)
            if all(pill_gap_any(p, q) >= CLEAR for q in c):
                c.append(p)
        if len(c) == NPLATE:
            return np.round(np.array(c), 4)
    sys.exit('could not place %d plates' % NPLATE)


def axes_for(seed, fn):
    """Radial plates first in a fixed random order: nested across fn."""
    order = np.random.default_rng(10000 + seed).permutation(NPLATE)
    nrad = int(round(fn / 100. * NPLATE))
    ax = [AX_Y.copy() for _ in range(NPLATE)]
    for k in order[:nrad]:
        ax[k] = AX_X.copy()
    return ax, nrad


def place(seed):
    rng = np.random.default_rng(seed)
    lo = np.array([BOX_X[0], BOX_Y[0], BOX_Z[0]])
    hi = np.array([BOX_X[1], BOX_Y[1], BOX_Z[1]])
    for attempt in range(1000):
        c = []
        tries = 0
        while len(c) < NPLATE and tries < 100000:
            tries += 1
            p = lo + (hi - lo) * rng.random(3)
            if all(pill_gap(p, q) >= CLEAR for q in c):
                c.append(p)
        if len(c) == NPLATE:
            return np.round(np.array(c), 4)
    sys.exit('could not place %d plates' % NPLATE)


def centres_from(deck):
    co = {}
    mode = False
    for l in open(deck):
        s = l.strip()
        if s.upper().startswith('*NODE'):
            mode = True
            continue
        if s.startswith('*'):
            if mode:
                break
            continue
        if mode and s:
            t = s.replace(',', ' ').split()
            co[int(t[0])] = np.array(list(map(float, t[1:4])))
    return np.array([(co[2 * k + 1] + co[2 * k + 2]) / 2.
                     for k in range(NPLATE)])


def mesh(centres, smin, smax, dmax, sin, algo3d, mseed, verbose,
         axes=None):
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber('General.Terminal', 1 if verbose else 0)
    gmsh.model.add('srad')
    occ = gmsh.model.occ
    box = occ.addBox(0, 0, 0, LX, LY, LZ)
    pills = []
    core = RAD - RIM
    for k, c in enumerate(centres):
        if axes is None:
            cyl = occ.addCylinder(c[0] - RIM, c[1], c[2], THK, 0, 0, core)
            tor = occ.addTorus(c[0], c[1], c[2], core, RIM, zAxis=[1, 0, 0])
        else:
            a = axes[k]
            b = c - RIM * a
            cyl = occ.addCylinder(b[0], b[1], b[2], THK * a[0], THK * a[1],
                                  THK * a[2], core)
            tor = occ.addTorus(c[0], c[1], c[2], core, RIM,
                               zAxis=list(a))
        out, _ = occ.fuse([(3, cyl)], [(3, tor)])
        pills.append(out[0][1])
    out, omap = occ.fragment([(3, box)], [(3, p) for p in pills])
    occ.synchronize()
    # omap[0] = pieces of the box, omap[1+k] = pieces of pill k
    ptag = []
    for k in range(NPLATE):
        vs = [t for d, t in omap[1 + k]]
        if len(vs) != 1:
            sys.exit('plate %d fragmented into %d volumes' % (k, len(vs)))
        ptag.append(vs[0])
    mtag = [t for d, t in omap[0] if t not in ptag]
    if len(mtag) != 1:
        sys.exit('matrix is %d volumes' % len(mtag))
    psurf = []
    for v in ptag:
        psurf.append([t for d, t in gmsh.model.getBoundary(
            [(3, v)], oriented=False)])
    allsurf = sorted(set(s for l in psurf for s in l))

    f = gmsh.model.mesh.field
    fd = f.add('Distance')
    f.setNumbers(fd, 'SurfacesList', allsurf)
    f.setNumber(fd, 'Sampling', 60)
    ft = f.add('Threshold')
    f.setNumber(ft, 'InField', fd)
    f.setNumber(ft, 'SizeMin', smin)
    f.setNumber(ft, 'SizeMax', smax)
    f.setNumber(ft, 'DistMin', 0.)
    f.setNumber(ft, 'DistMax', dmax)
    # the hydride is meshed finer than its own surface: s3rad's plate
    # tetrahedra have a median longest edge of 0.074 against a facet edge
    # of 0.056, i.e. more than one layer through the 0.09 thickness
    fc = f.add('Constant')
    f.setNumbers(fc, 'VolumesList', ptag)
    f.setNumber(fc, 'VIn', sin)
    f.setNumber(fc, 'VOut', 1.e22)
    f.setNumber(fc, 'IncludeBoundary', 0)
    fm = f.add('Min')
    f.setNumbers(fm, 'FieldsList', [ft, fc])
    f.setAsBackgroundMesh(fm)
    gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary', 0)
    gmsh.option.setNumber('Mesh.MeshSizeFromPoints', 0)
    gmsh.option.setNumber('Mesh.MeshSizeFromCurvature', 0)
    gmsh.option.setNumber('Mesh.Algorithm', 6)
    gmsh.option.setNumber('Mesh.Algorithm3D', algo3d)
    gmsh.option.setNumber('Mesh.Optimize', 1)
    gmsh.option.setNumber('Mesh.OptimizeNetgen', 1)
    gmsh.option.setNumber('Mesh.ElementOrder', 1)
    if mseed:
        gmsh.option.setNumber('Mesh.RandomSeed', mseed)
    gmsh.model.mesh.generate(3)

    tags, xyz, _ = gmsh.model.mesh.getNodes()
    xyz = xyz.reshape(-1, 3)
    X = {int(t): xyz[i] for i, t in enumerate(tags)}

    def tets(v):
        ty, et, en = gmsh.model.mesh.getElements(3, v)
        for t, n in zip(ty, en):
            if t == 4:
                return n.reshape(-1, 4).astype(int)
        return np.zeros((0, 4), int)

    def tris(s):
        ty, et, en = gmsh.model.mesh.getElements(2, s)
        for t, n in zip(ty, en):
            if t == 2:
                return n.reshape(-1, 3).astype(int)
        return np.zeros((0, 3), int)

    mat = tets(mtag[0])
    hyd = [tets(v) for v in ptag]
    ptri = [np.vstack([tris(s) for s in l]) for l in psurf]
    gmsh.finalize()
    return X, mat, hyd, ptri


def orient_tet(t, X):
    a, b, c, d = (X[i] for i in t)
    if np.dot(np.cross(b - a, c - a), d - a) < 0.:
        return [t[0], t[2], t[1], t[3]]
    return list(t)


def build(centres, X, mat, hyd, ptri, axes=None):
    # renumber: nodes 1..N in gmsh order, duplicates appended
    old = sorted(X)
    nid = {o: i + 1 for i, o in enumerate(old)}
    co = {nid[o]: X[o] for o in old}
    nnext = len(old) + 1
    elems = []          # (id, type, nodes)
    eid = 1
    matrix, plate, pplate = [], [], []
    for t in mat:
        elems.append((eid, 'C3D4', orient_tet([nid[i] for i in t], co)))
        matrix.append(eid)
        eid += 1
    facets = []
    seed_pick = int(np.argmin(centres[:, 0]))
    # the pill's -x extreme: the flat face of a radial plate (RIM from the
    # centre), the rim of a tangential one (RAD)
    seed_ext = RIM if (axes is None or axes[seed_pick][0] > 0.5) else RAD
    seed_target = centres[seed_pick] - np.array([seed_ext, 0., 0.])
    seed_best = (1e9, None)
    for k in range(NPLATE):
        surf = set(nid[i] for i in ptri[k].ravel())
        dup = {}
        for n in sorted(surf):
            dup[n] = nnext
            co[nnext] = co[n].copy()
            nnext += 1
        pc = centres[k]
        for t in hyd[k]:
            nodes = [nid[i] for i in t]
            nodes = [dup.get(n, n) for n in nodes]
            elems.append((eid, 'C3D4', orient_tet(nodes, co)))
            plate.append(eid)
            pplate.append(k)
            eid += 1
        for tri in ptri[k]:
            m = [nid[i] for i in tri]
            a, b, c = (co[i] for i in m)
            nrm = np.cross(b - a, c - a)
            cen = (a + b + c) / 3.
            # inward direction of the pill at cen: towards the core disc
            d = cen - pc
            ax = AX_X if axes is None else axes[k]
            h = d @ ax
            radv = d - h * ax
            rho = np.linalg.norm(radv)
            core = RAD - RIM
            if rho > core:
                q = pc + radv * core / rho
            else:
                q = pc + radv
            inward = q - cen
            if np.linalg.norm(inward) < 1e-12:
                inward = -np.sign(h) * ax
            if np.dot(nrm, inward) < 0.:
                m = [m[0], m[2], m[1]]
            facets.append([eid, m + [dup[n] for n in m], k, cen])
            if k == seed_pick:
                dist = np.linalg.norm(cen - seed_target)
                if dist < seed_best[0]:
                    seed_best = (dist, eid)
            eid += 1
    return co, elems, matrix, plate, facets, seed_best[1], pplate


def write(out, centres, co, elems, matrix, plate, facets, seed, src_tail,
          header, axes=None, pplate=None):
    def rows(ids, per=16):
        return '\n'.join(', '.join(str(i) for i in ids[j:j + per])
                         for j in range(0, len(ids), per))
    L = []
    L.append('** C3D4 two-phase disk with duplicated ZrH interface nodes '
             'and UC6 facets')
    for h in header:
        L.append('** ' + h)
    L.append('** Interface facets: %d' % len(facets))
    L.append('** Plate centres (x y z):')
    for k, c in enumerate(centres):
        if axes is None:
            L.append('**   %2d  %.4f %.4f %.4f' % (k, c[0], c[1], c[2]))
        else:
            L.append('**   %2d  %.4f %.4f %.4f  %s' % (
                k, c[0], c[1], c[2],
                'radial (axis x)' if axes[k][0] > 0.5 else
                'tangential (axis y)'))
    L.append('*Node')
    for n in sorted(co):
        x = co[n]
        L.append('%d, %.12E, %.12E, %.12E' % (n, x[0], x[1], x[2]))
    L.append('*User Element, Type=UC6, Nodes=6, Integration Points=3, '
             'MaxDof=3')
    L.append('*Element, Type=C3D4')
    for e, ty, nd in elems:
        L.append('%d, %s' % (e, ', '.join(map(str, nd))))
    L.append('*Element, Type=UC6')
    for f in facets:
        L.append('%d, %s' % (f[0], ', '.join(map(str, f[1]))))
    L.append('*Elset, Elset=MATRIX')
    L.append(rows(matrix))
    L.append('*Elset, Elset=PLATETANGENTIAL')
    L.append(rows(plate))
    if axes is not None:
        # the same elements as PLATETANGENTIAL (which carries the ZrH
        # section whatever the orientation), split for post-processing
        rad = [e for e, k in zip(plate, pplate) if axes[k][0] > 0.5]
        tan = [e for e, k in zip(plate, pplate) if axes[k][0] <= 0.5]
        if rad:
            L.append('*Elset, Elset=HYDRIDE_RADIAL')
            L.append(rows(rad))
        if tan:
            L.append('*Elset, Elset=HYDRIDE_TANGENTIAL')
            L.append(rows(tan))
    L.append('*Elset, Elset=INTERFACE_SEED')
    L.append(str(seed))
    L.append('*Elset, Elset=INTERFACE_REGULAR')
    L.append(rows([f[0] for f in facets if f[0] != seed]))
    L.append('*Elset, Elset=INTERFACE')
    L.append(rows([f[0] for f in facets]))
    # grip faces and the two points holding rigid motion: only nodes of
    # the matrix mesh (duplicates are never on the box faces)
    x0 = [n for n in sorted(co) if abs(co[n][0]) < 1e-9]
    xl = [n for n in sorted(co) if abs(co[n][0] - LX) < 1e-9]

    def at(p):
        n = min(co, key=lambda i: np.linalg.norm(co[i] - np.array(p)))
        if np.linalg.norm(co[n] - np.array(p)) > 1e-9:
            sys.exit('no node at %s' % (p,))
        return n
    L.append('*Nset, Nset=FACE_X0_NSET')
    L.append(rows(x0))
    L.append('*Nset, Nset=FACE_XL_NSET')
    L.append(rows(xl))
    L.append('*Nset, Nset=FIXPOINTA')
    L.append(str(at((0., 0., 0.))))
    L.append('*Nset, Nset=FIXPOINTB')
    L.append(str(at((0., LY, 0.))))
    open(out, 'w').write('\n'.join(L) + '\n' + src_tail)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--seed', type=int)
    ap.add_argument('--centres-from')
    ap.add_argument('-o', '--out', required=True)
    ap.add_argument('--smin', type=float, default=0.0592)
    ap.add_argument('--smax', type=float, default=0.24)
    ap.add_argument('--dmax', type=float, default=0.40)
    ap.add_argument('--sin', type=float, default=0.033)
    ap.add_argument('--algo3d', type=int, default=1)
    # same geometry, a different mesh: measures how much of a difference
    # between seeds a mesh alone can make (the deck is local)
    ap.add_argument('--mesh-seed', type=int, default=0)
    # Turon et al. 2007 eq. 7: K = alpha*E/t with alpha >= 50 keeps the
    # interface's added compliance under 2 percent.  The committed deck's
    # Kn = 2e5 N/mm^3 is alpha ~ 0.14 for a 0.09 thick ZrH plate (E=1.3e5):
    # plate plus its two interfaces ~ 12 percent of the plate's stiffness.
    # The dead facet keeps gmin*Kn in tension (cohesive_uc6.f, g =
    # max(gmin, 1-Dvis)).  Raising Kn 500x raises that residual 500x: at
    # Kn=1e8 a 'dead' facet opened 0.05 mm still carries ~50 MPa.  --gmin
    # rescales it; gmin*Kn = 2 N/mm^3 is the committed deck's value.
    ap.add_argument('--gmin', type=float, default=None,
                    help='residual stiffness fraction gmin for both user '
                         'sections; default keeps the deck value')
    ap.add_argument('--kn', type=float, default=None,
                    help='cohesive normal stiffness Kn [N/mm^3] for both '
                         'user sections; default keeps the deck value')
    ap.add_argument('--fn', type=float, default=None,
                    help='percent of plates radial (axis x); the rest are '
                         'tangential (axis y).  Needs --seed; see the '
                         'module docstring')
    ap.add_argument('-v', action='store_true')
    a = ap.parse_args()
    if (a.seed is None) == (a.centres_from is None):
        sys.exit('give exactly one of --seed, --centres-from')

    src = open(SRC, 'rb').read()
    if hashlib.sha256(src).hexdigest() != SRC_SHA:
        sys.exit('source deck hash mismatch')
    txt = src.decode('ascii')
    i = txt.index('*Material, Name=ZR')
    src_tail = txt[i:]
    if not src_tail.endswith('\n'):
        src_tail += '\n'
    kn_note = 'Kn as in the committed deck'
    if a.kn is not None or a.gmin is not None:
        # --kn replaces the FIRST constant (Kn, normal penalty stiffness,
        # cohesive_uc6.f prop 1) of every *User Section data line and
        # nothing else.  delta_f = 2 Gc / Tn0 does not depend on it, so
        # the fracture energy and the softening slope are unchanged; only
        # delta_0 = Tn0 / Kn and the interface compliance move.
        lines = src_tail.split('\n')
        nrep = 0
        for j, l in enumerate(lines):
            if l.strip().upper().startswith('*USER SECTION'):
                f = lines[j + 1].split(',')
                if a.kn is not None:
                    f[0] = '%.6e' % a.kn
                if a.gmin is not None:
                    f[4] = ' %.6e' % a.gmin
                lines[j + 1] = ','.join(f)
                nrep += 1
        if nrep != 2:
            sys.exit('expected two *User Section cards, found %d' % nrep)
        src_tail = '\n'.join(lines)
        kn_note = ('both *User Section cards: Kn %s, gmin %s'
                   % ('%.6e' % a.kn if a.kn is not None else 'as deck',
                      '%.6e' % a.gmin if a.gmin is not None else 'as deck'))

    axes = None
    if a.fn is not None:
        if a.seed is None:
            sys.exit('--fn needs --seed')
        if not 0. <= a.fn <= 100.:
            sys.exit('--fn is a percentage, 0..100')
        centres = place_any(a.seed)
        axes, nrad = axes_for(a.seed, a.fn)
        header = ['Generated by mkseeddeck.py --seed %d --fn %g (numpy '
                  'default_rng; centres valid for any orientation)'
                  % (a.seed, a.fn),
                  'Fn: %d of %d plates radial (axis x, the load), %d '
                  'tangential (axis y) = %.1f percent radial'
                  % (nrad, NPLATE, NPLATE - nrad, 100. * nrad / NPLATE)]
    elif a.seed is not None:
        centres = place(a.seed)
        header = ['Generated by mkseeddeck.py --seed %d (numpy default_rng)'
                  % a.seed]
    else:
        centres = centres_from(a.centres_from)
        header = ['Generated by mkseeddeck.py --centres-from %s'
                  % os.path.basename(a.centres_from)]
    header.append('Mesh: smin %.4g smax %.4g dmax %.4g sin %.4g mesh-seed %d; all cards after '
                  'the sets copied from m12_s3rad_gc24_w.inp sha256 %s'
                  % (a.smin, a.smax, a.dmax, a.sin, a.mesh_seed,
                     SRC_SHA[:16]))
    header.append(kn_note)
    X, mat, hyd, ptri = mesh(centres, a.smin, a.smax, a.dmax, a.sin,
                             a.algo3d, a.mesh_seed, a.v, axes)
    co, elems, matrix, plate, facets, seed, pplate = build(
        centres, X, mat, hyd, ptri, axes)
    write(a.out, centres, co, elems, matrix, plate, facets, seed, src_tail,
          header, axes, pplate)
    print('%s: nodes %d, MATRIX %d, PLATETANGENTIAL %d, UC6 %d, seed facet '
          '%d' % (a.out, len(co), len(matrix), len(plate), len(facets),
                  seed))


if __name__ == '__main__':
    main()
