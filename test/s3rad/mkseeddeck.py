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


def mesh(centres, smin, smax, dmax, sin, algo3d, mseed, verbose):
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber('General.Terminal', 1 if verbose else 0)
    gmsh.model.add('srad')
    occ = gmsh.model.occ
    box = occ.addBox(0, 0, 0, LX, LY, LZ)
    pills = []
    core = RAD - RIM
    for c in centres:
        cyl = occ.addCylinder(c[0] - RIM, c[1], c[2], THK, 0, 0, core)
        tor = occ.addTorus(c[0], c[1], c[2], core, RIM, zAxis=[1, 0, 0])
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


def build(centres, X, mat, hyd, ptri):
    # renumber: nodes 1..N in gmsh order, duplicates appended
    old = sorted(X)
    nid = {o: i + 1 for i, o in enumerate(old)}
    co = {nid[o]: X[o] for o in old}
    nnext = len(old) + 1
    elems = []          # (id, type, nodes)
    eid = 1
    matrix, plate = [], []
    for t in mat:
        elems.append((eid, 'C3D4', orient_tet([nid[i] for i in t], co)))
        matrix.append(eid)
        eid += 1
    facets = []
    seed_pick = int(np.argmin(centres[:, 0]))
    seed_target = centres[seed_pick] - np.array([RIM, 0., 0.])
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
            eid += 1
        for tri in ptri[k]:
            m = [nid[i] for i in tri]
            a, b, c = (co[i] for i in m)
            nrm = np.cross(b - a, c - a)
            cen = (a + b + c) / 3.
            # inward direction of the pill at cen: towards the core disc
            d = cen - pc
            rho = np.hypot(d[1], d[2])
            core = RAD - RIM
            if rho > core:
                q = pc + np.array([0., d[1], d[2]]) * core / rho
            else:
                q = pc + np.array([0., d[1], d[2]])
            inward = q - cen
            if np.linalg.norm(inward) < 1e-12:
                inward = -np.array([np.sign(d[0]), 0., 0.])
            if np.dot(nrm, inward) < 0.:
                m = [m[0], m[2], m[1]]
            facets.append([eid, m + [dup[n] for n in m], k, cen])
            if k == seed_pick:
                dist = np.linalg.norm(cen - seed_target)
                if dist < seed_best[0]:
                    seed_best = (dist, eid)
            eid += 1
    return co, elems, matrix, plate, facets, seed_best[1]


def write(out, centres, co, elems, matrix, plate, facets, seed, src_tail,
          header):
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
        L.append('**   %2d  %.4f %.4f %.4f' % (k, c[0], c[1], c[2]))
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

    if a.seed is not None:
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
    X, mat, hyd, ptri = mesh(centres, a.smin, a.smax, a.dmax, a.sin,
                             a.algo3d, a.mesh_seed, a.v)
    co, elems, matrix, plate, facets, seed = build(centres, X, mat, hyd,
                                                   ptri)
    write(a.out, centres, co, elems, matrix, plate, facets, seed, src_tail,
          header)
    print('%s: nodes %d, MATRIX %d, PLATETANGENTIAL %d, UC6 %d, seed facet '
          '%d' % (a.out, len(co), len(matrix), len(plate), len(facets),
                  seed))


if __name__ == '__main__':
    main()
