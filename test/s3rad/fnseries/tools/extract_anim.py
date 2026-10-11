#!/usr/bin/env python3
"""Pack one finished run into a compact .npz for animations and plots.

    extract_anim.py <run_dir> <out.npz> [fn_nominal]

Reads m.inp (mesh, phases, plates), m.frd (nodal U, S, PEEQ every output
frame), m.dat (grip force; cohesive state of the interface facets), m.damage
(exact deletion time of every element), m.de1.vtk (omega and D of the cells
left at the end) and m.de1stats (damage counts per increment).  The raw
files are 0.7 GB per run; the .npz keeps what a picture needs, in float16
where 3 digits are enough.

Keys (E tets, F interface facets, N nodes, K field frames, M increments):
  meta            json: run, fn_nominal, k_radial, A0, L0, grip, units
  nodes (N,3)     f4  initial coordinates, mm (x = load / T, y = R, z = L)
  tets (E,4)      i4  0-based node indices        tet_id (E,)  element numbers
  tet_phase (E,)  i1  0 matrix, 1 radial hydride, 2 tangential hydride
  tet_plate (E,)  i1  plate index 0..11, -1 for the matrix
  tet_tdel (E,)   f4  time the element was deleted (inf: never)
  tet_D_end, tet_omega_end (E,)  f2  at the last state (nan: deleted)
  plate_c (12,3) f4, plate_axis (12,) i1  0: normal along x (radial), 1: along y
  coh (F,6)       i4  0-based nodes, 3 on each face   coh_id, coh_seed, coh_plate
  coh_d0, coh_df (F,) f4  separation at the strength peak / at full debonding, mm
  coh_tdel (F,)   f4  deletion time of the facet (inf: never)
  t (K,)          f8  step time of each field frame (eps = grip * t / L0)
  U (K,N,3)       f2  displacement, mm
  S (K,N,6)       f2  stress xx yy zz xy xz yz, MPa (nodal average)
  PEEQ (K,N)      f2  equivalent plastic strain
  t_coh (K2,)     f8  times of the interface prints
  coh_D (K2,F)    f2  cohesive damage SDV3 (stiffness loss; mean of the 3 points;
                      it is ~0.5 already at twice coh_d0, so for pictures use
                      (coh_open - coh_d0) / (coh_df - coh_d0): the drop of traction)
  coh_open (K2,F) f2  largest separation reached, mm (mean of the 3 points)
  curve (K3,2)    f8  eps, sigma (MPa) from the grip force
  de1 (M,11)      f8  de1stats rows: step inc step_time total_time passes
                      active Dgt0.1 Dgt0.5 Dgt0.9 Dfull Dmax ... (see header)
"""
import json, os, re, sys
import numpy as np

A0, L0 = 3.6 * 1.6, 4.0


def fixed(seg, n_head, width, nv):
    """Parse fixed-width numeric columns: n_head chars of label, then nv
    fields of `width` chars, per line of equal length."""
    W = seg.index(b'\n') + 1
    n = len(seg) // W
    if n * W != len(seg):
        raise ValueError('ragged block')
    a = np.frombuffer(seg, np.uint8).reshape(n, W)
    v = np.full((n, nv, width + 1), ord(' '), np.uint8)
    v[:, :, :width] = a[:, n_head:n_head + width * nv].reshape(n, nv, width)
    vals = np.fromstring(v.tobytes().decode(), sep=' ').reshape(n, nv)
    lab = np.array(a[:, 3:n_head].tobytes().split(), dtype=np.int64)
    return lab, vals


def read_inp(path):
    nodes, nid = [], []
    tets, tid, coh, cid = [], [], [], []
    sets, plates = {}, []
    mode, cur = None, None
    grip = None
    for line in open(path, errors='replace'):
        if line.startswith('**'):
            m = re.match(r'\*\*\s+(\d+)\s+([-0-9.]+)\s+([-0-9.]+)\s+([-0-9.]+)\s+(radial|tangential)', line)
            if m:
                plates.append((float(m.group(2)), float(m.group(3)), float(m.group(4)),
                               0 if m.group(5) == 'radial' else 1))
            m = re.match(r'\*\*\s*Fn:\s*(\d+) of 12', line)
            if m:
                sets['__k'] = int(m.group(1))
            continue
        s = line.strip()
        if s.startswith('*'):
            low = s.lower().replace(' ', '')
            mode = None
            if low.startswith('*node') and not low.startswith('*nodefile') and not low.startswith('*nodeprint'):
                mode = 'node'
            elif low.startswith('*element,'):
                mode = 'c3d4' if 'type=c3d4' in low else ('uc6' if 'type=uc6' in low else None)
            elif low.startswith('*elset,') or low.startswith('*nset,'):
                cur = re.search(r'(?:elset|nset)=([^,]+)', low).group(1).upper()
                mode = 'set'; sets[cur] = []
            elif low.startswith('*boundary'):
                mode = 'bnd'
            elif low.startswith('*usersection'):
                cur = re.search(r'elset=([^,]+)', low).group(1).upper()
                mode = 'usec'
            continue
        if not s:
            continue
        v = [x for x in s.split(',') if x.strip()]
        if mode == 'node':
            nid.append(int(v[0])); nodes.append([float(x) for x in v[1:4]])
        elif mode == 'c3d4':
            tid.append(int(v[0])); tets.append([int(x) for x in v[1:5]])
        elif mode == 'uc6':
            cid.append(int(v[0])); coh.append([int(x) for x in v[1:7]])
        elif mode == 'set':
            sets[cur] += [int(x) for x in v]
        elif mode == 'usec':
            sets['__usec_' + cur] = [float(x) for x in v]
            mode = None
        elif mode == 'bnd' and v[0].strip().upper() == 'FACE_XL_NSET':
            grip = float(v[3])
    return (np.array(nid), np.array(nodes), np.array(tid), np.array(tets), np.array(cid),
            np.array(coh), sets, plates, grip)


def read_frd(path, want):
    """Nodal blocks by name: lists of (node ids, values) per output frame.
    Nodes left without elements drop out of later frames."""
    buf = open(path, 'rb').read()
    out = {w: [] for w in want}
    times = {w: [] for w in want}
    for m in re.finditer(rb'\n  100CL[^\n]*\n -4  (\w+) ', buf):
        name = m.group(1).decode()
        if name not in want:
            continue
        hdr = buf[m.start() + 1:buf.index(b'\n', m.start() + 1)]
        t = float(hdr[12:24])
        s = buf.index(b'\n -1', m.end()) + 1
        e = buf.index(b'\n -3', s) + 1
        out[name].append(fixed(buf[s:e], 13, 12, want[name]))
        times[name].append(t)
    return times, out


def read_dat(path):
    """Grip force per print; interface state variables per print."""
    buf = open(path, 'rb').read()
    tf, fx = [], []
    for m in re.finditer(rb'total force \(fx,fy,fz\) for set FACE_XL_NSET and time\s+(\S+)\s*\n\s*\n([^\n]*)', buf):
        tf.append(float(m.group(1))); fx.append(float(m.group(2).split()[0]))
    tc, rows = [], []
    for m in re.finditer(rb'internal state variables \(elem, integ\.pnt\.,values\) for set INTERFACE and time\s+(\S+)\s*\n\s*\n', buf):
        e = buf.find(b'\n\n', m.end())
        e = len(buf) if e < 0 else e
        a = np.fromstring(buf[m.end():e].decode(), sep=' ')
        tc.append(float(m.group(1))); rows.append(a)
    return np.array(tf), np.array(fx), np.array(tc), rows


def read_vtk_end(path):
    lines = open(path).read().split('\n')
    arrs = {}
    i = 0
    while i < len(lines):
        if lines[i].startswith('CELL_DATA'):
            nc = int(lines[i].split()[1])
        if lines[i].startswith('SCALARS'):
            name = lines[i].split()[1]
            vals = np.array(' '.join(lines[i + 2:i + 2 + nc]).split(), dtype=float)
            arrs[name] = vals
            i += 2 + nc
            continue
        i += 1
    return arrs


def main(run, out, fn_nominal=None):
    nid, xyz, tid, tets, cid, coh, sets, plates, grip = read_inp(os.path.join(run, 'm.inp'))
    nmap = np.full(nid.max() + 1, -1, np.int64); nmap[nid] = np.arange(len(nid))
    tets0 = nmap[tets].astype(np.int32)
    coh0 = nmap[coh].astype(np.int32)

    phase = np.zeros(len(tid), np.int8)
    is_r = np.isin(tid, sets.get('HYDRIDE_RADIAL', []))
    is_t = np.isin(tid, sets.get('HYDRIDE_TANGENTIAL', []))
    phase[is_r] = 1; phase[is_t] = 2
    pc = np.array([p[:3] for p in plates], np.float32)
    pax = np.array([p[3] for p in plates], np.int8)
    plate = np.full(len(tid), -1, np.int8)
    hyd = np.flatnonzero(phase > 0)
    if len(pc) and len(hyd):
        from scipy.sparse import coo_matrix
        from scipy.sparse.csgraph import connected_components
        # tets and nodes as one graph: hydride tets of a plate share nodes,
        # the duplicated interface nodes keep the plates apart from the matrix
        nn = len(nid)
        r = np.repeat(np.arange(len(hyd)), 4) + nn
        c = tets0[hyd].ravel()
        g = coo_matrix((np.ones(len(r)), (r, c)), shape=(nn + len(hyd),) * 2)
        _, lab = connected_components(g, directed=False)
        lab = lab[nn:]
        cen = xyz[tets0[hyd]].mean(1)
        for L in np.unique(lab):
            j = lab == L
            d = np.linalg.norm(pc - cen[j].mean(0), axis=1)
            d[pax != phase[hyd[j][0]] - 1] = np.inf
            plate[hyd[j]] = np.argmin(d)
    ccen = xyz[coh0].mean(1)
    coh_plate = np.argmin(np.linalg.norm(ccen[:, None, :] - pc[None], axis=2), 1).astype(np.int8) \
        if len(pc) else np.full(len(cid), -1, np.int8)
    coh_seed = np.isin(cid, sets.get('INTERFACE_SEED', []))
    coh_d0 = np.zeros(len(cid), np.float32); coh_df = np.zeros(len(cid), np.float32)
    for key, mask in (('INTERFACE_SEED', coh_seed), ('INTERFACE_REGULAR', ~coh_seed)):
        kn, tn0, ts0, gc = sets['__usec_' + key][:4]
        coh_d0[mask] = tn0 / kn; coh_df[mask] = 2 * gc / tn0    # cohesive_uc6.f

    # deletion history
    tdel = {}
    for line in open(os.path.join(run, 'm.damage')):
        if line.startswith('#'):
            continue
        v = line.split()
        tdel.setdefault(int(v[0]), float(v[4]))
    tet_tdel = np.array([tdel.get(int(e), np.inf) for e in tid], np.float32)
    coh_tdel = np.array([tdel.get(int(e), np.inf) for e in cid], np.float32)

    # fields
    ft, f = read_frd(os.path.join(run, 'm.frd'), {'DISP': 3, 'STRESS': 6, 'PE': 1})
    K = min(len(f['DISP']), len(f['STRESS']), len(f['PE']))
    U = np.full((K, len(nid), 3), np.nan, np.float16)
    S = np.full((K, len(nid), 6), np.nan, np.float16)
    PE = np.full((K, len(nid)), np.nan, np.float16)
    for k in range(K):
        for arr, (ids, vals) in ((U, f['DISP'][k]), (S, f['STRESS'][k]), (PE, f['PE'][k])):
            j = nmap[ids]
            if np.any(j < 0):
                raise ValueError('frd node not in the deck')
            arr[k, j] = vals if arr.ndim == 3 else vals[:, 0]
    t = np.array(ft['DISP'][:K])

    # grip force and interface
    tf, fx, tc, rows = read_dat(os.path.join(run, 'm.dat'))
    curve = np.c_[grip * tf / L0, fx / A0]
    cpos = np.full(cid.max() + 1, -1, np.int64); cpos[cid] = np.arange(len(cid))
    coh_D = np.zeros((len(tc), len(cid)), np.float16)
    coh_open = np.zeros((len(tc), len(cid)), np.float16)
    for k, a in enumerate(rows):
        ncol = None
        for nc in (2 + 13, 2 + 4):
            if len(a) % nc == 0 and np.all(np.isin(a.reshape(-1, nc)[:, 1], (1, 2, 3))):
                ncol = nc; break
        a = a.reshape(-1, ncol)
        j = cpos[a[:, 0].astype(np.int64)]
        dsum = np.bincount(j, a[:, 4], len(cid)); osum = np.bincount(j, a[:, 2], len(cid))
        cnt = np.bincount(j, None, len(cid))
        with np.errstate(invalid='ignore', divide='ignore'):
            coh_D[k] = np.where(cnt > 0, dsum / cnt, 1.)
            coh_open[k] = np.where(cnt > 0, osum / cnt, np.nan)

    # last state
    tet_D_end = np.full(len(tid), np.nan, np.float16)
    tet_om_end = np.full(len(tid), np.nan, np.float16)
    vp = os.path.join(run, 'm.de1.vtk')
    if os.path.exists(vp):
        v = read_vtk_end(vp)
        tpos = np.full(tid.max() + 1, -1, np.int64); tpos[tid] = np.arange(len(tid))
        j = tpos[v['ELEMENT_ID'].astype(np.int64)]
        tet_D_end[j] = v['DE1_D']; tet_om_end[j] = v['DUCT_IP']
    de1 = np.loadtxt(os.path.join(run, 'm.de1stats'), comments='#', ndmin=2)

    meta = dict(run=os.path.basename(os.path.normpath(run)), fn_nominal=fn_nominal,
                k_radial=sets.get('__k'), n_plates=len(pc), A0=A0, L0=L0, grip=grip,
                box=[float(x) for x in xyz.max(0) - xyz.min(0)],
                units='mm, MPa; time t in [0,1] scales the grip displacement',
                de1_cols='step inc step_time total_time passes active Dgt0.1 Dgt0.5 Dgt0.9 Dfull Dmax max_dD',
                coh_sdv='SDV1 max separation, SDV2 viscous damage, SDV3 damage, SDV4 failed')
    np.savez_compressed(out, meta=json.dumps(meta), nodes=xyz.astype(np.float32), node_id=nid.astype(np.int32),
                        tets=tets0, tet_id=tid.astype(np.int32), tet_phase=phase, tet_plate=plate,
                        tet_tdel=tet_tdel, tet_D_end=tet_D_end, tet_omega_end=tet_om_end,
                        plate_c=pc, plate_axis=pax, coh=coh0, coh_id=cid.astype(np.int32),
                        coh_seed=coh_seed, coh_plate=coh_plate, coh_d0=coh_d0, coh_df=coh_df, coh_tdel=coh_tdel,
                        t=t, U=U, S=S, PEEQ=PE, t_coh=tc, coh_D=coh_D, coh_open=coh_open,
                        curve=curve, de1=de1)
    print('%s: N %d, E %d, F %d, K %d, K2 %d, deleted tets %d, facets %d -> %.1f MB'
          % (meta['run'], len(nid), len(tid), len(cid), K, len(tc), np.isfinite(tet_tdel).sum(),
             np.isfinite(coh_tdel).sum(), os.path.getsize(out) / 1e6))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], float(sys.argv[3]) if len(sys.argv) > 3 else None)
