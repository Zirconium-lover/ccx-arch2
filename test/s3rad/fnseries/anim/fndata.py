"""The finished runs, ready for pictures.

The raw output of a run (0.7 GB of .frd / .dat) is packed by
tools/extract_anim.py into animdata/<run>.npz (9-14 MB); this module loads
it and turns it into what a scene draws:

    from fndata import Run, SERIES
    r = Run(0.4)                         # by Fn (the series run) or by run name
    r.eps_f10, r.peak                    # numbers of the final table
    r.at(0.05)                           # state at the grip strain 5 %:
        .U, .X (deformed nodes), .alive (tets not deleted),
        .peeq / .omega / .D / .mises / .triax / .sxx (nodal),
        .debond (facet traction drop 0..1), .sigma (point on the curve)
    r.slice('z', 0.8, state, 'omega')    # section polygons with nodal values
    r.surface(state.alive)               # faces of the remaining body: box / crack
    r.mechanisms(eps)                    # deleted volume: matrix, hydride; debonded area

Time: the solver's step time t in [0, 1] scales the grip displacement; the
grip strain is eps = grip * t / L0.  Fields exist at the output frames (every
10 increments, ~0.5 % of strain); between frames they are interpolated
linearly.  Deletion times are exact (every increment), so a crack grows
smoothly however finely it is sampled.

omega and D are not written per frame by the solver; they are rebuilt here
from the nodal PEEQ and stress history with the same law as calcdamage.f:
omega = sum d(eps_p) / eps_f(eta) up to 1, then D = l_e * d(eps_p) / u_f.
Check against the solver's own values at the last state: check_omega().
"""
import json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get('FN_ANIMDATA', os.path.join(HERE, '..', 'animdata'))

# the series: Fn -> run (SERIES_FINAL.txt)
SERIES = [(0.0, 'fn0'), (0.20, 'fn20'), (0.25, 'fn25v1e-3'), (0.33, 'fn33v1e-3'),
          (0.40, 'fn40'), (0.60, 'fn60'), (0.80, 'fn80v1e-3'), (1.0, 'fn100')]

# ductile damage of the decks: eps_f(eta) tables and u_f (Zr, ZrH)
EF_ETA = np.array([0.0, 0.33, 0.50, 0.67, 1.00, 1.50])
EF = {0: np.array([0.300, 0.120, 0.070, 0.0425, 0.025, 0.015]),     # matrix
      1: np.array([0.120, 0.070, 0.045, 0.028, 0.015, 0.010])}      # hydride
UF = {0: 0.015, 1: 0.004}

FACES = np.array([[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]])     # face i is opposite node i
EDGES = np.array([[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]])
AXIS = {'x': 0, 'y': 1, 'z': 2}


def mises(S):
    s = S.astype(np.float64)
    p = s[..., :3].mean(-1)
    d = s[..., :3] - p[..., None]
    return np.sqrt(1.5 * ((d ** 2).sum(-1) + 2 * (s[..., 3:] ** 2).sum(-1)))


def triax(S):
    s = S.astype(np.float64)
    q = mises(s)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(q > 1e-6, s[..., :3].mean(-1) / q, 0.)


class State:
    """Everything at one grip strain; built by Run.at()."""
    pass


class Run:
    def __init__(self, key):
        name = dict(SERIES).get(key, key) if not isinstance(key, str) else key
        d = np.load(os.path.join(DATA, name + '.npz'))
        self.meta = json.loads(str(d['meta']))
        self.name, self.fn = name, self.meta['fn_nominal']
        self.k_radial = self.meta['k_radial']
        self.L0, self.A0, self.grip = self.meta['L0'], self.meta['A0'], self.meta['grip']
        g = lambda k, t=None: d[k] if t is None else d[k].astype(t)
        self.X0 = g('nodes', np.float64)
        self.T = g('tets'); self.phase = g('tet_phase'); self.plate = g('tet_plate')
        self.tdel = g('tet_tdel', np.float64)
        self.D_end, self.omega_end = g('tet_D_end', np.float64), g('tet_omega_end', np.float64)
        self.plate_c, self.plate_axis = g('plate_c'), g('plate_axis')
        self.coh = g('coh'); self.coh_plate = g('coh_plate')
        self.coh_d0, self.coh_df = g('coh_d0', np.float64), g('coh_df', np.float64)
        self.coh_tdel = g('coh_tdel', np.float64)
        self.t = g('t'); self.eps = self.grip * self.t / self.L0
        self._U, self._S, self._P = d['U'], d['S'], d['PEEQ']
        self.t_coh = g('t_coh'); self._open = d['coh_open']
        self.curve = g('curve')
        self.de1 = g('de1')
        self.lo, self.hi = self.X0.min(0), self.X0.max(0)
        # nodes on both faces of the interface are duplicated: a geometric id
        # joins them (for surfaces and for "which phase does a node carry")
        q = np.round((self.X0 - self.lo) / 1e-6).astype(np.int64)
        _, self.geo = np.unique(q, axis=0, return_inverse=True)
        self.geo = self.geo.ravel()
        nph = np.zeros(len(self.X0), np.int8)
        nph[self.T[self.phase > 0].ravel()] = 1
        self.node_phase = nph                  # 1: the node belongs to hydride tets
        P = self.X0[self.T]
        self.vol = np.abs(np.einsum('ij,ij->i', P[:, 1] - P[:, 0],
                                    np.cross(P[:, 2] - P[:, 0], P[:, 3] - P[:, 0]))) / 6
        C = self.X0[self.coh[:, :3]]
        self.coh_area = 0.5 * np.linalg.norm(np.cross(C[:, 1] - C[:, 0], C[:, 2] - C[:, 0]), axis=1)
        self._omD = None
        c = self.curve
        i = int(np.argmax(c[:, 1]))
        self.peak, self.eps_peak = c[i, 1], c[i, 0]
        after = np.flatnonzero(c[i:, 1] <= 0.1 * self.peak)
        self.eps_f10 = c[i + after[0], 0] if len(after) else np.nan

    # ---------------- time ----------------
    def t_of(self, eps):
        return eps * self.L0 / self.grip

    def _w(self, tt, times):
        """Frame pair and weight for time tt (clamped to the frames)."""
        tt = float(np.clip(tt, times[0], times[-1]))
        k1 = int(np.searchsorted(times, tt))
        k1 = min(max(k1, 1), len(times) - 1)
        k0 = k1 - 1
        w = (tt - times[k0]) / (times[k1] - times[k0])
        if tt <= times[0]:
            return 0, 0, 0.
        return k0, k1, w

    def _lerp(self, A, tt, times):
        k0, k1, w = self._w(tt, times)
        a = A[k0].astype(np.float64)
        return a if w == 0 else (1 - w) * a + w * A[k1].astype(np.float64)

    def omega_D(self):
        """Nodal omega (0..1) and D (0..1) at every frame, rebuilt from PEEQ
        and triaxiality as in calcdamage.f (frame by frame, mid-point eta)."""
        if self._omD is not None:
            return self._omD
        K, N = self._P.shape
        P = np.nan_to_num(self._P.astype(np.float64))
        eta = np.nan_to_num(triax(self._S))
        ef = np.empty((K, N))
        uf = np.where(self.node_phase == 1, UF[1], UF[0])
        # element size: the solver's (6 V)^(1/3) of the tets around a node
        le = np.zeros(N); cnt = np.zeros(N)
        np.add.at(le, self.T.ravel(), np.repeat((6 * self.vol) ** (1 / 3), 4))
        np.add.at(cnt, self.T.ravel(), 1)
        le /= np.maximum(cnt, 1)
        om = np.zeros((K, N)); D = np.zeros((K, N))
        w = np.zeros(N); dd = np.zeros(N)
        for k in range(K):
            e_mid = eta[k] if k == 0 else 0.5 * (eta[k] + eta[k - 1])
            for ph in (0, 1):
                m = self.node_phase == ph
                ef[k, m] = np.interp(e_mid[m], EF_ETA, EF[ph])
            dp = np.maximum(P[k] - (P[k - 1] if k else 0.), 0.)
            dw = dp / ef[k]
            room = np.maximum(1 - w, 0.)
            used = np.minimum(dw, room)                     # the part of the increment before omega = 1
            post = np.where(dw > 0, (dw - used) / np.where(dw > 0, dw, 1) * dp, 0.)
            w = w + used
            dd = np.minimum(dd + le * post / uf, 1.)
            om[k], D[k] = w, dd
        self._omD = (om, D)
        return self._omD

    def at(self, eps, deform=1.0):
        """State at the grip strain eps (deform scales the displacements)."""
        s = State()
        s.eps = eps; s.tt = self.t_of(eps)
        s.U = np.nan_to_num(self._lerp(self._U, s.tt, self.t))
        s.X = self.X0 + deform * s.U
        s.alive = self.tdel > s.tt
        s.peeq = np.nan_to_num(self._lerp(self._P, s.tt, self.t))
        S = self._lerp(self._S, s.tt, self.t)
        s.mises, s.triax, s.sxx = mises(S), triax(S), np.nan_to_num(S[:, 0])
        om, D = self.omega_D()
        s.omega, s.D = self._lerp(om, s.tt, self.t), self._lerp(D, s.tt, self.t)
        op = self._lerp(self._open, s.tt, self.t_coh)
        s.debond = np.clip((op - self.coh_d0) / (self.coh_df - self.coh_d0), 0, 1)
        s.debond[self.coh_tdel <= s.tt] = 1.
        s.sigma = float(np.interp(eps, self.curve[:, 0], self.curve[:, 1]))
        return s

    # ---------------- geometry ----------------
    def surface(self, mask=None, X=None):
        """Boundary faces of the tets in mask (default: all).  Returns
        faces (M,3) node indices oriented outward, owner tet, kind
        (0: on the box, 1: crack / deleted region)."""
        idx = np.arange(len(self.T)) if mask is None else np.flatnonzero(mask)
        F = self.T[idx][:, FACES]                         # (n,4,3)
        opp = self.T[idx]                                 # node i is opposite face i
        key = np.sort(self.geo[F], axis=2).reshape(-1, 3)
        _, inv, cnt = np.unique(key, axis=0, return_inverse=True, return_counts=True)
        one = cnt[inv.ravel()] == 1
        faces = F.reshape(-1, 3)[one]
        owner = np.repeat(idx, 4)[one]
        far = opp.reshape(-1)[one]
        X = self.X0 if X is None else X
        a, b, c = X[faces[:, 0]], X[faces[:, 1]], X[faces[:, 2]]
        n = np.cross(b - a, c - a)
        flip = np.einsum('ij,ij->i', n, X[far] - a) > 0
        faces[flip] = faces[flip][:, [0, 2, 1]]
        P0 = self.X0[faces]
        tol = 1e-6
        on_box = np.zeros(len(faces), bool)
        for ax in range(3):
            for v in (self.lo[ax], self.hi[ax]):
                on_box |= np.all(np.abs(P0[:, :, ax] - v) < tol, axis=1)
        return faces, owner, np.where(on_box, 0, 1).astype(np.int8)

    def slice(self, axis, pos, state=None, field=None, mask=None):
        """Cut of the tets by the material plane x_axis = pos.  Returns a
        list of polygons (k,3) in current coordinates, the owner tets, and
        per-vertex values of `field` (a State attribute name or a nodal
        array) when given."""
        a = AXIS.get(axis, axis)
        X = self.X0 if state is None else state.X
        if mask is None:
            mask = np.ones(len(self.T), bool) if state is None else state.alive
        idx = np.flatnonzero(mask)
        T = self.T[idx]
        d = self.X0[T][:, :, a] - pos
        d[d == 0] = 1e-12
        s = d > 0
        n = s.sum(1)
        cut = (n > 0) & (n < 4)
        T, d, s, idx = T[cut], d[cut], s[cut], idx[cut]
        ea, eb = EDGES[:, 0], EDGES[:, 1]
        cross = s[:, ea] != s[:, eb]                       # (m,6)
        with np.errstate(divide='ignore', invalid='ignore'):
            f = np.where(cross, d[:, ea] / (d[:, ea] - d[:, eb]), 0.)
        pa, pb = X[T[:, ea]], X[T[:, eb]]
        P = pa + (pb - pa) * f[..., None]
        vals = None
        if field is not None:
            v = getattr(state, field) if isinstance(field, str) else field
            va, vb = v[T[:, ea]], v[T[:, eb]]
            vals = va + (vb - va) * f
        polys, owners, pv = [], [], []
        for m in (3, 4):
            sel = np.flatnonzero(cross.sum(1) == m)
            if not len(sel):
                continue
            order = np.argsort(~cross[sel], axis=1, kind='stable')[:, :m]
            Q = np.take_along_axis(P[sel], order[..., None], 1)
            if m == 4:                                    # order the quad around its centre
                c = Q.mean(1, keepdims=True)
                u, w_ = [i for i in range(3) if i != a]
                ang = np.arctan2(Q[:, :, w_] - c[:, :, w_], Q[:, :, u] - c[:, :, u])
                o2 = np.argsort(ang, axis=1)
                Q = np.take_along_axis(Q, o2[..., None], 1)
                order = np.take_along_axis(order, o2, 1)
            polys += list(Q); owners += list(idx[sel])
            if vals is not None:
                pv += list(np.take_along_axis(vals[sel], order, 1))
        return polys, np.array(owners, int), pv

    # ---------------- what fails ----------------
    def mechanisms(self, eps):
        """Fractions at grip strain eps: deleted volume of the matrix and of the
        hydrides (each of its own phase), debonded interface area (traction
        dropped to zero) of radial and of tangential plates."""
        tt = self.t_of(eps)
        dead = self.tdel <= tt
        out = {}
        m = self.phase == 0
        out['matrix'] = self.vol[m & dead].sum() / self.vol[m].sum()
        h = self.phase > 0
        out['hydride'] = self.vol[h & dead].sum() / self.vol[h].sum()
        s = self.at(eps)
        for ax, nm in ((0, 'debond_rad'), (1, 'debond_tan')):
            f = self.plate_axis[self.coh_plate] == ax
            out[nm] = (self.coh_area[f] * (s.debond[f] >= 1)).sum() / self.coh_area[f].sum() if f.any() else np.nan
        return out

    def check_omega(self):
        """Rebuilt omega + D against the solver's at the last state (tets
        left at the end): DUCT_IP of the solver runs 0..1 for omega, then
        1..2 as 1 + D."""
        om, D = self.omega_D()
        al = np.isfinite(self.omega_end)
        mine = (np.minimum(om[-1], 1) + D[-1])[self.T].mean(1)[al]
        ref = self.omega_end[al]
        return dict(corr=float(np.corrcoef(mine, ref)[0, 1]),
                    median_ratio=float(np.median(mine[ref > 0.05] / ref[ref > 0.05])),
                    agree_init=float(np.mean((mine >= 1) == (ref >= 1))))


def project(X, phi, theta, gamma=0., focal=None):
    """Manim ThreeDCamera view of points X (scene coordinates): returns
    (n,2) screen coordinates and depth (larger = nearer the camera)."""
    def rz(a):
        c, s = np.cos(a), np.sin(a)
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    c, s = np.cos(-phi), np.sin(-phi)
    rx = np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    R = rz(gamma) @ rx @ rz(-theta - np.pi / 2)
    P = X @ R.T
    if focal:
        f = focal / np.maximum(focal - P[:, 2], 1e-6)
        return P[:, :2] * f[:, None], P[:, 2]
    return P[:, :2], P[:, 2]
