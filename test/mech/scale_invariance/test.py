"""Scale invariance of the element-level laws, x1 against x1/25.

library notes/scaling-levers.md sec. 3: when geometry is shrunk by k the
PHYSICAL parameters (E, sigma_y, eps_f, G_f, Tn0, Gc) stay, the NUMERICAL ones
are rescaled to keep their groups (Kn t/E, mu v/delta0, ...).  A law that has
no hidden absolute length must then give the same dimensionless answer.

(a) Bulk law, one C3D8 under uniaxial stress (onelem.py), edge a = 0.5 and
    0.02.  With u_f scaled as the size (the artificial "G_f proportional to
    D") the deck is the same problem in other units: sigma(eps), omega and D
    must agree to the printed digits.  With the physical u_f kept, the
    crack-band law turns the same G_f into a softening strain u_f/L that is
    k times larger - the element-level size effect - while the work per crack
    area stays (V/L) G_f on both scales.
(b) UC6 facet (uc6_law/uc6deck.py), edge 1 and 1/25, same separation path:
    the law knows no length but the separations, so T(delta) and W/Gc agree to
    the printed digits.  With Kn -> k Kn (the rule that keeps alpha = Kn t/E)
    delta0 falls by k and W stays Gc (W_exact moves only through the gmin
    sliver).
(c) UC6 viscosity: the group mu v/delta0 decides the overshoot.  Constant
    opening rate v = 1, Kn = 1e5, Tn0 = 50, Gc = 0.3; exact Duvaut-Lions
    values from library notes/scripts/cohesive_ref.py (viscous_ramp, Ei
    solution): group 0.01 / 0.1 / 1 -> T_max/Tn0 1.0085 / 1.0772 / 1.5827,
    W/Gc 1.0028 / 1.0285 / 1.3597.  Then the two groups of the disk and
    the cell deck, 0.067 and 1.25: the cell's interfaces are ~1.7x stronger and
    tougher than their card says because mkdeck did not rescale mu.
    Tolerance 1 %: 2000 backward-Euler steps over the ramp (the librarian's
    200-step Euler value at group 1 is 0.7 % off the exact W).

Measured on ef52329: (a) sigma and dam identical to the last digit; rate
ratio 25.3 (the elastic share of the strain rate, 2 % allowed); (b) T to
8.8e-8 Tn0, W ratio 1.000000, Kn x25 W/W_exact 0.99928; (c) 0.01 / 0.1 / 1:
T_max 1.0074 / 1.0758 / 1.5809, W 1.0028 / 1.0286 / 1.3607 - and at the deck
groups: disk (0.067) T_max/Tn0 1.052, W/Gc 1.019; cell (1.25) 1.702 and
1.477.

(e) Elastic bar + one UC6 plane (mkcohesive.py, length 4, section 1),
    x1 and x1/10.  With Kn x10 (alpha kept), Gc/10 and the displacement /10
    every separation scales and the run is the same run: .sta identical,
    W/(Gc A) 1.00011 on both.  With the physics kept (Gc, Tn0; Kn x10 by the
    alpha rule) the interface still dissipates Gc per area (0.996) and the
    strength is the same, but the bar breaks at a 10x larger Delta/L.

RED PATH.  The identity checks are the same comparison as the third bulk
arm, which differs from x1 by the factor k it is asserted to; run with
MECH_RED_ENV=CCX_DAMAGE_CHARLEN=0 the bulk arms stop rc=201 (the legacy
width refuses C3D8) and the test is red.
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "bulk_law"))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "uc6_law"))
from mechlib import Check, ccx, totals, blocks, trapz, write  # noqa: E402
import onelem  # noqa: E402
import uc6deck  # noqa: E402

WHAT = "scale x1 vs x1/25: bulk law and UC6 facet invariant in their groups; UC6 overshoot vs mu v/delta0"
SOURCE = "library notes/scaling-levers.md sec. 3; cohesive_ref.py viscous_ramp (exact Duvaut-Lions)"

K = 25.0
E, NU, S0, EPS0 = 200000.0, 0.3, 400.0, 0.005


def bulk(work, name, a, uf, exe):
    mat = ["*MATERIAL, NAME=M", "*ELASTIC", "%g, %g" % (E, NU), "*PLASTIC", "%g, 0." % S0,
           "*DAMAGE INITIATION, CRITERION=RICETRACEY, EVOLUTION=DISPLACEMENT",
           "%g, 1.0, %.12g, 0." % (EPS0, uf)]
    d = os.path.join(work, name)
    write(os.path.join(d, "o.inp"), onelem.deck("C3D8", (1, None, None), 0.05, mat, 0.0025, a=a))
    rc, _ = ccx(d, "o", exe, {"CCX_DAMAGE_CHARLEN": "1"})
    dat = os.path.join(d, "o.dat")
    s, w = blocks(dat, "stresses"), blocks(dat, "ductile damage")
    hist = [(t, s[t][0][2], w[t][0][2]) for t in sorted(s) if s[t] and w.get(t)]
    return rc, hist


def facet(work, name, a, kn, path, dt, exe, mu=0.0, tn0=400.0, ts0=300.0, gc=2.0, gmin=1e-3):
    d = os.path.join(work, name)
    write(os.path.join(d, "u.inp"), uc6deck.deck(path, (kn, tn0, ts0, gc, gmin, mu), dt, a=a))
    rc, _ = ccx(d, "u", exe)
    area = a * a / 2
    f = totals(os.path.join(d, "u.dat"), "NPLUS", 0) if rc == 0 else []
    return rc, [(t, x / area) for t, x in f]


def work_of(path, tr):
    """work per area along a pure-normal path, trapezoidal in the printed points"""
    pts = [(0.0, 0.0)] + [(path[-1][1][0] * t / path[-1][0], T) for t, T in tr]
    return trapz(pts)


def run(work, exe):
    ck = []
    # (a) bulk law
    rc1, h1 = bulk(work, "bulk_x1", 0.5, 0.01, exe)
    rc2, h2 = bulk(work, "bulk_xk_uf_scaled", 0.5 / K, 0.01 / K, exe)
    rc3, h3 = bulk(work, "bulk_xk_uf_kept", 0.5 / K, 0.01, exe)
    ck.append(Check("bulk: runs (rc)", max(abs(rc1), abs(rc2), abs(rc3)), 0, tol=0))
    if not (rc1 or rc2):
        n = min(len(h1), len(h2))
        ds = max(abs(x[1] - y[1]) for x, y in zip(h1[:n], h2[:n])) / S0
        dd = max(abs(x[2] - y[2]) for x, y in zip(h1[:n], h2[:n]))
        ck.append(Check("bulk x1/25, u_f/25: max |sigma - sigma_x1| / sigma_0", ds, 0.0, tol=1e-5))
        ck.append(Check("bulk x1/25, u_f/25: max |dam - dam_x1| (omega, 1+D)", dd, 0.0, tol=1e-5))
    if not (rc1 or rc3):
        def slope(h):     # dD/dt after initiation, t = strain/0.05 (same path on both)
            post = [(t, w) for t, _, w in h if 1.005 < w < 1.97]
            return (post[-1][1] - post[0][1]) / (post[-1][0] - post[0][0])
        ck.append(Check("bulk x1/25, u_f kept: softening rate ratio x1 / x1/25 (= k)", slope(h1) / slope(h3), K,
                        rtol=0.02, why="D = L eps_p/u_f with L = a/k: the same G_f needs k times the strain"))
    # (b) UC6 facet, same path at two sizes, then Kn x k
    path = [(0, (0, 0, 0)), (1.2, (0.012, 0, 0))]
    rb1, t1 = facet(work, "uc6_x1", 1.0, 1e6, path, 0.01, exe)
    rb2, t2 = facet(work, "uc6_xk", 1.0 / K, 1e6, path, 0.01, exe)
    ck.append(Check("UC6: runs (rc)", max(abs(rb1), abs(rb2)), 0, tol=0))
    if not (rb1 or rb2):
        ck.append(Check("UC6 x1/25: max |T - T_x1| / Tn0", max(abs(x[1] - y[1]) for x, y in zip(t1, t2)) / 400.0,
                        0.0, tol=1e-5))
        ck.append(Check("UC6 x1/25: W / W_x1", work_of(path, t2) / work_of(path, t1), 1.0, tol=1e-5))
    rb3, t3 = facet(work, "uc6_xk_knk", 1.0 / K, 1e6 * K, path, 0.0025, exe)
    if rb3 == 0:
        # work to d_f of the bilinear law with the gmin floor, closed form (uc6_law.w_exact)
        kn, tn0, gc, gmin = 1e6 * K, 400.0, 2.0, 1e-3
        d0, df = tn0 / kn, 2 * gc / tn0
        sl = tn0 / (df - d0)
        dc = sl * df / (gmin * kn + sl)
        wex = kn * d0 ** 2 / 2 + (tn0 + sl * (df - dc)) / 2 * (dc - d0) + gmin * kn * (df ** 2 - dc ** 2) / 2
        pts = [(0.0, 0.0)] + [(0.01 * t, T) for t, T in t3 if 0.01 * t <= df + 1e-12]
        ck.append(Check("UC6 x1/25, Kn x25: W to d_f / W_exact", trapz(pts) / wex, 1.0, tol=2e-3,
                        why="d0 falls by k, Gc stays; the gmin floor (gmin Kn grew x25) is in W_exact"))
    # (c) viscosity group
    kn, tn0, gc = 1e5, 50.0, 0.3
    d0 = tn0 / kn
    exact = {0.01: (1.0085, 1.0028), 0.1: (1.0772, 1.0285), 1.0: (1.5827, 1.3597)}
    res = {}
    for g in (0.01, 0.1, 1.0, 0.067, 1.25):
        mu = g * d0                         # v = 1
        tend = 1.2e-2 * 1.5 + 10 * mu
        p = [(0, (0, 0, 0)), (tend, (tend, 0, 0))]
        rc, tr = facet(work, "visc_%g" % g, 1.0, kn, p, tend / 2000, exe, mu=mu, tn0=tn0, ts0=tn0, gc=gc,
                       gmin=1e-6)
        if rc:
            ck.append(Check("visc group %g: run rc" % g, rc, 0, tol=0))
            continue
        res[g] = (max(T for _, T in tr) / tn0, work_of(p, tr) / gc)
    for g, (tm, wr) in exact.items():
        if g in res:
            ck.append(Check("mu v/d0 = %g: T_max/Tn0" % g, res[g][0], tm, rtol=0.01))
            ck.append(Check("mu v/d0 = %g: W/Gc" % g, res[g][1], wr, rtol=0.01))
    for g, tag in ((0.067, "disk deck"), (1.25, "cell deck (mu not rescaled)")):
        if g in res:
            ck.append(Check("mu v/d0 = %g (%s): T_max/Tn0" % (g, tag), res[g][0], None, cmp="info"))
            ck.append(Check("mu v/d0 = %g (%s): W/Gc" % (g, tag), res[g][1], None, cmp="info"))
    # (e) elastic bar + one UC6 interface (test/pathfollow/mkcohesive.py), x1 and x1/10
    ck += bar_cohesive(work, exe)
    return ck


def bar_cohesive(work, exe):
    import filecmp
    import subprocess
    from mechlib import ROOT
    ck, r = [], {}
    for tag, k, kn, gc, ue in (("x1", 1, 1e6, 2.0, 0.02), ("xk_scaled", 10, 1e7, 0.2, 0.002),
                               ("xk_physics", 10, 1e7, 2.0, 0.02)):
        d = os.path.join(work, "coh_" + tag)
        os.makedirs(d, exist_ok=True)
        subprocess.check_call(["python3", os.path.join(ROOT, "test/pathfollow/mkcohesive.py"), "--nslice", "4",
                               "--isplit", "2", "--h", repr(1.0 / k), "--kn", repr(kn), "--gc", repr(gc),
                               "--uend", repr(ue), "--gmin", "1e-6", "-o", os.path.join(d, "t.inp")],
                              stdout=subprocess.DEVNULL)
        rc, _ = ccx(d, "t", exe)
        c = [(0.0, 0.0)] + [(ue * t, f) for t, f in totals(os.path.join(d, "t.dat"), "NRIGHT")]
        a, length = (1.0 / k) ** 2, 4.0 / k
        fm = max(f for _, f in c)
        ubreak = next((u for u, f in c if u > 0 and f < 0.01 * fm and c.index((u, f)) > 1), c[-1][0])
        r[tag] = dict(d=d, rc=rc, w=trapz(c) / (gc * a), fmax=fm / (400.0 * a), ub=ubreak / length)
    ck.append(Check("UC6 bar: runs (rc)", max(abs(x["rc"]) for x in r.values()), 0, tol=0))
    same = filecmp.cmp(os.path.join(r["x1"]["d"], "t.sta"), os.path.join(r["xk_scaled"]["d"], "t.sta"),
                       shallow=False)
    ck.append(Check("UC6 bar x1/10 with Kn x10, Gc/10, Delta/10: .sta identical", same, True, cmp="true"))
    ck.append(Check("UC6 bar x1/10 scaled: W/(Gc A) vs x1", r["xk_scaled"]["w"], r["x1"]["w"], rtol=1e-6))
    ck.append(Check("UC6 bar x1/10 physics kept (Gc, Kn x10): W/(Gc A)", r["xk_physics"]["w"], 1.0, rtol=0.01,
                    why="the interface dissipates Gc per area at any size"))
    ck.append(Check("UC6 bar x1/10 physics kept: strength F_max/(Tn0 A)", r["xk_physics"]["fmax"],
                    r["x1"]["fmax"], rtol=0.03, why="strength regime on both; 3 % for the sampling of the peak"))
    ck.append(Check("UC6 bar x1/10 physics kept: Delta_break/L x1/10 over x1 (= 10, size effect)",
                    r["xk_physics"]["ub"] / r["x1"]["ub"], 10.0, rtol=0.15,
                    why="the same delta_f on a bar 10x shorter"))
    return ck
