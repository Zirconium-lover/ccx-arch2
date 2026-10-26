"""Bulk ductile damage on ONE element: initiation, evolution, width, energy.

The element is strained homogeneously along a fixed path (onelem.py), so the
stress state and the triaxiality are exact and the same for a C3D4 or a C3D8,
aligned or turned 45 deg about z.  Material: E = 2e5, nu = 0.3, J2 with NO
hardening, sigma_0 = 400.

Initiation (src/calcdamage.f, damageupdatepoint):
    omega += d(eps_p) / eps_f(eta),   eta = p / sigma_vM of the effective stress
    RICETRACEY: eps_f = 1.65 eps0 exp(-1.5 eta)   (Rice and Tracey 1969, the
                constant and exponent as coded; eps0 = 0.005)
    DUCTILE, NPOINTS=n: eps_f linear in eta between the table points, clamped
Checked INCREMENT BY INCREMENT - d(omega) against d(eps_p)/eps_f(eta_n) with
eta_n from the printed stress - so a non-proportional path is checked as
exactly as a proportional one.  Paths: uniaxial (eta=1/3), equibiaxial (2/3),
pure shear (0), and eps_y = eps_x/2 (between table points).

Evolution:  D = L (eps_p - eps_p,init) / u_f, so dD/d(eps_p) = L/u_f, with
    CHARLEN=0: L = (6V)^(1/3)  (C3D4 only)
    CHARLEN=1: L = extent of the element along the major principal direction
               (Jirasek and Bauer 2012, sec. 7) - here along x.
Element edge a = 0.5: aligned C3D8 L = 0.5; C3D8 at 45 deg L = 0.5*sqrt(2);
C3D4 at 45 deg: legacy 0.5, projection 0.5*sqrt(2).

Energy: with sigma_eff = sigma_0 fixed while the point flows, the external
work from initiation to D = 1 is exactly  V sigma_0 u_f / (2L) = (V/L) G_f,
G_f = sigma_0 u_f / 2: one G_f per crack area V/L (Bazant and Oh 1983).
EVOLUTION=ENERGY reads G_f from the card and uses u_f = 2 G_f / sigma_0.

FOUND, kept as XFAIL and reported (71ed70d): the triaxiality is computed
from the 2nd Piola-Kirchhoff stress, not the Cauchy stress.  On the paths
where S stays proportional to sigma (uniaxial, equibiaxial) the two agree and
the Cauchy check passes to 3e-6; on pure shear and eps_y = eps_x/2 it is off
by 1e-3..1e-2 in d(omega), while the same increments match eps_f(eta_PK2) to
2e-6 - the PK2 checks pin what the code does.  The ENERGY sigma_0 is the
PEAK von Mises of the same PK2 stress before initiation, which a non-hardening
material reaches at FIRST YIELD (S = J sigma/lam^2 falls as it stretches):
u_f +0.37 % here, matched to 1e-6 by 2 G_f / max vM(S).  Both grow with the
strain at initiation, O(eps): small in this test, not in a necking tube.

RED PATH, run: MECH_RED_ENV=CCX_DAMAGE_CHARLEN=0 gives the C3D8 arms rc=201
(the legacy width refuses hexahedra) and the C3D4-at-45-deg CHARLEN=1 arm
L = 0.5 against 0.7071 and W +40 %: 13 checks red.  MECH_RED_ENV=
CCX_DAMAGE_VISCOSITY=0.05 (D lags behind the law) turns every energy check
red, W +23..32 %.

Hardening after initiation (library notes/bulk-damage.md sec. 2b): with
sigma_eff = sigma_i + H x during softening the work per crack area is
G_f(sigma_i) (1 + H u_f / (3 sigma_i L)), so the crack band is not objective
in L.  H = 4000, one C3D8 at a = L = 0.5 and 0.25: 1.0832 and 1.1443 against
1.0635 and 1.1269 (+1.9 / +1.5 %, the NLGEOM and elastic-effective-strain
bias, 3 % allowed); the L-dependent part itself, the difference of the two
excesses, 0.0611 against 0.0634 (15 % allowed).  A law without the term
(both excesses 0) fails that check by the whole 0.061.

Tolerances.  Initiation and slope: the .dat prints 7 digits, but the slope
is a difference of printed values: 2e-4.  Energy: the run is NLGEOM (the
damage path needs it), so the force carries the area change of a plastic
strain of ~1-3 % and the work integral differs from the small-strain identity
by that order: 3 %; measured 0.6-1.5 %.
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mechlib import Check, ccx, blocks, write  # noqa: E402
import onelem  # noqa: E402

WHAT = "bulk damage, one element: omega vs eps_f(eta) (RT, table), dD/deps_p = L/u_f (CHARLEN 0/1, 45 deg), W = (V/L) G_f"
SOURCE = "calcdamage.f law; Rice & Tracey 1969 form as coded; Jirasek & Bauer 2012 sec. 7; Bazant & Oh 1983"

E, NU, S0, EPS0, UF, A = 200000.0, 0.3, 400.0, 0.005, 0.01, 0.5
TABLE = ((0.0, 0.012), (1.0 / 3, 0.006), (2.0 / 3, 0.004))


def material(kind="RT", evolution="DISPLACEMENT", gf=None, hmod=0.0):
    m = ["*MATERIAL, NAME=M", "*ELASTIC", "%g, %g" % (E, NU), "*PLASTIC", "%g, 0." % S0]
    if hmod:
        m.append("%g, 1." % (S0 + hmod))
    if kind == "RT":
        m += ["*DAMAGE INITIATION, CRITERION=RICETRACEY, EVOLUTION=%s" % evolution,
              "%g, 1.0, %g, 0." % (EPS0, gf if gf else UF)]
    else:
        m += ["*DAMAGE INITIATION, CRITERION=DUCTILE, EVOLUTION=DISPLACEMENT, NPOINTS=%d" % len(TABLE),
              # 2 + 2n constants, eight per line, then the temperature
              "1.0, %g, " % UF + ", ".join("%.12g, %.12g" % p for p in TABLE), "0."]
    return m


def ef_table(eta):
    if eta <= TABLE[0][0]:
        return TABLE[0][1]
    for (e1, f1), (e2, f2) in zip(TABLE, TABLE[1:]):
        if eta <= e2:
            return f1 + (f2 - f1) * (eta - e1) / (e2 - e1)
    return TABLE[-1][1]


def history(d):
    dat = os.path.join(d, "o.dat")
    s, p, w = blocks(dat, "stresses"), blocks(dat, "equivalent plastic"), blocks(dat, "ductile damage")
    u, f = blocks(dat, "displacements"), blocks(dat, "forces")
    out = []
    for t in sorted(s):
        if not (s[t] and p.get(t) and w.get(t) and u.get(t) and f.get(t)):
            break                               # the element has been deleted
        sx, sy, sz, sxy, sxz, syz = s[t][0][2:8]
        pm = (sx + sy + sz) / 3
        vm = math.sqrt(1.5 * ((sx - pm) ** 2 + (sy - pm) ** 2 + (sz - pm) ** 2
                              + 2 * (sxy ** 2 + sxz ** 2 + syz ** 2)))
        # stretch: u_c = eps_c * x_c for every node, so node 1 (at the shift
        # point, every coordinate nonzero) gives lam_c = 1 + u_c / x_c
        x1 = onelem.nodes("C3D8")[0]
        lam = [1 + u[t][0][1 + c] / x1[c] for c in range(3)]
        out.append(dict(t=t, eta=pm / vm if vm > 1e-9 else 0.0, pe=p[t][0][2], dam=w[t][0][2],
                        sig=(sx, sy, sz), lam=lam,
                        u=[r[1:4] for r in u[t]], f=[r[1:4] for r in f[t]]))
    return out


def one(work, name, eltype, ratios, mat, exe, rotz=0.0, charlen="1", emax=0.05, dt=0.005, a=A):
    d = os.path.join(work, name)
    write(os.path.join(d, "o.inp"), onelem.deck(eltype, ratios, emax, mat, dt, a=a, rotz=rotz))
    rc, _ = ccx(d, "o", exe, {"CCX_DAMAGE_CHARLEN": charlen})
    return rc, (history(d) if rc == 0 else [])


def pk2(x):
    """triaxiality and von Mises of the 2nd Piola-Kirchhoff stress,
    S_i = J sigma_i / lam_i^2 for a pure stretch along the axes."""
    lam = x["lam"]
    j = lam[0] * lam[1] * lam[2]
    sv = [j * x["sig"][c] / lam[c] ** 2 for c in range(3)]
    m = sum(sv) / 3
    vm = math.sqrt(1.5 * sum((v - m) ** 2 for v in sv))
    return (m / vm if vm > 1e-9 else 0.0), vm


def init_check(name, h, eff, measure="cauchy", **kw):
    """max relative error of d(omega) against d(eps_p)/eps_f(eta), pre-initiation."""
    err, n = 0.0, 0
    for a, b in zip(h, h[1:]):
        dp = b["pe"] - a["pe"]
        if b["dam"] >= 1.0 or dp < 1e-5:
            continue
        eta = b["eta"] if measure == "cauchy" else pk2(b)[0]
        err = max(err, abs((b["dam"] - a["dam"]) / (dp / eff(eta)) - 1))
        n += 1
    return Check("%s: d(omega) = d(eps_p)/eps_f(eta_%s)  [%d inc]" % (name, measure, n),
                 err if n else None, 0.0, tol=2e-4, **kw)


def slope(h):
    post = [x for x in h if 1.02 < x["dam"] < 1.97]
    a, b = post[0], post[-1]
    return (b["dam"] - a["dam"]) / (b["pe"] - a["pe"])


def work_after_init(h):
    def w_upto(k):
        w = 0.0
        for a, b in zip(h[:k], h[1:k + 1]):
            for ua, ub, fa, fb in zip(a["u"], b["u"], a["f"], b["f"]):
                w += sum(0.5 * (fa[c] + fb[c]) * (ub[c] - ua[c]) for c in range(3))
        return w
    ki = next(i for i, x in enumerate(h) if x["dam"] >= 1.0)
    # D = 1, or the last state before deletion (at D_vis >= 0.999)
    kf = next((i for i, x in enumerate(h) if x["dam"] >= 2.0), len(h) - 1)
    # initiation inside increment ki: interpolate the work linearly in omega
    a, b = h[ki - 1], h[ki]
    wa, wb = w_upto(ki - 1), w_upto(ki)
    frac = (1.0 - a["dam"]) / max(b["dam"] - a["dam"], 1e-30)  # dam is omega up to 1, then 1+D
    frac = min(max(frac, 0.0), 1.0)
    return w_upto(kf) - (wa + frac * (wb - wa))


def run(work, exe):
    ck = []
    rt = lambda eta: 1.65 * EPS0 * math.exp(-1.5 * eta)
    paths = (("uniaxial", (1, None, None), 1 / 3), ("equibiaxial", (1, 1, None), 2 / 3),
             ("pure shear", (1, "iso", None), 0.0), ("eps_y=eps_x/2", (1, 0.5, None), None))
    for tag, ratios, eta0 in paths:
        for kind, eff in (("RT", rt), ("table", ef_table)):
            rc, h = one(work, "%s_%s" % (kind, tag.replace(" ", "_").replace("/", "")), "C3D8",
                        ratios, material(kind), exe)
            ck.append(Check("%s %s: run rc" % (kind, tag), rc, 0, tol=0))
            if rc:
                continue
            if eta0 is not None and kind == "RT":
                etas = [x["eta"] for x in h if x["pe"] > 0 and x["dam"] < 2]
                ck.append(Check("%s: triaxiality while flowing" % tag, max(abs(e - eta0) for e in etas),
                                0.0, tol=2e-3, why="the path fixes eta, up to the elastic volume change (1.2e-3 on pure shear)"))
            # eta of the Cauchy stress is the triaxiality; the code takes it from
            # the 2nd PK stress.  The two agree exactly on the uniaxial and
            # equibiaxial paths (S stays proportional to sigma) and differ at
            # O(strain) on the others: XFAIL there, reported; the PK2 check
            # pins what the code does meanwhile.
            exact = tag in ("uniaxial", "equibiaxial")
            ck.append(init_check("%s %s" % (kind, tag), h, eff, expect_fail=not exact,
                                 why="" if exact else "eta from the 2nd PK stress, not Cauchy (reported)"))
            if not exact:
                ck.append(init_check("%s %s" % (kind, tag), h, eff, measure="pk2"))
    # width and slope: C3D8 aligned / 45 deg, C3D4 45 deg legacy and projection
    cases = (("C3D8 aligned", "C3D8", 0.0, "1", A),
             ("C3D8 45deg", "C3D8", math.pi / 4, "1", A * math.sqrt(2)),
             ("C3D4 45deg CHARLEN=0", "C3D4", math.pi / 4, "0", A),
             ("C3D4 45deg CHARLEN=1", "C3D4", math.pi / 4, "1", A * math.sqrt(2)))
    for tag, et, rz, cl, lref in cases:
        rc, h = one(work, tag.replace(" ", "_").replace("=", ""), et, (1, None, None),
                    material(), exe, rotz=rz, charlen=cl, emax=0.06, dt=0.004)
        ck.append(Check("%s: run rc" % tag, rc, 0, tol=0))
        if rc:
            continue
        ck.append(Check("%s: L = u_f dD/d(eps_p)" % tag, UF * slope(h), lref, rtol=2e-4))
        vol = A ** 3 / (6.0 if et == "C3D4" else 1.0)
        ck.append(Check("%s: W after init / ((V/L) G_f)" % tag,
                        work_after_init(h) / (vol / lref * S0 * UF / 2), 1.0, rtol=0.03))
    # hardening AFTER initiation: sigma = (1-D)(sigma_i + H x), D = L x/u_f, so
    # the work per crack area is  G_f (1 + H u_f / (3 sigma_i L))  - a term
    # that depends on the band width L (library notes/bulk-damage.md sec. 2b)
    hm = 4000.0
    exc = {}
    for a in (0.5, 0.25):
        rc, h = one(work, "hardening_a%g" % a, "C3D8", (1, None, None), material(hmod=hm), exe,
                    emax=0.08 if a == 0.25 else 0.05, dt=0.0025, a=a)
        ck.append(Check("hardening L=%g: run rc" % a, rc, 0, tol=0))
        if rc:
            continue
        k = next(i for i, x in enumerate(h) if x["dam"] >= 1.0)
        x0, x1 = h[k - 1], h[k]
        epi = x0["pe"] + (1 - x0["dam"]) / (x1["dam"] - x0["dam"]) * (x1["pe"] - x0["pe"])
        si = S0 + hm * epi
        gfi = si * UF / 2
        pred = 1 + hm * UF / (3 * si * a)
        got = work_after_init(h) / (a * a * gfi)
        exc[a] = (got, pred)
        ck.append(Check("hardening L=%g: W/(A G_f(sigma_i)) vs 1 + H u_f/(3 sigma_i L)" % a, got, pred,
                        rtol=0.03, why="sigma_i = %.1f at initiation" % si))
    if len(exc) == 2:
        ck.append(Check("hardening: excess(L=0.25) - excess(L=0.5)", exc[0.25][0] - exc[0.5][0],
                        exc[0.25][1] - exc[0.5][1], rtol=0.15,
                        why="the L-dependent term itself; the NLGEOM bias common to both cancels"))
    # EVOLUTION=ENERGY: u_f = 2 G_f / sigma_0 at initiation
    gf = 3.0                                    # -> u_f = 0.015
    rc, h = one(work, "energy", "C3D8", (1, None, None), material(evolution="ENERGY", gf=gf), exe,
                emax=0.08, dt=0.004)
    ck.append(Check("ENERGY: run rc", rc, 0, tol=0))
    if rc == 0:
        ck.append(Check("ENERGY: u_f = L/(dD/d eps_p) vs 2 G_f/sigma_0", A / slope(h), 2 * gf / S0,
                        rtol=2e-4, expect_fail=True,
                        why="sigma_0 is the PEAK von Mises of the 2nd PK stress before initiation, "
                        "reached at first yield (reported)"))
        s_pk2 = max(pk2(x)[1] for x in h if x["dam"] < 1.0)
        ck.append(Check("ENERGY: u_f vs 2 G_f / max vM(S_PK2) pre-init (as coded)", A / slope(h),
                        2 * gf / s_pk2, rtol=2e-4))
        ck.append(Check("ENERGY: W after init / ((V/L) G_f)", work_after_init(h) / (A ** 2 * gf),
                        1.0, rtol=0.03))
    return ck
