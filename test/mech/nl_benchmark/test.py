"""Nonlocal model against the librarian's exact bar B1: where does W/G_f = 2 come from?

library (zirconium-lover/library, branch claude/ccx-arch2-damage-cohesion-
lsehkr-lib) notes/nonlocal-benchmarks.md sec. 1: this code's class - local
perfect plasticity, damage driven by the nonlocal plastic-strain increment,
L = l_d in the law (Bazant and Jirasek 2002 eqs. 61-63) - localises the plastic
flow into a set of zero measure, so in a bar with one weak element ALL the
plastic work after the peak goes into that element and the work to rupture is
    W / G_f = (l_d,h / l_d) * 1.0057     (INTEGRAL; 1.0057 = damage of the
                                          elastic neighbours, E=1e5, L=10)
with l_d,h/l_d = 1.768/1.764 at h = ell/4 (INTEGRAL), 2.266/2 (GRADIENT).
bar1d.py builds that bar: C3D8, one element across, 10 long, ell = 0.5,
h = ell/4, weak element sigma_y = 100 (bulk 102), u_f = 0.1, G_f = 5.

MEASURED on c4b0290, INTEGRAL (GRADIENT C3D8I: 1.377 against 1.140):
    family  W/G_f   F_max   plastic work / G_f   share in the weak element   elements flowing
    C3D8    1.985   101.96  1.82                 11 %                        33
    C3D8I   1.149   101.04  1.13                 33 %                         7
    expected 1.0075  100     1.00                 100 %                        1
The plastic flow does NOT localise into the weak element.  It cannot: a
uniaxial isochoric plastic stretch needs a lateral contraction of eps_p/2,
and an element whose faces are shared with elastic neighbours (nu = 0) is held
at their width - the force rises above sigma_y A (101-102 against 100) and the
neighbours yield.  Full-integration C3D8 adds volumetric locking on top, which
is why the incompatible-mode C3D8I halves the excess.  The 1-D result of the
reference assumes a uniaxial stress state the 3-D bar does not have, so the
deviation is the KINEMATICS (the specimen and the element family), not the
l_d calibration of the law: l_d,h/l_d would explain 0.2 %, not 99 %.

Consequences, recorded here because they decide how the numbers of
nl_objectivity are read (W_post/(G_f A) = 2.0 GRADIENT, 1.6 INTEGRAL on C3D4,
the family that locks worst): the dissipated energy converges in h, but to
G_f times a factor set by how far the flow spreads, which depends on the
element family and on the cross-section against l_d (C3D8I: 1.149 at ell =
0.5, 1.079 at ell = 1.0).

Checks.  The agreement with the reference is XFAIL on both families, with the
mechanism as a second XFAIL (the share of the plastic work in the weak
element, >= 90 % for the reference).  The PASS checks pin the measured
anatomy so a change in it is seen: rupture reached, C3D8I INTEGRAL at most
1.25 G_f, and C3D8I below C3D8 by at least 0.5 G_f.

RED PATH, run: MECH_RED_ENV=CCX_DAMAGE_NLWIDTH=0 (the element width in the
law instead of l_d) turns five checks red: no arm ruptures by the end of the
step, C3D8I INTEGRAL W/G_f 2.57 (ceiling 1.25), the C3D8-C3D8I gap 0.21.
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mechlib import Check, ccx, totals, blocks, trapz, write  # noqa: E402
import bar1d  # noqa: E402

WHAT = "nonlocal model vs exact 1-D bar B1 (librarian): W/G_f and where the plastic work goes"
SOURCE = "library notes/nonlocal-benchmarks.md sec. 1 (Bazant & Jirasek 2002 eqs. 61-63)"

NX, ELL, GF, UEND = 80, 0.5, 5.0, 0.15
REF = {"INTEGRAL": 1.768 / 1.764 * 1.0057, "GRADIENT": 2.266 / 2.0 * 1.0050}


def vm(r):
    sx, sy, sz, a, b, c = r[2:8]
    m = (sx + sy + sz) / 3
    return math.sqrt(1.5 * ((sx - m) ** 2 + (sy - m) ** 2 + (sz - m) ** 2 + 2 * (a * a + b * b + c * c)))


def arm(work, et, mode, exe):
    d = os.path.join(work, "%s_%s" % (et, mode))
    txt = bar1d.deck(NX, ell=ELL, uend=UEND).replace("TYPE=C3D8,", "TYPE=%s," % et)
    txt = txt.replace("*EL PRINT, ELSET=EWEAK\nPEEQ", "*EL PRINT, ELSET=EALL\nS, PEEQ")
    txt = txt.replace("*SOLID SECTION, ELSET=EBULK", "*ELSET, ELSET=EALL\nEBULK, EWEAK\n*SOLID SECTION, ELSET=EBULK")
    write(os.path.join(d, "t.inp"), txt)
    rc, _ = ccx(d, "t", exe, {"CCX_DAMAGE_CHARLEN": "1", "CCX_DAMAGE_NONLOCAL_MODE": mode})
    dat = os.path.join(d, "t.dat")
    c = [(0.0, 0.0)] + [(UEND * t, f) for t, f in totals(dat)]
    s, p = blocks(dat, "stresses"), blocks(dat, "equivalent plastic")
    ts = sorted(p)
    h = 10.0 / NX
    wp = {}
    # plastic work per element: (damaged) von Mises times d(peeq), 8 points of volume h/8
    for t0, t1 in zip(ts, ts[1:]):
        for r0, r1, q in zip(p[t0], p[t1], s[t1]):
            e = int(r1[0])
            wp[e] = wp.get(e, 0.0) + vm(q) * (r1[2] - r0[2]) * h / 8
    tot = sum(wp.values())
    return dict(rc=rc, w=trapz(c) / GF, fmax=max(f for _, f in c), fend=c[-1][1],
                wp=tot / GF, share=wp.get(NX // 2 + 1, 0.0) / tot if tot else 0.0,
                nflow=sum(1 for v in wp.values() if v > 0.01 * max(wp.values())) if wp else 0)


def run(work, exe):
    ck = []
    r = {(et, m): arm(work, et, m, exe) for et, m in (("C3D8", "INTEGRAL"), ("C3D8I", "INTEGRAL"),
                                                       ("C3D8I", "GRADIENT"))}
    for (et, m), x in r.items():
        tag = "%s %s" % (et, m)
        ck.append(Check("%s: ruptures (rc 0, F_end < 1 %% of peak)" % tag,
                        x["rc"] == 0 and abs(x["fend"]) < 0.01 * x["fmax"], True, cmp="true"))
        ck.append(Check("%s: W/G_f = (l_d,h/l_d) * elastic term" % tag, x["w"], REF[m], rtol=0.05,
                        expect_fail=True, why="plastic flow spreads past the weak element (3-D "
                        "lateral constraint%s)" % (", plus volumetric locking" if et == "C3D8" else "")))
        ck.append(Check("%s: share of plastic work in the weak element" % tag, x["share"], 0.9,
                        cmp="ge", expect_fail=True,
                        why="%d elements flow; plastic work %.2f G_f" % (x["nflow"], x["wp"])))
        ck.append(Check("%s: F_max / (sigma_y A)" % tag, x["fmax"] / 100.0, None, cmp="info"))
    ck.append(Check("C3D8I INTEGRAL: W/G_f ceiling (measured 1.149)", r[("C3D8I", "INTEGRAL")]["w"],
                    1.25, cmp="le", why="pins the measured excess so a growth of it is seen"))
    ck.append(Check("C3D8 - C3D8I (INTEGRAL): locking share of W/G_f",
                    r[("C3D8", "INTEGRAL")]["w"] - r[("C3D8I", "INTEGRAL")]["w"], 0.5, cmp="ge",
                    why="measured 0.84: the full-integration hexahedron locks under isochoric flow"))
    return ck
