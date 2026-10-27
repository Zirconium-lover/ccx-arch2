"""A fracturing bar at x1 and x1/25: what is the same problem in other units, and what is size effect.

The C3D4 bar of neck_dissipation (library B1: weak slice, H = 0, sigma_y =
100, G_f = sigma_y u_f/2), length 10 and section 1, scaled by 1/k, k = 25,
with the mesh (h = L/40 or L/80).  Viscosity 1e-3 and the separation end on
every arm; the step in pseudo-time is the same on both scales, which keeps the
viscous groups (eta d_dot/u_f) where the library's rules put them when u_f is
scaled with the size.

SAME PROBLEM IN OTHER UNITS ("G_f proportional to D", the library's check that
the numerical parameters are rescaled right): u_f -> u_f/k, ell -> ell/k,
end displacement -> /k.  The runs must be the same run: .sta byte-identical.
    LOCAL CHARLEN=1:   identical while the bar carries load (W/(G_f A) 5.5095
                       on both); the one increment after separation converges
                       in 3 iterations against 2 - absolute floors of the
                       convergence test at zero force, see loaded_sta()
    INTEGRAL ell=0.53: identical
A tie, found here and kept out on purpose: at ell = 0.5 the averaging radius
2 ell is exactly 8 h, so centroids sit ON the radius and rounding decides the
neighbour set differently at the two scales (82.7 against 81.8 neighbours per
element); the x1/25 run then stops rc=201 at theta = 0.32 where x1 finishes.
Not a hidden absolute length - at ell = 0.53 both are byte-identical - but
a radius commensurate with the mesh is a sensitivity worth knowing.

SIZE EFFECT (physics kept: u_f and G_f the same, bar 25x smaller) - recorded,
not "fixed", per library notes/scaling-levers.md sec. 1:
    strength F_max/(sigma_y A): the same (1.019 both: strength regime)
    LOCAL crack band: the softening strain u_f/h goes from 0.4 to 10 - the
        element must stretch tenfold before D = 1; W/(G_f A) 5.51 -> 2.04 and
        the bar is not severed at Delta = 2 L (2.6 % of the strength left)
    INTEGRAL, ell/k (numerical length): W/(G_f A) 1.389 -> 0.963,
        Delta_end/L 0.0097 -> 0.75: the same energy makes a 25x shorter bar
        ~75x more ductile in Delta/L
    INTEGRAL, ell kept (material length 0.53 against a bar 0.4 long): W/(G_f A)
        0.40 - the averaging covers the whole bar.

RED PATH, run: MECH_RED_ENV=CCX_DAMAGE_NONLOCAL=0.05 - an ABSOLUTE internal
length forced on every arm, the very thing a scale test exists to catch -
turns all four identity checks red (W/(G_f A) 1.41 against 5.52 on the local
pair, 0.88 against 5.53 on the integral pair).

Tolerances: identity is byte-level on the loaded increments; strength 1 %; the ductility gain >= 10x.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "nl_benchmark"))
from mechlib import Check, ccx, totals, trapz, write  # noqa: E402
import bar1d  # noqa: E402

WHAT = "fracturing C3D4 bar x1 vs x1/25: identical when G_f ~ D; size effect recorded when physics kept"
SOURCE = "library notes/scaling-levers.md sec. 1, 3, 5 (size effect is physics; numerical groups rescaled)"

K, SY, ELL = 25.0, 100.0, 0.53


def arm(work, tag, mode, nx, ell, k, uf, uend, exe):
    d = os.path.join(work, tag)
    dt = 0.00075 / 0.8 if mode == "LOCAL" else 0.0025
    write(os.path.join(d, "t.inp"), bar1d.deck(nx, length=10 / k, width=1 / k, ell=ell, uf=uf,
                                               eltype="C3D4", uend=uend, dt=dt))
    env = {"CCX_DAMAGE_CHARLEN": "1", "CCX_FRACTURE_TERMINATION": "NLEFT:NRIGHT",
           "CCX_DAMAGE_VISCOSITY": "1.e-3"}
    if mode != "LOCAL":
        env["CCX_DAMAGE_NONLOCAL_MODE"] = mode
    rc, _ = ccx(d, "t", exe, env)
    c = [(0.0, 0.0)] + [(uend * t, f) for t, f in totals(os.path.join(d, "t.dat"))]
    a = (1 / k) ** 2
    return dict(d=d, rc=rc, w=trapz(c) / (SY * uf / 2 * a), fmax=max(f for _, f in c) / (SY * a),
                uend=c[-1][0] / (10 / k), fend=c[-1][1] / (SY * a))


def loaded_sta(rec):
    """.sta lines of the increments that still carry > 1e-3 of the peak force.
    After separation the force is ~0 and the convergence test meets the
    solver's ABSOLUTE floors, which a force 625x smaller crosses differently
    (measured: the last LOCAL increment converges in 2 iterations at x1/25 and
    3 at x1, after rupture) - that is the solver's convergence criterion, not
    the mechanics."""
    f = dict(totals(os.path.join(rec["d"], "t.dat")))
    fmax = max(abs(v) for v in f.values())
    out = []
    for line in open(os.path.join(rec["d"], "t.sta")):
        x = line.split()
        if len(x) >= 5 and x[0].isdigit():
            t = float(x[4])
            if abs(f.get(t, fmax)) > 1e-3 * fmax:
                out.append(line)
    return out


def run(work, exe):
    ck = []
    r = {}
    r["L1"] = arm(work, "L1", "LOCAL", 40, None, 1, 0.1, 0.8, exe)
    r["L2"] = arm(work, "L2", "LOCAL", 40, None, K, 0.1 / K, 0.8 / K, exe)
    r["L3"] = arm(work, "L3", "LOCAL", 40, None, K, 0.1, 0.8, exe)
    r["I1"] = arm(work, "I1", "INTEGRAL", 80, ELL, 1, 0.1, 0.3, exe)
    r["I2"] = arm(work, "I2", "INTEGRAL", 80, ELL / K, K, 0.1 / K, 0.3 / K, exe)
    r["I3"] = arm(work, "I3", "INTEGRAL", 80, ELL / K, K, 0.1, 0.3, exe)
    r["I4"] = arm(work, "I4", "INTEGRAL", 80, ELL, K, 0.1, 0.3, exe)
    for a, b, tag in (("L1", "L2", "LOCAL"), ("I1", "I2", "INTEGRAL")):
        la, lb = loaded_sta(r[a]), loaded_sta(r[b])
        ck.append(Check("%s x1/25 with u_f, ell, Delta /25: .sta identical while loaded (%d inc)" % (tag, len(la)),
                        la == lb, True, cmp="true"))
        ck.append(Check("%s x1/25 scaled: W/(G_f A) vs x1" % tag, r[b]["w"], r[a]["w"], rtol=1e-6))
    ck.append(Check("x1 arms rupture (rc)", max(abs(r[t]["rc"]) for t in ("L1", "I1", "I3", "I4")), 0, tol=0))
    for a, b in (("L1", "L3"), ("I1", "I3")):
        ck.append(Check("%s -> %s (physics kept): strength F_max/(sigma_y A) unchanged" % (a, b),
                        r[b]["fmax"], r[a]["fmax"], rtol=0.01))
    ck.append(Check("INTEGRAL ell/k, physics kept: ductility Delta_end/L grows >= 10x (size effect)",
                    r["I3"]["uend"] / r["I1"]["uend"], 10.0, cmp="ge"))
    for t, lab in (("L1", "LOCAL x1"), ("L3", "LOCAL x1/25 physics kept (u_f/h = 10)"),
                   ("I1", "INTEGRAL x1"), ("I3", "INTEGRAL x1/25 ell/k"), ("I4", "INTEGRAL x1/25 ell kept")):
        ck.append(Check("%s: W/(G_f A)" % lab, r[t]["w"], None, cmp="info"))
        ck.append(Check("%s: Delta_end/L (F_end/(sigma_y A) = %.3f)" % (lab, r[t]["fend"]), r[t]["uend"], None,
                        cmp="info"))
    return ck
