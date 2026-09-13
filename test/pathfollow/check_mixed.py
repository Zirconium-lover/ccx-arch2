#!/usr/bin/env python3
"""Verify a mixed.inp run against the closed-form branch.

What is checked
---------------
1. UC6 kinematics: the printed local separations (dn,ds1,ds2) reproduce
   deff = sqrt(max(dn,0)^2 + beta*(ds1^2+ds2^2)) and the shear fraction
   the deck header predicts.  The tolerance is 1e-5 because the .dat
   file carries seven significant digits, not because the identity is
   approximate.  This is what makes the benchmark a
   MIXED-MODE test rather than a rescaled Mode-I one.
2. Constitutive branch: the reaction against
       sigma(deff) = (kappa/c)*Tn0*(df-deff)/(df-d0)
   which involves no assumption about the bulk at all.
3. Kinematic branch: the prescribed end displacement against
       u = D + a*L ,   D = deff/kappa ,   (1+a)*E*(a+a^2/2) = sigma
   i.e. the FINITE-STRAIN bulk response.  CalculiX runs UC6 only on the
   NLGEOM path (user elements are rejected on the linear path), so the
   bulk is St-Venant-Kirchhoff and the small-strain form u = D + sigma*L/E
   is wrong by O(eps^2) ~ 4e-3 relative - which is exactly the residual
   the first version of this script showed.
4. Snap-back: the end displacement must run BACKWARDS.
5. Continuation: how many consecutive post-peak increments were accepted
   and the largest constraint residual reported at an accepted state.

Usage:  check_mixed.py <rundir> [--jobname mixed]
"""

import argparse
import math
import os
import re
import sys


def parse_dat(path):
    """Return per-increment (F, u, [(elem,ip,dn,ds1,ds2,dvisc,deff,dmax)])."""
    inc, out = None, {}
    lines = open(path).read().split('\n')
    i = 0
    while i < len(lines):
        s = lines[i]
        m = re.match(r'\s*INCREMENT\s+(\d+)\s*$', s)
        if m:
            inc = int(m.group(1))
            out[inc] = {'F': None, 'u': None, 'e': []}
            i += 1
            continue
        if inc is None:
            i += 1
            continue
        if 'total force' in s:
            out[inc]['F'] = float(lines[i + 2].split()[0])
            i += 3
            continue
        if 'displacements' in s:
            out[inc]['u'] = float(lines[i + 2].split()[1])
            i += 3
            continue
        if s.lstrip().startswith('strains (elem'):
            i += 2
            while i < len(lines) and lines[i].strip() and \
                    not lines[i].lstrip().startswith(('INCREMENT', 'total',
                                                      'displacements',
                                                      'internal', 'strains')):
                v = lines[i].split()
                if len(v) >= 8:
                    out[inc]['e'].append(tuple(float(x) for x in v[2:8]))
                i += 1
            continue
        i += 1
    return out


def hdr(path):
    """Read the parameters the generator wrote into the deck header."""
    p = {}
    for l in open(path):
        if not l.startswith('**'):
            break
        for k, rx in [('E', r'E=([\d.eE+-]+)'), ('Kn', r'Kn=([\d.eE+-]+)'),
                      ('Tn0', r'Tn0=([\d.eE+-]+)'), ('Ts0', r'Ts0=([\d.eE+-]+)'),
                      ('Gc', r'Gc=([\d.eE+-]+)'), ('gmin', r'gmin=([\d.eE+-]+)'),
                      ('L', r'L=([\d.eE+-]+)'), ('c', r'n_x=([\d.eE+-]+)'),
                      ('kap2', r'kappa\^2=([\d.eE+-]+)'),
                      ('shear', r'shear carries ([\d.eE+-]+)')]:
            m = re.search(rx, l)
            if m and k not in p:
                p[k] = float(m.group(1))
    return p


def bulk_stretch(sig, E):
    """(1+a)*E*(a+a^2/2) = sigma, solved by Newton.  St-Venant-Kirchhoff
    uniaxial with the lateral strains held at zero (nu=0, u_y=u_z=0), so
    J = 1+a and the first Piola stress is (1+a)*S."""
    a = sig / E
    for _ in range(60):
        f = E * (a + 1.5 * a * a + 0.5 * a ** 3) - sig
        d = E * (1. + 3. * a + 1.5 * a * a)
        a -= f / d
    return a


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('rundir')
    ap.add_argument('--jobname', default='mixed')
    ap.add_argument('--ftol', type=float, default=1e-4)
    ap.add_argument('--utol', type=float, default=1e-5)
    ap.add_argument('--minpost', type=int, default=20)
    a = ap.parse_args()

    dat = os.path.join(a.rundir, a.jobname + '.dat')
    inp = os.path.join(a.rundir, a.jobname + '.inp')
    log = os.path.join(a.rundir, 'run.log')
    p = hdr(inp)
    E, Kn, Tn0, Ts0 = p['E'], p['Kn'], p['Tn0'], p['Ts0']
    Gc, gmin, L, c = p['Gc'], p['gmin'], p['L'], p['c']
    beta = (Ts0 / Tn0) ** 2
    kap = math.sqrt(p['kap2'])
    d0, df = Tn0 / Kn, 2. * Gc / Tn0
    rec = parse_dat(dat)
    bad = 0

    # ---- 1. kinematics ------------------------------------------------
    worst_k, worst_s = 0., None
    for inc, r in rec.items():
        for (dn, ds1, ds2, dvisc, deff, dmax) in r['e']:
            de = math.sqrt(max(dn, 0.) ** 2 + beta * (ds1 ** 2 + ds2 ** 2))
            worst_k = max(worst_k, abs(de - deff) / max(deff, 1e-30))
            if deff > d0:
                sf = beta * (ds1 ** 2 + ds2 ** 2) / deff ** 2
                d = abs(sf - p['shear'])
                worst_s = d if worst_s is None else max(worst_s, d)
    print("1 UC6 kinematics      max rel |deff_recomputed-deff| = %.3e  %s"
          % (worst_k, "ok" if worst_k < 1e-5 else "FAIL"))
    bad += worst_k >= 1e-5
    if worst_s is not None:
        print("  shear fraction      predicted %.4f, max deviation %.3e  %s"
              % (p['shear'], worst_s, "ok" if worst_s < 1e-5 else "FAIL"))
        bad += worst_s >= 1e-5

    # ---- 2/3. branch --------------------------------------------------
    def sigma(deff):
        g = 1. if deff <= d0 else (gmin if deff >= df else
                                   max(gmin, 1. - df * (deff - d0) /
                                       (deff * (df - d0))))
        return g * Kn * kap * deff / c

    ef, eu, npost, us = 0., 0., 0, []
    efwhere, euwhere = None, None
    for inc in sorted(rec):
        r = rec[inc]
        if (r['F'] is None) or (r['u'] is None) or (not r['e']):
            continue
        deff = r['e'][0][4]
        us.append(r['u'])
        if not (d0 < deff < df):
            continue
        npost += 1
        sg = sigma(deff)
        e = abs(r['F'] - sg) / max(abs(sg), 1e-30)
        if e > ef:
            ef, efwhere = e, (inc, deff)
        uex = deff / kap + bulk_stretch(sg, E) * L
        e = abs(r['u'] - uex) / max(abs(uex), 1e-30)
        if e > eu:
            eu, euwhere = e, (inc, deff)
    print("2 reaction branch     max rel error vs closed form = %.3e at "
          "inc %s (deff=%.6e)  %s"
          % (ef, efwhere[0] if efwhere else '-',
             efwhere[1] if efwhere else 0., "ok" if ef < a.ftol else "FAIL"))
    bad += ef >= a.ftol
    print("3 displacement branch max rel error vs closed form = %.3e at "
          "inc %s (deff=%.6e)  %s"
          % (eu, euwhere[0] if euwhere else '-',
             euwhere[1] if euwhere else 0., "ok" if eu < a.utol else "FAIL"))
    bad += eu >= a.utol

    # ---- 4. snap-back --------------------------------------------------
    back = min(us[i + 1] - us[i] for i in range(len(us) - 1)) if len(us) > 1 else 0.
    print("4 snap-back           u max=%.8f last=%.8f largest backward step "
          "%.3e  %s" % (max(us), us[-1], back, "ok" if back < 0. else "FAIL"))
    bad += not (back < 0.)

    # ---- 5. continuation ------------------------------------------------
    gmax, nacc = 0., 0
    if os.path.exists(log):
        for l in open(log):
            m = re.search(r'\[CRACKCTL\] inc=\d+ ACCEPTED .* g=([\d.eE+-]+)', l)
            if m:
                nacc += 1
                gmax = max(gmax, abs(float(m.group(1))))
    print("5 continuation        %d accepted increments, %d of them post-peak, "
          "max |g| at an accepted state = %.3e" % (nacc, npost, gmax))
    print("                      post-peak requirement >= %d  %s"
          % (a.minpost, "ok" if npost >= a.minpost else "FAIL"))
    bad += npost < a.minpost

    print("\n%s (%d failure(s))" % ("PASSED" if bad == 0 else "FAILED", bad))
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
