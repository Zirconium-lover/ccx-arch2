#!/usr/bin/env python3
"""Replace the material/interface/step cards of an mkseeddeck.py deck with the
E635 property set ("props v1"), leaving mesh, sets and boundary sets untouched.

    mkprops.py <in.inp> <out.inp>

Sources (see props_v1 notes in the session):
  ZR   E 100 GPa (Kim&Lee 2024 tab.2: Zr 95.5-103, Zry 99-115)
       sigma_y 490 MPa: Fedotov et al., fig.35, E635 ring samples, 293 K, unirradiated
       hardening sigma/sigma_y: Fedotov et al., fig.33, E635 at 293 K
       eps_f(eta): shape of the s3rad deck, scaled so that eps_f(1/3) = 0.97 =
       -ln(1-0.62), dS/S of E635 rings with Fn < 0.25 at 20 C (Egorochkin, tab.3.1)
  ZRH  E 130 GPa (Kim&Lee tab.2, DFT/ultrasonic cluster), sigma_y 1160 MPa
       (Chan et al. 2016 tab.1), eps_f ~1 % (cleavage after slight yield);
       u_f = 2 G_c / sigma_y with G_c = K_C^2 (1-nu^2)/E, K_C = 3.3 MPa m^0.5
  INTERFACE  G_c = (2.15/2.8)^2 G_c(ZrH) (Chan fig.3 fracture-load ratio)
  Fracture energies are multiplied by SCALE = model plate radius / real hydride
  half-length (0.3 mm / 0.03 mm): the LEFM similitude that keeps the plate as
  brittle as a ~60 um hydride (Plyasov et al. 2023: radial hydrides > 50 um).
  Step: grip displacement 2.0 mm (50 % nominal) with the SAME displacement
  increments as before (dtheta halved).  Field output every 10 increments.
"""
import sys

SCALE = 10.0
KC = 3.3 * 1000 ** 0.5                     # MPa sqrt(mm)
E_H, NU_H, SY_H = 130000., 0.33, 1160.
G_H = KC ** 2 * (1 - NU_H ** 2) / E_H      # N/mm, ~0.0746
G_I = (2.15 / 2.8) ** 2 * G_H              # ~0.044
G_H_M, G_I_M = SCALE * G_H, SCALE * G_I
UF_H = 2. * G_H_M / SY_H

SY_Z = 490.
HARD = [(0., 1.0), (0.02, 1.15), (0.05, 1.27), (0.10, 1.38), (0.20, 1.52),
        (0.30, 1.60), (0.50, 1.80), (0.75, 2.05), (1.00, 2.45)]
ETA = [0.00, 0.33, 0.50, 0.67, 1.00, 1.50]
EF_OLD_Z = [0.30, 0.12, 0.07, 0.0425, 0.025, 0.015]
EF_OLD_H = [0.120, 0.070, 0.045, 0.028, 0.015, 0.010]
EF_Z = [e * 0.97 / 0.12 for e in EF_OLD_Z]
EF_H = [e * 0.010 / 0.070 for e in EF_OLD_H]

TN, TS, TN_S, TS_S = 300.0, 260.0, 275.0, 240.0
KN, GMIN, VISC = 1.0e6, 2.0e-6, 1.e-4
UX = 2.0
FREQ = 10


def locus(uf, ef):
    pairs = ['%.2f, %.4g' % (a, b) for a, b in zip(ETA, ef)]
    return ('1.0, %.4g, ' % uf + ', '.join(pairs[:3]) + '\n' +
            ', '.join(pairs[3:]))


def tail():
    L = ['*Material, Name=ZR', '*Elastic', '100000., 0.335', '*Plastic']
    L += ['%.1f, %.3f' % (SY_Z * f, e) for e, f in HARD]
    L += ['*Damage Initiation, Criterion=Ductile, Evolution=Displacement, Npoints=6',
          locus(0.015, EF_Z)]
    L += ['*Material, Name=ZRH', '*Elastic', '%.0f., %.2f' % (E_H, NU_H), '*Plastic',
          '%.1f, 0.' % SY_H, '%.1f, 0.0050' % (SY_H + 20), '%.1f, 0.0200' % (SY_H + 40),
          '*Damage Initiation, Criterion=Ductile, Evolution=Displacement, Npoints=6',
          locus(UF_H, EF_H)]
    L += ['*Material, Name=COHESIVE', '*Elastic', '1., 0.', '*Depvar', '4',
          '*Solid Section, Elset=MATRIX, Material=ZR',
          '*Solid Section, Elset=PLATETANGENTIAL, Material=ZRH',
          '** Kn, Tn0, Ts0, Gc, residual fraction, viscosity',
          '*User Section, Elset=INTERFACE_SEED, Material=COHESIVE, Constants=6',
          '%.6e, %.1f, %.1f, %.4f, %.6e, %.0e' % (KN, TN_S, TS_S, G_I_M * 2.0 / 2.4, GMIN, VISC),
          '*User Section, Elset=INTERFACE_REGULAR, Material=COHESIVE, Constants=6',
          '%.6e, %.1f, %.1f, %.4f, %.6e, %.0e' % (KN, TN, TS, G_I_M, GMIN, VISC),
          '*Step, Nlgeom, Inc=200000', '*Static, Solver=Pardiso',
          '5.000000e-04, 1., 5.000000e-10, 1.000000e-03',
          '*Boundary', 'FACE_X0_NSET, 1, 1, 0.', 'FACE_XL_NSET, 1, 1, %.3f' % UX,
          'FIXPOINTA, 2, 3, 0.', 'FIXPOINTB, 3, 3, 0.',
          '*Output, Frequency=%d' % FREQ, '*Node File', 'U, RF', '*El File',
          'S, E, PEEQ, SDV', '*Node Print, Nset=FACE_XL_NSET, Totals=Yes', 'RF, U',
          '*El Print, Elset=INTERFACE, Frequency=%d' % FREQ, 'S, E, SDV', '*End Step']
    return '\n'.join(L) + '\n'


def main():
    src, out = sys.argv[1], sys.argv[2]
    lines = open(src).read().split('\n')
    k = lines.index('*Material, Name=ZR')
    head = ['** props v1 (E635: Fedotov, Kim&Lee, Chan, Egorochkin; energies x%g)' % SCALE]
    open(out, 'w', newline='\n').write('\n'.join(head + lines[:k]) + '\n' + tail())
    print('%s -> %s  G_H=%.4f G_I=%.4f (x%g: %.3f / %.3f N/mm)  uf_H=%.3e  '
          'delta_f/delta_0 = %.1f' % (src, out, G_H, G_I, SCALE, G_H_M, G_I_M, UF_H,
                                      (2 * G_I_M / TN) / (TN / KN)))


if __name__ == '__main__':
    main()
