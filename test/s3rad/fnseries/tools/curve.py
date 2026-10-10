#!/usr/bin/env python3
"""Stress-strain curve of a run and the series numbers, from m.dat only.

    curve.py <run_dir> [...]      -> <run_dir>/curve.csv and one summary line each

sigma = total RF_x on FACE_XL / initial cross-section (3.6 x 1.6 mm),
eps = u_x(grip) / 4.0.  The fracture strain is read where the force has fallen
to 10 % of its peak after the peak (eps_f10), so a run that stops before
FRACTURE COMPLETE but past that point still gives the number.
"""
import os, re, sys

A0, L0 = 3.6 * 1.6, 4.0


def grip(inp):
    """Prescribed x displacement of FACE_XL_NSET in the run's deck."""
    on = False
    for l in open(inp, errors='replace'):
        s = l.strip()
        if s.startswith('*'):
            on = s.lower().replace(' ', '').startswith('*boundary')
            continue
        v = [x.strip() for x in s.split(',')]
        if on and v and v[0].upper() == 'FACE_XL_NSET' and len(v) >= 4:
            return float(v[3])
    raise SystemExit('no FACE_XL_NSET displacement in ' + inp)


def read(dat, ux):
    eps, sig = [], []
    t_rf, rf = None, None
    mode = None
    for l in open(dat, errors='replace'):
        m = re.search(r'total force \(fx,fy,fz\) for set FACE_XL_NSET and time\s+([0-9.Ee+-]+)', l)
        if m:
            mode = 'rf'; t_rf = float(m.group(1)); continue
        if mode == 'rf':
            v = l.split()
            if len(v) >= 3:
                try:
                    rf = float(v[0]); mode = None
                    eps.append(t_rf * ux / L0); sig.append(rf / A0)
                except ValueError:
                    pass
    return eps, sig


def main():
    for d in sys.argv[1:]:
        if not os.path.isdir(d):
            continue
        eps, sig = read(d + '/m.dat', grip(d + '/m.inp'))
        if not eps:
            print(d, 'no data'); continue
        open(d + '/curve.csv', 'w').write('eps,sigma_MPa\n' + '\n'.join(
            '%.6f,%.3f' % p for p in zip(eps, sig)))
        k = max(range(len(sig)), key=sig.__getitem__)
        W = sum(0.5 * (sig[i] + sig[i - 1]) * (eps[i] - eps[i - 1]) for i in range(1, len(eps)))
        f10 = next((eps[i] for i in range(k, len(sig)) if sig[i] < 0.1 * sig[k]), None)
        print('%-28s peak %.0f MPa at eps %.2f %%  last eps %.2f %% sigma %.0f (%.0f %% of peak)  '
              'eps_f10 %s  W %.1f MJ/m3' % (d.split('/')[-1], sig[k], 100 * eps[k], 100 * eps[-1],
              sig[-1], 100 * sig[-1] / sig[k], '%.2f %%' % (100 * f10) if f10 else '-', W))


if __name__ == '__main__':
    main()
