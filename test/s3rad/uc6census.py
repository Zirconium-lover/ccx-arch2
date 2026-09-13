#!/usr/bin/env python3
"""Count the UC6 interface state from a CalculiX .dat SDV block.

The s3rad deck prints `SDV` for the whole INTERFACE element set, and for
a UC6 point the four state variables are

    1 dmax   maximum effective separation reached
    2 dvisc  viscous damage
    3 dback  inviscid backbone damage
    4 failed 0/1

(resultsmech_uc6.f).  With d0 = Tn0/Kn and df = 2*Gc/Tn0 per *USER
SECTION, that is enough to count initiated points, the process zone and
the failed points WITHOUT re-running anything - which is what makes a
stock run and a continuation run comparable on the same yardstick.

The deck has two interface materials, so the thresholds differ per
element set; --seed gives the SEED element range.
"""

import argparse
import re
import sys

# *User Section values from the committed deck
SEED = dict(kn=200000.0, tn0=275.0, ts0=240.0, gc=2.000)
REG = dict(kn=200000.0, tn0=300.0, ts0=260.0, gc=2.400)


def thresholds(p):
    return p['tn0'] / p['kn'], 2.0 * p['gc'] / p['tn0']


def elset(path, name):
    """Element numbers in *Elset, Elset=<name>."""
    out, on = set(), False
    for l in open(path):
        s = l.strip()
        if s.startswith('*'):
            on = s.lower().replace(' ', '').startswith('*elset,elset=' +
                                                       name.lower())
            continue
        if on and s and not s.startswith('**'):
            for t in s.split(','):
                t = t.strip()
                if t:
                    out.add(int(t))
    return out


def blocks(path):
    """Yield (increment, {elem: [dmax,...]}) for each SDV block."""
    inc, cur, on = None, None, False
    with open(path) as f:
        for l in f:
            m = re.match(r'\s*INCREMENT\s+(\d+)\s*$', l)
            if m:
                if cur is not None:
                    yield inc, cur
                inc, cur, on = int(m.group(1)), None, False
                continue
            if 'internal state variables' in l:
                cur, on = {}, True
                continue
            if on:
                v = l.split()
                if len(v) >= 6 and v[0].isdigit():
                    cur.setdefault(int(v[0]), []).append(
                        [float(x) for x in v[2:6]])
                elif l.strip() == '':
                    continue
                else:
                    on = False
    if cur is not None:
        yield inc, cur


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('dat')
    ap.add_argument('--inp', required=True)
    ap.add_argument('--last', type=int, default=1,
                    help='report only the last N blocks')
    a = ap.parse_args()

    seed = elset(a.inp, 'interface_seed')
    d0s, dfs = thresholds(SEED)
    d0r, dfr = thresholds(REG)
    keep = []
    for inc, b in blocks(a.dat):
        keep.append((inc, b))
        if len(keep) > a.last:
            keep.pop(0)
    if not keep:
        print('no SDV block found', file=sys.stderr)
        return 2
    print("SEED d0=%.6e df=%.6e ; REGULAR d0=%.6e df=%.6e"
          % (d0s, dfs, d0r, dfr))
    print("  inc     nip   init    zone    fail   dmax_max      "
          "dmax_mean_zone")
    for inc, b in keep:
        nip = ninit = nzone = nfail = 0
        mx, ssum = 0.0, 0.0
        for e, rows in b.items():
            d0, df = (d0s, dfs) if e in seed else (d0r, dfr)
            for r in rows:
                dmax = r[0]
                nip += 1
                mx = max(mx, dmax)
                if dmax > d0:
                    ninit += 1
                    if dmax < df:
                        nzone += 1
                        ssum += dmax
                    else:
                        nfail += 1
        print(" %5s %7d %6d %7d %7d  %.6e  %s"
              % (inc, nip, ninit, nzone, nfail, mx,
                 ("%.6e" % (ssum / nzone)) if nzone else "-"))
    return 0


if __name__ == '__main__':
    sys.exit(main())
