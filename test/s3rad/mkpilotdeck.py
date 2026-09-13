#!/usr/bin/env python3
"""Derive a reduced-output copy of the s3rad target deck.

Why a copy exists at all
------------------------
The committed deck asks for `*El Print, Elset=INTERFACE` and
`*Node Print` at the default frequency of one, i.e. SDV, S and E for
5400 UC6 elements on EVERY increment.  That is ~1.2 MB of .dat per
increment; a continuation run of several thousand increments would write
tens of gigabytes and stop on disk rather than on mechanics.

The ONLY difference this script introduces is output frequency.  Nothing
about the mesh, the materials, the interface, the boundary conditions,
the solver or the step control is touched, and the script prints the
exact diff and both SHA-256 hashes so the change is auditable.

    ./mkpilotdeck.py -o /path/to/pilot.inp [--every 200]
"""

import argparse
import difflib
import hashlib
import os
import re
import sys

SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   'm12_s3rad_gc24_w.inp')
SRC_SHA = ('2fb0cf4e3554282e1f85cf641c788939bc7a2fa5918dd842f54d38844df'
           'a5391')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-o', '--out', required=True)
    ap.add_argument('--every', type=int, default=200,
                    help='output frequency for the print and file cards')
    a = ap.parse_args()

    src = open(SRC, 'rb').read()
    got = hashlib.sha256(src).hexdigest()
    if got != SRC_SHA:
        print('source deck hash mismatch\n expected %s\n actual   %s'
              % (SRC_SHA, got), file=sys.stderr)
        return 2

    lines = src.decode('ascii').split('\n')
    out = []
    for l in lines:
        s = l.strip().lower()
        if s.startswith('*output'):
            l = '*Output, Frequency=%d' % a.every
        elif s.startswith('*node print') or s.startswith('*el print'):
            if 'frequency' not in s:
                l = l.rstrip() + ', Frequency=%d' % a.every
            else:
                l = re.sub(r'(?i)frequency\s*=\s*\d+',
                           'Frequency=%d' % a.every, l)
        out.append(l)
    txt = '\n'.join(out)
    open(a.out, 'w', newline='\n').write(txt)

    print('source     %s\n           sha256 %s' % (SRC, got))
    print('derived    %s\n           sha256 %s'
          % (a.out, hashlib.sha256(txt.encode('ascii')).hexdigest()))
    print('\nexact difference:')
    for d in difflib.unified_diff(lines, out, 'committed', 'pilot',
                                  n=0, lineterm=''):
        if d.startswith(('---', '+++', '@@', '-', '+')):
            print('  ' + d)
    return 0


if __name__ == '__main__':
    sys.exit(main())
