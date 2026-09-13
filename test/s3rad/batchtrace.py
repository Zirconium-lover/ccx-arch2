#!/usr/bin/env python3
"""Answer the repetition question from a CCX_DAMAGE_BATCH_TRACE log.

The claim under test is:

    a deletion batch forms, the same-load re-equilibration fails, the
    batch is rolled back, and the run returns to the SAME committed state
    and forms the SAME batch again.

Both halves are needed.  A committed state that repeats is ordinary - every
retry of an increment starts from it - and a batch that repeats across
different states is ordinary too.  What would prove the claim is the PAIR
(committed_state, batch) occurring more than once.

Usage: batchtrace.py <run.log> [--top N]
"""

import argparse
import collections
import re
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('log')
    ap.add_argument('--top', type=int, default=8)
    a = ap.parse_args()

    att = re.compile(r'\[BATCHTRACE\] attempt inc=(\d+) icutb=(\d+) '
                     r'theta=([\d.eE+-]+) dtheta=([\d.eE+-]+) '
                     r'committed_state=([0-9a-f]+)')
    bat = re.compile(r'\[BATCHTRACE\] batch inc=(\d+) icutb=(\d+) pass=(\d+) '
                     r'n=(\d+) committed_state=([0-9a-f]+) '
                     r'batch=([0-9a-f]+) sorted:(.*)')

    attempts, batches = [], []
    for l in open(a.log):
        m = att.search(l)
        if m:
            attempts.append(dict(inc=int(m.group(1)), icutb=int(m.group(2)),
                                 theta=float(m.group(3)),
                                 dtheta=float(m.group(4)), st=m.group(5)))
            continue
        m = bat.search(l)
        if m:
            batches.append(dict(inc=int(m.group(1)), icutb=int(m.group(2)),
                                npass=int(m.group(3)), n=int(m.group(4)),
                                st=m.group(5), b=m.group(6),
                                elems=m.group(7).split()))

    print("attempts traced : %d" % len(attempts))
    print("batches traced  : %d" % len(batches))
    if not attempts:
        print("no trace found - was CCX_DAMAGE_BATCH_TRACE=1 set?",
              file=sys.stderr)
        return 2

    st_att = collections.Counter(x['st'] for x in attempts)
    print("\ndistinct committed states seen at an attempt: %d" % len(st_att))
    print("most-retried committed states (state, attempts, increments):")
    for st, n in st_att.most_common(a.top):
        incs = sorted({x['inc'] for x in attempts if x['st'] == st})
        print("  %s  %3d attempt(s)  increments %s"
              % (st, n, incs if len(incs) < 8 else
                 "%s..%s (%d)" % (incs[0], incs[-1], len(incs))))

    if batches:
        b_c = collections.Counter(x['b'] for x in batches)
        print("\ndistinct batch compositions: %d" % len(b_c))
        print("most-repeated batch compositions:")
        for b, n in b_c.most_common(a.top):
            r = [x for x in batches if x['b'] == b]
            print("  %s  %2d time(s)  n=%d  increments %s"
                  % (b, n, r[0]['n'], sorted({x['inc'] for x in r})))
            if n > 1:
                print("      elements: %s" % " ".join(r[0]['elems'][:20]))

        pair = collections.Counter((x['st'], x['b']) for x in batches)
        rep = [(k, v) for k, v in pair.items() if v > 1]
        print("\n*** THE CLAIM: (same committed state, same batch) pairs "
              "occurring more than once: %d" % len(rep))
        for (st, b), n in sorted(rep, key=lambda kv: -kv[1])[:a.top]:
            r = [x for x in batches if x['st'] == st and x['b'] == b]
            print("  state %s batch %s  x%d  increments %s  n=%d"
                  % (st, b, n, sorted({x['inc'] for x in r}), r[0]['n']))
            print("      elements: %s" % " ".join(r[0]['elems'][:20]))
        if not rep:
            print("  none - no batch was ever re-formed from a state it had "
                  "already been formed from, so the batch/rollback loop did "
                  "not occur in this run.")
    else:
        print("\nno deletion batch was traced at all in this run.")
    return 0


if __name__ == '__main__':
    sys.exit(main())
