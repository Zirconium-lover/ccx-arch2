#!/usr/bin/env python3
"""Plot the s3rad A/B: legacy backtracking ladder against the fixed one.

Everything plotted is read from the two run logs, which come from the SAME
binary with one flag between them.

Usage: plot_ab.py <A-rundir> <B-rundir> -o out.png
"""

import argparse
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

C_A, C_B = '#b04a3a', '#2e5c8a'


def read(path):
    inc, theta, act, dele, batch = [], [], [], [], []
    ndel = 0
    for l in open(path):
        m = re.search(r'\[DAMAGE DE1.2 COMMIT\] inc=(\d+) time=([\d.eE+-]+) '
                      r'active=(\d+)', l)
        if m:
            inc.append(int(m.group(1)))
            theta.append(float(m.group(2)))
            act.append(int(m.group(3)))
            dele.append(ndel)
            continue
        m = re.search(r'\[DAMAGE COMMIT\] batch=(\d+) inc=(\d+) .* '
                      r'deleted=(\d+)', l)
        if m:
            ndel += int(m.group(3))
            batch.append((int(m.group(2)), int(m.group(1))))
    return dict(inc=inc, theta=theta, act=act, dele=dele, batch=batch)


def style(ax, xl, yl, t):
    ax.set_xlabel(xl); ax.set_ylabel(yl); ax.set_title(t, fontsize=10)
    ax.grid(alpha=.25, lw=.6)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.legend(fontsize=8, frameon=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('a'); ap.add_argument('b')
    ap.add_argument('-o', '--out', default='ab.png')
    p = ap.parse_args()
    A = read(p.a + '/run.log'); B = read(p.b + '/run.log')

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    for d, lab, c in ((A, 'legacy ladder {1.0, 0.5, 0.1}', C_A),
                      (B, 'fixed ladder, floor 1e-3, best rung', C_B)):
        ax[0].plot(d['inc'], d['theta'], color=c, lw=1.3, label=lab)
        ax[1].plot(d['inc'], d['dele'], color=c, lw=1.3, label=lab)
        ax[2].plot(d['inc'], d['act'], color=c, lw=1.3, label=lab)
    if A['theta']:
        for a in ax:
            a.axhline(0, lw=0)
        ax[0].axhline(A['theta'][-1], color=C_A, lw=.8, ls=':')
        ax[0].annotate('the recorded wall', (A['inc'][-1], A['theta'][-1]),
                       textcoords='offset points', xytext=(-90, 8),
                       fontsize=8, color=C_A)
    style(ax[0], 'accepted increment', 'load factor theta',
          's3rad, DEADALL=1.e-2: the wall is passed')
    style(ax[1], 'accepted increment', 'elements deleted (cumulative)',
          'front advance')
    style(ax[2], 'accepted increment', 'live bulk elements',
          'material remaining')
    fig.tight_layout(); fig.savefig(p.out, dpi=140)
    print('wrote', p.out)


if __name__ == '__main__':
    main()
