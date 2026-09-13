#!/usr/bin/env python3
"""Plot what a crack-control run actually did.

Three figures, one per kind of evidence:

  --modeI   <dir>   lambda and reaction against the controlled opening,
                    with the closed-form bilinear branch overlaid
  --mixed   <dir>   the same for the inclined-interface benchmark, with
                    the finite-strain closed form
  --s3rad   <dirs>  load factor, front census and BOTH residuals against
                    accepted increment, stock run against continuation

Output is a .png; nothing here is committed (see .gitignore).
"""

import argparse
import math
import os
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

C_STOCK, C_CTRL, C_EXACT = '#b04a3a', '#2e5c8a', '#3f3f3f'


def read_dat(path):
    F, U = [], []
    lines = open(path).read().split('\n')
    i = 0
    while i < len(lines):
        if 'total force' in lines[i]:
            F.append(float(lines[i + 2].split()[0])); i += 3
        elif 'displacements' in lines[i]:
            U.append(float(lines[i + 2].split()[1])); i += 3
        else:
            i += 1
    return F, U


def read_log(path, tag):
    out = []
    rx = re.compile(r'\[' + tag + r'\] inc=(\d+) ACCEPTED lambda=([\d.eE+-]+)')
    for l in open(path):
        m = rx.search(l)
        if not m:
            continue
        d = dict(inc=int(m.group(1)), lam=float(m.group(2)))
        for k, r in [('phi', r'phi=([\d.eE+-]+)'),
                     ('ach', r'achieved=([\d.eE+-]+)'),
                     ('g', r' g=([\d.eE+-]+)'), ('R', r'\|R\|=([\d.eE+-]+)'),
                     ('zone', r'zone=(\d+)'), ('load', r'load=(\d+)'),
                     ('fail', r'fail=(\d+)'), ('dead', r'dead=(\d+)'),
                     ('deff', r'deffmax=([\d.eE+-]+)'),
                     ('shear', r'shear=([\d.eE+-]+)')]:
            mm = re.search(r, l)
            if mm:
                d[k] = float(mm.group(1))
        out.append(d)
    return out


def style(ax, xl, yl, t=None):
    ax.set_xlabel(xl); ax.set_ylabel(yl)
    if t:
        ax.set_title(t, fontsize=10)
    ax.grid(alpha=.25, lw=.6)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)


def fig_bar(stock, ctrl, out, title, exact=None):
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for d, lab, c in ((stock, 'stock displacement control', C_STOCK),
                      (ctrl, 'crack control', C_CTRL)):
        if d is None:
            continue
        F, U = d
        ax[0].plot(U, F, color=c, lw=1.2, label=lab)
    if exact:
        ax[0].plot(exact[0], exact[1], color=C_EXACT, lw=1.0, ls='--',
                   label='closed form')
    style(ax[0], 'prescribed end displacement', 'reaction',
          title + ' - the branch')
    ax[0].legend(fontsize=8, frameon=False)

    for d, lab, c in ((stock, 'stock', C_STOCK), (ctrl, 'crack control', C_CTRL)):
        if d is None:
            continue
        F, U = d
        ax[1].plot(range(1, len(U) + 1), U, color=c, lw=1.2, label=lab)
    style(ax[1], 'accepted increment', 'prescribed end displacement',
          'the controlled variable must run backwards')
    ax[1].legend(fontsize=8, frameon=False)
    fig.tight_layout(); fig.savefig(out, dpi=140)
    print('wrote', out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--modeI', nargs=2, metavar=('STOCK', 'CTRL'))
    ap.add_argument('--mixed', nargs=2, metavar=('STOCK', 'CTRL'))
    ap.add_argument('--s3rad', nargs=1, metavar='RUNDIR',
                    help='one pilot run; it is split at the engagement '
                         'increment, because the armed-but-not-engaged '
                         'phase reproduces the stock trajectory '
                         'line-for-line')
    ap.add_argument('-o', '--out', default='out.png')
    a = ap.parse_args()

    if a.modeI:
        s = read_dat(os.path.join(a.modeI[0], 'cohesive.dat'))
        c = read_dat(os.path.join(a.modeI[1], 'cohesive.dat'))
        Tn0, d0, df, Kn = 400., 4e-4, 1e-2, 1e6
        xs = [d0 * i / 200. for i in range(201)] + \
             [d0 + (df - d0) * i / 400. for i in range(401)]
        ys = [Kn * x if x <= d0 else Tn0 * (df - x) / (df - d0) for x in xs]
        us = [x + y * 20. / 200000. for x, y in zip(xs, ys)]
        fig_bar(s, c, a.out, 'Mode I', exact=(us, ys))

    if a.mixed:
        s = read_dat(os.path.join(a.mixed[0], 'mixed.dat'))
        c = read_dat(os.path.join(a.mixed[1], 'mixed.dat'))
        Tn0, Ts0, Kn, Gc, E, L = 400., 300., 1e6, 2., 2e5, 20.
        d0, df = Tn0 / Kn, 2 * Gc / Tn0
        cx = 1. / math.sqrt(2.); beta = (Ts0 / Tn0) ** 2
        kap = math.sqrt(cx * cx + beta * (1 - cx * cx))

        def sg(d):
            return (Kn * kap * d / cx) if d <= d0 else \
                (kap / cx) * Tn0 * (df - d) / (df - d0)

        def st(s_):
            x = s_ / E
            for _ in range(60):
                x -= (E * (x + 1.5 * x * x + .5 * x ** 3) - s_) / \
                     (E * (1. + 3. * x + 1.5 * x * x))
            return x
        xs = [d0 * i / 200. for i in range(201)] + \
             [d0 + (df - d0) * i / 400. for i in range(401)]
        ys = [sg(x) for x in xs]
        us = [x / kap + st(y) * L for x, y in zip(xs, ys)]
        fig_bar(s, c, a.out, 'mixed mode, 45 deg / 36% shear',
                exact=(us, ys))

    if a.s3rad:
        d = read_log(os.path.join(a.s3rad[0], 'run.log'), 'CRACKCTL')
        eng = None
        for l in open(os.path.join(a.s3rad[0], 'run.log')):
            m = re.search(r'\[CRACKCTL\] engaged at inc=(\d+)', l)
            if m:
                eng = int(m.group(1))
        if eng is None:
            eng = d[-1]['inc'] + 1
        st = [x for x in d if x['inc'] <= eng]
        ct = [x for x in d if x['inc'] >= eng]
        fig, ax = plt.subplots(1, 4, figsize=(19, 4.3))
        for dd, lab, c in ((st, 'ordinary control (= stock)', C_STOCK),
                           (ct, 'crack control', C_CTRL)):
            if not dd:
                continue
            ax[0].plot([x['inc'] for x in dd], [x['lam'] for x in dd],
                       color=c, lw=1.2, label=lab)
            ax[1].plot([x['inc'] for x in dd],
                       [x.get('deff', 0) for x in dd], color=c, lw=1.2,
                       label=lab)
            ax[2].semilogy([x['inc'] for x in dd],
                           [max(abs(x.get('R', 0.)), 1e-16) for x in dd],
                           color=c, lw=.9, label=lab + r'  $\|R\|_\infty$')
            ax[2].semilogy([x['inc'] for x in dd],
                           [max(abs(x.get('g', 0.)), 1e-24) for x in dd],
                           color=c, lw=.9, ls='--', label=lab + '  |g|')
        # fourth panel: the engaged stretch on its own, because on the
        # scale of the whole run it is a few pixels wide
        if len(ct) > 1:
            ax[3].plot([x['inc'] for x in ct], [x['lam'] for x in ct],
                       color=C_CTRL, lw=1.4, marker='o', ms=3,
                       label='load factor')
            a2 = ax[3].twinx()
            a2.plot([x['inc'] for x in ct], [x.get('deff', 0) for x in ct],
                    color=C_STOCK, lw=1.4, marker='s', ms=3,
                    label='largest deff')
            a2.set_ylabel('largest UC6 deff', color=C_STOCK)
            a2.tick_params(axis='y', colors=C_STOCK)
            for sp in ('top',):
                a2.spines[sp].set_visible(False)
            ax[3].legend(fontsize=7, frameon=False, loc='upper left')
        for x in ax[:3]:
            x.axvline(eng, color='#888', lw=.8, ls=':')
        style(ax[0], 'accepted increment', 'load factor',
              's3rad load factor (dotted line = engagement)')
        style(ax[1], 'accepted increment', 'largest UC6 deff',
              's3rad front: largest effective separation')
        style(ax[2], 'accepted increment', 'residual',
              's3rad: equilibrium and constraint residuals')
        style(ax[3], 'accepted increment', 'load factor',
              'the engaged stretch: lambda turns, deff advances')
        for x in ax[:3]:
            x.legend(fontsize=7, frameon=False)
        fig.tight_layout(); fig.savefig(a.out, dpi=140)
        print('wrote', a.out)


if __name__ == '__main__':
    main()
