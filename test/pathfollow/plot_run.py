#!/usr/bin/env python3
"""Plot reaction vs end displacement, and the lambda history, for one run.

Reads the CalculiX .dat produced by snapback.inp (which prints the total
reaction on NRIGHT and the displacement of NRIGHT) plus, optionally, the
solver log so that the path-following load factor can be drawn against the
increment number.

Usage:  plot_run.py --dat snapback.dat [--log run.log] [--label stock]
                    [--out curve.png]
"""

import argparse
import re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_dat(path):
    """Return [(time, reaction_x, u_x)] in increment order."""
    times, rx, ux = [], [], []
    cur_t = None
    mode = None
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s.startswith("total force"):
                m = re.search(r"time\s+([0-9.E+-]+)", s)
                cur_t = float(m.group(1)) if m else None
                mode = "rf"
                continue
            if s.startswith("displacements"):
                mode = "u"
                continue
            if s.startswith("internal state") or s.startswith("stresses"):
                mode = None
                continue
            if not s:
                continue
            p = s.split()
            try:
                vals = [float(x) for x in p]
            except ValueError:
                mode = None
                continue
            if mode == "rf" and len(vals) == 3 and cur_t is not None:
                times.append(cur_t)
                rx.append(vals[0])
                mode = None
            elif mode == "u" and len(vals) == 4:
                # one line per node; the bar is uniform so any node will do
                if len(ux) < len(times):
                    ux.append(vals[1])
    n = min(len(times), len(rx), len(ux))
    return list(zip(times[:n], rx[:n], ux[:n]))


def read_lambda(path):
    lam, inc = [], []
    pat = re.compile(r"inc=(\d+) ACCEPTED lambda=([0-9.]+)")
    with open(path) as f:
        for line in f:
            m = pat.search(line)
            if m:
                inc.append(int(m.group(1)))
                lam.append(float(m.group(2)))
    return inc, lam


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dat", action="append", required=True,
                    help="path[:label]")
    ap.add_argument("--log", action="append", default=[],
                    help="path[:label]")
    ap.add_argument("--out", default="curve.png")
    ap.add_argument("--upeak", type=float, default=None)
    ap.add_argument("--ufail", type=float, default=None)
    a = ap.parse_args()

    ncol = 2 if a.log else 1
    fig, ax = plt.subplots(1, ncol, figsize=(6.2 * ncol, 4.4))
    if ncol == 1:
        ax = [ax]

    for spec in a.dat:
        path, _, label = spec.partition(":")
        rows = read_dat(path)
        if not rows:
            print("no data in", path)
            continue
        u = [r[2] for r in rows]
        f = [r[1] for r in rows]
        ax[0].plot(u, f, marker="o", ms=2.5, lw=1.1, label=label or path)
        print("%-28s %4d points, u in [%.6f, %.6f], peak F = %.4f"
              % (label or path, len(rows), min(u), max(u), max(f)))

    if a.upeak is not None:
        ax[0].axvline(a.upeak, ls="--", lw=0.9, color="0.5")
        ax[0].annotate("analytic peak\nu=%.5f" % a.upeak, (a.upeak, 0),
                       textcoords="offset points", xytext=(4, 14),
                       fontsize=8, color="0.35")
    if a.ufail is not None:
        ax[0].axvline(a.ufail, ls=":", lw=0.9, color="0.5")
        ax[0].annotate("analytic failure\nu=%.5f" % a.ufail, (a.ufail, 0),
                       textcoords="offset points", xytext=(4, 14),
                       fontsize=8, color="0.35")

    ax[0].set_xlabel("end displacement  u  (control variable)")
    ax[0].set_ylabel("total reaction  F")
    ax[0].set_title("reaction vs controlled displacement")
    ax[0].grid(alpha=0.3)
    ax[0].legend(fontsize=8)

    if a.log:
        for spec in a.log:
            path, _, label = spec.partition(":")
            inc, lam = read_lambda(path)
            if inc:
                ax[1].plot(inc, lam, lw=1.2, label=label or path)
                print("%-28s lambda in [%.6f, %.6f], last %.6f"
                      % (label or path, min(lam), max(lam), lam[-1]))
        ax[1].set_xlabel("accepted increment")
        ax[1].set_ylabel(r"load factor  $\lambda$")
        ax[1].set_title(r"$\lambda$ history (must decrease to snap back)")
        ax[1].grid(alpha=0.3)
        ax[1].legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(a.out, dpi=140)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
