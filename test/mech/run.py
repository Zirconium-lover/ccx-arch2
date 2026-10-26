#!/usr/bin/env python3
"""Run the mechanics tests: small decks checked against a reference.

    CCX_EXE=$PWD/src/ccx_2.23_pardiso python3 test/mech/run.py [name ...] [--work DIR] [-v]

Each test/mech/<name>/test.py generates its own decks, runs them and returns
checks (measured value, reference, tolerance, why).  The exit code is the
number of tests with a failing check, so a CI step or the regression gate can
use it directly.  XFAIL is a known limitation that is still failing as
documented; XPASS means it stopped failing and the test must be revisited, so
it counts as a failure.
"""
import argparse
import importlib.util
import os
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mechlib  # noqa: E402


def tests():
    return sorted(d for d in os.listdir(HERE)
                  if os.path.isfile(os.path.join(HERE, d, "test.py")))


def load(name):
    spec = importlib.util.spec_from_file_location("mech_" + name,
                                                  os.path.join(HERE, name, "test.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("names", nargs="*")
    ap.add_argument("--work", help="keep run directories here (default: a temp dir)")
    ap.add_argument("-l", "--list", action="store_true")
    a = ap.parse_args()
    names = a.names or tests()
    if a.list:
        for n in tests():
            print("%-22s %s" % (n, load(n).WHAT))
        return 0
    exe = os.environ.get("CCX_EXE", os.path.join(mechlib.ROOT, "src", "ccx_2.23_pardiso"))
    if not os.access(exe, os.X_OK):
        print("no solver at %s: set CCX_EXE" % exe)
        return 2
    work = a.work or tempfile.mkdtemp(prefix="mech_")
    summary, nbad = [], 0
    for n in names:
        if n not in tests():
            print("unknown test %s (have: %s)" % (n, ", ".join(tests())))
            return 2
        m = load(n)
        print("== %s: %s" % (n, m.WHAT))
        print("   reference: %s" % m.SOURCE)
        t0 = time.time()
        try:
            checks = m.run(os.path.join(work, n), exe)
            err = None
        except Exception as e:  # a crashed test is a failed test, not a skipped one
            checks, err = [], "%s: %s" % (type(e).__name__, e)
        dt = time.time() - t0
        for c in checks:
            print(c.line())
        ok = err is None and checks and all(c.counts_as_pass for c in checks)
        if err:
            print("  ERROR %s" % err)
        if dt > 60:
            print("  SLOW  %.0f s exceeds the 60 s budget of a mechanics test" % dt)
            ok = False
        nbad += 0 if ok else 1
        cnt = {}
        for c in checks:
            cnt[c.status] = cnt.get(c.status, 0) + 1
        summary.append((n, "pass" if ok else "FAIL", dt,
                        " ".join("%s=%d" % kv for kv in sorted(cnt.items()))))
        print()
    print("%-22s %-5s %7s  %s" % ("test", "", "time", "checks"))
    for n, st, dt, cnt in summary:
        print("%-22s %-5s %6.1fs  %s" % (n, st, dt, cnt))
    print("work dir: %s" % work)
    return nbad


if __name__ == "__main__":
    sys.exit(main())
