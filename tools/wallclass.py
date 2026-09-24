#!/usr/bin/env python3
"""Where did a run's wall come from - the step, or the state before it?

    tools/wallclass.py <run-dir> [<run-dir> ...]
    tools/wallclass.py --selftest

A run "hits a wall" when an increment fails on every attempt and the solver
gives up after cutting the step down.  Cutting the step is the remedy for a
residual the LOAD INCREMENT creates: halve the increment, and that residual
shrinks with it.  It cannot touch a residual that was already there when the
increment began - an out-of-balance force the committed state carried in, for
example the force an element held when it was deleted, or a damage increment
applied at the start of the step whatever the step's size.

The two are told apart by one number, read off the .cvg the solver already
writes: the first-iteration residual R1 of each attempt, against that
attempt's increment size.

    ratio = R1(smallest attempt) / R1(largest attempt)

  STEP-DRIVEN   ratio < 0.3: the residual fell by more than three when the
                step was cut, so the step made it and cutting was the right
                remedy, it simply did not reach far enough.
  CARRIED-IN    ratio >= 0.3: the residual stayed put while the step fell by
                orders of magnitude, so the increment started out of balance
                and no cutback can fix it.
  UNDETERMINED  the attempts did not span a factor of four in step size -
                typically a rescue retrying at the minimum step - so the
                ratio cannot say how the residual depends on the step.

The 0.3 is not tuned.  Every cutback here divides the step by at least four,
and a residual proportional to the step divides with it; one that did not
fall by three across all attempts is not proportional to the step.

The report also names the commit events of the increment before the wall,
because that is where a carried-in residual would come from - the gate's own
fast-wrapped wall follows an element deletion at the increment before it.

A run with no fatal increment is reported as such and is not an error.
"""
import os, re, sys, glob, tempfile, shutil

THRESH = 0.3
STA_ROW = re.compile(r'^\s+(\d+)\s+(\d+)\s+(\d+)(U?)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)')


def stem_of(d):
    for f in sorted(os.listdir(d)):
        if f.endswith('.sta'):
            return f[:-4]
    return None


def attempts(sta):
    """(inc, att, unconverged, inc_time) for every attempt line of the .sta"""
    out = []
    for ln in open(sta, errors='replace'):
        m = STA_ROW.match(ln)
        if m:
            out.append((int(m.group(2)), int(m.group(3)), m.group(4) == 'U',
                        float(m.group(8))))
    return out


def first_residuals(cvg, inc):
    """{attempt: R1} - the residual of iteration 1 of each attempt of inc"""
    r = {}
    for ln in open(cvg, errors='replace'):
        t = ln.split()
        if len(t) < 6 or t[1] != str(inc) or t[3] != '1':
            continue
        try:
            r.setdefault(int(t[2]), float(t[5]))
        except ValueError:
            pass
    return r


def commits_before(log, inc):
    """[DAMAGE ... COMMIT] lines of the increment just before the wall"""
    if not log or not os.path.exists(log):
        return []
    pat = re.compile(r'\[DAMAGE[^\]]*COMMIT\][^\n]*\binc=%d\b' % (inc - 1))
    return [ln.strip() for ln in open(log, errors='replace') if pat.search(ln)]


def classify(d):
    """dict describing the run's wall, or None when there is no wall"""
    stem = stem_of(d)
    if not stem:
        return {'dir': d, 'error': 'no .sta'}
    sta = os.path.join(d, stem + '.sta')
    cvg = os.path.join(d, stem + '.cvg')
    rows = attempts(sta)
    if not rows or not rows[-1][2]:
        return None
    inc = rows[-1][0]
    att = [r for r in rows if r[0] == inc]
    if not os.path.exists(cvg):
        return {'dir': d, 'inc': inc, 'error': 'no .cvg'}
    r1 = first_residuals(cvg, inc)
    pts = [(a[3], r1[a[1]]) for a in att if a[1] in r1]
    if len(pts) < 2:
        return {'dir': d, 'inc': inc, 'error': 'fewer than two attempts'}
    big = max(pts, key=lambda p: p[0])
    small = min(pts, key=lambda p: p[0])
    ratio = small[1] / big[1] if big[1] > 0 else float('inf')
    log = os.path.join(d, 'run.log')
    # The ratio only means something if the step actually changed.  A rescue
    # that retries at the minimum step gives several attempts at ONE step
    # size, and a constant residual across them says nothing about how it
    # depends on the step.  The first version of this tool counted those as
    # CARRIED-IN, and so did the count its docstring quotes.
    if big[0] < 4.0 * small[0]:
        kind = 'UNDETERMINED'
    else:
        kind = 'CARRIED-IN' if ratio >= THRESH else 'STEP-DRIVEN'
    return {'dir': d, 'inc': inc, 'n': len(pts),
            'dt': (big[0], small[0]), 'r1': (big[1], small[1]),
            'ratio': ratio, 'kind': kind,
            'commits': commits_before(log, inc)}


def report(c):
    if c is None:
        return "no wall: the last increment converged"
    if 'error' in c:
        return "cannot judge: %s" % c['error']
    s = ("%s at inc %d: R1 %.4g at dt %.3g -> %.4g at dt %.3g over %d "
         "attempts, ratio %.2f"
         % (c['kind'], c['inc'], c['r1'][0], c['dt'][0], c['r1'][1],
            c['dt'][1], c['n'], c['ratio']))
    for ln in c['commits']:
        s += "\n    before the wall: " + ln[:150]
    return s


def selftest():
    """two synthetic runs with known answers, one of each kind"""
    tmp = tempfile.mkdtemp(prefix='wallclass_')
    bad = 0
    try:
        def make(name, r1s, commit):
            d = os.path.join(tmp, name)
            os.makedirs(d)
            with open(os.path.join(d, 'm.sta'), 'w') as f:
                f.write("     1         9     1     3  0.1E+00  0.1E+00  0.1E-01\n")
                dt = 1e-2
                for k, _ in enumerate(r1s, 1):
                    f.write("     1        10     %dU    6  0.1E+00  0.1E+00  %.6E\n"
                            % (k, dt))
                    dt /= 4
            with open(os.path.join(d, 'm.cvg'), 'w') as f:
                for k, r in enumerate(r1s, 1):
                    f.write("     1    10     %d     1        0  %.4E  0.1E+03\n"
                            % (k, r))
                    f.write("     1    10     %d     2        0  %.4E  0.1E+03\n"
                            % (k, r * 3))
            with open(os.path.join(d, 'run.log'), 'w') as f:
                if commit:
                    f.write("[DAMAGE COMMIT] batch=1 inc=9 time=0.1 deleted=1\n")
            return d
        carried = make('carried', [294.0, 294.0, 294.0, 294.0], True)
        driven = make('driven', [800.0, 200.0, 50.0, 12.0], False)
        c = classify(carried)
        if not c or c.get('kind') != 'CARRIED-IN' or len(c['commits']) != 1:
            print("FAILED: constant residual not classified CARRIED-IN "
                  "with its commit: %r" % c)
            bad += 1
        c = classify(make('onestep', [294.0, 294.0, 294.0], False))
        # same residual AND same step on every attempt is not evidence
        # either way; the helper halves dt, so force equal steps by hand
        d1 = os.path.join(tmp, 'onestep')
        with open(os.path.join(d1, 'm.sta'), 'w') as f:
            f.write("     1         9     1     3  0.1E+00  0.1E+00  0.1E-01\n")
            for k in (1, 2, 3):
                f.write("     1        10     %dU    6  0.1E+00  0.1E+00  3.9E-06\n" % k)
        c = classify(d1)
        if not c or c.get('kind') != 'UNDETERMINED':
            print("FAILED: attempts at one step size were given a verdict: %r" % c)
            bad += 1
        c = classify(driven)
        if not c or c.get('kind') != 'STEP-DRIVEN' or c['commits']:
            print("FAILED: shrinking residual not classified STEP-DRIVEN: %r" % c)
            bad += 1
        # a converged run has no wall and must not be reported as one
        d = os.path.join(tmp, 'clean')
        os.makedirs(d)
        with open(os.path.join(d, 'm.sta'), 'w') as f:
            f.write("     1        10     1     3  0.1E+01  0.1E+01  0.1E-01\n")
        if classify(d) is not None:
            print("FAILED: a converged run was reported as a wall")
            bad += 1
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print("[WALLCLASS SELFTEST] %s" % ("PASSED" if bad == 0 else "FAILED"))
    return 1 if bad else 0


def main(argv):
    if '--selftest' in argv:
        return selftest()
    if not argv:
        print(__doc__.split('\n\n')[0])
        return 2
    rc = 0
    for d in argv:
        if not os.path.isdir(d):
            print("%s: not a directory" % d)
            rc = 2
            continue
        print("%s: %s" % (d, report(classify(d))))
    return rc


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
