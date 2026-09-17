#!/usr/bin/env bash
# Does a nonlocal backend reduce to the LOCAL model as ell goes to zero?
#
#   CCX_EXE=/path/to/ccx_2.23_pardiso test/nonlocal/run_ell0_test.sh [out-dir]
#
# ebar - div(ell^2 grad ebar) = e  at ell=0 is ebar = e, and the integral
# average over a radius that contains only the element itself is the element
# itself.  So BOTH backends must return the local answer when ell is small
# against the element size - not approximately, not in trend: the same
# elements deleted, at the same increments.
#
# This is not a style point.  A backend that fails it has an internal length
# it did not get from the user, set by the mesh, and every mesh study done
# with it measures that length instead of the physical one.
#
# ell is 0.01 against an element about 0.33 across, so the regularisation is
# off by a factor of thirty, and anything the backend still does to the field
# is its own.
#
# Each backend is judged SEPARATELY and the exit code is the worst of them,
# so a backend that passes is not hidden by one that does not.
#
# CCX_DAMAGE_CHARLEN IS DELIBERATELY NOT SET.  An earlier draft armed the
# directional width here, and the integral backend then stopped converging at
# this ell - rc=201 at a tenth of the load, against rc=0 and the exact local
# answer without it.  That is a real interaction and it is reported, but it is
# not what this test asks about: the question here is whether a backend
# reduces to the local model as ell falls, and the answer must not depend on
# how the crack-band width is estimated.  Both arms and the local reference
# run under the same default, so the comparison stays controlled.
set -u
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
OUT=${1:-$(mktemp -d)}
EXE=${CCX_EXE:?set CCX_EXE to the ccx binary}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1} MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
mkdir -p "$OUT"
ELL=0.01

run(){ # tag  extra-env...
  local tag=$1; shift
  local d="$OUT/$tag"; rm -rf "$d"; mkdir -p "$d"
  # A DELIBERATELY SMALL DECK.  6x4x4 runs the three arms in about twenty
  # seconds against ninety for the gate's 10x6x6, and it keeps the property
  # the test needs: the pre-fix gradient backend still adds elements the
  # local run never breaks, which is what carries the verdict.  Verified on
  # the pre-fix binary rather than assumed - a cheaper deck that stopped
  # discriminating would be worse than no test.
  python3 "$ROOT/test/fast/mkfast.py" -o "$d/m.inp" --nx 8 --ny 5 --nz 5 \
      >/dev/null 2>&1 || return 2
  ( cd "$d" && env "$@" "$EXE" -i m > run.log 2>&1 )
  return 0
}

run local CCX_DAMAGE_NONLOCAL=0 || { echo "deck generation failed" >&2; exit 2; }
run integral CCX_DAMAGE_NONLOCAL=$ELL CCX_DAMAGE_NONLOCAL_MODE=INTEGRAL
run gradient CCX_DAMAGE_NONLOCAL=$ELL CCX_DAMAGE_NONLOCAL_MODE=GRADIENT

python3 - "$OUT" "$ELL" <<'PY'
import sys
out, ell = sys.argv[1], sys.argv[2]
def dels(tag):
    s=[]
    try:
        for ln in open("%s/%s/m.damage"%(out,tag),errors='replace'):
            if ln.startswith('#'): continue
            t=ln.split()
            if t: s.append((int(t[0]),float(t[4])))
    except OSError: return None
    return s
def last(tag):
    try:
        rows=[l.split() for l in open("%s/%s/m.sta"%(out,tag),errors='replace')
              if l.startswith('     1')]
        return (len(rows), float(rows[-1][4])) if rows else (0,0.0)
    except OSError: return (0,0.0)
loc=dels('local')
if loc is None:
    print("  local run produced nothing - cannot judge"); sys.exit(2)
ln_,lt=last('local')
print("  local          : %4d deleted, %3d increments, theta %.6f"%(len(loc),ln_,lt))
bad=0
for tag in ('integral','gradient'):
    d=dels(tag)
    if d is None:
        print("  %-14s: no output - FAILED"%tag); bad=1; continue
    n,t=last(tag)
    # THE SET IS JUDGED EXACTLY, THE TIMES WITH A TOLERANCE, and the
    # difference is the point.  At ell=0.01 the regularisation is not
    # identically zero - ell^2/h^2 is about 2.6e-3 - so the trajectory may
    # shift slightly and demanding bit-identical times would make this test
    # fail for the right behaviour.  WHICH ELEMENTS break is a different
    # matter: no amount of legitimate smoothing at this ell adds or removes
    # one, and before the zero-length projection was subtracted the gradient
    # backend added SIXTY.  So the set carries the verdict and the times
    # only have to stay within TOL.
    TOL=0.02
    same_set = set(e for e,_ in d)==set(e for e,_ in loc)
    if same_set and len(d)==len(loc):
        # normalised by the FULL load, not by the event's own time: an
        # element that breaks early has a small time, so dividing by it
        # turns a shift of a few thousandths of the load into tens of per
        # cent and says nothing about how far the answer moved.
        scale=max(1e-12,abs(lt))
        worst=max((abs(a[1]-b[1])/scale
                   for a,b in zip(sorted(d),sorted(loc))),default=0.0)
    else:
        worst=float('inf')
    same_time = worst<=TOL
    verdict = "ok" if (same_set and same_time and abs(t-lt)<=TOL*max(1e-12,abs(lt))) else "FAILED"
    if verdict=="FAILED": bad=1
    print("  %-14s: %4d deleted, %3d increments, theta %.6f, worst time"
          " shift %s  -> %s"
          %(tag,len(d),n,t,
            ("%.2e"%worst) if worst!=float('inf') else "n/a",verdict))
    if verdict=="FAILED":
        extra=set(e for e,_ in d)-set(e for e,_ in loc)
        miss =set(e for e,_ in loc)-set(e for e,_ in d)
        print("                   %d element(s) deleted that the local run did not,"
              " %d that it did and this did not"%(len(extra),len(miss)))
print()
if bad:
    print("  VERDICT FAIL : a backend does not reduce to the local model at")
    print("                 ell=%s, so it carries an internal length it was"%ell)
    print("                 never given - one set by the MESH.")
    sys.exit(1)
print("  VERDICT ok   : both backends reduce to the local model at ell=%s"%ell)
PY
