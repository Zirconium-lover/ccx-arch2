#!/usr/bin/env bash
# Does CCX_DAMAGE_NLWIDTH reduce the ell dependence of the fracture energy?
#
#   CCX_EXE=/path/to/ccx_2.23_pardiso test/crackband/run_nlwidth_gate.sh <dir>
#
# WHY THIS COMPARES TWO RUNS AND NOT A NUMBER.  The quantity under test is
# the ratio of the largest to the smallest post-peak work across the
# nonlocal arms - 1.0 is a fracture energy independent of ell, which is what
# the model claims.  It is not there yet, so any fixed ceiling would have to
# be chosen to sit above whatever today's value happens to be, and a
# threshold chosen to pass is the exact mistake this directory already made
# once: a band width measured by a count over D>0.5 gave three verdicts for
# the three thresholds de1stats writes, and the published one was the middle.
#
# So the criterion is relative and has nothing to tune: the same sweep, the
# same binary, the same decks, run with the substitution off and on, and the
# dependence must come DOWN.  That is falsifiable - a substitution that does
# nothing, or makes it worse, fails - and it stays meaningful as the number
# approaches 1.
set -u
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
OUT=${1:-$(mktemp -d)}
: "${CCX_EXE:?set CCX_EXE to the ccx binary}"
mkdir -p "$OUT"
dep(){ # $1 = NLW value
  NLW=$1 bash "$ROOT/test/crackband/run_nlwidth_scaling.sh" "$OUT/nlw$1" \
      2>&1 | tee "$OUT/nlw$1.txt" | grep -o 'ELL_DEPENDENCE=[0-9.]*' \
      | cut -d= -f2
}
OFF=$(dep 0)
ON=$(dep 1)
echo
echo "  ell dependence of W_post across the nonlocal arms:"
echo "    CCX_DAMAGE_NLWIDTH=0   x$OFF"
echo "    CCX_DAMAGE_NLWIDTH=1   x$ON"
python3 - "$OFF" "$ON" <<'PY'
import sys
off,on=sys.argv[1],sys.argv[2]
if not off or not on:
    print("    could not measure both arms - no verdict"); sys.exit(2)
off,on=float(off),float(on)
print("    a fracture energy that does not depend on ell would be 1.0000")
if on < off:
    print("    PASS: the substitution reduced it, %.4f -> %.4f" % (off,on))
    sys.exit(0)
print("    FAIL: the substitution did not reduce it, %.4f -> %.4f" % (off,on))
sys.exit(1)
PY
