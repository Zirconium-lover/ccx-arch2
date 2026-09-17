#!/bin/bash
# EVOLUTION=ENERGY against EVOLUTION=DISPLACEMENT: do a G_f card and the
# equivalent u_f card describe the same material?
#
#   CCX_EXE=/path/to/ccx_2.23_pardiso test/crackband/run_energy_equiv.sh <dir>
#
# For the linear law, G_f = sigma_0*u_f/2 with sigma_0 the flow stress at
# initiation, so the deck writes the G_f equivalent to its own u_f and the two
# runs must agree.  The comparison is made at TWO strain levels, because the
# residual is expected to be a finite-strain effect: the conversion uses the
# stress in sti, and the work is done by the Cauchy stress.  Ten times smaller
# strain must shrink the disagreement, and that is the discriminating test -
# a fixed offset would not move.
#
# That criterion is CHECKED and sets the exit code.  It used to be a closing
# sentence telling the reader what a good result looks like, which left the
# script unable to fail however wrong the conversion became.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
OUT=${1:?usage: run_energy_equiv.sh <dir>}
EXE=${CCX_EXE:?set CCX_EXE to the ccx binary}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1} MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
mkdir -p "$OUT"

# label   E        eps0    uend
run_level(){
  local tag=$1 e=$2 eps0=$3 uend=$4
  for ev in DISPLACEMENT ENERGY; do
    local d="$OUT/${tag}_$ev"; rm -rf "$d"; mkdir -p "$d"
    python3 "$HERE/mkcross.py" --ncross 1 --hmod 10 --e "$e" --eps0 "$eps0" \
            --uf 0.02 --uend "$uend" --evolution "$ev" \
            -o "$d/t.inp" >/dev/null || return 2
    ( cd "$d" && "$EXE" t > run.log 2>&1 )
    echo "  $tag $ev exit=$?"
  done
  CB_HERE="$HERE" python3 - "$OUT" "$tag" "$uend" <<'PY'
import sys, os
sys.path.insert(0, os.environ["CB_HERE"])
from check_objectivity import curve, work
out, tag, uend = sys.argv[1], sys.argv[2], float(sys.argv[3])
a = curve(f"{out}/{tag}_DISPLACEMENT/t.dat", uend)
b = curve(f"{out}/{tag}_ENERGY/t.dat", uend)
wa, wb = work(a), work(b)
d = 100.0 * abs(wa - wb) / wa
print("  %-8s W(u_f)=%11.6f  W(G_f)=%11.6f  difference=%7.4f %%"
      % (tag, wa, wb, d))
open(f"{out}/{tag}.diff", "w").write("%.10f" % d)
PY
}

run_level strain1 200000  0.01  0.30
run_level strain01 2000000 0.001 0.10

# The criterion used to be a sentence telling the reader what to look for,
# which meant this script could not fail however wrong the conversion got.
# It is the RATIO that discriminates, not either value: a conversion wrong
# by a constant factor disagrees by about the same amount at both strain
# levels and so shows a ratio near 1, while an exact conversion leaves only
# the finite-strain residual, which shrinks with the strain.  Pinning the
# absolute difference instead would cement that residual, which is a
# documented limitation and not something to freeze.
echo
python3 - "$OUT" <<'PY'
import sys
out = sys.argv[1]
d1 = float(open(f"{out}/strain1.diff").read())
d01 = float(open(f"{out}/strain01.diff").read())
print("  difference at normal strain   : %7.4f %%" % d1)
print("  difference at one tenth of it : %7.4f %%" % d01)
FLOOR, NEED = 0.05, 5.0
if d01 <= FLOOR:
    print("  VERDICT ok   : %.4f %% is at the noise floor (<= %.2f %%), so"
          % (d01, FLOOR))
    print("                 the conversion is exact as far as this deck can")
    print("                 tell, and no ratio is needed to say so")
    sys.exit(0)
ratio = d1 / d01 if d01 > 0 else float("inf")
print("  ratio                         : %7.2f  (need >= %.1f)"
      % (ratio, NEED))
if ratio >= NEED:
    print("  VERDICT ok   : the disagreement shrinks with the strain, so it")
    print("                 is a finite-strain effect, not the conversion")
    sys.exit(0)
print("  VERDICT FAIL : the disagreement barely moved when the strain fell")
print("                 by ten.  That is the signature of a WRONG")
print("                 CONVERSION - a constant offset ignores the strain.")
sys.exit(1)
PY
