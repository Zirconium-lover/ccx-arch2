#!/bin/bash
# EVOLUTION=ENERGY against EVOLUTION=DISPLACEMENT: do a G_f card and the
# equivalent u_f card describe the same material?
#
#   CCX_EXE=/path/to/ccx_2.23_pardiso test/crackband/run_energy_equiv.sh <dir>
#
# For the linear law, G_f = sigma_0*u_f/2 with sigma_0 the flow stress at
# initiation, so the deck writes the G_f equivalent to its own u_f and the two
# runs must agree.  The comparison is made at TWO strain levels a factor of
# ten apart, and BOTH are bounded; the criterion sets the exit code rather
# than being a closing sentence for a reader to interpret.
#
# RETRACTED, and left here because the retraction is the useful part: this
# header used to say the remaining disagreement "is expected to be a
# finite-strain effect - the conversion uses the stress in sti, and the work
# is done by the Cauchy stress", and that ten times smaller strain shrinking
# it was "the discriminating test".  It was not discriminating.  The residual
# was sigma_0 sampled after the point had begun to unload, which scales with
# strain in exactly the same way, so the test could not tell the two apart -
# and the push-forward to Cauchy that the old explanation implies was built
# and measured and makes the answer WORSE.  See src/calcdamage.f, damcbufset.
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
#
# RETRACTED, third instance of the same thing, found by sweeping the CODE
# for a withdrawn claim rather than only the discussion of it.  This block
# used to argue that the ratio discriminates and that "pinning the absolute
# difference would cement that residual, which is a documented limitation
# and not something to freeze".  Both halves were wrong, and both followed
# from the same retracted diagnosis: the residual was not a limitation, it
# was a defect worth 36x, and pinning the absolute difference is exactly
# what catches it.  An argument for NOT checking something is the most
# expensive kind to leave standing after its premise goes.
echo
python3 - "$OUT" <<'PY'
import sys
out = sys.argv[1]
d1 = float(open(f"{out}/strain1.diff").read())
d01 = float(open(f"{out}/strain01.diff").read())
print("  difference at normal strain   : %7.4f %%" % d1)
print("  difference at one tenth of it : %7.4f %%" % d01)
# WHAT DECIDES, AND WHY IT IS NOT THE RATIO ANY MORE.
#
# This script used to judge by the RATIO between the two levels: an exact
# conversion was said to leave only a residual that shrinks with strain,
# so a large ratio meant "fine".  That reasoning named a cause - a
# finite-strain stress-measure effect - which has since been retracted.
# The residual was not that at all; it was sigma_0 sampled after the point
# had begun to unload.
#
# Agent 1's review then found what the retraction had left behind: the
# ratio branch still EXITED GREEN, and did so while asserting the
# retracted cause.  It was the only branch reachable once the disagreement
# rose above the floor, so precisely when something had gone wrong the
# script would have explained it away.
#
# So the ratio no longer decides.  It could not anyway: both failures this
# deck can show blow a ceiling long before the ratio notices, measured -
# a conversion wrong by a constant factor gives 12.60 % at working strain,
# and the freezing defect gave 2.98 %.  The ratio is still PRINTED,
# because how the disagreement scales is genuinely informative about what
# kind of residual is left, but information is not a verdict.
#
# Both levels are now bounded directly, each from measurement on 37d0d53
# with headroom: 0.0822 % against 0.5, and 0.0082 % against 0.05.
CEIL1, CEIL01 = 0.5, 0.05
bad = []
if d1 > CEIL1:
    bad.append("  %.4f %% at working strain exceeds the %.2f %% ceiling"
               % (d1, CEIL1))
if d01 > CEIL01:
    bad.append("  %.4f %% at one tenth strain exceeds the %.2f %% ceiling"
               % (d01, CEIL01))
if d01 > 0:
    print("  ratio                         : %7.2f  (reported, does not"
          " decide)" % (d1 / d01))
if bad:
    print("  VERDICT FAIL :")
    for line in bad:
        print("  " + line)
    print("                 The two cards do not describe the same material.")
    print("                 No cause is named here on purpose: above these")
    print("                 bounds we do not have one, and the last time a")
    print("                 cause was named from the way the residual")
    print("                 scaled, it was the wrong cause.")
    sys.exit(1)
print("  VERDICT ok   : both levels are within bounds set by measurement")
sys.exit(0)
ratio = d1 / d01 if d01 > 0 else float("inf")
print("  ratio                         : %7.2f  (need >= %.1f)"
      % (ratio, NEED))
if ratio < NEED:
    print("  VERDICT FAIL : the disagreement barely moved when the strain")
    print("                 fell by ten.  That is the signature of a WRONG")
    print("                 CONVERSION - a constant offset ignores strain.")
    bad = True
elif not bad:
    print("  VERDICT ok   : the disagreement shrinks with the strain and")
    print("                 stays under the ceiling")
sys.exit(1 if bad else 0)
PY
