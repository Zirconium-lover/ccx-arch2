#!/bin/bash
# Crack-band objectivity sweep: refine the cross-section only, on a fixed
# specimen, and see whether the dissipated energy moves.
#
#   CCX_EXE=/path/to/ccx_2.23_pardiso test/crackband/run_objectivity.sh <dir>
#
# The band is one slice thick in every run, so its width does not change.
# Only the legacy (6V)^(1/3) length changes, as ncross^(-2/3).
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
OUT=${1:?usage: run_objectivity.sh <dir>}
EXE=${CCX_EXE:?set CCX_EXE to the ccx binary}
UEND=0.30
NCROSS="1 2 3"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1} MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
mkdir -p "$OUT"

for nc in $NCROSS; do
  python3 "$HERE/mkcross.py" --ncross "$nc" --uend "$UEND" \
          -o "$OUT/cx_$nc.inp" || exit 2
  python3 "$HERE/mkcross.py" --ncross "$nc" --uend "$UEND" --eltype C3D8 \
          -o "$OUT/h8_$nc.inp" || exit 2
done

for nc in $NCROSS; do
  for arm in legacy proj; do
    d="$OUT/${nc}_$arm"; rm -rf "$d"; mkdir -p "$d"
    cp "$OUT/cx_$nc.inp" "$d/t.inp"
    if [ "$arm" = proj ]; then ex="CCX_DAMAGE_CHARLEN=1"; else ex=""; fi
    ( cd "$d" && env $ex "$EXE" t > run.log 2>&1 )
    echo "ncross=$nc C3D4 $arm exit=$?"
  done
  # C3D8 has ONE arm by construction: the legacy (6V)^(1/3) is defined for
  # tetrahedra only and refuses this family, so the projection is not merely
  # more objective here - it is the only way to run the model at all.
  d="$OUT/${nc}_hex"; rm -rf "$d"; mkdir -p "$d"
  cp "$OUT/h8_$nc.inp" "$d/t.inp"
  ( cd "$d" && env CCX_DAMAGE_CHARLEN=1 "$EXE" t > run.log 2>&1 )
  echo "ncross=$nc C3D8 proj exit=$?  (nonzero = branch not completed;"\
       "the window below is what all three meshes reached)"
done

# The legacy length per mesh, from the deck headers the generator wrote.
lof(){ grep "legacy    L=" "$OUT/cx_$1.inp" | sed 's/.*= *//'; }

for arm in legacy proj hex; do
  echo
  if [ "$arm" = hex ]; then echo "=== C3D8, projection (no legacy arm exists) ==="
  else echo "=== C3D4, arm=$arm ==="; fi
  args=""
  for nc in $NCROSS; do
    if [ "$arm" = legacy ]; then L=$(lof "$nc"); else L=1.0; fi
    args="$args ncross=$nc:$L:$OUT/${nc}_$arm/t.dat"
  done
  # shellcheck disable=SC2086
  python3 "$HERE/check_objectivity.py" "$UEND" $args
done
