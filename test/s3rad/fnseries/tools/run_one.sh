#!/usr/bin/env bash
# run_one.sh <deck.inp> <run_dir> [NAME=VALUE ...]
# One Fn-series run with the flag set of the disc series that reached
# FRACTURE COMPLETE (fractrure-tests/runs/run.sh, solver c75ad9b).
# THREADS (default 1) and TLIMIT (seconds, default 6 h) from the environment.
deck=$(readlink -f "$1"); D=$(readlink -m "$2"); shift 2
ROOT=/home/user/ccx-fn
mkdir -p "$(dirname "$D")"
S3RAD_DECK=$deck OMP_NUM_THREADS=${THREADS:-1} MKL_NUM_THREADS=${THREADS:-1} \
CCX_EXE=$ROOT/src/ccx_2.23_pardiso timeout ${TLIMIT:-21600} bash $ROOT/test/s3rad/run_s3rad.sh "$D" \
  CCX_DAMAGE_VISCOUS_DAMPING=2e-4 CCX_DAMAGE_FACET_DEBRIS=1 CCX_DAMAGE_GRADUAL_DELETE=8 \
  CCX_DAMAGE_DEADSOLE=1e-3 CCX_DAMAGE_DEADSOLE_LAW=1 \
  CCX_DAMAGE_BARE_MASK=0.95 CCX_DAMAGE_BARE_MASK_LAW=1 CCX_DAMAGE_HINGE=0.95 \
  CCX_DAMAGE_HINGE_PENDANT=2 CCX_DAMAGE_HINGE_CLUSTER=1 CCX_DAMAGE_HINGE_CLUSTER_BIG=1 \
  CCX_DAMAGE_TR_STICKY=1 CCX_DAMAGE_RESCUE_CUTBACKS=1 CCX_DAMAGE_TR_MAXEVAL=100000 \
  CCX_DAMAGE_TR_MAXFACT=100000 CCX_DAMAGE_TR_MAXARM=1000 "$@" > "$D.out" 2>&1
rc=$?
echo "$(basename "$D") rc=$rc $(tail -1 "$D/m.sta" 2>/dev/null | awk '{print "inc",$2,"theta",$6}') del=$(grep -vc '^#' "$D/m.damage" 2>/dev/null) $(grep -a -h 'FRACTURE COMPLETE\] inc' "$D/run.log" 2>/dev/null | head -1 | cut -c1-60)" >> "$(dirname "$D")/summary.txt"
exit $rc
