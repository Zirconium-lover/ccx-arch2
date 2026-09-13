#!/usr/bin/env bash
# Generate and run the fast fracture regression specimen.
#
#   CCX_EXE=/path/to/ccx_2.23_pardiso test/fast/run_fast.sh [run-dir] [NAME=VALUE ...]
#
# The deck is generated, not committed, so it can never drift from the
# generator; the generated deck's sha256 is recorded in provenance.txt by
# run_s3rad.sh, which this script reuses so that the solver environment is
# character for character the one the target deck runs under.
#
# FAST_VARIANT selects the specimen:
#   wrapped (default) - the eroding phase wrapped in cohesive facets over half
#                       the width.  Reaches a node with NO bulk support left,
#                       held by facets alone, at the standard AUTOSPC=1.e-3.
#   plain             - the original partial-cohesive-plane bar.  Completes at
#                       theta=1 and never masks a node.
set -u
SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$SCRIPT_DIR/../.." && pwd)
RUN_DIR=${1:-"$SCRIPT_DIR/_runs/$(date +%Y%m%d-%H%M%S)"}
if [ "$#" -gt 0 ]; then shift; fi
case "$RUN_DIR" in /*) ;; *) RUN_DIR="$PWD/$RUN_DIR" ;; esac

VARIANT=${FAST_VARIANT:-wrapped}
case "$VARIANT" in
    wrapped) GEN_ARGS="--incl 6,3,0,2,3,3" ;;
    plain)   GEN_ARGS="" ;;
    *) echo "unknown FAST_VARIANT: $VARIANT" >&2; exit 2 ;;
esac

mkdir -p "$RUN_DIR"
DECK="$RUN_DIR/fast.inp"
python3 "$SCRIPT_DIR/mkfast.py" -o "$DECK" --nx 10 --ny 6 --nz 6 \
        --from-deck "$ROOT/test/s3rad/m12_s3rad_gc24_w.inp" $GEN_ARGS || exit 2

S3RAD_DECK="$DECK" "$ROOT/test/s3rad/run_s3rad.sh" "$RUN_DIR" "$@"
rc=$?

L="$RUN_DIR/run.log"
if [ -f "$L" ]; then
    echo "--- load-path census (variant=$VARIANT) ---"
    for t in 1e-1 1e-2 1e-3; do
        printf "  nodes below %s of their own intact diagonal, max over the run: %s\n" \
            "$t" "$(grep -o "below_$t=[0-9]*" "$L" | sed 's/.*=//' | sort -n | tail -1)"
    done
    printf "  worst ratio reached: %s\n" \
        "$(grep -o 'worst=node_[0-9]*_at_[0-9.e+-]*' "$L" | sed 's/worst=node_//;s/_at_/ at /' | sort -k3 -g | head -1)"
    printf "  elements terminally deleted: %s\n" \
        "$(grep -vc '^#' "$RUN_DIR/m.damage" 2>/dev/null)"
fi
exit "$rc"
