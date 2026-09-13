#!/usr/bin/env bash
set -u

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$SCRIPT_DIR/../.." && pwd)
# S3RAD_DECK lets a pilot run use a DERIVED deck (see mkpilotdeck.py).
# The committed deck's hash is then not the right thing to check, so the
# check is replaced by recording the derived deck's own hash in
# provenance.txt together with the fact that it is not the committed one.
DECK=${S3RAD_DECK:-"$SCRIPT_DIR/m12_s3rad_gc24_w.inp"}
EXE=${CCX_EXE:-"$ROOT/src/ccx_2.23"}
RUN_DIR=${1:-"$SCRIPT_DIR/_runs/$(date +%Y%m%d-%H%M%S)"}
if [ "$#" -gt 0 ]; then shift; fi

case "$EXE" in
    /*) ;;
    *) EXE=$(CDPATH= cd -- "$(dirname -- "$EXE")" && pwd)/$(basename -- "$EXE") ;;
esac
case "$RUN_DIR" in
    /*) ;;
    *) RUN_DIR="$PWD/$RUN_DIR" ;;
esac

EXPECTED_DECK_SHA=2fb0cf4e3554282e1f85cf641c788939bc7a2fa5918dd842f54d38844dfa5391

if [ ! -f "$DECK" ]; then
    echo "s3rad deck not found: $DECK" >&2
    exit 2
fi
if [ ! -x "$EXE" ]; then
    echo "CalculiX executable is not executable: $EXE" >&2
    echo "Set CCX_EXE to a PARDISO-enabled ccx_2.23 binary." >&2
    exit 2
fi

ACTUAL_DECK_SHA=$(sha256sum "$DECK" | awk '{print $1}')
DECK_IS_COMMITTED=yes
if [ "$ACTUAL_DECK_SHA" != "$EXPECTED_DECK_SHA" ]; then
    if [ -n "${S3RAD_DECK:-}" ]; then
        DECK_IS_COMMITTED=no
        echo "NOTE: running a DERIVED deck, not the committed one." >&2
        echo "      $DECK" >&2
        echo "      sha256 $ACTUAL_DECK_SHA" >&2
    else
        echo "deck hash mismatch" >&2
        echo "expected: $EXPECTED_DECK_SHA" >&2
        echo "actual:   $ACTUAL_DECK_SHA" >&2
        exit 2
    fi
fi

export OMP_NUM_THREADS=${OMP_NUM_THREADS:-6}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-6}
export MKL_CBWR=${MKL_CBWR:-COMPATIBLE}
# The deck's configuration.  ${NAME:-default} rather than a bare
# assignment, because a bare one SILENTLY OVERWRITES a value the caller
# set in the environment, and that has now cost two measurements: the
# CCX_DAMAGE_AUTOSPC case that made a gate case run with the mask on for
# its whole life, and a CCX_DAMAGE_DEADALL arm on 2026-09-12 that spent
# nineteen minutes reproducing the baseline.  Positional NAME=VALUE
# overrides below still win over both, and the provenance line says which
# values came from where, so nothing is decided silently.
export CCX_DAMAGE_AUTOSPC=${CCX_DAMAGE_AUTOSPC:-1.e-3}
export CCX_DAMAGE_DEADALL=${CCX_DAMAGE_DEADALL:-1.e-2}
export CCX_DAMAGE_DELETE_MAT=${CCX_DAMAGE_DELETE_MAT:-ALL}
export CCX_DAMAGE_LINESEARCH=${CCX_DAMAGE_LINESEARCH:-ADAPTIVE}
export CCX_DAMAGE_REEQ_RESCUE2=${CCX_DAMAGE_REEQ_RESCUE2:-1}
export CCX_DAMAGE_REEQ_SCALE=${CCX_DAMAGE_REEQ_SCALE:-PHYSICAL}
export CCX_DAMAGE_TANGENT=${CCX_DAMAGE_TANGENT:-UNSYM}
export CCX_DAMAGE_TOPOLOGY=${CCX_DAMAGE_TOPOLOGY:-DEFERRED}
export CCX_DAMAGE_TR_DOGLEG=${CCX_DAMAGE_TR_DOGLEG:-1}
export CCX_DAMAGE_VISCOSITY=${CCX_DAMAGE_VISCOSITY:-1.e-4}
export CCX_FRACTURE_TERMINATION=${CCX_FRACTURE_TERMINATION:-FACE_X0_NSET:FACE_XL_NSET}
# Stop when the load path between the termination sets has narrowed to this
# fraction of its width at the first committed deletion batch.  The five
# topological rules all say CONNECTED at the end of a run whose grips are
# joined by one triangular face (research/09-SEVERANCE.md), so without a
# WIDTH the deck has no ending of its own and reports rc=201 - the solver
# giving up - as its result.
#
# 1e-4 is measured, not chosen.  At 1e-2 the criterion fires at theta
# 0.2427, BEFORE the wall, so it would end the run early and hide
# everything past it.  At the wall itself (theta 0.2550) the cut is 0.0031
# against a reference of 4.1547, i.e. ratio 7.5e-4, so any threshold above
# that also fires too early.  1e-4 is an order of magnitude tighter than
# "one face at the residual-stiffness floor" and fires at ratio 4.7e-05,
# theta 0.2588, with the grip reaction at 1.46 percent of peak:
# [FRACTURE COMPLETE] instead of rc=201 (research/14, 15, 17).
#
# CONFIRMED on this recipe with no overrides at all: the run ends
# [FRACTURE COMPLETE] inc=667 step_time=2.587500e-01, cut=1.962064e-04
# ratio=4.722478e-05, and its five output files are byte for byte the
# measured arm's (wall clock excepted).  The default IS the measurement.
export CCX_FRACTURE_CUT=${CCX_FRACTURE_CUT:-1e-4}
export CCX_PARDISO_REUSE_SYMBOLIC=${CCX_PARDISO_REUSE_SYMBOLIC:-1}



unset CCX_DISSIPATION_CONTROL
unset CCX_DISSIPATION_TARGET
unset CCX_FRACTURE_LINK
unset CCX_FRACTURE_DEADFACET

for kv in "$@"; do
    name=${kv%%=*}
    value=${kv#*=}
    case "$name" in
        ''|*[!A-Za-z0-9_]*)
            echo "invalid environment override: $kv" >&2
            exit 2
            ;;
    esac
    if [ "$value" = "$kv" ]; then
        echo "override must be NAME=VALUE: $kv" >&2
        exit 2
    elif [ -z "$value" ]; then
        unset "$name"
    else
        export "$kv"
    fi
done

mkdir -p "$RUN_DIR"
cp "$DECK" "$RUN_DIR/m.inp"

EXE_SHA=$(sha256sum "$EXE" | awk '{print $1}')
{
    echo "deck=$DECK"
    echo "deck_sha256=$ACTUAL_DECK_SHA"
    echo "deck_is_committed=$DECK_IS_COMMITTED"
    echo "executable=$EXE"
    echo "executable_sha256=$EXE_SHA"
    echo "started=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "overrides=$*"
    echo "--- environment ---"
    env | grep -E '^(CCX_|OMP_|MKL_)' | sort
} > "$RUN_DIR/provenance.txt"

echo "Run directory: $RUN_DIR"
echo "Deck SHA-256: $ACTUAL_DECK_SHA"
echo "Executable SHA-256: $EXE_SHA"
echo "The committed deck requests PARDISO."

start=$(date +%s)
(
    cd "$RUN_DIR" || exit 2
    "$EXE" -i m > run.log 2>&1
)
rc=$?
finish=$(date +%s)

{
    echo "return_code=$rc"
    echo "wall_seconds=$((finish-start))"
} >> "$RUN_DIR/provenance.txt"

echo "Return code: $rc"
echo "Wall seconds: $((finish-start))"
tail -n 5 "$RUN_DIR/m.sta" 2>/dev/null || true
exit "$rc"
