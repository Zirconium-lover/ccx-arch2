#!/usr/bin/env bash
# Build and run the damnlelmk self test.
#
#   test/nonlocal/run_elmk_test.sh [out-dir]
#
# Checks the element mass and conductivity of the implicit-gradient
# equation against closed forms and against three family-independent
# identities.  Links against the solver archive, so src/ must have been
# built.  Returns non-zero when any case disagrees; that path has been
# exercised by changing an expected value, which turns that case red and
# the exit code to 1.
set -u
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
OUT=${1:-$(mktemp -d)}
mkdir -p "$OUT"
A="$ROOT/src/ccx_2.23.a"
[ -f "$A" ] || { echo "build src/ first: $A is missing" >&2; exit 2; }
MKL="-Wl,--start-group /usr/lib/x86_64-linux-gnu/libmkl_gf_lp64.a
     /usr/lib/x86_64-linux-gnu/libmkl_gnu_thread.a
     /usr/lib/x86_64-linux-gnu/libmkl_core.a -Wl,--end-group"
# shellcheck disable=SC2086
gfortran -O0 -o "$OUT/elmktest" "$ROOT/test/nonlocal/elmk_test.f" \
    "$A" -lspooles -larpack $MKL -fopenmp -lpthread -lm -ldl || exit 2
"$OUT/elmktest"
