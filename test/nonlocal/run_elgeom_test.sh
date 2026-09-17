#!/usr/bin/env bash
# Build and run the damnlelgeom self test.
#
#   test/nonlocal/run_elgeom_test.sh
#
# Links against the solver archive, so src/ must have been built. Returns
# non-zero when any case disagrees with its hand-computed answer; that
# path has been exercised by changing an expected volume, which turns the
# first case red and the exit code to 1.
set -u
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
OUT=${1:-$(mktemp -d)}
A="$ROOT/src/ccx_2.23.a"
[ -f "$A" ] || { echo "build src/ first: $A is missing" >&2; exit 2; }
MKL="-Wl,--start-group /usr/lib/x86_64-linux-gnu/libmkl_gf_lp64.a
     /usr/lib/x86_64-linux-gnu/libmkl_gnu_thread.a
     /usr/lib/x86_64-linux-gnu/libmkl_core.a -Wl,--end-group"
# shellcheck disable=SC2086
gfortran -O0 -o "$OUT/elgeomtest" "$ROOT/test/nonlocal/elgeom_test.f" \
    "$A" -lspooles -larpack $MKL -fopenmp -lpthread -lm -ldl || exit 2
"$OUT/elgeomtest"
