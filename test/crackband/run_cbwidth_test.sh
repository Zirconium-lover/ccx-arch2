#!/bin/bash
# Unit test of the CB1 crack-band width against hand-computable geometry.
#
#   test/crackband/run_cbwidth_test.sh [src-dir]
#
# Links the driver against the built ccx archive and calls damcbwidth
# directly, so the geometry is checked without a finite element run.  Every
# case prints the measured width beside the expected one; the expected
# values are derived in the driver's comments, not fitted to the output.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
SRC=${1:-$HERE/../../src}
AR="$SRC/ccx_2.23.a"
[ -f "$AR" ] || { echo "no $AR - build the tree first" >&2; exit 2; }
BIN=$(mktemp -d)/cbwidth_test
gfortran -O2 -o "$BIN" "$HERE/cbwidth_test.f" "$AR" -lpthread -lm -ldl || exit 2
"$BIN"
