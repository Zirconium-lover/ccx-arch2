#!/usr/bin/env bash
# Out-of-tree PARDISO/MKL build.  Sources are symlinked, objects are not
# shared with the SPOOLES build, so the two binaries can never be mixed.
#
#   ./src/build_mkl.sh [build-dir]      default: <repo>/build-mkl
set -eu
SRC=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
OUT=${1:-"$SRC/../build-mkl"}
mkdir -p "$OUT"
OUT=$(CDPATH= cd -- "$OUT" && pwd)
cd "$OUT"
for f in "$SRC"/*.c "$SRC"/*.f "$SRC"/*.h "$SRC"/Makefile.inc \
         "$SRC"/Makefile.ubuntu2404.mkl; do
    ln -sf "$f" .
done
make -f Makefile.ubuntu2404.mkl "-j${NPROC:-4}"
ls -l "$OUT/ccx_2.23_pardiso"
