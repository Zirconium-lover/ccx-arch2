#!/bin/bash
# Prepare a Claude Code on the web container for this repository.
#
# Three things have to be true before any work here means anything:
# the Fortran/C toolchain and the solver libraries are present, the
# generated switch registry is not stale, and the tree actually builds.
# The third is what turns the first two from a claim into a check.
set -euo pipefail

# Local machines are left alone; this is only for the remote container.
if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

cd "${CLAUDE_PROJECT_DIR:-$(dirname "$0")/../..}"

PKGS="gfortran libblas-dev liblapack-dev libarpack2-dev libspooles-dev
      intel-mkl make gcc python3"

missing=""
for p in $PKGS; do
  dpkg -s "$p" >/dev/null 2>&1 || missing="$missing $p"
done

if [ -n "$missing" ]; then
  echo "installing:$missing"
  export DEBIAN_FRONTEND=noninteractive
  # intel-mkl lives in universe/multiverse; enable them before asking for it.
  apt-get update -qq || true
  apt-get install -y -qq software-properties-common >/dev/null 2>&1 || true
  add-apt-repository -y universe    >/dev/null 2>&1 || true
  add-apt-repository -y multiverse  >/dev/null 2>&1 || true
  apt-get update -qq
  # shellcheck disable=SC2086
  apt-get install -y -qq $missing
else
  echo "toolchain already present"
fi

# The switch registry and its documentation are generated from the sources.
# The regression gate's preflight fails when they are stale, and the file
# is not kept in the repository, so it has to be made here.
if [ -f tools/mkswitches.py ]; then
  python3 tools/mkswitches.py >/dev/null
  echo "switch registry regenerated"
fi

# Build once, so the container is cached with a working binary and the
# toolchain is proven rather than assumed.  make is idempotent; a later
# source edit rebuilds only what it has to.
if [ -f src/Makefile.ubuntu2404.mkl ]; then
  # Stock CalculiX compiles with warnings; they are noise here and would
  # bury a real failure.  The whole build log is kept and printed only if
  # the build fails, so nothing is hidden that matters.
  log=$(mktemp)
  if make -C src -f Makefile.ubuntu2404.mkl -j"$(nproc)" >"$log" 2>&1; then
    echo "built src/ccx_2.23_pardiso"
  else
    echo "BUILD FAILED:"
    cat "$log"
    rm -f "$log"
    exit 1
  fi
  rm -f "$log"
fi
