#!/bin/bash
# Prepare a Claude Code on the web container for this repository.
#
# Three things have to be true before any work here means anything:
# the Fortran/C toolchain and the solver libraries are present, the
# generated switch registry is not stale, and the tree actually builds.
# The third is what turns the first two from a claim into a check.
set -uo pipefail
# No -e: apt is allowed to warn about repositories that do not matter here;
# each step below decides for itself what counts as a failure.

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
  # The image carries third-party PPAs (deadsnakes, ondrej/php) the proxy
  # refuses with 403, and that fails apt-get update outright rather than
  # warning.  They have nothing to do with this build, so they are moved
  # aside.  universe and multiverse are already enabled in ubuntu.sources,
  # which is where intel-mkl lives, so add-apt-repository is not needed -
  # and it was the fragile part, because it talks to launchpad too.
  for f in /etc/apt/sources.list.d/*.sources /etc/apt/sources.list.d/*.list; do
    [ -f "$f" ] || continue
    grep -qs "ppa.launchpadcontent.net" "$f" && mv "$f" "$f.disabled"
  done

  export DEBIAN_FRONTEND=noninteractive
  apt-get update -qq || echo "apt-get update reported problems; continuing"

  echo "installing:$missing"
  log=$(mktemp)
  # shellcheck disable=SC2086
  if apt-get install -y -qq $missing >"$log" 2>&1; then
    echo "installed:$missing"
    rm -f "$log"
  else
    echo "INSTALL FAILED for:$missing"
    cat "$log"
    rm -f "$log"
    exit 1
  fi
else
  echo "toolchain already present"
fi

# The switch registry and its documentation are generated from the sources.
# The regression gate's preflight fails when they are stale, and the file
# is not kept in the repository, so it has to be made here.
if [ -f tools/mkswitches.py ]; then
  if python3 tools/mkswitches.py >/dev/null; then
    echo "switch registry regenerated"
  else
    echo "FAILED to regenerate the switch registry"
    exit 1
  fi
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
