#!/usr/bin/env bash
# Does the nonlocal model give MESH OBJECTIVITY - the property it exists for?
#
#   CCX_EXE=/path/to/ccx_2.23_pardiso test/crackband/run_nlobjectivity.sh <dir>
#
# WHY THIS SWEEP AND NOT run_objectivity.sh.  That one refines the CROSS
# SECTION, which holds the band one slice thick on every mesh - right for a
# crack-band width, and useless here, because a band whose thickness the deck
# fixes asks an internal length nothing at all.  This refines ALONG the axis
# with the specimen fixed, so the band must choose its own width, and that
# choice is what the internal length is supposed to govern.
#
# WHY A NOTCH AND NOT A WEAK SLICE.  A weak slice pins the band to itself.
# Peerlings et al. 1996 section 5 keeps the material uniform and reduces the
# cross-sectional area over a central segment; the band then forms where the
# stress is highest but sizes itself.
#
# WHY THE SAME VISCOSITY ON EVERY ARM.  Without it the ell=0.5 arm at nslice=12
# dies at u=0.062, before the peak, and the comparison is then made in the
# elastic range - a difference of STAGE masquerading as a difference of model.
# Raising it only on the nonlocal arms would introduce a second variable;
# raising it on ALL arms leaves one.  That distinction is the whole reason
# this sweep can be run at all.
#
# WHY ell=0.5.  2*ell = 1.0 spans 1, 2 and 4 slices on the three meshes, so
# the coarsest CANNOT resolve the band and is reported separately rather than
# averaged in.
set -u
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
OUT=${1:-$(mktemp -d)}
EXE=${CCX_EXE:?set CCX_EXE to the ccx binary}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1} MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
mkdir -p "$OUT"
ELL=${ELL:-0.5}
for ns in 6 12 24; do
  for arm in local nl; do
    d="$OUT/${ns}_${arm}"; rm -rf "$d"; mkdir -p "$d"
    python3 "$ROOT/test/crackband/mkcross.py" --nslice "$ns" --ltot 6.0 \
        --ncross 1 --alldamage --trigger 1.0 --notch 0.10 --uend 0.30 \
        -o "$d/t.inp" >/dev/null || exit 2
    [ "$arm" = nl ] && sed -i "s/EVOLUTION=DISPLACEMENT\$/EVOLUTION=DISPLACEMENT, NONLOCAL=$ELL/g" "$d/t.inp"
    ( cd "$d" && env CCX_DAMAGE_CHARLEN=1 CCX_DAMAGE_VISCOSITY=1.e-3 "$EXE" t > run.log 2>&1 )
    printf "  nslice=%-3s %-6s rc=%-4s theta=%s\n" "$ns" "$arm" "$?" \
           "$(awk '{t=$3}END{print t}' "$d/t.sta" 2>/dev/null)"
  done
done
echo
echo "  band width, elements with D>0.5, in LENGTH - the quantity the"
echo "  internal length is supposed to hold fixed:"
python3 - "$OUT" <<'PY'
import sys
out=sys.argv[1]
for arm in ("local","nl"):
    row=[]
    for ns in (6,12,24):
        best=0
        try:
            for i,l in enumerate(open("%s/%d_%s/t.de1stats"%(out,ns,arm))):
                if i<2 or l.startswith('#'): continue
                f=l.split()
                if len(f)>8: best=max(best,int(f[7]))
        except Exception: pass
        row.append(best/6.0*(6.0/ns))
    print("    %-6s %s" % (arm, "  ".join("%.3f"%v for v in row)))
PY
echo
for arm in local nl; do
  echo "  === $arm, dissipated work ==="
  python3 "$ROOT/test/crackband/check_objectivity.py" 0.30 \
     n6:1.0:"$OUT/6_$arm/t.dat" n12:1.0:"$OUT/12_$arm/t.dat" \
     n24:1.0:"$OUT/24_$arm/t.dat" 2>/dev/null | sed -n '1,6p'
done
