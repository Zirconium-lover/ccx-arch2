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
# WHY notchlen=1/3 AND A TENT.  The narrowing is imposed on nodes, so a step
# over a short segment is sampled differently by each mesh - over the middle
# 0.1 of this bar, nslice 6 and 12 catch one node and nslice 24 catches
# three, and the imperfection sharpens as the mesh is refined.  The sweep
# would then vary the specimen as well as the element size, and the band
# narrowing with refinement would be partly geometry with no way to divide
# it.  A tent over the middle 1/3 has its kinks at x = 2, 3 and 4, which are
# nodes at nslice 6, 12 and 24 alike, so all three decks describe the SAME
# solid to the last digit.  The preflight below asserts that rather than
# trusting it.
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
VISC=${VISC:-1.e-3}

# PREFLIGHT: one specimen on all three meshes.  Refining a mesh must not move
# the solid it discretises, and the only reason to state that as a check is
# that the obvious way of imposing the imperfection breaks it silently.  The
# check compares the COARSE deck's own surface, linearly interpolated, against
# every node of the finer ones, because node coincidence alone is not enough:
# a narrow notch can put its nodes in the right places and still be a taper of
# a different length on each mesh.  With --notchlen 0.1 this fails by 5.8e-2
# of the cross-section - which is what the first version of this sweep ran.
python3 - "$ROOT" "$OUT" <<'PY' || exit 3
import subprocess,sys,os
root,out=sys.argv[1],sys.argv[2]
nlen="0.33333333333333333"
def prof(ns):
    d=os.path.join(out,"pre%d"%ns); os.makedirs(d,exist_ok=True)
    subprocess.check_call(["python3",os.path.join(root,"test/crackband/mkcross.py"),
        "--nslice",str(ns),"--ltot","6.0","--ncross","1","--alldamage",
        "--trigger","1.0","--notch","0.10","--notchlen",nlen,
        "--uend","0.30","-o",os.path.join(d,"t.inp")],
        stdout=subprocess.DEVNULL)
    r,on={},False
    for l in open(os.path.join(d,"t.inp")):
        t=l.strip()
        if t.startswith('*'):
            on=t.upper().startswith('*NODE'); continue
        if on and t and not t.startswith('**'):
            f=[v.strip() for v in t.split(',')]
            if len(f)>=4:
                x=round(float(f[1]),9); y=float(f[2])
                r[x]=max(r.get(x,0.0),y)
    return r
def interp(r,x):
    xs=sorted(r)
    if x<=xs[0]: return r[xs[0]]
    for i in range(1,len(xs)):
        if xs[i]>=x:
            x0,x1=xs[i-1],xs[i]
            return r[x0]+(r[x1]-r[x0])*(x-x0)/(x1-x0)
    return r[xs[-1]]
base=prof(6); bad=0; worst=0.0
for ns in (12,24):
    fine=prof(ns)
    for x,y in sorted(fine.items()):
        d=abs(interp(base,x)-y)
        worst=max(worst,d)
        if d>1.0e-12:
            bad+=1
            if bad<=4: print("    PREFLIGHT FAIL nslice=%d x=%g: coarse"
                             " solid gives %.9f, this mesh %.9f"%(ns,x,interp(base,x),y))
if bad:
    print("  preflight FAILED: the three decks are not the same solid"
          " (%d of the finer nodes disagree, worst %.3e).  Refining the"
          " mesh moved the specimen, so the sweep would vary two things."
          % (bad,worst))
    sys.exit(1)
print("  preflight: one specimen on all three meshes (worst surface"
      " disagreement %.1e)" % worst)
PY
for ns in 6 12 24; do
  for arm in local nl; do
    d="$OUT/${ns}_${arm}"; rm -rf "$d"; mkdir -p "$d"
    python3 "$ROOT/test/crackband/mkcross.py" --nslice "$ns" --ltot 6.0 \
        --ncross 1 --alldamage --trigger 1.0 --notch 0.10 \
        --notchlen 0.33333333333333333 --uend 0.30 \
        -o "$d/t.inp" >/dev/null || exit 2
    [ "$arm" = nl ] && sed -i "s/EVOLUTION=DISPLACEMENT\$/EVOLUTION=DISPLACEMENT, NONLOCAL=$ELL/g" "$d/t.inp"
    ( cd "$d" && env CCX_DAMAGE_CHARLEN=1 CCX_DAMAGE_VISCOSITY="$VISC" "$EXE" t > run.log 2>&1 )
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
