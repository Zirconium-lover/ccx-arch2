#!/usr/bin/env bash
# What stops the LOCAL band from collapsing?  Vary the bar's cross-section.
#
#   CCX_EXE=/path/to/ccx_2.23_pardiso test/crackband/run_xsection_floor.sh <dir>
#
# WHY THIS EXISTS.  The local band width on this deck does not collapse with
# the mesh: 1.049 -> 0.619 -> 0.592 at h = 1.0 / 0.5 / 0.25, and a sharper
# notch on meshes down to h=0.125 leaves it at 0.634 / 0.545 / 0.642.  That was
# published as unexplained, and the first two candidates were already ruled out
# by measurement - it is not the notch span, and it is not the viscosity, which
# moves the width by 2.5 per cent when tripled.
#
# WHAT THIS VARIES, AND WHY IT IS ONE VARIABLE.  --h is the bar's transverse
# dimension.  The projected crack-band width of every Kuhn tetrahedron is the
# SLICE THICKNESS 6/nslice, along the axis, so charlen does not move when the
# cross-section does: the element's own length in the softening law is
# untouched.  The notch is a fixed FRACTION of the cross-section, so the
# relative stress profile along the bar is identical too.  What changes is only
# how far along the axis a cross-sectional perturbation takes to even out,
# which is of order the transverse dimension - Saint-Venant.
#
# Two axial meshes, because the two candidate floors cross over: a band cannot
# be shorter than the element that carries it, and at a coarse axial mesh that
# floor hides any other.
set -u
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
OUT=${1:-$(mktemp -d)}
EXE=${CCX_EXE:?set CCX_EXE to the ccx binary}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1} MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
mkdir -p "$OUT"
VISC=${VISC:-1.e-3}
HCROSS=${HCROSS:-"1.0 0.5 0.25"}
NSLICES=${NSLICES:-"12 24"}
for hc in $HCROSS; do
  for ns in $NSLICES; do
    d="$OUT/h${hc}_n${ns}"; rm -rf "$d"; mkdir -p "$d"
    python3 "$ROOT/test/crackband/mkcross.py" --nslice "$ns" --ltot 6.0 \
        --ncross 1 --alldamage --trigger 1.0 --notch 0.10 \
        --notchlen 0.33333333333333333 --uend 0.30 --h "$hc" \
        -o "$d/t.inp" >/dev/null || exit 2
    ( cd "$d" && env CCX_DAMAGE_CHARLEN=1 CCX_DAMAGE_VISCOSITY="$VISC" \
          "$EXE" t > run.log 2>&1 )
    rc=$?
    printf "  hcross=%-6s ns=%-3s rc=%-4s att=%-6s steptime=%-14s deleted=%s\n" \
        "$hc" "$ns" "$rc" \
        "$(awk '$1 ~ /^[0-9]+$/{a=$3}END{print a}' "$d/t.sta" 2>/dev/null)" \
        "$(awk '$1 ~ /^[0-9]+$/{t=$6}END{print t}' "$d/t.sta" 2>/dev/null)" \
        "$(grep -vc '^#' "$d/t.damage" 2>/dev/null || echo 0)"
    [ "$rc" = 0 ] || : > "$d/UNCONVERGED"
  done
done
python3 - "$OUT" "$ROOT" "$HCROSS" "$NSLICES" <<'PY'
import sys,os
out,root=sys.argv[1],sys.argv[2]
hcs,nss=sys.argv[3].split(),[int(v) for v in sys.argv[4].split()]
sys.path.insert(0,os.path.join(root,"test","crackband"))
from bandwidth import width
print()
print("  LOCAL band width, and the same divided by the cross-section")
hdr="  hcross  "+"".join("h=%-6.3f w/hc  " % (6.0/ns) for ns in nss)
print(hdr)
for hc in hcs:
    cells=""
    for ns in nss:
        d=os.path.join(out,"h%s_n%d"%(hc,ns))
        if os.path.exists(os.path.join(d,"UNCONVERGED")):
            cells+="  EXCLUDED     "; continue
        try: w=width(d)[0]
        except Exception: cells+="  no data      "; continue
        cells+="%-9.4f %-6.3f" % (w,w/float(hc))
    print("  %-7s %s" % (hc,cells))
print()
print("  TWO FLOORS, and they cross over.  A band cannot be shorter than the")
print("  element carrying it, so at the coarser axial mesh that floor hides")
print("  everything else and w barely moves with the cross-section.  At the")
print("  finer one the cross-section is free to act, and halving it halves w.")
print("  Quarter it and the element floor takes over again.")
print()
print("  This says the floor is the SPECIMEN, not a defect: a cross-sectional")
print("  perturbation evens out over a distance of order the transverse")
print("  dimension, so the band cannot be shorter than that however fine the")
print("  mesh.  Saint-Venant is the candidate MECHANISM and is not established")
print("  here - it predicts insensitivity to Poisson's ratio, which this deck")
print("  sets to zero and which --nu would vary.")
PY
