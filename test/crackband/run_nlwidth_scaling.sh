#!/usr/bin/env bash
# Does the internal length SET the band width, and does the dissipation
# follow the width it sets?
#
#   CCX_EXE=/path/to/ccx_2.23_pardiso test/crackband/run_nlwidth_scaling.sh <dir>
#
# WHY THIS AND NOT run_nlobjectivity.sh.  That sweep refines the mesh at a
# fixed ell and asks whether the answer stops moving.  It shows that the
# width stops moving and that the dissipated work does not, which is a
# result and not an explanation: holding one quantity fixed while the other
# drifts could be the averaging, the softening law, or the interaction of
# the two, and a refinement sweep cannot separate them because every
# candidate changes together with h.
#
# Here h is FIXED and only ell moves.  That makes the question a scaling
# law rather than a convergence test, and the two candidates predict
# different laws:
#
#   if the internal length sets the width          width  ~ ell
#   if each element inside the band still charges
#   the full G_f over its OWN projected size       W_post ~ ell / h
#
# The second is the double-count to look for.  The crack-band law advances
# damage as D += L*dEps_p/u_f with L the element's own width, so an element
# dissipates G_f = sigma_0*u_f/2 per unit area whatever its size - that is
# the property Bazant and Oh 1983 build it for, and it is correct exactly
# when the band is ONE element wide.  A nonlocal length that widens the
# band to n element layers then charges n*G_f for one crack, and the
# fracture energy of the model is no longer a material constant but a
# function of ell/h.  Jirasek and Bauer 2012 section 5 state the
# requirement the other way round: the width entering the softening law
# must be the width of the band that actually forms.
#
# A dissipation that grows with ell at fixed h is therefore not a bug in
# the averaging - the averaging is doing what it is asked - but a statement
# that the two regularisations are being combined without reconciling their
# lengths.
set -u
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
OUT=${1:-$(mktemp -d)}
EXE=${CCX_EXE:?set CCX_EXE to the ccx binary}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1} MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
mkdir -p "$OUT"
VISC=${VISC:-1.e-3}
NS=${NS:-24}
# ell = 0 is the local arm, the reference every ratio below is taken against.
ELLS=${ELLS:-"0 0.25 0.5 1.0"}

for ell in $ELLS; do
  d="$OUT/e$ell"; rm -rf "$d"; mkdir -p "$d"
  python3 "$ROOT/test/crackband/mkcross.py" --nslice "$NS" --ltot 6.0 \
      --ncross 1 --alldamage --trigger 1.0 --notch 0.10 \
      --notchlen 0.33333333333333333 --uend 0.30 \
      -o "$d/t.inp" >/dev/null || exit 2
  [ "$ell" != 0 ] && sed -i \
      "s/EVOLUTION=DISPLACEMENT\$/EVOLUTION=DISPLACEMENT, NONLOCAL=$ell/g" \
      "$d/t.inp"
  ( cd "$d" && env CCX_DAMAGE_CHARLEN=1 CCX_DAMAGE_VISCOSITY="$VISC" \
        "$EXE" t > run.log 2>&1 )
  printf "  ell=%-5s rc=%-4s theta=%s\n" "$ell" "$?" \
         "$(awk '{t=$3}END{print t}' "$d/t.sta" 2>/dev/null)"
done
echo
echo "  nslice=$NS, element size h=$(python3 -c "print(6.0/$NS)"), viscosity $VISC"
python3 - "$OUT" "$NS" "$ELLS" <<'PY'
import sys,os,re
out,ns,ells=sys.argv[1],int(sys.argv[2]),sys.argv[3].split()
sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               "test","crackband"))
h=6.0/ns

def width(d):
    """Largest simultaneously damaged volume, as a length along the bar.

    Taken as the maximum over the accepted history and not the final state,
    because an element reaching D=1 is deleted and drops out of the counts:
    the final state of a severed bar reports almost no damage at all.  Six
    Kuhn tetrahedra fill one slice, so count/6 is the number of slices and
    count/6*h is the length.
    """
    best=0
    for i,l in enumerate(open(os.path.join(d,"t.de1stats"))):
        if l.startswith('#'): continue
        f=l.split()
        if len(f)>8: best=max(best,int(f[7]))
    return best/6.0*h

def curve(path,uend=0.30):
    out,pend=[],None
    for line in open(path):
        m=re.search(r"total force .* and time\s+([-\dEe.+]+)",line)
        if m: pend=float(m.group(1)); continue
        if pend is not None:
            f=line.split()
            if len(f)==3:
                try: out.append((uend*pend,float(f[0])))
                except ValueError: pass
                pend=None
    return out

def work(c,umax):
    tot=0.0
    for i in range(1,len(c)):
        u0,s0=c[i-1]; u1,s1=c[i]
        if u0>=umax: break
        if u1>umax:
            s1=s0+(s1-s0)*(umax-u0)/(u1-u0); u1=umax
        tot+=0.5*(s0+s1)*(u1-u0)
    return tot

def gf(d):
    """G_f as the deck itself states it, not as this script assumes it."""
    for l in open(os.path.join(d,"t.inp")):
        if not l.startswith('**'): break
        m=re.search(r"G_f=([-\d.Ee+]+)",l)
        if m: return float(m.group(1))
    return None

rows=[]
for e in ells:
    d=os.path.join(out,"e"+e)
    try: c=curve(os.path.join(d,"t.dat"))
    except Exception: continue
    if not c: continue
    rows.append((float(e),width(d),c,gf(d)))
if len(rows)<2:
    print("  not enough arms completed to compare"); sys.exit(0)
usp=min(max(c,key=lambda p:p[1])[0] for _,_,c,_ in rows)
print("  common split at u=%.4f, the earliest peak of the set" % usp)
print("  ell      width    width/2ell   layers   W_pre     W_post   W_post/(loc)")
w0=None
for e,wd,c,_ in rows:
    wpre=work(c,usp); wpost=work(c,0.30)-wpre
    if w0 is None: w0=wpost
    print("  %-7s %7.4f  %10s %8.2f %8.4f %9.4f %9.3f"
          % (e,wd,("%.3f"%(wd/(2*e))) if e>0 else "   -",wd/h,wpre,wpost,
             wpost/w0))
print()

# THE PREDICTION, WITH NOTHING FITTED.  An element of size h reaches D=1
# when charlen*eps_p = u_f, so with charlen = h its plastic displacement at
# failure is u_f and it dissipates G_f per unit area - once.  A band that is
# w wide contains w/h such layers and therefore dissipates G_f*w/h, which is
# a material constant only when w = h.  Writing it out:
#
#     W_post ~ G_f * w / charlen        charlen = h  ->  G_f * w/h
#                                       charlen = w  ->  G_f, always
#
# so the same arithmetic that predicts the numbers below also says what
# would remove the dependence: when an internal length sets w, the length in
# the softening law must be w and not the element's own size.  G_f is taken
# from the deck header, and w and h are measured, so the column has no free
# parameter to absorb a disagreement.
g=rows[0][3]
if g:
    print("  W_post against G_f*w/h, PARAMETER FREE (G_f=%.4f from the deck):"
          % g)
    for e,wd,c,_ in rows:
        wpost=work(c,0.30)-work(c,usp)
        pred=g*wd/h
        print("    ell=%-6s w/h=%5.2f  predicted %9.4f  measured %9.4f"
              "  off %6.1f %%" % (e,wd/h,pred,wpost,
                                  100.0*(wpost-pred)/pred))
    print("    The layer count is a count of integration points over a")
    print("    threshold, so it is coarsest where the band is thinnest;")
    print("    read the arms with w/h >= 3 and treat the rest as the")
    print("    resolution of the measure rather than of the model.")
print()
print("  The two laws this set can tell apart, taken against the SMALLEST")
print("  nonlocal arm rather than the local one, because the local arm has")
print("  no ell to scale:")
nl=[r for r in rows if r[0]>0]
if len(nl)>=2:
    e0,wd0,c0,_=nl[0]
    p0=work(c0,0.30)-work(c0,usp)
    for e,wd,c,_ in nl[1:]:
        p=work(c,0.30)-work(c,usp)
        print("    ell %g -> %g : ell x%.2f   width x%.2f   W_post x%.2f"
              % (e0,e,e/e0,wd/wd0,p/p0))
    print("    width tracking ell means the averaging sets the width.")
    print("    W_post tracking ell too means each element layer inside the")
    print("    band still charges its own G_f, so the fracture energy of")
    print("    the model depends on ell/h - the lengths are not reconciled.")
PY
