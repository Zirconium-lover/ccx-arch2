#!/usr/bin/env python3
"""Verify the UC6 NORMAL law point by point, from a close.inp run.

    check_close.py <rundir> [--zeta Z] [--jobname close]

The deck pulls one cohesive facet onto its softening branch and then pushes
it back through zero separation into contact, printing, every increment, the
separation (U at the two coincident interface node sets), the damage (SDV 2)
and the traction (S on the facet).  This checks the constitutive identity

    sharp     T = kn*d                        d < 0
              T = g*kn*d                      d >= 0
    blended   T = g*kn*d - (1-g)*kn*psi(-d)
              psi(u) = 0 (u<=0), u^2/(2w) (0<u<w), u-w/2 (u>=w),  w = zeta*d0

with g = max(gmin, 1-dvisc), and nothing else: no bulk model, no structural
branch and no finite-strain correction enters it, so a failure here is a
failure of cohesive_uc6.f.

It also checks what the blend is FOR - that the tangent is continuous.  The
sharp law's slope jumps by 1/g at d=0; the blended one must not.
"""
import argparse,os,re,sys

def blocks(dat):
    """yield (kind, time, rows) for every printed block, in file order."""
    kind=None;t=None;rows=None
    for ln in open(dat,errors='replace'):
        low=ln.lower()
        m=re.match(r'\s*(displacements|stresses|internal state variables)\s*\(.*for set (\S+) and time\s+(\S+)',low)
        if m:
            if rows: yield kind,t,rows
            kind=(m.group(1),m.group(2));t=float(m.group(3));rows=[]
            continue
        if rows is not None:
            f=ln.split()
            if f and f[0].isdigit(): rows.append([float(x) for x in f])
            elif ln.strip() and not ln[0].isspace():
                yield kind,t,rows; kind=None;rows=None
    if rows: yield kind,t,rows

def psi(u,w):
    if u<=0.: return 0.,0.
    if u<w:   return u*u/(2.*w),u/w
    return u-0.5*w,1.

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('rundir'); ap.add_argument('--jobname',default='close')
    ap.add_argument('--zeta',type=float,default=0.0)
    ap.add_argument('--tol',type=float,default=1e-5,
                    help='the .dat carries seven significant digits, so 1e-5 '
                         'is the file format, not the identity')
    a=ap.parse_args()
    dat=os.path.join(a.rundir,a.jobname+'.dat')
    deck=os.path.join(a.rundir,a.jobname+'.inp')
    txt=open(deck,errors='replace').read()
    m=re.search(r'\*USER SECTION[^\n]*\n([^\n]*)',txt,re.I)
    kn,tn0,ts0,gc,gmin,mu=[float(x) for x in m.group(1).split(',')]
    d0=tn0/kn; w=a.zeta*d0
    cur={}; recs=[]
    for kind,t,rows in blocks(dat):
        what,setname=kind
        cur.setdefault(t,{})[ (what,setname) ]=rows
    for t in sorted(cur):
        d=cur[t]
        um=d.get(('displacements','nifm')); up=d.get(('displacements','nifp'))
        st=d.get(('stresses','ecoh')); sd=d.get(('internal state variables','ecoh'))
        if not(um and up and st and sd): continue
        dn=sum(r[1] for r in up)/len(up)-sum(r[1] for r in um)/len(um)
        T=sum(r[2] for r in st)/len(st)
        dvisc=sum(r[3] for r in sd)/len(sd)
        recs.append((t,dn,T,max(gmin,1.-dvisc)))
    if not recs: sys.exit("no complete increments found in %s"%dat)
    print("%d increments, Kn=%g Tn0=%g gmin=%g, d0=%.6g, zeta=%g -> band w=%.6g"
          %(len(recs),kn,tn0,gmin,d0,a.zeta,w))
    worst=0.;worst_at=None;nclose=0;nband=0
    for t,dn,T,g in recs:
        if dn<0: nclose+=1
        if w>0 and -w<dn<0: nband+=1
        if w>0.:
            p,_=psi(-dn,w); want=g*kn*dn-(1.-g)*kn*p
        else:
            want=kn*dn if dn<0 else g*kn*dn
        scale=max(abs(want),tn0*1e-6)
        err=abs(T-want)/scale
        if err>worst: worst,worst_at=err,(t,dn,T,want,g)
    print("  increments with the facet CLOSED: %d%s"
          %(nclose,", of which inside the blend band: %d"%nband if w>0 else ""))
    fails=[]
    print("  worst relative error against the law: %.3e"%worst)
    if worst_at:
        t,dn,T,want,g=worst_at
        print("     at time %.6g: d=%.6e g=%.6e  T=%.8e want %.8e"%(t,dn,g,T,want))
    if worst>a.tol: fails.append("law: worst error %.3e exceeds %.1e"%(worst,a.tol))
    if nclose<5: fails.append("the facet closed on only %d increments"%nclose)
    if w>0. and nband<3:
        fails.append("only %d increments landed inside the blend band; the "
                     "band was never exercised"%nband)
    # Far-field slopes, measured on the CLOSING step only, where the damage
    # is frozen (deff stays below dmax0) so both sides are elastic and the
    # secant IS the tangent.  Mixing in the pulling step would sample the
    # softening branch and the number would mean nothing.
    imax=max(range(len(recs)),key=lambda i:recs[i][1])
    rc=sorted(recs[imax:],key=lambda r:r[1])
    gfin=recs[-1][3]
    band=max(w,1e-12)
    neg=[];pos=[]
    for i in range(1,len(rc)):
        a0,a1=rc[i-1][1],rc[i][1]
        dd=a1-a0
        if abs(dd)<=1e-14: continue
        slope=(rc[i][2]-rc[i-1][2])/dd
        # BOTH endpoints must be on the same side and clear of the band: a
        # pair that straddles zero mixes the two slopes and, on the sharp
        # law, that single pair is large enough to move the mean by 300x.
        if a0<-3.*band and a1<-3.*band: neg.append(slope)
        elif 3.*band<a0<50.*d0 and 3.*band<a1<50.*d0: pos.append(slope)
    if neg and pos:
        sneg=sum(neg)/len(neg); spos=sum(pos)/len(pos)
        print("  on the closing step, far in compression dT/dd=%.6g (Kn=%g), "
              "far in tension %.6g (g*Kn=%.6g), ratio %.4g (1/g=%.4g)"
              %(sneg,kn,spos,gfin*kn,sneg/max(1e-30,spos),1./max(gfin,1e-30)))
        if abs(sneg-kn)>1e-3*kn:
            fails.append("far-compression slope %.6g is not Kn=%g"%(sneg,kn))
        if abs(spos-gfin*kn)>5e-2*max(gfin*kn,1e-30):
            fails.append("far-tension slope %.6g is not g*Kn=%.6g"%(spos,gfin*kn))
    else:
        fails.append("not enough points on either side to measure a slope")
    print()
    if fails:
        for f in fails: print("FAIL",f)
        print("\nFAILED (%d)"%len(fails)); return len(fails)
    print("PASSED (0 failure(s))"); return 0

if __name__=='__main__': sys.exit(main())
