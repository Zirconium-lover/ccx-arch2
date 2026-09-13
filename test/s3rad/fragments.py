#!/usr/bin/env python3
"""Connected components of the LIVE model, weighted by real facet state.

    fragments.py <rundir> [--deck DECK] [--grips SET:SET]

Answers two questions that the per-node diagonal cannot:

1. **Is anything adrift?**  A component of live bulk elements that reaches no
   grip has six rigid-body modes with no stiffness.  Its nodes still look
   healthy to a diagonal test - their own elements are alive - so AUTOSPC is
   blind to it by construction.  This finds it by topology instead.

2. **Is the specimen actually in two pieces?**  Count a cohesive facet as a
   load path only when it still retains more than a given fraction of its
   stiffness, and see at which fraction the model falls apart.  That is an
   independent check on the severance question, using connectivity rather
   than the reaction force.

Facet state comes from the `*El Print, Elset=INTERFACE` SDV block in the .dat
file (state variable 2 is the viscous damage, so g = max(gmin, 1-dvisc)); the
deleted bulk comes from m.damage.  Nothing is re-run.
"""
import argparse,collections,os,sys

def read_deck(path):
    mode=None;c3={};uc={};nsets=collections.defaultdict(set);cur=None
    for ln in open(path,errors='replace'):
        t=ln.strip()
        if t.startswith('**'): continue
        if t.startswith('*'):
            low=t.lower().replace(' ','')
            if low.startswith('*element,type=c3d4'): mode='c3'
            elif low.startswith('*element,type=uc6'): mode='uc'
            elif low.startswith('*nset'): mode='ns';cur=low.split('nset=')[1].split(',')[0]
            else: mode=None
            continue
        if mode=='c3':
            v=[int(x) for x in t.replace(',',' ').split()]; c3[v[0]]=v[1:5]
        elif mode=='uc':
            v=[int(x) for x in t.replace(',',' ').split()]; uc[v[0]]=v[1:7]
        elif mode=='ns' and cur:
            nsets[cur]|={int(x) for x in t.replace(',',' ').split() if x}
    return c3,uc,nsets

def read_facet_state(dat,gmin=1.e-5):
    """last INTERFACE SDV block; returns {element: best g over its points}."""
    rows=None;when=None;last_when=None
    for ln in open(dat,errors='replace'):
        low=ln.lower()
        if 'internal state variables' in low:
            rows=[];last_when=low.split('and time')[-1].strip();continue
        if rows is not None and ln.startswith(' '):
            f=ln.split()
            if len(f)>=6 and f[0].isdigit(): rows.append((int(f[0]),float(f[3])))
        if rows is not None and ln.strip() and not ln[0].isspace(): rows=None
    g=collections.defaultdict(float)
    for e,dv in (rows or []): g[e]=max(g[e],max(gmin,1.0-dv))
    return g,last_when

def components(c3,uc,dele,gmax,grips,gthr):
    par={}
    def find(x):
        while par.setdefault(x,x)!=x: par[x]=par[par[x]]; x=par[x]
        return x
    def uni(a,b):
        ra,rb=find(a),find(b)
        if ra!=rb: par[ra]=rb
    live=[e for e in c3 if e not in dele]
    for e in live:
        n=c3[e]
        for x in n[1:]: uni(n[0],x)
    for e,n in uc.items():
        if e in dele: continue
        if gmax.get(e,1.0)<=gthr: continue
        for i in range(3): uni(n[i],n[i+3])
    comp=collections.defaultdict(list)
    for e in live: comp[find(c3[e][0])].append(e)
    anchored=[];adrift=[]
    for r,els in comp.items():
        ns=set()
        for e in els: ns|=set(c3[e])
        (anchored if (ns&grips) else adrift).append(len(els))
    return sorted(anchored,reverse=True),sorted(adrift,reverse=True),len(live)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('rundir')
    ap.add_argument('--deck',default=None,help='defaults to <rundir>/m.inp')
    ap.add_argument('--grips',default='FACE_X0_NSET:FACE_XL_NSET')
    ap.add_argument('--gmin',type=float,default=1.e-5)
    a=ap.parse_args()
    deck=a.deck or os.path.join(a.rundir,'m.inp')
    c3,uc,nsets=read_deck(deck)
    gmax,when=read_facet_state(os.path.join(a.rundir,'m.dat'),a.gmin)
    dele=set()
    dpath=os.path.join(a.rundir,'m.damage')
    if os.path.exists(dpath):
        dele={int(l.split()[0]) for l in open(dpath) if not l.startswith('#') and l.strip()}
    grips=set()
    for s in a.grips.split(':'): grips|=nsets.get(s.lower(),set())
    if not grips: sys.exit("no grip nodes found for %s"%a.grips)
    print("deck   %s: %d C3D4, %d UC6"%(deck,len(c3),len(uc)))
    print("state  %d facets, from the last printed increment (time %s)"%(len(gmax),when))
    print("       %d bulk elements deleted, %d live"%(len(dele&set(c3)),len(c3)-len(dele&set(c3))))
    h=collections.Counter()
    for e,g in gmax.items():
        h['>0.5' if g>0.5 else ('>1e-2' if g>1.e-2 else ('>floor' if g>2*a.gmin else 'floor'))]+=1
    print("       facet stiffness: %d above 0.5, %d in (1e-2,0.5], %d just above the floor, %d AT the residual floor"
          %(h['>0.5'],h['>1e-2'],h['>floor'],h['floor']))
    print()
    print("  a facet counts as a load path only above g =")
    print("  %-10s %-12s %-10s %-22s %s"%("g","components","anchored","adrift (components)","elements adrift"))
    for thr in (0.,2*a.gmin,1.e-3,1.e-2,0.05,0.1,0.25,0.5,0.75):
        anch,adr,live=components(c3,uc,dele,gmax,grips,thr)
        print("  %-10g %-12d %-10d %-22d %d (%.2f%%)"
              %(thr,len(anch)+len(adr),len(anch),len(adr),sum(adr),100.*sum(adr)/max(1,live)))
    print()
    print("adrift>0 at a threshold near the residual floor means a piece of the")
    print("model is held only by failed facets: healthy diagonals, no load path.")
    return 0

if __name__=='__main__': sys.exit(main())
