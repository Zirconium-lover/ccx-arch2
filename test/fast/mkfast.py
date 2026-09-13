#!/usr/bin/env python3
"""Generate a small 3-D fracture specimen that reproduces the wall classes in
minutes instead of hours.

The target deck s3rad has 42807 elements and takes ~2.5 h to reach its walls,
so no refactor of the damage module can be validated cheaply.  The other
decks in this tree (cohesive.inp, mixed.inp, snapback.inp) are 1-D chains:
they exercise the cohesive law and the continuation, but a 1-D chain cannot
produce a node hanging by one tetrahedron, which is the mechanism behind the
second s3rad wall.

This builds the smallest thing that can: a 3-D bar of tetrahedra, split by a
cohesive plane, with a brittle band across it to localise the crack, pulled
in displacement control.  Fragments arise naturally in a tet mesh once
deletion starts.

    ./mkfast.py -o fast.inp [--nx 8 --ny 5 --nz 5]
"""
import argparse,sys

def build(nx,ny,nz,lx,ly,lz,seed,incl=None,bandback=0):
    """Bar of tetrahedra with a PARTIAL cohesive plane at mid-span.

    The plane is cohesive only over the first `seed` fraction of the y range;
    the rest stays continuous.  A full-section interface is the wrong model
    here - it fails first, the bulk unloads and never damages, and the run
    finishes at theta=1 with zero deletions (measured).  A partial one is a
    pre-crack: the remaining ligament has to tear through the bulk, which is
    what produces deletion, fragments and the wall classes."""
    dx,dy,dz=lx/nx,ly/ny,lz/nz
    mid=nx//2
    jseed=max(1,int(round(seed*ny)))          # cohesive rows 0..jseed-1
    def in_seed(j): return j<=jseed           # node row index
    nid={}; nodes=[]
    def add_raw(i,j,k,side):
        key=(i,j,k,side)
        if key in nid:
            return nid[key]
        nid[key]=len(nodes)+1
        nodes.append((i*dx,j*dy,k*dz))
        return nid[key]
    def add(i,j,k,side):
        # only nodes inside the seeded part of the mid plane are duplicated
        if not(i==mid and in_seed(j)): side=0
        return add_raw(i,j,k,side)
    # --- inclusion -----------------------------------------------------
    # A box of cells made of the eroding phase (ZRH), with EVERY one of its
    # nodes duplicated and its whole boundary wrapped in cohesive facets.
    # That is the arrangement the target deck has and this generator did not:
    # in s3rad the facets wrap the phase that erodes, so when it erodes the
    # nodes on the wrap keep no bulk support at all and their assembled
    # diagonal falls to the cohesive residual.  A plane that merely cuts the
    # section, as below, can never isolate anything.
    ibox=None
    if incl is not None:
        i0,j0,k0,ni,nj,nk=incl
        i1,j1,k1=i0+ni,j0+nj,k0+nk
        if not(0<=i0<i1<=nx and 0<=j0<j1<=ny and 0<=k0<k1<=nz):
            raise SystemExit("inclusion box is outside the bar")
        if i0<=mid<=i1:
            raise SystemExit("inclusion node range must not touch the mid "
                             "plane i=%d, whose nodes are already duplicated"%mid)
        ibox=(i0,j0,k0,i1,j1,k1)
    def in_box(i,j,k):
        if ibox is None: return False
        i0,j0,k0,i1,j1,k1=ibox
        return i0<=i<i1 and j0<=j<j1 and k0<=k<k1
    def addI(i,j,k): return add_raw(i,j,k,'I')
    def inside(i,j,k): return 0<=i<nx and 0<=j<ny and 0<=k<nz
    def N(i,j,k,elem_i):
        side = 0 if elem_i<mid else 1
        return add(i,j,k,side)
    TETS=[(0,1,3,7),(0,1,7,5),(0,5,7,4),(0,3,2,7),(0,6,4,7),(0,2,6,7)]
    bulk=[]
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                if in_box(i,j,k):
                    c=[addI(i+(m&1),j+((m>>1)&1),k+((m>>2)&1)) for m in range(8)]
                    mat=2
                else:
                    c=[N(i+(m&1),j+((m>>1)&1),k+((m>>2)&1),i) for m in range(8)]
                    brittle = (jseed-bandback <= j < jseed+max(1,ny//4)) and (abs(i-mid)<=1)
                    mat=2 if brittle else 1
                for t in TETS:
                    bulk.append(([c[t[0]],c[t[1]],c[t[2]],c[t[3]]],mat))
    coh=[]
    # wrap the inclusion.  For direction d the quad is chosen so that the
    # first (matrix) triangle's normal points INTO the box, i.e. from the
    # first triangle towards the second, which is the orientation the mid
    # plane below uses.
    WRAP=(((-1,0,0),(0,2,6,4)),((1,0,0),(1,5,7,3)),
          ((0,-1,0),(0,4,5,1)),((0,1,0),(2,3,7,6)),
          ((0,0,-1),(0,1,3,2)),((0,0,1),(4,6,7,5)))
    cohA=[]
    if ibox is not None:
        i0,j0,k0,i1,j1,k1=ibox
        for i in range(i0,i1):
            for j in range(j0,j1):
                for k in range(k0,k1):
                    for (d,q) in WRAP:
                        nb=(i+d[0],j+d[1],k+d[2])
                        if in_box(*nb): continue
                        # a box face on the outer surface of the bar has no
                        # matrix on the other side: the inclusion is simply
                        # exposed there and no facet exists.
                        if not inside(*nb): continue
                        cm=[(i+(m&1),j+((m>>1)&1),k+((m>>2)&1)) for m in q]
                        o=[add(*p,0) for p in cm]
                        n=[addI(*p) for p in cm]
                        cohA.append([o[0],o[1],o[2], n[0],n[1],n[2]])
                        cohA.append([o[0],o[2],o[3], n[0],n[2],n[3]])
    for j in range(ny):
        if not(in_seed(j) and in_seed(j+1)): continue
        for k in range(nz):
            l=[add(mid,j+a,k+b,0) for a,b in ((0,0),(1,0),(1,1),(0,1))]
            r=[add(mid,j+a,k+b,1) for a,b in ((0,0),(1,0),(1,1),(0,1))]
            coh.append([l[0],l[1],l[2], r[0],r[1],r[2]])
            coh.append([l[0],l[2],l[3], r[0],r[2],r[3]])
    return nodes,bulk,coh,cohA,nid,mid

def write(path,nodes,bulk,coh,cohA,nid,nx,ny,nz,mid,lx,materials,anchor):
    matrix=[i+1 for i,(n,m) in enumerate(bulk) if m==1]
    plate =[i+1 for i,(n,m) in enumerate(bulk) if m==2]
    off=len(bulk)
    x0=sorted({v for (i,j,k,s),v in nid.items() if i==0})
    xl=sorted({v for (i,j,k,s),v in nid.items() if i==nx})
    def lines(vals,per=8):
        return "\n".join(", ".join(str(v) for v in vals[i:i+per]) for i in range(0,len(vals),per))
    with open(path,'w') as f:
        f.write("** Fast fracture regression specimen - see mkfast.py\n*Node\n")
        for i,(x,y,z) in enumerate(nodes): f.write("%d, %.6f, %.6f, %.6f\n"%(i+1,x,y,z))
        f.write("*Element, Type=C3D4\n")
        for i,(n,m) in enumerate(bulk): f.write("%d, %d, %d, %d, %d\n"%(i+1,*n))
        f.write("*User Element, Type=UC6, Nodes=6, Integration Points=3, MaxDof=3\n")
        f.write("*Element, Type=UC6\n")
        for i,n in enumerate(coh+cohA): f.write("%d, %d, %d, %d, %d, %d, %d\n"%(off+i+1,*n))
        f.write("*Elset, Elset=MATRIX\n%s\n"%lines(matrix))
        f.write("*Elset, Elset=PLATETANGENTIAL\n%s\n"%lines(plate))
        # Without --anchor the wrap is the SAME interface as the rest, so the
        # generator states no cohesive constant of its own and the physics
        # cannot drift from the deck these were lifted from.  That is also the
        # arrangement that was measured to work; see README.
        nreg=len(coh)+(0 if anchor else len(cohA))
        f.write("*Elset, Elset=INTERFACE\n%s\n"%lines([off+i+1 for i in range(nreg)]))
        if anchor and cohA:
            f.write("*Elset, Elset=INTERFACE_ANCHOR\n%s\n"
                    %lines([off+len(coh)+i+1 for i in range(len(cohA))]))
        f.write("*Nset, Nset=FACE_X0_NSET\n%s\n"%lines(x0))
        f.write("*Nset, Nset=FACE_XL_NSET\n%s\n"%lines(xl))
        f.write("*Nset, Nset=FIXPOINTA\n%d\n"%x0[0])
        f.write("*Nset, Nset=FIXPOINTB\n%d\n"%x0[-1])
        f.write(materials)
        if anchor and cohA:
            # --anchor gives the wrap its own, stronger constants.  It was
            # tried and it is NOT what makes the specimen work: a wrap that
            # outlives the phase it wraps leaves the collapsed node held by
            # INTACT facets, which floors the diagonal ratio at ~2.5e-2 and
            # never reaches 1.e-3 (README, "what was tried").  Kept because it
            # is the only way to ask that question again.
            f.write("** wrap around the eroding phase: Kn, Tn0, Ts0, Gc, "
                    "residual fraction, viscosity\n")
            f.write("*User Section, Elset=INTERFACE_ANCHOR, Material=COHESIVE,"
                    " Constants=6\n%s\n"%anchor)
        f.write("""*Step, Nlgeom, Inc=20000
*Static, Solver=Pardiso
1.000000e-03, 1., 1.000000e-09, 2.000000e-03
*Boundary
FACE_X0_NSET, 1, 1, 0.
FACE_XL_NSET, 1, 1, %.4f
FIXPOINTA, 2, 3, 0.
FIXPOINTB, 3, 3, 0.
*Output, Frequency=50
*Node File
U, RF
*El File
S, E, PEEQ, SDV
*Node Print, Nset=FACE_XL_NSET, Totals=Yes, Frequency=50
RF, U
*End Step
"""%(0.25*lx))

def lift_materials(src):
    """Copy the material and section cards verbatim from the target deck, so
    this specimen can never drift from the physics it is meant to stand in
    for.  Only the element-set names on the *Solid/User Section cards are
    reused; the mesh above defines those sets itself."""
    txt=open(src,errors='replace').read()
    a=txt.index('*Material, Name=ZR')
    b=txt.index('*Step,')
    block=txt[a:b]
    # this specimen has a single INTERFACE set; the target deck has two, so
    # drop the seed section entirely (card plus its constants line) and point
    # the regular one at our set.
    out=[];skip=0
    for ln in block.split('\n'):
        if skip: skip-=1; continue
        if 'Elset=INTERFACE_SEED' in ln: skip=1; continue
        out.append(ln.replace('Elset=INTERFACE_REGULAR','Elset=INTERFACE'))
    return '\n'.join(out)

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('-o',required=True); ap.add_argument('--nx',type=int,default=8)
    ap.add_argument('--ny',type=int,default=5); ap.add_argument('--nz',type=int,default=5)
    ap.add_argument('--lx',type=float,default=4.0); ap.add_argument('--ly',type=float,default=2.0)
    ap.add_argument('--lz',type=float,default=2.0)
    ap.add_argument('--from-deck',default='test/s3rad/m12_s3rad_gc24_w.inp')
    ap.add_argument('--seed',type=float,default=0.5,help='fraction of the mid plane that is a pre-crack')
    ap.add_argument('--band-back',type=int,default=0,dest='bandback',
                    help='extend the brittle band this many cell rows BEHIND '
                         'the crack tip, so the eroding phase straddles the '
                         'tip instead of sitting only ahead of it')
    ap.add_argument('--incl',default=None,
                    help='ZRH inclusion wrapped in cohesive facets, as '
                         'i0,j0,k0,ni,nj,nk in CELL indices, or "auto" to '
                         'place a 2x1x2 box just ahead of the crack tip')
    ap.add_argument('--anchor',default=None,
                    help='give the wrap its own cohesive constants instead of '
                         'the interface the deck already defines; see README '
                         'for why the default is not to')
    a=ap.parse_args()
    incl=None
    if a.incl=='auto':
        mid=a.nx//2; jseed=max(1,int(round(a.seed*a.ny)))
        incl=(mid+1,jseed,max(1,a.nz//2-1),2,1,2)
    elif a.incl:
        incl=tuple(int(v) for v in a.incl.replace(',',' ').split())
        if len(incl)!=6: raise SystemExit('--incl needs i0,j0,k0,ni,nj,nk')
    nodes,bulk,coh,cohA,nid,mid=build(a.nx,a.ny,a.nz,a.lx,a.ly,a.lz,a.seed,incl,a.bandback)
    mats=lift_materials(a.from_deck)
    write(a.o,nodes,bulk,coh,cohA,nid,a.nx,a.ny,a.nz,mid,a.lx,mats,a.anchor)
    print("wrote %s: %d nodes, %d C3D4, %d UC6 (%d of them wrapping the "
          "inclusion %s)"%(a.o,len(nodes),len(bulk),len(coh)+len(cohA),
                           len(cohA),incl))
