!
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2025 Guido Dhondt
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     ==================================================================
!     Nonlocal (integral) averaging of the damage driving variable.
!
!     WHY THIS EXISTS.  Crack-band scaling makes the DISSIPATION per unit
!     crack area independent of element size and does that job correctly:
!     on three uniform bar meshes L*sigma_y/(u_f*E) scales exactly as L
!     and stays far below 1.  What it does NOT do is give the localisation
!     a WIDTH.  Once softening makes the tangent indefinite the boundary
!     value problem loses ellipticity and the band is as wide as the
!     discretisation makes it.  Measured (E-83): the load-displacement
!     curves of those three meshes agree to 0.94% up to the peak and
!     spread to 73-89% after it, and the post-peak curves are not even
!     ordered.
!
!     Viscosity does not fix it.  At eta=1.e-4 the same three meshes go
!     from failing at t=0.32-0.35 to running to t=1.0 while the curves
!     move by under a percent: solvability yes, objectivity no.  A rate
!     term cannot stand in for a length, and its help weakens under
!     refinement because the band localises into a smaller volume.
!
!     WHAT IT DOES.  The driving variable - the increment of equivalent
!     plastic strain - is replaced by its weighted average over a
!     neighbourhood of radius 2*ell:
!
!         dpeq_nl(e) = sum_f w(d) V_f dpeq(f) / sum_f w(d) V_f
!         w(d) = exp(-(d/ell)^2)
!
!     The element VOLUME sits in the weight, so cutting the same
!     neighbourhood into more elements does not change the average.
!     Leaving V out would smuggle the mesh dependence back in.
!
!     ell enters the formulation, so the band has a width of order ell on
!     any mesh and the problem is well posed.  This is the integral form
!     (Pijaudier-Cabot and Bazant).  The tangent is left LOCAL, which is
!     inconsistent and costs Newton iterations but avoids a dense
!     operator - the standard trade for this family.
!
!     THE MODULE IS DELIBERATELY PRIVATE TO THIS FILE.  Two places compute
!     damage and both must see the field: calcdamagebase, which nonlingeo
!     calls at the predictor and commit stages, and damageupdatepoint,
!     which resultsmech calls per integration point inside the stress
!     update.  The second is the one that drives the response - patching
!     only the first left the answer BIT-IDENTICAL to the local model,
!     which is how that was found.  Sharing through accessor subroutines
!     instead of a `use` in calcdamage.f keeps the module inside this
!     file, so no cross-file module dependency can break a parallel build.
!     ==================================================================
!
      module damnlmod
      implicit none
      real*8 :: ellsave=0.d0
      real*8, allocatable :: dpsave(:),cen(:,:),evol(:),wgt(:),dploc(:)
      integer, allocatable :: nbhead(:),nbnext(:),nblist(:),
     &     nbstart(:),nbcount(:)
      integer :: nbuilt=0,nesave=0
!
!     GRADIENT backend state.  imodenl: 0 = integral (the original), 1 =
!     implicit gradient.  gmass/gkpk/elnod cache the ELEMENT MATRICES of
!     the implicit-gradient equation - the lumped mass int N_a dV and the
!     conductivity int grad N_a . grad N_b dV - on the REFERENCE
!     configuration, for the same reason the neighbour list is built there:
!     the internal length is a material property and must not drift as the
!     mesh distorts.  nknl is the node count, derived from kon rather than
!     threaded through calcdamagebase and four call sites in nonlingeo.c.
!
!     These used to be hand-written linear-tetrahedron gradients in
!     bgrad(3,4,ne0), and every other family was skipped.  They now come
!     from damnlelmk, which uses ccx's own element library; see the
!     comment there for the two decisions that carries.
!
!     gkpk is PACKED: element i owns gnc(i)*(gnc(i)+1)/2 entries of the
!     upper triangle starting at gkst(i), because a tetrahedral mesh
!     should not pay a hexahedron's 36 slots for its 10.  The price is
!     that the cache is built in two passes - one to size it, one to fill
!     it - and damnlelmk is therefore called twice per element.  That is
!     once per ANALYSIS, not per iteration, and it buys back roughly a
!     factor of three on the dominant array for the mesh family the
!     module is most often used on.
!
      integer :: imodenl=0,gbuilt=0,nknl=0,fbuilt=0,nkpack=0
      integer :: ndnmax=0
!
!     TIME-CONSISTENT LAG.  dpsave is an INCREMENT of plastic strain,
!     averaged, measured over one increment of size dtref.  It is read by
!     damageupdatepoint on every Newton iteration of the NEXT attempt,
!     whose size is dtcur.  When the solver cuts the step, dtcur falls and
!     dpsave did not: the damage law was handed a full previous step's
!     damage increment at a step a thousand times smaller.  idtscale=1
!     hands out the lagged RATE times the current step instead, see
!     damnonlocalval.
!
      real*8 :: dtref=0.d0,dtcur=0.d0
      integer :: idtscale=1
      real*8 :: relmax=0.d0
      real*8, allocatable :: vele(:),ebar(:),
     &     gdiag(:),grhs(:),gp(:),gap(:),gz(:),gr(:),
     &     gmass(:,:),gkpk(:)
      integer, allocatable :: elnod(:,:),gnc(:),gkst(:)
      real*8, allocatable :: gmsum(:),ebar0(:)
!
!     PER-ELEMENT internal length.  ell2e(i) is what the two assembly
!     loops read; everything else in the backend is unchanged.  ellmf
!     holds a multiplier per material index and nlocn/nlocr the
!     localizing exponent and residual.  ellinit guards the one-time
!     environment read.
!
      real*8, allocatable :: ell2e(:)
      real*8 :: ellmf(16)=1.d0,nlocn=0.d0,nlocr=0.d0
      integer :: ellinit=0,nlocon=0,ellmn=0
!
!     W6.  ell as a MATERIAL constant, read from *DAMAGE INITIATION by
!     damnlellsetmat rather than from the environment.  ellcard(imat) is
!     the card's value, 0 meaning "this material did not give one"; ncard
!     counts how many did, so the run can say whether any card was seen
!     at all without scanning the array.
!
!     WHY A SETTER AND NOT dmcon.  There is no room in dmcon: nonlingeo.c
!     checks the constant count EXACTLY (nconst==4 for Rice-Tracey, and
!     (nconst-3) even for the tabulated locus), so one appended constant
!     turns damage_progressive_material into a silent no - the deck runs,
!     with the legacy hard-deletion model, and says nothing. Found by
!     Agent 2, confirmed here at nonlingeo.c:1002-1003 and at
!     calcdamage.f:1023, where npoints is derived from the same count.
!
!     WHY THAT IS SAFE HERE, WHEN THE SAME SHAPE WAS A BUG ELSEWHERE.
!     CB1's cbmode was module state too and was read before anything set
!     it, because its read was LAZY - the first reader triggered it, and
!     one reader ran before the writer ever did. This value is set by the
!     input reader, and input parsing finishes before any solver path
!     exists, so no reader can precede the write. The hazard is lazy
!     initialisation, not module state.
!
      real*8 :: ellcard(64)=0.d0
      integer :: ncard=0
!
!     Integration points per element, filled by the same pass that fills
!     evol.  A C3D4 has one and the old code hard-coded that; every other
!     family has more, and averaging only slot 1 would read one point and
!     call it the element.
!
      integer, allocatable :: nipe(:)
!
!     Element -> material map for the per-material internal length.
!     Filled by damnlmatmap, which the caller supplies because the
!     averaging routines do not receive ielmat.  matbuilt=0 means nobody
!     has supplied it, and then every element falls back to the single
!     global length - which is exactly the behaviour before W6, so a tree
!     in which no caller has been updated yet still answers as it did.
!
      integer, allocatable :: elmat(:)
      integer :: matbuilt=0
!
!     ellover records that the environment overrode a card, so the run
!     can say so once instead of leaving the answer depending on
!     something the input deck does not mention.
!
      integer :: ellover=0,ellsaid=0
!
!     T1 CHAIN FACTOR and PROBE MODE.
!
!     chainf(i) is d(ebar_i)/d(e_i): how much the REGULARISED driving
!     variable of element i moves when its own LOCAL one moves.  It is
!     what turns a local derivative into the diagonal of the nonlocal
!     one, and it is defined per backend - see damnlchain.  chnbuilt is
!     0 until a backend has filled it, so a backend that has no cheap
!     expression for it (GRADIENT) simply reports "not available" and
!     nothing downstream changes.
!
!     iprobee(i) is set while resultsmech perturbs the strain of element
!     i to measure dD/d(eps).  With it set, damnonlocalval reports "not
!     available" for THAT element, which makes damageupdatepoint fall
!     back to its LOCAL driving variable through the contract that
!     accessor already has.  That is the whole reason calcdamage.f needs
!     no change.
!
!     IT IS PER ELEMENT, NOT A SINGLE FLAG, AND THAT IS NOT DEFENSIVE.
!     results.c runs resultsmech on pthreads (pthread_create of
!     resultsmechmt), splitting the ELEMENT range across threads.  One
!     shared flag would let thread A, probing its own element, silently
!     switch thread B's element to the local driving variable and change
!     B's answer.  Indexing by element makes every write land in a slot
!     only its own thread touches, so the threads cannot see each other.
!     The gate runs at OMP_NUM_THREADS=1 and would never have shown it.
!
      real*8, allocatable :: chainf(:)
      integer, allocatable :: iprobee(:)
      integer :: chnbuilt=0,ichainon=0,iinner=0,ichainset=0
      end module damnlmod
!
!     ------------------------------------------------------------------
!
      subroutine damnlellinit()
      use damnlmod
      implicit none
      character*132 carg
      integer i,j,k
!
!     Read once.  Both controls are OFF unless set, and when both are off
!     ell2e is filled with the single global value, which reproduces the
!     previous binary exactly.
!
      if(ellinit.ne.0) return
      ellinit=1
!
!     ell multiplier per MATERIAL INDEX, comma separated, in deck order.
!     "1.0,0.2" gives material 1 the full ell and material 2 one fifth.
!     The material card is where this belongs eventually; the coefficient
!     it produces is identical either way.
!
      call getenv('CCX_DAMAGE_NONLOCAL_ELL_MAT',carg)
      if(carg(1:1).ne.' ') then
        j=1
        k=0
        do i=1,len_trim(carg)+1
          if((i.gt.len_trim(carg)).or.(carg(i:i).eq.',')) then
            if(i.gt.j) then
              k=k+1
              if(k.le.16) then
                read(carg(j:i-1),*,err=8300,end=8300) ellmf(k)
                if(ellmf(k).lt.0.d0) ellmf(k)=0.d0
              endif
            endif
            j=i+1
          endif
        enddo
 8300   continue
        ellmn=k
      endif
!
!     LOCALIZING gradient damage: "n,R".  g(D)=(1-R)*(1-D)**n+R, so the
!     interaction radius collapses inside the band and is unchanged in
!     the undamaged material.  R keeps the operator non-degenerate.
!
!     CCX_DAMAGE_NONLOCAL_CHAIN=1 arms the chain rule.  OFF BY DEFAULT,
!     and the default was changed from on to off in review: the factor is
!     the diagonal of the monolithic problem, while the scheme that runs
!     here freezes the driving variable for the whole Newton loop, whose
!     exact Jacobian has a zero on this path.  Until the refresh moves
!     inside the loop, arming it improves the Jacobian of a problem this
!     code does not solve.  The switch stays because that change is the
!     next step and this is the half of it that is already measured.
!
!     The environment value arrives through damnonlocalset, called from
!     nonlingeo.c only when CCX_DAMAGE_NONLOCAL is set, so a nonzero
!     ellsave at this point means the environment spoke.  Combined with a
!     card, that is an override and is announced by damnlellsay.
!
      if(ellsave.gt.0.d0) ellover=1
!
      call getenv('CCX_DAMAGE_NONLOCAL_CHAIN',carg)
      if(carg(1:1).ne.' ') then
        ichainset=1
        if(carg(1:1).eq.'1') ichainon=1
      endif
!
!     CCX_DAMAGE_NONLOCAL_INNER=1 refreshes the regularised driving
!     variable ONCE PER NEWTON ITERATION instead of once per increment.
!     That is what gives the chain factor an addressee: while the field is
!     frozen for the whole loop the exact Jacobian has a zero on that path
!     (see damnlchain), so a diagonal cannot help a scheme that does not
!     contain it.
!
!     IT ARMS THE CHAIN FACTOR TOO, unless CCX_DAMAGE_NONLOCAL_CHAIN says
!     otherwise explicitly.  The two are halves of one change and each
!     half is expected to look WORSE alone: refreshing inside Newton makes
!     the coupling Gauss-Seidel, which costs iterations unless the tangent
!     accounts for it, and the tangent term costs a little unless the
!     field it differentiates actually moves.  Measured separately, the
!     chain half alone cost +0.17 percent iterations.  Setting CHAIN
!     explicitly still wins, so the halves can still be measured apart.
!
      call getenv('CCX_DAMAGE_NONLOCAL_INNER',carg)
      if(carg(1:1).eq.'1') then
        iinner=1
        if(ichainset.eq.0) ichainon=1
      endif
!
      call getenv('CCX_DAMAGE_NONLOCAL_DTSCALE',carg)
      if(carg(1:1).eq.'0') idtscale=0
!
      call getenv('CCX_DAMAGE_NONLOCAL_LOCALIZING',carg)
      if(carg(1:1).ne.' ') then
        nlocn=2.d0
        nlocr=1.d-2
        j=index(carg,',')
        if(j.gt.1) then
          read(carg(1:j-1),*,err=8310,end=8310) nlocn
          read(carg(j+1:),*,err=8310,end=8310) nlocr
        else
          read(carg,*,err=8310,end=8310) nlocn
        endif
 8310   continue
        if(nlocn.lt.0.d0) nlocn=2.d0
        if(nlocr.lt.1.d-6) nlocr=1.d-6
        if(nlocr.gt.1.d0) nlocr=1.d0
        nlocon=1
      endif
      return
      end
!
!     ------------------------------------------------------------------
!
!     ==================================================================
!     Element geometry and integration-point count for the nonlocal
!     average, for every volume family the damage model can reach.
!
!     WHY THIS EXISTS.  Both the neighbour list and the local driving
!     variable were written for C3D4 and skipped everything else, which
!     was correct while DE1 itself refused non-tetrahedra.  CB1 lifted
!     that refusal, so a hexahedral deck now reaches the damage law with
!     no regularisation available to it - the two halves have to move
!     together or the internal length silently stops applying to exactly
!     the meshes CB1 just enabled.
!
!     CENTROID is the mean of the CORNER nodes.  Midside nodes are left
!     out on purpose: including them biases the centre towards the curved
!     faces of a distorted quadratic element, and the centre is only used
!     to place the element in the neighbour search, where a bias is a
!     systematic error in who is a neighbour of whom.
!
!     VOLUME is the corner hull, decomposed into tetrahedra: one for a
!     tet, three for a wedge, six for a hexahedron in the standard
!     diagonal split. For a quadratic element with curved edges this is
!     the straight-edged volume, which is what the weight needs - the
!     weight asks how much material sits near a point, and the corner
!     hull answers that to the accuracy the Gaussian kernel can use.
!
!     IOK=0 means the family is not one this routine knows, and the
!     caller must then leave the element out of the average entirely
!     rather than average it with a wrong weight.
!     ==================================================================
!
      subroutine damnlelgeom(lakonl,kon,co,indexe,cx,cy,cz,vol,nip,iok)
      implicit none
      character*8 lakonl
      integer kon(*),indexe,nip,iok,ncor,j,node,it,nt
      integer itet(4,6)
      real*8 co(3,*),cx,cy,cz,vol,xl(3,8),d6,a(3),b(3),c(3)
!
      iok=0
      cx=0.d0
      cy=0.d0
      cz=0.d0
      vol=0.d0
      nip=0
      if(lakonl(1:1).ne.'C') return
!
!     corner count and integration-point count per family.  The point
!     counts mirror calcdamage.f, which is the routine whose state this
!     average consumes: if the two ever disagree, the average would read
!     slots calcdamage never filled.
!
      if(lakonl(4:5).eq.'8R') then
        ncor=8
        nip=1
      elseif(lakonl(4:4).eq.'8') then
        ncor=8
        nip=8
      elseif(lakonl(4:5).eq.'20') then
        ncor=8
        if(lakonl(6:6).eq.'R') then
          nip=8
        else
          nip=27
        endif
      elseif(lakonl(4:5).eq.'10') then
        ncor=4
        nip=4
      elseif(lakonl(4:4).eq.'4') then
        ncor=4
        nip=1
      elseif(lakonl(4:5).eq.'15') then
        ncor=6
        nip=9
      elseif(lakonl(4:4).eq.'6') then
        ncor=6
        nip=2
      else
        return
      endif
!
      do j=1,ncor
        node=kon(indexe+j)
        if(node.le.0) return
        xl(1,j)=co(1,node)
        xl(2,j)=co(2,node)
        xl(3,j)=co(3,node)
        cx=cx+xl(1,j)
        cy=cy+xl(2,j)
        cz=cz+xl(3,j)
      enddo
      cx=cx/dble(ncor)
      cy=cy/dble(ncor)
      cz=cz/dble(ncor)
!
      if(ncor.eq.4) then
        nt=1
        itet(1,1)=1
        itet(2,1)=2
        itet(3,1)=3
        itet(4,1)=4
      elseif(ncor.eq.6) then
!
!       wedge 1-2-3 / 4-5-6 into three tetrahedra
!
        nt=3
        itet(1,1)=1
        itet(2,1)=2
        itet(3,1)=3
        itet(4,1)=5
        itet(1,2)=1
        itet(2,2)=3
        itet(3,2)=6
        itet(4,2)=5
        itet(1,3)=1
        itet(2,3)=6
        itet(3,3)=4
        itet(4,3)=5
      else
!
!       hexahedron into six tetrahedra sharing the 1-7 diagonal
!
        nt=6
        itet(1,1)=1
        itet(2,1)=2
        itet(3,1)=3
        itet(4,1)=7
        itet(1,2)=1
        itet(2,2)=3
        itet(3,2)=4
        itet(4,2)=7
        itet(1,3)=1
        itet(2,3)=4
        itet(3,3)=8
        itet(4,3)=7
        itet(1,4)=1
        itet(2,4)=8
        itet(3,4)=5
        itet(4,4)=7
        itet(1,5)=1
        itet(2,5)=5
        itet(3,5)=6
        itet(4,5)=7
        itet(1,6)=1
        itet(2,6)=6
        itet(3,6)=2
        itet(4,6)=7
      endif
!
      do it=1,nt
        do j=1,3
          a(j)=xl(j,itet(2,it))-xl(j,itet(1,it))
          b(j)=xl(j,itet(3,it))-xl(j,itet(1,it))
          c(j)=xl(j,itet(4,it))-xl(j,itet(1,it))
        enddo
        d6=a(1)*(b(2)*c(3)-b(3)*c(2))
     &    -a(2)*(b(1)*c(3)-b(3)*c(1))
     &    +a(3)*(b(1)*c(2)-b(2)*c(1))
        vol=vol+dabs(d6)/6.d0
      enddo
      if(vol.le.0.d0) return
      iok=1
      return
      end
!
!     ==================================================================
!     Element matrices of the implicit-gradient equation, for ANY volume
!     family - this is what lifts C3D4 out of the GRADIENT backend.
!
!     Peerlings et al. 1996 write the nonlocal equation as
!
!         ebar - div( ell^2 grad ebar ) = e
!
!     whose weak form on one element contributes
!
!         ml(a)    = int N_a dV            (row-sum lumped mass)
!         kel(a,b) = int grad N_a . grad N_b dV
!
!     and NOTHING else.  The backend used to hardwire both for a linear
!     tetrahedron, where ml = V/4 and grad N is constant, and skip every
!     other family.  Both integrals exist for every family; only the
!     quadrature and the shape functions change, and ccx already owns
!     those.
!
!     Two decisions, and the reasons they are decisions and not defaults:
!
!     1. The shape functions come from ccx's OWN shape4tet/shape8h/
!        shape6w and the quadrature from gauss.f, the same pair e_c3d.f
!        uses.  A private copy of the element library would be a second
!        source of truth for geometry that ccx already defines, and would
!        drift from it at the first change to either.
!
!     2. ebar is interpolated LINEARLY, on the CORNER nodes, for every
!        family including the quadratic ones.  That is the standard
!        choice for the implicit-gradient model (Peerlings et al. 1996):
!        the nonlocal field is taken one order below the displacement
!        field, which is both cheaper and better behaved - an equal-order
!        ebar on a quadratic element admits oscillations the averaging is
!        supposed to remove.  The consequence to state plainly: on a
!        quadratic element with CURVED edges the geometry seen here is
!        the straight-edged one.  For the internal length that is a
!        second-order error in a quantity that is itself a modelling
!        constant; for a straight-edged quadratic element it is exact.
!
!     The lumping is row-sum, which for these shape functions is
!     int N_a dV: positive for all three linear families, and exactly the
!     V/4 the old tetrahedral path used.  That equality is not a
!     coincidence to be trusted - it is checked in
!     test/nonlocal/elmk_test.f against the closed form.
!
      subroutine damnlelmk(lakonl,kon,co,indexe,nc,nod,ml,kel,iok)
      implicit none
      character*8 lakonl
      integer kon(*),indexe,nc,nod(8),iok,i,j,k,nip,ifam,node
      real*8 co(3,*),ml(8),kel(8,8),xl(3,8),shp(4,8),xsj,xi,et,ze,w
!
      include "gauss.f"
!
      iok=0
      nc=0
      nip=0
      ifam=0
      do i=1,8
        nod(i)=0
        ml(i)=0.d0
        do j=1,8
          kel(i,j)=0.d0
        enddo
      enddo
      if(lakonl(1:1).ne.'C') return
!
!     family by CORNER count.  The corner count is what ebar is
!     interpolated on, so C3D20 and C3D8 land in the same branch by
!     construction rather than by a list that has to be maintained.
!     ifam: 1 = tetrahedron, 2 = wedge, 3 = hexahedron.
!
      if((lakonl(4:5).eq.'20').or.(lakonl(4:4).eq.'8')) then
        nc=8
        ifam=3
        nip=8
      elseif((lakonl(4:5).eq.'15').or.(lakonl(4:4).eq.'6')) then
        nc=6
        ifam=2
        nip=6
      elseif((lakonl(4:5).eq.'10').or.(lakonl(4:4).eq.'4')) then
        nc=4
        ifam=1
        nip=4
      else
        nc=0
        return
      endif
!
      do j=1,nc
        node=kon(indexe+j)
        if(node.le.0) then
          nc=0
          return
        endif
        nod(j)=node
        xl(1,j)=co(1,node)
        xl(2,j)=co(2,node)
        xl(3,j)=co(3,node)
      enddo
!
      do k=1,nip
        if(ifam.eq.3) then
          xi=gauss3d2(1,k)
          et=gauss3d2(2,k)
          ze=gauss3d2(3,k)
          w=weight3d2(k)
          call shape8h(xi,et,ze,xl,xsj,shp,3)
        elseif(ifam.eq.2) then
          xi=gauss3d10(1,k)
          et=gauss3d10(2,k)
          ze=gauss3d10(3,k)
          w=weight3d10(k)
          call shape6w(xi,et,ze,xl,xsj,shp,3)
        else
          xi=gauss3d5(1,k)
          et=gauss3d5(2,k)
          ze=gauss3d5(3,k)
          w=weight3d5(k)
          call shape4tet(xi,et,ze,xl,xsj,shp,3)
        endif
!
!       a non-positive Jacobian is a degenerate or inverted element;
!       refuse the whole element rather than integrate part of it, so
!       the caller can keep it out of the assembly instead of getting a
!       matrix that is quietly wrong.
!
        if(xsj.le.0.d0) then
          nc=0
          return
        endif
        w=w*xsj
        do i=1,nc
          ml(i)=ml(i)+w*shp(4,i)
          do j=1,nc
            kel(i,j)=kel(i,j)+w*(shp(1,i)*shp(1,j)+shp(2,i)*shp(2,j)
     &           +shp(3,i)*shp(3,j))
          enddo
        enddo
      enddo
!
!     a lumped mass must be positive for every node, or the mass term
!     stops being a norm and the CG below loses its guarantee
!
      do i=1,nc
        if(ml(i).le.0.d0) then
          nc=0
          return
        endif
      enddo
      iok=1
      return
      end
!
!
!     ==================================================================
!     IS THE INTERNAL LENGTH RESOLVED BY THE MESH AT ALL?
!
!     One warning for both backends, because it is one disease with two
!     symptoms.  When the requested ell falls below what the mesh can
!     resolve, the model is ARMED BUT DEGENERATE: every element averages
!     essentially itself, so the user asked for regularisation and did
!     not get it - and that has been accepted in silence.
!
!     Measured, both on the reference deck:
!
!       INTEGRAL, ell=0.01, mean neighbours 0.98
!           the run DIVERGES - increments 4U, 5U, 6U in a row - while at
!           ell=0.2, mean neighbours 27.6, the same configuration
!           completes.  So it is not that a directional width is heavier;
!           it holds fine when the regularisation is real.
!       GRADIENT, ell=0.01
!           the run COMPLETED and the answer was WRONG: the field handed
!           to the damage law differed from the local one by up to 83 per
!           cent, entirely from the element->node->element projection.
!
!     A wrong answer that converges and a divergence with the right
!     answer look nothing alike, which is exactly why neither was
!     recognised as the same thing for as long as it was.
!
!     THE NUMBER IS PARTICIPATION, NOT ell/h.  How many elements actually
!     enter the average is what the question asks, it is already measured
!     for the integral backend, and it needs no threshold chosen by
!     taste: an element is always in its own neighbourhood, so a count
!     below TWO means not one neighbour joined it.  That is not a tuned
!     bound, it is the definition of averaging nothing.  The other agent
!     proposed this measure; the floor of two is where it becomes
!     arithmetic rather than judgement.
!
      subroutine damnlresolve(part,ell,iwhat)
      implicit none
      real*8 part,ell
      integer iwhat
!
      if(ell.le.0.d0) return
      if(part.ge.2.d0) return
      write(*,*) '*WARNING in damnonlocal: the internal length is NOT'
      write(*,*) '         resolved by this mesh.  ell=',ell
      if(iwhat.eq.0) then
        write(*,*) '         mean elements per neighbourhood=',part
      else
        write(*,*) '         est. elements within 2*ell=',part
      endif
      write(*,*) '         Below 2 the average has no neighbour to'
      write(*,*) '         average with, so the model is ARMED but'
      write(*,*) '         DEGENERATE: the regularisation asked for is'
      write(*,*) '         not being applied.  Measured consequences on'
      write(*,*) '         the reference deck: the integral backend'
      write(*,*) '         DIVERGES here, the gradient one COMPLETES'
      write(*,*) '         with a driving variable off by up to 83 per'
      write(*,*) '         cent.  Refine the mesh or raise ell.'
      return
      end
!
!
!     ==================================================================
!     THE LENGTH AN ELEMENT WAS ACTUALLY AVERAGED WITH.
!
!     Handed to the crack-band width so the damage law can use 2*ell in
!     place of the element size when the regularisation is real.  Agent 2
!     measured why this is needed: the nonlocal model widens the band to
!     about 2*ell and holds it against mesh refinement, but the law still
!     divides by a length of order h, so the dissipation follows the WIDTH
!     instead of staying put - W_post 8.87 -> 11.26 -> 16.55 -> 29.99 as
!     ell goes 0 -> 0.25 -> 0.5 -> 1.0 at fixed h, while W_pre agrees to
!     0.03 per cent.  The band obeys ell; the law does not know.
!
!     NOT damnlellmax, which is the largest length configured anywhere.
!     This is the element's OWN length, including the per-material
!     multiplier and the localizing g(D) - that is, the number standing in
!     the operator that produced this element's averaged value.  Handing
!     back anything else would substitute a length the element was not
!     averaged with.
!
!     iok=1 only when the mesh RESOLVES that length, by the same test the
!     shared guard uses: an element is always in its own neighbourhood, so
!     participation below two means no neighbour joined and the band is
!     not 2*ell but h.  Agent 2's scaling sweep found the same boundary
!     independently - w/2ell is 1.50 / 1.08 / 0.88 at ell = 0.25 / 0.5 /
!     1.0, and the outlier is exactly the row where 2*ell/h = 2.  The
!     caller therefore does not have to restate the criterion, and cannot
!     restate it differently by accident.
!
      subroutine damnlelleff(iel,ell,iok)
      use damnlmod
      implicit none
      integer iel,iok
      real*8 ell,part,vmean
!
      ell=0.d0
      iok=0
      if((iel.lt.1).or.(iel.gt.nesave)) return
!
      if(imodenl.eq.1) then
        if(.not.allocated(ell2e)) return
        if(.not.allocated(gnc)) return
        if(gnc(iel).le.0) return
        if(ell2e(iel).le.0.d0) return
        ell=dsqrt(ell2e(iel))
!
!       participation, estimated as in damgradient: the 2*ell ball over
!       this element's own volume, floored at one for the element itself
!
        if(.not.allocated(vele)) return
        if(vele(iel).le.0.d0) return
        part=4.18879020478639d0*(2.d0*ell)**3/vele(iel)
        if(part.lt.1.d0) part=1.d0
      else
        if(.not.allocated(nbcount)) return
        if(.not.allocated(evol)) return
        if(evol(iel).le.0.d0) return
        call damnlellel(iel,ell)
        if(ell.le.0.d0) return
!
!       counted, not estimated: the integral backend keeps the list
!
        part=dble(nbcount(iel))
      endif
      if(part.lt.2.d0) return
      iok=1
      return
      end
!
!
!     ==================================================================
!     Step sizes for the time-consistent lag.  damnldtbase is called by
!     calcdamagebase just before a backend rebuilds dpsave from the
!     converged increment, so both sizes are that increment's.
!     damnldtstep is called by results.c, before its threads start, with
!     the size of the attempt being solved - what the lagged rate is
!     multiplied by.
!
      subroutine damnldtbase(dt)
      use damnlmod
      implicit none
      real*8 dt
      dtref=dt
      dtcur=dt
      return
      end
!
      subroutine damnldtstep(dt)
      use damnlmod
      implicit none
      real*8 dt
      dtcur=dt
      return
      end
!
      subroutine damnonlocalset(ellin)
      use damnlmod
      implicit none
      real*8 ellin
      ellsave=ellin
      return
      end
!
!     W6 setter, called once per material from the *DAMAGE INITIATION
!     reader.  ell<=0 means "no NONLOCAL= on this card" and is recorded as
!     such rather than rejected: a deck may regularise one material and
!     not another.
!
!     ==================================================================
!     Refresh the regularised driving variable INSIDE the Newton loop.
!
!     WHERE THIS MUST BE CALLED FROM, and why it cannot be anywhere else.
!     The average needs xstate for every element of a neighbourhood at the
!     SAME iteration.  results.c runs resultsmech on num_cpus threads split
!     by element range and joins them afterwards, so before the join part
!     of xstate still belongs to the previous iteration: an element would
!     be averaged over a mixture of two iterations and the answer would
!     depend on the thread count.  No lock fixes that.  So the only sound
!     call site is AFTER the join - in practice, in the Newton loop of
!     nonlingeo.c, immediately after results() returns.
!
!     WHAT IT DOES AND DOES NOT DO.  It does not make the problem
!     monolithic.  It turns "ebar frozen for the increment" into "ebar
!     updated once per iteration", i.e. a Gauss-Seidel coupling inside
!     Newton.  Peerlings et al. 1996 is the monolithic alternative and
!     needs a nodal degree of freedom; this does not, and does not claim
!     that convergence rate.
!
!     Inert unless CCX_DAMAGE_NONLOCAL_INNER=1, so a binary built with it
!     is byte for byte the previous one until the switch is set.
!     ==================================================================
!
      subroutine damnlrefresh(ipkon,kon,lakon,co,ne0,mi,xstate,
     &     xstateini,nstate_,ielmat,dam)
      use damnlmod
      implicit none
      character*8 lakon(*)
      integer ipkon(*),kon(*),ne0,mi(*),nstate_,ielmat(mi(3),*)
      integer ionnl
      real*8 co(3,*),xstate(nstate_,mi(1),*),
     &     xstateini(nstate_,mi(1),*),dam(mi(1),*)
!
      call damnlactive(ionnl)
      if(ionnl.eq.0) return
      call damnlellinit()
      if(iinner.eq.0) return
!
!     Supply the element -> material map while ielmat is in scope, so the
!     per-material internal length applies on this path.  Idempotent.
!
      call damnlmatmap(ipkon,lakon,ielmat,ne0,mi)
!
!     The same dispatch calcdamagebase makes at the start of an
!     increment.  Kept identical on purpose: if the two ever disagree,
!     the field seen inside Newton would not be the field the increment
!     commits.
!
      if(imodenl.eq.2) then
        call damfrozen(ipkon,kon,lakon,co,ne0,mi,xstate,
     &       xstateini,nstate_)
      elseif(imodenl.eq.1) then
        call damgradient(ipkon,kon,lakon,co,ne0,mi,xstate,
     &       xstateini,nstate_,ielmat,dam)
      else
        call damnonlocal(ipkon,kon,lakon,co,ne0,mi,xstate,
     &       xstateini,nstate_)
      endif
!
!     the field just rebuilt was measured over the step being solved
!
      dtref=dtcur
      return
      end
!
      subroutine damnlinnerget(ion)
      use damnlmod
      implicit none
      integer ion
      call damnlellinit()
      ion=iinner
      return
      end
!
!     Supply the element -> material map.  Called by whoever has ielmat
!     in scope; it is idempotent and cheap, so calling it every increment
!     costs nothing and calling it never costs only the per-material
!     length, which then falls back to the global one.
!
      subroutine damnlmatmap(ipkon,lakon,ielmat,ne0,mi)
      use damnlmod
      implicit none
      character*8 lakon(*)
      integer ne0,mi(*),ipkon(*),ielmat(mi(3),*),i
      if(ncard.le.0) return
      if(allocated(elmat)) then
        if(size(elmat).lt.ne0) deallocate(elmat)
      endif
      if(.not.allocated(elmat)) allocate(elmat(ne0))
      do i=1,ne0
        elmat(i)=0
        if(lakon(i)(1:1).ne.'C') cycle
        elmat(i)=ielmat(1,i)
      enddo
      matbuilt=1
      return
      end
!
!     The internal length that applies to ONE element.
!
!     PRECEDENCE, agreed in review: the card sets the value and the
!     environment overrides it.  The override is announced once, naming
!     the material, what the card asked for and what is actually used -
!     otherwise the answer depends on something the deck does not
!     mention, which is the disease W6 exists to cure.
!
!     The largest internal length any material asks for, so the neighbour
!     search can size its radius to cover all of them.  Falls back to the
!     global value when no card carried NONLOCAL=.
!
!     Say once, in the log, that the environment overrode what the deck
!     asked for.  Precedence agreed in review: the card sets the value,
!     CCX_DAMAGE_NONLOCAL overrides it.  An override that is not printed
!     would put the answer back under the control of something the input
!     file does not mention, which is the disease W6 exists to cure, so
!     the announcement is not decoration.
!
      subroutine damnlellsay()
      use damnlmod
      implicit none
      integer i
      if(ellsaid.ne.0) return
      ellsaid=1
      if(ncard.le.0) return
      if(ellover.eq.0) then
        write(*,*) '[DAMAGE NONLOCAL] internal length from the',
     &       ' material card, per material:'
        do i=1,64
          if(ellcard(i).gt.0.d0) then
            write(*,'(a,i5,a,e13.6)')
     &        ' [DAMAGE NONLOCAL]   material ',i,'  ell=',ellcard(i)
          endif
        enddo
      else
        write(*,*) '[DAMAGE NONLOCAL] CCX_DAMAGE_NONLOCAL OVERRIDES',
     &       ' the NONLOCAL= parameter on the material card.'
        do i=1,64
          if(ellcard(i).gt.0.d0) then
            write(*,'(a,i5,a,e13.6,a,e13.6)')
     &        ' [DAMAGE NONLOCAL]   material ',i,'  card asked ',
     &        ellcard(i),'  applied ',ellsave
          endif
        enddo
      endif
      call flush(6)
      return
      end
!
!     Is ANY internal length configured - from the environment or from a
!     material card?
!
!     THIS EXISTS BECAUSE THE TWO SOURCES HAD DIFFERENT REACH.  ellsave
!     is written only by damnonlocalset, which nonlingeo.c calls only
!     when CCX_DAMAGE_NONLOCAL is set, so every guard phrased as
!     "ellsave > 0" meant "the environment spoke" and not "a length
!     exists".  A deck that put NONLOCAL= on its material card therefore
!     never entered the nonlocal block at all: the regularisation was
!     silently absent, and the run looked like a healthy local one.
!     Found by Agent 2 in review of W6.  Every guard now asks this
!     instead, so the card reaches as far as the environment does.
!
      subroutine damnlactive(ion)
      use damnlmod
      implicit none
      integer ion
      ion=0
      call damnlellinit()
      if(ellsave.gt.0.d0) ion=1
      if(ncard.gt.0) ion=1
      return
      end
!
      subroutine damnlellmax(ellmax)
      use damnlmod
      implicit none
      real*8 ellmax
      integer i
!
!     Initialise before reading.  ellover and the card table are set
!     by damnlellinit, and this routine used to read them without
!     calling it - so a caller that reached here first would see
!     ellover=0 and take a card value where the environment was
!     meant to override.  Not reachable today with a wrong ANSWER,
!     because the one early caller uses this only as a yes/no gate,
!     but it is the exact shape of the cbmode defect: a value read
!     by a path that did not set it.  Found by searching my own
!     files for that shape rather than by a run.
!
      call damnlellinit()
      ellmax=ellsave
      if(ncard.le.0) return
      if(ellover.eq.1) return
      do i=1,64
        if(ellcard(i).gt.ellmax) ellmax=ellcard(i)
      enddo
      return
      end
!
      subroutine damnlellel(iel,ell)
      use damnlmod
      implicit none
      integer iel,im
      real*8 ell,emx
!
!     Initialise before reading.  ellover and the card table are set
!     by damnlellinit, and this routine used to read them without
!     calling it - so a caller that reached here first would see
!     ellover=0 and take a card value where the environment was
!     meant to override.  Not reachable today with a wrong ANSWER,
!     because the one early caller uses this only as a yes/no gate,
!     but it is the exact shape of the cbmode defect: a value read
!     by a path that did not set it.  Found by searching my own
!     files for that shape rather than by a run.
!
      call damnlellinit()
      ell=ellsave
      if(ncard.le.0) return
!
!     A card-only run has ellsave=0, and handing that back as the default
!     would make the kernel width zero for any material that named no
!     length of its own.  The largest configured length is the honest
!     default there: generous rather than degenerate.
!
      if(ell.le.0.d0) then
        call damnlellmax(emx)
        ell=emx
      endif
      if(matbuilt.eq.0) return
      if(.not.allocated(elmat)) return
      if((iel.lt.1).or.(iel.gt.size(elmat))) return
      im=elmat(iel)
      if((im.lt.1).or.(im.gt.64)) return
      if(ellcard(im).le.0.d0) return
      if(ellover.eq.1) return
      ell=ellcard(im)
      return
      end
!
      subroutine damnlellsetmat(imat,ell)
      use damnlmod
      implicit none
      integer imat
      real*8 ell
      if((imat.lt.1).or.(imat.gt.64)) then
        write(*,*) '*ERROR in damnlellsetmat: material index out of'
        write(*,*) '       range for the nonlocal length table:',imat
        call exit(201)
      endif
      if(ell.le.0.d0) return
      ellcard(imat)=ell
      ncard=ncard+1
      return
      end
!
!     Hand back what the cards said.  iok=0 means no card carried
!     NONLOCAL=, which is how a caller tells "not regularised" from
!     "regularised with a length of zero", a distinction the single
!     global ellsave cannot make.
!
      subroutine damnlellgetmat(imat,ell,iok)
      use damnlmod
      implicit none
      integer imat,iok
      real*8 ell
!
!     Initialise before reading.  ellover and the card table are set
!     by damnlellinit, and this routine used to read them without
!     calling it - so a caller that reached here first would see
!     ellover=0 and take a card value where the environment was
!     meant to override.  Not reachable today with a wrong ANSWER,
!     because the one early caller uses this only as a yes/no gate,
!     but it is the exact shape of the cbmode defect: a value read
!     by a path that did not set it.  Found by searching my own
!     files for that shape rather than by a run.
!
      call damnlellinit()
      iok=0
      ell=0.d0
      if(ncard.le.0) return
      if((imat.lt.1).or.(imat.gt.64)) return
      if(ellcard(imat).le.0.d0) return
      ell=ellcard(imat)
      iok=1
      return
      end
!
      subroutine damnlellcardcount(n)
      use damnlmod
      implicit none
      integer n
      n=ncard
      return
      end
!
      subroutine damnonlocalget(ellout)
      use damnlmod
      implicit none
      real*8 ellout
      ellout=ellsave
      return
      end
!
!     Hand back the averaged driving variable for one element.  iok=0
!     means "not available" and the caller must keep its local value,
!     which is what happens before the first refresh.
!
      subroutine damnonlocalval(iel,val,iok)
      use damnlmod
      implicit none
      integer iel,iok
      integer ionnl
      real*8 val
      iok=0
      val=0.d0
      call damnlactive(ionnl)
      if(ionnl.eq.0) return
!
!     PROBE MODE.  While resultsmech is measuring dD/d(eps) the
!     regularised value must not be handed out: it is a stored field and
!     no strain perturbation can move it, so every difference quotient
!     taken against it is exactly zero.  Reporting "not available" here
!     makes the caller use its own local value, which is the quantity
!     whose derivative is actually wanted; damnlchain then supplies the
!     factor that turns that local derivative into the nonlocal one.
!
      if(.not.allocated(dpsave)) return
      if((iel.lt.1).or.(iel.gt.nesave)) return
      if(allocated(iprobee)) then
        if(iprobee(iel).ne.0) return
      endif
!
!     AN ELEMENT THE BACKEND LEFT OUT HAS NO REGULARISED VALUE, AND ZERO
!     IS NOT ONE.  dpsave is zero for an element the backend refused - a
!     family it cannot integrate, a non-positive Jacobian, or, in the
!     integral backend, one with no volume.  Handing that zero back with
!     iok=1 tells the caller "the regularised increment is zero", which is
!     a statement about the material, not about the backend: the element
!     would stop accumulating damage entirely.
!
!     Reporting iok=0 instead makes the caller fall back to its own local
!     value, which is what an unregularised element should use and what
!     the warning printed at build time already says happens.
!
!     HONEST SCOPE: this changes no answer on any deck in the gate.  The
!     elements either backend refuses on these decks are UC6 cohesive
!     facets, which never ask for this value.  It is a defect by
!     construction, found by reading the consumer after fixing the
!     producer, not by a run - and the gate confirms only that it costs
!     nothing, not that it was costing something.
!
      if(imodenl.eq.1) then
        if(allocated(gnc)) then
          if(gnc(iel).le.0) return
        endif
      else
        if(allocated(evol)) then
          if(evol(iel).le.0.d0) return
        endif
      endif
      val=dpsave(iel)
!
!     The lagged RATE times the current step, not the lagged increment.
!     Unscaled, a cutback left the damage law a full previous step's
!     increment at a step up to a thousand times smaller.  Measured on the
!     notched bar of mkcross.py (24 slices, CHARLEN=1, VISCOSITY=1e-3):
!       INTEGRAL l=0.375  unscaled: wall at t=0.3675, rc=201
!                         scaled:   rupture, rc=0, 124 inc, 276 iter
!       GRADIENT l=0.5    wall at t=0.2575 either way; its first Newton
!                         correction hits the cap at every step size, a
!                         different mechanism, still open
!     and the two hexahedral gate decks (hex-nonlocal-card, INTEGRAL, and
!     hex-gradient-card, GRADIENT): theta=0.215, rc=201 unscaled, the
!     end of the step, rc=0, scaled;
!     and on test/fast 10x6x6 GRADIENT l=0.4: 1150 iter against 1152,
!     same 90 deletions.  Refreshing inside Newton (INNER=1) rescues both
!     bar arms, NOT the hexahedral decks (theta=0.215 as unscaled), and
!     costs +30 percent iterations on that deck; refreshing only on cut
!     attempts rescued nothing this does not, so it was dropped.
!     The factor is 1 while the step does not change; when the solver
!     grows the step it grows the lagged increment with it, which is the
!     same consistency argument.
!
      if((idtscale.eq.1).and.(dtref.gt.0.d0).and.(dtcur.gt.0.d0))
     &     val=val*dtcur/dtref
      iok=1
      return
      end
!
!     ------------------------------------------------------------------
!
      subroutine damnlprobeset(iel,ion)
      use damnlmod
      implicit none
      integer iel,ion
      if(.not.allocated(iprobee)) return
      if((iel.lt.1).or.(iel.gt.nesave)) return
      iprobee(iel)=ion
      return
      end
!
!     ------------------------------------------------------------------
!
!     d(ebar_e)/d(e_e) for one element.  iok=0 means "no cheap expression
!     for this backend", and the caller must then leave its tangent
!     exactly as it was.
!
!     INTEGRAL (imodenl=0).  ebar_e = sum_f w_ef V_f e_f / sum_f w_ef V_f
!     with w(d)=exp(-(d/ell)^2), so the self term has d=0 and w_ee V_e is
!     simply V_e, and the factor is V_e / sum_f w_ef V_f.
!
!     WHAT THAT NUMBER IS, EXACTLY.  It is the diagonal of the MONOLITHIC
!     nonlocal problem - the one in which ebar responds to the trial
!     state.  It is NOT the Jacobian of the scheme this code runs, which
!     refreshes dpsave once per increment and therefore has a zero here
!     too, for the same reason spelled out under FROZEN below.  So this
!     is an approximation applied to a scheme it does not belong to, and
!     it pays nothing until the refresh moves inside the Newton loop:
!     measured, +0.17 percent iterations and 2 coefficients of 47400
!     moved in the structural probe.  That is why it is OFF by default.
!     Do not read "exact" into it - an earlier version of this comment
!     said exactly that and was wrong.
!
!     FROZEN (imodenl=2).  REFUSED, and the reason is the one that decides
!     this whole routine.  It is tempting to say that with no spatial
!     averaging ebar_e IS e_e, so the factor is 1.  That is wrong.
!     damfrozen writes dpsave once per increment, from calcdamagebase,
!     which nonlingeo.c calls only after the Newton loop has converged -
!     so during the iterations dpsave is a CONSTANT and
!
!         d(ebar_e)/d(e_e^trial) = 0 ,  not 1 .
!
!     A frozen value does not answer to the trial state; that is what the
!     word means.  The exact Jacobian of the scheme as implemented has a
!     zero on this path, so the original zero was CORRECT and a factor of
!     1 would insert a term the iteration does not contain.  Reported by
!     Agent 2 in review of T1, and confirmed at the three calcdamagebase
!     call sites in nonlingeo.c - all three are outside the Newton loop.
!
!     GRADIENT (imodenl=1).  NOT AVAILABLE, deliberately.  Here ebar
!     solves (I - ell^2 div grad) ebar = e, so the sensitivity is a row
!     of the INVERSE of that operator, i.e. a solve with a unit source -
!     not a ratio of weights.  Substituting the integral form's factor
!     would put a number unrelated to the operator into the tangent, so
!     this backend keeps its present behaviour until the diagonal is
!     computed from the same CG that produces ebar.
!
      subroutine damnlchain(iel,val,iok)
      use damnlmod
      implicit none
      integer iel,iok
      integer ionnl
      real*8 val
      iok=0
      val=1.d0
      call damnlactive(ionnl)
      if(ionnl.eq.0) return
!
!     The factor and the probe are one mechanism: reporting a factor that
!     the caller cannot pair with an open probe would scale a derivative
!     that is still zero.  So refuse unless the marker exists.
!
      call damnlellinit()
      if(ichainon.eq.0) return
      if(.not.allocated(iprobee)) return
      if((iel.lt.1).or.(iel.gt.nesave)) return
      if(imodenl.eq.2) return
      if(imodenl.ne.0) return
      if(chnbuilt.eq.0) return
      if(.not.allocated(chainf)) return
      if((iel.lt.1).or.(iel.gt.nesave)) return
      val=chainf(iel)
      if(val.le.0.d0) then
        val=1.d0
        return
      endif
      if(val.gt.1.d0) val=1.d0
      iok=1
      return
      end
!
!     ------------------------------------------------------------------
!
      subroutine damnonlocal(ipkon,kon,lakon,co,ne0,mi,xstate,
     &     xstateini,nstate_)
!
!     Refresh the average.  The neighbour list is built once, on the
!     REFERENCE coordinates: the neighbourhood is a material property,
!     and rebuilding it as the mesh distorts would let ell drift.
!
      use damnlmod
      implicit none
!
      character*8 lakon(*)
!
      integer ipkon(*),kon(*),ne0,mi(*),nstate_,
     &     i,j,k,m,indexe,node,nn,ip,jp,kp,ib,jb,kb,
     &     ncell,ix,iy,iz,icell,jcell,ifree,nipel,iokel,ip1,ionnl,
     &     nact
      real*8 co(3,*),xstate(nstate_,mi(1),*),
     &     xstateini(nstate_,mi(1),*),ell,
     &     xc,yc,zc,d2,w,swv,sv,det6,rmax,volel,elli,
     &     xmin,xmax,ymin,ymax,zmin,zmax,csize
!
      call damnlactive(ionnl)
      if(ionnl.eq.0) return
      call damnlellinit()
      ell=ellsave
      call damnlellsay()
!
!     W6.  The search radius must cover the LARGEST length any material
!     asks for, otherwise a material with a longer ell would silently get
!     a truncated neighbourhood - the average would still be computed and
!     would still look reasonable, which is the worst way for a length to
!     stop applying.  The per-element length itself enters the weight
!     below; only the radius is global, and it is generous rather than
!     wrong.
!
      call damnlellmax(rmax)
      if(rmax.le.0.d0) rmax=ell
      rmax=2.d0*rmax
!
      if((nbuilt.ne.0).and.(nesave.ne.ne0)) nbuilt=0
!
      if(nbuilt.eq.0) then
        if(allocated(cen)) deallocate(cen,evol,wgt,dploc,dpsave,
     &       nbhead,nbnext,nblist,nbstart,nbcount)
        if(allocated(chainf)) deallocate(chainf,iprobee)
        allocate(cen(3,ne0),evol(ne0),dploc(ne0),dpsave(ne0))
        allocate(chainf(ne0),iprobee(ne0))
        chnbuilt=0
        do i=1,ne0
          chainf(i)=1.d0
          iprobee(i)=0
        enddo
        do i=1,ne0
          dpsave(i)=0.d0
        enddo
!
        xmin=1.d30
        xmax=-1.d30
        ymin=1.d30
        ymax=-1.d30
        zmin=1.d30
        zmax=-1.d30
        if(allocated(nipe)) deallocate(nipe)
        allocate(nipe(ne0))
        do i=1,ne0
          evol(i)=0.d0
          cen(1,i)=0.d0
          cen(2,i)=0.d0
          cen(3,i)=0.d0
          nipe(i)=0
          if(ipkon(i).lt.0) then
            indexe=-ipkon(i)-2
          else
            indexe=ipkon(i)
          endif
          if(indexe.lt.0) cycle
!
!         Every volume family, not only C3D4.  CB1 lifted the DE1
!         refusal on non-tetrahedra, and a hexahedral deck that reaches
!         the damage law with no regularisation available to it is worse
!         than one that is refused: the internal length would silently
!         stop applying to exactly the meshes that were just enabled.
!
          call damnlelgeom(lakon(i),kon,co,indexe,xc,yc,zc,volel,
     &         nipel,iokel)
          if(iokel.eq.0) cycle
          cen(1,i)=xc
          cen(2,i)=yc
          cen(3,i)=zc
          nipe(i)=nipel
          evol(i)=volel
          xmin=dmin1(xmin,cen(1,i))
          xmax=dmax1(xmax,cen(1,i))
          ymin=dmin1(ymin,cen(2,i))
          ymax=dmax1(ymax,cen(2,i))
          zmin=dmin1(zmin,cen(3,i))
          zmax=dmax1(zmax,cen(3,i))
        enddo
!
!       bucket grid of cell size rmax, so a neighbourhood search only has
!       to look at the 27 cells around the element
!
        csize=rmax
        ip=max(1,int((xmax-xmin)/csize)+1)
        jp=max(1,int((ymax-ymin)/csize)+1)
        kp=max(1,int((zmax-zmin)/csize)+1)
        ncell=ip*jp*kp
        allocate(nbhead(ncell),nbnext(ne0))
        do i=1,ncell
          nbhead(i)=0
        enddo
        do i=1,ne0
          nbnext(i)=0
          if(evol(i).le.0.d0) cycle
          ix=min(ip-1,max(0,int((cen(1,i)-xmin)/csize)))
          iy=min(jp-1,max(0,int((cen(2,i)-ymin)/csize)))
          iz=min(kp-1,max(0,int((cen(3,i)-zmin)/csize)))
          icell=(iz*jp+iy)*ip+ix+1
          nbnext(i)=nbhead(icell)
          nbhead(icell)=i
        enddo
!
        allocate(nbcount(ne0),nbstart(ne0+1))
        do m=1,2
          ifree=0
          do i=1,ne0
            nn=0
            call damnlellel(i,elli)
!
!           With a length that comes only from a card, ellsave is zero,
!           so the fallback has to be the largest configured length and
!           not ellsave - otherwise this divides by zero in the weight
!           below for any element whose own material named none.
!
            if(elli.le.0.d0) elli=ell
            if(elli.le.0.d0) elli=0.5d0*rmax
            if(evol(i).gt.0.d0) then
              ix=min(ip-1,max(0,int((cen(1,i)-xmin)/csize)))
              iy=min(jp-1,max(0,int((cen(2,i)-ymin)/csize)))
              iz=min(kp-1,max(0,int((cen(3,i)-zmin)/csize)))
              do ib=max(0,ix-1),min(ip-1,ix+1)
                do jb=max(0,iy-1),min(jp-1,iy+1)
                  do kb=max(0,iz-1),min(kp-1,iz+1)
                    jcell=(kb*jp+jb)*ip+ib+1
                    j=nbhead(jcell)
                    do
                      if(j.eq.0) exit
                      d2=(cen(1,i)-cen(1,j))**2+(cen(2,i)-cen(2,j))**2
     &                  +(cen(3,i)-cen(3,j))**2
                      if(d2.le.rmax*rmax) then
                        nn=nn+1
                        if(m.eq.2) then
!
!                         The kernel width is the length of the element
!                         being averaged - i, the receiver - not of the
!                         neighbour: ell is a property of the material
!                         whose damage is being regularised.
!
                          nblist(ifree+nn)=j
                          wgt(ifree+nn)=dexp(-d2/(elli*elli))*evol(j)
                        endif
                      endif
                      j=nbnext(j)
                    enddo
                  enddo
                enddo
              enddo
            endif
            nbcount(i)=nn
            nbstart(i)=ifree+1
            ifree=ifree+nn
          enddo
          nbstart(ne0+1)=ifree+1
          if(m.eq.1) allocate(nblist(max(1,ifree)),wgt(max(1,ifree)))
        enddo
        nbuilt=1
        nesave=ne0
!
!       PER REGULARISED ELEMENT, NOT PER ELEMENT IN THE MODEL.  The
!       divisor used to be ne0, which counts every element the deck has -
!       including the UC6 cohesive facets, which carry no volume and are
!       in nobody's neighbourhood, their own included.  On the reference
!       deck that is 2160 tetrahedra against 2196 elements, and the
!       average came out 2160/2196 = 0.9836.
!
!       An element is always in its OWN neighbourhood, so this average
!       cannot honestly be below 1.  It read 0.98 - which is the shape of
!       defect this whole session has been about: a number that looks
!       like a mild result and is actually an impossible one.  Found by
!       the other agent reading a banner in a log of mine.
!
!       hex-nonlocal-card pins this at 1.0 and stays green: that deck is
!       all hexahedra, so the two divisors agree there and the case was
!       never able to see the difference.
!
        nact=0
        do i=1,ne0
          if(evol(i).gt.0.d0) nact=nact+1
        enddo
        write(*,*) '[DAMAGE NONLOCAL] ell=',ell,' radius=',rmax,
     &       ' mean neighbours=',dble(ifree)/dble(max(1,nact)),
     &       ' over ',nact,' regularised of ',ne0
        call damnlresolve(dble(ifree)/dble(max(1,nact)),ell,0)
        call flush(6)
      endif
!
!     Local driving variable.  mi(1) is the ALLOCATED integration-point
!     stride - the maximum over every element type, 3 as soon as a UC6
!     facet exists - while a C3D4 has exactly ONE integration point
!     (calcdamage.f sets mint3d=1 for lakonl(4:4)=='4').  Looping to
!     mi(1) would average in uninitialised slots, so the loop below goes
!     to nipe(i), which calcdamage's own table fills.
!
!     THIS COMMENT USED TO END: "Only C3D4 reaches here, evol being zero
!     elsewhere, so the count is 1 by construction."  That was true until
!     3368f12 lifted the C3D4 restriction from this backend, and false
!     from that commit on - hexahedra reach here and their count is 8 or
!     27.  It sat directly above the loop that was added BECAUSE the
!     count is no longer 1, arguing that the loop was unnecessary.
!
!     Kept visible rather than deleted because of what it is: not a wrong
!     number, which a test can catch, but an argument against checking,
!     which nothing can - the reader who believes it stops looking.  The
!     other agent named this class today and proposed the remedy that
!     found it: after a retraction, sweep the TREE for the retracted
!     claim instead of marking the discussion.
!
      do i=1,ne0
        dploc(i)=0.d0
        if(evol(i).le.0.d0) cycle
!
!       The element's driving variable is the mean over ITS OWN
!       integration points.  The old code read slot 1 and nothing else,
!       which is the whole element only for a C3D4; on a C3D8 it would
!       have called one of eight points the element, and on a C3D20 one
!       of twenty-seven.  nipe(i) comes from the same table calcdamage
!       uses, so no slot is read that calcdamage never wrote.
!
        nipel=nipe(i)
        if(nipel.lt.1) nipel=1
        if(nipel.gt.mi(1)) nipel=mi(1)
        do ip1=1,nipel
          dploc(i)=dploc(i)+xstate(1,ip1,i)-xstateini(1,ip1,i)
        enddo
        dploc(i)=dploc(i)/dble(nipel)
      enddo
!
      do i=1,ne0
        dpsave(i)=0.d0
        if(evol(i).le.0.d0) cycle
        swv=0.d0
        sv=0.d0
        do k=nbstart(i),nbstart(i)+nbcount(i)-1
          j=nblist(k)
!
!         a deleted element carries no material any more, so it must not
!         weigh in; the normalisation below absorbs its loss
!
          if(ipkon(j).lt.0) cycle
          w=wgt(k)
          swv=swv+w*dploc(j)
          sv=sv+w
        enddo
        if(sv.gt.0.d0) then
          dpsave(i)=swv/sv
!
!         The self weight is exp(0)*evol(i)=evol(i), and sv is the
!         normalisation formed just above from the SURVIVING neighbours,
!         so this is d(ebar_i)/d(e_i) for the integral average, exactly.
!         It is recorded here rather than recomputed because the loop
!         that has both numbers is this one.
!
          chainf(i)=evol(i)/sv
        else
          dpsave(i)=dploc(i)
          chainf(i)=1.d0
        endif
      enddo
      chnbuilt=1
!
      return
      end
!
!     ==================================================================
!     FROZEN-LOCAL control (imodenl = 2).  NOT a regularisation.
!
!     WHY IT EXISTS.  With ell > 0 the driving variable is refreshed once
!     per calcdamagebase call - which nonlingeo.c makes only at the
!     predictor and at the commit - and damageupdatepoint then reads that
!     same value on every Newton iteration.  So the damage driver stops
!     responding to the trial state: at ell = 0 damage is integrated
!     IMPLICITLY inside Newton, at ell > 0 it changes only between
!     increments.  That is a change of SCHEME, not a lag.
!
!     Consequence: every local-vs-nonlocal comparison in this tree moves
!     TWO factors at once - the internal length and the integration
!     scheme.  E-84, E-110, E-114 and E-118 are all in that position, and
!     E-84's headline (spread 89.20% -> 8.56% on three bar meshes) cannot
!     be attributed to the length until the two are separated.
!
!     This mode is the missing arm: the SAME staggered update, no spatial
!     averaging, no length.  Run the three E-83 bars as
!     local / FROZEN / INTEGRAL / GRADIENT.  If FROZEN already delivers
!     most of the objectivity, then the length is not what delivered it.
!
!     It regularises nothing and must never be adopted as a model.
!     ==================================================================
!
      subroutine damfrozen(ipkon,kon,lakon,co,ne0,mi,xstate,
     &     xstateini,nstate_)
      use damnlmod
      implicit none
!
      character*8 lakon(*)
      integer ipkon(*),kon(*),ne0,mi(*),nstate_,i,indexf,nipf,iokf,ipf
      integer ionnl
      real*8 co(3,*),xstate(nstate_,mi(1),*),
     &     xstateini(nstate_,mi(1),*),xcf,ycf,zcf,volf
!
      call damnlactive(ionnl)
      if(ionnl.eq.0) return
!
      if(allocated(dpsave).and.(nesave.ne.ne0)) then
        deallocate(dpsave)
        if(allocated(dploc)) deallocate(dploc)
        if(allocated(chainf)) deallocate(chainf,iprobee)
        fbuilt=0
        chnbuilt=0
      endif
      if(.not.allocated(dpsave)) allocate(dpsave(ne0))
      if(.not.allocated(dploc)) allocate(dploc(ne0))
      if(.not.allocated(iprobee)) then
        allocate(chainf(ne0),iprobee(ne0))
        do i=1,ne0
          chainf(i)=1.d0
          iprobee(i)=0
        enddo
      endif
      nesave=ne0
      if(fbuilt.eq.0) then
        fbuilt=1
        write(*,*) '[DAMAGE NONLOCAL] backend FROZEN: CONTROL',
     &       ' ONLY - staggered update, NO averaging, NO length.'
        call flush(6)
      endif
!
!     mi(1) is the allocated integration-point stride; a C3D4 has exactly
!     one point, so index 1 is the only valid one here - the same reason
!     damnonlocal loops to 1 and not to mi(1).
!
      do i=1,ne0
        dploc(i)=0.d0
        dpsave(i)=0.d0
        if(ipkon(i).lt.0) cycle
!
!       Same generalisation as the integral backend: every volume family,
!       and the mean over the element's own integration points.  This arm
!       is a control and must stay comparable with the one it controls,
!       so it cannot keep the C3D4 assumption after the other dropped it.
!
        indexf=ipkon(i)
        if(indexf.lt.0) cycle
        call damnlelgeom(lakon(i),kon,co,indexf,xcf,ycf,zcf,volf,
     &       nipf,iokf)
        if(iokf.eq.0) cycle
        if(nipf.gt.mi(1)) nipf=mi(1)
        do ipf=1,nipf
          dploc(i)=dploc(i)+xstate(1,ipf,i)-xstateini(1,ipf,i)
        enddo
        dploc(i)=dploc(i)/dble(nipf)
        dpsave(i)=dploc(i)
      enddo
!
      return
      end
!
!     ------------------------------------------------------------------
!     Select the regularisation backend.  0 = integral (Pijaudier-Cabot /
!     Bazant, the original), 1 = implicit gradient (Peerlings / Geers),
!     2 = FROZEN-LOCAL control, which is not a backend at all.
!     ------------------------------------------------------------------
!
      subroutine damnonlocalmode(imodein)
      use damnlmod
      implicit none
      integer imodein
      imodenl=imodein
      return
      end
!
      subroutine damnonlocalmodeget(imodeout)
      use damnlmod
      implicit none
      integer imodeout
      imodeout=imodenl
      return
      end
!
!     ==================================================================
!     IMPLICIT GRADIENT regularisation.
!
!         ebar - div( ell^2 grad ebar ) = e
!
!     WHY THIS EXISTS ALONGSIDE THE INTEGRAL FORM.  The integral form
!     stores an explicit neighbour list, so its memory is
!     O(N * (2*ell/h)^3).  That is affordable on a uniform lattice and
!     catastrophic on a graded gmsh mesh, where the neighbourhood of one
!     hydride element holds thousands of tiny elements: measured on
!     demo_realistic_clusters3, 6.5 GB at ell=0.15, 20.8 GB at 0.30 and
!     58.4 GB at 0.45, against 31.6 GB of machine.  E-87 concluded from
!     the resolution side that no admissible ell existed for those
!     meshes; the memory is the other half of the same obstruction.  The
!     PDE form has no neighbour list at all and costs O(nnz) whatever the
!     grading.
!
!     WHY IT IS CHEAP TO ADD HERE.  Everything downstream reads the
!     regularised variable through ONE accessor, damnonlocalval, which
!     hands back one scalar per element.  The constitutive law, the
!     rank-1 damage tangent, terminal deletion, rollback and the whole
!     topology machinery are untouched.
!
!     DISCRETISATION.  Nodal ebar, linear tetrahedra, LUMPED mass:
!
!         A_ab = delta_ab V/4  +  ell^2 V (grad N_a . grad N_b)
!         f_a  = V/4 * e_elem
!
!     A is symmetric positive definite and mass dominated, so a
!     Jacobi-preconditioned CG converges in tens of iterations and no
!     sparse structure has to be built or handed to PARDISO: the
!     matrix-vector product is a loop over elements.  Deleted elements
!     (ipkon<0) simply do not contribute, which is the right boundary
!     condition on a growing crack - the average must not diffuse across
!     a surface that is no longer there.
!
!     The tangent stays LOCAL, exactly as in the integral form.  That is
!     the standard trade for this family; it costs Newton iterations, not
!     correctness.
!
!     NOT YET: ell is one global value here.  E-85 rejected exactly that
!     on the ladder because the rungs differ 5x in scale, and the fix is
!     to read ell from the material card - in this form it enters as a
!     per-element coefficient inside K, so that is a change to how ell2
!     is fetched in the two loops below and nothing else.
!     ==================================================================
!
      subroutine damgradient(ipkon,kon,lakon,co,ne0,mi,xstate,
     &     xstateini,nstate_,ielmat,dam)
      use damnlmod
      implicit none
!
      character*8 lakon(*)
      integer ipkon(*),kon(*),ne0,mi(*),nstate_,ielmat(mi(3),*)
      real*8 co(3,*),xstate(nstate_,mi(1),*),
     &     xstateini(nstate_,mi(1),*),dam(mi(1),*)
!
      integer i,j,a,b,indexe,it,maxit,nd,im,nc,nod(8),iokel,iokg
      integer ionnl,nskip,nptot,ip,nipel,ib,ndneg,nact
      real*8 ell2,s,rz,rzold,pap,alpha,beta,rnorm,rnorm0,tol,fm,
     &     dloc,gloc,ml(8),kel(8,8),cx,cy,cz,volel,wsum,elli,ellrep,
     &     dnwrst,dgmx,dgden,vsum,part,elli2
!
      call damnlactive(ionnl)
      if(ionnl.eq.0) return
      call damnlellinit()
      call damnlmatmap(ipkon,lakon,ielmat,ne0,mi)
      ell2=ellsave*ellsave
!
!     ---------------- geometry cache, reference configuration ---------
!
      if((gbuilt.ne.0).and.(nesave.ne.ne0)) gbuilt=0
      nskip=0
      if(gbuilt.eq.0) then
        if(allocated(gmass)) deallocate(gmass,vele,elnod,gnc,gkst)
        allocate(gmass(8,ne0),vele(ne0),elnod(8,ne0),gnc(ne0),
     &       gkst(ne0))
        if(allocated(nipe)) deallocate(nipe)
        allocate(nipe(ne0))
        nknl=0
        nptot=0
!
!       PASS ONE: masses, connectivity, and the size of the packed
!       conductivity.  An element damnlelmk refuses - a family that is
!       not a volume element, an inverted Jacobian - is counted and left
!       out of the assembly with gnc=0, which every loop below tests.
!       It is NOT silently averaged as if it were regularised, which is
!       what the C3D4-only version did to every hexahedron in the mesh.
!
        do i=1,ne0
          vele(i)=0.d0
          gnc(i)=0
          gkst(i)=0
          nipe(i)=0
          do j=1,8
            elnod(j,i)=0
            gmass(j,i)=0.d0
          enddo
          if(ipkon(i).lt.0) cycle
          indexe=ipkon(i)
          call damnlelmk(lakon(i),kon,co,indexe,nc,nod,ml,kel,iokel)
          if(iokel.eq.0) then
            nskip=nskip+1
            cycle
          endif
!
!         nipe comes from damnlelgeom, the same table calcdamage fills
!         its state from: the driving variable below averages the
!         element's integration points, and reading a slot calcdamage
!         never filled would be an average over noise.
!
          call damnlelgeom(lakon(i),kon,co,indexe,cx,cy,cz,volel,
     &         nipel,iokg)
          if(iokg.eq.0) then
            nskip=nskip+1
            cycle
          endif
          gnc(i)=nc
          nipe(i)=nipel
          vele(i)=volel
          gkst(i)=nptot
          nptot=nptot+nc*(nc+1)/2
          do j=1,nc
            elnod(j,i)=nod(j)
            gmass(j,i)=ml(j)
            nknl=max(nknl,nod(j))
          enddo
        enddo
!
!       PASS TWO: the packed upper triangles.
!
        nkpack=nptot
        if(allocated(gkpk)) deallocate(gkpk)
        allocate(gkpk(max(nptot,1)))
        do i=1,max(nptot,1)
          gkpk(i)=0.d0
        enddo
        do i=1,ne0
          nc=gnc(i)
          if(nc.le.0) cycle
          indexe=ipkon(i)
          call damnlelmk(lakon(i),kon,co,indexe,nc,nod,ml,kel,iokel)
          if(iokel.eq.0) then
            gnc(i)=0
            cycle
          endif
          ib=gkst(i)
          do b=1,nc
            do a=1,b
              gkpk(ib+a+b*(b-1)/2)=kel(a,b)
            enddo
          enddo
        enddo
        if(allocated(ebar)) deallocate(ebar,gdiag,grhs,gp,gap,gz,gr)
        if(allocated(gmsum)) deallocate(gmsum,ebar0)
        allocate(ebar(nknl),gdiag(nknl),grhs(nknl),gp(nknl),
     &       gap(nknl),gz(nknl),gr(nknl),gmsum(nknl),ebar0(nknl))
        do i=1,nknl
          ebar(i)=0.d0
        enddo
        if(.not.allocated(dpsave)) allocate(dpsave(ne0))
        if(.not.allocated(dploc)) allocate(dploc(ne0))
!
!       The gradient backend has no weight ratio to report; damnlchain
!       refuses for imodenl=1 and this keeps the flag honest as well.
!
        chnbuilt=0
        gbuilt=1
        nesave=ne0
        call damnlellmax(ellrep)
        write(*,*) '[DAMAGE NONLOCAL] gradient backend, ell=',ellrep,
     &       ' nodes=',nknl,' packed=',nkpack
!
!       The same participation number, ESTIMATED rather than counted:
!       this backend keeps no neighbour list, so how many elements sit
!       within 2*ell is the volume of that ball over the mean element
!       volume, floored at one because an element always counts itself.
!       Checked against the integral backend, which counts exactly: at
!       ell=0.2 on the reference deck this gives about 36 where the count
!       is 27.6 - the kernel is Gaussian and truncated, so the estimate
!       is high by construction - and at ell=0.01 both give one.  It is
!       used only to decide whether ANY neighbour joins, where a factor
!       of 1.3 cannot change the answer.
!
        vsum=0.d0
        nact=0
        do i=1,ne0
          if(gnc(i).le.0) cycle
          vsum=vsum+vele(i)
          nact=nact+1
        enddo
        if((nact.gt.0).and.(vsum.gt.0.d0).and.(ellrep.gt.0.d0)) then
          part=4.18879020478639d0*(2.d0*ellrep)**3
     &         /(vsum/dble(nact))
          if(part.lt.1.d0) part=1.d0
          call damnlresolve(part,ellrep,1)
        endif
        if(nskip.gt.0) then
          write(*,*) '*WARNING in damgradient: ',nskip,
     &         ' element(s) REFUSED and are NOT regularised:'
          write(*,*) '         not a volume family, or a non-positive'
          write(*,*) '         Jacobian.  Their damage stays LOCAL.'
          write(*,*) '         This is no longer the C3D4 restriction,'
          write(*,*) '         which was lifted: hexahedra, wedges and'
          write(*,*) '         the quadratic families are regularised.'
        endif
        if(ellmn.gt.0) then
          write(*,*) '[DAMAGE NONLOCAL] ell multiplier per material:',
     &         (ellmf(i),i=1,ellmn)
        endif
        if(nlocon.eq.1) then
          write(*,*) '[DAMAGE NONLOCAL] LOCALIZING g(D)=(1-R)(1-D)^n+R',
     &         ' n=',nlocn,' R=',nlocr,
     &         ' - the interaction radius collapses inside the band'
        endif
        call flush(6)
      endif
!
!     ---------------- local driving variable --------------------------
!
      do i=1,ne0
        dploc(i)=0.d0
        if(ipkon(i).lt.0) cycle
        if(gnc(i).le.0) cycle
!
!       ONE INTEGRATION POINT IS NOT THE ELEMENT.  This read used to be
!       xstate(1,1,i) alone, which is the element only when the element
!       has one point.  On C3D4 it does, which is why the line survived
!       the whole time the backend was tetrahedra-only - and why it would
!       have become wrong in the same commit that let a hexahedron in,
!       silently, by regularising one eighth of it.  The integral backend
!       carries the identical loop.
!
        nipel=nipe(i)
        if(nipel.lt.1) nipel=1
        if(nipel.gt.mi(1)) nipel=mi(1)
        do ip=1,nipel
          dploc(i)=dploc(i)+xstate(1,ip,i)-xstateini(1,ip,i)
        enddo
        dploc(i)=dploc(i)/dble(nipel)
      enddo
!
!     ---------------- per-element internal length ---------------------
!
!     Filled every call: g depends on the damage, which moves.
!
!     THE LENGTH COMES FROM damnlellel, NOT FROM ellsave.  ellsave is
!     only what the ENVIRONMENT set; a deck that puts NONLOCAL= on its
!     *DAMAGE INITIATION card leaves it at zero.  This loop used to read
!     ellsave directly, so on a card-only deck every ell2e came out
!     ZERO - the gradient equation degenerated to ebar = e and the
!     backend regularised NOTHING, while its banner printed a backend
!     name and the run looked healthy.  That is the same seam that was
!     closed for the integral backend in 18cc358 and 4f9ebea; the
!     gradient half was never closed, and nothing pointed at it because
!     no gate case runs this backend at all.
!
!     Found by reading the banner of the first hexahedral run this
!     backend was ever able to do: it said ell=0.0 while the card said
!     0.3.
!
      if(.not.allocated(ell2e)) allocate(ell2e(ne0))
      do i=1,ne0
        ell2e(i)=ell2
        if(ipkon(i).lt.0) cycle
!
!       THE ELEMENT'S OWN LENGTH, AND IT MUST SURVIVE TO THE ASSIGNMENT
!       BELOW.  This used to set ell2e here and then be overwritten two
!       dozen lines down by ell2e(i)=ell2*fm*fm*gloc, which carries the
!       GLOBAL ellsave - the environment's value.  So a deck whose length
!       comes from its *DAMAGE INITIATION card ran the whole backend at
!       ell2e=0, i.e. with no regularisation, while the banner printed
!       the card's length because the banner asked damnlellmax and the
!       assembly did not.
!
!       That is the seam I closed for the ARMING of the model and left
!       open for its USE, in the same commit.  It hid because the gate
!       case pins gradell from the BANNER - the number the model reports,
!       not the number it uses - so the case was green on a backend doing
!       nothing.  Found by reading this loop while answering a question
!       about which length to hand the crack-band width, not by a run.
!
        call damnlellel(i,elli)
        elli2=ell2
        if(elli.gt.0.d0) elli2=elli*elli
        fm=1.d0
        if(ellmn.gt.0) then
          im=ielmat(1,i)
          if((im.ge.1).and.(im.le.16)) fm=ellmf(im)
        endif
        gloc=1.d0
        if(nlocon.eq.1) then
          dloc=dam(1,i)-1.d0
          if(dloc.lt.0.d0) dloc=0.d0
          if(dloc.gt.1.d0) dloc=1.d0
          gloc=(1.d0-nlocr)*(1.d0-dloc)**nlocn+nlocr
        endif
        ell2e(i)=elli2*fm*fm*gloc
      enddo
!
!     ---------------- assemble diag(A) and the right-hand side --------
!
      do i=1,nknl
        gdiag(i)=0.d0
        grhs(i)=0.d0
      enddo
      do i=1,nknl
        gmsum(i)=0.d0
      enddo
      do i=1,ne0
        if(ipkon(i).lt.0) cycle
        nc=gnc(i)
        if(nc.le.0) cycle
        ib=gkst(i)
        do a=1,nc
          nd=elnod(a,i)
          gdiag(nd)=gdiag(nd)+gmass(a,i)
     &         +ell2e(i)*gkpk(ib+a+a*(a-1)/2)
          grhs(nd)=grhs(nd)+gmass(a,i)*dploc(i)
          gmsum(nd)=gmsum(nd)+gmass(a,i)
        enddo
      enddo
!
!     THE ZERO-LENGTH PROJECTION.  ebar0 is what this backend returns when
!     ell=0: the system is then M ebar = M e with M diagonal, so the answer
!     is available without a solve - it is the mass-weighted mean of the
!     elements meeting at the node.
!
!     It is NOT the identity, and that is the whole point.  dploc lives on
!     elements; ebar lives on nodes; the round trip element -> node ->
!     element averages twice and smooths over a distance set by the MESH,
!     not by ell.  Measured on this deck: at ell=0.01, thirty times smaller
!     than an element, max|dpsave-dploc|/max|dploc| was 0.55 to 0.83, and
!     the band width scaled with h (ratio 2.02 across a halving, against
!     1.11 along the band).  So the backend carried an internal length
!     nobody gave it.
!
      do i=1,nknl
        ebar0(i)=0.d0
        if(gmsum(i).gt.0.d0) ebar0(i)=grhs(i)/gmsum(i)
      enddo
!
!     a node with no live support carries no equation; keep the row
!     regular so the solve stays well posed and leave its value at zero
!
      do i=1,nknl
        if(gdiag(i).le.0.d0) then
          gdiag(i)=1.d0
          grhs(i)=0.d0
        endif
      enddo
!
!     ---------------- Jacobi-preconditioned CG, matrix free -----------
!
      rnorm0=0.d0
      do i=1,nknl
        rnorm0=rnorm0+grhs(i)*grhs(i)
      enddo
      rnorm0=dsqrt(rnorm0)
      if(rnorm0.le.0.d0) then
        do i=1,ne0
          dpsave(i)=0.d0
        enddo
        return
      endif
      tol=1.d-10*rnorm0
      maxit=400
!
!     start from the previous solution: between two Newton iterations the
!     field barely moves, so this typically halves the iteration count
!
      call damgradmv(ipkon,ne0,ebar,gap)
      do i=1,nknl
        gr(i)=grhs(i)-gap(i)
        gz(i)=gr(i)/gdiag(i)
        gp(i)=gz(i)
      enddo
      rz=0.d0
      do i=1,nknl
        rz=rz+gr(i)*gz(i)
      enddo
      do it=1,maxit
        rnorm=0.d0
        do i=1,nknl
          rnorm=rnorm+gr(i)*gr(i)
        enddo
        rnorm=dsqrt(rnorm)
        if(rnorm.le.tol) exit
        call damgradmv(ipkon,ne0,gp,gap)
        pap=0.d0
        do i=1,nknl
          pap=pap+gp(i)*gap(i)
        enddo
        if(dabs(pap).lt.1.d-300) exit
        alpha=rz/pap
        do i=1,nknl
          ebar(i)=ebar(i)+alpha*gp(i)
          gr(i)=gr(i)-alpha*gap(i)
        enddo
        rzold=rz
        rz=0.d0
        do i=1,nknl
          gz(i)=gr(i)/gdiag(i)
          rz=rz+gr(i)*gz(i)
        enddo
        if(dabs(rzold).lt.1.d-300) exit
        beta=rz/rzold
        do i=1,nknl
          gp(i)=gz(i)+beta*gp(i)
        enddo
      enddo
!
!     ---------------- back to the element -----------------------------
!
!     Back to the element, as the DIFFERENCE of two projections.
!
!         dpsave = dploc + Q( ebar - ebar0 )
!
!     Q is the mass-weighted mean of the element's nodal values, which is
!     int ebar dV / V; on a linear tetrahedron every mass is V/4 and it is
!     the plain quarter-sum the backend used to take.
!
!     WHY THE DIFFERENCE AND NOT ebar ITSELF.  Q(ebar) alone contains the
!     element -> node -> element round trip, which smooths over a distance
!     set by h and NOT by ell: at ell=0 it is exactly [1 2 1]/4 on a
!     uniform chain, a parasitic length of 0.707 h.  Subtracting Q(ebar0),
!     the same round trip applied to the same field at ell=0, removes that
!     term identically.  At ell=0 the correction vanishes and dpsave =
!     dploc EXACTLY, which is what the model requires and what this
!     backend did not do.
!
!     THIS IS A CONSTRUCTION, NOT THE SCHEME OF Peerlings et al. 1996.
!     Subtracting the zero-length projection changes the operator at EVERY
!     ell, not only at zero: what is solved is still their PDE, but what
!     reaches the damage law is the ell-dependent PART of its discrete
!     solution.  The justification is that the removed part is an artefact
!     of the discretisation rather than physics - it has no ell in it at
!     all - but it should be read as a deliberate choice made here, not as
!     an implementation of that paper.
!
!     One property worth naming because it is not obvious: both operators
!     preserve the total, sum Q(.) m = sum dploc m, so their difference
!     sums to zero and sum dpsave = sum dploc EXACTLY at any ell.  For a
!     quantity that enters an energy balance that is worth having.
!
!     POSITIVITY IS NOT ASSUMED.  The correction is a difference of two
!     averaging operators and is therefore sign-changing, so dpsave can go
!     negative beside a peak even though dploc never does.  It is clipped
!     below at zero - but a clip HIDES, so the clipping is COUNTED and
!     reported rather than done quietly.
!
      ndneg=0
      dnwrst=0.d0
      do i=1,ne0
        dpsave(i)=0.d0
        if(ipkon(i).lt.0) cycle
        nc=gnc(i)
        if(nc.le.0) cycle
        s=0.d0
        wsum=0.d0
        do a=1,nc
          nd=elnod(a,i)
          s=s+gmass(a,i)*(ebar(nd)-ebar0(nd))
          wsum=wsum+gmass(a,i)
        enddo
        if(wsum.le.0.d0) cycle
        dpsave(i)=dploc(i)+s/wsum
        if(dpsave(i).lt.0.d0) then
          ndneg=ndneg+1
          dnwrst=dmin1(dnwrst,dpsave(i))
          dpsave(i)=0.d0
        endif
      enddo
!
!     HOW FAR THE REGULARISED FIELD IS FROM THE LOCAL ONE.  This is the
!     number that makes the backend's own health visible: at an ell well
!     below the element size it must be SMALL, because there is nothing
!     for the regularisation to do.  Before the zero-length projection was
!     subtracted it was 0.55 to 0.83 at ell=0.01 on the reference deck -
!     the field handed to the damage law differed from the local one by
!     most of its own magnitude while the run looked healthy.  It is now
!     about 1e-3 there, which is the size of ell^2/h^2 and therefore the
!     regularisation itself rather than an artefact.
!
!     Printed only when the running maximum GROWS, so a long run says it
!     a handful of times instead of once per call.  A diagnostic that
!     prints every time stops being read, and we have spent a day on
!     things that were printed and not read.
!
      dgmx=0.d0
      dgden=0.d0
      do i=1,ne0
        if(gnc(i).le.0) cycle
        dgmx=dmax1(dgmx,dabs(dpsave(i)-dploc(i)))
        dgden=dmax1(dgden,dabs(dploc(i)))
      enddo
      if(dgden.gt.0.d0) then
        if(dgmx/dgden.gt.relmax*1.05d0+1.d-12) then
          relmax=dgmx/dgden
          write(*,*) '[DAMAGE NONLOCAL] gradient: |ebar-e|/|e| =',
     &         relmax,' at ell=',ellsave
        endif
      endif
      if(ndneg.gt.0) then
        if(ndneg.gt.ndnmax) then
          ndnmax=ndneg
          write(*,*) '[DAMAGE NONLOCAL] gradient: ',ndneg,
     &         ' element(s) clipped at zero, worst ',dnwrst
        endif
      endif
!
      return
      end
!
!     A*p for the operator  M_lumped + ell^2 K, assembled on the fly.
!
      subroutine damgradmv(ipkon,ne0,p,ap)
      use damnlmod
      implicit none
      integer ipkon(*),ne0,i,a,b,nd,nc,ib,ik
      real*8 p(*),ap(*),pa(8),s
!
      do i=1,nknl
        ap(i)=0.d0
      enddo
      do i=1,ne0
        if(ipkon(i).lt.0) cycle
        nc=gnc(i)
        if(nc.le.0) cycle
        ib=gkst(i)
        do a=1,nc
          pa(a)=p(elnod(a,i))
        enddo
        do a=1,nc
          s=0.d0
          do b=1,nc
            if(a.le.b) then
              ik=ib+a+b*(b-1)/2
            else
              ik=ib+b+a*(a-1)/2
            endif
            s=s+gkpk(ik)*pa(b)
          enddo
          nd=elnod(a,i)
          ap(nd)=ap(nd)+gmass(a,i)*pa(a)+ell2e(i)*s
        enddo
      enddo
      return
      end
