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
!     implicit gradient.  bgrad/vele/elnod cache the linear-tetrahedron
!     shape-function gradients, volumes and connectivity on the REFERENCE
!     configuration, for the same reason the neighbour list is built there:
!     the internal length is a material property and must not drift as the
!     mesh distorts.  nknl is the node count, derived from kon rather than
!     threaded through calcdamagebase and four call sites in nonlingeo.c.
!
      integer :: imodenl=0,gbuilt=0,nknl=0,fbuilt=0
      real*8, allocatable :: bgrad(:,:,:),vele(:),ebar(:),
     &     gdiag(:),grhs(:),gp(:),gap(:),gz(:),gr(:)
      integer, allocatable :: elnod(:,:)
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
      subroutine damnonlocalset(ellin)
      use damnlmod
      implicit none
      real*8 ellin
      ellsave=ellin
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
      real*8 val
      iok=0
      val=0.d0
      if(ellsave.le.0.d0) return
      if(.not.allocated(dpsave)) return
      if((iel.lt.1).or.(iel.gt.nesave)) return
      val=dpsave(iel)
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
     &     ncell,ix,iy,iz,icell,jcell,ifree
      real*8 co(3,*),xstate(nstate_,mi(1),*),
     &     xstateini(nstate_,mi(1),*),ell,
     &     xc,yc,zc,d2,w,swv,sv,det6,rmax,
     &     xmin,xmax,ymin,ymax,zmin,zmax,csize
!
      ell=ellsave
      if(ell.le.0.d0) return
      rmax=2.d0*ell
!
      if((nbuilt.ne.0).and.(nesave.ne.ne0)) nbuilt=0
!
      if(nbuilt.eq.0) then
        if(allocated(cen)) deallocate(cen,evol,wgt,dploc,dpsave,
     &       nbhead,nbnext,nblist,nbstart,nbcount)
        allocate(cen(3,ne0),evol(ne0),dploc(ne0),dpsave(ne0))
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
        do i=1,ne0
          evol(i)=0.d0
          cen(1,i)=0.d0
          cen(2,i)=0.d0
          cen(3,i)=0.d0
          if(lakon(i)(1:4).ne.'C3D4') cycle
          if(ipkon(i).lt.0) then
            indexe=-ipkon(i)-2
          else
            indexe=ipkon(i)
          endif
          if(indexe.lt.0) cycle
          xc=0.d0
          yc=0.d0
          zc=0.d0
          do j=1,4
            node=kon(indexe+j)
            xc=xc+co(1,node)
            yc=yc+co(2,node)
            zc=zc+co(3,node)
          enddo
          cen(1,i)=xc/4.d0
          cen(2,i)=yc/4.d0
          cen(3,i)=zc/4.d0
          det6=
     &   (co(1,kon(indexe+2))-co(1,kon(indexe+1)))*
     &  ((co(2,kon(indexe+3))-co(2,kon(indexe+1)))*
     &   (co(3,kon(indexe+4))-co(3,kon(indexe+1)))-
     &   (co(3,kon(indexe+3))-co(3,kon(indexe+1)))*
     &   (co(2,kon(indexe+4))-co(2,kon(indexe+1))))
     & -(co(2,kon(indexe+2))-co(2,kon(indexe+1)))*
     &  ((co(1,kon(indexe+3))-co(1,kon(indexe+1)))*
     &   (co(3,kon(indexe+4))-co(3,kon(indexe+1)))-
     &   (co(3,kon(indexe+3))-co(3,kon(indexe+1)))*
     &   (co(1,kon(indexe+4))-co(1,kon(indexe+1))))
     & +(co(3,kon(indexe+2))-co(3,kon(indexe+1)))*
     &  ((co(1,kon(indexe+3))-co(1,kon(indexe+1)))*
     &   (co(2,kon(indexe+4))-co(2,kon(indexe+1)))-
     &   (co(2,kon(indexe+3))-co(2,kon(indexe+1)))*
     &   (co(1,kon(indexe+4))-co(1,kon(indexe+1))))
          evol(i)=dabs(det6)/6.d0
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
                          nblist(ifree+nn)=j
                          wgt(ifree+nn)=dexp(-d2/(ell*ell))*evol(j)
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
        write(*,*) '[DAMAGE NONLOCAL] ell=',ell,' radius=',rmax,
     &       ' mean neighbours=',dble(ifree)/dble(max(1,ne0))
        call flush(6)
      endif
!
!     Local driving variable.  mi(1) is the ALLOCATED integration-point
!     stride - the maximum over every element type, 3 as soon as a UC6
!     facet exists - while a C3D4 has exactly ONE integration point
!     (calcdamage.f sets mint3d=1 for lakonl(4:4)=='4').  Looping to
!     mi(1) would average in uninitialised slots.  Only C3D4 reaches
!     here, evol being zero elsewhere, so the count is 1 by construction.
!
      do i=1,ne0
        dploc(i)=0.d0
        if(evol(i).le.0.d0) cycle
        dploc(i)=xstate(1,1,i)-xstateini(1,1,i)
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
        else
          dpsave(i)=dploc(i)
        endif
      enddo
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
      integer ipkon(*),kon(*),ne0,mi(*),nstate_,i
      real*8 co(3,*),xstate(nstate_,mi(1),*),
     &     xstateini(nstate_,mi(1),*)
!
      if(ellsave.le.0.d0) return
!
      if(allocated(dpsave).and.(nesave.ne.ne0)) then
        deallocate(dpsave)
        if(allocated(dploc)) deallocate(dploc)
        fbuilt=0
      endif
      if(.not.allocated(dpsave)) allocate(dpsave(ne0))
      if(.not.allocated(dploc)) allocate(dploc(ne0))
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
        if(lakon(i)(1:4).ne.'C3D4') cycle
        dploc(i)=xstate(1,1,i)-xstateini(1,1,i)
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
      integer i,j,a,n1,n2,n3,n4,indexe,it,maxit,nd,im
      real*8 ell2,det,dv,x1(3),e1(3),e2(3),e3(3),ji(3,3),
     &     s,rz,rzold,pap,alpha,beta,rnorm,rnorm0,tol,fm,dloc,gloc
!
      if(ellsave.le.0.d0) return
      call damnlellinit()
      ell2=ellsave*ellsave
!
!     ---------------- geometry cache, reference configuration ---------
!
      if((gbuilt.ne.0).and.(nesave.ne.ne0)) gbuilt=0
      if(gbuilt.eq.0) then
        if(allocated(bgrad)) deallocate(bgrad,vele,elnod)
        allocate(bgrad(3,4,ne0),vele(ne0),elnod(4,ne0))
        nknl=0
        do i=1,ne0
          vele(i)=0.d0
          elnod(1,i)=0
          elnod(2,i)=0
          elnod(3,i)=0
          elnod(4,i)=0
          if(ipkon(i).lt.0) cycle
          if(lakon(i)(1:4).ne.'C3D4') cycle
          indexe=ipkon(i)
          n1=kon(indexe+1)
          n2=kon(indexe+2)
          n3=kon(indexe+3)
          n4=kon(indexe+4)
          elnod(1,i)=n1
          elnod(2,i)=n2
          elnod(3,i)=n3
          elnod(4,i)=n4
          nknl=max(nknl,n1,n2,n3,n4)
          do j=1,3
            x1(j)=co(j,n1)
            e1(j)=co(j,n2)-x1(j)
            e2(j)=co(j,n3)-x1(j)
            e3(j)=co(j,n4)-x1(j)
          enddo
          det=e1(1)*(e2(2)*e3(3)-e2(3)*e3(2))
     &       -e1(2)*(e2(1)*e3(3)-e2(3)*e3(1))
     &       +e1(3)*(e2(1)*e3(2)-e2(2)*e3(1))
          if(dabs(det).lt.1.d-30) cycle
          vele(i)=dabs(det)/6.d0
!
!         the rows of J^-1 are grad(xi), grad(eta), grad(zeta), which for
!         N1=1-xi-eta-zeta, N2=xi, N3=eta, N4=zeta are grad N_2,3,4
!
          ji(1,1)=(e2(2)*e3(3)-e2(3)*e3(2))/det
          ji(1,2)=(e1(3)*e3(2)-e1(2)*e3(3))/det
          ji(1,3)=(e1(2)*e2(3)-e1(3)*e2(2))/det
          ji(2,1)=(e2(3)*e3(1)-e2(1)*e3(3))/det
          ji(2,2)=(e1(1)*e3(3)-e1(3)*e3(1))/det
          ji(2,3)=(e1(3)*e2(1)-e1(1)*e2(3))/det
          ji(3,1)=(e2(1)*e3(2)-e2(2)*e3(1))/det
          ji(3,2)=(e1(2)*e3(1)-e1(1)*e3(2))/det
          ji(3,3)=(e1(1)*e2(2)-e1(2)*e2(1))/det
          do j=1,3
            bgrad(j,2,i)=ji(1,j)
            bgrad(j,3,i)=ji(2,j)
            bgrad(j,4,i)=ji(3,j)
            bgrad(j,1,i)=-(ji(1,j)+ji(2,j)+ji(3,j))
          enddo
        enddo
        if(allocated(ebar)) deallocate(ebar,gdiag,grhs,gp,gap,gz,gr)
        allocate(ebar(nknl),gdiag(nknl),grhs(nknl),gp(nknl),
     &       gap(nknl),gz(nknl),gr(nknl))
        do i=1,nknl
          ebar(i)=0.d0
        enddo
        if(.not.allocated(dpsave)) allocate(dpsave(ne0))
        if(.not.allocated(dploc)) allocate(dploc(ne0))
        gbuilt=1
        nesave=ne0
        write(*,*) '[DAMAGE NONLOCAL] gradient backend, ell=',ellsave,
     &       ' nodes=',nknl
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
        if(vele(i).le.0.d0) cycle
        dploc(i)=xstate(1,1,i)-xstateini(1,1,i)
      enddo
!
!     ---------------- per-element internal length ---------------------
!
!     Filled every call: g depends on the damage, which moves.  With no
!     environment set this is ell^2 for every element and the assembly
!     below is arithmetically what it was.
!
      if(.not.allocated(ell2e)) allocate(ell2e(ne0))
      do i=1,ne0
        ell2e(i)=ell2
        if(ipkon(i).lt.0) cycle
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
        ell2e(i)=ell2*fm*fm*gloc
      enddo
!
!     ---------------- assemble diag(A) and the right-hand side --------
!
      do i=1,nknl
        gdiag(i)=0.d0
        grhs(i)=0.d0
      enddo
      do i=1,ne0
        if(ipkon(i).lt.0) cycle
        if(vele(i).le.0.d0) cycle
        dv=vele(i)*0.25d0
        do a=1,4
          nd=elnod(a,i)
          s=bgrad(1,a,i)**2+bgrad(2,a,i)**2+bgrad(3,a,i)**2
          gdiag(nd)=gdiag(nd)+dv+ell2e(i)*vele(i)*s
          grhs(nd)=grhs(nd)+dv*dploc(i)
        enddo
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
      do i=1,ne0
        dpsave(i)=0.d0
        if(ipkon(i).lt.0) cycle
        if(vele(i).le.0.d0) cycle
        dpsave(i)=0.25d0*(ebar(elnod(1,i))+ebar(elnod(2,i))
     &       +ebar(elnod(3,i))+ebar(elnod(4,i)))
        if(dpsave(i).lt.0.d0) dpsave(i)=0.d0
      enddo
!
      return
      end
!
!     A*p for the operator  M_lumped + ell^2 K, assembled on the fly.
!
      subroutine damgradmv(ipkon,ne0,p,ap)
      use damnlmod
      implicit none
      integer ipkon(*),ne0,i,a,nd
      real*8 p(*),ap(*),g(3),dv,c
!
      do i=1,nknl
        ap(i)=0.d0
      enddo
      do i=1,ne0
        if(ipkon(i).lt.0) cycle
        if(vele(i).le.0.d0) cycle
        g(1)=0.d0
        g(2)=0.d0
        g(3)=0.d0
        do a=1,4
          c=p(elnod(a,i))
          g(1)=g(1)+bgrad(1,a,i)*c
          g(2)=g(2)+bgrad(2,a,i)*c
          g(3)=g(3)+bgrad(3,a,i)*c
        enddo
        dv=vele(i)*0.25d0
        do a=1,4
          nd=elnod(a,i)
          ap(nd)=ap(nd)+dv*p(nd)
     &         +ell2e(i)*vele(i)*(bgrad(1,a,i)*g(1)+bgrad(2,a,i)*g(2)
     &         +bgrad(3,a,i)*g(3))
        enddo
      enddo
      return
      end
