!
!     Self test for damnlelmk: the element mass and conductivity of the
!     implicit-gradient equation
!
!         ebar - div( ell^2 grad ebar ) = e      (Peerlings et al. 1996)
!
!     for every volume family.  Exists because the gradient backend used
!     to hardwire both for a linear tetrahedron and SKIP every other
!     family - loudly since the tree said so, but skip it all the same,
!     leaving that damage local.  Generalising the two integrals is the
!     lift; this is what makes the lift checkable without a solve.
!
!     Every expectation below is written down in closed form, not read
!     off a run.  Three of them are family-independent identities and
!     they are the ones that would catch a wrong quadrature:
!
!       sum_a ml(a)     = V          (the shape functions are a partition
!                                     of unity, so the lumped masses sum
!                                     to the volume)
!       sum_b kel(a,b)  = 0          (a CONSTANT ebar has zero gradient,
!                                     so the conductivity must annihilate
!                                     it - this is the single strongest
!                                     check here: it fails for almost any
!                                     error in shp, xsj or the weights)
!       kel(a,b)        = kel(b,a)   (symmetry of the bilinear form)
!
      program elmktest
      implicit none
      integer kon(64),nc,nod(8),iok,nbad,i,j,nbadx
      real*8 co(3,64),ml(8),kel(8,8),tol,s
      real*8 mlq(8),kelq(8,8)
      character*8 lak
!
      nbad=0
      nbadx=0
      tol=1.d-11
!
!     ---------------- hexahedron: the unit cube -----------------------
!
!     N_1 = (1-x)(1-y)(1-z) on the unit cube, so
!       int N_1 dV = 1/8
!       int |grad N_1|^2 dV = 3 * int (1-y)^2 (1-z)^2 = 3 * 1/3 * 1/3
!                           = 1/3
!
      call cube(co,kon,1.d0,1.d0,1.d0)
      lak='C3D8    '
      call damnlelmk(lak,kon,co,0,nc,nod,ml,kel,iok)
      call chk('C3D8 unit cube  ',nc,8,iok,1,ml,kel,1.d0,tol,nbad)
      call lin('  linear field  ',nc,nod,co,kel,1.d0,1.d-12,nbad)
      call one('  ml(1)=1/8     ',ml(1),0.125d0,tol,nbad)
      call one('  kel(1,1)=1/3  ',kel(1,1),1.d0/3.d0,tol,nbad)
!
!     ---------------- hexahedron: a brick 2 x 3 x 5 -------------------
!
!     only the identities are asserted here: they must hold on a shape
!     that is not a cube, which is where a quadrature that happens to be
!     right on the unit cube by symmetry would part company.
!
      call cube(co,kon,2.d0,3.d0,5.d0)
      call damnlelmk(lak,kon,co,0,nc,nod,ml,kel,iok)
      call chk('C3D8 brick 2x3x5',nc,8,iok,1,ml,kel,30.d0,tol,nbad)
      call lin('  linear field  ',nc,nod,co,kel,30.d0,1.d-12,nbad)
!
!     ---------------- C3D20 on the SAME corners -----------------------
!
!     ebar is interpolated on the corners for every family, so a
!     quadratic element with these corners must give the SAME matrices,
!     to the last bit.  This is the decision of damnlelmk stated as a
!     test: if C3D20 ever grew its own branch, this goes red.
!
      lak='C3D20   '
      call damnlelmk(lak,kon,co,0,nc,nod,mlq,kelq,iok)
      call chk('C3D20 same corn.',nc,8,iok,1,mlq,kelq,30.d0,tol,nbad)
      call lin('  linear field  ',nc,nod,co,kelq,30.d0,1.d-12,nbad)
      s=0.d0
      do i=1,8
        s=s+dabs(mlq(i)-ml(i))
        do j=1,8
          s=s+dabs(kelq(i,j)-kel(i,j))
        enddo
      enddo
      call one('  = C3D8 exactly',s,0.d0,tol,nbad)
!
!     ---------------- tetrahedron: the corner tet ---------------------
!
!     (0,0,0) (1,0,0) (0,1,0) (0,0,1): V = 1/6, and for the LINEAR
!     tetrahedron the gradients are constant,
!       grad N_2 = (1,0,0), grad N_3 = (0,1,0), grad N_4 = (0,0,1),
!       grad N_1 = -(1,1,1)
!     so ml(a) = V/4 = 1/24, kel(1,1) = 3V = 1/2, kel(2,2) = V = 1/6,
!     kel(1,2) = -V = -1/6.
!
!     These four numbers are also the OLD tetrahedral path written out:
!     it used ml = V/4 and kel = V grad N . grad N by hand.  Agreement
!     here is what says the generalisation did not move the family that
!     already worked.
!
      call tet(co,kon)
      lak='C3D4    '
      call damnlelmk(lak,kon,co,0,nc,nod,ml,kel,iok)
      call chk('C3D4 corner tet ',nc,4,iok,1,ml,kel,1.d0/6.d0,tol,nbad)
      call lin('  linear field  ',nc,nod,co,kel,1.d0/6.d0,1.d-12,nbad)
      call one('  ml(1)=1/24    ',ml(1),1.d0/24.d0,tol,nbad)
      call one('  kel(1,1)=1/2  ',kel(1,1),0.5d0,tol,nbad)
      call one('  kel(2,2)=1/6  ',kel(2,2),1.d0/6.d0,tol,nbad)
      call one('  kel(1,2)=-1/6 ',kel(1,2),-1.d0/6.d0,tol,nbad)
!
!     ---------------- a SKEWED tetrahedron ----------------------------
!
!     Every closed form above sits on a SYMMETRIC Jacobian - a cube, a
!     brick, a corner tetrahedron - and that is exactly the geometry on
!     which a gradient read off the wrong index of J^-1 still comes out
!     right.  This case has none of those symmetries, which is why it
!     is here and why it is the one that caught the defect described
!     below.  V = 0.5765 from det = 3.459, computed by hand.
!
      call skewtet(co,kon)
      lak='C3D4    '
      call damnlelmk(lak,kon,co,0,nc,nod,ml,kel,iok)
      call chk('C3D4 skewed tet ',nc,4,iok,1,ml,kel,0.5765d0,tol,nbad)
      call lin('  linear field  ',nc,nod,co,kel,0.5765d0,1.d-12,nbad)
!
!     THE OLD CONVENTION, KEPT AS A GUARD.  The tetrahedral path this
!     replaces read grad(xi),grad(eta),grad(zeta) off the ROWS of the
!     inverted edge matrix.  They are the COLUMNS: with the edges as the
!     rows of J, x = x1 + J^T . (xi,eta,zeta), so the isoparametric
!     coordinates carry (J^-1)^T and their gradients are the columns.
!
!     The two agree whenever J is symmetric - which a corner tetrahedron
!     of a cube is, and which is why every closed form above passes
!     under both and why the error survived.  On a general tetrahedron
!     they do not: measured here, kel(1,1) = 0.844305 against 0.456206,
!     85 per cent high.
!
!     What makes this checkable at all is the LINEAR FIELD test, not the
!     row-sum identity: sum_b kel(a,b) = 0 holds for the wrong gradients
!     too, because they are still built to sum to zero.  A test can hold
!     an invariant and still miss the defect, which is the whole reason
!     this case exists.
!
      call oldrows(co,mlq,kelq)
      call lin('  old rows: FAIL',4,nod,co,kelq,0.5765d0,1.d-12,nbadx)
      if(nbadx.eq.1) then
        write(*,*) '  old convention fails the linear field, as it must'
      else
        write(*,*) '  old convention PASSES - the guard is broken'
        nbad=nbad+1
      endif
!
!     ---------------- wedge -------------------------------------------
!
!     triangle (0,0) (1,0) (0,1) extruded by 1: V = 1/2.
!
      call wedge(co,kon)
      lak='C3D6    '
      call damnlelmk(lak,kon,co,0,nc,nod,ml,kel,iok)
      call chk('C3D6 unit wedge ',nc,6,iok,1,ml,kel,0.5d0,tol,nbad)
      call lin('  linear field  ',nc,nod,co,kel,0.5d0,1.d-12,nbad)
!
!     ---------------- refusals ----------------------------------------
!
!     an INVERTED hexahedron.  The Jacobian goes negative, and the
!     element must be refused rather than integrated: a negative xsj
!     would give a negative mass, and the CG that consumes these
!     matrices would lose the property it converges on.
!
      call cube(co,kon,1.d0,1.d0,1.d0)
      call swapz(co)
      lak='C3D8    '
      call damnlelmk(lak,kon,co,0,nc,nod,ml,kel,iok)
      call ref('inverted C3D8   ',iok,nbad)
!
!     a family that is not a volume element at all
!
      call cube(co,kon,1.d0,1.d0,1.d0)
      lak='C3D2    '
      call damnlelmk(lak,kon,co,0,nc,nod,ml,kel,iok)
      call ref('unknown family  ',iok,nbad)
!
      if(nbad.eq.0) then
        write(*,*) 'damnlelmk: all cases agree with the closed form'
        call exit(0)
      endif
      write(*,*) 'damnlelmk: ',nbad,' case(s) FAILED'
      call exit(1)
      end
!
!     ------------------------------------------------------------------
!
      subroutine cube(co,kon,ax,ay,az)
      implicit none
      real*8 co(3,*),ax,ay,az
      integer kon(*),i
      co(1,1)=0.d0
      co(2,1)=0.d0
      co(3,1)=0.d0
      co(1,2)=ax
      co(2,2)=0.d0
      co(3,2)=0.d0
      co(1,3)=ax
      co(2,3)=ay
      co(3,3)=0.d0
      co(1,4)=0.d0
      co(2,4)=ay
      co(3,4)=0.d0
      do i=1,4
        co(1,i+4)=co(1,i)
        co(2,i+4)=co(2,i)
        co(3,i+4)=az
      enddo
      do i=1,64
        kon(i)=i
      enddo
      return
      end
!
      subroutine tet(co,kon)
      implicit none
      real*8 co(3,*)
      integer kon(*),i,j
      do i=1,4
        do j=1,3
          co(j,i)=0.d0
        enddo
      enddo
      co(1,2)=1.d0
      co(2,3)=1.d0
      co(3,4)=1.d0
      do i=1,64
        kon(i)=i
      enddo
      return
      end
!
      subroutine wedge(co,kon)
      implicit none
      real*8 co(3,*)
      integer kon(*),i
      co(1,1)=0.d0
      co(2,1)=0.d0
      co(3,1)=0.d0
      co(1,2)=1.d0
      co(2,2)=0.d0
      co(3,2)=0.d0
      co(1,3)=0.d0
      co(2,3)=1.d0
      co(3,3)=0.d0
      do i=1,3
        co(1,i+3)=co(1,i)
        co(2,i+3)=co(2,i)
        co(3,i+3)=1.d0
      enddo
      do i=1,64
        kon(i)=i
      enddo
      return
      end
!
!     turn the cube inside out by swapping its two faces
!
      subroutine swapz(co)
      implicit none
      real*8 co(3,*),t
      integer i
      do i=1,4
        t=co(3,i)
        co(3,i)=co(3,i+4)
        co(3,i+4)=t
      enddo
      return
      end
!
!     ------------------------------------------------------------------
!     chk asserts the three family-independent identities plus the
!     corner count, the status and the volume the masses sum to.
!
!     a tetrahedron with no right angles and no unit edges: volume 0.35
!     by construction, det = 2.1
!
      subroutine skewtet(co,kon)
      implicit none
      real*8 co(3,*)
      integer kon(*),i
      co(1,1)=0.3d0
      co(2,1)=-0.2d0
      co(3,1)=0.7d0
      co(1,2)=1.9d0
      co(2,2)=0.1d0
      co(3,2)=0.4d0
      co(1,3)=0.5d0
      co(2,3)=1.4d0
      co(3,3)=-0.2d0
      co(1,4)=0.8d0
      co(2,4)=0.6d0
      co(3,4)=1.6d0
      do i=1,64
        kon(i)=i
      enddo
      return
      end
!
!     the OLD backend's arithmetic, written out here and nowhere else:
!     rows of J^-1 are grad(xi), grad(eta), grad(zeta), which for
!     N1 = 1-xi-eta-ze are grad N_2,3,4, and grad N_1 = -sum of them.
!
      subroutine oldrows(co,ml,kel)
      implicit none
      real*8 co(3,*),ml(8),kel(8,8),x1(3),e1(3),e2(3),e3(3),ji(3,3),
     &     bg(3,4),det,vol
      integer i,j,a,b
      do j=1,3
        x1(j)=co(j,1)
        e1(j)=co(j,2)-x1(j)
        e2(j)=co(j,3)-x1(j)
        e3(j)=co(j,4)-x1(j)
      enddo
      det=e1(1)*(e2(2)*e3(3)-e2(3)*e3(2))
     &   -e1(2)*(e2(1)*e3(3)-e2(3)*e3(1))
     &   +e1(3)*(e2(1)*e3(2)-e2(2)*e3(1))
      vol=dabs(det)/6.d0
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
        bg(j,2)=ji(1,j)
        bg(j,3)=ji(2,j)
        bg(j,4)=ji(3,j)
        bg(j,1)=-(ji(1,j)+ji(2,j)+ji(3,j))
      enddo
      do a=1,4
        ml(a)=vol*0.25d0
        do b=1,4
          kel(a,b)=vol*(bg(1,a)*bg(1,b)+bg(2,a)*bg(2,b)
     &         +bg(3,a)*bg(3,b))
        enddo
      enddo
      return
      end
!
!     LINEAR FIELD (patch) test on the conductivity alone.
!
!     For a nodal field f_a = a . x_a the discrete gradient is
!     sum_a f_a grad N_a, which equals the constant a for ANY element
!     that reproduces linear fields.  Then
!
!         f^T K f = int |grad f|^2 dV = |a|^2 V
!
!     so one matrix-vector product decides linear completeness without
!     the gradients being exposed.  Three independent directions are
!     used because a wrong gradient can be right along one of them.
!
      subroutine lin(name,nc,nod,co,kel,vol,tol,nbad)
      implicit none
      character*16 name
      integer nc,nod(8),nbad,i,j,k
      real*8 co(3,*),kel(8,8),vol,tol,f(8),s,want,emax,av(3)
      emax=0.d0
      do k=1,3
        av(1)=0.d0
        av(2)=0.d0
        av(3)=0.d0
        av(k)=1.d0
        if(k.eq.3) then
          av(1)=0.6d0
          av(2)=-0.8d0
          av(3)=0.5d0
        endif
        do i=1,nc
          f(i)=av(1)*co(1,nod(i))+av(2)*co(2,nod(i))
     &         +av(3)*co(3,nod(i))
        enddo
        s=0.d0
        do i=1,nc
          do j=1,nc
            s=s+f(i)*kel(i,j)*f(j)
          enddo
        enddo
        want=(av(1)*av(1)+av(2)*av(2)+av(3)*av(3))*vol
        emax=dmax1(emax,dabs(s-want)/dmax1(1.d0,dabs(want)))
      enddo
      if(emax.le.tol) then
        write(*,'(a,a,a,e11.4)') ' ',name,'  max rel err ',emax
      else
        write(*,'(a,a,a,e11.4)') ' ',name,'  FAILED, rel err ',emax
        nbad=nbad+1
      endif
      return
      end
!
      subroutine chk(name,nc,ncw,iok,iokw,ml,kel,volw,tol,nbad)
      implicit none
      character*16 name
      integer nc,ncw,iok,iokw,nbad,bad,i,j
      real*8 ml(8),kel(8,8),volw,tol,s,rmax,asym,mmin
      bad=0
      if(nc.ne.ncw) bad=1
      if(iok.ne.iokw) bad=1
      s=0.d0
      mmin=1.d30
      rmax=0.d0
      asym=0.d0
      if(bad.eq.0) then
        do i=1,nc
          s=s+ml(i)
          mmin=dmin1(mmin,ml(i))
        enddo
        do i=1,nc
          rmax=dmax1(rmax,dabs(sum(kel(i,1:nc))))
          do j=1,nc
            asym=dmax1(asym,dabs(kel(i,j)-kel(j,i)))
          enddo
        enddo
        if(dabs(s-volw).gt.tol*dmax1(1.d0,dabs(volw))) bad=1
        if(rmax.gt.tol*dmax1(1.d0,dabs(volw))) bad=1
        if(asym.gt.tol*dmax1(1.d0,dabs(volw))) bad=1
        if(mmin.le.0.d0) bad=1
      endif
      if(bad.eq.0) then
        write(*,'(a,a,f13.8,a,e10.3,a,e10.3)') ' ',name,
     &       s,'  rowsum ',rmax,'  asym ',asym
      else
        write(*,'(a,a,a,f13.8,a,f13.8)') ' ',name,
     &       ' FAILED  sum(ml)=',s,'  expected ',volw
        nbad=nbad+1
      endif
      return
      end
!
      subroutine one(name,got,want,tol,nbad)
      implicit none
      character*16 name
      real*8 got,want,tol
      integer nbad
      if(dabs(got-want).le.tol*dmax1(1.d0,dabs(want))) then
        write(*,'(a,a,f16.10)') ' ',name,got
      else
        write(*,'(a,a,f16.10,a,f16.10)') ' ',name,got,
     &       '  FAILED, expected ',want
        nbad=nbad+1
      endif
      return
      end
!
      subroutine ref(name,iok,nbad)
      implicit none
      character*16 name
      integer iok,nbad
      if(iok.eq.0) then
        write(*,'(a,a,a)') ' ',name,' refused, as it must'
      else
        write(*,'(a,a,a)') ' ',name,' ACCEPTED - FAILED'
        nbad=nbad+1
      endif
      return
      end
