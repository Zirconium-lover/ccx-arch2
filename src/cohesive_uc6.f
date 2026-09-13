!
!     Six-node zero-thickness triangular cohesive element: local law
!
!     Node convention:
!       1,2,3 : minus side
!       4,5,6 : plus side, paired with 1,2,3
!
!     The reference normal is (x2-x1) x (x3-x1).  Positive normal
!     separation is therefore controlled by the ordering of nodes 1..3.
!
!     User-section properties:
!       1  Kn       normal penalty stiffness
!       2  Tn0      nominal normal strength
!       3  Ts0      nominal shear strength (same in both directions)
!       4  Gc       effective fracture energy
!       5  gmin     residual tensile/shear stiffness fraction
!       6  mu       viscous relaxation time (0 = inviscid)
!
!     State variables:
!       1  maximum effective separation
!       2  viscous damage
!       3  inviscid backbone damage
!       4  failure status (0/1)
!
      subroutine cohesive_uc6(co,kon,indexe,v,mi,prop,indexp,dtime,
     &     xstateini,nstate_,nelem,mint,traction,ctan,area,rmat,
     &     shape,deltal,deff,dmax,dback,dvisc)
!
      implicit none
!
      integer kon(*),indexe,mi(*),indexp,nstate_,nelem,mint
      integer i,j,k,nminus,nplus
!
      real*8 co(3,*),v(0:mi(2),*),prop(*),dtime
      real*8 xstateini(nstate_,mi(1),*)
      real*8 traction(3),ctan(3,3),rmat(3,3),deltal(3)
      real*8 shape(3)
      real*8 edge1(3),edge2(3),crossv(3),jump(3),qvec(3),grad(3)
      real*8 area,norm1,normc,kn,tn0,ts0,gc,gmin,mu,beta
      real*8 d0,df,deff,dmax,dmax0,dback,dvisc,dvisc0,alpha,g
      real*8 ddback,ddvisc,tol,initialgap,scale
!
!     CCX_UC6_CONTACT_SMOOTH.  The normal law is continuous at
!     deltal(1)=0 but its SLOPE jumps from g*kn to kn there, a factor of
!     1/gmin.  Measured on the fast wrapped deck at its wall: seven UC6
!     integration points sit on that kink and flip category at EVERY
!     line-search rung down to eps=6.1e-5, while the residual grows
!     strictly linearly in eps and never falls below its base value.
!     This blends the two slopes over a penetration band, so the normal
!     law becomes C1 and Newton has a differentiable problem.  Zero (the
!     default) reproduces the sharp law bit for bit.
!
      real*8 zsmooth,wsm,upen,psism,dpsism,ctan11
      character*256 csm
      save zsmooth
!
!     Diagnostic-only local rescue.  Acts on an explicit list of
!     element numbers and on nothing else: no material card, no
!     global damage law, no deletion logic, no controller.  Its
!     single purpose is the causal A/B around one node.
!
      integer nresc,iresc,ilist(64),idump,ndump,jdump
      integer idlist(64),iinitr,irescinc,uc6inc
      common /uc6ctl/ uc6inc
      real*8 gresc
      character*256 cresc
      save nresc,ilist,gresc,ndump,idlist,iinitr,irescinc
      data iinitr /0/
!
      if(nstate_.lt.4) then
         write(*,*) '*ERROR in cohesive_uc6: UC6 requires at least'
         write(*,*) '       four state variables (*DEPVAR 4)'
         call exit(201)
      endif
      if((mint.lt.1).or.(mint.gt.3)) then
         write(*,*) '*ERROR in cohesive_uc6: invalid integration point'
         write(*,*) '       element=',nelem,' point=',mint
         call exit(201)
      endif
!
      kn=prop(indexp+1)
      tn0=prop(indexp+2)
      ts0=prop(indexp+3)
      gc=prop(indexp+4)
      gmin=prop(indexp+5)
      mu=prop(indexp+6)
!
      if((kn.le.0.d0).or.(tn0.le.0.d0).or.(ts0.le.0.d0).or.
     &     (gc.le.0.d0).or.(gmin.lt.0.d0).or.(gmin.ge.1.d0).or.
     &     (mu.lt.0.d0)) then
         write(*,*) '*ERROR in cohesive_uc6: invalid *USER SECTION'
         write(*,*) '       properties for element ',nelem
         call exit(201)
      endif
!
      d0=tn0/kn
      df=2.d0*gc/tn0
      if(df.le.d0*(1.d0+1.d-10)) then
         write(*,*) '*ERROR in cohesive_uc6: fracture energy is too'
         write(*,*) '       small: delta_f must exceed delta_0'
         write(*,*) '       element=',nelem,' delta_0=',d0,
     &        ' delta_f=',df
         call exit(201)
      endif
!
!     Reference local frame and area.
!
      do k=1,3
         edge1(k)=co(k,kon(indexe+2))-co(k,kon(indexe+1))
         edge2(k)=co(k,kon(indexe+3))-co(k,kon(indexe+1))
      enddo
      crossv(1)=edge1(2)*edge2(3)-edge1(3)*edge2(2)
      crossv(2)=edge1(3)*edge2(1)-edge1(1)*edge2(3)
      crossv(3)=edge1(1)*edge2(2)-edge1(2)*edge2(1)
      norm1=dsqrt(edge1(1)**2+edge1(2)**2+edge1(3)**2)
      normc=dsqrt(crossv(1)**2+crossv(2)**2+crossv(3)**2)
      if((norm1.le.1.d-30).or.(normc.le.1.d-30)) then
         write(*,*) '*ERROR in cohesive_uc6: degenerate reference face'
         write(*,*) '       in element ',nelem
         call exit(201)
      endif
      area=0.5d0*normc
!
!     Rows of rmat map global vectors to (normal,shear-1,shear-2).
!
      do k=1,3
         rmat(1,k)=crossv(k)/normc
         rmat(2,k)=edge1(k)/norm1
      enddo
      rmat(3,1)=rmat(1,2)*rmat(2,3)-rmat(1,3)*rmat(2,2)
      rmat(3,2)=rmat(1,3)*rmat(2,1)-rmat(1,1)*rmat(2,3)
      rmat(3,3)=rmat(1,1)*rmat(2,2)-rmat(1,2)*rmat(2,1)
!
!     Enforce the intended zero-thickness paired-node topology.  The
!     constitutive separation is displacement jump only; an initial gap
!     would otherwise be silently ignored.
!
      initialgap=0.d0
      scale=dsqrt(area)
      do i=1,3
         nminus=kon(indexe+i)
         nplus=kon(indexe+i+3)
         do k=1,3
            initialgap=max(initialgap,dabs(co(k,nplus)-co(k,nminus)))
         enddo
      enddo
      if(initialgap.gt.max(1.d-10,1.d-8*scale)) then
         write(*,*) '*ERROR in cohesive_uc6: paired reference nodes'
         write(*,*) '       are not coincident in element ',nelem
         call exit(201)
      endif
!
!     Three-point triangular rule.  It integrates the linear jump field
!     and its quadratic stiffness work exactly, avoiding the hourglass
!     modes of a one-point zero-thickness triangle.
!
      do i=1,3
         shape(i)=1.d0/6.d0
      enddo
      shape(mint)=2.d0/3.d0
      do k=1,3
         jump(k)=0.d0
         do i=1,3
            nminus=kon(indexe+i)
            nplus=kon(indexe+i+3)
            jump(k)=jump(k)+shape(i)*(v(k,nplus)-v(k,nminus))
         enddo
      enddo
      do i=1,3
         deltal(i)=0.d0
         do k=1,3
            deltal(i)=deltal(i)+rmat(i,k)*jump(k)
         enddo
      enddo
!
!     A potential-based mixed-mode metric.  beta gives the requested
!     pure-shear strength while preserving a symmetric material tangent.
!
      beta=(ts0/tn0)**2
      deff=dsqrt(max(deltal(1),0.d0)**2+
     &     beta*(deltal(2)**2+deltal(3)**2))
      dmax0=max(0.d0,xstateini(1,mint,nelem))
      dmax=max(dmax0,deff)
!
      if(dmax.le.d0) then
         dback=0.d0
      elseif(dmax.ge.df) then
         dback=1.d0
      else
         dback=df*(dmax-d0)/(dmax*(df-d0))
      endif
!
!     Backward-Euler viscous regularization from the committed baseline.
!
      dvisc0=max(0.d0,min(1.d0,xstateini(2,mint,nelem)))
      if(mu.gt.0.d0) then
         alpha=max(dtime,0.d0)/(mu+max(dtime,0.d0))
      else
         alpha=1.d0
      endif
      dvisc=dvisc0+alpha*(dback-dvisc0)
      dvisc=max(dvisc0,min(1.d0,dvisc))
      g=max(gmin,1.d0-dvisc)
!
      if(iinitr.eq.0) then
        iinitr=1
        nresc=0
        ndump=0
        gresc=0.d0
        irescinc=0
        zsmooth=0.d0
        call getenv('CCX_UC6_CONTACT_SMOOTH',csm)
        if(csm(1:1).ne.' ') read(csm,*,err=9103,end=9103) zsmooth
 9103   continue
        if(zsmooth.lt.0.d0) zsmooth=0.d0
        if(zsmooth.gt.0.d0) then
          write(*,'(a,e12.5,a)')
     &      '[UC6 CONTACT SMOOTH] the normal slope is blended from g*kn'
     &      //' to kn over a penetration band of ',zsmooth,
     &      ' x d0.  The law stays continuous and monotone; only its'
     &      //' derivative changes, and only within that band.'
        endif
        call getenv('CCX_UC6_RESCUE_INC',cresc)
        if(cresc(1:1).ne.' ') read(cresc,*,err=9102,end=9102) irescinc
 9102   continue
        call getenv('CCX_UC6_RESCUE_G',cresc)
        if(cresc(1:1).ne.' ') read(cresc,*,err=9101,end=9101) gresc
 9101   continue
        call getenv('CCX_UC6_RESCUE_ELEMS',cresc)
        if(cresc(1:1).ne.' ') then
          call uc6intlist(cresc,ilist,64,nresc)
        endif
        call getenv('CCX_UC6_DUMP_ELEMS',cresc)
        if(cresc(1:1).ne.' ') then
          call uc6intlist(cresc,idlist,64,ndump)
        endif
        if(nresc.gt.0) then
          write(*,'(a,i4,a,e12.5)')
     &      '[UC6 RESCUE] active on ',nresc,
     &      ' element(s), g_rescue = ',gresc
          write(*,'(a,64(1x,i7))') '[UC6 RESCUE] elements:',
     &      (ilist(iresc),iresc=1,nresc)
        endif
      endif
!
!     Gated by increment so the control and the intervention share
!     the same history.  The first attempt applied the rescue from
!     increment one, which changed the crack path from the start and
!     made the comparison meaningless.
!
      if((nresc.gt.0).and.(uc6inc.ge.irescinc)) then
        do iresc=1,nresc
          if(nelem.eq.ilist(iresc)) then
            g=max(g,gresc)
          endif
        enddo
      endif
!
!     Local traction.  Compression retains the full normal penalty;
!     shear and positive normal traction are degraded by the same damage.
!
      wsm=zsmooth*d0
      if(wsm.gt.0.d0) then
!
!        psi(u), u = -deltal(1) the penetration: 0 for u<=0, u^2/(2 wsm)
!        on [0,wsm], u-wsm/2 beyond.  psi and psi' are both continuous,
!        psi' rises monotonically from 0 to 1, so the normal slope goes
!        from g*kn to kn without ever exceeding kn.
!
         upen=-deltal(1)
         if(upen.le.0.d0) then
            psism=0.d0
            dpsism=0.d0
         elseif(upen.lt.wsm) then
            psism=upen*upen/(2.d0*wsm)
            dpsism=upen/wsm
         else
            psism=upen-0.5d0*wsm
            dpsism=1.d0
         endif
         traction(1)=g*kn*deltal(1)-(1.d0-g)*kn*psism
         ctan11=g*kn+(1.d0-g)*kn*dpsism
      elseif(deltal(1).lt.0.d0) then
         traction(1)=kn*deltal(1)
         ctan11=kn
      else
         traction(1)=g*kn*deltal(1)
         ctan11=g*kn
      endif
      traction(2)=g*kn*beta*deltal(2)
      traction(3)=g*kn*beta*deltal(3)
!
!     Elastic/unloading tangent.
!
      do i=1,3
         do j=1,3
            ctan(i,j)=0.d0
         enddo
      enddo
      ctan(1,1)=ctan11
      ctan(2,2)=g*kn*beta
      ctan(3,3)=g*kn*beta
!
!     Consistent softening term for an advancing inviscid maximum.
!
      tol=1.d-12*max(1.d0,dmax0)
      if((deff.ge.dmax0-tol).and.(deff.gt.d0).and.(deff.lt.df)
     &     .and.(deff.gt.1.d-30).and.(1.d0-dvisc.gt.gmin)) then
         ddback=df*d0/((df-d0)*deff**2)
         ddvisc=alpha*ddback
         qvec(1)=max(deltal(1),0.d0)
         qvec(2)=beta*deltal(2)
         qvec(3)=beta*deltal(3)
         do j=1,3
            grad(j)=qvec(j)/deff
         enddo
         do i=1,3
            do j=1,3
               ctan(i,j)=ctan(i,j)-kn*qvec(i)*ddvisc*grad(j)
            enddo
         enddo
      endif
!
      if(ndump.gt.0) then
        do jdump=1,ndump
          if(nelem.eq.idlist(jdump)) then
            write(*,'(a,i7,a,i2,a,e12.5,a,e12.5)')
     &        '[UC6 TANGENT] elem ',nelem,' ip ',mint,
     &        '  dvisc=',dvisc,'  g=',g
            do idump=1,3
              write(*,'(a,3(1x,e15.8))') '[UC6 TANGENT]   ',
     &          ctan(idump,1),ctan(idump,2),ctan(idump,3)
            enddo
            write(*,'(a,3(1x,e12.5))') '[UC6 TANGENT]   delta=',
     &        deltal(1),deltal(2),deltal(3)
          endif
        enddo
      endif
!
      return
      end
!
      subroutine uc6intlist(str,list,nmax,n)
!
!     Parse a comma or space separated list of element numbers.
!
      implicit none
      character*256 str
      integer list(*),nmax,n,i,ival,ic
      n=0
      ival=0
      ic=0
      do i=1,256
        if((str(i:i).ge.'0').and.(str(i:i).le.'9')) then
          ival=10*ival+(ichar(str(i:i))-ichar('0'))
          ic=1
        else
          if(ic.eq.1) then
            if(n.lt.nmax) then
              n=n+1
              list(n)=ival
            endif
          endif
          ival=0
          ic=0
        endif
      enddo
      if(ic.eq.1) then
        if(n.lt.nmax) then
          n=n+1
          list(n)=ival
        endif
      endif
      return
      end
!
      subroutine uc6setinc(i)
!
!     Publishes the current increment to the cohesive law so a
!     diagnostic intervention can be switched on part way through a run.
!
      implicit none
      integer i,uc6inc
      common /uc6ctl/ uc6inc
      uc6inc=i
      return
      end
