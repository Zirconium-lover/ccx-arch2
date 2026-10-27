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
      real*8 ddback,ddvisc,tol,initialgap,scale,gfl
!
!     B-K mixed mode (Turon et al. 2006, with the stiffness condition of
!     Turon et al. 2010): active when the section has 8 constants.
!
      integer ncon,ibk
      real*8 giic,etabk,bmix,lam0,xmb,hcur,hmax0,hmax,s2
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
!     CCX_UC6_FAILED_OPEN.  A facet whose degradation has reached the
!     floor (1-dvisc <= gmin) still carries gmin*Kn*delta in tension and
!     shear, a force that GROWS with the opening.  Over a large opening
!     that is real extra work, W/Gc of order gmin*delta^2/(delta_0*
!     delta_f) (library note notes/cohesive-uc6.md: +2.4 Gc at ten times
!     delta_f for Kn=1e5, Tn0=50, Gc=0.3, gmin=1e-3), and it is what let
!     fast-wrapped run on to theta=1 on the stretch of dead facets.  1
!     removes the floor: g = max(0,1-dvisc), so tension and shear
!     traction and tangent fall continuously to zero at delta_f and the
!     facet does exactly Gc; compression keeps the full-Kn contact.  0
!     (default) is the law as it was, bit for bit.
!
!     CONTINUOUS, AND THAT WAS MEASURED, NOT ASSUMED.  Two switched
!     versions came first.  Zeroing the traction at the floor but keeping
!     gmin*Kn in the tangent walled fast-wrapped at peak (theta 0.1075,
!     rc=201): the operator kept a stiffness the residual no longer had.
!     Zeroing both at the floor ran fast-wrapped to fracture, but walled
!     the notched plate of test/mech/size_effect before its peak (theta
!     0.302, gmin=1e-3): at the floor the traction jumps by gmin*Kn*delta,
!     about 1 per cent of Tn0, on many facets at once.  At gmin=1e-6 that
!     jump vanishes and the plate is the same with the switch as
!     without, which is what pointed at the jump.
!
      integer ifopen
      save ifopen
!
!     CCX_UC6_MIDSURFACE.  The local frame and the area are taken from
!     the REFERENCE minus face, so under a finite rotation of an open
!     facet the normal stays where the reference put it: a rigid 90
!     degree turn of an opened elastic facet reads its opening as shear
!     (test/mech/uc6_law: |RF| 56.25 where 100 is right, and a facet
!     with Ts0 > Tn0 is damaged by the turn alone).  1 takes both from
!     the CURRENT mid-surface x = X + (u+ + u-)/2, as Segurado and
!     LLorca 2004 (IJSS 41:2977, eqs. 6-12) do; like them the tangent
!     carries no derivative of the frame or of the area.  0 (default)
!     keeps the reference frame bit for bit.
!
      integer imid,imidinit
      save imid,imidinit
      data imidinit /0/
      character*8 cmid
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
!     Constants 7 and 8, G_IIc and eta, select the B-K law.  A section
!     with the six constants takes the original law below, statement for
!     statement, so every 6-constant deck is bit for bit as before.
!
      call uc6nconget(indexp,ncon)
      ibk=0
      if(ncon.ge.8) then
         giic=prop(indexp+7)
         etabk=prop(indexp+8)
         if((giic.le.0.d0).or.(etabk.le.0.d0)) then
            write(*,*) '*ERROR in cohesive_uc6: B-K needs G_IIc > 0'
            write(*,*) '       and eta > 0 (constants 7 and 8),'
            write(*,*) '       element ',nelem
            call exit(201)
         endif
         ibk=1
      endif
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
!     Reference local frame and area - or, with CCX_UC6_MIDSURFACE=1,
!     those of the current mid-surface.
!
      if(imidinit.eq.0) then
        imid=0
        call getenv('CCX_UC6_MIDSURFACE',cmid)
        if(cmid(1:1).eq.'1') imid=1
        imidinit=1
      endif
      if(imid.eq.0) then
      do k=1,3
         edge1(k)=co(k,kon(indexe+2))-co(k,kon(indexe+1))
         edge2(k)=co(k,kon(indexe+3))-co(k,kon(indexe+1))
      enddo
      else
      do k=1,3
         edge1(k)=co(k,kon(indexe+2))-co(k,kon(indexe+1))
     &     +0.5d0*(v(k,kon(indexe+2))+v(k,kon(indexe+5))
     &            -v(k,kon(indexe+1))-v(k,kon(indexe+4)))
         edge2(k)=co(k,kon(indexe+3))-co(k,kon(indexe+1))
     &     +0.5d0*(v(k,kon(indexe+3))+v(k,kon(indexe+6))
     &            -v(k,kon(indexe+1))-v(k,kon(indexe+4)))
      enddo
      endif
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
      if(ibk.eq.0) then
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
      else
!
!     B-K, as in the library note notes/cohesive-uc6.md (section 4) and
!     its reference notes/scripts/cohesive_ref.py, class BK:
!
!       kappa = Ks/Kn = (G_Ic/G_IIc)(Ts0/Tn0)^2      (Turon 2010: keeps the
!                                  pure-mode strengths AND toughnesses)
!       lambda = sqrt(<dn>^2 + kappa s^2),  B = kappa s^2 / lambda^2
!       lambda_0(B) = (Tn0/Kn) sqrt(1 + (G_IIc/G_Ic - 1) B^eta)
!       m = lambda_f/lambda_0 = 2 G_Ic Kn/Tn0^2 = df/d0, independent of B
!       h = max over the history of lambda/lambda_0(B)
!       D = m (h-1) / (h (m-1))
!
!     The history is the ratio h, not lambda, so D cannot fall when the
!     mode changes (Turon 2006 keep r = max lambda and evaluate the
!     thresholds at the current B, which lets D decrease).  State 1
!     stores h*d0: with G_IIc = G_Ic, lambda_0 = d0 and h*d0 = deff, so
!     the variable keeps its meaning, and the law reduces to the one
!     above to rounding.  beta carries kappa from here on, so the
!     traction and tangent below are shared with the original law.
!
      beta=(gc/giic)*(ts0/tn0)**2
      s2=deltal(2)**2+deltal(3)**2
      deff=dsqrt(max(deltal(1),0.d0)**2+beta*s2)
      if(deff.gt.1.d-300) then
         bmix=beta*s2/deff**2
      else
         bmix=0.d0
      endif
      lam0=d0*dsqrt(1.d0+(giic/gc-1.d0)*bmix**etabk)
      xmb=df/d0
      hcur=deff/lam0
      dmax0=max(0.d0,xstateini(1,mint,nelem))
      hmax0=dmax0/d0
      hmax=max(hmax0,hcur)
      dmax=hmax*d0
      if(hmax.le.1.d0) then
         dback=0.d0
      elseif(hmax.ge.xmb) then
         dback=1.d0
      else
         dback=xmb*(hmax-1.d0)/(hmax*(xmb-1.d0))
      endif
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
        ifopen=0
        call getenv('CCX_UC6_FAILED_OPEN',csm)
        if(csm(1:1).eq.'1') ifopen=1
        if(ifopen.eq.1) then
          write(*,'(a)')
     &      '[UC6 FAILED OPEN] no gmin floor: tension and shear fall'
     &      //' continuously to zero at delta_f, a failed facet carries'
     &      //' nothing; compression keeps the full-Kn contact.  This'
     &      //' CHANGES THE ANSWER after facets fail'
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
!     CCX_UC6_FAILED_OPEN: no floor.  g follows 1-dvisc to zero, so the
!     traction and the tangent of tension and shear fall CONTINUOUSLY to
!     zero at delta_f; compression keeps kn below.  gfl is the floor the
!     softening tangent below is gated on.
!
      gfl=gmin
      if(ifopen.eq.1) then
         g=max(0.d0,1.d0-dvisc)
         gfl=0.d0
      endif!
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
!
!     B-K: the same rank-one term with B frozen over the increment
!     (Turon 2006, eq. 49): dD/d(delta) = dD/dh * q/(lambda lambda_0),
!     so ddback below carries dD/dh / lambda_0 and grad stays q/lambda.
!
      if(ibk.eq.1) then
         if((hcur.ge.hmax0-1.d-12*max(1.d0,hmax0)).and.(hcur.gt.1.d0)
     &     .and.(hcur.lt.xmb).and.(deff.gt.1.d-30)
     &     .and.(1.d0-dvisc.gt.gfl)) then
            ddback=xmb/((xmb-1.d0)*hcur**2)/lam0
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
      elseif((deff.ge.dmax0-tol).and.(deff.gt.d0).and.(deff.lt.df)
     &     .and.(deff.gt.1.d-30).and.(1.d0-dvisc.gt.gfl)) then
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
!
!     ------------------------------------------------------------------
!
!     How many *USER SECTION constants a UC6 section has.  usersections.f
!     knows the count and the element knows only the start index into
!     prop, so the reader records (start, count) here and the law asks.
!     A start nobody recorded answers 6, the original law.
!
      subroutine uc6nconset(istart,n)
      implicit none
      integer istart,n,nsec,ist(4096),nct(4096)
      common /uc6ncon/ nsec,ist,nct
      data nsec /0/
      if(nsec.ge.4096) then
         write(*,*) '*ERROR in uc6nconset: more than 4096 sections'
         call exit(201)
      endif
      nsec=nsec+1
      ist(nsec)=istart
      nct(nsec)=n
      return
      end
!
      subroutine uc6nconget(istart,n)
      implicit none
      integer istart,n,nsec,ist(4096),nct(4096),i
      common /uc6ncon/ nsec,ist,nct
      n=6
      do i=1,nsec
         if(ist(i).eq.istart) then
            n=nct(i)
            return
         endif
      enddo
      return
      end
