!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2025 Guido Dhondt
!     
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
!     
!     
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of 
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
!     GNU General Public License for more details.
!     
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!     
      subroutine calcdamage(ipkon,lakon,kon,co,mi,
     &     thicke,ielmat,ielprop,prop,ne0,ndmat_,ntmat_,
     &     ndmcon,dmcon,dam,dtime,sti,ithermal,t1,xstate,
     &     xstateini,nstate_,vold,idamage,imode,alphaevent)
!
!     Compatibility wrapper preserving the stock CalculiX 2.23 ABI.
!     Existing callers keep the original accumulated-damage behaviour.
!     nonlingeo.c in the active-set patch calls calcdamagebase directly
!     and supplies the immutable physical-increment damage baseline.
!
      implicit none
      character*8 lakon(*)
      integer ipkon(*),kon(*),mi(*),ielmat(mi(3),*),ielprop(*),ne0,
     &     ndmcon(2,*),ndmat_,ntmat_,ithermal(*),nstate_,idamage,imode
      real*8 co(3,*),thicke(mi(3),*),prop(*),
     &     dmcon(0:ndmat_,ntmat_,*),dam(mi(1),*),dtime,
     &     sti(6,mi(1),*),t1(*),xstate(nstate_,mi(1),*),
     &     xstateini(nstate_,mi(1),*),vold(0:mi(2),*),alphaevent
!
      call calcdamagebase(ipkon,lakon,kon,co,mi,
     &     thicke,ielmat,ielprop,prop,ne0,ndmat_,ntmat_,
     &     ndmcon,dmcon,dam,dam,dtime,sti,ithermal,t1,xstate,
     &     xstateini,nstate_,vold,idamage,imode,alphaevent)
      return
      end
!
      subroutine calcdamagebase(ipkon,lakon,kon,co,mi,
     &     thicke,ielmat,ielprop,prop,ne0,ndmat_,ntmat_,
     &     ndmcon,dmcon,dam,dambase,dtime,sti,ithermal,t1,xstate,
     &     xstateini,nstate_,vold,idamage,imode,alphaevent)
!     
!     calculates the damage due to plastic strain
!     
      implicit none
!     
      character*8 lakon(*),lakonl
!     
      integer ipkon(*),kon(*),mi(*),nope,indexe,i,j,k,ii,
     &     konl(20),mint3d,jj,iflag,ki,kl,ilayer,nlayer,kk,
     &     nopes,ielmat(mi(3),*),mint2d,null,ielprop(*),ne0,
     &     ndmcon(2,*),ndmat_,ntmat_,i1,nopered,ithermal(*),
     &     nstate_,id,imat,ndmconst,idamage,imode,ielemcross,
     &     ide1,imodenlv
!     
      real*8 co(3,*),prop(*),xl(3,20),xi,et,ze,xsj,shp(4,20),weight,
     &     a,gs(8,4),dlayer(4),tlayer(4),thickness,skl(3,3),s(3,3),
     &     thicke(mi(3),*),xlayer(mi(3),4),shp2(7,8),xs2(3,7),xsj2(3),
     &     xl2(3,8),eps0RT,xlimit,d1,d2,d3,d4,d5,Tmelt,Ttrans,eps0p,shy,
     &     svm,triax,t1l,dmcon(0:ndmat_,ntmat_,*),dam(mi(1),*),
     &     dambase(mi(1),*),dpeq,dpeqdt,dtime,ef,
     &     sti(6,mi(1),*),t1(*),vold(0:mi(2),*),
     &     xstate(nstate_,mi(1),*),xstateini(nstate_,mi(1),*),
     &     dmconloc(ndmat_),That,alphaevent,damold,ddam,damtrial,
     &     alphai,ufail,charlen,det6v,ax,ay,az,bx,by,bz,cx,cy,cz,
     &     damagebaseval,damagecur,damageD,damagetarget,dpeqpost,
     &     de1tol,ellnl
      integer iokv
!
      real*8 dpnlv
!     
      include "gauss.f"
!
!     imode = 0: predictor only. The routine does not modify dam/ipkon.
!                idamage is the number of elements predicted to cross
!                the damage limit and alphaevent is the earliest estimated
!                crossing fraction within the current physical increment.
!
!     imode = 1: recompute the damage increment from dambase and perform
!                the usual hard element deletion.
!
!     dambase is the immutable damage state at the beginning of the
!     physical increment. Using it here is essential during same-load
!     re-equilibration: xstate-xstateini spans the complete physical
!     increment, so adding that increment to already-updated dam would
!     double count plastic strain.
!
!     The predictor uses a linear interpolation of the accumulated damage
!     increment. It is only an event locator; the physical increment is
!     re-solved by nonlingeo at the predicted load level before deletion.
!
!
!     Nonlocal average of the damage driving variable (E-83, E-84).
!     Inert unless CCX_DAMAGE_NONLOCAL has set an internal length, in
!     which case dpnl(i) replaces the LOCAL plastic-strain increment
!     further down.  The stress update stays local: only what DRIVES
!     damage is averaged.  That is the integral nonlocal DAMAGE model,
!     not nonlocal plasticity, and it is the whole change needed to put
!     a length into the formulation.
!
      call damnonlocalget(ellnl)
      if(ellnl.gt.0.d0) then
        call damnonlocalmodeget(imodenlv)
        if(imodenlv.eq.2) then
!
!         FROZEN-LOCAL control: the staggered update with NO averaging.
!         Separates the internal length from the integration scheme.
!
          call damfrozen(ipkon,kon,lakon,co,ne0,mi,xstate,
     &         xstateini,nstate_)
        elseif(imodenlv.eq.1) then
          call damgradient(ipkon,kon,lakon,co,ne0,mi,xstate,
     &         xstateini,nstate_,ielmat,dam)
        else
          call damnonlocal(ipkon,kon,lakon,co,ne0,mi,xstate,
     &         xstateini,nstate_)
        endif
      endif
      idamage=0
      alphaevent=2.d0
      de1tol=1.d-6
!     
      do i=1,ne0
        ielemcross=0
!     
!     element must exist and be a volume element
!     
        if((ipkon(i).lt.0).or.(lakon(i)(1:1).ne.'C')) cycle
!     
        lakonl=lakon(i)
        indexe=ipkon(i)
!     
        if(lakonl(1:5).eq.'C3D8I') then
          nope=11
        elseif(lakonl(4:4).eq.'2') then
          nope=20
        elseif(lakonl(4:4).eq.'8') then
          nope=8
        elseif(lakonl(4:5).eq.'10') then
          nope=10
        elseif(lakonl(4:4).eq.'4') then
          nope=4
        elseif(lakonl(4:5).eq.'15') then
          nope=15
        elseif(lakonl(4:5).eq.'6') then
          nope=6
        else
          cycle
        endif
!     
!     material
!     
        if(lakonl(7:8).ne.'LC') then
!     
!         no composite material: one material per element, all
!         integration points correspond to the same material
!     
          imat=ielmat(1,i)
          if(ndmcon(2,imat).eq.0) cycle
!     
!     determining the model for this element
!     
          if(int(dmcon(1,1,imat)).eq.1) then
!     
!     Rice-Tracey model
!     
            eps0RT=dmcon(2,1,imat)
            xlimit=dmcon(3,1,imat)
          elseif(int(dmcon(1,1,imat)).eq.2) then
!     
!     Johnson-Cook model
!     
            d1=dmcon(2,1,imat)
            d2=dmcon(3,1,imat)
            d3=dmcon(4,1,imat)
            d4=dmcon(5,1,imat)
            d5=dmcon(6,1,imat)
            Tmelt=dmcon(7,1,imat)
            Ttrans=dmcon(8,1,imat)
            eps0p=dmcon(9,1,imat)
            xlimit=dmcon(10,1,imat)
          elseif(int(dmcon(1,1,imat)).eq.3) then
!
!     DM2.0 tabulated ductile model
!
            xlimit=dmcon(2,1,imat)
          endif
        else
!     
!     composite materials
!     
!     determining the number of layers
!     
          nlayer=0
          do k=1,mi(3)
            if(ielmat(k,i).ne.0) then
              nlayer=nlayer+1
            endif
          enddo
!     
!     the thickness of the composite layers is only needed for
!     models requiring the temperature at the integration points
!     (so far only for the Johnson-Cook model)
!     
          if(int(dmcon(1,1,imat)).eq.2) then
            if(lakonl(4:4).eq.'2') then
              mint2d=4
              nopes=8
!     
!     determining the layer thickness and global thickness
!     at the shell integration points
!     
              iflag=1
              indexe=ipkon(i)
              do kk=1,mint2d
                xi=gauss3d2(1,kk)
                et=gauss3d2(2,kk)
                call shape8q(xi,et,xl2,xsj2,xs2,shp2,iflag)
                tlayer(kk)=0.d0
                do ii=1,nlayer
                  thickness=0.d0
                  do j=1,nopes
                    thickness=thickness+thicke(ii,indexe+j)*shp2(4,j)
                  enddo
                  tlayer(kk)=tlayer(kk)+thickness
                  xlayer(ii,kk)=thickness
                enddo
              enddo
              iflag=2
!     
              ilayer=0
              do ii=1,4
                dlayer(ii)=0.d0
              enddo
            elseif(lakonl(4:5).eq.'15') then
              mint2d=3
              nopes=6
!     
!     determining the layer thickness and global thickness
!     at the shell integration points
!     
              iflag=1
              indexe=ipkon(i)
              do kk=1,mint2d
                xi=gauss3d10(1,kk)
                et=gauss3d10(2,kk)
                call shape6tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
                tlayer(kk)=0.d0
                do ii=1,nlayer
                  thickness=0.d0
                  do j=1,nopes
                    thickness=thickness+thicke(ii,indexe+j)*shp2(4,j)
                  enddo
                  tlayer(kk)=tlayer(kk)+thickness
                  xlayer(ii,kk)=thickness
                enddo
              enddo
              iflag=2
!     
              ilayer=0
              do ii=1,3
                dlayer(ii)=0.d0
              enddo
            endif
          endif
!     
        endif
!     
        do j=1,nope
          konl(j)=kon(indexe+j)
          do k=1,3
            xl(k,j)=co(k,konl(j))
          enddo
        enddo
!
!       DE1 characteristic length for a linear tetrahedron.  The
!       coordinates in co are the reference mesh coordinates, so the
!       regularization length does not change during NLGEOM iterations.
!       det6v is six times the tetrahedral reference volume.
!
        charlen=0.d0
        if(lakonl(4:4).eq.'4') then
          ax=xl(1,2)-xl(1,1)
          ay=xl(2,2)-xl(2,1)
          az=xl(3,2)-xl(3,1)
          bx=xl(1,3)-xl(1,1)
          by=xl(2,3)-xl(2,1)
          bz=xl(3,3)-xl(3,1)
          cx=xl(1,4)-xl(1,1)
          cy=xl(2,4)-xl(2,1)
          cz=xl(3,4)-xl(3,1)
          det6v=dabs(ax*(by*cz-bz*cy)-ay*(bx*cz-bz*cx)
     &         +az*(bx*cy-by*cx))
          if(det6v.gt.1.d-30) charlen=det6v**(1.d0/3.d0)
        endif
!     
        if(lakonl(4:5).eq.'8R') then
          mint3d=1
        elseif(lakonl(4:7).eq.'20RB') then
          if((lakonl(8:8).eq.'R').or.(lakonl(8:8).eq.'C')) then
            mint3d=50
          else
            call beamintscheme(lakonl,mint3d,ielprop(i),prop,
     &           null,xi,et,ze,weight)
          endif
        elseif((lakonl(4:4).eq.'8').or.
     &         (lakonl(4:6).eq.'20R')) then
          if(lakonl(7:8).eq.'LC') then
            mint3d=8*nlayer
          else
            mint3d=8
          endif
        elseif(lakonl(4:4).eq.'2') then
          mint3d=27
        elseif(lakonl(4:5).eq.'10') then
          mint3d=4
        elseif(lakonl(4:4).eq.'4') then
          mint3d=1
        elseif(lakonl(4:5).eq.'15') then
          if(lakonl(7:8).eq.'LC') then
            mint3d=6*nlayer
          else
            mint3d=9
          endif
        elseif(lakonl(4:5).eq.'6') then
          mint3d=2
        else
          cycle
        endif
!     
        do jj=1,mint3d
!     
          if(lakonl(4:5).eq.'8R') then
            xi=gauss3d1(1,jj)
            et=gauss3d1(2,jj)
            ze=gauss3d1(3,jj)
            weight=weight3d1(jj)
          elseif(lakonl(4:7).eq.'20RB') then
            if((lakonl(8:8).eq.'R').or.(lakonl(8:8).eq.'C')) then
              xi=gauss3d13(1,jj)
              et=gauss3d13(2,jj)
              ze=gauss3d13(3,jj)
              weight=weight3d13(jj)
            else
              call beamintscheme(lakonl,mint3d,ielprop(i),prop,
     &             kk,xi,et,ze,weight)
            endif
          elseif((lakonl(4:4).eq.'8').or.
     &           (lakonl(4:6).eq.'20R'))
     &           then
            if(lakonl(7:8).ne.'LC') then
              xi=gauss3d2(1,jj)
              et=gauss3d2(2,jj)
              ze=gauss3d2(3,jj)
              weight=weight3d2(jj)
            else
              kl=mod(jj,8)
              if(kl.eq.0) kl=8
!     
              xi=gauss3d2(1,kl)
              et=gauss3d2(2,kl)
              ze=gauss3d2(3,kl)
              weight=weight3d2(kl)
!     
              ki=mod(jj,4)
              if(ki.eq.0) ki=4
!     
              if(kl.eq.1) then
                ilayer=ilayer+1
                if(ilayer.gt.1) then
                  do ii=1,4
                    dlayer(ii)=dlayer(ii)+xlayer(ilayer-1,ii)
                  enddo
                endif
              endif
              ze=2.d0*(dlayer(ki)+(ze+1.d0)/2.d0*xlayer(ilayer,ki))/
     &             tlayer(ki)-1.d0
              weight=weight*xlayer(ilayer,ki)/tlayer(ki)
              imat=ielmat(ilayer,i)
              if(ndmcon(2,imat).eq.0) cycle
            endif
          elseif(lakonl(4:4).eq.'2') then
            xi=gauss3d3(1,jj)
            et=gauss3d3(2,jj)
            ze=gauss3d3(3,jj)
            weight=weight3d3(jj)
          elseif(lakonl(4:5).eq.'10') then
            xi=gauss3d5(1,jj)
            et=gauss3d5(2,jj)
            ze=gauss3d5(3,jj)
            weight=weight3d5(jj)
          elseif(lakonl(4:4).eq.'4') then
            xi=gauss3d4(1,jj)
            et=gauss3d4(2,jj)
            ze=gauss3d4(3,jj)
            weight=weight3d4(jj)
          elseif(lakonl(4:5).eq.'15') then
            if(lakonl(7:8).ne.'LC') then
              xi=gauss3d8(1,jj)
              et=gauss3d8(2,jj)
              ze=gauss3d8(3,jj)
              weight=weight3d8(jj)
            else
              kl=mod(jj,6)
              if(kl.eq.0) kl=6
!     
              xi=gauss3d10(1,kl)
              et=gauss3d10(2,kl)
              ze=gauss3d10(3,kl)
              weight=weight3d10(kl)
!     
              ki=mod(jj,3)
              if(ki.eq.0) ki=3
!     
              if(kl.eq.1) then
                ilayer=ilayer+1
                if(ilayer.gt.1) then
                  do ii=1,3
                    dlayer(ii)=dlayer(ii)+xlayer(ilayer-1,ii)
                  enddo
                endif
              endif
              ze=2.d0*(dlayer(ki)+(ze+1.d0)/2.d0*xlayer(ilayer,ki))/
     &             tlayer(ki)-1.d0
              weight=weight*xlayer(ilayer,ki)/tlayer(ki)
              imat=ielmat(ilayer,i)
              if(ndmcon(2,imat).eq.0) cycle
            endif
          else
            xi=gauss3d7(1,jj)
            et=gauss3d7(2,jj)
            ze=gauss3d7(3,jj)
            weight=weight3d7(jj)
          endif
!     
!     determining the material model for this layer
!     (for composites only)
!     
          if(lakonl(7:8).eq.'LC') then
            if(int(dmcon(1,1,imat)).eq.1) then
!     
!     Rice-Tracey model
!     
              eps0RT=dmcon(2,1,imat)
              xlimit=dmcon(3,1,imat)
            elseif(int(dmcon(1,1,imat)).eq.2) then
!     
!     Johnson-Cook model
!     
              d1=dmcon(2,1,imat)
              d2=dmcon(3,1,imat)
              d3=dmcon(4,1,imat)
              d4=dmcon(5,1,imat)
              d5=dmcon(6,1,imat)
              Tmelt=dmcon(7,1,imat)
              Ttrans=dmcon(8,1,imat)
              eps0p=dmcon(9,1,imat)
              xlimit=dmcon(10,1,imat)
            elseif(int(dmcon(1,1,imat)).eq.3) then
!
!     DM2.0 tabulated ductile model
!
              xlimit=dmcon(2,1,imat)
            endif
          endif
!
!         DE1 is opt-in and is identified by the fourth Rice-Tracey
!         constant (u_f).  The first implementation is deliberately
!         restricted to C3D4 to keep the characteristic-length definition
!         unambiguous and to avoid changing stock behavior for other
!         elements.
!
          ide1=0
          if(((int(dmcon(1,1,imat)).eq.1).and.
     &       (ndmcon(1,imat).ge.4)).or.
     &       (int(dmcon(1,1,imat)).eq.3)) then
            ide1=1
            if(int(dmcon(1,1,imat)).eq.1) then
              ufail=dmcon(4,1,imat)
            else
              ufail=dmcon(3,1,imat)
            endif
            if(lakonl(4:4).ne.'4') then
              write(*,*) '*ERROR in calcdamage: DE1 displacement'
              write(*,*) '       evolution is presently implemented'
              write(*,*) '       for C3D4 elements only. Element: ',i
              call exit(201)
            endif
            if((ufail.le.0.d0).or.(charlen.le.0.d0)) then
              write(*,*) '*ERROR in calcdamage: invalid DE1 data'
              write(*,*) '       element=',i,' u_f=',ufail,
     &                   ' L=',charlen
              call exit(201)
            endif
          endif
!     
!     shape functions need only be determined if the     
!     temperature is needed, i.e. for the Johnson-Cook model
!     
          if(int(dmcon(1,1,imat)).eq.2) then
            iflag=1
            if(lakonl(1:5).eq.'C3D8R') then
              call shape8hr(xl,xsj,shp,gs,a)
            elseif(lakonl(1:5).eq.'C3D8I') then
              call shape8hu(xi,et,ze,xl,xsj,shp,iflag)
            elseif(nope.eq.20) then
              call shape20h(xi,et,ze,xl,xsj,shp,iflag)
            elseif(nope.eq.8) then
              call shape8h(xi,et,ze,xl,xsj,shp,iflag)
            elseif(nope.eq.10) then
              call shape10tet(xi,et,ze,xl,xsj,shp,iflag)
            elseif(nope.eq.4) then
              call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
            elseif(nope.eq.15) then
              call shape15w(xi,et,ze,xl,xsj,shp,iflag)
            else
              call shape6w(xi,et,ze,xl,xsj,shp,iflag)
            endif
          endif
!     
!     stress at the integration point
!     
          skl(1,1)=sti(1,jj,i)
          skl(2,2)=sti(2,jj,i)
          skl(3,3)=sti(3,jj,i)
          skl(1,2)=sti(4,jj,i)
          skl(1,3)=sti(5,jj,i)
          skl(2,3)=sti(6,jj,i)
!     
!     hydrostatic stress
!     
          shy=(skl(1,1)+skl(2,2)+skl(3,3))/3.d0
!     
!     deviatoric stress tensor
!     
          s(1,1)=skl(1,1)-shy
          s(2,2)=skl(2,2)-shy
          s(3,3)=skl(3,3)-shy
          s(1,2)=skl(1,2)
          s(1,3)=skl(1,3)
          s(2,3)=skl(2,3)
!     
!     von Mises stress
!     
          svm=dsqrt(3.d0/2.d0*(
     &         s(1,1)*s(1,1)+s(2,2)*s(2,2)+s(3,3)*s(3,3)+
     &         2.d0*(s(1,2)*s(1,2)+s(1,3)*s(1,3)+s(2,3)*s(2,3))))
!     
!     triaxiality
!     
          if(svm.gt.1.d-20) then
            triax=shy/svm
          else
            triax=0.d0
          endif
!     
!         calculate the temperature (only needed for the Johnson-Cook
!         model
!     
          if(int(dmcon(1,1,imat)).eq.2) then
            t1l=Ttrans
            if(ithermal(1).ge.1) then
              t1l=0.d0
              if(ithermal(1).eq.1) then
                if((lakonl(4:5).eq.'8 ').or.
     &               (lakonl(4:5).eq.'8I')) then
                  do i1=1,8
                    t1l=t1l+t1(konl(i1))/8.d0
                  enddo
                elseif(lakonl(4:6).eq.'20 ') then
                  nopered=20
                  call lintemp(t1,konl,nopered,jj,t1l)
                elseif(lakonl(4:6).eq.'10T') then
                  call linscal10(t1,konl,t1l,null,shp)
                else
                  do i1=1,nope
                    t1l=t1l+shp(4,i1)*t1(konl(i1))
                  enddo
                endif
              elseif(ithermal(1).ge.2) then
                if((lakonl(4:5).eq.'8 ').or.
     &               (lakonl(4:5).eq.'8I')) then
                  do i1=1,8
                    t1l=t1l+vold(0,konl(i1))/8.d0
                  enddo
                elseif(lakonl(4:6).eq.'20 ') then
                  nopered=20
                  call lintemp_th1(vold,konl,nopered,jj,t1l,mi)
                elseif(lakonl(4:6).eq.'10T') then
                  call linscal10(vold,konl,t1l,mi(2),shp)
                else
                  do i1=1,nope
                    t1l=t1l+shp(4,i1)*vold(0,konl(i1))
                  enddo
                endif
              endif
            endif
          endif
!     
!     interpolating the material data
!     (for models with temperature dependent parameters;
!     no such model is implemented so far; a model
!     number of 3 or higher is assumed)          
!     
          if(int(dmcon(1,1,imat)).gt.3) then
!     
!           number of constants in this model    
!     
            ndmconst=ndmcon(1,imat)
            if(ithermal(1).eq.0) then
              do k=1,ndmconst
                dmconloc(k)=dmcon(k,1,imat)
              enddo
            else
              call ident2(dmcon(0,1,imat),t1l,ndmcon(2,imat),ndmat_+1,
     &             id)
              if(ndmcon(2,imat).eq.1) then
                do k=1,ndmconst
                  dmconloc(k)=dmcon(k,1,imat)
                enddo
              elseif(id.eq.0) then
                do k=1,ndmconst
                  dmconloc(k)=dmcon(k,1,imat)
                enddo
              elseif(id.eq.ndmcon(2,imat)) then
                do k=1,ndmconst
                  dmconloc(k)=dmcon(k,id,imat)
                enddo
              else
                do k=1,ndmconst
                  dmconloc(k)=dmcon(k,id,imat)+
     &                 (dmcon(k,id+1,imat)-dmcon(k,id,imat))*
     &                 (t1l-dmcon(0,id,imat))/
     &                 (dmcon(0,id+1,imat)-dmcon(0,id,imat))
                enddo
              endif
            endif
          endif          
!     
!     change in equivalent plastic strain
!     
          dpeq=xstate(1,jj,i)-xstateini(1,jj,i)
          if(ellnl.gt.0.d0) then
            call damnonlocalval(i,dpnlv,iokv)
            if(iokv.eq.1) dpeq=dpnlv
          endif
          damagecur=dam(jj,i)
          damagebaseval=dambase(jj,i)

          if(dtime.gt.1.d-30) then
            dpeqdt=dpeq/dtime
          else
            dpeqdt=0.d0
          endif
!
!         ------------------------------------------------------------
!         DE1/DM2: initiation + displacement-based evolution
!         ------------------------------------------------------------
!
!         Storage convention (DE1 only, xlimit is required to be 1):
!           0 <= dam < 1 : Rice-Tracey initiation variable omega_D
!           1 <= dam <=2 : 1 + scalar degradation D
!
!         Thus existing DUCT output remains monotone and no new global
!         history array is required.  The physical degradation variable is
!         D=max(0,dam-1).  Fully damaged elements are not removed from the
!         topology in DE1; resultsmech applies a residual stiffness floor.
!
          if(ide1.eq.1) then

!           DE1.2 mode dispatch:
!             imode=0: legacy/A3 predictor only; DE1 is skipped.
!             imode=1: legacy/A3 apply only; DE1 is skipped because its
!                      trial state is already integrated inside Newton.
!             imode=2: DE1.2 trial update from the immutable physical-
!                      increment baseline. Legacy damage is skipped below.

            if(imode.ne.2) then
              cycle
            endif

!           Recompute from the immutable physical-increment baseline.
!           This makes DE1 transactional under Newton retries and A3
!           rollback, exactly like the A2/A3 initiation history.

            if(damagebaseval.ge.xlimit) then

!             Damage evolution was already active at increment start.

              damageD=damagebaseval-xlimit
              if(damageD.lt.0.d0) damageD=0.d0
              if(damageD.gt.1.d0) damageD=1.d0
              dpeqpost=dmax1(0.d0,dpeq)
              damageD=damageD+charlen*dpeqpost/ufail
              if(damageD.gt.1.d0) damageD=1.d0
              damagetarget=xlimit+damageD

            else

!             Still in the initiation stage.  The selected criterion
!             determines the fraction at which omega_D reaches
!             one.  Only the remaining plastic strain drives evolution.

              damagetarget=damagebaseval
              if(dpeq.gt.0.d0) then
                if(int(dmcon(1,1,imat)).eq.1) then
                  ef=1.65d0*eps0RT*dexp(-3.d0*triax/2.d0)
                elseif(int(dmcon(1,1,imat)).eq.3) then
                  call dm2failurestrain(dmcon,ndmat_,ntmat_,imat,
     &                 ndmcon,triax,ef)
                else
                  ef=-1.d0
                endif
                if(ef.gt.1.d-14) then
                  ddam=dpeq/ef
                  damtrial=damagebaseval+ddam
                  if(damtrial.lt.xlimit) then
                    damagetarget=damtrial
                  else
                    alphai=(xlimit-damagebaseval)/ddam
                    if(alphai.lt.0.d0) alphai=0.d0
                    if(alphai.gt.1.d0) alphai=1.d0
                    dpeqpost=(1.d0-alphai)*dpeq
                    damageD=charlen*dpeqpost/ufail
                    if(damageD.lt.0.d0) damageD=0.d0
                    if(damageD.gt.1.d0) damageD=1.d0
                    damagetarget=xlimit+damageD
                  endif
                else
                  write(*,*) '*WARNING in calcdamage: non-positive',
     &                       ' failure strain in element ',i,
     &                       ', integration point ',jj
                endif
              endif
            endif

            dam(jj,i)=damagetarget

!           idamage means "the constitutive damage state changed enough
!           to require another same-load equilibrium", not deletion, for
!           DE1.  Count an element only once even if it has several IPs.

            if(((damagetarget.gt.xlimit+de1tol).or.
     &          (damagecur.gt.xlimit+de1tol)).and.
     &          (dabs(damagetarget-damagecur).gt.de1tol)) then
              if(ielemcross.eq.0) idamage=idamage+1
              ielemcross=1
            endif

            cycle
          endif
!
!         ------------------------------------------------------------
!         Legacy CalculiX/A3 damage path (unchanged hard deletion)
!         ------------------------------------------------------------
!
!         imode=2 is reserved for the DE1.2 Newton-trial update.  Legacy
!         damage is evaluated only by the original predictor/apply modes.
!
          if(imode.eq.2) cycle
!
!         In apply/recompute mode the result at every surviving damage
!         integration point must be reconstructed from dambase.
!
          if(imode.ne.0) dam(jj,i)=dambase(jj,i)

          if(dpeq.gt.0.d0) then
!     
            if(int(dmcon(1,1,imat)).eq.1) then
!     
!     Rice-Tracey model
!     
              ef=1.65d0*eps0RT*dexp(-3.d0*triax/2.d0)
!     
            elseif(int(dmcon(1,1,imat)).eq.2) then
!     
!     Johnson-Cook model
!
              if(t1l.lt.Ttrans) then
                That=0.d0
              elseif(t1l.gt.Tmelt) then
                That=1.d0
              else
                That=(t1l-Ttrans)/(Tmelt-Ttrans)
              endif
!
              if(dpeqdt.gt.eps0p) then
                ef=(d1+d2*dexp(-d3*triax))*(1.d0+d4*dlog(dpeqdt/eps0p))*
     &               (1.d0+d5*That)
              else
!
!               replacing the logarithmic function by a linear function
!
                ef=(d1+d2*dexp(-d3*triax))*
     &               (1.d0+d4*(dpeqdt/eps0p-1.d0))*(1.d0+d5*That)
              endif
            endif
!     
!     damage increment
!     
            if(ef.gt.1.d-14) then
              damold=dambase(jj,i)
              ddam=dpeq/ef
              damtrial=damold+ddam
            else
              write(*,*) '*WARNING in calcdamage: non-positive',
     &                   ' failure strain in element ',i,
     &                   ', integration point ',jj
              cycle
            endif
!
!     predictor mode: determine whether the element crosses the damage
!     limit during this physical increment, without changing dam/ipkon.
!
            if(imode.eq.0) then
              if(damold.ge.xlimit) then
                alphai=0.d0
                if(ielemcross.eq.0) idamage=idamage+1
                ielemcross=1
                if(alphai.lt.alphaevent) alphaevent=alphai
              elseif((ddam.gt.1.d-30).and.
     &               (damtrial.ge.xlimit)) then
                alphai=(xlimit-damold)/ddam
                if(alphai.lt.0.d0) alphai=0.d0
                if(alphai.gt.1.d0) alphai=1.d0
                if(ielemcross.eq.0) idamage=idamage+1
                ielemcross=1
                if(alphai.lt.alphaevent) alphaevent=alphai
              endif
!
!     apply mode: preserve the original CalculiX hard-deletion law.
!
            else
              dam(jj,i)=damtrial
              if(dam(jj,i).ge.xlimit) then
                ipkon(i)=-ipkon(i)-2
                idamage=idamage+1
                exit
              endif
            endif
          endif
        enddo
      enddo
!     
      return
      end

!
!     BK1 integration-point update for progressive DE1/DM2 damage.
!     This routine is called immediately after the effective constitutive
!     update in resultsmech.  It deliberately implements only the existing
!     C3D4 displacement-evolution path; legacy Rice-Tracey/Johnson-Cook
!     hard deletion continues to use calcdamagebase.
!
      subroutine damageupdatepoint(iel,iint,ipkon,lakon,kon,co,mi,
     &     ielmat,ne0,ndmat_,ntmat_,ndmcon,dmcon,dam,dambase,
     &     stre,xstate,xstateini,nstate_)
      implicit none
      real*8 ellnlp,dpnlp
      integer ioknl
!
      character*8 lakon(*)
      integer iel,iint,ipkon(*),kon(*),mi(*),ielmat(mi(3),*),ne0,
     &     ndmat_,ntmat_,ndmcon(2,*),nstate_,imat,itype,indexe,
     &     n1,n2,n3,n4
      real*8 co(3,*),dmcon(0:ndmat_,ntmat_,*),dam(mi(1),*),
     &     dambase(mi(1),*),stre(6),xstate(nstate_,mi(1),*),
     &     xstateini(nstate_,mi(1),*),shy,s1,s2,s3,svm,triax,
     &     dpeq,xlimit,ufail,ef,damagebaseval,damagetarget,
     &     damageD,dpeqpost,ddam,damtrial,alphai,charlen,
     &     ax,ay,az,bx,by,bz,cx,cy,cz,det6v
!
      if((iel.lt.1).or.(iel.gt.ne0)) return
      if(ipkon(iel).lt.0) return
      if(lakon(iel)(1:1).ne.'C') return
      if(lakon(iel)(7:8).eq.'LC') then
        write(*,*) '*ERROR in damageupdatepoint: progressive damage'
        write(*,*) '       is not supported for composite elements.'
        call exit(201)
      endif
!
      imat=ielmat(1,iel)
      if(imat.le.0) return
      if(ndmcon(2,imat).eq.0) return
      itype=int(dmcon(1,1,imat))
      if(itype.eq.1) then
        if(ndmcon(1,imat).lt.4) return
        xlimit=dmcon(3,1,imat)
        ufail=dmcon(4,1,imat)
      elseif(itype.eq.3) then
        xlimit=dmcon(2,1,imat)
        ufail=dmcon(3,1,imat)
      else
        return
      endif
      if(lakon(iel)(4:4).ne.'4') then
        write(*,*) '*ERROR in damageupdatepoint: DE1 displacement'
        write(*,*) '       evolution is presently implemented for'
        write(*,*) '       C3D4 elements only. Element: ',iel
        call exit(201)
      endif
      if(ufail.le.0.d0) return
!
!     Reference characteristic length L=(6 V0)^(1/3).
!
      indexe=ipkon(iel)
      n1=kon(indexe+1)
      n2=kon(indexe+2)
      n3=kon(indexe+3)
      n4=kon(indexe+4)
      ax=co(1,n2)-co(1,n1)
      ay=co(2,n2)-co(2,n1)
      az=co(3,n2)-co(3,n1)
      bx=co(1,n3)-co(1,n1)
      by=co(2,n3)-co(2,n1)
      bz=co(3,n3)-co(3,n1)
      cx=co(1,n4)-co(1,n1)
      cy=co(2,n4)-co(2,n1)
      cz=co(3,n4)-co(3,n1)
      det6v=dabs(ax*(by*cz-bz*cy)-ay*(bx*cz-bz*cx)
     &     +az*(bx*cy-by*cx))
      if(det6v.le.1.d-30) return
      charlen=det6v**(1.d0/3.d0)
!
!     Triaxiality from the effective stress (scalar degradation would
!     leave it invariant, but using the effective stress is unambiguous).
!
      shy=(stre(1)+stre(2)+stre(3))/3.d0
      s1=stre(1)-shy
      s2=stre(2)-shy
      s3=stre(3)-shy
      svm=dsqrt(1.5d0*(s1*s1+s2*s2+s3*s3+
     &     2.d0*(stre(4)*stre(4)+stre(5)*stre(5)+
     &     stre(6)*stre(6))))
      if(svm.gt.1.d-20) then
        triax=shy/svm
      else
        triax=0.d0
      endif
!
!     Rebuild the complete trial state from the committed baseline.
!
      dpeq=xstate(1,iint,iel)-xstateini(1,iint,iel)
!
!     Nonlocal substitution.  THIS is the path that drives the response:
!     resultsmech calls damageupdatepoint per integration point inside
!     the stress update, and patching only calcdamagebase left the answer
!     bit-identical to the local model.
!
!     The field is refreshed once per calcdamagebase call, so what is read
!     here is the last refresh rather than the current trial state - the
!     explicit (staggered) nonlocal scheme.  The lag is one damage update
!     and the tangent is local either way, so nothing is given up that
!     was not already given up.
!
      call damnonlocalget(ellnlp)
      if(ellnlp.gt.0.d0) then
        call damnonlocalval(iel,dpnlp,ioknl)
        if(ioknl.eq.1) dpeq=dpnlp
      endif
      damagebaseval=dambase(iint,iel)
      damagetarget=damagebaseval
!
      if(damagebaseval.ge.xlimit) then
        damageD=damagebaseval-xlimit
        if(damageD.lt.0.d0) damageD=0.d0
        if(damageD.gt.1.d0) damageD=1.d0
        dpeqpost=dmax1(0.d0,dpeq)
        damageD=damageD+charlen*dpeqpost/ufail
        if(damageD.gt.1.d0) damageD=1.d0
        damagetarget=xlimit+damageD
      elseif(dpeq.gt.0.d0) then
        if(itype.eq.1) then
          ef=1.65d0*dmcon(2,1,imat)*dexp(-1.5d0*triax)
        else
          call dm2failurestrain(dmcon,ndmat_,ntmat_,imat,
     &         ndmcon,triax,ef)
        endif
        if(ef.gt.1.d-14) then
          ddam=dpeq/ef
          damtrial=damagebaseval+ddam
          if(damtrial.lt.xlimit) then
            damagetarget=damtrial
          else
            alphai=(xlimit-damagebaseval)/ddam
            if(alphai.lt.0.d0) alphai=0.d0
            if(alphai.gt.1.d0) alphai=1.d0
            dpeqpost=(1.d0-alphai)*dpeq
            damageD=charlen*dpeqpost/ufail
            if(damageD.lt.0.d0) damageD=0.d0
            if(damageD.gt.1.d0) damageD=1.d0
            damagetarget=xlimit+damageD
          endif
        endif
      endif
!
      dam(iint,iel)=damagetarget
      return
      end

!
!     DM2.0: linearly interpolate the tabulated fracture locus
!     eps_f(eta).  The table is stored as
!       dmcon(2)=xlimit, dmcon(3)=u_f,
!       dmcon(4:)=eta_1,eps_f_1,eta_2,eps_f_2,...
!     Values outside the tabulated eta range are clamped to the end points.
!
      subroutine dm2failurestrain(dmcon,ndmat_,ntmat_,imat,ndmcon,
     &     triax,ef)
      implicit none
      integer ndmat_,ntmat_,imat,ndmcon(2,*),npoints,ipt,
     &     idxeta,idxeps
      real*8 dmcon(0:ndmat_,ntmat_,*),triax,ef,eta1,eta2,
     &     eps1,eps2,fac
!
      npoints=(ndmcon(1,imat)-3)/2
      if(npoints.lt.2) then
        ef=-1.d0
        return
      endif
!
      eta1=dmcon(4,1,imat)
      eps1=dmcon(5,1,imat)
      if(triax.le.eta1) then
        ef=eps1
        return
      endif
!
      idxeta=2*npoints+2
      idxeps=idxeta+1
      eta2=dmcon(idxeta,1,imat)
      eps2=dmcon(idxeps,1,imat)
      if(triax.ge.eta2) then
        ef=eps2
        return
      endif
!
      do ipt=1,npoints-1
        idxeta=2*ipt+2
        idxeps=idxeta+1
        eta1=dmcon(idxeta,1,imat)
        eps1=dmcon(idxeps,1,imat)
        eta2=dmcon(idxeta+2,1,imat)
        eps2=dmcon(idxeps+2,1,imat)
        if((triax.ge.eta1).and.(triax.le.eta2)) then
          fac=(triax-eta1)/(eta2-eta1)
          ef=eps1+fac*(eps2-eps1)
          return
        endif
      enddo
!
      ef=-1.d0
      return
      end
