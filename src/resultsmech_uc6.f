!
!     Internal force, output and trial state for the UC6 cohesive element.
!
      subroutine resultsmech_uc6(co,kon,ipkon,ne,v,stx,xstateini,
     &     xstate,mi,nstate_,fn,qa,nal,calcul_fn,calcul_qa,iout,eei,
     &     dtime,nelem,ielprop,prop)
!
      implicit none
!
      integer kon(*),ipkon(*),ne,mi(*),nstate_,nal,calcul_fn
      integer calcul_qa,iout,nelem,ielprop(*)
      integer indexe,indexp,i,j,k,node,signi,mint
!
      real*8 co(3,*),v(0:mi(2),*),stx(6,mi(1),*)
      real*8 xstateini(nstate_,mi(1),*),xstate(nstate_,mi(1),*)
      real*8 fn(0:mi(2),*),qa(*),eei(6,mi(1),*),prop(*),dtime
      real*8 traction(3),ctan(3,3),area,rmat(3,3),deltal(3)
      real*8 deff,dmax,dback,dvisc,tglobal(3),q(3,6),force
      real*8 shape(3)
!
      indexe=ipkon(nelem)
      indexp=ielprop(nelem)
      if(indexp.lt.0) then
         write(*,*) '*ERROR in resultsmech_uc6: missing *USER SECTION'
         write(*,*) '       for element ',nelem
         call exit(201)
      endif
!
      if(calcul_qa.eq.1) then
         do i=1,6
            node=kon(indexe+i)
            do k=1,3
               q(k,i)=fn(k,node)
            enddo
         enddo
      endif
!
      do mint=1,3
         call cohesive_uc6(co,kon,indexe,v,mi,prop,indexp,dtime,
     &        xstateini,nstate_,nelem,mint,traction,ctan,area,rmat,
     &        shape,deltal,deff,dmax,dback,dvisc)
!
         xstate(1,mint,nelem)=dmax
         xstate(2,mint,nelem)=dvisc
         xstate(3,mint,nelem)=dback
         if(dback.ge.1.d0-1.d-12) then
            xstate(4,mint,nelem)=1.d0
         else
            xstate(4,mint,nelem)=0.d0
         endif
!
!        Standard slots are cohesive diagnostics: three tractions,
!        damage, effective separation and maximum separation.
!
         stx(1,mint,nelem)=traction(1)
         stx(2,mint,nelem)=traction(2)
         stx(3,mint,nelem)=traction(3)
         stx(4,mint,nelem)=dvisc
         stx(5,mint,nelem)=deff
         stx(6,mint,nelem)=dmax
!
         if(iout.gt.0) then
            eei(1,mint,nelem)=deltal(1)
            eei(2,mint,nelem)=deltal(2)
            eei(3,mint,nelem)=deltal(3)
            eei(4,mint,nelem)=dvisc
            eei(5,mint,nelem)=deff
            eei(6,mint,nelem)=dmax
         endif
!
         do k=1,3
            tglobal(k)=0.d0
            do j=1,3
               tglobal(k)=tglobal(k)+rmat(j,k)*traction(j)
            enddo
         enddo
!
         if(calcul_fn.eq.1) then
            do i=1,6
               if(i.le.3) then
                  signi=-1
               else
                  signi=1
               endif
               node=kon(indexe+i)
               do k=1,3
                  force=signi*(area/3.d0)*shape(mod(i-1,3)+1)*
     &                 tglobal(k)
                  fn(k,node)=fn(k,node)+force
               enddo
            enddo
         endif
      enddo
!
      if(calcul_qa.eq.1) then
         do i=1,6
            node=kon(indexe+i)
            do k=1,3
               qa(1)=qa(1)+dabs(fn(k,node)-q(k,i))
            enddo
         enddo
         nal=nal+18
      endif
!
      return
      end
