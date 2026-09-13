!
!     Element matrix for the UC6 zero-thickness cohesive triangle.
!
      subroutine e_c3d_uc6(co,kon,ipkon,s,sm,ff,nelem,vold,mi,
     &     ielprop,prop,dtime,xstateini,nstate_,mass,stiffness,rhsi,
     &     nmethod)
!
      implicit none
!
      integer kon(*),ipkon(*),nelem,mi(*),ielprop(*),nstate_
      integer mass(*),stiffness,rhsi,nmethod,indexe,indexp
      integer i,j,a,b,k,signi,mint
!
      real*8 co(3,*),s(60,60),sm(60,60),ff(60)
      real*8 vold(0:mi(2),*),prop(*),dtime
      real*8 xstateini(nstate_,mi(1),*)
      real*8 traction(3),ctan(3,3),area,rmat(3,3),deltal(3)
      real*8 deff,dmax,dback,dvisc,bmat(3,18),value,shape(3)
!
      if((nmethod.ne.1).and.(nmethod.ne.4)) then
         write(*,*) '*ERROR in e_c3d_uc6: UC6 currently supports'
         write(*,*) '       static and implicit dynamic procedures only'
         call exit(201)
      endif
!
      indexe=ipkon(nelem)
      indexp=ielprop(nelem)
      if(indexp.lt.0) then
         write(*,*) '*ERROR in e_c3d_uc6: missing *USER SECTION'
         write(*,*) '       for element ',nelem
         call exit(201)
      endif
!
      do i=1,18
         ff(i)=0.d0
         do j=1,18
            s(i,j)=0.d0
            sm(i,j)=0.d0
         enddo
      enddo
!
!     Three integration points, each with weight area/3.
!
      do mint=1,3
         call cohesive_uc6(co,kon,indexe,vold,mi,prop,indexp,dtime,
     &        xstateini,nstate_,nelem,mint,traction,ctan,area,rmat,
     &        shape,deltal,deff,dmax,dback,dvisc)
!
         do i=1,6
            if(i.le.3) then
               signi=-1
            else
               signi=1
            endif
            do k=1,3
               j=3*(i-1)+k
               do a=1,3
                  bmat(a,j)=signi*shape(mod(i-1,3)+1)*rmat(a,k)
               enddo
            enddo
         enddo
!
         if(stiffness.eq.1) then
            do i=1,18
               do j=1,18
                  value=0.d0
                  do a=1,3
                     do b=1,3
                        value=value+bmat(a,i)*ctan(a,b)*bmat(b,j)
                     enddo
                  enddo
                  s(i,j)=s(i,j)+(area/3.d0)*value
               enddo
            enddo
         endif
      enddo
!
!     No mass or distributed load contribution is defined for UC6.
!
      return
      end
