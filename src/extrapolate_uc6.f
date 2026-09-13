!
!     Linear extrapolation from the three triangular integration points.
!
      subroutine extrapolate_uc6(yi,yn,ipkon,inum,kon,nfield,mi,
     &     ndim,iorienloc,i)
!
      implicit none
!
      integer ipkon(*),inum(*),kon(*),nfield,mi(*),ndim,iorienloc,i
      integer indexe,j,k,node,m
      real*8 value
      real*8 yi(ndim,mi(1),*),yn(nfield,*)
!
      if(iorienloc.ne.0) then
         write(*,*) '*ERROR in extrapolate_uc6: local output is not'
         write(*,*) '       supported for UC6'
         call exit(201)
      endif
!
      indexe=ipkon(i)
      do j=1,6
         node=kon(indexe+j)
         do k=1,nfield
            value=0.d0
            do m=1,3
               if(m.eq.mod(j-1,3)+1) then
                  value=value+(5.d0/3.d0)*yi(k,m,i)
               else
                  value=value-(1.d0/3.d0)*yi(k,m,i)
               endif
            enddo
            yn(k,node)=yn(k,node)+value
         enddo
         inum(node)=inum(node)+1
      enddo
!
      return
      end
