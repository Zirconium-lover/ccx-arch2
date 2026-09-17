!
!     Self test for damnlelgeom: the centroid, volume and
!     integration-point count the nonlocal average uses for every volume
!     family.  Exists because the average is a weighted sum in which the
!     weight is the element VOLUME - a volume wrong by a constant factor
!     is invisible in a uniform mesh and silently reweights a graded one,
!     which is the mesh the whole regularisation is for.
!
!     Each case has an answer that can be written down without running
!     anything.
!
      program elgeomtest
      implicit none
      integer kon(64),nip,iok,nbad,nc
      real*8 co(3,64),cx,cy,cz,vol,tol
      character*8 lak
!
      nbad=0
      tol=1.d-12
!
!     unit cube as C3D8: volume 1, centre (0.5,0.5,0.5), 8 points
!
      call cube(co,kon,1.d0,1.d0,1.d0)
      lak='C3D8    '
      call damnlelgeom(lak,kon,co,0,cx,cy,cz,vol,nip,iok)
      call chk('C3D8 unit cube      ',vol,1.d0,cx,0.5d0,nip,8,iok,1,
     &     tol,nbad)
!
!     the same cube called C3D8R: same geometry, ONE point
!
      lak='C3D8R   '
      call damnlelgeom(lak,kon,co,0,cx,cy,cz,vol,nip,iok)
      call chk('C3D8R unit cube     ',vol,1.d0,cx,0.5d0,nip,1,iok,1,
     &     tol,nbad)
!
!     a brick 2 x 3 x 5: volume 30, centre (1,1.5,2.5)
!
      call cube(co,kon,2.d0,3.d0,5.d0)
      lak='C3D8    '
      call damnlelgeom(lak,kon,co,0,cx,cy,cz,vol,nip,iok)
      call chk('C3D8 brick 2x3x5    ',vol,30.d0,cy,1.5d0,nip,8,iok,1,
     &     tol,nbad)
!
!     C3D20 on the same eight corners: the corner hull is the same
!     volume, and the centre must NOT move - midside nodes are excluded
!     on purpose, and this is the case that catches it if they are not
!
      lak='C3D20   '
      call damnlelgeom(lak,kon,co,0,cx,cy,cz,vol,nip,iok)
      call chk('C3D20 same corners  ',vol,30.d0,cy,1.5d0,nip,27,iok,1,
     &     tol,nbad)
!
!     corner tetrahedron of the unit cube: volume 1/6
!
      call cube(co,kon,1.d0,1.d0,1.d0)
      kon(1)=1
      kon(2)=2
      kon(3)=4
      kon(4)=5
      lak='C3D4    '
      call damnlelgeom(lak,kon,co,0,cx,cy,cz,vol,nip,iok)
      call chk('C3D4 corner tet     ',vol,1.d0/6.d0,cx,0.25d0,nip,1,
     &     iok,1,tol,nbad)
!
!     wedge: half the unit cube cut on the x-y diagonal, volume 1/2
!
      call cube(co,kon,1.d0,1.d0,1.d0)
      kon(1)=1
      kon(2)=2
      kon(3)=3
      kon(4)=5
      kon(5)=6
      kon(6)=7
      lak='C3D6    '
      call damnlelgeom(lak,kon,co,0,cx,cy,cz,vol,nip,iok)
      call chk('C3D6 half cube      ',vol,0.5d0,cz,0.5d0,nip,2,iok,1,
     &     tol,nbad)
!
!     a family the routine does not know must REFUSE, not guess
!
      lak='C3D99   '
      call damnlelgeom(lak,kon,co,0,cx,cy,cz,vol,nip,iok)
      if(iok.ne.0) then
        write(*,*) ' unknown family      FAILED: iok=',iok,' expected 0'
        nbad=nbad+1
      else
        write(*,*) ' unknown family       ok   (refused, as it must)'
      endif
!
      write(*,*)
      if(nbad.eq.0) then
        write(*,*) '[ELGEOM SELFTEST] PASSED'
      else
        write(*,*) '[ELGEOM SELFTEST] ',nbad,' failure(s)'
        call exit(1)
      endif
      end
!
      subroutine cube(co,kon,a,b,c)
      implicit none
      real*8 co(3,*),a,b,c
      integer kon(*),i
      co(1,1)=0.d0
      co(2,1)=0.d0
      co(3,1)=0.d0
      co(1,2)=a
      co(2,2)=0.d0
      co(3,2)=0.d0
      co(1,3)=a
      co(2,3)=b
      co(3,3)=0.d0
      co(1,4)=0.d0
      co(2,4)=b
      co(3,4)=0.d0
      co(1,5)=0.d0
      co(2,5)=0.d0
      co(3,5)=c
      co(1,6)=a
      co(2,6)=0.d0
      co(3,6)=c
      co(1,7)=a
      co(2,7)=b
      co(3,7)=c
      co(1,8)=0.d0
      co(2,8)=b
      co(3,8)=c
      do i=1,8
        kon(i)=i
      enddo
      return
      end
!
      subroutine chk(name,vol,volw,ctr,ctrw,nip,nipw,iok,iokw,tol,nbad)
      implicit none
      character*20 name
      real*8 vol,volw,ctr,ctrw,tol
      integer nip,nipw,iok,iokw,nbad,bad
      bad=0
      if(dabs(vol-volw).gt.tol*dmax1(1.d0,dabs(volw))) bad=1
      if(dabs(ctr-ctrw).gt.tol*dmax1(1.d0,dabs(ctrw))) bad=1
      if(nip.ne.nipw) bad=1
      if(iok.ne.iokw) bad=1
      if(bad.eq.0) then
        write(*,'(a,a,f12.7,a,f10.6,a,i3)') ' ',name,vol,
     &       '  centre ',ctr,'  points ',nip
      else
        write(*,'(a,a,f12.7,a,f12.7,a)') ' ',name,vol,
     &       '  FAILED, expected volume ',volw,' (or centre/points)'
        nbad=nbad+1
      endif
      return
      end
