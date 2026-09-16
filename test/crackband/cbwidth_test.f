!     Unit test of the CB1 projected crack-band width.
!     Each case has a hand-computable answer, printed next to the
!     measured one, so the geometry is checked rather than asserted.
      program cbtest
      implicit none
      real*8 xl(3,20),stre(6),w,r2,r3,vol,legacy
      integer nope,iok,i
      character*8 lakonl
      r2=dsqrt(2.d0)
      r3=dsqrt(3.d0)
!
!     ---- case 1: unit cube C3D8, uniaxial tension along x.
!     n=(1,0,0), so L is the x extent = 1.
!
      call cube(xl,1.d0,1.d0,1.d0)
      nope=8
      lakonl='C3D8    '
      call zero(stre)
      stre(1)=100.d0
      call damcbwidth(stre,xl,nope,lakonl,w,iok)
      write(*,'(a,f12.7,a,f12.7,a,i2)')
     &  ' cube 1x1x1, tension x        L=',w,'  expect=',1.d0,'  iok=',
     &  iok
!
!     ---- case 2: same cube, tension along (1,1,0)/sqrt2.
!     This is the diagonal band Jirasek & Bauer describe: the correct
!     width is sqrt(2) times the aligned one, and a volume estimate
!     cannot see it (it would still say 1).
!
      call zero(stre)
      stre(1)=100.d0
      stre(2)=100.d0
      stre(4)=100.d0
      call damcbwidth(stre,xl,nope,lakonl,w,iok)
      write(*,'(a,f12.7,a,f12.7,a,i2)')
     &  ' cube 1x1x1, tension 45 deg   L=',w,'  expect=',r2,'  iok=',iok
!
!     ---- case 3: elongated brick 1 x 0.2 x 1, tension along x.
!     L=1, while (V)^(1/3)=0.2^(1/3)=0.5848: the elongation error.
!
      call cube(xl,1.d0,0.2d0,1.d0)
      call zero(stre)
      stre(1)=100.d0
      call damcbwidth(stre,xl,nope,lakonl,w,iok)
      legacy=(1.d0*0.2d0*1.d0)**(1.d0/3.d0)
      write(*,'(a,f12.7,a,f12.7,a,f12.7)')
     &  ' brick 1x0.2x1, tension x     L=',w,'  expect=',1.d0,
     &  '  V^(1/3)=',legacy
!
!     ---- case 4: same brick, tension along y (across the thin side).
!     L=0.2.  Same element, same volume, five times the width: this is
!     the direction dependence the volume estimate does not have.
!
      call zero(stre)
      stre(2)=100.d0
      call damcbwidth(stre,xl,nope,lakonl,w,iok)
      write(*,'(a,f12.7,a,f12.7,a,f12.7)')
     &  ' brick 1x0.2x1, tension y     L=',w,'  expect=',0.2d0,
     &  '  V^(1/3)=',legacy
!
!     ---- case 5: C3D4, the corner tetrahedron of a unit cube,
!     tension along x.  Nodes (0,0,0),(1,0,0),(0,1,0),(0,0,1):
!     the x extent is 1, while 6V=1 so the legacy length is also 1.
!     The two agree here, which is why C3D4 on a regular mesh never
!     exposed the defect.
!
      do i=1,20
        xl(1,i)=0.d0
        xl(2,i)=0.d0
        xl(3,i)=0.d0
      enddo
      xl(1,2)=1.d0
      xl(2,3)=1.d0
      xl(3,4)=1.d0
      nope=4
      lakonl='C3D4    '
      call zero(stre)
      stre(1)=100.d0
      call damcbwidth(stre,xl,nope,lakonl,w,iok)
      vol=1.d0
      write(*,'(a,f12.7,a,f12.7,a,f12.7)')
     &  ' tet corner, tension x        L=',w,'  expect=',1.d0,
     &  ' (6V)^(1/3)=',vol**(1.d0/3.d0)
!
!     ---- case 6: hydrostatic tension.  Every eigenvalue is equal, the
!     band normal is not determined by the stress, and the routine must
!     REFUSE (iok=0) rather than invent a direction.
!
      call cube(xl,1.d0,1.d0,1.d0)
      nope=8
      lakonl='C3D8    '
      call zero(stre)
      stre(1)=100.d0
      stre(2)=100.d0
      stre(3)=100.d0
      call damcbwidth(stre,xl,nope,lakonl,w,iok)
      write(*,'(a,i2,a,f12.7,a,f12.7)')
     &  ' hydrostatic cube: iok=',iok,'  L=',w,'  expect=',1.d0
!
!     ---- case 7: hydrostatic on the ELONGATED brick.  No band normal, so
!     the direction-free mean of the three extents: (1+0.2+1)/3 = 0.7333.
!     A finite, shape-dependent width for a family that has no (6V)^(1/3),
!     which is what stops the two damage entry points disagreeing.
!
      call cube(xl,1.d0,0.2d0,1.d0)
      call damcbwidth(stre,xl,nope,lakonl,w,iok)
      write(*,'(a,i2,a,f12.7,a,f12.7)')
     &  ' hydrostatic brick: iok=',iok,' L=',w,'  expect=',
     &  (1.d0+0.2d0+1.d0)/3.d0
!
!     ---- case 8: hydrostatic on a C3D8 is now answered rather than
!     refused, so a hexahedron near hydrostatic no longer leaves charlen
!     at zero.  Degeneracy is reserved for a degenerate ELEMENT:
!
      call cube(xl,0.d0,0.d0,0.d0)
      call damcbwidth(stre,xl,nope,lakonl,w,iok)
      write(*,'(a,i2,a)') ' degenerate element: iok=',iok,
     &  '  expect=0 (the only real data error)'
!
      end
!
      subroutine cube(xl,a,b,c)
      implicit none
      real*8 xl(3,20),a,b,c
      integer i
      do i=1,20
        xl(1,i)=0.d0
        xl(2,i)=0.d0
        xl(3,i)=0.d0
      enddo
      xl(1,2)=a
      xl(1,3)=a
      xl(2,3)=b
      xl(2,4)=b
      xl(3,5)=c
      xl(1,6)=a
      xl(3,6)=c
      xl(1,7)=a
      xl(2,7)=b
      xl(3,7)=c
      xl(2,8)=b
      xl(3,8)=c
      return
      end
!
      subroutine zero(s)
      implicit none
      real*8 s(6)
      integer i
      do i=1,6
        s(i)=0.d0
      enddo
      return
      end
