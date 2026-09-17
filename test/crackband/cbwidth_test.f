!     Unit test of the CB1 projected crack-band width.
!     Each case has a hand-computable answer, derived in the comment
!     beside it rather than fitted to the output, and each is CHECKED
!     rather than merely printed: a disagreement sets the exit code.
!
!     It used to only print measured next to expected and exit 0 no
!     matter what, which made it unusable as a gate check and made a
!     claim of mine about it false.  The verdict column and the exit
!     code below are what that claim should have been describing.
      program cbtest
      implicit none
      real*8 xl(3,20),stre(6),w,r2,vol,legacy
      integer nope,iok,i,nfail
      character*8 lakonl
      nfail=0
      r2=dsqrt(2.d0)
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
      call chk('cube 1x1x1, tension x       ',w,1.d0,nfail)
      call chki('  and it is determined      ',iok,1,nfail)
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
      call chk('cube 1x1x1, tension 45 deg  ',w,r2,nfail)
!
!     ---- case 3: elongated brick 1 x 0.2 x 1, tension along x.
!     The band runs across the LONG direction, so L=1, while the
!     volume estimate says (0.2)^(1/3)=0.5848 for both this case and
!     the next one - it cannot tell them apart at all.
!
      call cube(xl,1.d0,0.2d0,1.d0)
      call zero(stre)
      stre(1)=100.d0
      call damcbwidth(stre,xl,nope,lakonl,w,iok)
      legacy=(1.d0*0.2d0*1.d0)**(1.d0/3.d0)
      call chk('brick 1x0.2x1, tension x    ',w,1.d0,nfail)
!
!     ---- case 4: same brick, tension along y (across the thin side).
!     L=0.2 now.  Same element, same volume, width five times smaller:
!     this is the pair that shows the volume estimate is not a width.
!
      call zero(stre)
      stre(2)=100.d0
      call damcbwidth(stre,xl,nope,lakonl,w,iok)
      call chk('brick 1x0.2x1, tension y    ',w,0.2d0,nfail)
      write(*,'(a,f12.7,a)')
     &  '   both of the above have V^(1/3)=',legacy,
     &  ', which is neither'
!
!     ---- case 5: C3D4, the corner tetrahedron of a unit cube,
!     tension along x.  Nodes (0,0,0),(1,0,0),(0,1,0),(0,0,1):
!     the x extent is 1, while 6V=1 so the legacy length is also 1.
!     The two agree here, which is why C3D4 on a regular mesh never
!     exposed the defect.
!
      call tet(xl)
      nope=4
      lakonl='C3D4    '
      call zero(stre)
      stre(1)=100.d0
      call damcbwidth(stre,xl,nope,lakonl,w,iok)
      vol=1.d0
      call chk('tet corner, tension x       ',w,1.d0,nfail)
      write(*,'(a,f12.7,a)')
     &  '   here (6V)^(1/3)=',vol**(1.d0/3.d0),' agrees, by accident'
!
!     ---- case 6: C3D10 on the SAME four corners, but with a midside
!     node pushed far outside the element, to x=5.
!
!     Two things are checked at once.  The width must be unchanged at
!     1, which it can only be if damcbcorner really does exclude the
!     midside nodes - a straight-edge midside would hide the bug,
!     because it lies inside the corner hull and changes no extent.
!     And this is the case that triggers damcbwarnonce, so the
!     higher-order warning is exercised rather than assumed.
!
!     Jirasek and Bauer 2012 section 6.5 is why that warning exists:
!     on quadratic elements the band localizes into one layer of GAUSS
!     POINTS, so the whole-element projection is too wide, releases
!     less than G_f, and makes the response too brittle.
!
      call tet(xl)
      xl(1,5)=5.d0
      nope=10
      lakonl='C3D10   '
      call damcbwidth(stre,xl,nope,lakonl,w,iok)
      call chk('C3D10, midside node at x=5  ',w,1.d0,nfail)
!
!     ---- case 7: hydrostatic tension.  Every eigenvalue is equal and
!     the band normal is not determined by the stress, so the routine
!     must not invent a direction: it falls back on the direction-free
!     mean width, which for a unit cube is 1.
!
      call cube(xl,1.d0,1.d0,1.d0)
      nope=8
      lakonl='C3D8    '
      call zero(stre)
      stre(1)=100.d0
      stre(2)=100.d0
      stre(3)=100.d0
      call damcbwidth(stre,xl,nope,lakonl,w,iok)
      call chk('hydrostatic cube            ',w,1.d0,nfail)
!
!     ---- case 8: hydrostatic on the ELONGATED brick.  No band normal,
!     so the direction-free mean of the three extents: (1+0.2+1)/3.
!     A finite, shape-dependent width for a family that has no
!     (6V)^(1/3), which is what stops the two damage entry points
!     disagreeing.
!
      call cube(xl,1.d0,0.2d0,1.d0)
      call damcbwidth(stre,xl,nope,lakonl,w,iok)
      call chk('hydrostatic brick           ',w,
     &  (1.d0+0.2d0+1.d0)/3.d0,nfail)
!
!     ---- case 9: a degenerate ELEMENT is the one case that is a data
!     error, and the only one that may refuse.
!
      call cube(xl,0.d0,0.d0,0.d0)
      call damcbwidth(stre,xl,nope,lakonl,w,iok)
      call chki('degenerate element refuses  ',iok,0,nfail)
!
      write(*,*)
      if(nfail.gt.0) then
        write(*,'(a,i3,a)') ' cbwidth_test: ',nfail,' check(s) FAILED'
        call exit(1)
      endif
      write(*,*) 'cbwidth_test: all checks ok'
      end
!
      subroutine chk(lab,got,want,nfail)
!
!     One checked value.  The tolerance is absolute and loose enough
!     for a square root and tight enough that a wrong extent, which is
!     the failure this test is for, can never pass it.
!
      implicit none
      character*(*) lab
      real*8 got,want,tol
      integer nfail
      character*4 verd
      tol=1.d-9
      if(dabs(got-want).le.tol) then
        verd=' ok '
      else
        verd='FAIL'
        nfail=nfail+1
      endif
      write(*,'(1x,a,a,f12.7,a,f12.7,a,a)')
     &  lab,' L=',got,'  expect=',want,'   ',verd
      return
      end
!
      subroutine chki(lab,got,want,nfail)
      implicit none
      character*(*) lab
      integer got,want,nfail
      character*4 verd
      if(got.eq.want) then
        verd=' ok '
      else
        verd='FAIL'
        nfail=nfail+1
      endif
      write(*,'(1x,a,a,i2,a,i2,a,a)')
     &  lab,' iok=',got,'  expect=',want,'             ',verd
      return
      end
!
      subroutine tet(xl)
!
!     The corner tetrahedron of a unit cube, with all midside slots
!     zeroed so a caller can place one deliberately.
!
      implicit none
      real*8 xl(3,20)
      integer i
      do i=1,20
        xl(1,i)=0.d0
        xl(2,i)=0.d0
        xl(3,i)=0.d0
      enddo
      xl(1,2)=1.d0
      xl(2,3)=1.d0
      xl(3,4)=1.d0
      return
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
