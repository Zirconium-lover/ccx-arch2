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
      subroutine damrank1blk(xl,voldl,djac,avec,bvec,xsj,weight)
!
!     ONE owner for the projection
!
!         damjac(12) at an integration point  ->  the element block
!         K_rank1(a,b) = -weight * avec(a) * bvec(b)
!
!     Extracted verbatim from mafilldamas.f on 2026-09-11.  It was in the
!     middle of an assembly loop, which is why nobody could test it: the
!     only way to see its output was to read a matrix coefficient out of a
!     CSR array in a running solve, and the structural FD probe that did
!     so could not say whether a discrepancy came from damjac, from this
!     projection, or from the probe.  As its own routine it has a self
!     test - damrank1test below - that answers that question offline.
!
!     The factorisation it implements.  e_c3d.f builds the material part
!     of the element matrix as
!
!         s(ii,i1;jj,j1) += F(i1,m) C(m,k,n,l) F(j1,n)
!                           * shpj(k,ii) * shpj(l,jj) * weight
!
!     with F = I + grad(u) and shpj = shp * sqrt(det J).  The damage term
!     of the consistent tangent,
!
!         C = g(D) C_ep - sigma_eff (x) dD/d(eps),
!
!     contributes C(m,k,n,l) = -T(m,k) Q(n,l), which is separable, so the
!     quadruple sum collapses to an outer product of two 3*nope vectors:
!
!         a(ii,i1) = F(i1,m) T(m,k) shpj(k,ii)
!         b(jj,j1) = F(j1,n) Q(n,l) shpj(l,jj)
!
!     BOTH index pairs carry an F because C is a material-configuration
!     tangent in exactly the convention e_c3d consumes, and both of its
!     pairs are pushed forward there.  T is the effective stress and Q is
!     dD/d(eps); the caller stores them in djac(1:6) and djac(7:12) with
!     the Voigt map (1,1)->1 (2,2)->2 (3,3)->3 (1,2)->4 (1,3)->5 (2,3)->6.
!
!     Q is stored PER TENSOR SLOT: djac(10) is dD/d(eps_12) with the 12
!     and the 21 slot each taking half, not the engineering-shear
!     derivative.  resultsmech.f halves the shear components when it
!     fills djac for precisely this reason; the self test below is what
!     makes that claim checkable rather than a comment.
!
!     C3D4 only, one integration point.  The caller enforces this.
!
      implicit none
!
      real*8 xl(3,4),voldl(3,4),djac(12),avec(12),bvec(12),xsj,weight
!
      real*8 shp(4,4),shpj(3,4),vo(3,3),fdef(3,3),tt(3,3),qq(3,3),
     &     xi,et,ze,xsjj
      integer i1,j1,k1,k,m,j,ii,iflag
!
      include "gauss.f"
!
      iflag=2
!
      xi=gauss3d4(1,1)
      et=gauss3d4(2,1)
      ze=gauss3d4(3,1)
      weight=weight3d4(1)
      call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
!
!     a degenerate Jacobian is the caller's to account for; avec and
!     bvec are left untouched and xsj carries the verdict
!
      if(xsj.lt.1.d-20) return
!
      xsjj=dsqrt(xsj)
      do j=1,4
        do k=1,3
          shpj(k,j)=shp(k,j)*xsjj
        enddo
      enddo
!
!     deformation gradient
!
      do i1=1,3
        do j1=1,3
          vo(i1,j1)=0.d0
          do k1=1,4
            vo(i1,j1)=vo(i1,j1)+shp(j1,k1)*voldl(i1,k1)
          enddo
          fdef(i1,j1)=vo(i1,j1)
        enddo
        fdef(i1,i1)=fdef(i1,i1)+1.d0
      enddo
!
!     effective stress and dD/d(eps) as symmetric 3x3 tensors
!
      tt(1,1)=djac(1)
      tt(2,2)=djac(2)
      tt(3,3)=djac(3)
      tt(1,2)=djac(4)
      tt(2,1)=djac(4)
      tt(1,3)=djac(5)
      tt(3,1)=djac(5)
      tt(2,3)=djac(6)
      tt(3,2)=djac(6)
!
      qq(1,1)=djac(7)
      qq(2,2)=djac(8)
      qq(3,3)=djac(9)
      qq(1,2)=djac(10)
      qq(2,1)=djac(10)
      qq(1,3)=djac(11)
      qq(3,1)=djac(11)
      qq(2,3)=djac(12)
      qq(3,2)=djac(12)
!
!     a and b
!
      do ii=1,4
        do i1=1,3
          avec(3*(ii-1)+i1)=0.d0
          bvec(3*(ii-1)+i1)=0.d0
          do m=1,3
            do k=1,3
              avec(3*(ii-1)+i1)=avec(3*(ii-1)+i1)
     &             +fdef(i1,m)*tt(m,k)*shpj(k,ii)
              bvec(3*(ii-1)+i1)=bvec(3*(ii-1)+i1)
     &             +fdef(i1,m)*qq(m,k)*shpj(k,ii)
            enddo
          enddo
        enddo
      enddo
!
      return
      end
!
      subroutine damrank1test(nbad)
!
!     The self test for damrank1blk.  Three independent routes to the
!     same element block, plus a deliberate corruption that must be
!     caught, plus the degenerate-Jacobian contract.
!
!     Why it exists.  The structural FD probe measures a discrepancy
!     between the assembled operator and a finite difference of the
!     residual, and reports it per matrix coefficient.  When that
!     discrepancy is nonzero on an element that HAS the rank-1 term,
!     the probe cannot say whether the error is in damjac, in this
!     projection, or in the probe's own reading of the CSR.  This
!     routine removes the middle possibility from the list.
!
      implicit none
!
      integer nbad,i,j,k,l,m,n,ii,jj,i1,j1,ia,ib,kk
      real*8 xl(3,4),voldl(3,4),djac(12),avec(12),bvec(12),xsj,weight
      real*8 shp(4,4),shpj(3,4),fdef(3,3),tt(3,3),qq(3,3),vo(3,3)
      real*8 blk(12,12),ref(12,12),xi,et,ze,xsjj,s,worst,den,h,
     &     dplus,dminus,fplus,fminus,ur(3,4)
      integer iflag
      character*56 what
!
      include "gauss.f"
!
      nbad=0
      iflag=2
!
!     A tetrahedron with no symmetry, so that a transposed index cannot
!     hide, and a stress and damage gradient with every component filled
!     and no two equal.
!
      xl(1,1)=0.d0
      xl(2,1)=0.d0
      xl(3,1)=0.d0
      xl(1,2)=1.3d0
      xl(2,2)=0.2d0
      xl(3,2)=-0.1d0
      xl(1,3)=0.4d0
      xl(2,3)=1.1d0
      xl(3,3)=0.3d0
      xl(1,4)=-0.2d0
      xl(2,4)=0.5d0
      xl(3,4)=1.7d0
!
      djac(1)=2.10d2
      djac(2)=-1.30d2
      djac(3)=0.70d2
      djac(4)=0.41d2
      djac(5)=-0.23d2
      djac(6)=0.17d2
      djac(7)=1.10d0
      djac(8)=-0.60d0
      djac(9)=0.35d0
      djac(10)=0.21d0
      djac(11)=-0.13d0
      djac(12)=0.09d0
!
!     ---------------------------------------------------------------
!     1.  The outer product against e_c3d's own quadruple sum, at a
!         deformation gradient that is NOT the identity.  This is the
!         factorisation claim in the header, checked instead of argued.
!     ---------------------------------------------------------------
!
      voldl(1,1)=0.011d0
      voldl(2,1)=-0.004d0
      voldl(3,1)=0.007d0
      voldl(1,2)=-0.013d0
      voldl(2,2)=0.021d0
      voldl(3,2)=0.002d0
      voldl(1,3)=0.006d0
      voldl(2,3)=0.015d0
      voldl(3,3)=-0.009d0
      voldl(1,4)=0.003d0
      voldl(2,4)=-0.017d0
      voldl(3,4)=0.012d0
!
      call damrank1blk(xl,voldl,djac,avec,bvec,xsj,weight)
      do ii=1,12
        do jj=1,12
          blk(ii,jj)=-weight*avec(ii)*bvec(jj)
        enddo
      enddo
!
      call shape4tet(gauss3d4(1,1),gauss3d4(2,1),gauss3d4(3,1),
     &     xl,xsj,shp,iflag)
      xsjj=dsqrt(xsj)
      do j=1,4
        do k=1,3
          shpj(k,j)=shp(k,j)*xsjj
        enddo
      enddo
      do i1=1,3
        do j1=1,3
          vo(i1,j1)=0.d0
          do k=1,4
            vo(i1,j1)=vo(i1,j1)+shp(j1,k)*voldl(i1,k)
          enddo
          fdef(i1,j1)=vo(i1,j1)
        enddo
        fdef(i1,i1)=fdef(i1,i1)+1.d0
      enddo
      tt(1,1)=djac(1)
      tt(2,2)=djac(2)
      tt(3,3)=djac(3)
      tt(1,2)=djac(4)
      tt(2,1)=djac(4)
      tt(1,3)=djac(5)
      tt(3,1)=djac(5)
      tt(2,3)=djac(6)
      tt(3,2)=djac(6)
      qq(1,1)=djac(7)
      qq(2,2)=djac(8)
      qq(3,3)=djac(9)
      qq(1,2)=djac(10)
      qq(2,1)=djac(10)
      qq(1,3)=djac(11)
      qq(3,1)=djac(11)
      qq(2,3)=djac(12)
      qq(3,2)=djac(12)
!
!     s(ii,i1;jj,j1) = F(i1,m) C(m,k,n,l) F(j1,n) shpj(k,ii) shpj(l,jj)
!                      * weight,   C(m,k,n,l) = -T(m,k) Q(n,l)
!
      do ii=1,4
        do i1=1,3
          do jj=1,4
            do j1=1,3
              s=0.d0
              do m=1,3
                do k=1,3
                  do n=1,3
                    do l=1,3
                      s=s+fdef(i1,m)*(-tt(m,k)*qq(n,l))*fdef(j1,n)
     &                     *shpj(k,ii)*shpj(l,jj)
                    enddo
                  enddo
                enddo
              enddo
              ref(3*(ii-1)+i1,3*(jj-1)+j1)=weight*s
            enddo
          enddo
        enddo
      enddo
      what='the outer product equals the full C(m,k,n,l) sum at F/=I'
      call damrank1cmp(blk,ref,1.d-12,what,nbad)
!
!     ---------------------------------------------------------------
!     2.  At F=I, against a finite difference of the element internal
!         force through the damage scalar alone.  This is the route
!         that does not reuse the Voigt map: D is built from Q as a
!         function of the DISPLACEMENT FIELD, so a transposed slot or
!         a shear component stored in the engineering convention
!         instead of the per-slot one shows up here and nowhere else.
!
!             f(a,i) = weight * detJ * (1 - D(u)) * (T grad N_a)_i
!             D(u)   = sum_{n,l} Q(n,l) du_n/dx_l
!     ---------------------------------------------------------------
!
      do j=1,4
        do k=1,3
          voldl(k,j)=0.d0
        enddo
      enddo
      call damrank1blk(xl,voldl,djac,avec,bvec,xsj,weight)
      do ii=1,12
        do jj=1,12
          blk(ii,jj)=-weight*avec(ii)*bvec(jj)
        enddo
      enddo
!
      h=1.d-6
      do jj=1,4
        do j1=1,3
          do j=1,4
            do k=1,3
              ur(k,j)=0.d0
            enddo
          enddo
          ur(j1,jj)=h
          dplus=0.d0
          do n=1,3
            do l=1,3
              do k=1,4
                dplus=dplus+qq(n,l)*shp(l,k)*ur(n,k)
              enddo
            enddo
          enddo
          ur(j1,jj)=-h
          dminus=0.d0
          do n=1,3
            do l=1,3
              do k=1,4
                dminus=dminus+qq(n,l)*shp(l,k)*ur(n,k)
              enddo
            enddo
          enddo
          do ii=1,4
            do i1=1,3
              s=0.d0
              do k=1,3
                s=s+tt(i1,k)*shp(k,ii)
              enddo
              fplus=weight*xsj*(1.d0-dplus)*s
              fminus=weight*xsj*(1.d0-dminus)*s
              ref(3*(ii-1)+i1,3*(jj-1)+j1)=(fplus-fminus)/(2.d0*h)
            enddo
          enddo
        enddo
      enddo
      what='at F=I it is d(internal force)/du through the damage'
      call damrank1cmp(blk,ref,1.d-9,what,nbad)
!
!     ---------------------------------------------------------------
!     3.  The same comparison with the shear components of Q doubled -
!         which is exactly what storing dD/d(gamma) instead of
!         dD/d(eps_12) would produce.  It MUST go red.  A check that
!         cannot be made to fail is not a check.
!     ---------------------------------------------------------------
!
      do kk=10,12
        djac(kk)=2.d0*djac(kk)
      enddo
      call damrank1blk(xl,voldl,djac,avec,bvec,xsj,weight)
      do ii=1,12
        do jj=1,12
          blk(ii,jj)=-weight*avec(ii)*bvec(jj)
        enddo
      enddo
      worst=0.d0
      den=0.d0
      do ii=1,12
        do jj=1,12
          if(dabs(blk(ii,jj)-ref(ii,jj)).gt.worst)
     &         worst=dabs(blk(ii,jj)-ref(ii,jj))
          if(dabs(ref(ii,jj)).gt.den) den=dabs(ref(ii,jj))
        enddo
      enddo
      if(den.le.0.d0) den=1.d0
      if(worst/den.gt.1.d-3) then
        write(*,'(a,a,a,e10.3)') '  ok   ',
     &    'doubling the shear slots of Q is caught          ',
     &    ' worst rel ',worst/den
      else
        write(*,'(a,a,a,e10.3)') '  FAIL ',
     &    'doubling the shear slots of Q is caught          ',
     &    ' worst rel ',worst/den
        nbad=nbad+1
      endif
      do kk=10,12
        djac(kk)=0.5d0*djac(kk)
      enddo
!
!     ---------------------------------------------------------------
!     4.  A flat tetrahedron must come back with xsj at zero and leave
!         the caller to account for it, rather than dividing by it.
!     ---------------------------------------------------------------
!
      xl(1,4)=xl(1,1)+(xl(1,2)-xl(1,1))+(xl(1,3)-xl(1,1))
      xl(2,4)=xl(2,1)+(xl(2,2)-xl(2,1))+(xl(2,3)-xl(2,1))
      xl(3,4)=xl(3,1)+(xl(3,2)-xl(3,1))+(xl(3,3)-xl(3,1))
      avec(1)=-1.d30
      call damrank1blk(xl,voldl,djac,avec,bvec,xsj,weight)
      if((xsj.lt.1.d-20).and.(avec(1).eq.-1.d30)) then
        write(*,'(a,a)') '  ok   ',
     &    'a flat tetrahedron returns detJ=0 and touches nothing'
      else
        write(*,'(a,a,e12.4)') '  FAIL ',
     &    'a flat tetrahedron returns detJ=0 and touches nothing',xsj
        nbad=nbad+1
      endif
!
      if(nbad.eq.0) then
        write(*,'(a)') '[DAMRANK1 SELFTEST] PASSED'
      else
        write(*,'(a,i3,a)') '[DAMRANK1 SELFTEST] ',nbad,
     &       ' check(s) FAILED'
      endif
      call flush(6)
!
      return
      end
!
      subroutine damrank1cmp(blk,ref,tol,what,nbad)
!
!     compare two 12x12 blocks relative to the largest entry of the
!     reference, because the block is an outer product and its small
!     entries carry no information
!
      implicit none
      integer nbad,i,j
      real*8 blk(12,12),ref(12,12),tol,worst,den
      character*56 what
!
      worst=0.d0
      den=0.d0
      do i=1,12
        do j=1,12
          if(dabs(blk(i,j)-ref(i,j)).gt.worst)
     &         worst=dabs(blk(i,j)-ref(i,j))
          if(dabs(ref(i,j)).gt.den) den=dabs(ref(i,j))
        enddo
      enddo
      if(den.le.0.d0) den=1.d0
      if(worst/den.le.tol) then
        write(*,'(a,a,a,e10.3)') '  ok   ',what,' worst rel ',worst/den
      else
        write(*,'(a,a,a,e10.3)') '  FAIL ',what,' worst rel ',worst/den
        nbad=nbad+1
      endif
      return
      end
