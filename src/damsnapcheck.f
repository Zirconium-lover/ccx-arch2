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
      subroutine damsnapcheck(co,kon,ipkon,lakon,ne0,ielmat,mi,
     &     ndmcon,dmcon,ndmat_,ntmat_,elcon,ncmat_,plicon,nplicon,
     &     npmat_,ratmax,ielmax,nbad)
!
!     Crack-band well-posedness check for the DE1/DM2 displacement
!     evolution law.
!
!     After initiation the law degrades an element over a plastic
!     displacement u_f, so the softening modulus seen by the element is
!     roughly sigma_y * L / u_f.  Compared with the elastic unloading
!     modulus E this gives the dimensionless indicator
!
!         r = L * sigma_y / (u_f * E)
!
!     r > 1 means the element's own softening branch is steeper than the
!     elastic release it can supply, i.e. local snap-back.
!
!     This is an INDICATOR, not a proof of anything about the global
!     problem.  A structure containing such elements may still reach
!     equilibrium by redistributing into stiffer neighbours, and whether
!     it does depends on the surrounding mesh, the boundary conditions
!     and the load path - none of which this scalar sees.  What r > 1
!     does reliably predict is that a consistent tangent will expose an
!     indefinite element operator there, so convergence may be poor
!     where the stock positive-definite secant tangent appeared to cope.
!
!     The routine only measures and reports; it never changes the model
!     and never stops the run.
!
      implicit none
!
      character*8 lakon(*)
!
      integer mi(*)
!
      integer kon(*),ipkon(*),ne0,ielmat(mi(3),*),ndmcon(2,*),
     &     ndmat_,ntmat_,ncmat_,nplicon(0:ntmat_,*),npmat_,
     &     ielmax,nbad,i,imat,itype,nconst,indexe,n1,n2,n3,n4
!
      real*8 co(3,*),dmcon(0:ndmat_,ntmat_,*),elcon(0:ncmat_,ntmat_,*),
     &     plicon(0:2*npmat_,ntmat_,*),ratmax,
     &     ax,ay,az,bx,by,bz,cx,cy,cz,det6v,charlen,ufail,ee,sy,ratio
!
      ratmax=0.d0
      ielmax=0
      nbad=0
!
      do i=1,ne0
!
        if(ipkon(i).lt.0) cycle
        if(lakon(i)(1:1).ne.'C') cycle
        if(lakon(i)(4:4).ne.'4') cycle
        if(lakon(i)(7:8).eq.'LC') cycle
!
        imat=ielmat(1,i)
        if(imat.lt.1) cycle
        nconst=ndmcon(1,imat)
        if(ndmcon(2,imat).eq.0) cycle
        itype=int(dmcon(1,1,imat))
        if(itype.eq.1) then
          if(nconst.ne.4) cycle
          ufail=dmcon(4,1,imat)
        elseif(itype.eq.3) then
          if(nconst.lt.7) cycle
          ufail=dmcon(3,1,imat)
        else
          cycle
        endif
        if(ufail.le.0.d0) cycle
!
        ee=elcon(1,1,imat)
        if(ee.le.0.d0) cycle
        if(nplicon(1,imat).lt.1) cycle
        sy=plicon(1,1,imat)
        if(sy.le.0.d0) cycle
!
!       characteristic length of the reference tetrahedron,
!       the same L=(6 V0)^(1/3) that calcdamage.f uses
!
        indexe=ipkon(i)
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
     &       +az*(bx*cy-by*cx))
        if(det6v.le.1.d-30) cycle
        charlen=det6v**(1.d0/3.d0)
!
        ratio=charlen*sy/(ufail*ee)
        if(ratio.gt.ratmax) then
          ratmax=ratio
          ielmax=i
        endif
        if(ratio.gt.1.d0) nbad=nbad+1
!
      enddo
!
      return
      end
