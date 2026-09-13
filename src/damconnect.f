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
      subroutine damconnect(ipkon,kon,lakon,ne,nk,nodesa,nna,
     &     nodesb,nnb,iconn,nreach,ifacdead)
!
!     Load-path test between two node sets over the surviving elements.
!
!     A crack is complete when no chain of still-existing elements links
!     the two boundaries any more:
!
!         tension     Face_X0       <-> Face_XL
!         Lame        InnerPressure <-> OuterRadius
!
!     iconn = 1  a path still exists
!     iconn = 0  the sets are disconnected, i.e. the crack is through
!     nreach     nodes reachable from set A, for diagnostics
!
!     Deleted elements carry ipkon<0 and are skipped, so the test sees
!     exactly the load-bearing topology the assembly sees.  Elements that
!     are merely degraded still conduct: a nonzero stiffness is a load
!     path, and calling it broken would be a modelling decision rather
!     than a topological fact.  Cohesive facets conduct for the same
!     reason - damfloat removes the ones that have lost their bulk.
!
      implicit none
!
      character*8 lakon(*)
!
      integer ipkon(*),kon(*),ne,nk,nodesa(*),nna,nodesb(*),nnb,
     &     iconn,nreach,ndof,ifacdead(*),
     &     i,j,k,node,nope,indexe,index,nstack,ifree,iel
!
      integer, allocatable :: mark(:),iponoel(:),inoel(:,:),istack(:)
!
      iconn=1
      nreach=0
      if(nk.le.0) return
!
      ifree=0
      do i=1,ne
        if(ipkon(i).lt.0) cycle
        if(ifacdead(i).ne.0) cycle
        call damfloatnope(lakon(i),nope,ndof)
        if(nope.le.0) cycle
        ifree=ifree+nope
      enddo
      if(ifree.le.0) then
        iconn=0
        return
      endif
!
      allocate(mark(nk))
      allocate(iponoel(nk))
      allocate(inoel(2,ifree))
      allocate(istack(nk))
!
      do i=1,nk
        mark(i)=0
        iponoel(i)=0
      enddo
!
      ifree=0
      do i=1,ne
        if(ipkon(i).lt.0) cycle
        if(ifacdead(i).ne.0) cycle
        call damfloatnope(lakon(i),nope,ndof)
        if(nope.le.0) cycle
        indexe=ipkon(i)
        do j=1,nope
          node=kon(indexe+j)
          if((node.lt.1).or.(node.gt.nk)) cycle
          ifree=ifree+1
          inoel(1,ifree)=i
          inoel(2,ifree)=iponoel(node)
          iponoel(node)=ifree
        enddo
      enddo
!
      nstack=0
      do i=1,nna
        node=nodesa(i)
        if((node.lt.1).or.(node.gt.nk)) cycle
        if(mark(node).ne.0) cycle
        if(iponoel(node).eq.0) cycle
        mark(node)=1
        nstack=nstack+1
        istack(nstack)=node
      enddo
!
      k=0
      do
        if(k.ge.nstack) exit
        k=k+1
        node=istack(k)
        index=iponoel(node)
        do
          if(index.eq.0) exit
          iel=inoel(1,index)
          call damfloatnope(lakon(iel),nope,ndof)
          indexe=ipkon(iel)
          do j=1,nope
            i=kon(indexe+j)
            if((i.lt.1).or.(i.gt.nk)) cycle
            if(mark(i).ne.0) cycle
            mark(i)=1
            nstack=nstack+1
            istack(nstack)=i
          enddo
          index=inoel(2,index)
        enddo
      enddo
      nreach=nstack
!
      iconn=0
      do i=1,nnb
        node=nodesb(i)
        if((node.lt.1).or.(node.gt.nk)) cycle
        if(mark(node).ne.0) then
          iconn=1
          exit
        endif
      enddo
!
      deallocate(istack,inoel,iponoel,mark)
      return
      end
!
!     ------------------------------------------------------------------
!
      subroutine damconnectface(ipkon,kon,lakon,ne,nk,nodesa,nna,
     &     nodesb,nnb,iconn,nreach,ifacdead)
!
!     The same load-path test as damconnect, but two elements conduct
!     only when they share a whole FACE, not merely a node.
!
!     WHY.  On run_m12_spc3 the specimen finished in two pieces of 12500
!     and 13300 elements that share no node with each other at all; they
!     were joined through ONE tetrahedron touching each of them at a
!     single vertex.  A point contact in a tetrahedral mesh transmits
!     nothing - three zero-energy rotations sit at that vertex - but the
!     node-level sweep reads it as a load path, so [FRACTURE COMPLETE]
!     never fired and the run spent a further 3.2% of the load history,
!     and about a sixth of its wall time, driving a severed specimen to
!     zero force (E-75).
!
!     This is NOT damfloatface (E-57, E-72).  That one DELETED material
!     it judged unreachable and cost R0, R1, R2 and the UC6 disk their
!     PASS.  Termination deletes nothing: it changes no equation, no
!     stiffness and no solution, only the moment the run is allowed to
!     stop.  Selected by CCX_FRACTURE_LINK=FACE; the default stays NODE
!     so every stored baseline keeps its meaning.
!
!     A zero-thickness cohesive facet shares its whole node triple with
!     the bulk on each side, so the >=3 rule joins facet to bulk without
!     a special case.  Two facets meeting along an edge share only two
!     nodes and do not join each other, which is right: each of them
!     reaches the other through the bulk it is glued to.
!
!     LIMITATION, and it is deliberately conservative.  The face rule is
!     written for linear tetrahedra and for the 6-node cohesive element.
!     Any other element type falls back to node conduction and prints a
!     warning once, because a termination test that stops a run too early
!     destroys the result, while one that stops too late merely wastes
!     time - the failure modes are not symmetric.
!
      implicit none
!
      character*8 lakon(*)
!
      integer ipkon(*),kon(*),ne,nk,nodesa(*),nna,nodesb(*),nnb,
     &     iconn,nreach,ndof,ifacdead(*),
     &     i,j,k,m,node,nope,nopj,indexe,indexj,index,nstack,ifree,
     &     iel,jel,nface,nshare,iwarn,ia,ib,ic
      integer nodface(3,4)
!
      integer, allocatable :: emark(:),iponoel(:),inoel(:,:),estack(:),
     &     nmark(:)
!
      iconn=1
      nreach=0
      iwarn=0
      if((nk.le.0).or.(ne.le.0)) return
!
      ifree=0
      do i=1,ne
        if(ipkon(i).lt.0) cycle
        if(ifacdead(i).ne.0) cycle
        call damfloatnope(lakon(i),nope,ndof)
        if(nope.le.0) cycle
        ifree=ifree+nope
      enddo
      if(ifree.le.0) then
        iconn=0
        return
      endif
!
      allocate(emark(ne))
      allocate(nmark(nk))
      allocate(iponoel(nk))
      allocate(inoel(2,ifree))
      allocate(estack(ne))
!
      do i=1,ne
        emark(i)=0
      enddo
      do i=1,nk
        iponoel(i)=0
        nmark(i)=0
      enddo
!
      ifree=0
      do i=1,ne
        if(ipkon(i).lt.0) cycle
        if(ifacdead(i).ne.0) cycle
        call damfloatnope(lakon(i),nope,ndof)
        if(nope.le.0) cycle
        indexe=ipkon(i)
        do j=1,nope
          node=kon(indexe+j)
          if((node.lt.1).or.(node.gt.nk)) cycle
          ifree=ifree+1
          inoel(1,ifree)=i
          inoel(2,ifree)=iponoel(node)
          iponoel(node)=ifree
        enddo
      enddo
!
!     seed: every live element that carries a node of set A
!
      nstack=0
      do i=1,nna
        node=nodesa(i)
        if((node.lt.1).or.(node.gt.nk)) cycle
        index=iponoel(node)
        do
          if(index.eq.0) exit
          iel=inoel(1,index)
          if(emark(iel).eq.0) then
            emark(iel)=1
            nstack=nstack+1
            estack(nstack)=iel
          endif
          index=inoel(2,index)
        enddo
      enddo
!
      k=0
      do
        if(k.ge.nstack) exit
        k=k+1
        iel=estack(k)
        call damfloatnope(lakon(iel),nope,ndof)
        indexe=ipkon(iel)
!
!       the faces of this element, as node triples
!
        nface=0
        if(nope.eq.4) then
          do j=1,4
            nface=nface+1
            m=0
            do i=1,4
              if(i.eq.j) cycle
              m=m+1
              nodface(m,nface)=kon(indexe+i)
            enddo
          enddo
        elseif((nope.eq.6).and.(lakon(iel)(1:1).eq.'U')) then
          nface=2
          do i=1,3
            nodface(i,1)=kon(indexe+i)
            nodface(i,2)=kon(indexe+3+i)
          enddo
        else
!
!         unknown type: fall back to node conduction so the run is never
!         stopped early on an element this routine does not understand
!
          if(iwarn.eq.0) then
            iwarn=1
            write(*,*) '*WARNING in damconnectface: element type ',
     &           lakon(iel)(1:8),' is not covered by the face rule;'
            write(*,*) '         it conducts through single nodes, '//
     &           'which can only DELAY termination.'
          endif
          do j=1,nope
            node=kon(indexe+j)
            if((node.lt.1).or.(node.gt.nk)) cycle
            index=iponoel(node)
            do
              if(index.eq.0) exit
              jel=inoel(1,index)
              if(emark(jel).eq.0) then
                emark(jel)=1
                nstack=nstack+1
                estack(nstack)=jel
              endif
              index=inoel(2,index)
            enddo
          enddo
          cycle
        endif
!
!       a neighbour across face (ia,ib,ic) must appear in the element
!       list of ia and must also contain ib and ic.  Walking one node's
!       list and testing the other two costs far less than comparing
!       every pair of elements that share any node.
!
        do j=1,nface
          ia=nodface(1,j)
          ib=nodface(2,j)
          ic=nodface(3,j)
          if((ia.lt.1).or.(ia.gt.nk)) cycle
          index=iponoel(ia)
          do
            if(index.eq.0) exit
            jel=inoel(1,index)
            index=inoel(2,index)
            if(jel.eq.iel) cycle
            if(emark(jel).ne.0) cycle
            call damfloatnope(lakon(jel),nopj,ndof)
            if(nopj.le.0) cycle
            indexj=ipkon(jel)
            nshare=0
            do m=1,nopj
              node=kon(indexj+m)
              if((node.eq.ib).or.(node.eq.ic)) nshare=nshare+1
            enddo
            if(nshare.lt.2) cycle
            emark(jel)=1
            nstack=nstack+1
            estack(nstack)=jel
          enddo
        enddo
      enddo
!
!     nodes carried by the reached elements, for the same diagnostic the
!     node rule reports
!
      do i=1,ne
        if(emark(i).eq.0) cycle
        call damfloatnope(lakon(i),nope,ndof)
        if(nope.le.0) cycle
        indexe=ipkon(i)
        do j=1,nope
          node=kon(indexe+j)
          if((node.ge.1).and.(node.le.nk)) nmark(node)=1
        enddo
      enddo
      nreach=0
      do i=1,nk
        if(nmark(i).ne.0) nreach=nreach+1
      enddo
!
      iconn=0
      do i=1,nnb
        node=nodesb(i)
        if((node.lt.1).or.(node.gt.nk)) cycle
        if(nmark(node).ne.0) then
          iconn=1
          exit
        endif
      enddo
!
      deallocate(estack,inoel,iponoel,nmark,emark)
      return
      end
