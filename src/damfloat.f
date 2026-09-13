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
      subroutine damfloat(ipkon,kon,lakon,ne,ne0,nk,nodeboun,nboun,
     &     ipompc,nodempc,nmpc,batchmax,nfloat,nreach,nfcoh,nfisl)
!
!     Detects material that has become detached from every support and
!     marks it terminal.
!
!     nfcoh / nfisl SPLIT THE REPORTING OF nfloat AND NOTHING ELSE.
!     nfloat stays exactly what it was: the single running total that
!     every batch test below compares against batchmax, and the value
!     returned to the caller.  The two new counters are pure observers
!     written after the fact - no loop is reordered, no budget is split,
!     no comparison reads them.  They exist because [DAMAGE DETACHED]
!     printed one number for two physically different removals - a
!     cohesive facet that lost its bulk on one side, and a bulk island
!     that lost every load path - and the logs could not tell which had
!     happened.  A cohesive facet goes at g=gmin*Kn, an island element
!     goes at whatever it was carrying; conflating them hid that.
!
!     BK4 only checks whether the nodes of the elements just deleted have
!     lost all their elements.  A region that is internally connected but
!     no longer attached to anything - an island cut out by the crack -
!     passes that test and leaves the global operator singular in its
!     rigid-body modes.  The symptom is unmistakable: a Newton correction
!     of the order of the specimen size, in a node nowhere near the crack
!     front, with the residual growing instead of contracting.
!
!     The test is a breadth-first sweep over the surviving elements
!     seeded from every node that is held by a boundary condition or
!     appears in an MPC.  Anything not reached has no load path to
!     ground.  Elements are all-or-nothing here: as soon as one of their
!     nodes is reached the sweep marks the rest, so an unreached element
!     has no reached node at all.
!
!     User elements are traversed as well.  An intact cohesive facet is a
!     genuine load path, and ignoring it would declare a perfectly
!     supported hydride island detached.
!
!     Marking uses exactly the DE1.3 convention, ipkon -> -ipkon-2, so the
!     caller's existing tentative/commit/rollback bookkeeping picks the
!     elements up unchanged.
!
      implicit none
!
      character*8 lakon(*)
!
      integer ipkon(*),kon(*),ne,ne0,nk,nodeboun(*),nboun,
     &     ipompc(*),nodempc(3,*),nmpc,batchmax,nfloat,nreach,
     &     nfcoh,nfisl,
     &     i,j,k,node,nope,indexe,index,nstack,ifree,iel,ndof
!
      integer, allocatable :: mark(:),iponoel(:),inoel(:,:),istack(:)
!
      integer damfaceinit,damfaceon
      character*132 damfaceenv
      save damfaceinit,damfaceon
      data damfaceinit /0/
!
      nfloat=0
      nreach=0
      nfcoh=0
      nfisl=0
      if(nk.le.0) return
!
!     First: cohesive facets that have lost the material on one side.
!
!     A zero-thickness cohesive element exists to tie two pieces of bulk
!     together.  Once every bulk element around one of its node triples
!     has been deleted, the facet connects that triple to nothing, and the
!     triple survives on the residual cohesive stiffness alone - gmin*Kn,
!     which is deliberately negligible.  The result is a nearly singular
!     3x3 block per node and a Newton correction of the order of the
!     specimen size, far from the crack front.
!
!     This is not the same thing as a detached island and the reachability
!     sweep below cannot see it, because the sweep traverses cohesive
!     elements too - correctly so, since an intact facet is a real load
!     path.  It has to be removed first.
!
      call damfloatcoh(ipkon,kon,lakon,ne,nk,batchmax,nfloat)
      nfcoh=nfloat
!
!     CCX_DAMAGE_FLOAT_FACE=1 switches the reachability test below from
!     node connectivity to FACE connectivity.
!
!     Why: the sweep below marks an element reached as soon as ONE of its
!     nodes is reached, so a tetrahedron that touches the structure at a
!     single node counts as supported.  It is not.  Such an element is a
!     free body carrying near-rigid rotational modes, and E-56 measured
!     exactly that: on m12_field at its stall the softest mode of the
!     assembled operator had 1/sigma_min >= 230 and was carried by TWO
!     nodes out of 9555, both belonging to ONE tetrahedron that shared no
!     face with the body and had no surviving cohesive facet.  72 of the
!     73 face-connected components were floating and 70 of them still
!     shared a node, hence invisible to the node-level test.
!
!     Face reachability is strictly stronger - anything reachable by faces
!     is reachable by nodes - so the face sweep replaces the node sweep
!     rather than supplementing it.
!
!     DEFAULT OFF.  damdangle entered this code as a plausible mechanism
!     backed by a real observation, was adopted on one model, and proved
!     net-negative on the full ladder (E-22).  This one is off until it
!     has been accepted on the whole ladder, so every existing baseline
!     stays valid and the A/B is exact.
!
      if(damfaceinit.eq.0) then
        damfaceinit=1
        damfaceon=0
        call getenv('CCX_DAMAGE_FLOAT_FACE',damfaceenv)
        if(damfaceenv(1:1).eq.'1') damfaceon=1
        if(damfaceon.eq.1) then
          write(*,*) '[DAMAGE FLOAT] face connectivity enabled'
        endif
      endif
!
      if(damfaceon.eq.1) then
        call damfloatface(ipkon,kon,lakon,ne,ne0,nk,nodeboun,nboun,
     &       ipompc,nodempc,nmpc,batchmax,nfloat,nreach)
        nfisl=nfloat-nfcoh
        return
      endif
!
!     node -> element lists over the surviving elements
!
      ifree=0
      do i=1,ne
        if(ipkon(i).lt.0) cycle
        call damfloatnope(lakon(i),nope,ndof)
        if(nope.le.0) cycle
        ifree=ifree+nope
      enddo
      if(ifree.le.0) return
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
!     seeds: every constrained node
!
      nstack=0
      do i=1,nboun
        node=nodeboun(i)
        if((node.lt.1).or.(node.gt.nk)) cycle
        if(mark(node).ne.0) cycle
        mark(node)=1
        nstack=nstack+1
        istack(nstack)=node
      enddo
      do i=1,nmpc
        index=ipompc(i)
        do
          if(index.eq.0) exit
          node=nodempc(1,index)
          if((node.ge.1).and.(node.le.nk)) then
            if(mark(node).eq.0) then
              mark(node)=1
              nstack=nstack+1
              istack(nstack)=node
            endif
          endif
          index=nodempc(3,index)
        enddo
      enddo
      if(nstack.eq.0) then
        deallocate(istack,inoel,iponoel,mark)
        return
      endif
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
!     any surviving bulk element with no reached node is detached
!
      do i=1,ne0
        if(nfloat.ge.batchmax) exit
        if(ipkon(i).lt.0) cycle
        if(lakon(i)(1:1).ne.'C') cycle
        call damfloatnope(lakon(i),nope,ndof)
        if(nope.le.0) cycle
        indexe=ipkon(i)
        do j=1,nope
          node=kon(indexe+j)
          if((node.lt.1).or.(node.gt.nk)) cycle
          if(mark(node).ne.0) goto 100
        enddo
        ipkon(i)=-ipkon(i)-2
        nfloat=nfloat+1
 100    continue
      enddo
!
      nfisl=nfloat-nfcoh
      deallocate(istack,inoel,iponoel,mark)
      return
      end
!
!     ------------------------------------------------------------------
!
      subroutine damfloatnope(lakonl,nope,ndof)
      implicit none
      character*8 lakonl
      integer nope,ndof
!
      nope=0
      ndof=3
      if(lakonl(1:1).eq.'F') return
      if(lakonl(1:1).eq.'U') then
!
!       user element: the node count is encoded in position 8, the same
!       way mafillsm.f reads it
!
        ndof=ichar(lakonl(7:7))
        nope=ichar(lakonl(8:8))
        if((nope.lt.1).or.(nope.gt.20)) nope=0
        return
      endif
      if(lakonl(1:1).ne.'C') return
      if(lakonl(4:4).eq.'4') then
        nope=4
      elseif(lakonl(4:5).eq.'10') then
        nope=10
      elseif(lakonl(4:4).eq.'8') then
        nope=8
      elseif(lakonl(4:4).eq.'2') then
        nope=20
      elseif(lakonl(4:5).eq.'15') then
        nope=15
      elseif(lakonl(4:4).eq.'6') then
        nope=6
      endif
      return
      end
!
!     ------------------------------------------------------------------
!     Removes cohesive facets that no longer have bulk material on both
!     sides.  Nodes 1..nope/2 are the minus surface, the rest the plus
!     surface; a facet is kept only while both halves still touch at
!     least one surviving bulk element.
!
      subroutine damfloatcoh(ipkon,kon,lakon,ne,nk,batchmax,nfloat)
      implicit none
!
      character*8 lakon(*)
!
      integer ipkon(*),kon(*),ne,nk,batchmax,nfloat,
     &     i,j,node,nope,ndof,indexe,nhalf,iminus,iplus
!
      integer, allocatable :: hasbulk(:)
!
      if(nk.le.0) return
      allocate(hasbulk(nk))
      do i=1,nk
        hasbulk(i)=0
      enddo
!
      do i=1,ne
        if(ipkon(i).lt.0) cycle
        if(lakon(i)(1:1).ne.'C') cycle
        call damfloatnope(lakon(i),nope,ndof)
        if(nope.le.0) cycle
        indexe=ipkon(i)
        do j=1,nope
          node=kon(indexe+j)
          if((node.ge.1).and.(node.le.nk)) hasbulk(node)=1
        enddo
      enddo
!
      do i=1,ne
        if(nfloat.ge.batchmax) exit
        if(ipkon(i).lt.0) cycle
        if(lakon(i)(1:1).ne.'U') cycle
        call damfloatnope(lakon(i),nope,ndof)
        if(nope.lt.2) cycle
        nhalf=nope/2
        if(2*nhalf.ne.nope) cycle
        indexe=ipkon(i)
        iminus=0
        iplus=0
        do j=1,nhalf
          node=kon(indexe+j)
          if((node.ge.1).and.(node.le.nk)) then
            if(hasbulk(node).ne.0) iminus=1
          endif
        enddo
        do j=nhalf+1,nope
          node=kon(indexe+j)
          if((node.ge.1).and.(node.le.nk)) then
            if(hasbulk(node).ne.0) iplus=1
          endif
        enddo
        if((iminus.eq.0).or.(iplus.eq.0)) then
          ipkon(i)=-ipkon(i)-2
          nfloat=nfloat+1
        endif
      enddo
!
      deallocate(hasbulk)
      return
      end
!
!     ------------------------------------------------------------------
!     Reachability from the supports by FACE connectivity.
!
!     Two elements are neighbours when they share at least three nodes.
!     For a tetrahedron that is a face; for a hexahedron a face has four;
!     for a zero-thickness cohesive element the shared triple is exactly
!     the surface it is glued to.  Sharing an edge (two nodes) or a vertex
!     (one) is NOT a load path in any useful sense - it leaves the piece
!     free to rotate about that edge or point - and that is the whole
!     difference from the node-level sweep.
!
!     LIMITATION: the >=3 rule is correct for LINEAR elements.  On
!     quadratic tetrahedra two elements sharing only an edge also share
!     three nodes (two corners plus the mid-side node), so the rule would
!     over-connect.  Every model this branch is used on is C3D4; a warning
!     is printed if a quadratic volume element is present.
!
      subroutine damfloatface(ipkon,kon,lakon,ne,ne0,nk,nodeboun,nboun,
     &     ipompc,nodempc,nmpc,batchmax,nfloat,nreach)
      implicit none
!
      character*8 lakon(*)
!
      integer ipkon(*),kon(*),ne,ne0,nk,nodeboun(*),nboun,
     &     ipompc(*),nodempc(3,*),nmpc,batchmax,nfloat,nreach,
     &     i,j,k,node,nope,nopej,indexe,indexj,index,nstack,ifree,
     &     iel,jel,ndof,ncom,iwarn
!
      integer, allocatable :: marke(:),iponoel(:),inoel(:,:),istack(:),
     &     nodeflag(:)
!
      nreach=0
      if(nk.le.0) return
!
      ifree=0
      iwarn=0
      do i=1,ne
        if(ipkon(i).lt.0) cycle
        call damfloatnope(lakon(i),nope,ndof)
        if(nope.le.0) cycle
        if((lakon(i)(1:1).eq.'C').and.(nope.gt.8)) iwarn=1
        ifree=ifree+nope
      enddo
      if(ifree.le.0) return
      if(iwarn.eq.1) then
        write(*,*) '[DAMAGE FLOAT] *WARNING: quadratic volume elements'
        write(*,*) '               present; the >=3 shared-node rule'
        write(*,*) '               over-connects on those.'
      endif
!
      allocate(marke(ne))
      allocate(iponoel(nk))
      allocate(inoel(2,ifree))
      allocate(istack(ne))
      allocate(nodeflag(nk))
!
      do i=1,ne
        marke(i)=0
      enddo
      do i=1,nk
        iponoel(i)=0
        nodeflag(i)=0
      enddo
!
      ifree=0
      do i=1,ne
        if(ipkon(i).lt.0) cycle
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
!     seeds: every element touching a constrained node.  A support is a
!     nodal condition, so the element carrying it is grounded whatever
!     its face connectivity.
!
      nstack=0
      do i=1,nboun
        node=nodeboun(i)
        if((node.lt.1).or.(node.gt.nk)) cycle
        index=iponoel(node)
        do
          if(index.eq.0) exit
          iel=inoel(1,index)
          if(marke(iel).eq.0) then
            marke(iel)=1
            nstack=nstack+1
            istack(nstack)=iel
          endif
          index=inoel(2,index)
        enddo
      enddo
      do i=1,nmpc
        index=ipompc(i)
        do
          if(index.eq.0) exit
          node=nodempc(1,index)
          if((node.ge.1).and.(node.le.nk)) then
            k=iponoel(node)
            do
              if(k.eq.0) exit
              iel=inoel(1,k)
              if(marke(iel).eq.0) then
                marke(iel)=1
                nstack=nstack+1
                istack(nstack)=iel
              endif
              k=inoel(2,k)
            enddo
          endif
          index=nodempc(3,index)
        enddo
      enddo
      if(nstack.eq.0) then
        deallocate(nodeflag,istack,inoel,iponoel,marke)
        return
      endif
!
      k=0
      do
        if(k.ge.nstack) exit
        k=k+1
        iel=istack(k)
        call damfloatnope(lakon(iel),nope,ndof)
        indexe=ipkon(iel)
!
!       flag this element's nodes, so a candidate's shared-node count is
!       one pass over the candidate instead of a nope x nope comparison
!
        do j=1,nope
          node=kon(indexe+j)
          if((node.ge.1).and.(node.le.nk)) nodeflag(node)=1
        enddo
        do j=1,nope
          node=kon(indexe+j)
          if((node.lt.1).or.(node.gt.nk)) cycle
          index=iponoel(node)
          do
            if(index.eq.0) exit
            jel=inoel(1,index)
            if(marke(jel).ne.0) then
              index=inoel(2,index)
              cycle
            endif
            call damfloatnope(lakon(jel),nopej,ndof)
            indexj=ipkon(jel)
            ncom=0
            do i=1,nopej
              node=kon(indexj+i)
              if((node.ge.1).and.(node.le.nk)) then
                if(nodeflag(node).ne.0) ncom=ncom+1
              endif
            enddo
            if(ncom.ge.3) then
              marke(jel)=1
              nstack=nstack+1
              istack(nstack)=jel
            endif
            index=inoel(2,index)
          enddo
        enddo
        do j=1,nope
          node=kon(indexe+j)
          if((node.ge.1).and.(node.le.nk)) nodeflag(node)=0
        enddo
      enddo
      nreach=nstack
!
!     any surviving bulk element not reached by faces is detached
!
      do i=1,ne0
        if(nfloat.ge.batchmax) exit
        if(ipkon(i).lt.0) cycle
        if(lakon(i)(1:1).ne.'C') cycle
        call damfloatnope(lakon(i),nope,ndof)
        if(nope.le.0) cycle
        if(marke(i).ne.0) cycle
        ipkon(i)=-ipkon(i)-2
        nfloat=nfloat+1
      enddo
!
      deallocate(nodeflag,istack,inoel,iponoel,marke)
      return
      end
!
!     ------------------------------------------------------------------
!     Flags the nodes of pieces that have come loose by FACE connectivity,
!     so the caller can put a small stiffness on their free modes instead
!     of deleting them.
!
!     Why not delete: E-57 measured that.  Removing floating pieces turns
!     R0, R1, R2 and the UC6 PASS disk from PASS to FAIL, because SP1
!     reaches rc=0 while carrying nine of them, each removal is another
!     topology event needing re-equilibration, and it changes the physics
!     to fix what is a conditioning problem.
!
!     What is actually wrong is narrower.  E-56 and the follow-up located
!     three independent pointers at the m12_field wall - the softest mode
!     of the operator (1/sigma_min = 230 against 3.7 in the elastic phase),
!     the dof whose Newton correction equals the whole increment, and the
!     node carrying the largest residual - and ALL THREE sit on nodes that
!     carry a one-element face-disconnected piece.  Such a piece touches
!     the structure at a node or an edge, so it retains free rotational
!     modes; the operator is singular in those directions and Newton
!     solves for them.
!
!     A node is flagged only when EVERY one of its surviving elements
!     belongs to a floating piece.  A node shared between the body and a
!     dangling element - node 6177 at the m12 wall - is NOT flagged: it is
!     properly supported, and it is the dangling element's own free nodes
!     that need holding.  This keeps the perturbation off any dof that
!     carries load.
!
      subroutine damfloatstab(ipkon,kon,lakon,ne,nk,nodeboun,nboun,
     &     ipompc,nodempc,nmpc,nodestab,nstab)
      implicit none
!
      character*8 lakon(*)
!
      integer ipkon(*),kon(*),ne,nk,nodeboun(*),nboun,
     &     ipompc(*),nodempc(3,*),nmpc,nodestab(*),nstab,
     &     i,j,k,node,nope,nopej,indexe,indexj,index,nstack,ifree,
     &     iel,jel,ndof,ncom
!
      integer, allocatable :: marke(:),iponoel(:),inoel(:,:),istack(:),
     &     nodeflag(:)
!
      nstab=0
      if(nk.le.0) return
      do i=1,nk
        nodestab(i)=0
      enddo
!
      ifree=0
      do i=1,ne
        if(ipkon(i).lt.0) cycle
        call damfloatnope(lakon(i),nope,ndof)
        if(nope.le.0) cycle
        ifree=ifree+nope
      enddo
      if(ifree.le.0) return
!
      allocate(marke(ne))
      allocate(iponoel(nk))
      allocate(inoel(2,ifree))
      allocate(istack(ne))
      allocate(nodeflag(nk))
!
      do i=1,ne
        marke(i)=0
      enddo
      do i=1,nk
        iponoel(i)=0
        nodeflag(i)=0
      enddo
!
      ifree=0
      do i=1,ne
        if(ipkon(i).lt.0) cycle
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
      do i=1,nboun
        node=nodeboun(i)
        if((node.lt.1).or.(node.gt.nk)) cycle
        index=iponoel(node)
        do
          if(index.eq.0) exit
          iel=inoel(1,index)
          if(marke(iel).eq.0) then
            marke(iel)=1
            nstack=nstack+1
            istack(nstack)=iel
          endif
          index=inoel(2,index)
        enddo
      enddo
      do i=1,nmpc
        index=ipompc(i)
        do
          if(index.eq.0) exit
          node=nodempc(1,index)
          if((node.ge.1).and.(node.le.nk)) then
            k=iponoel(node)
            do
              if(k.eq.0) exit
              iel=inoel(1,k)
              if(marke(iel).eq.0) then
                marke(iel)=1
                nstack=nstack+1
                istack(nstack)=iel
              endif
              k=inoel(2,k)
            enddo
          endif
          index=nodempc(3,index)
        enddo
      enddo
      if(nstack.eq.0) then
        deallocate(nodeflag,istack,inoel,iponoel,marke)
        return
      endif
!
      k=0
      do
        if(k.ge.nstack) exit
        k=k+1
        iel=istack(k)
        call damfloatnope(lakon(iel),nope,ndof)
        indexe=ipkon(iel)
        do j=1,nope
          node=kon(indexe+j)
          if((node.ge.1).and.(node.le.nk)) nodeflag(node)=1
        enddo
        do j=1,nope
          node=kon(indexe+j)
          if((node.lt.1).or.(node.gt.nk)) cycle
          index=iponoel(node)
          do
            if(index.eq.0) exit
            jel=inoel(1,index)
            if(marke(jel).ne.0) then
              index=inoel(2,index)
              cycle
            endif
            call damfloatnope(lakon(jel),nopej,ndof)
            indexj=ipkon(jel)
            ncom=0
            do i=1,nopej
              node=kon(indexj+i)
              if((node.ge.1).and.(node.le.nk)) then
                if(nodeflag(node).ne.0) ncom=ncom+1
              endif
            enddo
            if(ncom.ge.3) then
              marke(jel)=1
              nstack=nstack+1
              istack(nstack)=jel
            endif
            index=inoel(2,index)
          enddo
        enddo
        do j=1,nope
          node=kon(indexe+j)
          if((node.ge.1).and.(node.le.nk)) nodeflag(node)=0
        enddo
      enddo
!
!     nodeflag is reused: 1 = the node has at least one REACHED element
!
      do i=1,nk
        nodeflag(i)=0
      enddo
      do i=1,ne
        if(ipkon(i).lt.0) cycle
        if(marke(i).eq.0) cycle
        call damfloatnope(lakon(i),nope,ndof)
        if(nope.le.0) cycle
        indexe=ipkon(i)
        do j=1,nope
          node=kon(indexe+j)
          if((node.ge.1).and.(node.le.nk)) nodeflag(node)=1
        enddo
      enddo
!
!     a node is flagged when it has surviving elements but not one of
!     them was reached
!
      do i=1,ne
        if(ipkon(i).lt.0) cycle
        if(marke(i).ne.0) cycle
        call damfloatnope(lakon(i),nope,ndof)
        if(nope.le.0) cycle
        indexe=ipkon(i)
        do j=1,nope
          node=kon(indexe+j)
          if((node.lt.1).or.(node.gt.nk)) cycle
          if(nodeflag(node).ne.0) cycle
          if(nodestab(node).eq.0) then
            nodestab(node)=1
            nstab=nstab+1
          endif
        enddo
      enddo
!
      deallocate(nodeflag,istack,inoel,iponoel,marke)
      return
      end
!
!     Remove a cohesive facet whose cohesive damage is COMPLETE.
!
!     WHY THIS EXISTS.  The bulk phases create a displacement discontinuity
!     when they fail - the element is deleted and the crack is a real free
!     surface.  The interface never does: `g = max(gmin, 1-dvisc)` floors the
!     facet stiffness at gmin, so a fully debonded facet keeps conducting
!     gmin*Kn for ever and the two sides are never allowed to separate.
!     `mh_facets.py` measures the consequence directly - on the radial twin
!     45.3% of backed facets are at Dmin >= 0.999 and ZERO have been removed
!     (E-113) - and E-111 item 1 named it as the remaining gap after the
!     internal length (E-114) and the penalty stiffness (E-113) were both
!     tested and excluded.  E-119/E-120 sharpened it: the route that will not
!     resolve is the DEBONDING one, and this is the only part of that route
!     that cannot create a crack.
!
!     `xstate(4,.,.)` is the per-integration-point flag `resultsmech_uc6.f`
!     already writes when dback >= 1-1e-12.  `CCX_FRACTURE_DEADFACET` reads it
!     to keep such a facet out of the TERMINATION test; this routine acts on
!     the same flag in the EQUATIONS.
!
!     UC6 carries exactly three integration points (cohesive_uc6.f, the
!     three-point triangular rule).  mi(1) is the allocated stride and looping
!     to it would test uninitialised slots - the defect E-84 found in the
!     nonlocal average and the C-side dead-facet scan already guards against.
!
!     KNOWN APPROXIMATION, stated because nothing in the output reports it:
!     deleting the facet removes its COMPRESSIVE penalty too, so the two faces
!     may afterwards interpenetrate.  The cohesive law keeps a full penalty in
!     compression (V5, E-41) and this discards it.  That is the same trade
!     every deleted bulk element already makes, and `damfloatcoh` already
!     removes facets on a different criterion, but it is a physical
!     approximation and not a free one.
!
!     Deletion uses the branch's own convention, ipkon -> -ipkon-2, so the
!     facet travels the unchanged transactional path and is restored by the
!     existing rollback (E-64) exactly as damfloatcoh's removals are.
!
      subroutine damfaildead(ipkon,lakon,ne,xstate,nstate_,mi,
     &     batchmax,ndead)
      implicit none
!
      character*8 lakon(*)
!
      integer ipkon(*),ne,nstate_,mi(*),batchmax,ndead,i,j,nip,iok
!
      real*8 xstate(nstate_,mi(1),*)
!
      if(nstate_.lt.4) return
      nip=3
      if(nip.gt.mi(1)) nip=mi(1)
      if(nip.lt.1) return
!
      do i=1,ne
        if(ndead.ge.batchmax) exit
        if(ipkon(i).lt.0) cycle
        if(lakon(i)(1:1).ne.'U') cycle
        iok=1
        do j=1,nip
          if(xstate(4,j,i).lt.0.5d0) then
            iok=0
            exit
          endif
        enddo
        if(iok.eq.1) then
          ipkon(i)=-ipkon(i)-2
          ndead=ndead+1
        endif
      enddo
!
      return
      end
