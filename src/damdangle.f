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
      subroutine damdangle(ipkon,kon,lakon,ne,nk,batchmax,ndanglemax,
     &     ihazin,ndangle,nweak,nbare)
!
!     Removes material left hanging on a single element after a deletion
!     batch, and marks it terminal.
!
!     BK4 asks whether a node has lost ALL of its elements, and damfloat
!     asks whether a region has lost every load path to a support.  A node
!     that goes from twenty-two elements down to one passes both tests: it
!     still has an element, and that element is still attached to ground
!     through its other three nodes.  Nothing in the chain notices, and the
!     next solve throws the node most of a specimen length.
!
!     This was measured, not guessed.  In the notched single-hydride model
!     the run stalls at lambda=0.3579 in a loop: increments 401, 403, 405
!     and 407 all collect the same three elements - 7007, 7038, 7101 - all
!     roll back, and all fail the same way.  Node 2003 belongs to all three
!     and drops from four surviving elements to one.  It is the node the
!     solver then moves by 0.90 while the load step is 2.9e-6.
!
!     Cutting the step cannot help, because the cascade is driven by the
!     deletion and not by the load: one element crosses the threshold, and
!     removing it takes the other two to D=1 exactly.  Halving the step
!     took the batch from six elements to four to three and no further.
!     That is also why the dissipation target made no difference - 0.30 and
!     0.10 both stop at the same lambda to four figures - and why a
!     hundredfold change in the residual stiffness floor changed nothing.
!
!     A node held by one tetrahedron is a sliver of material that no longer
!     takes part in carrying load, so deleting it is the same physical
!     statement damfloat already makes about detached islands, applied to a
!     case its reachability sweep cannot see.
!
!     The closure is small and it terminates: on the stuck batch it adds
!     exactly one element in one sweep, and applied to the whole committed
!     deletion history of that model it adds none at all.  It therefore
!     cannot change any state that already equilibrated.
!
!     Support is counted over every surviving element type, cohesive ones
!     included, for the reason damfloat gives: an intact facet is a real
!     load path, and ignoring it would condemn a properly held node.
!
!     Only nodes that LOST support are eligible.  A mesh can legitimately
!     contain a vertex on one element - this one has three of them from the
!     start - and those must be left alone.
!
!     Marking uses the DE1.3 convention, ipkon -> -ipkon-2, so the caller's
!     tentative/commit/rollback bookkeeping picks the elements up unchanged.
!
      implicit none
!
      character*8 lakon(*)
!
      integer ipkon(*),kon(*),ne,nk,batchmax,ndanglemax,ndangle,nweak,
     &     ihazin(*),nbare,i,j,node,nope,indexe,nsweep,nnew,ndof
!
      integer, allocatable :: nlive(:),norig(:),ihaz(:),nbulk(:)
!
      ndangle=0
      nweak=0
      nbare=0
      if(nk.le.0) return
!
      allocate(nlive(nk))
      allocate(norig(nk))
      allocate(ihaz(nk))
      allocate(nbulk(nk))
!
!     Repeat until closure: removing the last element of a dangling node
!     can in principle strand another one.  Bounded so a pathological mesh
!     cannot spin here; in practice one sweep is enough.
!
      do nsweep=1,10
!
        do i=1,nk
          nlive(i)=0
          norig(i)=0
          ihaz(i)=0
          nbulk(i)=0
        enddo
!
!       nlive counts elements still assembled, norig counts those the mesh
!       started with.  A deleted element carries ipkon=-indexe-2, so it is
!       distinguishable from one that never existed (ipkon=-1).
!
        do i=1,ne
          if(ipkon(i).eq.-1) cycle
          call damfloatnope(lakon(i),nope,ndof)
          if(nope.le.0) cycle
          if(ipkon(i).ge.0) then
            indexe=ipkon(i)
          else
            indexe=-ipkon(i)-2
          endif
          if(indexe.lt.0) cycle
          do j=1,nope
            node=kon(indexe+j)
            if((node.lt.1).or.(node.gt.nk)) cycle
            norig(node)=norig(node)+1
            if(ipkon(i).ge.0) then
              nlive(node)=nlive(node)+1
              if(lakon(i)(1:1).eq.'C') nbulk(node)=nbulk(node)+1
            endif
          enddo
        enddo
!
!       a node is a hazard when it still has an element but has lost its
!       way down to the threshold from a support that used to be larger
!
        nnew=0
        do i=1,nk
          if(nlive(i).lt.1) cycle
!
!         Two independent ways to be stranded.  The topological
!         one - support collapsed to ndanglemax from more - is what
!         caught node 2003.  It cannot see node 158, which keeps
!         plenty of elements but whose assembled diagonal has
!         fallen to 1.1e-4 of its intact value; that node stops
!         contracting and its Newton correction equals the whole
!         displacement increment, so the correction criterion can
!         never be met however far the step is cut.  The caller
!         supplies that judgement in ihazin because only it has the
!         assembled operator.
!
          if(ihazin(i).eq.1) then
            ihaz(i)=1
            nnew=nnew+1
            cycle
          endif
          if(ndanglemax.lt.1) cycle
          if(nlive(i).gt.ndanglemax) cycle
          if(norig(i).le.ndanglemax) cycle
          ihaz(i)=1
          nnew=nnew+1
        enddo
!
!       Count, but do not act on, nodes that have lost every bulk
!       element and hang on cohesive facets alone.  A debonded
!       facet carries g*kn with g at the residual floor, so such a
!       node is very nearly free - and it is invisible to all three
!       existing tests: BK4 wants zero elements of ANY type,
!       damfloatcoh works on facet triples rather than single
!       nodes, and the support count above deliberately treats a
!       cohesive element as support, which is right for an intact
!       facet and wrong for a debonded one.
!
        if(nsweep.eq.1) then
          do i=1,nk
            if((nlive(i).ge.1).and.(nbulk(i).eq.0).and.
     &           (norig(i).gt.nlive(i))) nbare=nbare+1
          enddo
        endif
        if(nnew.eq.0) exit
        nweak=nweak+nnew
!
!       mark what is left at those nodes
!
        nnew=0
        do i=1,ne
          if(ipkon(i).lt.0) cycle
          if(ndangle.ge.batchmax) exit
          call damfloatnope(lakon(i),nope,ndof)
          if(nope.le.0) cycle
          indexe=ipkon(i)
          do j=1,nope
            node=kon(indexe+j)
            if((node.lt.1).or.(node.gt.nk)) cycle
            if(ihaz(node).eq.0) cycle
            ipkon(i)=-ipkon(i)-2
            ndangle=ndangle+1
            nnew=nnew+1
            exit
          enddo
        enddo
        if(nnew.eq.0) exit
        if(ndangle.ge.batchmax) exit
      enddo
!
      deallocate(nbulk)
      deallocate(ihaz)
      deallocate(norig)
      deallocate(nlive)
!
      return
      end
