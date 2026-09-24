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
!     ==================================================================
!     CB1: directional crack-band width for the DE1/DM2 softening law.
!
!     WHY THIS EXISTS.  The DE1 evolution law spreads the fracture energy
!     over a band whose width is the characteristic length L:
!
!         D = D_base + L * dEps_p / u_f
!
!     so L is what makes the DISSIPATION per unit crack area independent
!     of element size.  What was here is L=(6V)^(1/3), the cube root of
!     six tetrahedral volumes, for C3D4 and nothing else.
!
!     That is the volume-based estimate, and it is the one the literature
!     rejects outright.  Jirasek and Bauer, "Numerical aspects of the
!     crack band approach", Comput. Struct. 110-111 (2012) 60-78, section
!     5.3.1 and the conclusions in section 7: criteria based "exclusively
!     on the element area or volume can lead to large errors and thus
!     cannot be recommended", the error being "comparable to a
!     misprediction of the fracture energy by 50 % or even more".  The
!     error is not confined to distorted meshes: a band running along the
!     diagonals of a regular mesh is wider than an aligned one by a
!     factor of sqrt(2) (their section 5.3.2), and a volume estimate
!     cannot see the difference because it does not know the direction.
!
!     Crack band theory asks for the width measured ACROSS the band
!     (Bazant and Oh, "Crack band theory for fracture of concrete",
!     Mater. Struct. 16 (1983) 155-177), so the width has to depend on
!     the band normal.  The method implemented here is the one Jirasek and
!     Bauer recommend in section 7: project the element onto the major
!     principal direction, taking that direction at the ELEMENT CENTRE
!     rather than at each integration point separately, which in their
!     notation is Pm:
!
!         L = max_I (x_I . n) - min_I (x_I . n)
!
!     over the CORNER nodes I.  Two of their findings are respected by
!     omission: the direction is never taken per integration point,
!     because that "may result into excessively large or even infinite
!     estimates", and no orientation factor is applied, because "there is
!     no need for the additional orientation factor".
!
!     A projection is defined for any element shape, so the C3D4
!     restriction disappears as a by-product rather than as a separate
!     feature: C3D8, C3D10, C3D20, C3D6 and C3D15 get a width by the same
!     formula.  Their conclusions do warn that "higher-order elements are
!     not suitable for crack band simulations", so the quadratic families
!     are accepted with a warning rather than silently.
!
!     WHY THE DIRECTION IS FROZEN, AND ON WHAT.  n is computed from the
!     COMMITTED stress sti at the start of a physical increment, never
!     from the trial state inside Newton, and it stops being refreshed
!     once initiation is committed (dambase >= xlimit).  Three reasons,
!     all raised by Agent 1 in review and all load-bearing:
!
!       (a) L multiplies dEps_p in the law.  A direction that followed the
!           trial state would make L a function of the Newton iterate.
!       (b) L would then enter the consistent tangent through dL/dEps,
!           a term the rank-1 correction in resultsmech does not have.
!       (c) DE1 is transactional: every quantity in the law is rebuilt
!           from the immutable baseline dambase so that a Newton retry or
!           an A3 rollback cannot leave a partially updated state.  A
!           length derived from the trial state would be the only
!           quantity in the law that is not.
!
!     Deriving n from committed data alone answers all three at once: it
!     is constant across the Newton iterations of an increment AND across
!     rollback retries of that increment, because its inputs are.  This
!     also matches the intent already recorded for the legacy length -
!     "the coordinates in co are the reference mesh coordinates, so the
!     regularization length does not change during NLGEOM iterations".
!
!     Freezing at initiation is the physical statement that the band
!     orientation is fixed when the crack forms and does not then follow
!     the field; a width that drifted during softening would not dissipate
!     G_f, which is the whole point of the band.
!
!     THE MODULE IS PRIVATE TO THIS FILE, for the reason damnonlocal.f
!     states for its own: sharing through accessor subroutines instead of
!     a `use' across files keeps a parallel build from breaking on a
!     module dependency.  The cache is written only from calcdamagebase,
!     where the stress is committed, and only read from
!     damageupdatepoint, which runs inside the stress update.
!
!     OFF BY DEFAULT.  CCX_DAMAGE_CHARLEN is unset or 0 reproduces
!     (6V)^(1/3) bit for bit, so the regression gate does not move; 1
!     selects the projection.
!     ==================================================================
!
      module damcbmod
      implicit none
!
!     cbmode: -1 not yet read, 0 legacy (6V)^(1/3), 1 directional Pm.
!     cbfroz(ip,element) is the frozen width, 0 meaning "not set yet".
!
      integer :: cbmode=-1,cbne=0,cbnip=0,cbwarn=0,cbmeanw=0
      integer :: cbnlw=-1,cbnlwarn=0
      real*8, allocatable :: cbfroz(:,:)
!
!     W2 ENERGY evolution.  cbevol(imat) is the evolution kind read from
!     the material card: 1 = DISPLACEMENT, the historical u_f, 2 = ENERGY,
!     where the same slot carries G_f instead.  cbuf(ip,element) is the
!     effective u_f derived from G_f, frozen at initiation beside the width
!     and for the same reason.
!
      integer, allocatable :: cbevol(:)
      real*8, allocatable :: cbuf(:,:)
      integer :: cbnmat=0,cbune=0,cbunip=0
      end module damcbmod
!
!     ------------------------------------------------------------------
!
      subroutine damcbevolset(imat,kind)
!
!     Record the evolution kind of a material.  Called from the reader.
!
!     WHY NOT dmcon.  The obvious place is one more constant on the card,
!     and that would silently switch progressive damage OFF:
!     damage_progressive_material in nonlingeo.c tests the constant COUNT
!     exactly - nconst==4 for Rice-Tracey, and a parity test on (nconst-3)
!     for the tabulated locus - so an extra constant fails the test, the
!     deck still runs, and it runs on the legacy hard-deletion path
!     instead.  dm2failurestrain derives its point count from the same
!     total and would break too.  Keeping the count fixed and carrying the
!     kind here leaves every one of those readers correct.
!
      use damcbmod
      implicit none
      integer imat,kind,n,i
      integer, allocatable :: tmp(:)
!
      if(imat.lt.1) return
      if(.not.allocated(cbevol)) then
        cbnmat=max(64,imat)
        allocate(cbevol(cbnmat))
        do i=1,cbnmat
          cbevol(i)=1
        enddo
      elseif(imat.gt.cbnmat) then
        n=max(2*cbnmat,imat)
        allocate(tmp(n))
        do i=1,n
          tmp(i)=1
        enddo
        do i=1,cbnmat
          tmp(i)=cbevol(i)
        enddo
        deallocate(cbevol)
        allocate(cbevol(n))
        do i=1,n
          cbevol(i)=tmp(i)
        enddo
        deallocate(tmp)
        cbnmat=n
      endif
      cbevol(imat)=kind
      return
      end
!
      subroutine damcbevolget(imat,kind)
      use damcbmod
      implicit none
      integer imat,kind
      kind=1
      if(.not.allocated(cbevol)) return
      if((imat.lt.1).or.(imat.gt.cbnmat)) return
      kind=cbevol(imat)
      return
      end
!
      subroutine damcbufail(imat,itype,dmcon,ndmat_,ntmat_,iel,iint,
     &     ufail,iok)
!
!     THE ONLY way the damage path should obtain u_f.
!
!     WHY IT EXISTS.  Under EVOLUTION=ENERGY the card slot that
!     EVOLUTION=DISPLACEMENT reads as u_f carries G_f instead, and which of
!     the two it is cannot be stored in dmcon: the constant count is tested
!     exactly in nonlingeo.c, so one more constant switches progressive
!     damage off silently.  The kind therefore lives outside dmcon, and
!     Agent 1's objection to that is correct and worth restating - a reader
!     that takes dmcon(4) or dmcon(3) and does not ask for the kind will
!     take an energy for a displacement, without a word.  That is the same
!     shape of defect as imode=2, as cbmode, and as the constant count
!     itself: a value living somewhere other than the structure that
!     describes it.
!
!     The remedy is not to encode the kind somewhere clever.  It is to stop
!     offering the raw slot: value and meaning leave this routine together,
!     and a caller cannot obtain one without the other.  Encoding the
!     combination in the type code instead - itype=4 for DUCTILE+ENERGY -
!     was the alternative, and it was declined because ENERGY is orthogonal
!     to the criterion, so the type space would grow with the product of
!     the options rather than their sum, and a RICETRACEY+ENERGY code would
!     already be the next one needed.
!
!     iok=0 means the effective u_f is not available yet, which happens
!     only under ENERGY before any equivalent stress has been committed for
!     this point, and means the caller must not evolve damage.
!
      implicit none
      integer imat,itype,ndmat_,ntmat_,iel,iint,iok,kind
      real*8 dmcon(0:ndmat_,ntmat_,*),ufail,cbufv
!
      iok=0
      ufail=0.d0
!
!     The raw slot.  Its MEANING is decided by damcbevolget, never here.
!
      if(itype.eq.1) then
        ufail=dmcon(4,1,imat)
      elseif(itype.eq.3) then
        ufail=dmcon(3,1,imat)
      else
        return
      endif
      if(ufail.le.0.d0) return
!
      call damcbevolget(imat,kind)
      if(kind.eq.2) then
        call damcbufget(iel,iint,cbufv,iok)
        if(iok.eq.0) return
        ufail=cbufv
      endif
      iok=1
      return
      end
!
      subroutine damcbufget(iel,iint,val,iok)
!
!     The effective u_f frozen at initiation.  iok=0 means "not available"
!     and the caller must not evolve damage this increment.
!
      use damcbmod
      implicit none
      integer iel,iint,iok
      real*8 val
      iok=0
      val=0.d0
      if(.not.allocated(cbuf)) return
      if((iel.lt.1).or.(iel.gt.cbune)) return
      if((iint.lt.1).or.(iint.gt.cbunip)) return
      if(cbuf(iint,iel).le.0.d0) return
      val=cbuf(iint,iel)
      iok=1
      return
      end
!
      subroutine damcbufset(iel,iint,val)
!
!     Keep the SMALLEST effective u_f the point has offered, which is the
!     LARGEST equivalent stress it has carried, because u_f = 2 G_f / s0
!     and G_f is a material constant.
!
!     WHY NOT THE LATEST.  This cache is refreshed on every increment for
!     which initiation has not yet committed, and taking the latest looks
!     like "the value at initiation" but is not: by the last increment
!     before the threshold is crossed, the neighbouring points of the same
!     band have often initiated already and are softening, and they unload
!     this one.  The stress sampled there is on the way DOWN.
!
!     MEASURED on the equivalence deck, element 19, point 1, whose card
!     says sigma_0 = 400.1:
!
!       dambase=0.931604   vM=400.0982   u_f=0.0200001
!       dambase=0.981067   vM=367.4603   u_f=0.0217765   <- was frozen
!
!     One increment before the crossing the point carries exactly the
!     flow stress the card describes, to four decimals.  The increment
!     that used to win carries 8 per cent less, and that error goes
!     straight into u_f and from there into the dissipated energy.
!
!     The peak is the right rule and not merely the better one: for a
!     hardening J2 material under monotonic loading the largest
!     equivalent stress before initiation IS the flow stress at
!     initiation, which is the sigma_0 that G_f = sigma_0 u_f / 2 is
!     written in terms of (Bazant and Oh 1983 for the relation itself).
!
!     LIMITATION, stated rather than hidden.  That identity needs the
!     loading to be monotonic and proportional.  What sigma_0 has to be
!     is the flow stress at the moment the initiation criterion is met,
!     and the peak equals it only because a hardening material at yield
!     carries its largest equivalent stress then.  Under NON-proportional
!     loading the largest value the point ever carried can come from a
!     different stress state - which matters here, because RICETRACEY is
!     triaxiality dependent and a point can meet its criterion in a state
!     it did not peak in.
!
!     The refinement, if a case ever needs it, is to take the equivalent
!     stress at the last increment in which the point actually FLOWED,
!     dpeq = xstate(1,jj,i)-xstateini(1,jj,i) > 0, rather than the largest
!     it ever carried.  That is what "the flow stress at initiation"
!     literally means, it coincides with the peak under monotonic
!     proportional loading, and it is strictly better outside it.
!
!     Not done, and the reason is not that it is hard: dpeq is not yet
!     computed where this cache is written (it is formed around line 1518,
!     well below), so taking it would mean reordering the routine.  That
!     is a change with no measurement behind it - no deck here loads
!     non-proportionally, so both rules would give the same numbers and
!     the reorder could not be shown to be an improvement or shown to be
!     safe.  It is written down instead, with the expression to use.
!
      use damcbmod
      implicit none
      integer iel,iint
      real*8 val
      if(.not.allocated(cbuf)) return
      if((iel.lt.1).or.(iel.gt.cbune)) return
      if((iint.lt.1).or.(iint.gt.cbunip)) return
      if(val.le.0.d0) return
      if(cbuf(iint,iel).gt.0.d0) then
        if(val.ge.cbuf(iint,iel)) return
      endif
      cbuf(iint,iel)=val
      return
      end
!
!     ------------------------------------------------------------------
!
      subroutine damcbread()
!
!     One-time environment read, and NOTHING else.
!
!     This is deliberately separate from damcbinit.  The cache has exactly
!     one writer, calcdamagebase, and that invariant is worth keeping - but
!     folding the environment read into the writer left the MODE unreadable
!     on the path that consults it FIRST.  damageupdatepoint runs inside the
!     stress update, on the first Newton iteration of the first increment,
!     before nonlingeo has ever reached calcdamagebase; there cbmode was
!     still -1, damcbmodeget mapped that to 0, and the guard below it
!     rejected every non-tetrahedron no matter what the switch said.  So
!     CCX_DAMAGE_CHARLEN=1 did not lift the C3D4 restriction at all.
!
!     Measured on a C3D8 deck: rc=201 "C3D4 elements only" with the switch
!     unset AND set, on the first increment, in both the uniaxial and the
!     hydrostatic variant.  Found in review by Agent 1 and reproduced here
!     (forum 2026-09-17 01:20).  Reading the environment from a routine
!     both paths call is the fix; the single-writer invariant is untouched.
!
      use damcbmod
      implicit none
      character*132 carg
!
      if(cbmode.ge.0) return
      cbmode=0
      call getenv('CCX_DAMAGE_CHARLEN',carg)
      if(carg(1:1).ne.' ') then
        if(carg(1:1).eq.'1') then
          cbmode=1
        elseif(carg(1:4).eq.'PROJ') then
          cbmode=1
        endif
      endif
      if(cbmode.eq.1) then
        write(*,*)
        write(*,*) '*INFO in calcdamage: CCX_DAMAGE_CHARLEN=1,'
        write(*,*) '      the crack-band width is the element'
        write(*,*) '      projection onto the major principal'
        write(*,*) '      direction taken at the element centre'
        write(*,*) '      (Jirasek and Bauer 2012, section 7),'
        write(*,*) '      frozen at damage initiation.'
        write(*,*)
      endif
      return
      end
!
      subroutine damcbinit(ne0,nip)
!
!     Size the frozen-width cache.  Called from calcdamagebase, which
!     remains its ONLY writer.
!
      use damcbmod
      implicit none
      integer ne0,nip,i,j
!
      call damcbread()
      if(cbmode.ne.1) return
!
      if(allocated(cbfroz)) then
        if((cbne.ge.ne0).and.(cbnip.ge.nip)) return
        deallocate(cbfroz)
      endif
      cbne=ne0
      cbnip=nip
      allocate(cbfroz(cbnip,cbne))
      do i=1,cbne
        do j=1,cbnip
          cbfroz(j,i)=0.d0
        enddo
      enddo
      return
      end
!
      subroutine damcbufinit(ne0,nip)
!
!     Size the frozen effective-u_f cache.  Separate from damcbinit because
!     ENERGY evolution is independent of the crack-band switch: a deck may
!     ask for G_f while keeping the legacy (6V)^(1/3) width.
!
      use damcbmod
      implicit none
      integer ne0,nip,i,j
!
      if(allocated(cbuf)) then
        if((cbune.ge.ne0).and.(cbunip.ge.nip)) return
        deallocate(cbuf)
      endif
      cbune=ne0
      cbunip=nip
      allocate(cbuf(cbunip,cbune))
      do i=1,cbune
        do j=1,cbunip
          cbuf(j,i)=0.d0
        enddo
      enddo
      return
      end
!
      subroutine damcbmodeget(mode)
!
!     Readable from either damage entry point, including the one that runs
!     before calcdamagebase has been called for the first time.
!
      use damcbmod
      implicit none
      integer mode
      call damcbread()
      mode=cbmode
      if(mode.lt.0) mode=0
      return
      end
!
      subroutine damcbget(iel,iint,val,iok)
!
!     Read the frozen width.  iok=0 means "not available" and the caller
!     must keep its own value - the same contract damnonlocalval uses.
!
      use damcbmod
      implicit none
      integer iel,iint,iok
      real*8 val
      iok=0
      val=0.d0
      if(cbmode.ne.1) return
      if(.not.allocated(cbfroz)) return
      if((iel.lt.1).or.(iel.gt.cbne)) return
      if((iint.lt.1).or.(iint.gt.cbnip)) return
      if(cbfroz(iint,iel).le.0.d0) return
      val=cbfroz(iint,iel)
      iok=1
      return
      end
!
      subroutine damcbset(iel,iint,val)
      use damcbmod
      implicit none
      integer iel,iint
      real*8 val
      if(cbmode.ne.1) return
      if(.not.allocated(cbfroz)) return
      if((iel.lt.1).or.(iel.gt.cbne)) return
      if((iint.lt.1).or.(iint.gt.cbnip)) return
      if(val.le.0.d0) return
      cbfroz(iint,iel)=val
      return
      end
!
      subroutine damcbnlwread()
!
!     CB-NL: read CCX_DAMAGE_NLWIDTH once.  OFF BY DEFAULT, so every deck
!     that ran before this switch existed keeps its answer to the bit.
!
      use damcbmod
      implicit none
      character*132 carg
!
      if(cbnlw.ge.0) return
      cbnlw=0
      call getenv('CCX_DAMAGE_NLWIDTH',carg)
      if(carg(1:1).eq.'1') cbnlw=1
      return
      end
!
      subroutine damcbnlwidth(iel,charlen)
!
!     WHEN THE DAMAGE DRIVER IS A NONLOCAL AVERAGE, THE WIDTH IN THE
!     SOFTENING LAW MUST BE THE WIDTH OF THE BAND THAT FORMS.
!
!     D advances as charlen*dEps_p/u_f, so an element reaches D=1 once the
!     plastic displacement ACROSS it reaches u_f, and it therefore
!     dissipates G_f = sigma_0 u_f / 2 per unit area whatever its size.
!     That is the property Bazant and Oh 1983 build the crack band for,
!     and it is exact for one reason only: the band is assumed to be ONE
!     element wide.  Charge it to every element of a band that is n layers
!     wide and the model dissipates n*G_f for one crack.
!
!     An internal length abolishes that assumption on purpose - it exists
!     to fix the band's width in the material rather than in the mesh - so
!     the two regularisations cannot both keep their own length.  Measured
!     on this tree (test/crackband/run_nlwidth_scaling.sh, forum
!     2026-09-17): at h=0.25 the band width approaches 2*ell from ABOVE as
!     ell/h grows - w/(2*ell) of 1.65, 1.26, 1.08 at ell = 0.25, 0.5, 1.0 -
!     and W_post follows it, growing by 2.66 between ell=0.25 and ell=1.0
!     while the pre-peak work agrees to 0.03 per cent.  The dissipation is
!     proportional to the number of ELEMENT LAYERS, which is ell/h, and so
!     is a material constant only when the band is one element wide.
!
!     An earlier version of this comment said the width follows 2*ell to
!     within +8/-12 per cent.  That came from measuring the width as a count
!     of integration points over D>0.5, which was retracted: the same runs
!     give three different verdicts for the three thresholds de1stats
!     writes.  The numbers above are the threshold-free measure,
!     sum(D*V)/A.  The retraction was grepped through test/crackband and not
!     through src/, which is how it survived here - after a retraction the
!     tree is what needs walking, not the directory the work happened in.
!
!     Jirasek and Bauer 2012, section 5, state the requirement directly:
!     the width entering the softening law must be the width of the band
!     that actually forms.
!
!     TWO GUARDS, both of which make this a no-op rather than a guess:
!
!       iok=0 from damnlelleff - the mesh does not resolve the averaging
!       length, so the band the averaging would form is not represented
!       and 2*ell is not the width of anything.  The criterion lives in
!       damnonlocal.f and is NOT repeated here: repeating it would let the
!       two copies disagree.
!
!       2*ell below the element's own width - a band cannot be narrower
!       than the one element that carries it, so the element size remains
!       the floor.  dmax1 rather than a replacement for that reason.
!
      use damcbmod
      implicit none
      integer iel,iokv
      real*8 charlen,ellv,wnl
!
      call damcbnlwread()
      if(cbnlw.ne.1) return
      if(charlen.le.0.d0) return
      call damnlelleff(iel,ellv,iokv)
      if(iokv.ne.1) return
      if(ellv.le.0.d0) return
      wnl=2.d0*ellv
      charlen=dmax1(charlen,wnl)
      if(cbnlwarn.eq.0) then
        cbnlwarn=1
        write(*,*)
        write(*,*) '*INFO in calcdamage: CCX_DAMAGE_NLWIDTH=1, the'
        write(*,*) '      crack-band width in the softening law is'
        write(*,*) '      the band width the nonlocal averaging'
        write(*,*) '      forms, 2*ell, where the mesh resolves it,'
        write(*,*) '      and the element width elsewhere.  Without'
        write(*,*) '      this the fracture energy of the combined'
        write(*,*) '      model scales with ell/h.'
        write(*,*)
      endif
      return
      end
!
      subroutine damcbwarnonce(lakonl)
!
!     Jirasek and Bauer 2012, sections 6.5 and 7: on quadratic elements
!     the band localizes into one layer of GAUSS POINTS rather than one
!     layer of elements, so the element is no longer the band.  The
!     projection then returns a width that is too LARGE, which releases
!     less than G_f and makes the response too brittle - their Fig. 31a
!     shows the whole-element curves falling furthest below the
!     reference, and the halved-width curves sitting closest to it.
!
!     Their section 7 does not regard this as a correctable detail: it
!     concludes that higher-order elements are not suitable for crack
!     band simulations at all, because the corrected factors work only
!     when the band happens to align with the mesh lines.  We warn once
!     rather than refuse, because refusing would lock out C3D10, which
!     is what most tetrahedral meshers produce.
!
      use damcbmod
      implicit none
      character*8 lakonl
      if(cbwarn.ne.0) return
      cbwarn=1
      write(*,*)
      write(*,*) '*WARNING in calcdamage: the crack-band projection'
      write(*,*) '         is being used on a higher-order element ('
     &     //lakonl(1:6)//').'
      write(*,*) '         Jirasek and Bauer 2012, section 6.5: on'
      write(*,*) '         quadratic elements the band localizes into'
      write(*,*) '         one layer of GAUSS POINTS, not one layer'
      write(*,*) '         of elements, so projecting the WHOLE'
      write(*,*) '         element overestimates the band width.'
      write(*,*) '         Less than G_f is then released and the'
      write(*,*) '         response comes out too BRITTLE, not too'
      write(*,*) '         ductile (their Fig. 31a).  In 2D they'
      write(*,*) '         report a true width of h/2 for 2x2'
      write(*,*) '         integration and 13h/18 for 3x3; no such'
      write(*,*) '         factor is established for these elements,'
      write(*,*) '         and their section 7 concludes higher-order'
      write(*,*) '         elements are not suitable for crack-band'
      write(*,*) '         work at all.  Prefer C3D4, C3D6 or C3D8.'
      write(*,*)
      return
      end
!
!     ------------------------------------------------------------------
!
      subroutine damcbcorner(nope,ncor)
!
!     Number of CORNER nodes, which are the first ncor local nodes in the
!     CalculiX ordering for every volume family.  Only corners take part
!     in the projection: midside nodes of a straight edge add nothing to
!     the extent, and on a curved edge they would make the width depend
!     on the curvature rather than on the band.
!
      implicit none
      integer nope,ncor
      if((nope.eq.4).or.(nope.eq.10)) then
        ncor=4
      elseif((nope.eq.8).or.(nope.eq.11).or.(nope.eq.20)) then
        ncor=8
      elseif((nope.eq.6).or.(nope.eq.15)) then
        ncor=6
      else
        ncor=0
      endif
      return
      end
!
      subroutine damcbevec(c,lam,v,iok)
!
!     Unit eigenvector of the symmetric 3x3 matrix c for the eigenvalue
!     lam.  The null space of c-lam*I is spanned by the cross product of
!     any two independent rows; all three pairs are tried and the largest
!     is kept, which is the numerically stable choice.
!
!     iok=0 signals a repeated largest eigenvalue: every pair is then
!     degenerate, the band normal is genuinely not determined by the
!     stress state, and the caller falls back to the volume estimate
!     rather than inventing a direction.
!
      implicit none
      integer iok,i,j,k1,k2
      real*8 c(3,3),lam,v(3),a(3,3),w(3),wn,best,scal
!
      iok=0
      v(1)=0.d0
      v(2)=0.d0
      v(3)=0.d0
!
      scal=0.d0
      do i=1,3
        do j=1,3
          a(i,j)=c(i,j)
          if(dabs(c(i,j)).gt.scal) scal=dabs(c(i,j))
        enddo
      enddo
      if(dabs(lam).gt.scal) scal=dabs(lam)
      if(scal.le.1.d-30) return
      do i=1,3
        a(i,i)=a(i,i)-lam
      enddo
!
      best=0.d0
      do i=1,3
        k1=i
        k2=i+1
        if(k2.gt.3) k2=1
        w(1)=a(k1,2)*a(k2,3)-a(k1,3)*a(k2,2)
        w(2)=a(k1,3)*a(k2,1)-a(k1,1)*a(k2,3)
        w(3)=a(k1,1)*a(k2,2)-a(k1,2)*a(k2,1)
        wn=dsqrt(w(1)*w(1)+w(2)*w(2)+w(3)*w(3))
        if(wn.gt.best) then
          best=wn
          v(1)=w(1)
          v(2)=w(2)
          v(3)=w(3)
        endif
      enddo
!
!     The cross product of two rows of a matrix scaled by scal is of
!     order scal**2, so that is what the test is measured against.
!
      if(best.le.1.d-8*scal*scal) return
      v(1)=v(1)/best
      v(2)=v(2)/best
      v(3)=v(3)/best
      iok=1
      return
      end
!
      subroutine damcbmean(nope,xl,width,iok)
!
!     DIRECTION-FREE width: the mean of the element extents along the three
!     global axes, over the corner nodes.
!
!     WHY THIS EXISTS.  damcbwidth refuses when the largest principal stress
!     is repeated, because then the band normal is genuinely not determined
!     by the stress state.  For C3D4 the caller could fall back on
!     (6V)^(1/3), but that formula was only ever written for tetrahedra, so
!     lifting the C3D4 restriction removed a safety net that no other family
!     ever had: a hexahedron whose stress went near-hydrostatic would leave
!     charlen at zero, and the two damage entry points then disagreed - one
!     stopped the run with "invalid DE1 data", the other returned silently
!     and left the point un-updated.  Raised in review by Agent 1 (forum
!     2026-09-17 00:05), and the disagreement is what this removes.
!
!     It is not a refinement of the projection and does not pretend to be.
!     When there is no preferred direction there is no width to project
!     onto one, and a mean over directions is the honest substitute.  For a
!     cube it returns the edge exactly, so the aligned case is unaffected;
!     for a 1 x 0.2 x 1 brick it returns 0.7333, between the two extents.
!
!     LIMITATION, stated rather than hidden: a mean over the GLOBAL axes is
!     not invariant to how the element is oriented in space.  There is no
!     cheap invariant substitute - averaging over the element's own
!     principal axes is not invariant either, because for a cube those axes
!     are arbitrary and a rotated triple gives larger extents.  The true
!     invariant is the mean caliper width over all directions, which is not
!     worth its cost for a state the model cannot orient anyway.  This path
!     is a fallback, it says so in the log, and it is never the default.
!
      implicit none
      integer nope,iok,ncor,i,k
      real*8 xl(3,20),width,pmin,pmax,tot
!
      iok=0
      width=0.d0
      call damcbcorner(nope,ncor)
      if(ncor.eq.0) return
!
      tot=0.d0
      do k=1,3
        pmin=1.d30
        pmax=-1.d30
        do i=1,ncor
          if(xl(k,i).lt.pmin) pmin=xl(k,i)
          if(xl(k,i).gt.pmax) pmax=xl(k,i)
        enddo
        tot=tot+(pmax-pmin)
      enddo
      width=tot/3.d0
      if(width.le.1.d-30) return
      iok=1
      return
      end
!
      subroutine damcbmeanonce()
      use damcbmod
      implicit none
      if(cbmeanw.ne.0) return
      cbmeanw=1
      write(*,*)
      write(*,*) '*WARNING in calcdamage: at least one integration'
      write(*,*) '         point reached equal largest principal'
      write(*,*) '         stresses, where the crack-band normal is'
      write(*,*) '         not determined by the stress state.  The'
      write(*,*) '         direction-free mean element width is used'
      write(*,*) '         there instead of the projection.'
      write(*,*)
      return
      end
!
      subroutine damcbwidth(stre,xl,nope,lakonl,width,iok)
!
!     Projected crack-band width.  stre is the COMMITTED stress at the
!     integration point in CalculiX order (11,22,33,12,13,23); xl holds
!     the REFERENCE coordinates of the element nodes.
!
!     The band is taken normal to the major principal direction.  For an
!     isotropic material, damaged isotropically, the principal directions
!     of stress and of strain coincide, so this is the direction Jirasek
!     and Bauer describe as "the major principal strain axis"; using the
!     stress avoids reconstructing the strain here.
!
      implicit none
      character*8 lakonl
      integer nope,iok,ncor,i
      real*8 stre(6),xl(3,20),width,c(3,3),al(3),v(3),p,pmin,pmax
!
      iok=0
      width=0.d0
!
      call damcbcorner(nope,ncor)
      if(ncor.eq.0) return
!
      c(1,1)=stre(1)
      c(2,2)=stre(2)
      c(3,3)=stre(3)
      c(1,2)=stre(4)
      c(2,1)=stre(4)
      c(1,3)=stre(5)
      c(3,1)=stre(5)
      c(2,3)=stre(6)
      c(3,2)=stre(6)
!
!     calceigenvalues is the stock CalculiX routine and returns the three
!     eigenvalues in INCREASING order, so al(3) is the major one.
!
      call calceigenvalues(c,al)
      call damcbevec(c,al(3),v,iok)
      if(iok.eq.0) then
!
!       No determined band normal.  Fall back on the direction-free width
!       rather than on nothing: iok=0 out of here has to mean "degenerate
!       element", the one case that IS a data error, so that the two damage
!       entry points can act on it the same way.
!
        call damcbmean(nope,xl,width,iok)
        if(iok.eq.1) call damcbmeanonce()
        return
      endif
!
      pmin=1.d30
      pmax=-1.d30
      do i=1,ncor
        p=xl(1,i)*v(1)+xl(2,i)*v(2)+xl(3,i)*v(3)
        if(p.lt.pmin) pmin=p
        if(p.gt.pmax) pmax=p
      enddo
      width=pmax-pmin
      if(width.le.1.d-30) then
        iok=0
        return
      endif
!
      if((nope.eq.10).or.(nope.eq.20).or.(nope.eq.15))
     &     call damcbwarnonce(lakonl)
      return
      end
!
      subroutine calcdamage(ipkon,lakon,kon,co,mi,
     &     thicke,ielmat,ielprop,prop,ne0,ndmat_,ntmat_,
     &     ndmcon,dmcon,dam,dtime,sti,ithermal,t1,xstate,
     &     xstateini,nstate_,vold,idamage,imode,alphaevent)
!
!     Compatibility wrapper preserving the stock CalculiX 2.23 ABI.
!     Existing callers keep the original accumulated-damage behaviour.
!     nonlingeo.c in the active-set patch calls calcdamagebase directly
!     and supplies the immutable physical-increment damage baseline.
!
      implicit none
      character*8 lakon(*)
      integer ipkon(*),kon(*),mi(*),ielmat(mi(3),*),ielprop(*),ne0,
     &     ndmcon(2,*),ndmat_,ntmat_,ithermal(*),nstate_,idamage,imode
      real*8 co(3,*),thicke(mi(3),*),prop(*),
     &     dmcon(0:ndmat_,ntmat_,*),dam(mi(1),*),dtime,
     &     sti(6,mi(1),*),t1(*),xstate(nstate_,mi(1),*),
     &     xstateini(nstate_,mi(1),*),vold(0:mi(2),*),alphaevent
!
      call calcdamagebase(ipkon,lakon,kon,co,mi,
     &     thicke,ielmat,ielprop,prop,ne0,ndmat_,ntmat_,
     &     ndmcon,dmcon,dam,dam,dtime,sti,ithermal,t1,xstate,
     &     xstateini,nstate_,vold,idamage,imode,alphaevent)
      return
      end
!
      subroutine calcdamagebase(ipkon,lakon,kon,co,mi,
     &     thicke,ielmat,ielprop,prop,ne0,ndmat_,ntmat_,
     &     ndmcon,dmcon,dam,dambase,dtime,sti,ithermal,t1,xstate,
     &     xstateini,nstate_,vold,idamage,imode,alphaevent)
!     
!     calculates the damage due to plastic strain
!     
      implicit none
!     
      character*8 lakon(*),lakonl
!     
      integer ipkon(*),kon(*),mi(*),nope,indexe,i,j,k,ii,
     &     konl(20),mint3d,jj,iflag,ki,kl,ilayer,nlayer,kk,
     &     nopes,ielmat(mi(3),*),mint2d,null,ielprop(*),ne0,
     &     ndmcon(2,*),ndmat_,ntmat_,i1,nopered,ithermal(*),
     &     nstate_,id,imat,ndmconst,idamage,imode,ielemcross,
     &     ide1,imodenlv
!     
      real*8 co(3,*),prop(*),xl(3,20),xi,et,ze,xsj,shp(4,20),weight,
     &     a,gs(8,4),dlayer(4),tlayer(4),thickness,skl(3,3),s(3,3),
     &     thicke(mi(3),*),xlayer(mi(3),4),shp2(7,8),xs2(3,7),xsj2(3),
     &     xl2(3,8),eps0RT,xlimit,d1,d2,d3,d4,d5,Tmelt,Ttrans,eps0p,shy,
     &     svm,triax,t1l,dmcon(0:ndmat_,ntmat_,*),dam(mi(1),*),
     &     dambase(mi(1),*),dpeq,dpeqdt,dtime,ef,
     &     sti(6,mi(1),*),t1(*),vold(0:mi(2),*),
     &     xstate(nstate_,mi(1),*),xstateini(nstate_,mi(1),*),
     &     dmconloc(ndmat_),That,alphaevent,damold,ddam,damtrial,
     &     alphai,ufail,charlen,det6v,ax,ay,az,bx,by,bz,cx,cy,cz,
     &     cbwid,cbstre(6),
     &     damagebaseval,damagecur,damageD,damagetarget,dpeqpost,
     &     de1tol,ellnl
      integer iokv
!
!     CB1 crack-band width.  cbmodev is the switch, cbiok the
!     availability flag of the frozen value.
!
      integer cbmodev,cbiok,cbset,cbkind
      real*8 cbsvm,cbsh,cbufv
!
      real*8 dpnlv
!     
      include "gauss.f"
!
!     imode = 0: predictor only. The routine does not modify dam/ipkon.
!                idamage is the number of elements predicted to cross
!                the damage limit and alphaevent is the earliest estimated
!                crossing fraction within the current physical increment.
!
!     imode = 1: recompute the damage increment from dambase and perform
!                the usual hard element deletion.
!
!     dambase is the immutable damage state at the beginning of the
!     physical increment. Using it here is essential during same-load
!     re-equilibration: xstate-xstateini spans the complete physical
!     increment, so adding that increment to already-updated dam would
!     double count plastic strain.
!
!     The predictor uses a linear interpolation of the accumulated damage
!     increment. It is only an event locator; the physical increment is
!     re-solved by nonlingeo at the predicted load level before deletion.
!
!
!     Nonlocal average of the damage driving variable (E-83, E-84).
!     Inert unless CCX_DAMAGE_NONLOCAL has set an internal length, in
!     which case dpnl(i) replaces the LOCAL plastic-strain increment
!     further down.  The stress update stays local: only what DRIVES
!     damage is averaged.  That is the integral nonlocal DAMAGE model,
!     not nonlocal plasticity, and it is the whole change needed to put
!     a length into the formulation.
!
!
!     CB1: read the crack-band switch once and size the frozen-width
!     cache.  calcdamagebase is the only writer of that cache, because it
!     is the only one of the two damage entry points that sees a
!     COMMITTED stress; damageupdatepoint runs inside the stress update
!     and only reads it.
!
      call damcbinit(ne0,mi(1))
      call damcbufinit(ne0,mi(1))
      call damcbmodeget(cbmodev)
!
!
!     ENTRY GATE for the nonlocal model.  It must ask for the largest length
!     ANY source asked for, not just the environment's.
!
!     damnonlocalget returns ellsave, and ellsave is written only by
!     damnonlocalset, which nonlingeo.c calls with the environment value.  A
!     deck that carries NONLOCAL= on its material card and sets no
!     environment variable therefore left ellnl at zero and skipped this
!     block entirely: W6's whole point - that the .inp is a complete
!     description of the problem - was inert unless the environment ALSO
!     spoke.  MEASURED on a C3D8 deck: card only gave zero
!     "[DAMAGE NONLOCAL]" banners, the same deck with the environment set as
!     well gave four.
!
!     damnlellmax is Agent 1's accessor for exactly this question - the
!     maximum over ellsave and every card - so the gate asks it instead.
!
      call damnlellmax(ellnl)
      if(ellnl.gt.0.d0) then
!
!       The element-to-material map, for Agent 1's per-material internal
!       length (forum 2026-09-17 07:20).  damnonlocal.f does not receive
!       ielmat, and this is the main path into it - the other two callers
!       are inside that file.  Without this call the card's length is still
!       applied, but globally rather than per material, which is the
!       pre-W6 behaviour and would be a silent difference rather than a
!       visible one.  It returns immediately when no card carried
!       NONLOCAL=.
!
        call damnlmatmap(ipkon,lakon,ielmat,ne0,mi)
        call damnldtbase(dtime)
        call damnonlocalmodeget(imodenlv)
        if(imodenlv.eq.2) then
!
!         FROZEN-LOCAL control: the staggered update with NO averaging.
!         Separates the internal length from the integration scheme.
!
          call damfrozen(ipkon,kon,lakon,co,ne0,mi,xstate,
     &         xstateini,nstate_)
        elseif(imodenlv.eq.1) then
          call damgradient(ipkon,kon,lakon,co,ne0,mi,xstate,
     &         xstateini,nstate_,ielmat,dam)
        else
          call damnonlocal(ipkon,kon,lakon,co,ne0,mi,xstate,
     &         xstateini,nstate_)
        endif
      endif
      idamage=0
      alphaevent=2.d0
      de1tol=1.d-6
!     
      do i=1,ne0
        ielemcross=0
!     
!     element must exist and be a volume element
!     
        if((ipkon(i).lt.0).or.(lakon(i)(1:1).ne.'C')) cycle
!     
        lakonl=lakon(i)
        indexe=ipkon(i)
!     
        if(lakonl(1:5).eq.'C3D8I') then
          nope=11
        elseif(lakonl(4:4).eq.'2') then
          nope=20
        elseif(lakonl(4:4).eq.'8') then
          nope=8
        elseif(lakonl(4:5).eq.'10') then
          nope=10
        elseif(lakonl(4:4).eq.'4') then
          nope=4
        elseif(lakonl(4:5).eq.'15') then
          nope=15
        elseif(lakonl(4:5).eq.'6') then
          nope=6
        else
          cycle
        endif
!     
!     material
!     
        if(lakonl(7:8).ne.'LC') then
!     
!         no composite material: one material per element, all
!         integration points correspond to the same material
!     
          imat=ielmat(1,i)
          if(ndmcon(2,imat).eq.0) cycle
!     
!     determining the model for this element
!     
          if(int(dmcon(1,1,imat)).eq.1) then
!     
!     Rice-Tracey model
!     
            eps0RT=dmcon(2,1,imat)
            xlimit=dmcon(3,1,imat)
          elseif(int(dmcon(1,1,imat)).eq.2) then
!     
!     Johnson-Cook model
!     
            d1=dmcon(2,1,imat)
            d2=dmcon(3,1,imat)
            d3=dmcon(4,1,imat)
            d4=dmcon(5,1,imat)
            d5=dmcon(6,1,imat)
            Tmelt=dmcon(7,1,imat)
            Ttrans=dmcon(8,1,imat)
            eps0p=dmcon(9,1,imat)
            xlimit=dmcon(10,1,imat)
          elseif(int(dmcon(1,1,imat)).eq.3) then
!
!     DM2.0 tabulated ductile model
!
            xlimit=dmcon(2,1,imat)
          endif
        else
!     
!     composite materials
!     
!     determining the number of layers
!     
          nlayer=0
          do k=1,mi(3)
            if(ielmat(k,i).ne.0) then
              nlayer=nlayer+1
            endif
          enddo
!     
!     the thickness of the composite layers is only needed for
!     models requiring the temperature at the integration points
!     (so far only for the Johnson-Cook model)
!     
          if(int(dmcon(1,1,imat)).eq.2) then
            if(lakonl(4:4).eq.'2') then
              mint2d=4
              nopes=8
!     
!     determining the layer thickness and global thickness
!     at the shell integration points
!     
              iflag=1
              indexe=ipkon(i)
              do kk=1,mint2d
                xi=gauss3d2(1,kk)
                et=gauss3d2(2,kk)
                call shape8q(xi,et,xl2,xsj2,xs2,shp2,iflag)
                tlayer(kk)=0.d0
                do ii=1,nlayer
                  thickness=0.d0
                  do j=1,nopes
                    thickness=thickness+thicke(ii,indexe+j)*shp2(4,j)
                  enddo
                  tlayer(kk)=tlayer(kk)+thickness
                  xlayer(ii,kk)=thickness
                enddo
              enddo
              iflag=2
!     
              ilayer=0
              do ii=1,4
                dlayer(ii)=0.d0
              enddo
            elseif(lakonl(4:5).eq.'15') then
              mint2d=3
              nopes=6
!     
!     determining the layer thickness and global thickness
!     at the shell integration points
!     
              iflag=1
              indexe=ipkon(i)
              do kk=1,mint2d
                xi=gauss3d10(1,kk)
                et=gauss3d10(2,kk)
                call shape6tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
                tlayer(kk)=0.d0
                do ii=1,nlayer
                  thickness=0.d0
                  do j=1,nopes
                    thickness=thickness+thicke(ii,indexe+j)*shp2(4,j)
                  enddo
                  tlayer(kk)=tlayer(kk)+thickness
                  xlayer(ii,kk)=thickness
                enddo
              enddo
              iflag=2
!     
              ilayer=0
              do ii=1,3
                dlayer(ii)=0.d0
              enddo
            endif
          endif
!     
        endif
!     
        do j=1,nope
          konl(j)=kon(indexe+j)
          do k=1,3
            xl(k,j)=co(k,konl(j))
          enddo
        enddo
!
!       DE1 characteristic length for a linear tetrahedron.  The
!       coordinates in co are the reference mesh coordinates, so the
!       regularization length does not change during NLGEOM iterations.
!       det6v is six times the tetrahedral reference volume.
!
        charlen=0.d0
        if(lakonl(4:4).eq.'4') then
          ax=xl(1,2)-xl(1,1)
          ay=xl(2,2)-xl(2,1)
          az=xl(3,2)-xl(3,1)
          bx=xl(1,3)-xl(1,1)
          by=xl(2,3)-xl(2,1)
          bz=xl(3,3)-xl(3,1)
          cx=xl(1,4)-xl(1,1)
          cy=xl(2,4)-xl(2,1)
          cz=xl(3,4)-xl(3,1)
          det6v=dabs(ax*(by*cz-bz*cy)-ay*(bx*cz-bz*cx)
     &         +az*(bx*cy-by*cx))
          if(det6v.gt.1.d-30) charlen=det6v**(1.d0/3.d0)
        endif
!     
        if(lakonl(4:5).eq.'8R') then
          mint3d=1
        elseif(lakonl(4:7).eq.'20RB') then
          if((lakonl(8:8).eq.'R').or.(lakonl(8:8).eq.'C')) then
            mint3d=50
          else
            call beamintscheme(lakonl,mint3d,ielprop(i),prop,
     &           null,xi,et,ze,weight)
          endif
        elseif((lakonl(4:4).eq.'8').or.
     &         (lakonl(4:6).eq.'20R')) then
          if(lakonl(7:8).eq.'LC') then
            mint3d=8*nlayer
          else
            mint3d=8
          endif
        elseif(lakonl(4:4).eq.'2') then
          mint3d=27
        elseif(lakonl(4:5).eq.'10') then
          mint3d=4
        elseif(lakonl(4:4).eq.'4') then
          mint3d=1
        elseif(lakonl(4:5).eq.'15') then
          if(lakonl(7:8).eq.'LC') then
            mint3d=6*nlayer
          else
            mint3d=9
          endif
        elseif(lakonl(4:5).eq.'6') then
          mint3d=2
        else
          cycle
        endif
!     
        do jj=1,mint3d
!     
          if(lakonl(4:5).eq.'8R') then
            xi=gauss3d1(1,jj)
            et=gauss3d1(2,jj)
            ze=gauss3d1(3,jj)
            weight=weight3d1(jj)
          elseif(lakonl(4:7).eq.'20RB') then
            if((lakonl(8:8).eq.'R').or.(lakonl(8:8).eq.'C')) then
              xi=gauss3d13(1,jj)
              et=gauss3d13(2,jj)
              ze=gauss3d13(3,jj)
              weight=weight3d13(jj)
            else
              call beamintscheme(lakonl,mint3d,ielprop(i),prop,
     &             kk,xi,et,ze,weight)
            endif
          elseif((lakonl(4:4).eq.'8').or.
     &           (lakonl(4:6).eq.'20R'))
     &           then
            if(lakonl(7:8).ne.'LC') then
              xi=gauss3d2(1,jj)
              et=gauss3d2(2,jj)
              ze=gauss3d2(3,jj)
              weight=weight3d2(jj)
            else
              kl=mod(jj,8)
              if(kl.eq.0) kl=8
!     
              xi=gauss3d2(1,kl)
              et=gauss3d2(2,kl)
              ze=gauss3d2(3,kl)
              weight=weight3d2(kl)
!     
              ki=mod(jj,4)
              if(ki.eq.0) ki=4
!     
              if(kl.eq.1) then
                ilayer=ilayer+1
                if(ilayer.gt.1) then
                  do ii=1,4
                    dlayer(ii)=dlayer(ii)+xlayer(ilayer-1,ii)
                  enddo
                endif
              endif
              ze=2.d0*(dlayer(ki)+(ze+1.d0)/2.d0*xlayer(ilayer,ki))/
     &             tlayer(ki)-1.d0
              weight=weight*xlayer(ilayer,ki)/tlayer(ki)
              imat=ielmat(ilayer,i)
              if(ndmcon(2,imat).eq.0) cycle
            endif
          elseif(lakonl(4:4).eq.'2') then
            xi=gauss3d3(1,jj)
            et=gauss3d3(2,jj)
            ze=gauss3d3(3,jj)
            weight=weight3d3(jj)
          elseif(lakonl(4:5).eq.'10') then
            xi=gauss3d5(1,jj)
            et=gauss3d5(2,jj)
            ze=gauss3d5(3,jj)
            weight=weight3d5(jj)
          elseif(lakonl(4:4).eq.'4') then
            xi=gauss3d4(1,jj)
            et=gauss3d4(2,jj)
            ze=gauss3d4(3,jj)
            weight=weight3d4(jj)
          elseif(lakonl(4:5).eq.'15') then
            if(lakonl(7:8).ne.'LC') then
              xi=gauss3d8(1,jj)
              et=gauss3d8(2,jj)
              ze=gauss3d8(3,jj)
              weight=weight3d8(jj)
            else
              kl=mod(jj,6)
              if(kl.eq.0) kl=6
!     
              xi=gauss3d10(1,kl)
              et=gauss3d10(2,kl)
              ze=gauss3d10(3,kl)
              weight=weight3d10(kl)
!     
              ki=mod(jj,3)
              if(ki.eq.0) ki=3
!     
              if(kl.eq.1) then
                ilayer=ilayer+1
                if(ilayer.gt.1) then
                  do ii=1,3
                    dlayer(ii)=dlayer(ii)+xlayer(ilayer-1,ii)
                  enddo
                endif
              endif
              ze=2.d0*(dlayer(ki)+(ze+1.d0)/2.d0*xlayer(ilayer,ki))/
     &             tlayer(ki)-1.d0
              weight=weight*xlayer(ilayer,ki)/tlayer(ki)
              imat=ielmat(ilayer,i)
              if(ndmcon(2,imat).eq.0) cycle
            endif
          else
            xi=gauss3d7(1,jj)
            et=gauss3d7(2,jj)
            ze=gauss3d7(3,jj)
            weight=weight3d7(jj)
          endif
!     
!     determining the material model for this layer
!     (for composites only)
!     
          if(lakonl(7:8).eq.'LC') then
            if(int(dmcon(1,1,imat)).eq.1) then
!     
!     Rice-Tracey model
!     
              eps0RT=dmcon(2,1,imat)
              xlimit=dmcon(3,1,imat)
            elseif(int(dmcon(1,1,imat)).eq.2) then
!     
!     Johnson-Cook model
!     
              d1=dmcon(2,1,imat)
              d2=dmcon(3,1,imat)
              d3=dmcon(4,1,imat)
              d4=dmcon(5,1,imat)
              d5=dmcon(6,1,imat)
              Tmelt=dmcon(7,1,imat)
              Ttrans=dmcon(8,1,imat)
              eps0p=dmcon(9,1,imat)
              xlimit=dmcon(10,1,imat)
            elseif(int(dmcon(1,1,imat)).eq.3) then
!
!     DM2.0 tabulated ductile model
!
              xlimit=dmcon(2,1,imat)
            endif
          endif
!
!         DE1 is opt-in and is identified by the fourth Rice-Tracey
!         constant (u_f).  The first implementation is deliberately
!         restricted to C3D4 to keep the characteristic-length definition
!         unambiguous and to avoid changing stock behavior for other
!         elements.
!
          ide1=0
          if(((int(dmcon(1,1,imat)).eq.1).and.
     &       (ndmcon(1,imat).ge.4)).or.
     &       (int(dmcon(1,1,imat)).eq.3)) then
            ide1=1
!
!           The RAW slot, deliberately: this is the WRITER side.  Under
!           EVOLUTION=ENERGY it holds G_f, and the conversion to u_f a few
!           lines below is what turns it into a length.  Every READER must
!           go through damcbufail instead, which hands back value and
!           meaning together - see the comment on that routine.
!
            if(int(dmcon(1,1,imat)).eq.1) then
              ufail=dmcon(4,1,imat)
            else
              ufail=dmcon(3,1,imat)
            endif
!
!           The C3D4 restriction belongs to the VOLUME estimate, which is
!           the only thing that could not be defined for other shapes.
!           With the projection on, any volume family has a width.
!
            if((lakonl(4:4).ne.'4').and.(cbmodev.ne.1)) then
              write(*,*) '*ERROR in calcdamage: DE1 displacement'
              write(*,*) '       evolution with the legacy (6V)^(1/3)'
              write(*,*) '       characteristic length is implemented'
              write(*,*) '       for C3D4 elements only. Element: ',i
              write(*,*) '       Set CCX_DAMAGE_CHARLEN=1 for the'
              write(*,*) '       directional width, which is defined'
              write(*,*) '       for every volume element.'
              call exit(201)
            endif
!
!           CB1.  Refresh the projected width from the COMMITTED stress
!           while initiation is not yet committed, and stop refreshing
!           once it is: that freezes n and L at the increment in which
!           the crack forms.  Both inputs - sti and dambase - are
!           committed quantities, so the value is identical across the
!           Newton iterations of an increment and across rollback
!           retries of it, which is what keeps the law transactional.
!
!           THE SAME FREEZING STRUCTURE COST 8 PER CENT IN sigma_0, AND
!           IT COSTS 0.007 PER CENT HERE.  Measured, not assumed, because
!           the two caches freeze the same way and the question had to be
!           asked of both.  The u_f cache took the LAST value before the
!           threshold, and by then the neighbouring points of the band
!           have initiated and are unloading this one, so it sampled a
!           stress on the way down - see damcbufset.  Element 19 of the
!           equivalence deck, the last two refreshes:
!
!             dambase=0.931604   L=1.0000000
!             dambase=0.981067   L=0.9999312
!
!           The same unloading moves the width by seven thousandths of a
!           per cent.  The asymmetry is not luck: sigma_0 is a MAGNITUDE
!           and drops the moment the point leaves the yield surface,
!           while L depends on the principal DIRECTION, which is set by
!           how the element is held by its neighbours and barely turns.
!
!           So this cache needs no peak rule, and that is now a
!           measurement rather than a plausible argument.
!
            if(cbmodev.eq.1) then
              cbset=0
              if(dambase(jj,i).lt.xlimit) then
                cbset=1
              else
                call damcbget(i,jj,cbwid,cbiok)
                if(cbiok.eq.0) cbset=1
              endif
              if(cbset.eq.1) then
                cbstre(1)=sti(1,jj,i)
                cbstre(2)=sti(2,jj,i)
                cbstre(3)=sti(3,jj,i)
                cbstre(4)=sti(4,jj,i)
                cbstre(5)=sti(5,jj,i)
                cbstre(6)=sti(6,jj,i)
                call damcbwidth(cbstre,xl,nope,lakonl,cbwid,cbiok)
                if(cbiok.eq.1) call damcbset(i,jj,cbwid)
              endif
              call damcbget(i,jj,cbwid,cbiok)
              if(cbiok.eq.1) charlen=cbwid
            endif
!
!           CB-NL: a nonlocal driver needs the band's width, not the
!           element's.  Inert unless CCX_DAMAGE_NLWIDTH=1.
!
            call damcbnlwidth(i,charlen)
!
!           W2: ENERGY evolution.  The card slot that DISPLACEMENT reads as
!           u_f is G_f here, and the law needs u_f, so it is derived:
!
!             G_f = int_0^{u_f} (1 - d/u_f) sigma_0 dd = sigma_0 u_f / 2
!             =>  u_f = 2 G_f / sigma_0
!
!           exactly, for the linear-in-displacement law actually
!           implemented.  sigma_0 is the equivalent stress AT INITIATION at
!           THIS point, which under hardening differs from point to point -
!           that is the whole reason G_f is the material constant and u_f is
!           not.  One u_f card prescribes a different dissipated energy at
!           every point of a hardening model; one G_f card prescribes the
!           same energy and a different u_f.
!
!           Frozen exactly like the width, on the same committed inputs and
!           for the same reasons: refreshed while initiation has not been
!           committed, fixed once it has, so the last refresh is the
!           increment in which the crack forms.  sti is committed, so the
!           value cannot move with the Newton iterate or with a rollback.
!
            call damcbevolget(imat,cbkind)
            if(cbkind.eq.2) then
              if(dambase(jj,i).lt.xlimit) then
                cbsh=(sti(1,jj,i)+sti(2,jj,i)+sti(3,jj,i))/3.d0
                cbsvm=dsqrt(1.5d0*(
     &               (sti(1,jj,i)-cbsh)**2+(sti(2,jj,i)-cbsh)**2+
     &               (sti(3,jj,i)-cbsh)**2+2.d0*(sti(4,jj,i)**2+
     &               sti(5,jj,i)**2+sti(6,jj,i)**2)))
!
!               imode=1 ONLY.  The predictor pass is documented as not
!               modifying dam or ipkon, and the cache belongs in that set
!               for the same reason: nonlingeo's event controller can cut
!               the increment AFTER the predictor has run (it sets
!               icutb=1 together with damage_event_cut=1), and the load
!               level the predictor saw is then never committed.
!
!               This did not matter while the cache took the LATEST value,
!               because the redone increment simply overwrote it.  It
!               matters now that the cache keeps the PEAK: a cut increment
!               sits at a higher load than the redo, so its stress is
!               higher, so its u_f is smaller, and the smaller value would
!               be latched permanently by a pass that was thrown away.
!
!               Found by asking what the peak rule does to the file's own
!               transactional claim - "a Newton retry or a rollback cannot
!               change it" - rather than by a test, and no deck here shows
!               it: the equivalence numbers are unchanged to four decimals.
!
                if((cbsvm.gt.1.d-10).and.(imode.eq.1)) then
                  call damcbufset(i,jj,2.d0*ufail/cbsvm)
                endif
              endif
              call damcbufget(i,jj,cbufv,cbiok)
              if(cbiok.eq.1) then
                ufail=cbufv
              else
!
!               No equivalent stress has been seen yet, so G_f cannot be
!               converted.  Leave the point alone this increment rather
!               than guess: D is still at its baseline here.
!
                cycle
              endif
            endif
!
            if((ufail.le.0.d0).or.(charlen.le.0.d0)) then
              write(*,*) '*ERROR in calcdamage: invalid DE1 data'
              write(*,*) '       element=',i,' u_f=',ufail,
     &                   ' L=',charlen
              call exit(201)
            endif
          endif
!     
!     shape functions need only be determined if the     
!     temperature is needed, i.e. for the Johnson-Cook model
!     
          if(int(dmcon(1,1,imat)).eq.2) then
            iflag=1
            if(lakonl(1:5).eq.'C3D8R') then
              call shape8hr(xl,xsj,shp,gs,a)
            elseif(lakonl(1:5).eq.'C3D8I') then
              call shape8hu(xi,et,ze,xl,xsj,shp,iflag)
            elseif(nope.eq.20) then
              call shape20h(xi,et,ze,xl,xsj,shp,iflag)
            elseif(nope.eq.8) then
              call shape8h(xi,et,ze,xl,xsj,shp,iflag)
            elseif(nope.eq.10) then
              call shape10tet(xi,et,ze,xl,xsj,shp,iflag)
            elseif(nope.eq.4) then
              call shape4tet(xi,et,ze,xl,xsj,shp,iflag)
            elseif(nope.eq.15) then
              call shape15w(xi,et,ze,xl,xsj,shp,iflag)
            else
              call shape6w(xi,et,ze,xl,xsj,shp,iflag)
            endif
          endif
!     
!     stress at the integration point
!     
          skl(1,1)=sti(1,jj,i)
          skl(2,2)=sti(2,jj,i)
          skl(3,3)=sti(3,jj,i)
          skl(1,2)=sti(4,jj,i)
          skl(1,3)=sti(5,jj,i)
          skl(2,3)=sti(6,jj,i)
!     
!     hydrostatic stress
!     
          shy=(skl(1,1)+skl(2,2)+skl(3,3))/3.d0
!     
!     deviatoric stress tensor
!     
          s(1,1)=skl(1,1)-shy
          s(2,2)=skl(2,2)-shy
          s(3,3)=skl(3,3)-shy
          s(1,2)=skl(1,2)
          s(1,3)=skl(1,3)
          s(2,3)=skl(2,3)
!     
!     von Mises stress
!     
          svm=dsqrt(3.d0/2.d0*(
     &         s(1,1)*s(1,1)+s(2,2)*s(2,2)+s(3,3)*s(3,3)+
     &         2.d0*(s(1,2)*s(1,2)+s(1,3)*s(1,3)+s(2,3)*s(2,3))))
!     
!     triaxiality
!     
          if(svm.gt.1.d-20) then
            triax=shy/svm
          else
            triax=0.d0
          endif
!     
!         calculate the temperature (only needed for the Johnson-Cook
!         model
!     
          if(int(dmcon(1,1,imat)).eq.2) then
            t1l=Ttrans
            if(ithermal(1).ge.1) then
              t1l=0.d0
              if(ithermal(1).eq.1) then
                if((lakonl(4:5).eq.'8 ').or.
     &               (lakonl(4:5).eq.'8I')) then
                  do i1=1,8
                    t1l=t1l+t1(konl(i1))/8.d0
                  enddo
                elseif(lakonl(4:6).eq.'20 ') then
                  nopered=20
                  call lintemp(t1,konl,nopered,jj,t1l)
                elseif(lakonl(4:6).eq.'10T') then
                  call linscal10(t1,konl,t1l,null,shp)
                else
                  do i1=1,nope
                    t1l=t1l+shp(4,i1)*t1(konl(i1))
                  enddo
                endif
              elseif(ithermal(1).ge.2) then
                if((lakonl(4:5).eq.'8 ').or.
     &               (lakonl(4:5).eq.'8I')) then
                  do i1=1,8
                    t1l=t1l+vold(0,konl(i1))/8.d0
                  enddo
                elseif(lakonl(4:6).eq.'20 ') then
                  nopered=20
                  call lintemp_th1(vold,konl,nopered,jj,t1l,mi)
                elseif(lakonl(4:6).eq.'10T') then
                  call linscal10(vold,konl,t1l,mi(2),shp)
                else
                  do i1=1,nope
                    t1l=t1l+shp(4,i1)*vold(0,konl(i1))
                  enddo
                endif
              endif
            endif
          endif
!     
!     interpolating the material data
!     (for models with temperature dependent parameters;
!     no such model is implemented so far; a model
!     number of 3 or higher is assumed)          
!     
          if(int(dmcon(1,1,imat)).gt.3) then
!     
!           number of constants in this model    
!     
            ndmconst=ndmcon(1,imat)
            if(ithermal(1).eq.0) then
              do k=1,ndmconst
                dmconloc(k)=dmcon(k,1,imat)
              enddo
            else
              call ident2(dmcon(0,1,imat),t1l,ndmcon(2,imat),ndmat_+1,
     &             id)
              if(ndmcon(2,imat).eq.1) then
                do k=1,ndmconst
                  dmconloc(k)=dmcon(k,1,imat)
                enddo
              elseif(id.eq.0) then
                do k=1,ndmconst
                  dmconloc(k)=dmcon(k,1,imat)
                enddo
              elseif(id.eq.ndmcon(2,imat)) then
                do k=1,ndmconst
                  dmconloc(k)=dmcon(k,id,imat)
                enddo
              else
                do k=1,ndmconst
                  dmconloc(k)=dmcon(k,id,imat)+
     &                 (dmcon(k,id+1,imat)-dmcon(k,id,imat))*
     &                 (t1l-dmcon(0,id,imat))/
     &                 (dmcon(0,id+1,imat)-dmcon(0,id,imat))
                enddo
              endif
            endif
          endif          
!     
!     change in equivalent plastic strain
!     
          dpeq=xstate(1,jj,i)-xstateini(1,jj,i)
          if(ellnl.gt.0.d0) then
            call damnonlocalval(i,dpnlv,iokv)
            if(iokv.eq.1) dpeq=dpnlv
          endif
          damagecur=dam(jj,i)
          damagebaseval=dambase(jj,i)

          if(dtime.gt.1.d-30) then
            dpeqdt=dpeq/dtime
          else
            dpeqdt=0.d0
          endif
!
!         ------------------------------------------------------------
!         DE1/DM2: initiation + displacement-based evolution
!         ------------------------------------------------------------
!
!         Storage convention (DE1 only, xlimit is required to be 1):
!           0 <= dam < 1 : Rice-Tracey initiation variable omega_D
!           1 <= dam <=2 : 1 + scalar degradation D
!
!         Thus existing DUCT output remains monotone and no new global
!         history array is required.  The physical degradation variable is
!         D=max(0,dam-1).  Fully damaged elements are not removed from the
!         topology in DE1; resultsmech applies a residual stiffness floor.
!
          if(ide1.eq.1) then

!           DE1.2 mode dispatch:
!             imode=0: legacy/A3 predictor only; DE1 is skipped.
!             imode=1: legacy/A3 apply only; DE1 is skipped because its
!                      trial state is already integrated inside Newton.
!             imode=2: DE1.2 trial update from the immutable physical-
!                      increment baseline. Legacy damage is skipped below.

            if(imode.ne.2) then
              cycle
            endif

!           Recompute from the immutable physical-increment baseline.
!           This makes DE1 transactional under Newton retries and A3
!           rollback, exactly like the A2/A3 initiation history.

            if(damagebaseval.ge.xlimit) then

!             Damage evolution was already active at increment start.

              damageD=damagebaseval-xlimit
              if(damageD.lt.0.d0) damageD=0.d0
              if(damageD.gt.1.d0) damageD=1.d0
              dpeqpost=dmax1(0.d0,dpeq)
              damageD=damageD+charlen*dpeqpost/ufail
              if(damageD.gt.1.d0) damageD=1.d0
              damagetarget=xlimit+damageD

            else

!             Still in the initiation stage.  The selected criterion
!             determines the fraction at which omega_D reaches
!             one.  Only the remaining plastic strain drives evolution.

              damagetarget=damagebaseval
              if(dpeq.gt.0.d0) then
                if(int(dmcon(1,1,imat)).eq.1) then
                  ef=1.65d0*eps0RT*dexp(-3.d0*triax/2.d0)
                elseif(int(dmcon(1,1,imat)).eq.3) then
                  call dm2failurestrain(dmcon,ndmat_,ntmat_,imat,
     &                 ndmcon,triax,ef)
                else
                  ef=-1.d0
                endif
                if(ef.gt.1.d-14) then
                  ddam=dpeq/ef
                  damtrial=damagebaseval+ddam
                  if(damtrial.lt.xlimit) then
                    damagetarget=damtrial
                  else
                    alphai=(xlimit-damagebaseval)/ddam
                    if(alphai.lt.0.d0) alphai=0.d0
                    if(alphai.gt.1.d0) alphai=1.d0
                    dpeqpost=(1.d0-alphai)*dpeq
                    damageD=charlen*dpeqpost/ufail
                    if(damageD.lt.0.d0) damageD=0.d0
                    if(damageD.gt.1.d0) damageD=1.d0
                    damagetarget=xlimit+damageD
                  endif
                else
                  write(*,*) '*WARNING in calcdamage: non-positive',
     &                       ' failure strain in element ',i,
     &                       ', integration point ',jj
                endif
              endif
            endif

            dam(jj,i)=damagetarget

!           idamage means "the constitutive damage state changed enough
!           to require another same-load equilibrium", not deletion, for
!           DE1.  Count an element only once even if it has several IPs.

            if(((damagetarget.gt.xlimit+de1tol).or.
     &          (damagecur.gt.xlimit+de1tol)).and.
     &          (dabs(damagetarget-damagecur).gt.de1tol)) then
              if(ielemcross.eq.0) idamage=idamage+1
              ielemcross=1
            endif

            cycle
          endif
!
!         ------------------------------------------------------------
!         Legacy CalculiX/A3 damage path (unchanged hard deletion)
!         ------------------------------------------------------------
!
!         imode=2 is reserved for the DE1.2 Newton-trial update.  Legacy
!         damage is evaluated only by the original predictor/apply modes.
!
          if(imode.eq.2) cycle
!
!         In apply/recompute mode the result at every surviving damage
!         integration point must be reconstructed from dambase.
!
          if(imode.ne.0) dam(jj,i)=dambase(jj,i)

          if(dpeq.gt.0.d0) then
!     
            if(int(dmcon(1,1,imat)).eq.1) then
!     
!     Rice-Tracey model
!     
              ef=1.65d0*eps0RT*dexp(-3.d0*triax/2.d0)
!     
            elseif(int(dmcon(1,1,imat)).eq.2) then
!     
!     Johnson-Cook model
!
              if(t1l.lt.Ttrans) then
                That=0.d0
              elseif(t1l.gt.Tmelt) then
                That=1.d0
              else
                That=(t1l-Ttrans)/(Tmelt-Ttrans)
              endif
!
              if(dpeqdt.gt.eps0p) then
                ef=(d1+d2*dexp(-d3*triax))*(1.d0+d4*dlog(dpeqdt/eps0p))*
     &               (1.d0+d5*That)
              else
!
!               replacing the logarithmic function by a linear function
!
                ef=(d1+d2*dexp(-d3*triax))*
     &               (1.d0+d4*(dpeqdt/eps0p-1.d0))*(1.d0+d5*That)
              endif
            endif
!     
!     damage increment
!     
            if(ef.gt.1.d-14) then
              damold=dambase(jj,i)
              ddam=dpeq/ef
              damtrial=damold+ddam
            else
              write(*,*) '*WARNING in calcdamage: non-positive',
     &                   ' failure strain in element ',i,
     &                   ', integration point ',jj
              cycle
            endif
!
!     predictor mode: determine whether the element crosses the damage
!     limit during this physical increment, without changing dam/ipkon.
!
            if(imode.eq.0) then
              if(damold.ge.xlimit) then
                alphai=0.d0
                if(ielemcross.eq.0) idamage=idamage+1
                ielemcross=1
                if(alphai.lt.alphaevent) alphaevent=alphai
              elseif((ddam.gt.1.d-30).and.
     &               (damtrial.ge.xlimit)) then
                alphai=(xlimit-damold)/ddam
                if(alphai.lt.0.d0) alphai=0.d0
                if(alphai.gt.1.d0) alphai=1.d0
                if(ielemcross.eq.0) idamage=idamage+1
                ielemcross=1
                if(alphai.lt.alphaevent) alphaevent=alphai
              endif
!
!     apply mode: preserve the original CalculiX hard-deletion law.
!
            else
              dam(jj,i)=damtrial
              if(dam(jj,i).ge.xlimit) then
                ipkon(i)=-ipkon(i)-2
                idamage=idamage+1
                exit
              endif
            endif
          endif
        enddo
      enddo
!     
      return
      end

!
!     BK1 integration-point update for progressive DE1/DM2 damage.
!     This routine is called immediately after the effective constitutive
!     update in resultsmech.  It deliberately implements only the existing
!     C3D4 displacement-evolution path; legacy Rice-Tracey/Johnson-Cook
!     hard deletion continues to use calcdamagebase.
!
      subroutine damageupdatepoint(iel,iint,ipkon,lakon,kon,co,mi,
     &     ielmat,ne0,ndmat_,ntmat_,ndmcon,dmcon,dam,dambase,
     &     stre,xstate,xstateini,nstate_)
      implicit none
      real*8 ellnlp,dpnlp
      integer ioknl
!
      character*8 lakon(*)
      integer iel,iint,ipkon(*),kon(*),mi(*),ielmat(mi(3),*),ne0,
     &     ndmat_,ntmat_,ndmcon(2,*),nstate_,imat,itype,indexe,
     &     n1,n2,n3,n4
      real*8 co(3,*),dmcon(0:ndmat_,ntmat_,*),dam(mi(1),*),
     &     dambase(mi(1),*),stre(6),xstate(nstate_,mi(1),*),
     &     xstateini(nstate_,mi(1),*),shy,s1,s2,s3,svm,triax,
     &     dpeq,xlimit,ufail,ef,damagebaseval,damagetarget,
     &     damageD,dpeqpost,ddam,damtrial,alphai,charlen,
     &     ax,ay,az,bx,by,bz,cx,cy,cz,det6v,cbwid,xlp(3,20),cbufv
      integer cbmodev,cbiok,nope4,ncor4,i4,k4,cbkind
!
      if((iel.lt.1).or.(iel.gt.ne0)) return
      if(ipkon(iel).lt.0) return
      if(lakon(iel)(1:1).ne.'C') return
      if(lakon(iel)(7:8).eq.'LC') then
        write(*,*) '*ERROR in damageupdatepoint: progressive damage'
        write(*,*) '       is not supported for composite elements.'
        call exit(201)
      endif
!
      imat=ielmat(1,iel)
      if(imat.le.0) return
      if(ndmcon(2,imat).eq.0) return
      itype=int(dmcon(1,1,imat))
      if(itype.eq.1) then
        if(ndmcon(1,imat).lt.4) return
        xlimit=dmcon(3,1,imat)
      elseif(itype.eq.3) then
        xlimit=dmcon(2,1,imat)
      else
        return
      endif
!
!     u_f comes from damcbufail and from nowhere else: under
!     EVOLUTION=ENERGY the slot holds G_f, and value and meaning have to
!     leave one routine together or a reader will take one for the other.
!     iok=0 means it is not available yet - only possible under ENERGY
!     before an equivalent stress has been committed here - and then this
!     point must not evolve.
!
      call damcbufail(imat,itype,dmcon,ndmat_,ntmat_,iel,iint,ufail,
     &     cbiok)
      if(cbiok.eq.0) return
      call damcbmodeget(cbmodev)
!
!     Node count of this element, needed for the corner count.
!
      if(lakon(iel)(1:5).eq.'C3D8I') then
        nope4=11
      elseif(lakon(iel)(4:4).eq.'2') then
        nope4=20
      elseif(lakon(iel)(4:4).eq.'8') then
        nope4=8
      elseif(lakon(iel)(4:5).eq.'10') then
        nope4=10
      elseif(lakon(iel)(4:4).eq.'4') then
        nope4=4
      elseif(lakon(iel)(4:5).eq.'15') then
        nope4=15
      elseif(lakon(iel)(4:4).eq.'6') then
        nope4=6
      else
        nope4=0
      endif
!
      if((lakon(iel)(4:4).ne.'4').and.(cbmodev.ne.1)) then
        write(*,*) '*ERROR in damageupdatepoint: DE1 displacement'
        write(*,*) '       evolution with the legacy (6V)^(1/3)'
        write(*,*) '       characteristic length is implemented for'
        write(*,*) '       C3D4 elements only. Element: ',iel
        call exit(201)
      endif
      if(ufail.le.0.d0) return
!     CB1: the projected width, when it is armed, comes from the cache
!     that calcdamagebase filled from the COMMITTED stress.  This routine
!     runs inside the stress update, where the only stress available is
!     the trial one, so it never computes the direction itself - that is
!     what keeps the width out of the tangent and out of the Newton path.
!
!     A point that initiates before calcdamagebase has ever seen it has
!     no cached value yet; it uses the volume estimate for that one
!     increment and the projection from the next one on.  The difference
!     is confined to the increment in which D leaves zero.
!
      charlen=0.d0
      if(cbmodev.eq.1) then
        call damcbget(iel,iint,cbwid,cbiok)
        if(cbiok.eq.1) then
          charlen=cbwid
          goto 100
        endif
!
!       Not cached yet.  This happens on the increment in which a point
!       initiates, because calcdamagebase - the only writer - runs after
!       the increment converges, while this routine runs inside it.
!
!       The width used here is the DIRECTION-FREE one, from geometry alone.
!       Using the projection would mean taking a direction from stre, which
!       is the TRIAL stress: that would make the width a function of the
!       Newton iterate and undo exactly what freezing it on committed data
!       was for.  The geometric width depends on nothing that a Newton
!       retry or a rollback can change, so the law stays transactional, and
!       it is defined for every volume family - which is what lets the
!       C3D4 restriction be lifted here as well as in calcdamagebase.
!       Previously this path returned silently for a non-tetrahedron while
!       calcdamagebase stopped the run on the same condition; the two now
!       agree.
!
        call damcbcorner(nope4,ncor4)
        if(ncor4.gt.0) then
          indexe=ipkon(iel)
          do i4=1,ncor4
            do k4=1,3
              xlp(k4,i4)=co(k4,kon(indexe+i4))
            enddo
          enddo
          call damcbmean(nope4,xlp,cbwid,cbiok)
          if(cbiok.eq.1) then
            charlen=cbwid
            goto 100
          endif
        endif
        if(lakon(iel)(4:4).ne.'4') then
          write(*,*) '*ERROR in calcdamage: no crack-band width could'
          write(*,*) '       be formed for element ',iel
          call exit(201)
        endif
      endif
!
!     Reference characteristic length L=(6 V0)^(1/3).
!
      indexe=ipkon(iel)
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
     &     +az*(bx*cy-by*cx))
      if(det6v.le.1.d-30) return
      charlen=det6v**(1.d0/3.d0)
!
  100 continue
!
!     CB-NL: the same substitution as in calcdamagebase, applied at the
!     one point both paths have a width.  Inert unless the switch is set.
!
      call damcbnlwidth(iel,charlen)
!
!     Triaxiality from the effective stress (scalar degradation would
!     leave it invariant, but using the effective stress is unambiguous).
!
      shy=(stre(1)+stre(2)+stre(3))/3.d0
      s1=stre(1)-shy
      s2=stre(2)-shy
      s3=stre(3)-shy
      svm=dsqrt(1.5d0*(s1*s1+s2*s2+s3*s3+
     &     2.d0*(stre(4)*stre(4)+stre(5)*stre(5)+
     &     stre(6)*stre(6))))
      if(svm.gt.1.d-20) then
        triax=shy/svm
      else
        triax=0.d0
      endif
!
!     Rebuild the complete trial state from the committed baseline.
!
      dpeq=xstate(1,iint,iel)-xstateini(1,iint,iel)
!
!     Nonlocal substitution.  THIS is the path that drives the response:
!     resultsmech calls damageupdatepoint per integration point inside
!     the stress update, and patching only calcdamagebase left the answer
!     bit-identical to the local model.
!
!     The field is refreshed once per calcdamagebase call, so what is read
!     here is the last refresh rather than the current trial state - the
!     explicit (staggered) nonlocal scheme.  The lag is one damage update
!     and the tangent is local either way, so nothing is given up that
!     was not already given up.
!
!     The same entry gate as in calcdamagebase, and for the same reason: a
!     length that came from the card alone must arm this path too.
      call damnlellmax(ellnlp)
      if(ellnlp.gt.0.d0) then
        call damnonlocalval(iel,dpnlp,ioknl)
        if(ioknl.eq.1) dpeq=dpnlp
      endif
      damagebaseval=dambase(iint,iel)
      damagetarget=damagebaseval
!
      if(damagebaseval.ge.xlimit) then
        damageD=damagebaseval-xlimit
        if(damageD.lt.0.d0) damageD=0.d0
        if(damageD.gt.1.d0) damageD=1.d0
        dpeqpost=dmax1(0.d0,dpeq)
        damageD=damageD+charlen*dpeqpost/ufail
        if(damageD.gt.1.d0) damageD=1.d0
        damagetarget=xlimit+damageD
      elseif(dpeq.gt.0.d0) then
        if(itype.eq.1) then
          ef=1.65d0*dmcon(2,1,imat)*dexp(-1.5d0*triax)
        else
          call dm2failurestrain(dmcon,ndmat_,ntmat_,imat,
     &         ndmcon,triax,ef)
        endif
        if(ef.gt.1.d-14) then
          ddam=dpeq/ef
          damtrial=damagebaseval+ddam
          if(damtrial.lt.xlimit) then
            damagetarget=damtrial
          else
            alphai=(xlimit-damagebaseval)/ddam
            if(alphai.lt.0.d0) alphai=0.d0
            if(alphai.gt.1.d0) alphai=1.d0
            dpeqpost=(1.d0-alphai)*dpeq
            damageD=charlen*dpeqpost/ufail
            if(damageD.lt.0.d0) damageD=0.d0
            if(damageD.gt.1.d0) damageD=1.d0
            damagetarget=xlimit+damageD
          endif
        endif
      endif
!
      dam(iint,iel)=damagetarget
      return
      end

!
!     DM2.0: linearly interpolate the tabulated fracture locus
!     eps_f(eta).  The table is stored as
!       dmcon(2)=xlimit, dmcon(3)=u_f,
!       dmcon(4:)=eta_1,eps_f_1,eta_2,eps_f_2,...
!     Values outside the tabulated eta range are clamped to the end points.
!
      subroutine dm2failurestrain(dmcon,ndmat_,ntmat_,imat,ndmcon,
     &     triax,ef)
      implicit none
      integer ndmat_,ntmat_,imat,ndmcon(2,*),npoints,ipt,
     &     idxeta,idxeps
      real*8 dmcon(0:ndmat_,ntmat_,*),triax,ef,eta1,eta2,
     &     eps1,eps2,fac
!
      npoints=(ndmcon(1,imat)-3)/2
      if(npoints.lt.2) then
        ef=-1.d0
        return
      endif
!
      eta1=dmcon(4,1,imat)
      eps1=dmcon(5,1,imat)
      if(triax.le.eta1) then
        ef=eps1
        return
      endif
!
      idxeta=2*npoints+2
      idxeps=idxeta+1
      eta2=dmcon(idxeta,1,imat)
      eps2=dmcon(idxeps,1,imat)
      if(triax.ge.eta2) then
        ef=eps2
        return
      endif
!
      do ipt=1,npoints-1
        idxeta=2*ipt+2
        idxeps=idxeta+1
        eta1=dmcon(idxeta,1,imat)
        eps1=dmcon(idxeps,1,imat)
        eta2=dmcon(idxeta+2,1,imat)
        eps2=dmcon(idxeps+2,1,imat)
        if((triax.ge.eta1).and.(triax.le.eta2)) then
          fac=(triax-eta1)/(eta2-eta1)
          ef=eps1+fac*(eps2-eps1)
          return
        endif
      enddo
!
      ef=-1.d0
      return
      end
