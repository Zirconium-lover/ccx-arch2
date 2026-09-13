!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2025 Guido Dhondt
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
      subroutine damageinitiations(inpc,textpart,matname,nmat,
     &  irstrt,istep,istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,
     &  ier,dmcon,ndmcon,ntmat_,ndmat_)
!
!     reading the input deck: *DAMAGE INITIATION
!
      implicit none
!
      character*1 inpc(*)
      character*80 matname(*)
      character*132 textpart(16)
!
      integer nmat,istep,istat,n,key,i,irstrt(*),iline,ipol,inl,
     &     ipoinp(2,*),inp(3,*),ipoinpc(0:*),ier,ndmcon(2,*),
     &     ntmat_,j,nconstants,isum,imax,ntmat,ndmat_,itype,
     &     ievolution,npoints,ipt,idxeta,idxeps
!
      real*8 dmcon(0:ndmat_,ntmat_,*),etaold
!
      ntmat=0
      itype=0
      ievolution=0
      nconstants=0
      npoints=0
!
      if((istep.gt.0).and.(irstrt(1).ge.0)) then
         write(*,*) '*ERROR reading *DAMAGE INITIATION:'
         write(*,*) '       *DAMAGE INITIATION should be placed'
         write(*,*) '       before all step definitions'
         ier=1
         return
      endif
!
      do i=2,n
        if(textpart(i)(1:10).eq.'CRITERION=') then
          if(textpart(i)(11:20).eq.'RICETRACEY') then
            itype=1
          elseif(textpart(i)(11:21).eq.'JOHNSONCOOK') then
            itype=2
          elseif(textpart(i)(11:17).eq.'DUCTILE') then
            itype=3
          else
            write(*,*)
     &           '*ERROR reading *DAMAGE INITIATION: criterion'
            write(*,*) '       is not known.'
            write(*,*) '       criterion:',textpart(i)(11:90)
            ier=1
            return
          endif
        elseif(textpart(i)(1:10).eq.'EVOLUTION=') then
          if(textpart(i)(11:22).eq.'DISPLACEMENT') then
            ievolution=1
          else
            write(*,*)
     &           '*ERROR reading *DAMAGE INITIATION: evolution'
            write(*,*) '       is not known.'
            write(*,*) '       evolution:',textpart(i)(11:90)
            ier=1
            return
          endif
        elseif(textpart(i)(1:8).eq.'NPOINTS=') then
          read(textpart(i)(9:18),'(i10)',iostat=istat) npoints
          if(istat.ne.0) then
            write(*,*) '*ERROR reading *DAMAGE INITIATION:'
            write(*,*) '       invalid NPOINTS value.'
            ier=1
            return
          endif
        else
          write(*,*) '*WARNING reading *DAMAGE INITIATION:'
          write(*,*) '         parameter not recognized:'
          write(*,*) '         ',
     &         textpart(i)(1:index(textpart(i),' ')-1)
          call inputwarning(inpc,ipoinpc,iline,
     &         "*DAMAGE INITIATION%")
        endif
      enddo
!
      if(itype.eq.1) then
        if(ievolution.eq.1) then
          nconstants=4
        else
          nconstants=3
        endif
      elseif(itype.eq.2) then
        if(ievolution.eq.1) then
          write(*,*) '*ERROR reading *DAMAGE INITIATION:'
          write(*,*) '       DE1 displacement evolution is presently'
          write(*,*) '       not implemented for Johnson-Cook.'
          ier=1
          return
        endif
        nconstants=10
      elseif(itype.eq.3) then
        if(ievolution.ne.1) then
          write(*,*) '*ERROR reading *DAMAGE INITIATION:'
          write(*,*) '       DM2 DUCTILE requires'
          write(*,*) '       EVOLUTION=DISPLACEMENT.'
          ier=1
          return
        endif
        if(npoints.lt.2) then
          write(*,*) '*ERROR reading *DAMAGE INITIATION:'
          write(*,*) '       DM2 DUCTILE requires NPOINTS >= 2.'
          ier=1
          return
        endif
        nconstants=3+2*npoints
      else
        write(*,*) '*ERROR reading *DAMAGE INITIATION:'
        write(*,*) '       CRITERION parameter is required.'
        ier=1
        return
      endif
!
!     For DE1 Rice-Tracey:
!       dmcon(2) = eps0RT
!       dmcon(3) = initiation limit omega_D
!       dmcon(4) = failure plastic displacement u_f
!
!     For DM2.0 DUCTILE with NPOINTS=N:
!       dmcon(2) = initiation limit omega_D (must be 1)
!       dmcon(3) = failure plastic displacement u_f
!       dmcon(4),dmcon(5) = eta_1, eps_f_1
!       dmcon(6),dmcon(7) = eta_2, eps_f_2, etc.
!       eps_f(eta) is linearly interpolated and clamped at the ends.
!
!     the damage initiation is stored as a mechanical user material
!
      ndmcon(1,nmat)=nconstants
!
      do
        do j=1,(nconstants-1)/8+1
          if(j.eq.1) then
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &           inl,ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) return
            ntmat=ntmat+1
            ndmcon(2,nmat)=ntmat
            if((itype.eq.3).and.(ntmat.gt.1)) then
              write(*,*) '*ERROR reading *DAMAGE INITIATION:'
              write(*,*) '       DM2.0 does not yet support'
              write(*,*) '       temperature-dependent locus tables.'
              ier=1
              return
            endif
            if(ntmat.gt.ntmat_) then
              write(*,*) 
     &             '*ERROR reading *DAMAGE INITIATION: increase ntmat_'
              ier=1
              return
            endif
            isum=1
            dmcon(1,ntmat,nmat)=itype+0.5d0
          else
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &           inl,ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) then
              write(*,*) 
     &          '*ERROR reading *DAMAGE INITIATION: damage definition'
              write(*,*) '  is not complete. '
              call inputerror(inpc,ipoinpc,iline,
     &             "*DAMAGE INITIATION%",ier)
              return
            endif
          endif
          imax=8
          if(isum+imax.gt.nconstants+1) then
            imax=nconstants-isum+1
          endif
          do i=1,imax
            if(isum+i.le.nconstants) then
              read(textpart(i)(1:20),'(f20.0)',iostat=istat) 
     &             dmcon(isum+i,ntmat,nmat)
            else
              read(textpart(i)(1:20),'(f20.0)',iostat=istat) 
     &             dmcon(0,ntmat,nmat)
            endif
            if(istat.gt.0) then
              call inputerror(inpc,ipoinpc,iline,
     &             "*DAMAGE INITIATION%",ier)
              return
            endif
          enddo
          isum=isum+imax
!
!         check constants for Johnson-Cook
!
          if(isum.eq.nconstants+1) then
            if((itype.eq.1).and.(ievolution.eq.1)) then
              if(dabs(dmcon(3,ntmat,nmat)-1.d0).gt.1.d-12) then
                write(*,*) '*ERROR reading *DAMAGE INITIATION:'
                write(*,*) '       DE1 presently requires the'
                write(*,*) '       Rice-Tracey initiation limit = 1.'
                ier=1
                return
              endif
              if(dmcon(4,ntmat,nmat).le.0.d0) then
                write(*,*) '*ERROR reading *DAMAGE INITIATION:'
                write(*,*) '       failure plastic displacement u_f'
                write(*,*) '       must be strictly positive for DE1.'
                ier=1
                return
              endif
            endif
            if(itype.eq.3) then
              if(dabs(dmcon(2,ntmat,nmat)-1.d0).gt.1.d-12) then
                write(*,*) '*ERROR reading *DAMAGE INITIATION:'
                write(*,*) '       DM2 DUCTILE presently requires'
                write(*,*) '       initiation limit omega_D = 1.'
                ier=1
                return
              endif
              if(dmcon(3,ntmat,nmat).le.0.d0) then
                write(*,*) '*ERROR reading *DAMAGE INITIATION:'
                write(*,*) '       DM2 failure displacement u_f'
                write(*,*) '       must be strictly positive.'
                ier=1
                return
              endif
              etaold=-1.d300
              do ipt=1,npoints
                idxeta=2*ipt+2
                idxeps=idxeta+1
                if(dmcon(idxeps,ntmat,nmat).le.0.d0) then
                  write(*,*) '*ERROR reading *DAMAGE INITIATION:'
                  write(*,*) '       DM2 eps_f must be positive.'
                  ier=1
                  return
                endif
                if((ipt.gt.1).and.
     &             (dmcon(idxeta,ntmat,nmat).le.etaold)) then
                  write(*,*) '*ERROR reading *DAMAGE INITIATION:'
                  write(*,*) '       DM2 eta values must be strictly'
                  write(*,*) '       increasing.'
                  ier=1
                  return
                endif
                etaold=dmcon(idxeta,ntmat,nmat)
              enddo
            endif
            if(itype.eq.2) then
              if(dmcon(7,ntmat,nmat).le.dmcon(8,ntmat,nmat)) then
                write(*,*) '*ERROR reading *DAMAGE INITIATION'
                write(*,*) '       melt temperarature for material',
     &               matname(nmat)
                write(*,*)
     &               '       does not exceed transition temperature'
                ier=1
                return
              endif
            endif
          endif
!     
        enddo
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end

