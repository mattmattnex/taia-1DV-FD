!     
!     
!***********************************************************************
!***********************************************************************
      subroutine BIODEF 
      Use comrunmod
      Use comdynmod
      USE mod_eco3m
      use mod_varphy_coupl
      
      Integer::ivar,ii
      character(200) :: clpipo
      real(8) :: vitsedfile

!     Lecture des parametres du modele biologique
      open(90,file=cdir(1:ndir)//'/inputs/bioparam.i',status='old')

!     Photo - phyto
      read(90,'(a)') clpipo
      read(90,*) teneurmin
      read(90,'(a)') clpipo
      read(90,*) xslopet
      read(90,'(a)') clpipo
      read(90,*) xzelight
      read(90,'(a)') clpipo
      read(90,*) flimlm0
      read(90,'(a)') clpipo
      read(90,*) chlaphea
      read(90,'(a)') clpipo
      read(90,*) cchlmin
      read(90,'(a)') clpipo
      read(90,*) cchlmax
      read(90,'(a)') clpipo
      read(90,*) xpurmx
      read(90,'(a)') clpipo
      read(90,*) ftimenut0
      read(90,'(a)') clpipo
      read(90,*) ftimenut1

!     Sedimentation
      read(90,'(a)') clpipo
      write(*,*) trim(clpipo)
      VITSEDfile = 0.0
      vitsed0 = 0.0
      read(90, *) vitsedfile
      write(*,*) 'vitsedfile =',vitsedfile
      do ivar = 1,nbvar
         if (VAR(ivar)%comp == 'mop' .or. VAR(ivar)%comp == 'det') then 
            if (VAR(ivar)%scomp == 'MOPS' .or. VAR(ivar)%scomp=='DETS') then 
                   vitsed0(ivar) = VITSEDfile
            elseif (VAR(ivar)%scomp=='MOPL'.or.VAR(ivar)%scomp=='DETL') then
                   vitsed0(ivar) = VITSEDfile * 25.d0  !100
!                   vitsed0(ivar) = VITSEDfile * 50.d0  !100
!                   vitsed0(ivar) = VITSEDfile * 12.5d0  !100
            endif
!                   vitsed0(ivar) = VITSEDfile
         endif
      enddo
!      #IFDEF DMIGR
!	    print*,"THE MIGR KEY WORKED!"
!      #ENDIF

!     Passage de m/jour en SI (m/sec)
! MB correction bug      do jtr=1,nbrprono
      do jtr=1,nbvar
         VITSED(jtr)=VITSED0(jtr)/(day*hour)
         write(*,*) 'jtr,vitsed',jtr,vitsed(jtr)
      enddo

!     Rappels NO3 DOR - Initialisation     
      read(90,'(a)') clpipo
      read(90,'(a)') coptique
      read(90,'(a)') clpipo
      read(90,*) mrapbio
      read(90,'(a)') clpipo
      read(90,'(a)') crapbio
      read(90,'(a)') clpipo
      read(90,'(a)')  cbioinit
!     
      close(90)
!     
!     Clha: variabile pronostique ?
!     
      nchlapro=0
      do jt=1,nbrprono
         if (jt.eq.nicchl) nchlapro=1
      enddo
!     
!     Eta = uptake ammonium / uptake DOM labile par bacteries
!     
!     MB       eta=REDCN(nidom)/REDCN(nibac)-1.
!     
      return
      end
!***********************************************************************
#if defined MIG_MZOO      

      subroutine CALC_VMIG
      Use comrunmod
      Use comdynmod
      USE mod_eco3m
      use mod_varphy_coupl

      Integer::ivar,ii

      Do ivar=1,nbvar
         if (VAR(ivar)%comp == 'zoo') then
            if (VAR(ivar)%scomp == 'COP') then
! par convention une vitesse negative va vers la surface et positive
! vers le fond                    
                   vitsed0(ivar) = -17.d0 *  cos((tps/3600.0d0 -1.0d0)*&
                   2 *pi / 24.0d0)
                   VITSED(ivar)=VITSED0(ivar)/(day*hour)
           endif
         endif
      Enddo
      End subroutine CALC_VMIG
#endif /* MIG_MZOO*/
!***********************************************************************
