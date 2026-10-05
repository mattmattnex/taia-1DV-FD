  
!***********************************************************************
!***********************************************************************
      subroutine DYNINIT
      Use comrunmod
      Use comdynmod
      use mod_varphy_coupl
      
!     Ouverture des flux de forcage, de T  S initiaux et des rappels (profondeur et surface)     
      open(94, file=cdir(1:ndir) // '/inputs/' // cflux, &
          status='old', err=101)
      write(*,'(a)') cdir(1:ndir) // '/inputs/' // cflux

!     barrier.n --- modified 2015-09-15
!     barrier.n --- commented out (else, "file open in other unit" error)
!     open(994,file=trim(cdir)//'/'//cflux,status='old',err=101)
      open(91, file=cdir(1:ndir) // '/inputs/' // ctsinit,& 
          status='old', err=101)
      
      open(95, file=cdir(1:ndir) // '/inputs/' // crappel, &
          status='old', err=101)

!     Lecture champs initiaux     
      do jk=1,nzt
         read(91,*) TN(jk),SN(jk)
         TN(jk) = TN(jk) 
         TN(jk)=TN(jk)*TMASK(jk)
         SN(jk)=SN(jk)*TMASK(jk)
         TB(jk)=TN(jk)
         SB(jk)=SN(jk)
         TA(jk)=TN(jk)
         SA(jk)=SN(jk)
      enddo
      close(91)

!     Initialisation dynamique
      do jk=1,nzt
         UB(jk)=0.
         VB(jk)=0.
         UN(jk)=0.
         VN(jk)=0.
         UA(jk)=0.
         VA(jk)=0.
         EB(jk)=emin0*TMASK(jk)
         EN(jk)=emin0*TMASK(jk)
         EA(jk)=emin0*TMASK(jk)
      enddo

!     Initialisation sauvegarde
      do jk=1,nzt
         TMOY(jk)=0.
         SMOY(jk)=0.
         DENSITE(jk)=0.
         TKEMOY(jk)=0.
         DIFFT(jk)=0.
         TBN(jk)=0.
         TSH(jk)=0.
         TMLDF(jk)=0.
         TSHTKE(jk)=0.
         TBNTKE(jk)=0.
         TDFTKE(jk)=0.
         TDSTKE(jk)=0.
         TRAPT(jk)=0.
         TRAPS(jk)=0.
         TDIFT(jk)=0.
         TDIFS(jk)=0.
         TPENT(jk)=0.
      enddo
      do jt=1,jptemps
         DEPTHML(jt)=0.
         FLUXSOL(jt)=0.
         FLUXNSOL(jt)=0.
         FLUXADVT(jt)=0.
         FLUXSAL(jt)=0.
         FLUXADVS(jt)=0.
         FLUXVENTX(jt)=0.
         FLUXVENTY(jt)=0.
         FLUXVENT(jt)=0.
         FLUXSURT(jt)=0.
         FLUXSURS(jt)=0.
         FLUXSURTKE(jt)=0.
      enddo
      nind=1

!     Lecture rappels surfaceC     
      if (msurft.eq.1) then
         open(92, file=cdir(1:ndir) // '/inputs/' // crapsurft, &
             status='old',err=101)
         ndsurft=0
 201     continue
         read(92,*,end=202) ztt,z1,z2,zxx
         ndsurft=ndsurft+1
         TTSURF(ndsurft)=ztt
         TSURF(ndsurft)=zxx
         RAPTPAS(ndsurft)=z1
         RAPTFUT(ndsurft)=z2

         goto 201
 202     continue
         close(92)
      endif
     

      if (msurfs.eq.1) then
         open(92, file=cdir(1:ndir) // '/inputs/' // crapsurfs, &
      status='old',err=101)
         ndsurfs=0
 211     continue
         read(92,*,end=212) ztt,z1,z2,zxx
         ndsurfs=ndsurfs+1
         TSSURF(ndsurfs)=ztt
         SSURF(ndsurfs)=zxx
         RAPSPAS(ndsurfs)=z1
         RAPSFUT(ndsurfs)=z2
         goto 211
 212     continue
         close(92)
      endif

!     Calcul des moyennes journalieres flux de chaleur > 0 (assimilation)     
      call DYNSOLDAY

!     Calage temporel initial des flux et rappels     
      call DYNFLUXINIT

!     Ouverture des fichiers sortie    
#if defined SAVE_PHY
      call DYNOPENFIC
#endif
      return
     
 101  write(*,*) 'pb avec fichier dyn'
      mflagtke=1
     
      return
      end
