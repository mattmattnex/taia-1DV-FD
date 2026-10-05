!
! ----------
! biologie.f
! ----------
! 27/02/97: Version 1 provenant de biotke.f (lm)
!
! Appele a chaque pas de temps par main.f.
!
! Programme en quatre parties:
!   1. Calcul des profils optiques PAR et PUR, ainsi que des moyennes
!           dans couche melangee, et couche euphotique
!            (subroutines dans optique.f)
!           PARPLUS = calcul de PAR0+ par ciel clair
!           PURPROF = calculs des profils de PAR et de PUR
!   2. Recherche et calage temporel des profils verticaux de rappels pour
!          les variables biogeochimiques (dissoutes) = BIORAPPEL
!   3. Calcul des flux a l'interface air - mer = BIOSURFACE
!   4. Estimation des termes non conservatifs biogeochimiques =
!                BIOMASSE
!
! ======================================================================

       subroutine BIOLOGIE
       Use comrunmod
       Use comdynmod
       Use mod_eco3m
       USE mod_varphy_coupl

! Calcul du PAR0+ et PAR0-
!       call OPTPARPLUS
!	write(*,*) "biologie xpar0p ",xpar0p

! Calcul du profil vertical de PAR et de PUR    
!       npurday0 = npurday0 + 1
!       if (xpar0p.gt.0.) then
!         call OPTPURPROF
!         DEPTHZE(nind)=DEPTHZE(nind)+seuph
!         XNEBUL(nind)=XNEBUL(nind)+trpar
!         XPAR0PLUS(nind)=XPAR0PLUS(nind)+xpar0p
!         nind0=nind0+1
!         npurd0=npurd0+1
!         fzede0=fzede0+seuph
!         fmlde0=fmlde0+xdml
!         fmlpur0=fmlpur0+fpurml
!        fzepur0=fzepur0+fpurze
!       endif

! Flux a la surface
       call BIOSURFACE

! Calcul des flux biogeochimiques
!       call BIOMASSE
       call eco3m_conc_update
       
       if (npurday0.eq.npurday) then
         npurday0=0
         npurd0=0
       endif

       return
       end
