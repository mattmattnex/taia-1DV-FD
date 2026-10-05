!
! -----------
! dynamique.f
! -----------
! 27/02/97: Version 1 provenant de biotke.f (lm)
!
! Appele a chaque pas de temps par main.f
!
! Programme en quatre parties:
!   1. Diagnostique sur stratification, longueur de melange et diffusion
!              verticale.
!   2. Equation pronostique de TKE
!   3. Equation pronostique de U et V (couche d'Ekman)
!   4. Equation pronostique pour T et S
!
! ======================================================================
!
       subroutine DYNAMIQUE
       Use comrunmod
       Use comdynmod
!
! 1111111111111111111111111111111111111111
!
! Calcul de la diffusion verticale
!
       call DYNKZ
!
! Calcul de la profondeur de la couche melangee
!
       call DYNMIXLAY
       DEPTHML(nind)=DEPTHML(nind)+xdml
!
! 1111111111111111111111111111111111111111
!
! 2222222222222222222222222222222222222222
!
! Lecture des forcages atmospheriques
!
       call DYNFLUXSUR
       FLUXVENTX(nind)=FLUXVENTX(nind)+windx
       FLUXVENTY(nind)=FLUXVENTY(nind)+windy
       FLUXVENT(nind)=FLUXVENT(nind)+wind
!
! Integration de TKE
!
       call DYNTKE
!
! 2222222222222222222222222222222222222222
!
! 3333333333333333333333333333333333333333
!
! Integration des vitesses
!
       call DYNVIT
!
! 3333333333333333333333333333333333333333
!
! 4444444444444444444444444444444444444444
!
! Integration temperature, salinite
!
       call DYNTS
!
! 4444444444444444444444444444444444444444
       
       return
       end
