!
!
!***********************************************************************
!***********************************************************************
       subroutine RUNFINAL 
!
! Appele par main.f a la fin de la simulation: sauvegarde des champs 1d
!   d'evolution temporelle et ecruture finale de statistiques sur les
!   valeurs sous seuil des variables biologiques.
!
! ======================================================================
!
! ----------------------------------------------
! Sauvegarde et ecriture finales a la fin du run
! ----------------------------------------------
!
       Use comrunmod
       Use comdynmod
!       include 'combio'
!
! Sauvegarde des champs 2D dynamiques sur nundyn1d
!
       write(nundyn1d,1001) (float(jt),jt=1,nind-1)
       write(nundyn1d,1001) (-DEPTHML(jt),jt=1,nind-1)
       write(nundyn1d,1002) (FLUXVENTX(jt),jt=1,nind-1)
       write(nundyn1d,1002) (FLUXVENTY(jt),jt=1,nind-1)
       write(nundyn1d,1002) (FLUXVENT(jt),jt=1,nind-1)
       write(nundyn1d,1003) (FLUXSOL(jt),jt=1,nind-1)
       write(nundyn1d,1003) (FLUXNSOL(jt),jt=1,nind-1)
       write(nundyn1d,1003) (FLUXADVT(jt),jt=1,nind-1)
       write(nundyn1d,1003) (FLUXSURT(jt),jt=1,nind-1)
       write(nundyn1d,1003) (FLUXSAL(jt),jt=1,nind-1)
       write(nundyn1d,1003) (FLUXADVS(jt),jt=1,nind-1)
       write(nundyn1d,1003) (FLUXSURS(jt),jt=1,nind-1)
       write(nundyn1d,1003) (FLUXSURTKE(jt),jt=1,nind-1)
 1001  format(0p,10f12.1)
 1002  format(0p,10f12.3)
 1003  format(1p,10e12.4)
!
! Sauvegarde des champs 2D biologiques sur nunbio1d
!
! Lignes commentees par MB
! MB       write(nunbio1d,1001) (float(jt),jt=1,nind-1)
! MB       write(nunbio1d,1001) (-DEPTHZE(jt),jt=1,nind-1)
! MB       write(nunbio1d,1002) (XNEBUL(jt),jt=1,nind-1)
! MB       write(nunbio1d,1003) (XPAR0PLUS(jt),jt=1,nind-1)
!       do jtr=1,nbrbio
!         if (MSTOCK(jtr).eq.1) then
!           write(89,1003) (FLUXAIR(jt,jtr),jt=1,nind-1)
!         endif
!       enddo
!
! Ecriture des statistiques finales sur les parametres biologiques < 0
! durant l'integration
!
       call WRUNFIN
!
       return
       end
