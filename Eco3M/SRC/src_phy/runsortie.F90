!
!
!***********************************************************************
!***********************************************************************
       subroutine RUNSORTIE
!
! Appele par main.f a chaque pas de temps: Moyenne et sauvegarde des
!   champs dynamiques et biologiques.
!   N'est effectif que tous les nsave pas de temps.
!
! ======================================================================
!
! ----------
! Sauvegarde
! ----------
!
! s'efectue tous les nsave pas de temps, ie nsa = nsave
!
       Use comrunmod
       Use comdynmod
       USE mod_varphy_coupl
!       include 'combio'
!
       nsa=nsa+1
       if (nsa.lt.nsave) goto 101
!
! Champs 1D
!
       call DYN1DUP
       call BIO1DUP
!
       nind=nind+1
!
! Champs 3D traceurs
!
       call DYN3DSAVE
! MB   call BIO3DSAVE
!
! Ecriture tous les nwrite*nsave pas de temps
!
       nwrite0=nwrite0+1
       if (mflagwrite.eq.1.and.nwrite0.eq.nwrite) then
         call WDYNPRINT
         call WBIOPRINT
         nwrite0=0
       endif
!
! Initialisation champs 3D traceurs apres ecriture
!
       call DYN3DUP
       call BIO3DUP
!
! Champs 3D tendances traceurs
!
       call DYNTEND
! MB       call BIOTEND
!
! Initialisation des indices apres sauvegarde
!
       nsa=0
       nind0=0
!
 101   continue
       return
       end
