!
!
! **********************************************************************
! **********************************************************************
       subroutine BIODYN 
!
! Appele par biologie.f a chaque passage: calcul des tendances dynamiques
! des variables biologiques et update des variables
!
! ======================================================================
!
       Use comrunmod
       Use comdynmod
       USE mod_varphy_coupl


       do jtr=1,jptract  ! barrier.n -> replacement of nbrprono by jptract
         workb=BIOAIR(jtr)*dts
         do jz=1,nzt
           WORK(jz)=TENEUR(jz,jtr)
           WORK1(jz)=TENEUR(jz,jtr)
         enddo
         WORK(1)=WORK(1)+workb
         WORK1(1)=WORK1(1)+workb

         call BIOINVERSE

         do jz=1,nzt
           TENEUR(jz,jtr)=WORK(jz)
         enddo
         FLUXAIR(nind,jtr)=FLUXAIR(nind,jtr)+BIOAIR(jtr)
         do jz=1,nzt
           DIFFBIO(jz,jtr)=DIFFBIO(jz,jtr)+(WORK(jz)-WORK1(jz))/dts
         enddo
       enddo

!
! Diagnostique sur rapport !/Chl
!
! MB       if (nchlapro.eq.1) then
! MB         do jz=1,nzt
! MB           if (TENEUR(jz,niparchl).gt.teneurmin) then
! MB             TENEUR(jz,nicchl)=carbmass*REDCN(niphy)*TENEUR(jz,niphy)/
! MB     &            TENEUR(jz,niparchl)
! MB           else
! MB             TENEUR(jz,nicchl)=cchlmin
! MB           endif
! MB         enddo
! MB       endif
!
! Stockage en carbone
!
!       VPCO2(nind)=VPCO2(nind)+xpco2
!       VITTRANS(nind)=VITTRANS(nind)+vitgaz
!
       return
       end
