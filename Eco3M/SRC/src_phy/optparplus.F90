!
! ======================================================================
!
       subroutine OPTPARPLUS
       Use comrunmod
       Use comdynmod
       Use mod_varphy_coupl

! Si changement de journee solaire, calcul du calendrier, et des
!   donnees astronomiques (eclairement)
       if (timeminu.gt.-dt/2.and.timeminu.lt.dt/2.) then
         call OPTDAILY
         call OPTASTRO
       endif

! Calcul de PAR0
       ztime=timeminu/hourm
       zdy0=(day-daylength)/2.
       zdy1=(day+daylength)/2.
       if (ztime.lt.zdy0.or.ztime.gt.zdy1) then
         xpar0p=0.
         xpar0m=0.
         do jk=1,nzt
           VPAR(jk)=0.
           VPUR(jk)=0.
         enddo          
         fparml=0.
         fparze=0.
         fpurml=0.
         fpurze=0.
       else
         call OPTPAR0
       endif

       return
       end
