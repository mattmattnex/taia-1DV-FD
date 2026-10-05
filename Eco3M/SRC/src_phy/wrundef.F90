!
!=======================================================================
!
       subroutine WRUNDEF
       Use comrunmod
       Use comdynmod
!
!       if (mflagrun.eq.1) then
!         write(99,1001) nbrbio,nbrprono
!         print*, ' ==================='
!         print*, ' STOP : mflagrun = 1'
!         print*, ' ==================='
!         stop
!       endif
! 1001  format(/,2x,'Nbr de variables biogeochimiques = ',i3,
!     &   'avec  nbr de variables pronostiques = ',i3)
!
       if (mflagrun.eq.11) then
         write(99,1002) nzt,jpzt       
         print*, ' ===================='
         print*, ' STOP : mflagrun = 11'
         print*, ' ===================='
         stop
       endif
 1002  format(/,2x,'Nbr de niveaux verticaux = ',i3,&
                    ' avec parameter jpzt = ',i3)
!
       if (mflagrun.eq.12) then
         write(99,1003) nbrprono,nbrdiag,nbrbio,jptract 
         print*, ' ===================='
         print*, ' STOP : mflagrun = 12'
         print*, ' ===================='
         stop
       endif
 1003  format(/,2x,'Nbr de traceurs biologiques prono = ',i3,5x,&
       'diag = ',i3,/,2x,'===> total = ',i3,&
          ' avec parameter jptract = ',i3)
!
       if (mflagrun.eq.13) then
         write(99,1004) nbrflux,jpfluxt
         print*, ' ===================='
         print*, ' STOP : mflagrun = 13'
         print*, ' ===================='
         stop
       endif
 1004  format(/,2x,'Nbr de flux biologiques = ',i3,&
                    ' avec parameter jpfluxt = ',i3)
!
       if (mflagrun.eq.2) then
         write(99,1005) dt,timesave,nsave
         print*, ' ==================='
         print*, ' STOP : mflagrun = 2'
         print*, ' ==================='
         stop
       endif
 1005  format(/,2x,'Dt = ',f6.0,' mn',3x,'Sauvegarde = ',f6.2,' h',&
         3x,' ===> Pas de temps de sauvegarde = ',i2)
!
       if (mflagrun.eq.4) then
         jsd=int(tfin0*day/timesave+0.001)
         write(99,1008) tfin0,timesave,jsd,jptemps
         print*, ' ==================='
         print*, ' STOP : mflagrun = 4'
         print*, ' ==================='
         stop
       endif
 1008  format(/,2x,'Dureee du run = ',f8.2,' jours',3x,  &
         'Sauvegarde = ',f8.2,' heures',/,2x,            &
         'Longueur vecteurs = ',i8,                      &
         '  >>> jptemps = ',i8,/,2x,                     &
         'Parameter jptemps dans comrun a modifier')
!
       return
       end
