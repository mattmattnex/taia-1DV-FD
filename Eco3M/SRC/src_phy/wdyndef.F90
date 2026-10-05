!
!=======================================================================
!
       subroutine WDYNDEF
       Use comrunmod
       Use comdynmod
!
       write(99,1001) xlat,emin,emin0,avmb,avtb,xpdl,fave,drag
 1001  format(/,2x,'PARAMETRES TKE',/,2x,'Latitude = ',f8.3,/,2x,    &
         'TKE min = ',1pe10.2,3x,'min surface = ',e10.2,/,2x,        &
         'Diffusion background: vitesse = ',e10.2,5x,'traceurs = ',  &
          e10.2,/,2x,                                                &
          'Prandtl = ',0pf5.2,2x,'fave = ',f5.2,5x,'drag = ',f9.5)
       write(99,1002) crappel
 1002  format(/,2x,'Fichier rappel profils T et S:',/,2x,a80)
       if (msurf.eq.1) then
         write(99,1003) crapsurf
         write(99,1004) rapmax,rapprop,alfcte,alfmin,alfmax
       endif
 1003  format(/,2x,'Fichier rappel surface T et S:',/,2x,a80)
 1004  format(2x,'rapmax = ',f6.2,' jours',5x,                        & 
         'Intensification temps de rappel = ',f6.2,/,2x,'alfcte = ', &
          f6.2,3x,'alfmin = ',f6.2,3x,'alfmax = ',f6.2)
!
       return
       end
