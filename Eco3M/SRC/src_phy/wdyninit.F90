!
!=======================================================================
!
       subroutine WDYNINIT
       Use comrunmod
       Use comdynmod
!
       if (mflagtke.eq.1) then
         write(99,1011) cflux,ctsinit,crappel,msurf,crapsurf
         print*, ' ==================='
         print*, ' STOP : mflagtke = 1'
         print*, ' ==================='
         stop
       endif
 1011  format(/,2x,'Fichier forcages atmospheriques:',/,2x,a80,/,     & 
         2x,'Fichier profils initiaux T et S:',/,2x,a80,/,            &
         2x,'Fichier rappel profils T et S:',/,2x,a80,/,              &
         2x,'msurf= ',i2,2x,'Fichier rappel surface T et S:',/,2x,a80,&
         /,2x,'UN DE CES FICHIERS AU MOINS N EXISTE PAS')
!
       if (mflagtke.eq.2) then
         write(99,1007) jpant,njourflux
         print*, ' ==================='
         print*, ' STOP : mflagtke = 2'
         print*, ' ==================='
         stop
       endif
 1007  format(/,2x,'Parameter comdyn jpant = ',i4,5x,&
        'Input tkeparam.i njourflux = ',i4,/,2x,     &
         'Probleme de dimensions de tableaux')
!
       if (mflagflux.eq.1) then
         write(99,1012) tdebut0,xperflux0,xperfluxt,xperflux0+xperfluxt
         print*, ' ==================='
         print*, ' STOP : mflagflux = 1'
         print*, ' ==================='
         stop
       endif
 1012  format(/,2x,'Debut run = ',f6.2,/,2x,'Debut forcage = ',f6.2,    &
           2x,/,                                                        &
         'Temps forcage total = ',f6.2,3x,'Fin forcage = ',f6.2,        &
          ' jours',/,                                                   &
         2x,'Debut run DOIT ETRE ENTRE debut forcages et fin forcages'  &
            ,/,                                                         &
         2x,'Voir les entrees dans runparam.i et le fichier forcages')
!
       if (mflagflux.eq.2) then
         write(99,1013) tdebut0,timeday,crappel
         print*, ' ==================='
         print*, ' STOP : mflagflux = 2'
         print*, ' ==================='
         stop
       endif
 1013  format(/,2x,'Debut run = ',f6.2,' Jours, ie jour dans annee = ', &
         f6.2,/,2x,                                                     &
         'Les profils de rappel en T et S doivent avoir le premier',    &
         ' intervalle de temps commencant a 1',/,2x,                    &
         'et le dernier a 365: voir fichier :',/,2x,a80)
!
       if (mflagflux.eq.3) then
         write(99,1014) tdebut0,xperflux0,xperflux,xperfluxt,nday0,nday1
         print*, ' ==================='
         print*, ' STOP : mflagflux = 3'
         print*, ' ==================='
         stop
       endif
 1014  format(/,2x,'Debut run = ',f6.2, ' jours - Debut forcage = ',  &
        f6.2,' jours',/,2x,'periode forcage = ',f6.2,                 &
        ' heures - duree forcage = ',f6.2,' jours',/,2x,              &
         ' ===> nday0 = ',i5,2x,'nday1 = ',i5,/,2x,                   &
          'Il y a un probleme initialisation flux dans la definition',&
          ' du jour')
!      
       write(99,1001) ctsinit
 1001  format(/,2x,'PROFILS INITIAUX T et S: Fichier',/,2x,a80)
       if (mflagwrite.eq.1) then
         do jk=1,nzt-1,2
           write(99,1002) jk,DEPT(jk),TN(jk),SN(jk)
         enddo
       endif
 1002  format(i5,f8.1,2f8.3)
!
       return
       end
