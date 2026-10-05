!
!=======================================================================
!
       subroutine WRUNPRINT
       Use comrunmod
       Use comdynmod
!       include 'combio'
!
       write(99,1000) 
 1000  format(/,2x,40('='),/,2x, &
         'DEFINITION DU RUN: directory de travail')
       write(99,'(a)') cdir(1:ndir)
       write(99,1098) nroot,croot
 1098  format(2x,'=====> Fichiers binaire: character length = ',i3,&
        /,2x,'ROOT = ',a60)
       write(99,1099) ctit,csbtit
 1099  format(2(/,2x,a80))
       write(99,1001) dt,tdebut0,tfin0,npurday
 1001  format(/,2x,'Pas de temps = ',f5.1,' mn',/,2x,             & 
          'Instant initial = ',f6.2,' jour',5x,'Duree du run = ', &
          f8.2,' jours',/,2x,'Nbr pas de temps / jour = ',i4)
       write(99,1002) nzt,nbrbio,nbrflux
 1002  format(2x,'Nbr de niveaux = ',i3,5x,', de traceurs bio = ',i3,5x,&
          'de flux bio = ',i3)
       write(99,1003) timesave,nsave,nwrite,&
         nwrite*timesave/day
 1003  format(2x,'Periode de sauvegarde = ',f6.2,2x, &
         'Nbr de pas de temps sauvegarde = ',i4,     &
         /,2x,'Ecriture: nbr pas de temps save = ',  &
          i5,5x,'ie ',f8.1,' jours')
       write(99,1011) (CNMTRA(jtr),jtr=1,nbrbio)
 1011  format(/,2x,'Variables biogeochimiques',10(/,2x,10a7))
       write(99,1014) (MSTOCK(jtr),jtr=1,nbrbio)
       write(99,1111) nbrprono
 1111  format(2x,'Nbr de traceurs pronostiques = ',i3)
       write(99,1112) (CNMTRA(jtr),jtr=1,nbrprono)
 1112  format(5(2x,12a7,/))
       write(99,1012) (MBIODIF(jtr),jtr=1,nbrprono)
! barrier.n --- modified 2015-09-15
! format changed from 10i7 to f10.7
 1012  format(2x,'Sauvegarde diffusion',10(/,2x,f10.7))
       write(99,1031) mbiosed
! barrier.n --- modified 2015-09-15
! format changed from i2 to f5.1
 1031  format(2x,'Sauvegarde de la sedimentation: flag = ',f5.1)
       write(99,1013) (CNMFLX(jtr),jtr=1,nbrflux)
 1013  format(/,2x,'Flux biogeochimiques',10(/,2x,10a7))
       write(99,1014) (MFLUX(jtr),jtr=1,nbrflux)
! barrier.n --- modified 2015-09-15
! format changed from 10i7 to f10.7
 1014  format(2x,'Sauvegarde',10(/,2x,f10.7))
       write(99,1004) cgrille
 1004  format(/,2x,'Fichier grille:',/,2x,a80)
!
       return
       end
