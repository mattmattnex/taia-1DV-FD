!
!
!***********************************************************************
!***********************************************************************
       subroutine DYN3DUP 
!
! ------------------------------------------------------
! Initialisation des moyennes apres sauvegarde effectuee
! ------------------------------------------------------
! Traceurs dynamiques
!
       Use comrunmod
       Use comdynmod
!       include 'combio'
!
       do jk=1,nzt
         TMOY(jk)=0.
       enddo
       do jk=1,nzt
         SMOY(jk)=0.
         DENSITE(jk)=0.
       enddo
       do jk=1,nzt
         TKEMOY(jk)=0.
       enddo
       do jk=1,nzt
         DIFFT(jk)=0.
       enddo
!
       return
       end
