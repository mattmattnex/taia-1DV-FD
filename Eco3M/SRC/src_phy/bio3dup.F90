!
!
!***********************************************************************
!***********************************************************************
       subroutine BIO3DUP 
!
! ------------------------------------------------------
! Initialisation des moyennes apres sauvegarde effectuee
! ------------------------------------------------------
! Traceurs biologiques
!
!       include 'comrun'
!       include 'comdyn'
!       include 'combio'
       Use comrunmod
       Use comdynmod
       Use mod_varphy_coupl

!
       do jtr=1,nbrbio
         do jk=1,nzt
           BIOAVER(jk,jtr)=0.
         enddo
       enddo
!
       return
       end
