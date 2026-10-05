!
!
!***********************************************************************
!***********************************************************************
       subroutine BIO3DSAVE
!
! ------------------------------------------------
! Sauvegarde sur nunbio des champs biogeochimiques
! ------------------------------------------------
!
!       include 'comrun'
!       include 'comdyn'
!       include 'combio'
       Use comrunmod
       Use comdynmod
       Use mod_varphy_coupl
!
! Moyenne
!
       do jtr=1,nbrbio
         do jk=1,nzt
           BIOAVER(jk,jtr)=BIOAVER(jk,jtr)/nsave
         enddo
       enddo
!CCCC!!!!!!!!!!!!!!!!!!!!!!!!
! PUR et PAR: moyenne sur periode eclairee de la journee
!
        do jk=1,nzt
          BIOAVER(jk,13)=BIOAVER(jk,13)*nsave/nind0
          BIOAVER(jk,14)=BIOAVER(jk,14)*nsave/nind0
        enddo
!CCCC!!!!!!!!!!!!!!!!!!!!!!!!
!
! Sauvegarde sur nunbio
!
       do jtr=1,nbrbio
         if (MSTOCK(jtr).eq.1) then
           write(nunbio,*) (BIOAVER(jk,jtr),jk=nzt,1,-1)
         endif
       enddo
!
       return
       end
