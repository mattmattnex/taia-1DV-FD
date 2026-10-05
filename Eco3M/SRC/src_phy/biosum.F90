!
!
!***********************************************************************
!***********************************************************************
       subroutine BIOSUM
!       include 'comrun'
!       include 'comdyn'
!       include 'combio'
       Use comrunmod
       Use comdynmod
       Use mod_varphy_coupl
!
       do jk=1,nzt
         do jtr=1,nbrbio
           BIOAVER(jk,jtr)=BIOAVER(jk,jtr)+TENEUR(jk,jtr)
         enddo
       enddo
!
       return
       end
