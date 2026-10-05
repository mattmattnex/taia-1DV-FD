!
!
!***********************************************************************
!***********************************************************************
       subroutine BIOSURFACE
              
       Use comrunmod
       Use comdynmod
       USE mod_varphy_coupl
       IMPLICIT None
       Integer:: jtr
       Real(8)::netflux
!
       do jtr=1,jptract
         BIOAIR(jtr)=0.
       enddo
!
       call ATMOCEOXY(jtr,netflux)
       BIOAIR(jtr) = netflux ! en mol O2/m2/s
       BIOAIR(jtr) = BIOAIR(jtr) / DE3T(1) ! en mol O2/m3/s
!       call BIOFLUXPCO2
!
       return
       end
