!
! ======================================================================
!
       subroutine OPTDAILY     
       Use comrunmod
       Use comdynmod
       Use mod_varphy_coupl

       dimension JMONTH(13)
       data JMONTH /-1,31,59,90,120,151,181,212,243,273,304,334,365/

       do jm=1,12
         if (ntimeday.gt.JMONTH(jm).and.ntimeday.le.JMONTH(jm+1)) goto 1
       enddo

 1     nmonth=jm
       njour=ntimeday-JMONTH(nmonth)

       return
       end
