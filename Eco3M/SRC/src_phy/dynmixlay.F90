!
!
!***********************************************************************
!***********************************************************************
       subroutine DYNMIXLAY    
       Use comrunmod
       Use comdynmod
!
!       do jk=2,nzt-1
!         if (AVT(jk).lt.AVT(jk-1)) goto 11
!       enddo
!       xdml=DEPW(2)
!       return
!
! 11    jk0=jk-1
!       do jk=jk0,nzt-1
!         if (AVT(jk).le.xlmkz) goto 13
!       enddo
! 13    if (jk.eq.jk0) then
!         xdml=DEPW(2)
!       else
!         xdml=DEPW(jk)
!       endif
!
        do jk=2,nzt-1
          if (EN(jk).lt.2.*emin) goto 11
        enddo
 11     xdml=DEPW(jk)
!
       return
       end
