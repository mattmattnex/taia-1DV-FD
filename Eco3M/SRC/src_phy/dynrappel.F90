!
!
!***********************************************************************
!***********************************************************************
       subroutine DYNRAPPEL
       Use comrunmod
       Use comdynmod
!
       if (dayflux0.ge.dayts0.and.dayflux0.lt.dayts1) return
 102   continue
       read(95,*,end=101) dayts0,dayts1
       do jk=1,nzt
         read(95,*) XTS(jk),TRP(jk),SRP(jk)
         TRP(jk) = TRP(jk) 
       enddo
       return
!
 101   rewind(95)
       goto 102
!
       end
