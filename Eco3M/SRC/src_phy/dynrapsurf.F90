!
!
!***********************************************************************
!***********************************************************************
       subroutine DYNRAPSURF
       Use comrunmod
       Use comdynmod
!
       rindfutt=1.
       rindpast=1.
       rindfuts=1.
       rindpass=1.
!
       if (msurft.eq.0) then
         rapfutt=0.
         rappast=0.
         tsurfut=0.
         tsurpas=0.
         rindfutt=0.
         rindpast=0.
         goto 199
       endif
!
       zday=dayflux0+1.
!
! Temperature
       if (zday.lt.TTSURF(1)) goto 101
       if (zday.ge.TTSURF(ndsurft)) goto 102
       do jd=2,ndsurft
         if (zday.ge.TTSURF(jd-1).and.zday.lt.TTSURF(jd)) goto 111
       enddo
!
101    rapfutt=TTSURF(1)-zday
       tsurfut=TSURF(1)
       if (rapfutt.gt.RAPTPAS(1)) rindfutt=0.
       rappast=xperfluxt-TTSURF(ndsurft)+zday
       tsurpas=TSURF(ndsurft)
       if (rappast.gt.RAPTFUT(ndsurft)) rindpast=0.
       goto 199
!
102    rapfutt=xperfluxt-zday+TTSURF(1)
       tsurfut=TSURF(1)
       if (rapfutt.gt.RAPTPAS(1)) rindfutt=0.
       rappast=zday-TTSURF(ndsurft)
       tsurpas=TSURF(ndsurft)
       if (rappast.gt.RAPTFUT(ndsurft)) rindpast=0.
       goto 199
!
111    rapfutt=TTSURF(jd)-zday
       tsurfut=TSURF(jd)
       if (rapfutt.gt.RAPTPAS(jd)) rindfutt=0.
       rappast=zday-TTSURF(jd-1)
       tsurpas=TSURF(jd-1)
       if (rappast.gt.RAPTFUT(jd-1)) rindpast=0.
       goto 199
!
! Salinite
 199   continue
       if (msurfs.eq.0) then
         rapfuts=0.
         rappass=0.
         ssurfut=0.
         ssurpas=0.
         rindfuts=0.
         rindpass=0.
         goto 299
       endif
!
       if (zday.lt.TSSURF(1)) goto 121
       if (zday.ge.TSSURF(ndsurfs)) goto 122
       do jd=2,ndsurfs
         if (zday.ge.TSSURF(jd-1).and.zday.lt.TSSURF(jd)) goto 131
       enddo
!
121    rapfuts=TSSURF(1)-zday
       ssurfut=SSURF(1)
       if (rapfuts.gt.RAPSPAS(1)) rindfuts=0.
       rappass=xperfluxt-TSSURF(ndsurfs)+zday
       ssurpas=SSURF(ndsurfs)
       if (rappass.gt.RAPSFUT(ndsurfs)) rindpass=0.
       goto 299
!
122    rapfuts=xperfluxt-zday+TSSURF(1)
       ssurfut=SSURF(1)
       if (rapfuts.gt.RAPSPAS(1)) rindfuts=0.
       rappass=zday-TSSURF(ndsurfs)
       ssurpas=SSURF(ndsurfs)
       if (rappass.gt.RAPSFUT(ndsurfs)) rindpass=0.
       goto 299
!
131    rapfuts=TSSURF(jd)-zday
       ssurfut=SSURF(jd)
       if (rapfuts.gt.RAPSPAS(jd)) rindfuts=0.
       rappass=zday-TSSURF(jd-1)
       ssurpas=SSURF(jd-1)
       if (rappass.gt.RAPSFUT(jd-1)) rindpass=0.
       goto 299
!
 299   continue
!
       return
       end
