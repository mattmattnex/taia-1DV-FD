!
!
!***********************************************************************
!***********************************************************************
       subroutine DYNTS
!
! Integration de Temperature et Salinite
!
       Use comrunmod
       Use comdynmod
       dimension ZWORK(jpzt)
!
! Initialisation - Rappel sur colonne
       call DYNRAPPEL
       do jk=1,nzt-1
         zexp1=rabs*EXP(-DEPW(jk)/xsi1)+(1.-rabs)*EXP(-DEPW(jk)/xsi2)
         zexp2=rabs*EXP(-DEPW(jk+1)/xsi1)+(1.-rabs)*EXP(-DEPW(jk+1)/&
           xsi2)
         zpen=qsr*(zexp1-zexp2)
         TPENT(jk)=TPENT(jk)+zpen/E3T(jk)
         TA(jk)=zpen/(raubase*xcp)/E3T(jk)
         zrapt=xcp*E3T(jk)*raubase*XTS(jk)*(TRP(jk)-TN(jk))
         TRAPT(jk)=TRAPT(jk)+zrapt/E3T(jk)
         zraps=XTS(jk)*(SRP(jk)-SN(jk))
         TRAPS(jk)=TRAPS(jk)+zraps
         TA(jk)=TA(jk)+XTS(jk)*(TRP(jk)-TN(jk))
         SA(jk)=XTS(jk)*(SRP(jk)-SN(jk))
         FLUXADVT(nind)=FLUXADVT(nind)+zrapt
         FLUXADVS(nind)=FLUXADVS(nind)+zraps
       enddo
       FLUXSOL(nind)=FLUXSOL(nind)+qsr
       FLUXNSOL(nind)=FLUXNSOL(nind)+q
       FLUXSAL(nind)=FLUXSAL(nind)+SN(1)*ep/(tcs*E3T(1))
!
! Rappel en surface 
!       if (rapfut.eq.0.) goto 501
       call DYNRAPSURF
! Rappel donnee future
       zrap0=AMAX1(rapmax,rapfutt)
       zrap=rindfutt*rapprop/(zrap0*day*hour)
       zrapt=xcp*E3T(1)*raubase*zrap*(tsurfut-TN(1))
       TRAPT(1)=TRAPT(1)+zrapt/E3T(1)
       FLUXSURT(nind)=FLUXSURT(nind)+zrapt
       TA(1)=TA(1)+zrap*(tsurfut-TN(1))
       zrap0=AMAX1(rapmax,rapfuts)
       zrap=rindfuts*rapprop/(zrap0*day*hour)
       zraps=zrap*(ssurfut-SN(1))
       TRAPS(1)=TRAPS(1)+zraps
       FLUXSURS(nind)=FLUXSURS(nind)+zraps
       SA(1)=SA(1)+zrap*(ssurfut-SN(1))
! Rappel donnee passee
! 501   if (rappas.eq.0.) goto 502
       zrap0=AMAX1(rapmax,rappast)
       zrap=rindpast*rapprop/(zrap0*day*hour)
       zrapt=xcp*E3T(1)*raubase*zrap*(tsurpas-TN(1))
       TRAPT(1)=TRAPT(1)+zrapt/E3T(1)
       FLUXSURT(nind)=FLUXSURT(nind)+zrapt
       TA(1)=TA(1)+zrap*(tsurpas-TN(1))
       zrap0=AMAX1(rapmax,rappass)
       zrap=rindpass*rapprop/(zrap0*day*hour)
       zraps=zrap*(ssurpas-SN(1))
       TRAPS(1)=TRAPS(1)+zraps
       FLUXSURS(nind)=FLUXSURS(nind)+zraps
       SA(1)=SA(1)+zrap*(ssurpas-SN(1))
!502    continue
!
! Initialisation T
       do jk=1,nzt-1
         DZW(jk)=-dtsd*DAVT(jk)/DE3T(jk)/DE3W(jk)
         DZY(jk)=-dtsd*DAVT(jk+1)/DE3T(jk)/DE3W(jk+1)
         DZX(jk)=dun-DZW(jk)-DZY(jk)
         DZZ(jk)=DBLE(TN(jk))+dtsd*DBLE(TA(jk))
       enddo
! Surface
       DZW(1)=0.
       DZX(1)=dun-DZY(1)
! MB: dans dynfluxsur.f:   q=(npasflux0*zq1+(npasflux-npasflux0)*zq0)/npasflux
       DZZ(1)=DZZ(1)+dtsd/DE3T(1)*DBLE(q)/DBLE(raubase)/DBLE(xcp)
       do jk=1,nzt-1
         ZWORK(jk)=SNGL(DZZ(jk))
       enddo
! Fond
       DZY(nzt-1)=0.
       DZX(nzt-1)=dun-DZW(nzt-1)
! Inversion
       DZR(1)=DZX(1)
       do jk=2,nzt-1
         DZR(jk)=DZX(jk)-DZW(jk)*DZY(jk-1)/DZR(jk-1)
       enddo
       DZX(1)=DZZ(1)
       do jk=2,nzt-1
         DZX(jk)=DZZ(jk)-DZW(jk)/DZR(jk-1)*DZX(jk-1)
       enddo
       DZW(nzt-1)=DZX(nzt-1)/DZR(nzt-1)
       do jk=nzt-2,1,-1
         DZW(jk)=(DZX(jk)-DZY(jk)*DZW(jk+1))/DZR(jk)
       enddo
! Champ after
       do jk=1,nzt-1
         TA(jk)=SNGL(DZW(jk))
         ZWORK(jk)=TA(jk)-ZWORK(jk)
         TDIFT(jk)=TDIFT(jk)+xcp*raubase*ZWORK(jk)/dts
       enddo
! Initialisation S
       do jk=1,nzt-1
         DZW(jk)=-dtsd*DAVT(jk)/DE3T(jk)/DE3W(jk)
         DZY(jk)=-dtsd*DAVT(jk+1)/DE3T(jk)/DE3W(jk+1)
         DZX(jk)=dun-DZW(jk)-DZY(jk)
         DZZ(jk)=DBLE(SN(jk))+dtsd*DBLE(SA(jk))
       enddo
! Surface
       DZW(1)=0.
       DZX(1)=dun-DZY(1)
       DZZ(1)=DZZ(1)+dtsd/DE3T(1)*DBLE(SN(1))*DBLE(ep)/dtcs
       do jk=1,nzt-1
         ZWORK(jk)=SNGL(DZZ(jk))
       enddo
! Fond
       DZY(nzt-1)=0.
       DZX(nzt-1)=dun-DZW(nzt-1)
! Inversion
       DZR(1)=DZX(1)
       do jk=2,nzt-1
         DZR(jk)=DZX(jk)-DZW(jk)*DZY(jk-1)/DZR(jk-1)
       enddo
       DZX(1)=DZZ(1)
       do jk=2,nzt-1
         DZX(jk)=DZZ(jk)-DZW(jk)/DZR(jk-1)*DZX(jk-1)
       enddo
       DZW(nzt-1)=DZX(nzt-1)/DZR(nzt-1)
       do jk=nzt-2,1,-1
         DZW(jk)=(DZX(jk)-DZY(jk)*DZW(jk+1))/DZR(jk)
       enddo
! Champ after
       do jk=1,nzt-1
         SA(jk)=SNGL(DZW(jk))
         ZWORK(jk)=SA(jk)-ZWORK(jk)
         TDIFS(jk)=TDIFS(jk)+ZWORK(jk)/dts
       enddo
! Filtre temporel et swap des tableaux
       do jk=1,nzt
         TN(jk)=TA(jk)*TMASK(jk)
         SN(jk)=SA(jk)*TMASK(jk)
       enddo
!       if (mt.eq.1) then
         do jk=1,nzt
           TMOY(jk)=TMOY(jk)+TN(jk)
         enddo
!       endif
!       if (ms.eq.1) then
         do jk=1,nzt
           SMOY(jk)=SMOY(jk)+SN(jk)
           DENSITE(jk)=DENSITE(jk)+RAUT(jk)-raubase
         enddo
!       endif
!
       return
       end
