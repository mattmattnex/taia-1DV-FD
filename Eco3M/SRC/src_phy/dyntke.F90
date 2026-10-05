!
!
!***********************************************************************
!***********************************************************************
       subroutine DYNTKE
!
! Integration de TKE
!
       Use comrunmod
       Use comdynmod
       Implicit none
       Integer :: jk
       real(8):: dzxxl1,dzzl1,dzzv1,dzzb1
       real(8):: zdtke,zsur,zttke
       real(8):: ZWW0(jpzt),ZXXD0(jpzt),ZXXL0(jpzt),ZYY0(jpzt)
       real(8):: ZZL0(jpzt),ZZV0(jpzt),ZZB0(jpzt)
!
       EN(1)=AMAX1(emin0,bb*SQRT(taux**2+tauy**2)/raubase)
       FLUXSURTKE(nind)=FLUXSURTKE(nind)+(EN(1)-EB(1))
       EA(1)=EN(1)
       do jk=2,nzt-1
         DZW(jk)=-0.5*dtsd*(DAVE(jk-1)+DAVE(jk))/DE3T(jk-1)/DE3W(jk)
         DZY(jk)=-0.5*dtsd*(DAVE(jk)+DAVE(jk+1))/DE3T(jk)/DE3W(jk)
         dzxxl1=1.5*dtsd*dxdiss/DBLE(XMLDS(jk))*DSQRT(DBLE(EN(jk)))
         DZX(jk)=dun-DZW(jk)-DZY(jk)+dzxxl1
         dzzl1=0.5*dtsd*dxdiss/DBLE(XMLDS(jk))*DSQRT(DBLE(EN(jk)))*&
              DBLE(EN(jk))
         dzzv1=dtsd*DAVM(jk)*DBLE(SH(jk))
         dzzb1=-dtsd*DAVT(jk)*DBLE(BN(jk))
         DZZ(jk)=DBLE(EN(jk))+dzzl1+dzzv1+dzzb1
! Sauvegarde
         ZWW0(jk)=SNGL(DZW(jk))
         ZYY0(jk)=SNGL(DZY(jk))
         ZXXD0(jk)=ZWW0(jk)+ZYY0(jk)
         ZXXL0(jk)=SNGL(dzxxl1)
         ZZL0(jk)=SNGL(dzzl1)
         ZZV0(jk)=SNGL(dzzv1)
         ZZB0(jk)=SNGL(dzzb1)
       enddo
! Traitement de surface
       DZW(2)=-0.5*dtsd*((2.*DE3T(2)+DE3T(1))*DAVE(2)-&
             DE3T(1)*DAVE(3))
       DZW(2)=DZW(2)/DE3T(2)/DE3W(2)/DE3T(1)
       DZX(2)=dun+dtsd*DAVE(2)*(DE3T(1)+DE3T(2))/DE3T(1)/DE3T(2)/DE3W(2)

       DZX(2)=DZX(2)+1.5*dtsd*dxdiss/DBLE(XMLDS(2))*DSQRT(DBLE(EN(2)))

       ZWW0(2)=SNGL(DZW(2))
       ZXXD0(2)=ZWW0(2)+ZYY0(2)
       DZZ(2)=DZZ(2)-DZW(2)*DBLE(EN(1))
       zsur=SNGL(DZW(2))*EN(1)
! Traitement du fond
       DZZ(nzt-1)=DZZ(nzt-1)-DZY(nzt-1)*DBLE(EN(nzt))
! Inversion
       DZR(2)=DZX(2)
       do jk=3,nzt-1
         DZR(jk)=DZX(jk)-DZW(jk)*DZY(jk-1)/DZR(jk-1)
       enddo
       DZX(2)=DZZ(2)
       do jk=3,nzt-1
         DZX(jk)=DZZ(jk)-DZW(jk)/DZR(jk-1)*DZX(jk-1)
       enddo
       DZW(nzt-1)=DZX(nzt-1)/DZR(nzt-1)
       do jk=nzt-2,2,-1
         DZW(jk)=(DZX(jk)-DZY(jk)*DZW(jk+1))/DZR(jk)
       enddo
! Champ after
       do jk=2,nzt
         EA(jk)=TMASK(jk)*AMAX1(emin,SNGL(DZW(jk)))
       enddo
! Swap des tableaux
       do jk=1,nzt
         EB(jk)=EN(jk)
         EN(jk)=EA(jk)
       enddo
! Sauvegarde tendances
!       if (mtke.eq.1) then
         do jk=1,nzt
           TKEMOY(jk)=TKEMOY(jk)+EN(jk)
         enddo
!       endif
!       if (mlength.eq.1) then
         do jk=1,nzt
           TBN(jk)=TBN(jk)+BN(jk)
           TSH(jk)=TSH(jk)+SH(jk)
           TMLDF(jk)=TMLDF(jk)+XMLDF(jk)
         enddo
!       endif
!       if (mdifft.eq.1) then
         do jk=1,nzt
           DIFFT(jk)=DIFFT(jk)+AVT(jk)
         enddo
!       endif
!       if (mtket.eq.1) then
!         do jk=1,nzt
         do jk=2,nzt-1
           zttke=-ZWW0(jk)*EN(jk-1)-ZYY0(jk)*EN(jk+1)+ZXXD0(jk)*EN(jk)
           TDFTKE(jk)=TDFTKE(jk)+zttke/dts
           zdtke=ZZL0(jk)-ZXXL0(jk)*EN(jk)
           TDSTKE(jk)=TDSTKE(jk)+zdtke/dts
           TBNTKE(jk)=TBNTKE(jk)+ZZB0(jk)/dts
           TSHTKE(jk)=TSHTKE(jk)+ZZV0(jk)/dts
         enddo
         TDFTKE(2)=TDFTKE(2)-zsur
!       endif
!
        return
        end
