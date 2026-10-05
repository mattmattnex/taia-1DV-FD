!
!
!***********************************************************************
!***********************************************************************
       subroutine DYNKZ
       Use comrunmod
       Use comdynmod
       Implicit none
       integer :: jk,jkk
       real::zden,zden0,zden1,zdiscr
       real::zxld,zdxld,zdxld0,zdxld1,zdxlu,zdxlu0,zdxlu1
       real::zecart,zlope,zrhs0,zss,zxdis,zydis,zxlu,zxmld2
       real::zsum,zsum0,zsum1,zsum2
       real:: ZY(jpzt),ZZ(jpzt)
       real,external::SIGMA0
!
! Calcul de la densite potentielle, Brunt Vaissalla, cisaillement

       do jk=1,nzt-1
         RAUT(jk)=raubase+SIGMA0(0.,TN(jk),SN(jk))
       enddo
       do jk=2,nzt-1
         zlope=(RAUT(jk-1)-RAUT(jk))/(DEPT(jk-1)-DEPT(jk))
         RAUW(jk)=RAUT(jk-1)+zlope*(DEPW(jk)-DEPT(jk-1))
         BN(jk)=-g/raubase*(RAUT(jk-1)-RAUT(jk))/DE3W(jk)
       enddo
       do jk=2,nzt-1
         SH(jk)=((UN(jk-1)-UN(jk))/DE3W(jk))**2+((VN(jk-1)-VN(jk))/&
           DE3W(jk))**2
       enddo
!
! Calcul des longueurs de melange (eq. de 2ieme degre)
!
       do jk=2,nzt-1
         zrhs0=EN(jk)*raubase/g
         zden0=RAUW(jk)
         zden=RAUT(jk-1)
         zden1=zden0
         zxdis=0.
         zydis=DEPW(jk)-DEPT(jk-1)
         zsum0=(zden0-zden)*zydis/2.
         zsum=zsum0
         zecart=zsum-zrhs0
         if (zecart.ge.0.) goto 1001
         if (jk.eq.2) goto 1002
         zxdis=zydis
         do jkk=jk-2,1,-1
           zden1=zden
           zden=RAUT(jkk)
           zydis=DE3W(jkk+1)
           zsum0=zydis*(2*zden0-zden-zden1)/2.
           zsum=zsum+zsum0
           zecart=zsum-zrhs0
           if (zecart.ge.0.) goto 1001
           zxdis=zxdis+zydis
         enddo
 1002    zxdis=DEPW(jk)
         zecart=0
 1001    zxlu=zxdis
         if (ABS(zecart).le.ecartmx) goto 1011
         zsum=zecart-zsum0
         zsum1=(zden0-zden1)*zydis
         zsum2=-zydis*(zden-zden1)/2.
         zss=ABS(zsum2/zydis)
         if (zss.le.xlims) then
           zdxlu0=-zsum/zsum1
           zdxlu1=-1.
           zdxlu=zdxlu0
         else
           zdiscr=zsum1*zsum1-4.*zsum*zsum2
           zdiscr=SQRT(zdiscr)
           zdxlu0=(-zsum1+zdiscr)/(2.*zsum2)
           zdxlu1=(-zsum1-zdiscr)/(2.*zsum2)
           zdxlu=zdxlu0
           if (zdxlu.lt.-xlimd.or.zdxlu.gt.1.+xlimd) zdxlu=zdxlu1
         endif
         zxlu=zxlu+zdxlu*zydis
 1011    ZY(jk)=zxlu
         zden=RAUT(jk)
         zden1=zden0
         zxdis=0.
         zydis=DEPT(jk)-DEPW(jk)
         zsum0=zydis*(zden0-zden)/2.
         zsum=zsum0
         zecart=zsum+zrhs0
         if (zecart.lt.0.) goto 1021
         if (jk.eq.nzt-1) goto 1022
         zxdis=zydis
         do jkk=jk+1,nzt-1
           zden1=zden
           zden=RAUT(jkk)
           zydis=DE3W(jkk)
           zsum0=zydis*(2.*zden0-zden-zden1)/2.
           zsum=zsum+zsum0
           zecart=zsum+zrhs0
           if (zecart.le.0.) goto 1021
           zxdis=zxdis+zydis
         enddo
 1022     zxdis=DEPW(nzt)-DEPW(jk)
         zecart=0.
 1021    zxld=zxdis
         if (ABS(zecart).le.ecartmx) goto 1031
         zsum=zecart-zsum0
         zsum1=(zden0-zden1)*zydis
         zsum2=-zydis*(zden-zden1)/2.
         zss=ABS(zsum2/zydis)
         if (zss.le.xlims) then
           zdxld0=-zsum/zsum1
           zdxld1=-1.
           zdxld=zdxld0
         else
           zdiscr=zsum1*zsum1-4.*zsum*zsum2
           zdiscr=SQRT(zdiscr)
           zdxld0=(-zsum1+zdiscr)/(2.*zsum2)
           zdxld1=(-zsum1-zdiscr)/(2.*zsum2)
           zdxld=zdxld0
           if (zdxld.lt.-xlimd.or.zdxld.gt.1.+xlimd) zdxld=zdxld1
         endif
         zxld=zxld+zdxld*zydis
 1031    ZZ(jk)=zxld
       enddo
!
       do jk=2,nzt-1
         XMLDF(jk)=AMIN1(ZY(jk),ZZ(jk))
         XMLDF(jk)=AMAX1(1.E-8,XMLDF(jk))
         zxmld2=AMAX1(ZY(jk)*ZZ(jk),1.E-16)
         XMLDS(jk)=SQRT(zxmld2)
       enddo
!
! Calcul des coefficients de diffusion
!
       do jk=1,nzt-1
         AVM(jk)=AMAX1(avmb,xdiff*XMLDF(jk)*SQRT(EN(jk)))*TMASK(jk)
! MB : coefficient de turbulence pour la chaleur = AVM * Prandtl :
         AVT(jk)=xpdl*AVM(jk)
         AVE(jk)=fave*AVM(jk)
       enddo
       AVT(nzt)=0.
       AVM(nzt)=0.

!
       do jk=1,nzt
         DAVE(jk)=DBLE(AVE(jk))
         DAVM(jk)=DBLE(AVM(jk))
         DAVT(jk)=DBLE(AVT(jk))
       enddo
!
       return
       end
