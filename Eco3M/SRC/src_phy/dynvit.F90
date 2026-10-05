!
!
!***********************************************************************
!***********************************************************************
       subroutine DYNVIT
!
! Integration des vitesses
!
       Use comrunmod
       Use comdynmod
!
! Coriolis
       do jk=1,nzt
         UA(jk)=-f*VN(jk)
         VA(jk)=f*UN(jk)
       enddo
! initialisation U
       dfact=2.*dtsd
       do jk=1,nzt-1
         DZW(jk)=-dfact*DAVM(jk)/DE3T(jk)/DE3W(jk)
         DZY(jk)=-dfact*DAVM(jk+1)/DE3T(jk)/DE3W(jk+1)
         DZX(jk)=dun-DZW(jk)-DZY(jk)
         DZZ(jk)=DBLE(UB(jk))+dfact*DBLE(UA(jk))
       enddo
! Surface
       DZW(1)=0.
       DZX(1)=dun-DZY(1)
       DZZ(1)=DZZ(1)+dfact/DE3T(1)*DBLE(taux)/DBLE(raubase)
! Fond
       DZY(nzt-1)=0.
       DZX(nzt-1)=dun-DZW(nzt-1)+dfact*2.*DAVM(nzt-1)/ &
            DE3W(nzt)/DE3T(nzt)
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
         UA(jk)=SNGL(DZW(jk))
       enddo
!
! initialisation V
       do jk=1,nzt-1
         DZW(jk)=-dfact*DAVM(jk)/DE3T(jk)/DE3W(jk)
         DZY(jk)=-dfact*DAVM(jk+1)/DE3T(jk)/DE3W(jk+1)
         DZX(jk)=dun-DZW(jk)-DZY(jk)
         DZZ(jk)=DBLE(VB(jk))+dfact*DBLE(VA(jk))
       enddo
! Surface
       DZW(1)=0.
       DZX(1)=dun-DZY(1)
       DZZ(1)=DZZ(1)+dfact/DE3T(1)*DBLE(tauy)/DBLE(raubase)
! Fond
       DZY(nzt-1)=0.
       DZX(nzt-1)=dun-DZW(nzt-1)+dfact*2.*DAVM(nzt-1)/ &
                              DE3W(nzt)/DE3T(nzt)
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
         VA(jk)=SNGL(DZW(jk))
       enddo
! Filtre temporel et swap des tableaux
       do jk=1,nzt
         UN(jk)=UN(jk)+tgamma*(UB(jk)+UA(jk)-2.*UN(jk))
         VN(jk)=VN(jk)+tgamma*(VB(jk)+VA(jk)-2.*VN(jk))
       enddo
       do jk=1,nzt
         UB(jk)=UN(jk)
         UN(jk)=UA(jk)*TMASK(jk)
         VB(jk)=VN(jk)
         VN(jk)=VA(jk)*TMASK(jk)
       enddo
!
       return
       end
