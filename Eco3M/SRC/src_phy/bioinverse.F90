!
!
! **********************************************************************
! **********************************************************************
       subroutine BIOINVERSE
!
! Appele par BIODYN: inversion d'une matrive triangulaire (resolution
!    implicate de la diffusion verticale)
!
! ======================================================================
!
       Use comrunmod
       Use comdynmod
       USE mod_varphy_coupl
       Implicit none
       Integer:: jk
!
! Initialisation pour le traceur
!
       do jk=1,nzt-1
         DZW(jk)=-dts*DAVT(jk)/DE3T(jk)/DE3W(jk)
         DZY(jk)=-dts*DAVT(jk+1)/DE3T(jk)/DE3W(jk+1)
         DZX(jk)=dun-DZW(jk)-DZY(jk)
       enddo
! 
! Traitement de surface: flux deja pris en compte
!
       DZW(1)=0.
       DZX(1)=dun-DZY(1)
!
! Traitement du fond: pas de flux
!
       DZY(nzt-1)=0.
       DZX(nzt-1)=dun-DZW(nzt-1)
!
! Inversion matrice tridiagonale
!
       DZR(1)=DZX(1)
       do jk=2,nzt-1
         DZR(jk)=DZX(jk)-DZW(jk)*DZY(jk-1)/DZR(jk-1)
       enddo

       DZX(1)=DBLE(WORK1(1))
       do jk=2,nzt-1
         DZX(jk)=DBLE(WORK1(jk))-DZW(jk)/DZR(jk-1)*DZX(jk-1)
       enddo
       DZW(nzt-1)=DZX(nzt-1)/DZR(nzt-1)
       do jk=nzt-2,1,-1
         DZW(jk)=(DZX(jk)-DZY(jk)*DZW(jk+1))/DZR(jk)
       enddo
!
! Update
!
       do jk=1,nzt-1
!MB         WORK(jk)=SNGL(DZW(jk))
         WORK(jk)=DZW(jk)
       enddo
!
       return
       end
