!
!
! **********************************************************************
! **********************************************************************
       subroutine BIOSEDIM
!
! Appele par biologie.f a chaque passage: calcul des flux de sedimenta-
! tion
!
! ======================================================================
!
        Use comrunmod
        Use comdynmod
        Use mod_varphy_coupl

!        write(*,*) "Sed speed ", vitsed
!
       do jtr=1,jptract  ! barrier.n nbrprono -> jptract
         SEDFLUX(1,jtr)=-VITSED(jtr)*TENEUR(1,jtr)/DE3T(1)
         EXPORT(1,jtr)= 0.0d0
!         do jk=2,nzt
         do jk=2,nzt-2
          SEDFLUX(jk,jtr)=VITSED(jtr)*(TENEUR(jk-1,jtr)-TENEUR(jk,jtr))&
                                                               /DE3T(jk)
          EXPORT(jk,jtr)=-VITSED(jtr)*TENEUR(jk-1,jtr)  & 
      - DAVT(jk) *(TENEUR(jk-1,jtr) - TENEUR(jk,jtr))/DE3T(jk-1)
         enddo
           SEDFLUX(nzt-1,jtr)=VITSED(jtr)*(TENEUR(nzt-2,jtr)) &
                                                    /DE3T(nzt-1)
          EXPORT(nzt-1,jtr)=-VITSED(jtr)*(TENEUR(nzt-2,jtr)) &
      -DAVT(nzt-1) *(TENEUR(nzt-2,jtr) - TENEUR(nzt-1,jtr))/DE3T(nzt-2)
           SEDFLUX(nzt,jtr)=0.0d0
           EXPORT(nzt,jtr)=0.0d0
       enddo
       
!
       do jtr=1,jptract  ! barrier.n nbrprono -> jptract
         do jk=1,nzt
           SEDBIO(jk,jtr)=SEDBIO(jk,jtr)+SEDFLUX(jk,jtr)
         enddo
       enddo
!
       do jtr=1,jptract  ! barrier.n nbrprono -> jptract
         do jk=1,nzt-1
           TENEUR(jk,jtr)=TENEUR(jk,jtr)+dts*SEDFLUX(jk,jtr)
         enddo
       enddo
!
       return
       end
