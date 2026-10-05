!     
!
!***********************************************************************
!***********************************************************************
       subroutine BIO1DUP
!
! ------------------------------------------------
! Moyenne et stockage des champs 1D de la biologie
! ------------------------------------------------
! nind = indice temporel
!
!       include 'comrun'
!       include 'comdyn'
        Use comrunmod
        Use comdynmod
        Use mod_varphy_coupl
!       include 'combio'
!

       if (nind0.ne.0) then
         DEPTHZE(nind)=DEPTHZE(nind)/nind0
         XNEBUL(nind)=XNEBUL(nind)/nind0
         XPAR0PLUS(nind)=XPAR0PLUS(nind)/nind0
       endif
       do jtr=1,jptract  ! barrier.n nbrprono -> jptract
         FLUXAIR(nind,jtr)=FLUXAIR(nind,jtr)/nsave
       enddo
!       VPCO2(nind)=VPCO2(nind)/nsave
!       VITTRANS(nind)=100.*hour*VITTRANS(nind)/nsave
!
       return
       end
