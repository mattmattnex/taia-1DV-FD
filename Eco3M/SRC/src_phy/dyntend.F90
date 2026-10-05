!
!
!***********************************************************************
!***********************************************************************
       subroutine DYNTEND 
       Use comrunmod
       Use comdynmod
!       include 'combio'
!
! ----------------------------------------------------------------
! Moyenne, sauvegarde et initialisation sur nundynfl des tendances
! ----------------------------------------------------------------
! Traceurs dynamiques
!
!
! Longueur melange, Brunt Vaissaila, Cisaillement
!
       do jk=1,nzt
         TMLDF(jk)=TMLDF(jk)/nsave
         TBN(jk)=TBN(jk)/nsave
         TSH(jk)=TSH(jk)/nsave
       enddo
#if defined SAVE_PHY       
       if (mlength.eq.1) then
          write(nundynfl) (TMLDF(jk),jk=nzt,1,-1)
          write(nundynfl) (TBN(jk),jk=nzt,1,-1)
          write(nundynfl) (TSH(jk),jk=nzt,1,-1)
       endif
#endif
       do jk=1,nzt
         TMLDF(jk)=0.                   
         TBN(jk)=0.                   
         TSH(jk)=0.                   
       enddo
!
! Tendances equation TKE
!
       do jk=1,nzt
         TSHTKE(jk)=TSHTKE(jk)/nsave
         TBNTKE(jk)=TBNTKE(jk)/nsave
         TDFTKE(jk)=TDFTKE(jk)/nsave
         TDSTKE(jk)=TDSTKE(jk)/nsave
       enddo
#if defined SAVE_PHY       
       if (mtket.eq.1) then
         write(nundynfl) (TSHTKE(jk),jk=nzt,1,-1)
         write(nundynfl) (TBNTKE(jk),jk=nzt,1,-1)
         write(nundynfl) (TDFTKE(jk),jk=nzt,1,-1)
         write(nundynfl) (TDSTKE(jk),jk=nzt,1,-1)
       endif
#endif
       do jk=1,nzt
         TSHTKE(jk)=0.
         TBNTKE(jk)=0.
         TDFTKE(jk)=0.
         TDSTKE(jk)=0.
       enddo
!
! Tendances equation Temperature
!
       do jk=1,nzt
         TRAPT(jk)=TRAPT(jk)/nsave
         TDIFT(jk)=TDIFT(jk)/nsave
         TPENT(jk)=TPENT(jk)/nsave
       enddo
#if defined SAVE_PHY       
       if (mtt.eq.1) then
         write(nundynfl) (TRAPT(jk),jk=nzt,1,-1)
         write(nundynfl) (TDIFT(jk),jk=nzt,1,-1)
         write(nundynfl) (TPENT(jk),jk=nzt,1,-1)
       endif
#endif       
       do jk=1,nzt
         TRAPT(jk)=0.
         TDIFT(jk)=0.
         TPENT(jk)=0.
       enddo
!
! Tendances equation Salinite
!
       do jk=1,nzt
         TRAPS(jk)=TRAPS(jk)/nsave
         TDIFS(jk)=TDIFS(jk)/nsave
       enddo
#if defined SAVE_PHY       
       if (mst.eq.1) then
         write(nundynfl) (TRAPS(jk),jk=nzt,1,-1)
         write(nundynfl) (TDIFS(jk),jk=nzt,1,-1)
       endif
#endif
       do jk=1,nzt-1
         TRAPS(jk)=0.
         TDIFS(jk)=0.
       enddo
!
       return
       end
