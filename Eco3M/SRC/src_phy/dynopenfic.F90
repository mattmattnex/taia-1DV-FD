!
!
!***********************************************************************
!***********************************************************************
       subroutine DYNOPENFIC
       Use comrunmod
       Use comdynmod
       Implicit None
       CHARACTER(LEN=100)::fic
       Integer:: nwldynfl
!
       nwdyn=mt+ms+mtke+mdifft
       if (nwdyn.gt.0) then
         fic = trim(croot(1:nroot))//'dyntrac'
         open(nundyn,file=fic,form='unformatted')
         write(nundyn) nwdyn,mt,ms,mtke,mdifft,nzt
       endif
       nwdynfl=mlength+mtket+mtt+mst
       nwldynfl=nlengthl*mlength+ntketl*mtket+nttl*mtt+nstl*mst
       if (nwdynfl.gt.0) then
         fic = trim(croot(1:nroot))//'dynflux'
        open(nundynfl,file=fic,form='unformatted')
         write(nundynfl) nwdynfl,mlength,mtket,mtt,mst,nzt
         write(nundynfl) nwldynfl,nlengthl,ntketl,nttl,nstl
       endif
       fic = trim(croot(1:nroot))//'dyn.1d'
       open(nundyn1d,file=fic)
!
       return
       end
