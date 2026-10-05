!
!
!***********************************************************************
!***********************************************************************
       subroutine DYNFLUXSUR
       Use comrunmod
       Use comdynmod
       USE mod_varphy_coupl
       Implicit none
       Integer::iy,im,id,ih
       real::zx,zq0,zq1,zqsr0,zqsr1
!
    
       npasflux0=npasflux0+1
       dayflux=dayflux+dts/(hour*day)
       dayflux0=dayflux0+dts/(hour*day)
       if (dayflux0.gt.xperfluxt) dayflux0=dayflux0-xperfluxt
       
       if (npasflux0.ge.npasflux) then
         fse0=fse1
         fle0=fle1
         rsw0=rsw1
         rlw0=rlw1
         ustr0=ustr1
         vstr0=vstr1
         nday0=nday1
         nday1=INT(dayflux+xperflux/day)
         npasflux0=0
          read(94,*,end=101) iy,im,id,ih,fse1,fle1,rsw1,rlw1,zx,zx,&
                                                        ustr1,vstr1
         goto 102
 101     rewind(94)
         write(99,1001) nstep,time,timeday,timeyr
          read(94,*) iy,im,id,ih,fse1,fle1,rsw1,rlw1,zx,zx, &
                                                ustr1,vstr1
       endif

       if (dayflux.gt.xperfluxt) then
         dayflux=dayflux-xperfluxt
         nday1=0
       endif
 1001  format(/,2x,'&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&',/,  &
        2x,'REWIND fichier forcages',/,2x,'nstep = ',i8,3x,'time = ', &
        f10.1,' mn',3x,'timeday = ',f7.1,' jours',3x,'timeyr = ',f5.1,&
               /,2x,'&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&') 
 102   continue
       zq0=fse0+fle0+rlw0
       zq1=fse1+fle1+rlw1
       q=(npasflux0*zq1+(npasflux-npasflux0)*zq0)/npasflux
       if (rsw0.eq.0.) then
         zqsr0=rsw0
       else
         zqsr0=alfcte*qsolmean
         if (QSOLDAY(nday0+1).gt.0.) then
           zqsr0=zqsr0+alfmax*(rsw0-qsolmean)
         else
           zqsr0=zqsr0+alfmin*(rsw0-qsolmean)
         endif
       endif  
       if (rsw1.eq.0.) then
         zqsr1=rsw1
       else
         zqsr1=alfcte*qsolmean
         if (QSOLDAY(nday1+1).gt.0.) then
           zqsr1=zqsr1+alfmax*(rsw1-qsolmean)
         else
           zqsr1=zqsr1+alfmin*(rsw1-qsolmean)
         endif
       endif  
       qsr=(npasflux0*zqsr1+(npasflux-npasflux0)*zqsr0)/npasflux
!-- 27/05/2016
       if (qsr < 0) qsr = 0.0
!------------------------------------
       qsr=qsr*(1.-albedo_phy)
       taux=(npasflux0*ustr1+(npasflux-npasflux0)*ustr0)/npasflux
       tauy=(npasflux0*vstr1+(npasflux-npasflux0)*vstr0)/npasflux
       ep=0.
!
       wind=SQRT(SQRT(taux*taux+tauy*tauy)/(drag*rauair))
       windx=taux/(rauair*drag*wind)
       windy=tauy/(rauair*drag*wind)
!
       return
       end
