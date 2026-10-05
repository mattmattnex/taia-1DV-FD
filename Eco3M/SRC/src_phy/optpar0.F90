!
! ======================================================================

       subroutine OPTPAR0
       Use comrunmod
       Use comdynmod
       Use mod_varphy_coupl

! Classe de vent
       jwind=1
       if (wind.gt.2.) jwind=2
       if (wind.gt.7.) jwind=3
       if (wind.gt.13.) jwind=4

! Angle solaire du soleil en deg.
       ahsun=pprtd*(timeminu/hourm-day/2.)*ppi/(day/2.)

! Hauteur du soleil
       zsin=SIN(xlat*ppdtr)*SIN(dcln*ppdtr)+ &
       COS(xlat*ppdtr)*COS(ahsun*ppdtr)*COS(dcln*ppdtr)
       zsin=ABS(zsin)
       zsin0=zsin
       zsin=AMAX1(zsin,1.E-4)
       hsun=pprtd*ASIN(zsin)

! Angle zenithal en deg,
       azsun=90.-hsun

! Discretisation de l'angle zenithal (1-19)
       jaz1=1+int(azsun/5.)
       if (azsun.ge.90.) jaz1=18
       jaz2=jaz1+1
       zaz1=5.*(jaz1-1.)
       zfct=1.

! Calcul de PAR0+ et PAR0- en ciel clair
       xpar0p=0.
       xpar0m=0.       
       do jl=1,nzt,nfrog
         zfcl=1.
         if ((jl-1)*(nzt-jl).eq.0) zfcl=0.5
         zrp1=SRDIR(jl,jaz1)+SRDIF(jl,jaz1)
         zrp2=SRDIR(jl,jaz2)+SRDIF(jl,jaz2)
         zrm1=SRDIR(jl,jaz1)*(1.-RDIR(jaz1,jwind))+ &
             SRDIF(jl,jaz1)*(1-rdif)
         zrm2=SRDIR(jl,jaz2)*(1.-RDIR(jaz2,jwind))+ &
             SRDIF(jl,jaz2)*(1-rdif)
         SR0P(jl)=(zrp1+(zrp2-zrp1)*(azsun-zaz1)/5.)*excen
         xpar0p=xpar0p+SR0P(jl)*zfct*zfcl*xdwl
         SR0M(jl)=(zrm1+(zrm2-zrm1)*(azsun-zaz1)/5.)*excen
         xpar0m=xpar0m+SR0M(jl)*zfct*zfcl*xdwl
       enddo

! Prise en compte de la nebulosite
!      ztrsw=1.-0.38*xneb-0.38*xneb*xneb
       zw0p043=xpar0p/0.43
       zqsr=qsr/(1.-albedo_phy)
       if (zw0p043.ne.0.) then
         ztrsw=zqsr/zw0p043
       else
         ztrsw=0.3
       endif
       if (ztrsw.gt.1.) ztrsw=1.
       zrdsw=1.-ztrsw
       zrdpar=0.75*zrdsw
       trpar=1.-zrdpar

       do jl=1,nzt,nfrog
         SR0P(jl)=SR0P(jl)*trpar
         SR0M(jl)=SR0M(jl)*trpar
       enddo

       xpar0m=xpar0m*trpar
       xpar0p=xpar0p*trpar

       return
       end
