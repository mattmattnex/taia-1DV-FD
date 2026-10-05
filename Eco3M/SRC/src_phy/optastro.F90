!
! ======================================================================
!
       subroutine OPTASTRO     
       Use comrunmod
       Use comdynmod
       Use mod_varphy_coupl

       isun=1

! Excentricite
       zjr=(ntimeday-2)*2.*ppi/xyear
       excen=1.000110+0.034221*COS(zjr)+0.00128*SIN(zjr)+ &
                   0.000719*COS(2.*zjr)+0.000077*SIN(2.*zjr)

! Declinaison
      zjr=ntimeday*2.*ppi/(xyear+1.)
      dcln=0.32281-22.984*COS(zjr)-0.3499*COS(2.*zjr)-0.1398*COS(3.*zjr) &
         + 3.7878*SIN(zjr)+0.03205*SIN(2.*zjr)+0.07187*SIN(3.*zjr)

! Angle horaire du couche du coleil et duree du jour
      zcos=-TAN(xlat*ppdtr)*TAN(dcln*ppdtr)
      if (ABS(zcos).ge.1.) then
        if (dcln*xlat.gt.0.) then
          daylength=24.
        else
          isun=0
        endif
      else
        ahsun=pprtd*ACOS(zcos)
        daylength=day*ahsun*ppdtr/ppi
      endif

      return
      end      
