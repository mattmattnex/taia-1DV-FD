!     
!     
!***********************************************************************
!***********************************************************************
      subroutine DYNDEF 
      Use comrunmod
      Use comdynmod
!     
      character*80 clpipo
!     
!     Lecture des parameters TKE et du run
!     
      open(90,file=cdir(1:ndir)//'/inputs/dynparam.i',status='old')
      read(90,'(a)') clpipo
      read(90,'(a)')  cflux
      read(90,'(a)') clpipo
      read(90,*) xperflux
      read(90,'(a)') clpipo
      read(90,*) xperfluxt
      read(90,'(a)') clpipo
      read(90,*) xperflux0
      read(90,'(a)') clpipo
      read(90,'(a)')  ctsinit
      read(90,'(a)') clpipo
      read(90,'(a)') crappel
      read(90,'(a)') clpipo
      read(90,*) msurft
      read(90,'(a)') clpipo
      read(90,'(a)') crapsurft
      read(90,'(a)') clpipo
      read(90,*) msurfs
      read(90,'(a)') clpipo
      read(90,'(a)') crapsurfs
      read(90,'(a)') clpipo
      read(90,*) xlat
      read(90,'(a)') clpipo
      read(90,*) g
      read(90,'(a)') clpipo
      read(90,*) xdiff
      read(90,'(a)') clpipo
      read(90,*) xdiss
      read(90,'(a)') clpipo
      read(90,*) avmb
      read(90,'(a)') clpipo
      read(90,*) avtb
      read(90,'(a)') clpipo
      read(90,*) emin0
      read(90,'(a)') clpipo
      read(90,*) emin
      read(90,'(a)') clpipo
      read(90,*) ecartmx
      read(90,'(a)') clpipo
      read(90,*) drag
      read(90,'(a)') clpipo
      read(90,*) rabs
      read(90,'(a)') clpipo
      read(90,*) xsi1
      read(90,'(a)') clpipo
      read(90,*) xsi2
      read(90,'(a)') clpipo
      read(90,*) raubase
      read(90,'(a)') clpipo
      read(90,*) rauair
      read(90,'(a)') clpipo
      read(90,*) xcp
      read(90,'(a)') clpipo
      read(90,*) tcs
      read(90,'(a)') clpipo
      read(90,*) xpdl
      read(90,'(a)') clpipo
      read(90,*) fave
      read(90,'(a)') clpipo
      read(90,*) bb
      read(90,'(a)') clpipo
      read(90,*) tgamma
      read(90,'(a)') clpipo
      read(90,*) albedo_phy
      read(90,'(a)') clpipo
      read(90,*) emiss
      read(90,'(a)') clpipo
      read(90,*) stefan
      read(90,'(a)') clpipo
      read(90,*) rapmax
      read(90,'(a)') clpipo
      read(90,*) rapprop
      read(90,'(a)') clpipo
      read(90,*) alfcte
      read(90,'(a)') clpipo
      read(90,*) alfmin
      read(90,'(a)') clpipo
      read(90,*) alfmax
      close(90)
!     
      f=2.*omega*SIN(xlat*(ppi/180.))
!     
      npasflux=int(xperflux*hourm/dt+0.001)
      njourflux=int(xperfluxt+0.001)
      nfreqflux=int(day/xperflux+0.001)
      ztot=xperfluxt*day/xperflux
      nfluxt=int(ztot+0.001)
!     
      if (njourflux.gt.jpant) then
         mflagtke=2
         return
      endif
!     
      dtsd=DBLE(dts)
      dfact=2.*dtsd
      dxdiss=DBLE(xdiss)
      dtcs=DBLE(tcs)

      return
      end
