!     
!     
!***********************************************************************
!***********************************************************************
      subroutine DYNFLUXINIT
      Use comrunmod
      Use comdynmod
      mflagflux=0
   
!     Initialisation flux surface    
      if (tdebut0.lt.xperflux0.or.tdebut0.gt.xperflux0+xperfluxt) then
         mflagflux=1
         return
      endif

      do jt=1,nfluxt
         zt=xperflux0+(jt-1)*xperflux/day
         if (tdebut0.lt.zt) goto 11
      enddo
 11   nflux0=jt-1
      idebut=int(tdebut/dt+0.0001)
      izt=int((zt*day-xperflux)*hourm/dt+0.0001)
      npasflux0=idebut-izt
      do jt=1,nflux0
!        barrier.n --- 994 changed to 94
         read(94,*) iy,im,id,ih,fse0,fle0,rsw0,rlw0,zx,zx,ustr0,vstr0
      enddo
      zday0=tdebut0-xperflux0
      zday1=zday0+xperflux/day

      nday0=INT(zday0)
      nday1=INT(zday1)
      dayflux=zday0+dt/(2.*hourm*day)
      dayflux0=dayflux+xperflux0
      if (dayflux0.lt.0.) dayflux0=dayflux0+xperfluxt
      if (dayflux0.gt.xperfluxt) dayflux0=dayflux0-xperfluxt
      if (nday1.gt.INT(xperfluxt)) then
         if (nday0.lt.INT(xperfluxt)) then
            mflagflux=3
            return
         else
            nday1=0
            if (dayflux.gt.xperfluxt) then
               dayflux=dayflux-xperfluxt
            endif
         endif
      endif
!     barrier.n --- modified barrier.n
!     994 changed to 94 in the next few lines
      read(94,*,end=201) iy,im,id,ih,fse1,fle1,rsw1,rlw1,zx,zx,&
                                                    ustr1,vstr1
      goto 202
 201  rewind(94)
      read(94,*,end=201) iy,im,id,ih,fse1,fle1,rsw1,rlw1,zx,zx,&
                                                    ustr1,vstr1
 202  continue

!     Initialisation rappel colonne T et S   
      rewind(95)
 101  continue
      read(95,*,end=109) dayts0,dayts1
      do jk=1,nzt
         read(95,*) XTS(jk),TRP(jk),SRP(jk)
      enddo
      if (dayflux0.ge.dayts0.and.dayflux0.lt.dayts1) goto 111
      goto 101
     
 109  mflagflux=2
      return
     
 111  continue
     
!     Ecriture initialisation flux.     
!      call WDYNFLUX

      return
      end
