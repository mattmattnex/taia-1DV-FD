!     
!     
!***********************************************************************
!***********************************************************************
      subroutine DYNSOLDAY
      Use comrunmod
      Use comdynmod
      integer:: jt, ji,ii ,ifreq

!     Calcul des moyennes journalieres flux de chaleur (assimilation)     
      do jt=1,jpant
         QSOLDAY(jt)=0.
      enddo

 301  continue
      ifluxt=0
      zfluxt=0.
      ifreq = nfreqflux
      ii=0
      do jt=1,njourflux
         iflux=0
         zflday=0.
         do ji=1,ifreq
!
            ii=ii+1
!           barrier.n --- modified 2015-09-15
!           994 has been replaced by 94
            read(94,*) iy,im,id,ih,zfse,zfle,zsw,zlw
            if (zsw.ne.0.) then
               iflux=iflux+1
               ifluxt=ifluxt+1
               zflday=zflday+zsw
               zfluxt=zfluxt+zsw
            endif
         enddo
         QSOLDAY(jt)=zflday/iflux
      enddo
      qsolmean=zfluxt/ifluxt
      rewind(94)
      return
      end
