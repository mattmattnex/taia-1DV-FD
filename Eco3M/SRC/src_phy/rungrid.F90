!     
!     
!***********************************************************************
!***********************************************************************
      subroutine RUNGRID
!     
      Use comrunmod
      Use comdynmod
      USE mod_varphy_coupl
!     
      open(93,file=cdir(1:ndir)//'/'//cgrille,status='old')
      do jk=1,nzt
         read(93,*) DEPW(jk)
         write(*,*)'jk,depw(jk)',jk,depw(jk)
      enddo
      do jk=1,nzt-1
         TMASK(jk)=1.
         DEPT(jk)=0.5*(DEPW(jk)+DEPW(jk+1))
! MB         E3T(jk)=DEPW(jk+1)-DEPW(jk)
         DE3T(jk)=DEPW(jk+1)-DEPW(jk)
      enddo
      TMASK(nzt)=0.
      DEPT(nzt)=DEPW(nzt)
      do jk=2,nzt-1
! MB         E3W(jk)=DEPT(jk)-DEPT(jk-1)
         DE3W(jk)=DEPT(jk)-DEPT(jk-1)
      enddo
      DE3T(nzt)=DE3T(nzt-1)
      DE3W(nzt)=DE3T(nzt)
      DE3W(1)=DE3T(1)
      close(93)
!     
      do jk=1,nzt
         E3T(jk)=DE3T(jk)
         E3W(jk)=DE3W(jk)
!        DE3T(jk)=DBLE(E3T(jk))
!        DE3W(jk)=DBLE(E3W(jk))
      enddo
!     
      return
      end
