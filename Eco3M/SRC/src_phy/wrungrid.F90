!     
!=======================================================================
!     
      subroutine WRUNGRID
!     
      Use comrunmod
      Use comdynmod
!     
      do jk=1,nzt-1,2
         write(99,1005) jk,DEPW(jk),DEPT(jk),E3W(jk),E3T(jk)
      enddo
!-- Ajout MB
!      do jk=1,nzt-1
!         write(199,1006) jk,DEPW(jk),DEPT(jk),E3W(jk),E3T(jk),DE3T(jk)
!      enddo
 1005 format(i5,4f8.1)  
!-- Ajout MB
! 1006 format(i5,2X,4(f10.5,1X),E16.10)
!     
      return
      end
