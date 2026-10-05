!
!=======================================================================
!
       subroutine WBIOPRINT
       Use comrunmod
       Use comdynmod
       Use mod_varphy_coupl
!       include 'combio'
!
       write(99,1001) DEPTHZE(nind-1),XPAR0PLUS(nind-1),XNEBUL(nind-1)
 1001  format(/,2x,'ZE depth = ',f8.1,' m',3x,'PAR0+ = ',f8.2,' W/m2',&
        3x,'NEBUL = ',f8.3)
!       zflux=FLUXAIR(nind-1,5)*DEPT(1)
!       write(99,1002) VPCO2(nind-1),VITTRANS(nind-1),zflux
! 1002  format(2x,'pCO2 = ',0pf9.2,' ppm',3x,'Vit. Transfert = ',f9.2,
!     &  ' cm/h',3x,'Flux = ',1pe12.4,' mmolC/m2/sec')
!
       write(99,1003) CNMTRA(1),CNMTRA(2),CNMTRA(3),CNMTRA(11),&
          CNMTRA(12),CNMTRA(4),CNMTRA(5),CNMTRA(6),CNMTRA(10)
       do jk=1,nzt-1,2
         write(99,1004) jk,DEPT(jk),BIOAVER(jk,1),BIOAVER(jk,2),&
             BIOAVER(jk,3),BIOAVER(jk,11),BIOAVER(jk,12),BIOAVER(jk,4),&
             BIOAVER(jk,5),BIOAVER(jk,6),BIOAVER(jk,10)
       enddo
 1004  format(2x,i3,f6.1,9f9.4)
 1003  format(16x,9(a3,6x))
!
       return
       end
