!
!=======================================================================
!
       subroutine WDYNPRINT
       Use comrunmod
       Use comdynmod
!       include 'combio'
!
       ztt=timep/(hourm*day)
       write(99,1001) nstep,ztt,DEPTHML(nind-1)
       write(99,1006) FLUXVENTX(nind-1),FLUXVENTY(nind-1), &
              FLUXVENT(nind-1)                               
      write(99,1003) FLUXSOL(nind-1),FLUXNSOL(nind-1),      &
        FLUXADVT(nind-1)                                     
      write(99,1004) FLUXSAL(nind-1),FLUXADVS(nind-1)        
      ztott=FLUXSOL(nind-1)+FLUXNSOL(nind-1)+               &
        FLUXADVT(nind-1)                                     
      ztots=FLUXSAL(nind-1)+FLUXADVS(nind-1)                 
      write(99,1005) ztott,ztots                             
      do jk=1,nzt-1,2                                        
        write(99,1002) jk,DEPT(jk),TMOY(jk),SMOY(jk),DENSITE(jk),&
            TKEMOY(jk),DIFFT(jk)
       enddo
 1001  format(/,2x,'====> nstep = ',i10,3x,'Jour : ',f10.3,5x, &
         'ML depth = ',f8.1,' m')
 1002  format(2x,i3,f6.1,3f8.3,2f8.3)
 1003  format(2x,'Forcages T: SOL = ',1pe12.4,2x,'NSOL = ',e12.4,2x,&
         'ADV = ',e12.4)
 1004  format(2x,'Forcages S:',21x,'SURF = ',1pe12.4,2x,'ADV = ',e12.4)
 1005  format(2x,'Bilan forcages: T = ',1pe12.4,5x,'S = ',e12.4)
 1006  format(2x,'Vent X = ',f8.3,'  Y = ',f8.3,'  Norm = ',f8.3,&
                  ' m/sec')
!
       return
       end
