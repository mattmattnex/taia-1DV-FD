!
!=======================================================================
!
       subroutine WRUNFIN 
       Use comrunmod
       Use comdynmod
       Use mod_varphy_coupl
!       include 'combio'
!
       write(99,1001) nstep,timep,timep/(hourm*day)
       write(99,1011) float(nstepneg)/(nzt-1), &
                  float(nstepneg)/(nstep*(nzt-1))
       zn=nstepneg
       if (zn.eq.0) zn=1
       write(99,1012) float(nstepneg0)/(nzt-1),float(nstepneg0)/zn
 1001  format(/,2x,50('$'),//,2x,'nstep final = ',i10,2x,'timep = ',&
           1pe12.4,' mn = ',0pf8.2,' jours')
 1011  format(2x,'Nbr de pas de temps avec < 0: ',f13.3,3x,'Rapport = ',&
              f8.3)
 1012  format(2x,'Nbr de passage dans boucle < 0 en moyenne ',f13.3,&
       3x,'Rapport vs nbr pas de temps < 0 ',f8.3)
       write(99,1002) (CNMTRA(jtr),jtr=1,nbrprono)
 1002  format(/,2x,'Proportion de valeurs negatives %',/,9x,&
          10(4x,a3,2x))
!
       iwr=0
       do jk=1,nzt-1
         do jtr=1,nbrprono
           if (MNEGB(jk,jtr).ne.0) iwr=1
           goto 101
         enddo
       enddo
 101   continue
       do jk=1,nzt-1
         write(99,1003) DEPT(jk),(MNEGB(jk,jtr)*100./nstep,&
                jtr=1,nbrprono)
       enddo
 1003  format(0pf8.1,10f9.3)
!
       if (iwr.eq.0) goto 201
       write(99,1004) (CNMTRA(jtr),jtr=1,nbrprono)
 1004  format(/,2x,'Ecart de valeurs negatives totales vs base %',/,9x,&
         10(4x,a3,2x))
       do jk=1,nzt-1
         isum=0
         do jtr=1,nbrprono
           isum=isum+MNEG(jk,jtr)
         enddo
         if (isum.ne.0) then
           write(99,1003) DEPT(jk),((MNEG(jk,jtr)-MNEGB(jk,jtr))*100./&
           AMAX1(1.,FLOAT(MNEGB(jk,jtr))),jtr=1,nbrprono)
         endif
       enddo
 201   continue
!
       return
       end
