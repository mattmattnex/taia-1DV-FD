!
!=======================================================================
!
       subroutine WDYNFLUX 
       Use comrunmod
       Use comdynmod
!       include 'combio'
!
       write(99,1006) cflux,xperflux,xperfluxt,xperflux0
 1006  format(/,2x,'Fichier forcages atmospheriques:',/,2x,a80,/,2x,   &
       'Periode forcages = ',f5.2,' heures',/,2x,'Longueur fichier = ' &
       ,f8.2,' jours',5x,'Enregistrement 1: ',f8.3,' jours')
       if (mflagwrite.eq.0) return
!
       write(99,1001) nfluxt,nflux0,dayflux
 1001  format(/,2x,'Initialisation Flux: nbr enregistrements = ',i5,/, &
       2x,'Enregistrement debut = ',i3,5x,'Jour initial flux = ',f6.2)
       write(99,1002) fse0,fle0,rsw0,rlw0,ustr0,vstr0
       write(99,1003) fse1,fle1,rsw1,rlw1,ustr1,vstr1
 1002  format(2x,'Before: sensible, latent,solaire, IR, U*, V*',/,5x,&
        1p,6e12.4)
 1003  format(2x,' After: sensible, latent,solaire, IR, U*, V*',/,5x,&
        1p,6e12.4)
       write(99,1004) timeday,dayts0,dayts1
 1004  format(/,2x,'Initialisation rappel profils T et S',/,2x,&
         'Jour : ',f6.2,' entre ',f6.2,' et ',f6.2)
       do jk=1,nzt-1,2
         write(99,1005) jk,DEPT(jk),TRP(jk),SRP(jk),XTS(jk)
       enddo
 1005  format(0p,i5,f8.1,3f8.3,1p,e12.4)
!
       return
       end
