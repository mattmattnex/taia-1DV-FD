!
!
!***********************************************************************
!***********************************************************************

       subroutine BIOINIT 
!***********************************************************************

       Use comrunmod
       Use comdynmod
       use mod_varphy_coupl
       
! Initialisation photosynthese (Rapport !/Chl)
       mflagbio=0
       npurday0=0
       npurd0=0
       fmlde0=0.
       fzede0=0.
       fmlpur0=0.
       fzepur0=0.

       fzede=DEPT(nzt/2)
       fmlde=DEPT(nzt/2)
       fmlpur=2.*xpurmx
       fzepur=2.*xpurmx
       
! Initialisation Flux interface ocean atmosphere
       do jtr=1,jptract
          BIOAIR(jtr)=0.
          do jt=1,jptemps
             FLUXAIR(jt,jtr)=0.
          enddo
       enddo

! initialisation moyenne sauvegarde
       do jtr=1,jptract
          do jk=1,nzt
             BIOAVER(jk,jtr)=0.
             DIFFBIO(jk,jtr)=0.
             SEDBIO(jk,jtr)=0.
          enddo
       enddo
       

       do jtr=1,jpfluxt
          do jk=1,nzt
             FLUXBIO(jk,jtr)=0.
          enddo
       enddo

! Initialisation champ 1D
       do jt=1,jptemps
          DEPTHZE(jt)=0.
          XNEBUL(jt)=0.
          XPAR0PLUS(jt)=0.
       enddo
       nind0=0

!     Initialisation indices
       do jtr=1,jptract
         do jk=1,jpzt
            MNEG(jk,jtr)=0
            MNEGB(jk,jtr)=0
         enddo
      enddo

! Initialisation parametres bio-optiques
! barrier.n: commented out for the dyfamed config
!      call OPTDATA

! Initialisation parametres astronomiques     
      call OPTDAILY
      call OPTASTRO

      return
      end
