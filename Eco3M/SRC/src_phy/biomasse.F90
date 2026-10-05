!
!
! **********************************************************************
! **********************************************************************
       subroutine BIOMASSE
!
! Appele par biologie.f a chaque passage: TOUTES LES SUBROUTINES APPELEES
!    PAR biomasse.f SE TROUVENT DANS biomodel.f
!
! Calculs des termes non conservatifs: boucle principale en z.
!
! Programme en parties:
!   1. Echanges des variables de formes vectorielles a forme scalaire
!         au niveau nk = BIOVARIABLE
!   2. Calculs des differents termes non conservatifs = BIOMODEL 
!             subroutines dans biomodel.f
!   3. Calculs des rappels dynamiques = BIORAP
!   4. Integration avant et verification de la positivite des champs
!   5. Inverse de 1: retour de l'information scalaire sous forme
!         vectorielle = BIOELBAIRAV
!
!            de traceurs updates = BIOAVANT
!
! ======================================================================
!
       Use comrunmod
       Use comdynmod
       USE mod_varphy_coupl

! Debut boucle en z
! 	barrier.n
!	No reason why the BIOAVANT function should be called for each z-level
!	In BIOAVANT, it is already done in 3D; i.e including the z-dimension
!       do jk=1,nzt
!          write(*,*) "+++++++++++++++++++++++++++++++jk"
!         nk=jk
!         call BIOAVANT
!       enddo

      return
      end
