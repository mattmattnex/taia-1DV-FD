!***************************************************************************
!***************************************************************************
!Copyright or © or Copr. CNRS/IRD/Université de la Méditerranée
!contributor(s) : Melika BAKLOUTI & Vincent FAURE (10/10/2006)
!
!m.baklouti@univmed.fr; vincent.faure@univmed.fr
!
!This software (Eco3M) is a computer program whose purpose is to perform 
!biogeochemical or coupled physical-biogeochemical modelling.
!
!This software is governed by the CeCILL license under French law and
!abiding by the rules of distribution of free software. You can  use, 
!modify and/ or redistribute the software under the terms of the CeCILL
!license as circulated by CEA, CNRS and INRIA at the following URL
!"http://www.cecill.info". 
!
!As a counterpart to the access to the source code and  rights to copy,
!modify and redistribute granted by the license, users are provided only
!with a limited warranty  and the software''s author,  the holder of the
!economic rights,  and the successive licensors  have only  limited
!liability. 
!
!In this respect, the user''s attention is drawn to the risks associated
!with loading,  using,  modifying and/or developing or reproducing the
!software by the user in light of its specific status of free software,
!that may mean  that it is complicated to manipulate,  and  that  also
!therefore means  that it is reserved for developers  and  experienced
!professionals having in-depth computer knowledge. Users are therefore
!encouraged to load and test the software''s suitability as regards their
!requirements in conditions enabling the security of their systems and/or 
!data to be ensured and,  more generally, to use and operate it in the 
!same conditions as regards security. 
!
!The fact that you are presently reading this means that you have had
!knowledge of the CeCILL license and that you accept its terms.
!***************************************************************************
!***************************************************************************
!--------------------------------------------------------------------
!
    function f_mort_quadratique_cell_2026(ili,jcol,iX,iCell,kmort1)
!
!--------------------------------------------------------------------
! This function calculates the mortality rate according to a QUADRATIC law due 
! to predation not explicitly represented in the model (closure term).
!
!
! Variables and parameters and their units:
!--------------------------------------- --
!
! kmort: mortality rate constant ((mol.X .l-1 s)^-1)
!
!
! Melika Baklouti 09/10/2007, 09/06/2017, 07/2025
!--------------------------------------------------------------------

    Implicit None

    !-- fonction:
    Real(8)             :: f_mort_quadratique_cell_2026(nx_min:nx_max,ny_min:ny_max,nz_max)

    !-- arguments
    Real(8)             :: kmort1, iCell, iX !-- cste spécifique de mortalité, (mol.X .m-3 s)^-1
    Integer             :: ili,jcol,iiCell,iiX,i,j,k
    
    !--variables locales:
    !logical :: isnan

    iiX = nint(iX)
    iiCell = nint(iCell)


    !-- Vitesse de mortalite  (mol.X .m-3 . s-1):
    do k=1,nz_max
        do j=ny_min,ny_max
            do i=nx_min,nx_max
                if (var(iiCell)%conc(i,j,k) > 1.1 *  kmort1*(var(iiCell)%conc(i,j,k))**2 * dt_bio) then 
                    f_mort_quadratique_cell_2026(i,j,k) = &
                        kmort1*(var(iiCell)%conc(i,j,k))**2 * (var(iiX)%conc(i,j,k)/var(iiCell)%conc(i,j,k)+1d-80)
                else
                    f_mort_quadratique_cell_2026(i,j,k) = 0.0d0
                endif
            enddo
        enddo
    enddo
End function f_mort_quadratique_cell_2026
!--------------------------------------------------------------------
