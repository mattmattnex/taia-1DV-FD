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
   function f_upt_monod_cell(ili,jcol,ibio,inut,k_demi,Vmax)
!
!--------------------------------------------------------------------

 This function calculates the uptake rate of a nutrient
! using the Monod-Michaelis-Menten model
!
! Variables and parameters and their units:
!------------------------------------ -----
! ibio: index of one of the concentrations (C or N or P...)
!         of the organism consuming the nutrient
! inut: parameter giving the index of the concentration
!         of the limiting nutrient (specified because it may differ from jcol)
! k_demi: half-saturation constant, (mol X)/l
! Vmax  : SPECIFIC maximum uptake rate, in  molX.cell-1.s^{-1}
!
! f_upt_monod: uptake rate, in (mol X).l-1 s^{-1}
!
! CAUTION:
!    * Vmax must be in mol X.cell-1.s-1
!
!
! Ref: Monod, J. (1950) Annales de l'Institut Pasteur, 79:390-410
!      Michaelis, L. & Menten, M. L. (1913), 49:333-369
!
! Melika Baklouti 01/07/2011, 01/2026
!--------------------------------------------------------------------
    !
    !-- Global variables 
    ! Use VAR_GLOBAL

    Implicit None

    Real(8) :: f_upt_monod_cell(nx_min:nx_max,ny_min:ny_max,nz_max)

    !-- arguments
    Real(8) :: ibio,inut,k_demi,Vmax
    Integer :: ili,jcol

    !-- variables locales
    Integer :: iibio,iinut,iorg_plancton,iX,iCell
    Character(L_VAR)  :: cell
    Integer::i,j,k
    !logical :: isnan

    !-- nutrient concentration index
    iinut=nint(inut)

    !-- Planktonic biomass organism:
    iibio = nint(ibio)
    iorg_plancton = var(iibio)%idorg 


    !-- cell concentration (abundance) of plankton:
    cell = 'cell'
    iCell= f_idorg2id(iorg_plancton,cell)

    do k=1,nz_max
        do j=ny_min,ny_max
            do i=nx_min,nx_max
                if (var(iinut)%conc(i,j,k) > Vmax*var(iinut)%conc(i,j,k)/(var(iinut)%conc(i,j,k) + k_demi) &
                                                              * var(iCell)%conc(i,j,k) * dt_bio) then
                  f_upt_monod_cell(i,j,k) = &
                     (Vmax * var(iinut)%conc(i,j,k)/(var(iinut)%conc(i,j,k) + k_demi)) * var(iCell)%conc(i,j,k)
                else 
                  f_upt_monod_cell(i,j,k) = 0.d0
                endif
           enddo
      enddo
  enddo
      
end function f_upt_monod_cell
!--------------------------------------------------------------------
