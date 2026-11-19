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
   function f_upt_monod_cell_multX(ili,jcol,ibio,inut,K_sat,Vmax,inut1,K1,inut2,K2,inut3,K3)
!
!--------------------------------------------------------------------
!
! This function calculates the absorption/assimilation (uptake) of a nutrient    
! associated with the inut index using the Monod-Michaelis-Menten model,
! taking into account other sources of the same biogenic element.   
! inut1 is the nutrient consumed first, followed by inut2, etc.
!
! Variables and parameters and their units:   
!----------------------------------------- 
! ibio  : consumer concentration index (in mol C or N or P, etc. per liter)
! inut  : parameter giving the index of the limiting nutrient 
!         concentration (specified because it may differ from jcol) 
! K     : half saturation constant, (mol X)/l
! Vmax  : maximum specific uptake rate, in  molX.cell-1.s^{-1}
!
! f_upt_monod : uptake rate, in (mol X).l-1 s^{-1}
!
! ATTENTION: 
!    * Vmax must be in mol X.cell-1.s-1
!
! Ref: Monod, J. (1950) Annales de l'Institut Pasteur, 79:390-410
!      Michaelis, L. & Menten, M. L. (1913), 49:333-369
!
! Melika Baklouti 08/2025
!--------------------------------------------------------------------
    !
    !-- Global variables 
    ! Use VAR_GLOBAL

    implicit None

    real(8) :: f_upt_monod_cell_multX(nx_min:nx_max,ny_min:ny_max,nz_max)

    !-- arguments
    real(8)          :: ibio,inut,K_sat,Vmax
    real(8),optional :: inut1,K1,inut2,K2,inut3,K3
    integer          :: ili,jcol

    !-- local variables
    integer           :: iibio,iinut,iorg_plancton,iX,iCell
    character(L_VAR)  :: cell
    integer           :: i,j,k
    integer           :: iinut1,iinut2,iinut3
    real(8)           :: flim,fluX1,fluX1_2

    !--  conversion of nutrient concentration indices into integers
    iinut=nint(inut)
    if (present (inut1)) iinut1=nint(inut1)
    if (present (inut2)) iinut2=nint(inut2)
    if (present (inut3)) iinut3=nint(inut3)

    !-- index of the corresponding plankton organism :
    iibio = nint(ibio)
    iorg_plancton = var(iibio)%idorg 


    !-- index of cell concentration (abundance):
    cell = 'cell'
    iCell= f_idorg2id(iorg_plancton,cell)

    ! initialization
    f_upt_monod_cell_multX = 0.d0

  if ( ((.not. present (inut1)) .and. (.not. present (inut2)) .and. (.not. present (inut3))) .or. (iinut==iinut1)) then 
     do k=1,nz_max
        do j=ny_min,ny_max
            do i=nx_min,nx_max
                if (var(iinut)%conc(i,j,k) > 1.1 * Vmax*var(iinut)%conc(i,j,k)/(var(iinut)%conc(i,j,k) + K_sat) &
                                                              * var(iCell)%conc(i,j,k) * dt_bio) then
                  f_upt_monod_cell_multX(i,j,k) = &
                     (Vmax * var(iinut)%conc(i,j,k)/(var(iinut)%conc(i,j,k) + K_sat)) * var(iCell)%conc(i,j,k)
		endif
           enddo
      enddo
    enddo
   elseif (iinut==iinut2) then
     do k=1,nz_max
        do j=ny_min,ny_max
            do i=nx_min,nx_max
               fluX1 =  1.0d0 - var(iinut1)%conc(i,j,k)/(var(iinut1)%conc(i,j,k) + K1)   
	       fluX1 = max(0.0d0,fluX1)
               flim=min(var(iinut2)%conc(i,j,k)/(var(iinut2)%conc(i,j,k) + K2),fluX1 )
               if (var(iinut)%conc(i,j,k) > 1.1 * Vmax*flim*var(iCell)%conc(i,j,k) * dt_bio) then
                  f_upt_monod_cell_multX(i,j,k) = Vmax * flim * var(iCell)%conc(i,j,k)
	       endif
           enddo
      enddo
    enddo
 elseif ((present (inut3)) .and. iinut==iinut3) then
     do k=1,nz_max
       do j=ny_min,ny_max
         do i=nx_min,nx_max
           fluX1_2 =  1.0d0 -var(iinut1)%conc(i,j,k)/(var(iinut1)%conc(i,j,k) + K1) & 
	                                 -var(iinut2)%conc(i,j,k)/(var(iinut2)%conc(i,j,k)+K2)
	   fluX1_2 = max(0.0d0,fluX1_2)
           flim=min(var(iinut3)%conc(i,j,k)/(var(iinut3)%conc(i,j,k) + K3),fluX1_2)
           if (var(iinut)%conc(i,j,k) > 1.1 * Vmax*flim*var(iCell)%conc(i,j,k) * dt_bio) then
              f_upt_monod_cell_multX(i,j,k) = Vmax * flim * var(iCell)%conc(i,j,k)
	   endif
         enddo
      enddo
    enddo
 endif 



      
end function f_upt_monod_cell_multX
!--------------------------------------------------------------------
