!***************************************************************************
!***************************************************************************
!Copyright or © or Copr. CNRS/IRD/Université de la Méditerranée
!contributor(s) : Melika BAKLOUTI & Vincent FAURE (10/10/2006)
!
!melika.baklouti@univ-amu.fr; 
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
!> This file contains the subroutines dedicated to the calculation of source 
!! minus sink term (SMS) for each variable of the biogeochemical model, which are stored in the TEND array.
!! Used only in CALC mode (key CALC)
!! It also contains the subroutine dedicated to the update of biological concentrations (only in NCOUPL mod)
!! \date 24/08/2017 

#if defined  CALC 
!-----------------------------------------------------------------------------------------------
    subroutine eco3m_sms_init

!! \author  M Baklouti, V Faure
!-----------------------------------------------------------------------------------------------
  use mod_eco3m
 !$ use omp_lib   ! when openmp activated, loading of the library !MB: A VOIR
  implicit none
!
  integer:: istat

  Allocate (TEND(nx_min:nx_max, ny_min:ny_max, 1:nz_max,1:nbvar),STAT=istat)
  if (istat /= 0) write(*,*) 'pb with the allocation of TEND'


    end subroutine eco3m_sms_init
!-----------------------------------------------------------------------------------------------
    subroutine eco3m_compute_sms

!> Computes the source minus sink term for each state variable. If first computes the flux values
!! by calling the "call.inc" file, then it updates the TEND array.
!! \author  M Baklouti, V Faure
!-----------------------------------------------------------------------------------------------
  use mod_eco3m_vartypes
  use mod_eco3m
  use mod_eco3m_fprocess
  implicit none
!        
  integer :: ivar, jvar, kvar ! Loop indexes
  integer :: i,j,k

        ! Calculation of the flux values
        include "call.inc"

        TEND = 0.d0

        ! Sum of the sources/sinks for each variables.
        do ivar = 1, nbvar

            ! Handling of sink terms for variable ivar
            if (ivar < nbvar ) then
                do jvar = ivar + 1, nbvar
                    if (allocated (FLUX_VAL(ivar, jvar)%idproc)) then
                      do k = 1,nz_max
                       do j=ny_min,ny_max
                        do i=nx_min,nx_max
                           TEND(i,j,k, ivar) = TEND(i,j,k, ivar) - &
                                  FLUX_VAL(ivar, jvar)%val(i,j,k) 
                          enddo
                         enddo
                        enddo
                    endif

                enddo
            endif

            ! Handling of source terms for variable ivar
            do kvar = 1, ivar 
                if (allocated (FLUX_VAL(kvar, ivar)%idproc)) then
                    TEND(:, :, :, ivar) = TEND(:, :, :, ivar) + &
                        FLUX_VAL(kvar, ivar)%val(:, :, :) 
                endif
            enddo
        enddo

#if !defined ECO3M_SUB
       TEND(:,:,:,:) = TEND(:,:,:,:) * dt_bio
#endif

    end subroutine eco3m_compute_sms

!-----------------------------------------------------------------------------------------------

    subroutine eco3m_conc_update
    
!> Subroutine which updates the state variable concentration by using an explicit first order euler method:
!! \f[
!! \C_{t+\delta t} = C_{t} + dt_{bio}\times F(t,C)
!! \f]
!! \author  M Baklouti, V Faure
!! \todo Propose new numerical methods
!-----------------------------------------------------------------------------------------------
  use mod_eco3m_vartypes
  use mod_eco3m
  implicit none

  integer :: ivar,i,j,k

    ! If not coupled
    ! time integration using the explicit Euler method 
    do ivar = 1, nbvar
      do k = 1,nz_max
        do j=ny_min,ny_max
         do i=nx_min,nx_max
        VAR(ivar)%conc(i,j,k) =  VAR(ivar)%conc(i,j,k) + TEND(i,j,k, ivar) 
         enddo
       enddo
      enddo
    enddo

  end subroutine eco3m_conc_update

!-----------------------------------------------------------------------------------------------
#endif
