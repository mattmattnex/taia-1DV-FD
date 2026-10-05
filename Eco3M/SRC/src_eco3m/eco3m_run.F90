!***************************************************************************
!***************************************************************************
!Copyright or © or Copr. CNRS/IRD/Université de la Méditerranée (now
!Aix-Marseille University)
!contributor(s) : Melika BAKLOUTI & Vincent FAURE (10/10/2006)
!
!melika.baklouti@univ-amu.fr; vincent.faure@univmed.fr
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
!---------------------------------------------------------------------------------------
! Main program of the Eco3M modelling tool when used: (i) in uncoupled mode (key NCOUPL) 
! or (ii) in COUPL mode when the coupling is handled by Eco3M 
!
!! The default coupled configuration uses the TKE 1DV physical model 
!! (key key_tke_eco3m). Track this key to see what should be changed for the
!! coupling with another physical model.
!!
!! Other coupled configurations exist (for which the coupling is handled by the
!! physical model):
!! - with the MARS3D physical model
!! - with the NEMO physical model
!!
!! Authors: M. Baklouti, V. Faure, N. Barrier
!!
!! Last revision of the Eco3M platform: M Baklouti, August 2017, April 2022
!---------------------------------------------------------------------------------------
#ifndef key_nemo_eco3m 

!---------------------------------------------------------------------------------------
#if !defined ECO3M_SUB
     program eco3m_main
#else
     subroutine eco3m_main
#endif
!---------------------------------------------------------------------------------------
   use eco3m_string
   use mod_eco3m_vartypes
   use mod_eco3m_files
   use mod_eco3m 
   use mod_eco3m_fprocess 
   use mod_eco3m_irrad
#ifdef key_tke_eco3m
   use comrunmod
   use comdynmod
   use mod_varphy_coupl
#endif 
   implicit none
   real :: t1, t2  ! CPU time

    call CPU_TIME(t1)


#if defined ECO3M_SUB
#ifdef CALC
       ! computes EPAR_z at time tps
        call eco3m_compute_eparz
       ! computes the Chl:C ratio if provided by a mathematical function
       if (CHL_C_BOOL .and. CHL_C_FUNC) then
            call compute_chl_c_ratio
       endif
       ! computes the eco3m source minus sinks (trends)
        call eco3m_compute_sms
#endif

#else /* not ECO3M_SUB */
    ! Initialisation of the Eco3M configuration (INI mode)
    ! In Calc mode, it also initializes the variables
      call eco3m_init_config
    !-- from here, the program is executed only in CALC mode
#ifdef CALC
#ifdef COUPL
    ! In CALC/COUPL mode
    ! Initialisation of the physical model
      call init_phys_model
    ! Time-steppping in coupled mode
      call TRANSFERT_ECO3M_PHYS  
#endif
    ! write depth locations
#if defined M3D_NCOUPL  || (defined COUPL && (! defined ECO3M_SUB))
      call eco3m_write_deptharray
#endif
#ifdef NCOUPL
    !  Main loop on time in uncoupled mode
    Do while (tps <= run_duration*24.d0*3600.d0) 
#endif

#ifdef COUPL
    !  Main loop on time in coupled mode
    Do while (time <= tfin) 
#endif
            
        if (mod(tps, 100*3600.) < 1d-5 )then
           write(*,*)
           write(*,*)'-------------------------------------------'
           write(*,*)'           tps(h) =',tps/3600
           write(*,*)'-------------------------------------------'
           write(*,*)
        endif


#ifdef COUPL
        !  Calling the dynamics of the physical model (the name of the subroutine
        !  depends on the physical model)
#ifdef key_tke_eco3m

        call DYNAMIQUE
#endif

        ! Transferring the outputs of the physical model to Eco3M
        CALL TRANSFERT_PHYS_ECO3M

#endif /* COUPL */
        ! We update the irradiance value (either from the lecture of a file or from a function when not provided by the physical code)
        if (fichirrad(5:) /= "CODEPHYS") then
          if (mod(tps, dt_irrad)==0) then
              call eco3m_irrad_update
          endif
        endif
        ! Computation of Eco3M EPAR_z at time tps
        CALL eco3m_compute_eparz
       
       ! computes the Chl:C ratio if provided by a mathematical function
       if (CHL_C_BOOL .and. CHL_C_FUNC) then
            call compute_chl_c_ratio
       endif

       ! We compute the eco3m source minus sinks (trends)
        call eco3m_compute_sms
       

       ! We update concentration values
#if defined NCOUPL 
       call eco3m_conc_update
#endif        
       ! save outputs     
#if ! defined ECO3M_SUB
        if (mod(tps, dt_save_bio) == 0) then
#ifdef INDIC
       ! calculation of user-defined indicators
             call calc_indicators
#endif 
             call eco3m_write_outputs
#endif
#ifdef SAVE_FLUX
             call eco3m_write_fluxes
#endif
        end if

! transfer of biological variables after the integration
! Necessary because the PREVIOUS function updates the Eco3M trends,
! and the FOLLOWING functions update the PHY trends
#ifdef COUPL
       CALL TRANSFERT_ECO3M_PHYS
#ifdef key_tke_eco3m
       ! Temporal integration of the equation taking into account only source-sink terms due to 
       ! biological processes (explicit euler scheme)
       CALL BIOLOGIE
       ! transfer of updated biological variables to physics 
       CALL TRANSFERT_Eco3M_PHYS


      ! Diffusion, sedimentation and saving of the biological tracer concentrations
      ! is done by the physical model (and the name of the associated routines depend
      ! on it):
    ! Diffusion of biological variables
       call BIODYN
    ! Sedimentation
       call BIOSEDIM 

    ! Sum in order to save the tracer fields
       call BIOSUM

    ! Saving
       call RUNSORTIE
#endif
#endif /* COUPL */
!
#ifdef NCOUPL
       ! time update in uncoupled mode
       tps = tps + dt_bio ! tps = time in s since the beginning of the simulation
#endif

#ifdef COUPL
       ! time update in coupled mode
#ifdef key_tke_eco3m
       time = time + dt  ! time in minutes since the beginning of the simulation
       nstep = nstep+1
       timep = timep + dt  ! absolut time = time + tdebut
       tps = time * 60.0  !  time in seconds
       timeyr = float(int(timep / (xyear*day*hourm))) + 1.
       timeday = timep/(hourm*day)- (timeyr - 1.)*xyear + 1.
       ntimeday = int(timeday + 0.001)
       timeminu = (timeday - ntimeday)*day*hourm
       ntimeminu = int(timeminu + 0.001)
#endif
#endif /* COUPL */
  ENDDO   ! End of the iteration loop

#ifdef NCOUPL
        call eco3m_write_varindex
#endif
    call CPU_TIME(t2)

    write(*,*) "End of program. CPU time is ", (t2-t1)/(60.)

#endif /*-- end of CALC key */
#endif /* ECO3M_SUB */

! deallocation of all global variables.
#if defined CALC
#if !defined ECO3M_SUB 
    call deallocate_globvar_arrays
#endif
#endif


#if !defined ECO3M_SUB
end program eco3m_main
#else
end subroutine eco3m_main
#endif
!------------------------------------------------------------------------------
#ifdef COUPL
#ifdef CALC
    SUBROUTINE TRANSFERT_ECO3M_PHYS

    !> Transfer of the biological variables 
    !! from Eco3M to the TKE model for transport calculations
    !! (sedimentation and sinks)
    !!
    !! - VAR(ii)\%conc(x,y,z)--> TENEUR(z,ii)
    !!
    !! \author Melika Baklouti, Vincent Faure
!------------------------------------------------------------------------------
   Use mod_eco3m
#ifdef key_tke_eco3m
   use mod_varphy_coupl
#endif
   Implicit None
       
 ! local variables
   integer :: ii
   integer :: l, n, m

#ifdef key_tke_eco3m
    do l = 1, nz_max
       do ii = 1, nbvar
            TENEUR(l, ii) = VAR(ii)%conc(1, 1, l)
       enddo
    enddo

#endif
    End SUBROUTINE TRANSFERT_ECO3M_PHYS
#endif
#endif
!---------------------------------------------------------------------------------
#ifdef COUPL
#ifdef CALC
    SUBROUTINE TRANSFERT_PHYS_ECO3M

    !> Transfer of the physical and biological
    !! variables from the TKE model to Eco3M, for biological
    !! processes calculations.
    !!
    !! - TENEUR(z, ii) --> VAR(ii)\%conc(x,y,z)
    !! - TN(z) --> TEMP_BIO(x,y,z)
    !! - SN(z) --> SALT_BIO(x,y,z)
    !! - irrad(x,y) = qsr
    !!
    !! \author Melika Baklouti, Vincent Faure 
!---------------------------------------------------------------------------------
     use mod_eco3m
     use mod_eco3m_irrad, only: irrad,fichirrad
#ifdef key_tke_eco3m
     use mod_varphy_coupl
     use comdynmod
     Implicit None

   ! local variables
     integer :: ii, l, n, m  ! Loop indexes

        do ii = 1,nbvar
            do l = 1, nz_max 
                        VAR(ii)%conc(1, 1, l) = TENEUR(l, ii) 
                        TEMP_BIO(1, 1, l) = TN(l) 
                        SAL_BIO(1, 1, l) = SN(l)
            enddo
        enddo


        ! barrier.n 
        ! In coupled mode, we always take the irradiance from
        ! the physical model
        if (fichirrad(5:) .eq. "CODEPHYS") Then
          do m = ny_min, ny_max
            do n = nx_min, nx_max
                irrad(n,m) = qsr 
            enddo
          enddo
        endif
#endif

    End SUBROUTINE TRANSFERT_PHYS_ECO3M
#endif
#endif
!------------------------------------------------------------------------------
#endif
