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
!!MB  potentiellement prioblematique de mettre une telle routine sur git
!--------------------------------------------------------------------------
!> Module that contains functions relative to \f$CO_2\f$
!> This foutine and the underlying model come from PISCES and from James Orr's work
!! \author Rémi Pagès, Nicolas Barrier 

               Module eco3m_co2
!--------------------------------------------------------------------------

#ifdef CO2

    use trc
    use mod_eco3m
    use dom_oce

    use eco3m_files
    use in_out_manager ! I/O manager
    use lib_mpp
    use sms_eco3m

    implicit none

    Integer :: jpdic=-1, jpcal=-1, jptal=-1, jpoxy=-1, jppo4=-1, jpsil=-1

    REAL(wp) :: sco2
    REAL(wp) :: alka0
    REAL(wp) :: oxyg0
    REAL(wp) :: bioma0,sil0   !init silice rp

    REAL(wp) ::  t_oce_co2_flx      !: Total ocean carbon flux
    REAL(wp) ::  t_oce_co2_flx_cum  !: Cumulative Total ocean carbon flux
    REAL(wp) ::  t_atm_co2_flx      !: global mean of atmospheric pco2

    REAL(wp) ::   po4r              !: ???

    !!* Variable for chemistry of the CO2 cycle
    REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ak13       !: ???
    REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ak23       !: ???
    REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   aksp       !: ???
    REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   hi         !: ???
    REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   excess     !: ???
    REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   aphscale   !: 

    INTEGER ::   numonp      = -1           !! Logical unit for namelist pisces output
    LOGICAL :: ln_atint              !: read the AT profile for init pages.r
    CHARACTER(len=34) ::  clname2     !: filename of At init values pages.r

contains

!--------------------------------------------------------------------------
    subroutine co2_find_index

        !! Initialisation of the "PISCES" index from the names 
        !! provided in the TOP namelist. If the names
        !! for the variables are not found, program stop

!--------------------------------------------------------------------------
        integer :: ivar
        CHARACTER (len=20)           :: cltra

        do ivar=1, jp_eco3m
            cltra = TRIM(ctrcnm(ivar))
            if (cltra .eq. "Oxy") then
                jpoxy = ivar
            else if (cltra .eq. "Alk") then
                jptal = ivar
            else if (cltra .eq. "DIC") then
                jpdic = ivar
            else if (cltra .eq. "CaCO3") then
                jpcal = ivar
            else if (cltra .eq. "MID-PO4-P") then
                jppo4 = ivar
            else if (cltra .eq. "Sil") then
                jpsil = ivar
            end if
        end do

        ! check that the index of usefull variables 
        ! are found
        if(jpoxy==-1) then
            write(*,*) "The 'Oxy' var has not been found!!!"
            stop
        end if

        if(jptal==-1) then
            write(*,*) "The 'Alk' var has not been found!!!"
            stop
        end if

        if(jpdic==-1) then
            write(*,*) "The 'DIC' var has not been found!!!"
            stop
        end if

        if(jpcal==-1) then
            write(*,*) "The 'CaCO3' var has not been found!!!"
            stop
        end if
        
        if(jppo4==-1) then
            write(*,*) "The 'PO4' var has not been found!!!"
            stop
        end if
        
        if(jpsil==-1) then
            write(*,*) "The 'Sil' var has not been found!!!"
            stop
        end if
        
    end subroutine co2_find_index

!--------------------------------------------------------------------------

    SUBROUTINE co2_init_norst

        !> initialisation of CO2 related variables
        !! when no restart. Overwrites the initialisation of Eco3M
        !! initialisation.
!--------------------------------------------------------------------------

        NAMELIST/namco2init/sco2, alka0, oxyg0, bioma0,ln_atint,clname2,sil0
        INTEGER :: ios

        REWIND(namelist_ecophy_id)
        READ  (namelist_ecophy_id, namco2init, IOSTAT = ios, ERR = 902 )
        902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisext in eco3m namelist', lwp )

        IF(lwm) WRITE (numonp, namco2init)

       ! trn(:, :, :, jpdic) = sco2
       ! trn(:, :, :, jptal) = alka0
        trn(:, :, :, jpoxy) = oxyg0
        trn(:, :, :, jpcal) = bioma0
        trn(:,:,:,jpsil)    = sil0  
        IF(lwp) THEN                         ! control print
            WRITE(numout,*) ' '
            WRITE(numout,*) ' Namelist parameters for init AT ,namatinit'
            WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
            WRITE(numout,*) '    Choice for reading constant value, ln_atint =',ln_atint
            WRITE(numout,*) ' '
        ENDIF
         IF( .NOT.ln_atint ) THEN
            IF(lwp) THEN                         ! control print
                WRITE(numout,*) '    Constant  value  AT=', alka0
                WRITE(numout,*) ' '
            ENDIF
            trn(:, :, :, jptal) = alka0      ! Initialisation of atmosphericpco2
           WRITE(numout,*) ' '
         ENDIF
END SUBROUTINE co2_init_norst

!--------------------------------------------------------------------------
    integer function sms_eco3m_alloc !* Variable for chemistry of the CO2 cycle

!----------------------------------------------------------------------
        USE lib_mpp , ONLY: ctl_warn
        INTEGER ::   ierr        ! Local variables
        
        ierr = 0

        !* Variable for chemistry of the CO2 cycle
        ALLOCATE( ak13  (jpi,jpj,jpk) ,                              &
            &      ak23(jpi,jpj,jpk)    , aksp  (jpi,jpj,jpk) ,      &
            &      hi  (jpi,jpj,jpk)    , excess(jpi,jpj,jpk) ,       &
            &      aphscale(jpi,jpj,jpk),        STAT=sms_eco3m_alloc )

        IF( sms_eco3m_alloc /= 0 )   CALL ctl_warn('sms_eco3m_alloc: failed to allocate arrays')

    END FUNCTION sms_eco3m_alloc

#else
#endif

end module eco3m_co2
