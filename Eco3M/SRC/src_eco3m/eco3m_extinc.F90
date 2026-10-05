!***************************************************************************
!***************************************************************************
!Copyright or © or Copr. CNRS/IRD/Université de la Méditerranée (now
!Aix-Marseille University)
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
!> Subroutines for the calculation of light irradiance and extinction 
! \author Melika Baklouti, Vincent Faure, Nicolas Barrier

!--------------------------------------------------------------------------------------------------------------    
    subroutine eco3m_irrad_init

    !> Initialization of the light-related variables. It first reads the
    !! extinction parameters provided in the namelist (i.e. if irradiance is provided by a file 
    !! or computed by a function). It also allocates the extinction-related variables.
    !!
    !!  - If RGB: it initialises the RGB module.
    !!  - If not RGB: in INI mode, it writes the calc_extinc.inc file, which contains the function
    !!  for the calculation of light extinction coefficients. 
    !! \author Melika Baklouti, Vincent Faure
    !! \author Nicolas Barrier

        use mod_eco3m
        use mod_eco3m_irrad
        use mod_eco3m_files
        implicit none
        
        ! local variables
#if !defined M0D_NCOUPL && !defined RGB
        integer :: k  ! Loop index
#endif 
        logical :: file_exist
        integer :: istat
        namelist/nam_eco3m_irrad/fichirrad,filename_irrad, irrad_max, dt_irrad, irr2par, &
    &            albedo, irr_param

        ! Reading the namelist for irradiance parameters
        REWIND(namelist_eco3m_id)
        READ  (namelist_eco3m_id, nam_eco3m_irrad)
#if defined CALC &&  !defined M0D_NCOUPL       
        if ((pos_nzmax.ne."BOTT") .and. (pos_nzmax.ne."SURF")) then
            write(*,*) "*********************************"
            write(*,*) "The pos_nzmax variable must either be BOTT or SURF"
            write(*,*) "Currently, it is ", pos_nzmax
            write(*,*) "Please correct the namelist"
            write(*,*) "*********************************"
            stop
        endif
#endif
#ifdef RGB
        ! If RGB key is used, initialisation of the rgb_ex_coef 
        ! matric read from the files
        call eco3m_rgb_init 
#else
        ! If no RGB key,
        ! reads extinction-related settings in the config.ini files
#if !defined M0D_NCOUPL        
        ! No light extinction by Chl in 0D mode
        call fscanf_extinc        
#ifdef INI
        ! Writes in calc_extinc.inc  file used for the 
        ! light extinction coef. calculation
        call fwritef_extinc_inc
#endif /* INI */
#endif /* !M0D_NCOUPL */
#endif /* !RGB  */

#ifdef CALC
        ! Allocation of 2D irradiance and EPAR0+
        Allocate (irrad(nx_min:nx_max,ny_min:ny_max),stat=istat)
        if (istat /=0) write(*,*) "problem with the allocation of irrad"

        ! Allocation of the dz array
        Allocate(dz(nx_min:nx_max,ny_min:ny_max,1:nz_max),stat=istat)
        if (istat /=0) write(*,*) "problem with the allocation of dz"

        ! Surface EPAR
        Allocate(E_PAR(nx_min:nx_max,ny_min:ny_max),stat=istat)
        if (istat /=0) write(*,*) "problem with the allocation of e_par"

        ! Vertical profile of EPAR
        Allocate (E_PARZ(nx_min:nx_max,ny_min:ny_max,1:nz_max),stat=istat)
        if (istat /=0) write(*,*) "problem with the allocation of e_parz"
        
        ! Allocation of the extinction coefficient
        Allocate(kextinc(nx_min:nx_max,ny_min:ny_max,1:nz_max),stat=istat)
        if (istat /=0) write(*,*) "problem with the allocation of kextinc"

        ! If irradiance has to be read from a file
        if (fichirrad(5:) .eq. "FILE") then
           filename_irrad=trim(adjustl(eco3m_root_dir))//trim(adjustl(filename_irrad))
           inquire(file=filename_irrad, exist=file_exist)
           if (.not.(file_exist)) then
              write(*,*) "****************************************"
              write(*,*) "The ", trim(filename_irrad), " file does not exist"
              write(*,*) "This program will stop"
              write(*,*) "****************************************"
              stop
           end if

          open(file_irrad_id, file=filename_irrad)
        end if
#endif

        ! WRITING LOG
#ifdef CALC
        write(file_CR_id,*) "================== Eco3M: Initialization of the light extinction module"
        write(file_CR_id,*) "Surface irradiance calculation mode: ", trim(fichirrad)
        select case(trim(fichirrad))
        case ("PAR_FONCTION")
            write(file_CR_id,*) "PAR computed by a function with maximum PAR = ", irrad_max
        case ("IRR_FONCTION")
            write(file_CR_id,'(A,F10.5)') "IRR computed by a function with maximum IRR = ", irrad_max
            write(file_CR_id,"(A, F10.5)") "Irradiance to PAR conversion factor: ", irr2par
#ifdef COUPL
        case("IRR_CODEPHYS")
            write(file_CR_id,*) "IRR provided by the physical model"
            write(file_CR_id,"(A, F10.5)") "Irradiance to PAR conversion factor: ", irr2par
        case("PAR_CODEPHYS")
            write(file_CR_id,*) "PAR provided by the physical model"
#endif
        case ("PAR_FILE")
            write(file_CR_id,*) ", PAR provided by the ", trim(fichirrad_long), " file"
        case ("IRR_FILE")
            write(file_CR_id,*) ", IRR provided by the ", trim(fichirrad_long), " file"
            write(file_CR_id,"(A, F10.5)") "Irradiance to PAR conversion factor: ", irr2par
        case default
            write(file_CR_id,*) "Bad setting -> Set to IRR_FUNCTION"
            fichirrad="IRR_FUNCTION"
            irrad_max = 450 ! W/m2
            write(file_CR_id,"(A,F10.5)") ", IRR computed by with maximum irradiance = ", irrad_max
        end select
!
        write(file_CR_id,*) 
        write(file_CR_id,'(A, F10.5)') "Albedo: ", albedo
        write(file_CR_id,*)
        ! light extinction is not relevant in 0D mode
#if !defined M0D_NCOUPL
#ifdef RGB
        write(file_CR_id,*) "Extinction calculation using RGB formulation"
        write(file_CR_id,*) 
#else
        write(file_CR_id,*) "Extinction calculation: ", trim(proc_mod(irr_par%idproc)%nomsub), " function"
        write(file_CR_id,*) "Parameters: "
        do k = 1, size(proc_mod(irr_par%idproc)%nompar)
            write(file_CR_id,"(A,A,F10.5)") proc_mod(irr_par%idproc)%nompar(k), &
                " = ", irr_par%valpar(k) 
        end do
        write(file_CR_id,*) 
#endif /*RGB*/
#endif /*!defined M0D_NCOUPL*/
#endif /*CALC*/
    end subroutine eco3m_irrad_init

#ifdef CALC
!-------------------------------------------------------------------------------------------------------------------
    SUBROUTINE eco3m_irrad_update

    !> **Only if NCOUPL key is on**. Update of surface irradiance values every 
    !! *dt_irrad*.
    !!
    !!  - Either computed by a mathematical function written in this subroutine 
    !!  - Either read from the *filename_irrad* file 
    !!   (the file, which ID is defined in the *mod_eco3m_files* module, is opened when
    !!   *nbcallbio == 0*)
    !! 
    !! \date 2008-06-26
    !! \author Melika Baklouti
!-------------------------------------------------------------------------------------------------------------------
        use mod_eco3m, only: tps,pi
        use mod_eco3m_irrad
        use mod_eco3m_files, only:filename_irrad,file_irrad_id
        implicit none

        ! Local variables
        Integer :: errlec
        real :: irradscal

        if (fichirrad(5:) /= 'FONCTION') Then ! irr is read in a file which name is given in the Eco3M namelist

            ! Reads irradiance
            Read(file_irrad_id, *, iostat=errlec) irradscal
            if (irradscal < 0.d0) irradscal = 0.d0
            if (errlec /= 0) then 
            ! case where the end of file is reached
              rewind (file_irrad_id)
              Read(file_irrad_id, *, iostat=errlec) irradscal
              if (errlec /= 0) then
                write(*,*) 'Problem with the reading of the  ', trim(filename_irrad), ' file.'
                write(*,*) "This program will be stopped"
                stop
              endif
            endif
        else ! surface irradiance or par is provided by a periodic function (irrad_MAX provided in the namelist) 
            irradscal = irrad_MAX*exp(3.7*(cos(2.0*pi*tps/24/3600.)-1))
            ! periodic function for surface irradiance with seasonal variation
            !irradscal = max(1.d0,max(50.,ABS(800.*cos(0.45*pi*((tps/3600.-4000.)/4100.)))) * sin(pi*tps/12./3600.))
        endif

        ! The surface irradiance is used everywhere on the domain
        irrad = irradscal

    End Subroutine eco3m_irrad_update
!----------------------------------------------------------------------------------------------------------------------
#endif

#ifdef CALC
!-------------------------------------------------------------------------------------------------------------------
    Subroutine eco3m_compute_eparz
 
    !> **Only if the M0D_NCOUPL key is not activated.** Calculation of the vertical 
    !! profiles of Photosynthetic Available Radiation (PAR).
    !!
    !!  - If the *irrad* array is irradiance, it first converts irradiance
    !!  into \f$E_{PAR}^{0+}\f$ by using \f$E_{PAR}^{0+}= irrad\times irr2par\f$
    !!  - It then converts \f$E_{PAR}^{0+}\f$ into \f$E_{PAR}^{0-}\f$ by using
    !!  \f$E_{PAR}^{0-}=(1-albedo) \times E_{PAR}^{0+}\f$
    !!  - It calculates the total chlorophyll by using the *compute_chl_tot* subroutine (*eco3m_chl* module)
    !!  - If RGB
    !!      * It calls the eco3m_rgb::compute_eparz_rgb() subroutine 
    !!  - If no RGB
    !!      * It computes the light extinction coefficients using the function defined
    !!      in the config.ini file (through a call to the **calc_extinc.inc** file)
    !!      * It calculates PAR values as a function of depth:  E_PARZ 
    !! \author Melika Baklouti, Vincent Faure
    !! \date 2008-06-26
!-------------------------------------------------------------------------------------------------------------------
        use mod_eco3m
        use mod_eco3m_irrad
        use mod_eco3m_fprocess
        implicit none

        ! Local variables

        Integer    :: i, j, k  ! Loop indexes
        character(3) :: debfichirrad
#if !defined M0D_NCOUPL
        Real(8):: Chl_tot(nx_min:nx_max,ny_min:ny_max,1:nz_max)
#endif
        ! Calculation of the EPAR0+, except if the input is already EPAR (instead of irrad)
        fichirrad = adjustl(fichirrad)
        debfichirrad= fichirrad(1:3)

        ! Conversion from irr to PAR if necessary using the irr2par scale factor  (0.43 usually)
        if (debfichirrad /= 'PAR') then 
            E_PAR(:,:) = irr2par * irrad(:,:) 
        else
            E_PAR(:,:) =  irrad(:,:)    
        endif

#if defined NCOUPL
        ! E_PAR(0-) is computed by taking into account the albedo and 
        ! the irr_param parameter
        E_PAR(:,:) = irr_param * (1.-albedo) * E_PAR(:,:)  
#else
        ! Case where the irradiance is already corrected from albedo
        ! See dynfluxsur.f: qsr=(1-albedo_phy)*qsr
        E_PAR(:,:) = irr_param * E_PAR(:,:) 
#endif 

#if !defined M0D_NCOUPL
        ! Calculation of total Chl, only used for light extinction
        call compute_chl_tot(Chl_tot)
         

!====  extinction coefficients        
#ifdef RGB
         call compute_rgb_kextinc(Chl_tot)
#else
        ! the formulation for the calculation of extinction is provided in the file "calc_extinc.inc" 
        ! which is automatically generated (see the subroutine fwritef_extinc_inc)
        include "calc_extinc.inc"
#endif


!===== extinction due to Suspended Material (if SM is provided)
#ifdef MUSTANG
         call compute_chl_mes_kextinct() 
#endif

!===== calculation of EPARZ 
#ifdef RGB
        call compute_eparz_rgb(Chl_tot)
#else
        ! Calculation of EPARZ if nz_max is at the surface
        if (pos_nzmax == 'SURF') then 
            do k = nz_max, 1, -1
                do j = ny_min, ny_max
                    do i = nx_min, nx_max
                        if (k == nz_max) then
                            E_PARZ(i, j, k) = E_PAR(i, j)*exp(-0.5*kextinc(i, j, k)*abs(dz(i, j, k)))
                        else
                            E_PARZ(i, j, k) = E_PARZ(i,j,k+1)*exp(-0.5*(kextinc(i, j, k+1)*abs(dz(i, j, k+1))+&
                                kextinc(i, j, k)*abs(dz(i, j, k))))
                        endif

                    enddo
                enddo
            enddo

            ! Calculation of EPARZ if nz_max is at the bottom
        elseif (pos_nzmax == 'BOTT') then
            do k = 1, nz_max
                do j = ny_min, ny_max
                    do i = nx_min, nx_max
                        if (k==1) then 
                            E_PARZ(i,j,k) = E_PAR(i,j)*exp(-0.5*kextinc(i,j,k)*abs(dz(i,j,k)))
                        else 
                            E_PARZ(i,j,k) = E_PARZ(i,j,k-1)*exp(-0.5*(kextinc(i,j,k-1)*abs(dz(i,j,k-1)) + &
                                kextinc(i,j,k)*abs(dz(i,j,k))))  
                        endif
                    enddo
                enddo
            enddo
        endif
#endif

#else
        ! In OD, there is no light extinction (EPARZ == EPAR at the surface)
        do k=1,nz_max
            E_PARZ(:,:,k) = E_PAR(:,:)
        enddo
#endif /* !defined M0D_NCOUPL */


    End Subroutine eco3m_compute_eparz
!------------------------------------------------------------------------------------------------------------------
#endif   

#ifdef INI
!------------------------------------------------------------------------------------------------------------------
    Subroutine fwritef_extinc_inc

    !> **Only if the INI key is on and the RGB key is off**. Automatic filling of the .inc file dedicated to the
    !! calculation of light extinction within the Eco3M model.
    !!
    !! \author Melika Baklouti
    !! \author Vincent Faure
    !! \date 31/08/2017
!------------------------------------------------------------------------------------------------------------------
    use mod_eco3m
    use mod_eco3m_irrad
    use mod_eco3m_files, only:file_extincinc_id
    implicit none

    !-- Local variables 

    Integer             :: idtemp
    integer             :: kk,dim_param
    character(L_VAR_LG)  :: nomsb

    write(file_extincinc_id,*) "! This file is automatically generated by the fwritef_extinc_inc subroutine"
    write(file_extincinc_id,*) "! It contains the include file for the calculation of light extinction coefficients"
    write(file_extincinc_id,*) "! It is generated by the settings in the config.ini file"
    write(file_extincinc_id,*) " " 
    write(file_extincinc_id,"(A)")"kextinc = &"

    idtemp = IRR_PAR%idproc
    write(nomsb,"(A)") PROC_MOD(idtemp)%nomsub
    nomsb=trim(adjustl(nomsb))
    dim_param = PROC_MOD(idtemp)%nbpar

    write(file_extincinc_id,*) nomsb(1:len_trim(nomsb)),'(& '
    do kk = 1, dim_param 
        write(file_extincinc_id,*)'     IRR_PAR%valpar(',  kk  ,'),  & '
    enddo
    write(file_extincinc_id,*)'Chl_tot) '

    End Subroutine fwritef_extinc_inc
!---------------------------------------------------------------------------------------------------
#endif
!---------------------------------------------------------------------------------------------------
!  No light absorption by chlorophyll in 0D configurations
#ifndef M0D_NCOUPL

    subroutine fscanf_extinc

    !> **Only if the RGB key is not set**. Reads the settings of of the biological 
    !! process dedicated to the calcultion of light extinction coefficient in the *config.ini* file.
    !! \author Melika Baklouti
!--------------------------------------------------------------------------------------------------
! global variables
    use mod_eco3m         
    use mod_eco3m_irrad
    use mod_eco3m_files, only:file_config_id
    use mod_eco3m_id_extract, only : f_proc2id
    implicit none

! local variables
    character(l_chain) :: chaine, tempo, tempo2, tempopo
    integer :: errlec,istat
    integer :: idtemp
    integer :: k
    integer :: nb_args

        do 
            Read(file_config_id,*,iostat=errlec) chaine
            if (errlec /= 0) stop 'pb with the reading of the light extinction filename'
            if (chaine(1:1) == '#') cycle

            ! name of the process
            tempo = f_chain(chaine, 1, '(')
            tempo = trim(adjustl(tempo))

            ! tempo2 contains the string within the brackets
            tempo2 = chaine(len_trim(tempo)+2:len_trim(chaine)-1)

            ! idtemp is the number associated with the "tempo" function 
            idtemp = f_proc2id(tempo)
            IRR_PAR%idproc = idtemp 
            nb_args= PROC_MOD(idtemp)%nbpar

            IF (nb_args >0) THEN   !MB2021
            Allocate(IRR_PAR%valpar(1:nb_args),stat=istat)
            if (istat/=0) write(*,*)"pb with the allocation of IRR_PAR%valpar"

            ! The parameter values are read from the file
            ! and put into the IRR_PAR%valpar(k) variable
            do k= 1, nb_args 
                 tempopo = f_chain(tempo2, k, '>')
                 Read(tempopo,*) IRR_PAR%valpar(k)
            enddo
            ENDIF !MB2021
            exit
        enddo

    end subroutine fscanf_extinc
#endif
!----------------------------------------------------------------------------------------------
!!MB: under construction 
!     to implement when ready (coupled with CROCO and MUSTANG)
!#ifdef MUSTANG
!    subroutine compute_chl_mes_kextinc()
!
!    ! use module du coupleur : cmes_3dmgl        
!#ifdef RGB
!    ! calcul du kextinct pour chaque longeur d'onde 
!    ! boucle sur les 4 indices (i,j,k,bandes)
!    ! do irgb = 1, rgb_nbands
!    ! kextinc_rgb(irgb, i, j, k)=kextinc_rgb(irgb, i, j, k)+cmes_3dmgl*coeff d'extinct
!#else
!    ! kextinc(i, j, k)=kextinc(i, j, k)+cmes_3dmgl*coeff d'extinct
!#endif
!
!    end subroutine compute_chl_mes_kextinc()
!#endif
!--------------------------------------------------------------------------------------------------

