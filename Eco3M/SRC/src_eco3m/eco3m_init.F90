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
!-------------------------------------------------------------------------------------------------
    subroutine eco3m_init_config

    !> Initialises the  Eco3M simulation:
    !!
    !! - Opens Eco3M input files (namelist, model.def, config.ini, etc): 
    !! eco3m_files::eco3m_files_init()
    !! - Initializes run parameters (namelist): eco3m_runpar::eco3m_runpar_init()
    !! - Reads modele.def file and initializes process definition: eco3m_mod::eco3m_mod_init()
    !! - Initializes extinction (allocation of tables, reading of config.ini
    !! file, eventually writing into the calc_extinc.inc file): eco3m_irrad::eco3m_irrad_init
    !! - Initializes Eco3M state variables (allocation of tables, reading of 
    !! config.ini file, couting of elements, allocation and 
    !! definition of index arrays): eco3m_var::eco3m_var_init()
    !! - Initializes fluxes (allocation of tables, reading of
    !! config.ini file, eventually writting into the call.inc file): eco3m_flux::eco3m_flux_init()
    !! - Initializes chlorophyll variables (reading of the config.ini file, allocation
    !! of arrays, eventually writting of the calc_ChlC.inc) file: eco3m_chl::eco3m_chl_init()
    !! - Initializes  \f$CO_2\f$ related variables (if key
    !! **MODULE_CO2**)
    !! - Initializes the time-stepping variables *tps*  
    !!  (if **CALC** key): eco3m_step::eco3m_step_init()
    !! - Initializes the vertical profiles for variable concentration (if **CALC** key) and
    !!  depth array (if **NCOUPL**): eco3m_init_values()
    !! - Initializes the output files (reading into the namelist for output parameters, opening 
    !!  output files) (if **CALC** key): eco3m_outputs::eco3m_outputs_init()
    !! - Closing of input files which are not used anymore after 
    !! that (namelist, etc): eco3m_files::eco3m_files_close()
    !!
    !! \author Melika Baklouti, Vincent Faure, Nicolas Barrier
   
!-- global variables 
      use mod_eco3m_vartypes
      use mod_eco3m
      use mod_eco3m_files
      use mod_eco3m_irrad
#ifdef CALC
      use mod_eco3m_outputs
#endif
      implicit none

      character(124) :: filename

        ! opens the configuration and the .log files 
        call eco3m_files_init  
#ifdef MODTEST
        write(*,*) "eco3m_files_init OK"
#endif

        ! reads the run parameters in the namelist
#ifdef CALC
        call eco3m_runpar_init
#ifdef MODTEST
        write(*,*) "eco3m_runpar_init OK"
#endif
        
        filename = trim(adjustl(eco3m_root_dir))//"OUTPUTS/"//trim(adjustl(run_name))//".log" 
        open(file_CR_id, file=trim(filename), status='replace')
        Write(file_CR_id,*) "===================== Eco3M configuration used: " 
        Write(file_CR_id,*) "Name of the configuration of the biogeochemical model: ", bioconf_name
        Write(file_CR_id,*) " "
#endif /* CALC */        
 ! reads the modele.def file; initializes the PROC_MOD matrix
        call eco3m_mod_init  
#ifdef MODTEST
        write(*,*) "eco3m_mod_init OK"
#endif


 ! reads the variable-related parameters in the config.ini file; allocates the VAR array
 ! and of the VAR%conc arrays for each variable
        call eco3m_var_init   
#ifdef MODTEST
        write(*,*) "eco3m_var_init OK"
#endif
 
 ! reads namelist in INI (for SAVE_FLUX);
 ! opens all the output files (flux, variables, physical variables)
 ! before flux_init call 
        call eco3m_outputs_init
#ifdef MODTEST
        write(*,*) "eco3m_outputs_init OK"
#endif

 ! initializes the FLUX_VAL and FLUX_PAR arrays by reading the
 ! config.ini file (dynamical allocation, etc.)
        call eco3m_flux_init
#ifdef MODTEST
        write(*,*) "eco3m_flux_init OK"
#endif

#if defined INI && defined INDIC  
       call eco3m_indic_init
#ifdef MODTEST
        write(*,*) "eco3m_flux_init OK"
#endif
#endif

 ! initializes chlorohyll variables in the config.ini file;
 ! dynamical allocation
        call eco3m_chl_init
#ifdef MODTEST
        write(*,*) "eco3m_chl_init OK"
#endif
        
 ! reads the parameters related to irrad in the namelist, 
 ! and those related to the extinction function in the config.ini
 ! file; dynamical allocation of the extinction related variables
        call eco3m_irrad_init   
#ifdef MODTEST
        write(*,*) "eco3m_irrad_init OK"
#endif

#if !defined ECO3M_SUB || defined ECO3M_1DV
#ifdef MODULE_CO2
     ! initializes CO2-related variables
        call eco3m_co2_init
#endif
#ifdef CALC        
     ! Assignment of initial conditions 
        call eco3m_init_values
#ifdef MODTEST
        write(*,*) "eco3m_init_values OK"
#endif
#endif /* CALC */
#endif /* ECO3M_SUB */

#ifdef CALC        
        ! Initialisation of source minus sinks (TEND variable)
        call eco3m_sms_init
#ifdef MODTEST
        write(*,*) "eco3m_sms_init OK"
#endif
#endif

#if (defined CALC) 
        write(file_CR_id,*) "=============== Eco3M: Initialisation of the",&
            &" run configuration ==============="
        write(file_CR_id,*) "Run name: ", run_name
#if (defined NCOUPL) || (defined COUPL && !defined ECO3M_SUB)
        write(file_CR_id,*) "X-dimensions : ", nx_min, nx_max
        write(file_CR_id,*) "Y-dimensions : ", ny_min, ny_max
        write(file_CR_id,*) "Z-dimensions : ", 1, nz_max
        write(file_CR_id,*) "Run duration (days): ", run_duration
        write(file_CR_id,*) "Biological time step (s): ", dt_bio
        write(file_CR_id,*) "Biological saving time step (s): ", dt_save_bio
#endif
#if (!defined M0D_NCOUPL)
        write(file_CR_id,*) "Position of the pos_nzmax level: ", pos_nzmax
#endif
        write(file_CR_id,*) 
#endif

        ! Closing of all the input files
        call eco3m_files_close

    end subroutine eco3m_init_config

!-------------------------------------------------------------------------------------------------------
    subroutine eco3m_files_init

 !> Subroutine that opens all the configuration files used in Eco3M (namelist(s),  config.ini and
 !! model.def files), but also the .inc files in INI mode.
 !! \note Output files are opened in the eco3m_outputs::eco3m_outputs_init() function
 !! \date 2015-12-23
 !! \author Melika Baklouti, Vincent Faure, Nicolas Barrier
!-------------------------------------------------------------------------------------------------------
        use mod_eco3m 
        use mod_eco3m_files
        implicit none

        namelist/nam_eco3m_config/eco3m_root_dir,filename_config,filename_confmod

!  local variables
        integer :: errlec
        logical :: file_exist
        character(len=L_CHAIN) :: chaine
        character(len=lfl) ::  filename

! opens the Namelist specific to Eco3M
#if defined COUPL && defined ECO3M_SUB
        inquire(file="namelist_eco3m", exist=file_exist)
#else
        inquire(file="../CONFIG_ECO3M/namelist_eco3m", exist=file_exist)
#endif

        if (.not.(file_exist)) then
            write(*,*) "****************************************"
            write(*,*) "The namelist_eco3m file does not exist"
            write(*,*) "This program will stop"
            write(*,*) "****************************************"
            stop
        end if

#if defined COUPL && defined ECO3M_SUB
        open(namelist_eco3m_id, file="namelist_eco3m", status='old', &
            form='formatted', access="sequential")
#else
        open(namelist_eco3m_id, file="../CONFIG_ECO3M/namelist_eco3m", status='old', &
            form='formatted', access="sequential")
#endif

#ifdef key_nemo_eco3m
! opens an additional namelist necessary for the coupling between NEMO and Eco3M
        inquire(file="namelist_nemo_eco3m", exist=file_exist)
        if (.not.(file_exist)) then
            write(*,*) "****************************************"
            write(*,*) "The namelist_nemo_eco3m file does not exist"
            write(*,*) "This program will stop"
            write(*,*) "****************************************"
            stop
        end if

        open(namelist_ecophy_id, file="namelist_nemo_eco3m", status='old', &
            form='formatted', access="sequential")
#endif
  ! Reads information on the configuration of the biogeochemical model 
        rewind(namelist_eco3m_id)
        read(namelist_eco3m_id, nam_eco3m_config)
#if defined INI
    eco3m_root_dir   = "../"
#endif



  ! opens the config.ini file
        filename_config=trim(adjustl(eco3m_root_dir))//trim(adjustl(filename_config))
        inquire(file=filename_config, exist=file_exist)
        if (.not.(file_exist)) then
            write(*,*) "****************************************"
            write(*,*) "The ", trim(filename_config), " file does not exist"
            write(*,*) "This program will stop"
            write(*,*) "****************************************"
            stop
        end if

        open(file_config_id, file=filename_config)

! reads the name of the configuration of the biogeochemical model in the config.ini file
        errlec = 0
        do 
            Read(file_config_id,*,iostat=errlec) chaine
            if (errlec /=0 ) stop 'pb with the reading of the number of compartments'
            if (chaine(1:1) == '#') cycle
            exit
        enddo
        bioconf_name = trim(adjustl(chaine))
        bioconf_name="Eco3M_"//bioconf_name

 ! opens the log file for the biogeochemical model 
#ifdef INI
        filename = trim(adjustl(eco3m_root_dir))//"OUTPUTS/"//trim(bioconf_name)//".log" 
        open(file_CR_id, file=trim(filename), status='replace')
        Write(file_CR_id,*) "===================== Eco3M configuration used: " 
        Write(file_CR_id,*) "Name of the configuration of the biogeochemical model: ", bioconf_name
        Write(file_CR_id,*) " "
        ! In INI mode, opens the .inc files
        open(file_extincinc_id, file="calc_extinc.inc", status="replace")
#endif 

    end subroutine eco3m_files_init

!-----------------------------------------------------------------------------------------
    subroutine eco3m_runpar_init
    !> Initialisation of Eco3M run parameters. They are read from the
    !! namelist_eco3m file. 
    !! \author Nicolas Barrier
!-----------------------------------------------------------------------------------------
#if ! defined M3D_NCOUPL
    use mod_eco3m_files, only: namelist_eco3m_id
#else
    use mod_eco3m_files, only:namelist_eco3m_id,filename_depth
#endif
    use mod_eco3m
    implicit none

#if (defined NCOUPL) && (defined M0D_NCOUPL) 
   namelist/nam_eco3m_run/ run_name, run_duration, dt_bio, dt_save_bio, &
            nx_min, nx_max, ny_min, ny_max, nz_max
#else /* ALL NOT 0D CONFIGURATIONS */
   namelist/nam_eco3m_run/ run_name,dt_bio,pos_nzmax
#if   (defined NCOUPL && defined M3D_NCOUPL)  
   namelist/nam_eco3m_run_main/run_duration, dt_bio, dt_save_bio, &
            nx_min, nx_max, ny_min, ny_max, nz_max, filename_depth
#elif (defined COUPL && ! defined ECO3M_SUB )
   namelist/nam_eco3m_run_main/run_duration, dt_bio, dt_save_bio, &
            nx_min, nx_max, ny_min, ny_max, nz_max
#endif
#endif /* end of the test on OD or NOT OD configuration */
        ! reads the namelist run parameters
        REWIND(namelist_eco3m_id)
        READ  (namelist_eco3m_id, nam_eco3m_run)
#if (defined NCOUPL && defined M3D_NCOUPL) || (defined COUPL && ! defined ECO3M_SUB)
        READ  (namelist_eco3m_id, nam_eco3m_run_main)
#endif
        !-- conversion in seconds:
        dt_save_bio = dt_save_bio * 60.0

#ifdef NCOUPL       
    ! time and loop counter initialisation 
    tps = 0
#endif
    end subroutine eco3m_runpar_init

!-----------------------------------------------------------------------------------------
#ifdef CALC
    subroutine eco3m_init_values

    !> Initialization of the vertical profiles. 
    !!
    !! - Reads the namelist parameters relative to initial conditions (*ic_input_directory*, 
    !! *default_conc*, etc)
    !! - Loops over all state variables:
    !!      - if the input concentration file exists for the variable, the data are read
    !!      - if the input concentration file does not exist, the concentration data are 
    !!        initialised to the *default_conc* namelist value
    !! - If **NCOUPL** key is on:
    !!     - Initializes temperature and salinity to the *default_temp* and 
    !!        *default_salt* namelist values, respectively
    !!     - Reads the input depth file (*depth_z* array).
    !      - Initializes the *bathymetry* array (total depth) depth as 
    !!        the maximum of the *depth_z* array
    !!
    !! \author Melika Baklouti, Vincent Faure, Nicolas Barrier
    !! \todo Possibility to read a vertical salinity/temperature profile 
!-----------------------------------------------------------------------------------------
!
! global variables
   use mod_eco3m
   use mod_eco3m_files
   implicit none


! local variables
   real(8) :: default_conc  ! default concentration if no input file is provided
   real(8) :: default_temp  ! default temperature in uncoupled mode (i.e. key COUPL not activated)
   real(8) :: default_salt  ! default salinity in uncoupled mode (i.e. key COUPL not activated)
   character(256) :: ic_input_directory  ! directory for initial conditions
   logical :: file_exist ! boolean to check the existence of input file
   integer :: ivar  ! index of the variable
   character(300) :: filename  ! name of the file
   
   namelist/nam_eco3m_ic/ic_input_directory, default_conc, default_temp, default_salt
#if ! defined M3D_NCOUPL
   namelist/nam_eco3m_run/run_name, run_duration, dt_bio, dt_save_bio, &
            nx_min, nx_max, ny_min, ny_max, nz_max, pos_nzmax
#else
   namelist/nam_eco3m_run/run_name, run_duration, dt_bio, dt_save_bio, &
            nx_min, nx_max, ny_min, ny_max, nz_max, pos_nzmax,filename_depth   
#endif

        ! Reads the namelist parameters relative to initial conditions
        rewind(namelist_eco3m_id)
        read(namelist_eco3m_id, nam_eco3m_ic)
        
        write(file_CR_id,*) "================ Eco3M:  Initialisation of the Eco3M model concentrations "

        ! initializes the state variable concentrations
        ! loops over the variables, checking the existence of the file for the
        ! initialisation of each state variable 
        do ivar = 1, nbvar

            VAR(ivar)%conc = default_conc  ! initialisation of the concentration to the default one

            ! construction of the name of the initial condition file
            filename = trim(adjustl(eco3m_root_dir))//trim(adjustl(ic_input_directory))&
                //trim(VAR(ivar)%comp) &
                // '_' // trim(VAR(ivar)%scomp) // '_' // trim(VAR(ivar)%elmt) // ".DAT"
            ! we inquire for the existence of the file. 
            inquire(file=trim(filename), exist=file_exist)
            
            if(file_exist) then
                ! if it exists, we call the subroutine that reads the file
                write(file_CR_id,"(A, I3, A, A)") "Initial conditions for variable ", ivar, &
                    " are read from the file ", trim(filename)
                call init_prof(var(ivar)%conc, filename)
            else
                write(file_CR_id,"(A, A, A, I3, A, A ,E15.3)") "The file ", trim(filename), &
                    " does not exist. Concentration for variable ",ivar, " will be",&
                    &" initialised by using the default concentration ", default_conc
            end if 
        end do
        write(file_CR_id,*)

#ifdef NCOUPL
        ! in uncoupled mode, the temperature and salinity values are constant
        ! and initialised through a namelist parameter
        TEMP_BIO = default_temp
        SAL_BIO = default_salt
#ifdef M3D_NCOUPL
        ! we inquire for the existence of the file with cell depths
        filename_depth=trim(adjustl(eco3m_root_dir))//trim(adjustl(filename_depth))
        inquire(file=trim(filename_depth), exist=file_exist)
        if(.not.(file_exist)) then
            write(*,*) "****************************************************"
            write(*,*) "The ", trim(filename_depth), " file does not exist"
            write(*,*) "This program will stop"
            write(*,*) "****************************************************"
        end if
        
        ! in 3D uncoupled mode,  depth_z is initialized by reading an input file.
        ! And the prof array is set to 0
        write(*,"(2A)") "Depth at W levels are read in the file ", trim(filename_depth)

        ! reads of the depths file and filling the depth_z array (depth of W levels)
        call init_prof(depth_z, filename_depth)

        ! computes the bathymetry and dz, moving depth_z to T levels (overwritting depth_z)
        call init_dz_array
#endif /* M3D_NCOUPL */
#endif /* NCOUPL */

    end subroutine eco3m_init_values
#endif

!-----------------------------------------------------------------------------------------
#ifdef CALC

    subroutine init_prof(val, filename)

    !> Initialisation of the vertical profile of a given variable by reading an
    !! input file. The latter must have *nz_max* lines
    !! and *2* columns: the first one with the k-index (integers), the second one with the input values (real).
    !! \author Nicolas Barrier
    !! \warning Contrary to the previous version, the data are not interpolated in the program! This must
    !! be done by the use of an external software (excel, python, matlab, etc)
!-----------------------------------------------------------------------------------------
        ! global variables
        use mod_eco3m
        use mod_eco3m_files, only:file_input_id,lfl
        implicit none

        ! Declaration of local variables
        real(8), intent(inout) :: val(nx_min:nx_max, ny_min:ny_max, nz_max)  !< Variable to read in the file
        character(lfl), intent(in) :: filename  !< Name of the file to read

        integer :: rstatus ! End of file flag
        integer :: i,j,k  
        integer :: errlec
        real(8) :: conc  ! Concentration read from the file
        integer :: numoflines  ! Number of lines within the file
        logical :: file_exist  ! Boolean for checking of input file
        character(120)::test

        ! initialisation of the rstatus value
        ! at 0 allows entering the loop
        rstatus = 0 

        ! inquires for the existence of the input file. 
        inquire(file=trim(filename), exist=file_exist)
        if(.not.(file_exist)) then
            write(*,*) "**********************************************"
            write(*,*) "The file ", trim(filename), " does not exist"
            write(*,*) "This program will be stopped"
            write(*,*) "**********************************************"
            stop
        end if

        ! opens the input file 
        open(file_input_id, file=trim(filename))

        ! first reading to count the number of lines
        rstatus = 0
        numoflines = 0
        do while (rstatus == 0)
            read(file_input_id, '(A)', iostat=rstatus) test
            numoflines = numoflines + 1
        end do

        ! correction of the number of lines
        numoflines = numoflines - 1

        ! checks for file consistency with nz_max
        if (numoflines < nz_max) then
            write(*,*) "**********************************************"
            write(*,"(A, A, A, I3, A, I3, A)") " The number of lines of the ", &
                trim(filename), " file has ", &
                numoflines, " lines, while ", nz_max, " were expected"
            write(*,*) "This program will be stopped"
            write(*,*) "**********************************************"
            stop
        end if

        ! goes back to the start of the file
        rewind(file_input_id)

        ! loops into the file until the rstatus has changed
        ! (i.e the enf of file is reached)

        errlec = 0
        do  while (errlec==0)
            Read(file_input_id,*,iostat=errlec) i,j,k , conc
            if (pos_nzmax=='SURF' .and. k <= nz_max) then 
                    val(i,j,nz_max-k+1) = conc
            elseif (pos_nzmax=='BOTT' .and. k<= nz_max) then
                    val(i,j,k) = conc
            else
                    val(i,j,k) = conc
            end if

        enddo
#if defined ECO3M_1DV
! some 1DV consider several horizontal grid meshes
        if ((nx_max - nx_min) /=0 .or. (ny_max - ny_min) /=0) THEN
           do k=1,nz_max
              val(:,:,k) = val(1,1,k) 
           enddo
        endif
#endif
        close(file_input_id)

    end subroutine init_prof
#endif

!-----------------------------------------------------------------------------------------
#ifdef CALC
    subroutine init_dz_array

    !> Initialises the dz global array in uncoupled mode, following
    !! the initialisation of the TKE model. It also converts the depth_z array
    !! from W to T levels.
    !! \author Nicolas Barrier
    !! \date 2016-04-26
!-----------------------------------------------------------------------------------------
 ! global variables and subroutine arguments
   use mod_eco3m
#if defined M3D_NCOUPL
   use mod_eco3m_files, only: filename_depth
#endif

 ! local variables
        integer :: i, j, k  ! Loop indexes
        real :: depth_temp(nx_min:nx_max, ny_min:ny_max, nz_max)

        ! copies the content of the depth_z variable in the depth_temp variable 
        depth_temp = depth_z

        if (pos_nzmax .eq. "BOTT") then  

            if(depth_temp(1, 1, 1) > depth_temp(1, 1, 2)) then 
                write(*,*) "**********************************************"
                write(*,*) "If nz_max = BOTT, the depth of the 1st level must be "
                write(*,*) "greater than that of the second layer "
                write(*,*) "Correct your the pos_nzmax value in the namelist"
                write(*,*) "**********************************************"
                stop
            endif

            ! Definition of T-grid width as the difference of 
            ! levels at W levels (read from file)
            do k = 1, nz_max-1
                do j = ny_min, ny_max
                    do i = nx_min, nx_max
                        dz(i, j, k) = abs(depth_temp(i, j, k+1) - depth_temp(i, j, k))
                        depth_z(i,j,k) = 0.5*(depth_temp(i,j,k)+depth_temp(i,j,k+1))
                    enddo
                enddo
            enddo

            ! dz at the last level
              dz(:,:,nz_max) = dz(i,j,nz_max-1)

            ! depth of the last level
              depth_z(:, :, nz_max) = depth_temp(:, :, nz_max)

            ! defines the bathymetry array
              bathymetry(:, :) = abs(depth_temp(:, :, nz_max))


        else if (pos_nzmax == "SURF") then

            if(depth_temp(1, 1, nz_max) > depth_temp(1, 1, nz_max-1)) then 
                write(*,*) "**********************************************"
                write(*,*) "If nz_max = SURF, the depth of the last level must be "
                write(*,*) "greater than that of the previous  layer "
                write(*,*) "Correct your the pos_nzmax value in the namelist"
                write(*,*) "**********************************************"
                stop
            endif

            ! Definition of T-grid width as the difference of 
            ! levels at W levels (read from file)
            do k = nz_max, 2, -1
                do j = ny_min, ny_max
                    do i = nx_min, nx_max
                        dz(i, j, k) = abs(depth_temp(i, j, k-1) - depth_temp(i, j, k))
                        depth_z(i,j,k) = 0.5*(depth_temp(i,j,k)+depth_temp(i,j,k-1))
                    enddo
                enddo
            enddo

            ! dz at the last level
              dz(i, j, 1) = dz(i, j, 2)

            ! depth of the last level
              depth_z(i, j, 1) = depth_temp(i, j, 2)
            
            ! defines the bathymetry array
              bathymetry(:, :) = abs(depth_temp(:, :, 1))

        end if


    end subroutine init_dz_array
#endif


!-----------------------------------------------------------------------------------------
    subroutine eco3m_files_close

!> Subroutine which closes all the files which are not read anymore after the model initialisation
!! (i.e. namelist, .inc files, kRGB file, config.ini, model.def)
!-----------------------------------------------------------------------------------------
    use mod_eco3m_files
    implicit none
    
    close(namelist_eco3m_id)
#ifdef key_nemo_eco3m
    close(namelist_ecophy_id)
#endif
    close(file_config_id)

#ifdef INI
! In INI mode, we close all the include files
   close(file_call_id)
   close(file_extincinc_id)
#ifdef SAVE_FLUX
   close (file_saveflx_id)
#endif
#endif 

! Closing the log file
  close(file_CR_id)

  end subroutine eco3m_files_close
!-----------------------------------------------------------------------------------------




