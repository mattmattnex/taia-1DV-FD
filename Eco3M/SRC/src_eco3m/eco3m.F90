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
!---------------------------------------------------------------------------
      Module mod_eco3m

! Module with the general Eco3m variables
!---------------------------------------------------------------------------

    use mod_eco3m_vartypes  ! Module for variable types

implicit none
!
    ! ==================== Chl related variables  ==================== 
!    real(8), allocatable :: Chl_tot(:,:,:) !< Total Chlorophyl concentration
    logical :: CHL_C_BOOL  !< False if the Chl:C ratio is constant, True if the ratio is variable
    logical :: CHL_C_FUNC  !< True if Chl:C ratio is computed by using an analytical function
    real(8) :: CHL_C0  !< value of the constant Chl:C ratio if needed
    integer, allocatable :: iphy_Chl(:)  !< Index of the phytoplankton Chl elements (if any)
    integer, allocatable :: iphy_C(:)  !< Index of the phytoplankton C elements (if any)
    TYPE(PARAM) :: CHL_C_PAR  !< Parameters of the Chl:C ratio function
    TYPE(VAR_GLOB_USER), Allocatable :: CHL_C(:) !< Chl:C ratios for sub-compartments containing chlorophyll

    ! =================== Light extinction related variables ==================== 
    REAL(8), Allocatable :: E_PAR(:,:) !< Surface Photosynthetic Available Radiation (2D)
    REAL(8), Allocatable :: E_PARZ(:,:,:)  !< Photosynthetic Available Radiation (3D)
    REAL(8), Allocatable :: kextinc(:,:,:)  !< Light extinction coefficient 
#ifdef RGB
    integer :: rgb_ncoeff  !< Number of coefficients (i.e lines in the file)
    integer :: rgb_nbands !< Number of columns in the file (filled in the namelist)
    real(8), allocatable :: rgb_coef_table(:, :)  !< Coefficients for light extinction
    REAL(8), Allocatable :: kextinc_rgb(:, :, :, :)  !< Light extinction coefficient for RGB
    REAL(8), Allocatable :: E_PARZ_RGB(:, :, :, :)  !< Photosynthetically Available Radiation for RGB 
#endif
    ! ======================== Flux related  variables calculation ====================
    TYPE(FLUX), Allocatable :: FLUX_VAL(:,:)  !< Array of biological fluxes between the variables
    TYPE(PARAM), Allocatable :: FLUX_PAR(:)  !< Array of subprocesses parameters

    ! ====================  Model related  variables ====================
    CHARACTER(20)::bioconf_name !< Name of the configuration of the biogeochemical model implemented in Eco3M
    TYPE(PROC), Allocatable :: PROC_MOD(:)  !< Array of model processes
    Integer :: nbvar        !< Number of state variables
    integer :: nb_proc      !< Number of process functions
    integer :: nbfunc_max   !< Maximum number of sub-processes within a flux
    integer :: nb_bio_flux  !< Number of direct mass fluxes in the biogeochemical model

    !====================  Physical variables within Eco3M ====================
    REAL(8), pointer :: TEMP_BIO(:,:,:)   !< Temperature used in the Eco3M model
    REAL(8), pointer :: SAL_BIO(:,:,:)    !< Salinity used in the Eco3M model
    REAL(8), pointer :: bathymetry(:,:)   !< Total water depth as function of longitude/latitudes
    REAL(8), pointer :: depth_z(:,:,:)    !< Depth of T-points 
    REAL(8), pointer :: dz(:,:,:)         !< Thickness of vertical grid cell (Tracer points)

    ! ================== Run related variables ====================
    ! Parameters defined in the config.ini or in the namelist_eco3m file
    Integer :: nx_min  !< Min spatial dimensions along X
    Integer :: nx_max  !< Max spatial dimensions along X
    Integer :: ny_min  !< Min spatial dimensions along Y
    Integer :: ny_max  !< Max spatial dimensions along Y
    Integer :: nz_max  !< Max number of vertical levels
    real(8) :: run_duration    !< Run duration (in days) **(only used if COUPL is not set)**
    real(8) :: dt_bio          !< Biological time-step (in seconds) 
    real(8) :: dt_save_bio     !< Biogeochemical saving time-step (in minutes) 
    Character(4) :: pos_nzmax  !< Position of the nz_max point (SURF or BOTT)
    Character(15) :: run_name  !< Name of the Eco3M run

    ! =============== Variable related variables ================================
    integer :: nbcomp    !< Number of compartments
    integer :: nbscomp   !< Number of sub-compartments
    integer :: nscp_phy  !< Number of phytoplankton PFTs
    integer :: nscp_zoo  !< Number of zooplankton PFTs 
    integer :: nscp_bac  !< Number of bacterial PFTs
    integer :: nphy_diaz !< Number of diazotroph phytoplankton PFTs
    integer :: nscp_mod  !< Number of dissolved organic sub-compartments 
    integer :: nscp_cell !< Number of cell elements

    Integer, Allocatable :: nb_elmt(:)   !< Number of elements per organism 
    Integer, Allocatable :: iscp_phy(:)  !< Index of phytoplankton PFTs
    Integer, Allocatable :: iscp_zoo(:)  !< Index of zooplankton PFTs
    Integer, Allocatable :: iscp_bac(:)  !< Index of bacterial PFTs
    Integer, Allocatable :: iphy_diaz(:) !< Index of diazotroph phytoplankton PFTs
    Integer, Allocatable :: iscp_mod(:)  !< Index of dissolved organic sub-compartments
    Integer, Allocatable :: iscp_cell(:) !< Index of cell elements

    ! declaration of biogeochemical state variables and fluxes
    TYPE(VAR_ETAT), Allocatable :: VAR(:) !< Biogeochemical state variables
    REAL(8), pointer :: TEND(:,:,:,:)     !< Biogeochemical trends

    TYPE(VAR_GLOB_USER), Allocatable :: mu_PPB(:)   !< specific gross primary production rate 
    TYPE(VAR_GLOB_USER), Allocatable :: PPB_NR(:)   !< gross primary production rate in nutrient replete cond.
    TYPE(VAR_GLOB_USER), Allocatable :: mu_PPB_NR(:)!< specific gross primary production rate in nut. replete cond.
    TYPE(VAR_GLOB_USER), Allocatable :: mu_resp(:)  !< specific respiration rate
    TYPE(VAR_GLOB_USER), Allocatable :: mu_nit(:)   !< specific nitrification rate
    TYPE(VAR_GLOB_USER), Allocatable :: PBB_NR(:)   !< bacterial production rate
    TYPE(VAR_GLOB_USER), Allocatable :: remin(:)    !< mineralisation rate


    real(8) :: tps  !< Simulation time (in seconds)
    
    REAL(8), parameter :: pi=3.14159265  !< Pi value
    logical :: dosave  !< Boolean that defines whether outputs are to be written

CONTAINS

End Module mod_eco3m
!--------------------------------------------------------------------------------
                     Module mod_eco3m_irrad

! Module specfic to light irradiance and extinction
!--------------------------------------------------------------------------------
    use mod_eco3m_vartypes  ! Module for variable types

implicit none
    real(8) :: dt_irrad  !< Time-step for the irradiance reading/calculation 
    real(8) :: irr2par   !< Conversion factor from irradiance to par (set to 1 if PAR is used instead of irradiance) 
    real(8) :: irr_param  !<  Parameter used to convert PAR(0+) into PAR(0-)
    Character(15) :: fichirrad  !< Definition of irradiance calculation.
    !! ( among IRR_FONCTION, IRR_filename, IRR_CODEPHYS or PAR_CODEPHYS)
    Character(40) :: fichirrad_long  !< Name of the irradiance file (if
    !! irradiance is read from a file)
    REAL(8) :: irrad_MAX  !< Maximum irradience (used if fichirrad=IRR_FONCTION)
    TYPE(PARAM) :: IRR_PAR  !< Parameters of the irradiance function
    REAL(8), pointer :: irrad(:,:)  !< 2D irradiance
    real(8) :: albedo  !< Albedo

CONTAINS

End Module mod_eco3m_irrad
!--------------------------------------------------------------------------------
                   Module mod_eco3m_outputs

! Module specific to Eco3M outputs
!--------------------------------------------------------------------------------
    
implicit none

    character(124) :: output_dir  !< Output directory (namelist)
    integer :: num_outphy  !< Number of Physical variables that will be saved (namelist)
    integer :: num_outflux !< Number of points at which biogeochemical fluxes will be saved (namelist)
    
    type outflux_type  !< Type associated with the saving of fluxes (namelist)
        integer :: indi  !< i-index 
        integer :: indj  !< j-index
        integer :: indk  !< k-index
    end type
    
    type outphy_type  !< Type associated with the saving of physical variables (namelist)
        character(124) :: outvarname    !< Name of the variable to output 
    end type

    type(outflux_type), allocatable :: outflux(:) !< Array of coordinates where to save fluxes
    type(outphy_type), allocatable :: outphy(:)   !< Array of physical variables to save

    integer :: filecount !< File counter (iterated after each writing in a file)

CONTAINS

end Module mod_eco3m_outputs
