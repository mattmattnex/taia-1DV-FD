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

  module mod_eco3m_files

!> Module that contains the variables relative to the files that are used in Eco3M. 
!! \author Nicolas Barrier, M. Baklouti
!\date 23/08/2017
!----------------------------------------------------------------------------
    INTEGER,parameter::lfl=300 ! long filenames length
!-- Root directory for Eco3M:
    Character(100)           :: eco3m_root_dir
!-- Configuration files   
    Character(lfl)          :: filename_config !< "config.ini" file defining the
                                               ! biogeochemical fluxes between state variables 
    Integer                 :: file_config_id = 1002 
    Character(lfl)          :: filename_confmod !< "modele.def" file defining the
                                               ! process functions uses to
                                               ! calculate the fluxes between
                                               ! state variables 
    Integer                 :: file_confmod_id = 1054 
!-- Eco3M namelist
  integer :: namelist_eco3m_id = 1000 !< ID of the Eco3M Namelist file
#ifdef key_nemo_eco3m
  integer :: namelist_ecophy_id = 1001 !< ID of the Eco3M Namelist file
#endif

!--  files automatically generated during the run in the INI mode (key INI)
#ifdef INI
! call.inc file 
  integer :: file_call_id = 1004  !< ID of the call.inc file containing the
                                  !< equations of the biogeochemical model
  integer :: file_extincinc_id = 1005  !< ID of the calc_extinc.inc file
#ifdef SAVE_FLUX
  integer :: file_saveflx_id = 1006 !< ID of the call_save_flux.inc file
                                    !< allowing to save the biogeochemical fluxes
                                    !< at given positions
#endif
#ifdef INDIC
  integer :: file_saveind_id = 1007 !< ID of the call_save_indic.inc file
                                    !< allowing to save the biogeochemical indicators
#endif
#endif
!-- Position of the vertical nodes (in M3D_NCOUPL)
#ifdef M3D_NCOUPL
  character(lfl) :: filename_depth !< File containing the depth of the grid cells (Eco3M used alone)
#endif

!-- Eco3M initial condition files
  integer :: file_input_id = 1888  !< ID for all the Eco3M initial condition files files

!-- Files for light irradiance and extinction
! KRGB file
#ifdef RGB
  character(lfl) :: filename_krgb !< Name of the KRGB file
    !! This file is closed after a call 
    !! to the rgb_coef_init function.  
#endif

 integer :: file_irrad_id=1013 !< ID of the irradiance file
 character(lfl) :: filename_irrad !< Name of the irradiance file


  integer :: file_CR_id = 8001   !< ID of the Eco3M configuration summary file
  integer :: file_output_id = 8002  !< ID of the first output file
  integer :: file_budget_id = 8000  !< Id of the file containing C/N/P budgets 


end module mod_eco3m_files

