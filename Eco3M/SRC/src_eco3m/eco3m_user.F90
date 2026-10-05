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
      module mod_eco3m_user
 !>   contains:
 !>   - declarations of additional user variables
 !>   - allocate_globvar_arrays_user
 !>   \author Camille Mazoyer
 !>   \date mars 2020
      
      use mod_eco3m_vartypes  ! Module for variable types
      use eco3m_string ! Module for string manipulation
      use mod_eco3m_files  ! Module for file variables
      use mod_eco3m

      implicit none
      public


      TYPE(VAR_GLOB_USER), Allocatable :: mu_graz(:,:)  !< specific grazing rate
      TYPE(VAR_GLOB_USER), Allocatable :: mu_growth(:)  !< specific growth rate
      ! to complete with other user variables

#if defined INDIC
      !GPP = gross primary production
      !NPP = net primary production
      REAL(8),Allocatable ::GPPS(:,:,:),GPPL(:,:,:),NPPS(:,:,:),NPPL(:,:,:)
      !GSP = gross secondary production
      !NSP = net   secondary production
      REAL(8),Allocatable ::GSP(:,:,:),NSP(:,:,:)
      !GBP = gross bacterial production
      REAL(8),Allocatable ::GBP(:,:,:)
      !IND = intermediate variable for the calculation of indicators
      REAL(8),Allocatable ::IND(:,:,:)
#endif

      contains

!---------------------------------------------------------------------------

      subroutine allocate_globvar_arrays_user
!---------------------------------------------------------------------------

              
#ifdef CALC

 !> Subroutine dedicated to the allocation of additional global arrays 
 !> only for user variables used by the biogeochemical model.
 !> \author Camille Mazoyer
 !> \date mars 2020

 implicit none
 integer :: i,j  ! Dummy loop index
 integer :: nn
 integer :: istat  ! Allocation status


 ! to complete
#ifdef INDIC
! gross and net primary productions for small and large phytoplankton
Allocate(GPPS(nx_min:nx_max, ny_min:ny_max, nz_max))
Allocate(GPPL(nx_min:nx_max, ny_min:ny_max, nz_max))
Allocate(NPPS(nx_min:nx_max, ny_min:ny_max, nz_max))
Allocate(NPPL(nx_min:nx_max, ny_min:ny_max, nz_max))
Allocate(GSP(nx_min:nx_max, ny_min:ny_max, nz_max))
Allocate(NSP(nx_min:nx_max, ny_min:ny_max, nz_max))
Allocate(GBP(nx_min:nx_max, ny_min:ny_max, nz_max))
Allocate(IND(nx_min:nx_max, ny_min:ny_max, nz_max))
#endif
!----------------------------------------------------------
!-- specific growth  rate 
!----------------------------------------------------------
  istat=0
  if (nscp_cell /=0) allocate (mu_growth(nscp_cell),STAT=istat)
  if (istat/=0) write(*,*) 'pb with the allocation of  mu_growth'

  if (nscp_cell /= 0 .and. allocated(mu_growth)) then
    do i=1,nscp_cell
!          NULLIFY(mu_growth(i)%val)
          Allocate(mu_growth(i)%val(nx_min:nx_max, &
                  ny_min:ny_max,1:nz_max),stat = istat)
          if (istat /=0) write(*,*) " problem with the ", &
                     "allocation of mu_growth(", i,", ",j, ")"
              mu_growth(i)%val(:,:,:)=0.d0
              mu_growth(i)%idorg=iscp_cell(i)
    enddo
  endif
!----------------------------------------------------------
!-- specific grazing rate 
!----------------------------------------------------------
  istat=0
  nn = nscp_phy + nscp_zoo + nscp_bac
  if (nscp_zoo /=0) allocate (mu_graz(nscp_zoo,nn),STAT=istat)
  if (istat/=0) write(*,*) 'pb with the allocation of  mu_graz'

  if (nscp_zoo /= 0 .and. allocated(mu_graz)) then
    do i=1,nscp_zoo
       do j=1,nn
!          NULLIFY(mu_graz(i,j)%val)
          Allocate(mu_graz(i,j)%val(nx_min:nx_max, &
                  ny_min:ny_max,1:nz_max),stat = istat)
          if (istat /=0) write(*,*) " problem with the ", &
                     "allocation of mu_graz(", i,", ",j, ")"
              mu_graz(i,j)%val(:,:,:)=0.d0
              write(*,*)mu_graz(i,j)%val(1,1,1)
          enddo
      enddo
  endif


#endif

  end subroutine allocate_globvar_arrays_user
!---------------------------------------------------------------------------

  subroutine deallocate_globvar_arrays_user
!---------------------------------------------------------------------------
              
#ifdef CALC

 !> Subroutine dedicated to the deallocation of global var arrays 
 !> only for user variables used by the biogeochemical model.
 !> \author Camille Mazoyer
 !> \date mars 2021

 implicit none
 integer :: istat

 ! to complete
!----------------------------------------------------------
!-- specific grazing rate 
!----------------------------------------------------------
  istat=-1
  if (nscp_zoo /=0) deallocate (mu_graz,STAT=istat)
  if (istat/=0) write(*,*) 'pb with the deallocation of  mu_graz'

  !deAllocate(mu_graz(i,j)%val)
#endif

  end subroutine deallocate_globvar_arrays_user
!---------------------------------------------------------------------------
end module
