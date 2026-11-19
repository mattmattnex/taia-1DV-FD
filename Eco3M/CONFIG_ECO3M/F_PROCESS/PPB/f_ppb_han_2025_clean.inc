!***************************************************************************
!***************************************************************************
!Copyright or © or Copr. CNRS/IRD/Université de la Méditerranée
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
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!
function f_ppb_han_2025_clean(ili,jcol,slopePE,tau,&
        sig_ps2,k_d,k_r,iphyC)
!
! This function calculates the GROSS primary production rate of the state variable
! whose concentration is given by the index variable ili
!
! Parameter variables and their units:
!--------------------------------------
! iphyC     : concentration of phytoplankton in carbon (mol C/l). If absent, ili used by default
! a_chl     : specific absorption coefficient of Chl, m^2 / (g Chl)
! phi_max_C : maximum quantum yield, (mol C)/ J
! slopePE   : initial slope  of P vs E curve = a_chl x phi_max_C
! tau       : electron turnover time, s
! sig_ps2   : effective cross-section of PSII, m^2 / Joules
! kd_Han    : PSII damage constant, adim in Han model
! kr        : PSII repair constant, s^-1
! f_ppb_han_slope_2026: in molC/l/s
!
! this function also fills in the muphy_NR table

! muphy_NR: specific primary production rate in the absence of nutrient limitations, in s^-1
!
! Ref:  Han, B. P. (2002) J. Theor. Biol. 214:519-527
!
! Melika Baklouti 6/07/2007, last revision: 22/07/2025
!-------------------------------------------------------------------------

    Implicit None

    !-- arguments
    Real(8) :: slopePE,tau,sig_ps2,k_d,k_r
    Real(8),optional:: iphyC
    Real(8) :: f_ppb_han_2025_clean(nx_min:nx_max,ny_min:ny_max,nz_max)
    Integer :: ili,jcol, ii

    !-- local variables 
    Integer:: iorg,iichl,iiC,iiphyC
    Real(8)  :: po(nx_min:nx_max,ny_min:ny_max,nz_max),Chl(nx_min:nx_max,ny_min:ny_max,nz_max)
    Character(L_VAR)::elmt

    !-- proportion of open PSII 
    po = 1.e0/(1.d0 + sig_ps2 &
        * E_PARZ * tau + k_d/k_r * (sig_ps2 * E_PARZ)**2 *tau)

    !-- organism:
    if (present(iphyC)) then
      iiphyC = nint(iphyC)
      iorg = var(iiphyC)%idorg
    else
      iorg = var(ili)%idorg
      iiphyC = ili
    endif

    !-- case of a constant  CHL:C ratio
    if (.not. CHL_C_BOOL) then
        Chl = CHL_C0 * var(iiphyC)%conc
    elseif ((CHL_C_BOOL) .and. (.not. CHL_C_FUNC)) then

        !-- Phytoplankton Chlorophyll Concentration Index 
        elmt ='Chl'
        iichl = f_idorg2id (iorg,elmt)
        Chl = var(iichl)%conc
    else !-- Case in which CHL:C is calculated using an analytical function
        elmt ='C'
        iiC = f_idorg2id (iorg,elmt)
        Chl = CHL_C(1)%val * var(iiC)%conc
    endif

    !-- Primary production rate (in molC/m3/s):
    f_ppb_han_2025_clean =  slopePE * po * E_PARZ * Chl


    !-- Save the specific growth rate under
    !-- nutrient replete conditions, in s^-1.
    ! Important: check that the dynamic allocation of mu_PPB_NR
    ! in the ALLOC_VAR_Eco3M subroutine of the coupler.F90 file
    ! is operational
    !
    do ii = 1, nscp_phy
        if (mu_PPB_NR(ii)%idorg == iorg) then
            mu_PPB_NR(ii)%val = &
                slopePE*po*E_PARZ*Chl/(var(iiphyC)%conc+1.0d-80)
            exit
        endif
    enddo

end function f_ppb_han_2025_clean
