!***************************************************************************
!***************************************************************************
!Copyright or � or Copr. CNRS/IRD/Universit� de la M�diterran�e
!contributor(s) : Melika BAKLOUTI & Vincent FAURE (10/10/2006)
!
!melika.baklouti@univ-amu.fr
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
!encouraged to load and testhe software''s suitability as regards their
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
!
   function f_regul_Chl_bak2026(ili,jcol,alpha_star,CHL_C_MAX,Q_MIN_PHY_C,&
	   Q_MAX_PHY_C,CN_min,CN_max)
!
!--------------------------------------------------------------------
! This function calculates the rate of chlorophyll production according to the formulation
! presented in Baklouti et al. (2021), which itself is adapted from Geider et al. (1998)
!
! variables, parameters and their units:
!----------------------------------------
! thetaN_max  = Chl:N maximum, gChl/molN
! alpha_star = initial slope of the P vs E curve, molC.m^2/(gChl.J)
! rho_chl    = efficacite de photosynthese, gChl/molN
! inut       = indice de la concentration du nutriment azot� absorb
! K          = constante de demi saturation de l'absorption de nutriment, molX m^{-3}
! ili        = contient l'indice de la concentration en chlorophylle de l'autotrophe
!
! f_regul_Chl_bak2026 = Vitesse de production de chlorophylle par mol de N absorbee, gChl/molN
!
!-- ordres de grandeur: thetaN_max= 0.3 * 14 gChl/molN (Geider, 1998)
!
! Refs:   Baklouti et al. (2006) Progr. in Oceanogr., 71:1-33
!         Geider et al. (1998) Limnol. Oceanog., 43(4):679-694
!
! Revisions:
! Melika Baklouti 27/06/2008
! Clement ALDEBERT 09/07/2012
! Remi Pages, Melika Baklouti 20/04/2018, 05/2025, 07/2025, 06/2026
!------------------------------------------------------------------
!-- Variables globales
    USE mod_eco3m_user

    Implicit None

    !-- arguments
    Real(8) :: alpha_star,CHL_C_MAX,Q_MIN_PHY_C,Q_MAX_PHY_C,CN_min,CN_max
    Real(8) :: f_regul_Chl_bak2026(nx_min:nx_max,ny_min:ny_max,nz_max)
    Integer :: ili,jcol

    !-- local variables 
    Character(L_CHAIN) :: chaine
    integer :: iorg
    integer :: i,j,k, ii,iiphy,ll
    real    :: nn
    integer :: iiC,iiChl,iiN, iicell,nsurf
    real(8) :: QC,fQC,CHL_C_MAX_VAR,CHL_C_MAX_LOC, thetaN_max,Q_MAX_PHY_N,NCmax,NCmin,QNC
    real(8) :: fCHL_C,fQNC,fregul,fCHL_N
    real(8) :: rhochl(nx_min:nx_max,ny_min:ny_max,nz_max), theta(nx_min:nx_max,ny_min:ny_max,nz_max)
    real(8) :: thetaN(nx_min:nx_max,ny_min:ny_max,nz_max)
    real :: nn_PHYS,nn_PHYL

    !-- Initialisation
    nn_PHYS =  0.06
    nn_PHYL =  4!0.06
    if (pos_nzmax.eq."BOTT") then
       nsurf = 1
    elseif  (pos_nzmax.eq."SURF") then
       nsurf = nz_max
    endif
   !-- conversion of C:X into X:C
     NCmin = 1.d0/CN_max
     NCmax = 1.d0/CN_min

     CHL_C_MAX_LOC = CHL_C_MAX
     thetaN_max = CHL_C_MAX_LOC * CN_min

    !-- The max Chl:C ratios is a function of available light (photoacclimation)
    
    if  (var(ili)%scomp == "phyl" .or. var(ili)%scomp == "PHYL" .or.  var(ili)%scomp == "Phyl") then
       CHL_C_MAX_VAR = CHL_C_MAX - 0.6      ! (C:Chl) must be > 10  gC/g Chl to synthetise new Chl
        nn=nn_PHYL
    elseif (var(ili)%scomp == "phys" .or. var(ili)%scomp == "PHYS" .or.  var(ili)%scomp == "Phys") then
       CHL_C_MAX_VAR = CHL_C_MAX - 0.6
        nn=nn_PHYS
    elseif (var(ili)%scomp == "phy" .or. var(ili)%scomp == "PHY" .or.  var(ili)%scomp == "Phy") then
       CHL_C_MAX_VAR = CHL_C_MAX - 0.6
       nn=nn_PHYS
    else
        STOP "problem in function f_regul_Chl_bak2026"
    endif

    !-- Recherche du numero d'organisme de l'autotrophe:
    iorg = var(ili)%idorg

    !-- de l'indice de sa conc en C :
    chaine = 'C'
    iiC = f_idorg2id (iorg,chaine)

    !-- de l'indice de sa conc en N :
    chaine = 'N'
    iiN = f_idorg2id (iorg,chaine)

    !-- et de celle de sa conc en Chl :
    chaine = 'Chl'
    iiChl = f_idorg2id (iorg,chaine)
    !-- et de celle de sa conc en Chl :
    chaine = 'cell'
    iicell = f_idorg2id (iorg,chaine)


    !--rapport CHL:C (g Chl/mol C):
    if ((CHL_C_BOOL .EQV. .TRUE.) .AND. (CHL_C_FUNC .EQV. .FALSE.) ) then
      theta = VAR(iiChl)%conc / (VAR(iiC)%conc + 1d-80)
      thetaN =  VAR(iiChl)%conc / (VAR(iiN)%conc + 1d-80)
    elseif ((CHL_C_BOOL .EQV. .TRUE.) .AND. (CHL_C_FUNC .EQV. .TRUE.) ) then
      theta = CHL_C(1)%val
      thetaN =  VAR(iiChl)%conc / (VAR(iiN)%conc+1d-80)
    elseif ((CHL_C_BOOL .EQV. .TRUE.) .AND. (CHL_C_FUNC .EQV. .TRUE.)) then
      theta = CHL_C(1)%val
      thetaN = CHL_C(1)%val * VAR(iiC)%conc / (VAR(iiN)%conc+1d-80)
    else
      theta = CHL_C0
      thetaN = CHL_C0 * VAR(iiC)%conc / (VAR(iiN)%conc+1d-80)
    endif

!--  Find the index of phytoplankton among all phytoplancton PFTs
    do ii=1,nscp_phy
      if (mu_PPB_NR(ii)%idorg == iorg) then
            iiphy = ii
            exit
      endif
    enddo
 f_regul_Chl_bak2026 = 0.0d0

!-- Chlorophyll production rate (in gChl/molN):

 if ((CHL_C_BOOL .EQV. .TRUE.) .AND. (CHL_C_FUNC .EQV. .FALSE.)) then
    do k=1,nz_max
     do j=ny_min,ny_max
       do i=nx_min,nx_max
            E_PARZ(i,j,k) = max(E_PARZ(i,j,k),0.0d0)

             CHL_C_MAX_LOC =  CHL_C_MAX_VAR * exp(-0.3 * E_PARZ(i,j,k)) + 0.6 
             thetaN_max = CHL_C_MAX_LOC * CN_min

      !-- Calcul of rho_chl coefficient as in  de Geider et al. (1998)
             rhochl(i,j,k) = mu_PPB_NR(iiphy)%val(i,j,k) / (alpha_star * E_PARZ(i,j,k) * theta(i,j,k) +1d-80) !-- egal a po finalement
             rhochl(i,j,k)  = min(1.0d0, rhochl(i,j,k) )
             rhochl(i,j,k)  = max(0.0d0, rhochl(i,j,k) )
    
        if (theta(i,j,k) <= CHL_C_MAX_LOC) then 
         
	!-- Regulation by C quota in phytoplankton
             QC = VAR(iiC)%conc(i,j,k) / (VAR(iicell)%conc(i,j,k) + 1.d-80)
             fQC= (Q_MAX_PHY_C-QC)/(Q_MAX_PHY_C-Q_MIN_PHY_C)
             fQC = min(1.0d0,fQC)
             fQC = max(0.0d0,fQC)

         !-- Regulation by Chl:N ratio in phytoplankton
            fCHL_N = ((thetaN_max - thetaN(i,j,k) )/(thetaN_max - 0.42))
            fCHL_N = min(1.0d0,fCHL_N)
            fCHL_N = max(0.0d0,fCHL_N)

         !-- Regulation by Chl:C ratio in phytoplankton
            fCHL_C = ((CHL_C_MAX_LOC - theta(i,j,k) )/(CHL_C_MAX_LOC - 0.12))
            fCHL_C = min(1.0d0,fCHL_C)
            fCHL_C = max(0.0d0,fCHL_C)

        !-- Regulation by the NC  ratio:
    
           QNC = var(iiN)%conc(i,j,k) / (var(iiC)%conc(i,j,k)+ 1.d-80) ! N:C ratio
           fQNC = (QNC  - NCmin) / (NCmax - NCmin)
           fQNC = max(0.d0,fQNC)
           fQNC = min(1.d0,fQNC)
                            
           fregul= fQC**0.06 * fQNC**nn *  fCHL_N**0.9 
        
	   f_regul_Chl_bak2026(i,j,k) =  thetaN_max *  rhochl(i,j,k) * fregul 
       endif
            enddo
        enddo
    enddo
 endif
 if ((CHL_C_BOOL .EQV. .FALSE.) .AND. (CHL_C_FUNC .EQV. .FALSE.)) then
      VAR(iiChl)%conc = CHL_C0 * VAR(iiC)%conc
  endif


end function f_regul_Chl_bak2026
