!***************************************************************************
!***************************************************************************
!Copyright or ï¿œ or Copr. CNRS/IRD/Universite de la Mediterranee
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
!--------------------------------------------------------------------
!
    function f_graz_KTW_revis(ili,jcol,ipred,iprey,nbprey,Beta,&
    iCell1,Imax1,F1,pref1,iCell2,Imax2,F2,pref2,icell3,Imax3,F3,pref3)
!
!--------------------------------------------------------------------
! This function calculates the predation rate of a predator
! that can feed on SEVERAL potential prey species; 
! here, the maximum number of preys is 3, but the code can easily 
! be adapted to increase this number
!
! Functional response: 

!     * class 3 multiple functional responses (cf. Gentleman et al. 2003 DSR 50: 2847-2875)
!     * takes into account interference between prey
!     * weighted preferences, i.e. variable preferences which are dependent on prey concentration (allows prey coexistence)
!     * maximum feeding     
!     * active prey switching
!
! Warning:
! ----------
!
! * this function applies to a FLUX involving
!   a predator whose concentration is given by VAR(ipred)%conc and
!   a prey whose concentration is given by VAR(iprey)%conc
!
! * The calculation is based on the C concentrations of the predator and prey,
!   and we obtain a GROSS Graz grazing rate in [mol C] m-3 s^{-1}
! which does not take into account losses due to egestion, and sloppy feeding
! 
!
! Variables, parameters and their units:
!---------------------------------------
!
! pref_i: fixed preference of the predator for prey i, (without unit [0,1])
! nbprey: number of prey likely to be grazed by the predator (between 1 and 4)
! Im_i: maximum rate of ingestion (in s-1)
! beta: constant that determines the form of the feeding probability Q = 2 in the reference paper
! beta = 1: type II response; beta = 2: type III response
! f_graz_KTW_revis: consumption rate (in preys/pred/s) of the prey iprey (in ind m-3)
! by the predator ipred (in ind m-3)
!
!
! Ref: Vallina et al. Progress in Oceanography 120 (2014),93-109
!
! Frederic Diaz, Melika Baklouti  02/05/2014, 31/07/2025
!--------------------------------------------------------------------

    Implicit None

    !-- fonction:
    Real(8)             :: f_graz_KTW_revis(nx_min:nx_max,ny_min:ny_max,nz_max)

    !-- arguments
    Real(8)             :: ipred,iprey,iprey_X,nbprey,iCell1,pref1,Beta,F1,Imax1
    Integer             :: ili,jcol
    Real(8), optional   :: iCell2,iCell3
    Real(8), optional   :: pref2,pref3
    Real(8), optional   :: F2,F3,Imax2,Imax3 
    !--local variables :
    Integer              :: int_nbprey,iiprey,iipred,pos_iprey,iscp_pred,iscp_prey
    Integer              :: iCell_pred,iCell_prey,iorg_prey,iorg_pred,iorg_prov
    Integer              :: i,j,k,ii,mm
    Integer              :: indCell(nint(nbprey))  ! index of prey abundance concentration
    Integer              :: indC(nint(nbprey))     ! index of prey carbon concentration
    Real(8)              :: pref(nint(nbprey)),vfilt(nint(nbprey)),Imax(nint(nbprey)),AM(nint(nbprey))
    Real(8)              :: ress,ress_v,Kad
    Real(8)              :: pref_v(nint(nbprey))
    Real(8)              :: sum_selec,Q,delta
    Character(L_CHAIN)   :: chaine

    !-- conversion in integer :
    int_nbprey = nint(nbprey)
    iipred = nint(ipred)
    iiprey = nint(iprey)

    !-- initialisation:
    indCell = 0
    indC = 0
    pref = 0.d0
    vfilt= 0.0
    Imax = 0.0
    AM   = 0.0
    f_graz_KTW_revis=0.d0

   
    !-- First prey (required):
    indCell(1) = nint(iCell1)
    iorg_prov  = VAR(indCell(1))%idorg
    chaine     = 'C'
    indC(1)    = f_idorg2id(iorg_prov,chaine)
    pref(1)    = pref1
    vfilt(1)   = F1
    Imax(1)    = Imax1
    AM(1)      = Imax(1)/vfilt(1)

    Select case (int_nbprey)
    case (2)
        indCell(2) = nint(iCell2)
        iorg_prov  = VAR(indCell(2))%idorg
        chaine     = 'C'
        indC(2)    = f_idorg2id(iorg_prov,chaine)
        pref(2)    = pref2
        vfilt(2)   = F2
        Imax(2)    = Imax2
	AM(2)      = Imax(2)/vfilt(2)

    case (3)
        indCell(2)  = nint(iCell2)
        iorg_prov   = VAR(indCell(2))%idorg
        chaine      = 'C'
        indC(2)     = f_idorg2id(iorg_prov,chaine)
        pref(2)     = pref2 
        vfilt(2)    = F2
        Imax(2)     = Imax2
	AM(2)       = Imax(2)/vfilt(2)

        indCell(3)  = nint(iCell3)
        iorg_prov   = VAR(indCell(3))%idorg
        chaine      = 'C'
        indC(3)     = f_idorg2id(iorg_prov,chaine)
        pref(3)     = pref3
        vfilt(3)    = F3
        Imax(3)     = Imax3
	AM(3)       = Imax(3)/vfilt(3)

    end select
    

    !-- index of the organisms associated with preys and predators:
     iorg_prey = VAR(iiprey)%idorg
     iorg_pred = VAR(iipred)%idorg

    !-- index of cell concentrations of preys and predators:
     chaine = 'cell'
     iCell_prey = f_idorg2id(iorg_prey,chaine)
     iCell_pred = f_idorg2id(iorg_pred,chaine)


    !-- initialisations:
     iscp_pred = 0
     iscp_prey = 0

    !-- search for the organism index of the predator (zoo) and that of its prey among zooplankton:
    do mm = 1,nscp_zoo
        if (iscp_zoo(mm)==VAR(iipred)%idorg) then
            iscp_pred = mm
        elseif (iscp_zoo(mm)==VAR(iiprey)%idorg) then
            iscp_prey = mm
        endif
    enddo


    !-- if not found in the zoo, search for prey organism index in vector iscp_phy then iscp_bact:

    if( iscp_prey == 0) then
        if (allocated (iscp_phy)) then
            do mm = 1,nscp_phy
                if (iscp_phy(mm)==VAR(iiprey)%idorg) then
                    iscp_prey = mm + nscp_zoo
                    exit
                endif
            enddo
        endif
    endif
    if (iscp_prey ==0) then
        if (allocated (iscp_bac)) then
            do mm = 1,nscp_bac
                if (iscp_bac(mm)==VAR(iiprey)%idorg) then
                    iscp_prey = mm +nscp_zoo + nscp_phy
                    exit
                endif
            enddo
        endif
    endif


  ! search for position pos_iprey (between 1 and nbprey) relative to
  ! the cell conc of prey iprey in vector indCell
    pos_iprey=0

    do i=1,int_nbprey
        if (indCell(i)==iCell_prey) then
            pos_iprey = i
            exit
        endif
    enddo
    
    if (pos_iprey ==0) then
        write(*,*) 'probleme d''indices de proies dans l appel de la fonction f_graz_KTW_revis',ili,jcol
        read(*,*)
    endif
    

    do k=1,nz_max
      do j=ny_min,ny_max
        do i=nx_min,nx_max
          !-- Initialisation of the vector for variable preferences (selectivites):
          pref_v=0.d0
          !-- Initialisation of the total food using variable and fixed preferences: 
          ress = 0.d0
          ress_v = 0.d0
          do ii = 1, int_nbprey
              ! total food in terms of normalized abundances using fixed  preferences  :
                  ress    = ress +  pref(ii) * (var(indCell(ii))%conc(i,j,k)/AM(ii))
          enddo

          do ii=1,int_nbprey
              ! variable preferences (dimensionless), calculated with carbon biomasses
                  pref_v(ii) = pref(ii) * (var(indCell(ii))%conc(i,j,k)/AM(ii)) / (ress+1.d-80)

              ! total food in terms of normalized abundances  using variable  preferences  :
                  ress_v  = ress_v +  pref_v(ii) * (var(indCell(ii))%conc(i,j,k)/AM(ii))
          enddo
    !--Dimensionless saturation constant (ksat*F_phi/F_rho in Vallina et al. 2014, with ksat = 1 here):
           Kad =  ress_v/(ress + 1.d-80)

   !-- predators’ probability of feeding as a saturating function of the available prey biomass
	  Q = (ress_v)**beta / ( Kad**beta + ress_v**beta + 1.d-80 )
            
          delta = pref_v(pos_iprey) * (var(iCell_prey)%conc(i,j,k)/AM(pos_iprey)) /(ress_v+1d-80)

   !-- specific grazing rate of prey  iscp_prey by predator  iscp_pred (nb prey /nb pred /s): 
           mu_graz(iscp_pred,iscp_prey)%val(i,j,k) =  Imax(pos_iprey) * delta * Q 

   !-- Gross grazing rate  in  (nb preys)l-1 s-1
        if (var(iCell_prey)%conc(i,j,k)  > 1.1 * mu_graz(iscp_pred,iscp_prey)%val(i,j,k) * var(iCell_pred)%conc(i,j,k) * dt_bio) then
              f_graz_KTW_revis(i,j,k) = mu_graz(iscp_pred,iscp_prey)%val(i,j,k) * var(iCell_pred)%conc(i,j,k)                  ! mu x Z_Cell
        else
          f_graz_KTW_revis(i,j,k) = 0.d0
          mu_graz(iscp_pred,iscp_prey)%val(i,j,k) = 0.d0
        endif
       enddo
     enddo
    enddo

    ! write(*,*)ili,jcol, mu_graz(iscp_pred,iscp_prey)%val(1,1,1) ,   f_graz_KTW_revis(1,1,1)
    ! write(*,*) 'end of function f_graz_KTW_revis'

End function f_graz_KTW_revis
!--------------------------------------------------------------------
