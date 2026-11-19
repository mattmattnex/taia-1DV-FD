!***************************************************************************
!***************************************************************************
!Copyright or © or Copr. CNRS/IRD/Université de la Méditerranée
!contributor(s) : Melika BAKLOUTI & Vincent FAURE (10/10/2006)
!
!m.baklouti@univmed.fr; vincent.faure@univmed.fr
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
function f_graz_hol2_mpreys_X_2025(ili,jcol,ipred,iprey)
    !
    !--------------------------------------------------------------------
    ! Cette fonction est à utiliser parallèlement à f_graz_hol2_mpreys
    ! 
    ! Elle permet d'exprimer la consommation d'une proie exprimee en un élément X par un prédateur
    ! dont la concentration est aussi exprimée en X, en presence de PLUSIEURS PROIES potentielles
    !
    ! Ce flux de proies consommé exprimé en cell a déjà été calculé par la fonction f_graz_hol2_mpreys
    ! (ce qui a permis de remplir le tableau mu_graz)
    ! 
    ! La fonction f_graz_hol2_mpreys_X_2025 se borne donc à convertir ce flux en une unité de biomasse 
    !
    ! Reponse fonctionnelle: 
    !     * de type Holling II
    !     * tient compte des interferences entre proies
    !     * preferences fixes
    !
    ! Attention:
    ! ----------
    !
    ! * Le calcul se base sur les concentrations en C du prédateur et des proies,
    !  et on obtient une vitesse BRUTE de grazing  en [mol X] l-1 s^{-1}
    !
    ! Variables et parametres et leurs unites:
    !----------------------------------------- 
    ! ipred    : indice de la variable d'état représentant la concentration X du prédateur
    ! iprey    : indice de la variable d'etat représentant la concentration X de la proie
    ! mu_graz  : ce tableau est une variable globale. Il contient les vitesses specifiques (en s^-1) de grazing des prédateurs.
    !            Ce tableau est rempli lors de l'appel de la fonction  f_graz_hol2_mpreys.
    !            l'élément (i,j) de ce tableau représente la vitesse de grazing du prédateur VAR(i) sur la proie VAR(j)
    !
    ! f_graz_hol2_mpreys_X_2025  : vitesse de consommation (en mol X m^-3 s-1) de la proie iprey par le prédateur ipred
    !
    !
    ! Ref: Holling, C. (1959)  Canadian Entomologist, 91:385-398
    !
    ! Melika Baklouti 6/07/2007
    !--------------------------------------------------------------------

    Implicit None

    !-- fonction:
    Real(8)             :: f_graz_hol2_mpreys_X_2025(nx_min:nx_max,ny_min:ny_max,nz_max)

    !-- arguments
    Integer             :: ili,jcol
    Real(8)             :: ipred,iprey
    !--variables locales:
    Integer             :: iipred,iiprey,iCell_pred,iCell_prey,iscp_prey,iscp_pred,nn
    Integer             :: iorg_prey,iorg_pred,iorg_ili,iorg_jcol
    ! Character(L_VAR)    :: elmt
    Character(L_VAR)    :: chain

    iipred=nint(ipred)
    iiprey=nint(iprey)

    !--recherche des numeros d'organisme de la proie et du predateur :
    iorg_prey= var(iiprey)%idorg
    iorg_pred= var(iipred)%idorg
    iorg_ili=var(ili)%idorg

    ! barrier.n
    ! if (jcol /=0)  iorg_jcol=var(jcol)%idorg  ! if vector is not a self
    if (jcol /= ili) iorg_jcol=var(jcol)%idorg

    !-- concentration en cellule de la proie et du predateur:
    chain='cell'
    iCell_prey = f_idorg2id(iorg_prey,chain)
    iCell_pred = f_idorg2id(iorg_pred,chain)


    !-- Rappel: mu_graz = vitesse specifique de grazing du predateur (Gamma x G) en s^-1
    !          calculé pour chaque couple pred/proie dans la fonction f_graz_hol2_mpreys

    !-- indice de l organisme du predateur (zoo):
    iscp_pred = 0
    iscp_prey = 0

    do nn = 1,nscp_zoo
        if (iscp_zoo(nn)==iorg_pred) then
            iscp_pred = nn
        elseif (iscp_zoo(nn)==iorg_prey) then
            iscp_prey = nn
        endif
    enddo
    !-- indice de l organisme de la proie dans le vecteur iscp_phy ou iscp_bact:
    if (iscp_prey ==0) then
        if (allocated (iscp_phy)) then
            do nn = 1,nscp_phy
                if (iscp_phy(nn)==iorg_prey) then
                    iscp_prey = nn + nscp_zoo
                    exit
                endif
            enddo
        endif
    endif

    if (iscp_prey ==0) then
        if (allocated (iscp_bac)) then
            do nn = 1,nscp_bac
                if (iscp_bac(nn)==iorg_prey) then
                    iscp_prey = nn +nscp_zoo + nscp_phy
                    exit
                endif
            enddo
        endif
    endif

    !-- vitesse de grazing en mol X m^-3 s-1

        f_graz_hol2_mpreys_X_2025 = mu_graz(iscp_pred,iscp_prey)%val * var(iCell_pred)%conc &
            * var(iiprey)%conc / (var(iCell_prey)%conc + 1d-80)           !mu x Z_C x (X/C)_prey


End function f_graz_hol2_mpreys_X_2025
