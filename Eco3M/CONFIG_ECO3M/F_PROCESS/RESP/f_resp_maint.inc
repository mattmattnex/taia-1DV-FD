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
!
  function f_resp_maint(ili,jcol,costm,QminC,QmaxC,iconc)
!
!--------------------------------------------------------------------------------
! Cette fonction calcule la respiration de maintenance en fonction de l'etat nutritionnel
!
! Variables et parametres et leurs unites:
!----------------------------------------- 
! costm     = cout de la maintenance
!  
! Ref: Thingstad,1987: Utilization of N;P and Organic C by heterotrophic bacteria. I. Outline of a chemostat 
!       theory with a consistent concept of 'maintenance metabolism'. MEPS Vol.35 99-109

!Authors: R. Mauriac, M. Baklouti
!--------------------------------------------------------------------


    Implicit None

    !-- fonction:
    Real(8) :: f_resp_maint(nx_min:nx_max,ny_min:ny_max,nz_max)

    !-- arguments
    Real(8)          :: costm,QminC,QmaxC
    Real(8), optional:: iconc
    Integer          :: ili,jcol

    !--variables locales:
    Integer            :: i,j,k
    Integer            :: iorg,iiCell,iiC,iiconc
    Real(8)            :: QC,reserve
    Character(L_CHAIN) :: chain1,chain2

    !-- Calcul du quota intracelulaire en carbone
    chain1 = 'cell'
    chain2 = 'C'
    if (present(iconc)) then 
      iiconc = nint(iconc)
      iorg = VAR(iiconc)%idorg
    else
      iorg = VAR(ili)%idorg
    endif
    iiCell = f_idorg2id(iorg,chain1)
    iiC    = f_idorg2id(iorg,chain2)


    do k = 1,nz_max
        do j = ny_min,ny_max
            do i = nx_min,nx_max
                QC = VAR(iiC)%conc(i,j,k)/(VAR(iiCell)%conc(i,j,k)+1d-80)
                reserve  = QC - QminC
      
               if (reserve > 0.0d0 .and. VAR(iiC)%conc(i,j,k) > costm*(1.d0+reserve/(QC+1.d-80))*VAR(iiC)%conc(i,j,k) * dt_bio) then
                    f_resp_maint(i,j,k) = costm*(1.d0+reserve/(QC + 1.d-80))*VAR(iiC)%conc(i,j,k) 
               else
                  if (VAR(iiC)%conc(i,j,k) > costm*VAR(iiC)%conc(i,j,k)* dt_bio) then 
                    f_resp_maint(i,j,k) = costm*VAR(iiC)%conc(i,j,k) 
                  else
                    f_resp_maint(i,j,k) = 0.d0
                  endif
               endif

            enddo
        enddo
    enddo


end function f_resp_maint
