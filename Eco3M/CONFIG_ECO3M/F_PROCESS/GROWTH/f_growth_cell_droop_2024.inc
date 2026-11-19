!***************************************************************************
!***************************************************************************
!Copyright or © or Copr. CNRS/IRD/Université de la Méditerranée
!contributor(s) : Melika BAKLOUTI & Vincent FAURE (10/10/2006)
!
! melika.baklouti@univ-amu.fr; vincent.faure@univmed.fr
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
!
    function f_growth_cell_droop_2024(ili,jcol,mumax,indcell,indX1,Qmin1,Qmax1, &
        &     indX2,Qmin2,Qmax2,indX3,Qmin3,Qmax3)
!
!---------------------------------------------------------------------------
! This function calculates the net specific cell division rate taking into account 
! the intracellular nutritional state of organisms through Droop quota functions
!
!
! Variables, parameters and their units :
!-----------------------------------------
! indX: index of the intracellular  conc. in element  X (X among C, N, P, ...)
! Qmin: minimum value of the cell quota in element X (mol X (cell)^-1)
!
!
! Refs: Droop 
! Melika Baklouti: last revision on 09/07/2024
!-------------------------------------------------------------------------
    USE mod_eco3m_user
    Implicit None

    !-- fonction:
    Real(8):: f_growth_cell_droop_2024(nx_min:nx_max,ny_min:ny_max,nz_max)

    !-- arguments
    Real(8) :: Qmin1,Qmax1,indX1,mumax,indcell
    Real(8),optional :: indX2,Qmin2,Qmax2,Qmin3,Qmax3,indX3
    Integer :: ili,jcol

    !-- local variables:
    Integer :: iorg1,iorg2,iorg3,iorgCell
    Integer :: iindX1,iindX2,iindx3,iiCell
    Integer :: nbnut
    Integer :: i,j,k,l,lab
    Real(8) :: Q1,fQ1,Q2,fQ2
    Real(8) :: Q3,fQ3
    Real(8) :: fQ

    !-- Cell concentration index:
    iiCell = nint(indcell)
    iorgCell =var(iiCell)%idorg

    !-- Concentration index of the first element (required)
    iindX1 = nint(indX1)
    iorg1=var(iindX1)%idorg
    nbnut = 1

    if (present(indx2))then
        iindx2 = nint(indX2)
        iorg2 = var(iindX2)%idorg
        nbnut = 2
        if (present(indx3))then
            iindx3 = nint(indX3)
            iorg3 = var(iindX3)%idorg
            nbnut = 3
        endif
    endif

    !-- verification test
    Select case (nbnut)
    case(1)
        if (iorgCell /= iorg1 ) then
            write(*,*) 'pb avec f_growth_cell_droop_2024'
            STOP
        endif
    case(2)
        if (iorgCell /= iorg1 .or. iorgCell/= iorg2) then
            write(*,*) 'pb avec f_growth_cell_droop_2024'
            STOP
        endif
    case(3)
        if (iorgCell /= iorg1 .or. iorgCell /= iorg2  .or. iorgCell /= iorg3 ) then
            write(*,*) 'pb avec f_growth_cell_droop_2024'
            STOP
        endif
    end select
  !-- Save the specific cell division rate of each organism : 
         do l = 1,nscp_cell
            if (mu_growth(l)%idorg == iiCell) then
                 lab = l 
                 exit    
            endif
        enddo
    
!-- Calcul of the limitation function (when nbnut == 1 ) :
if (nbnut == 1) then
    do k=1,nz_max
        do j=ny_min,ny_max
            do i=nx_min,nx_max

                fQ1 = 0.d0

                Q1 = var(iindX1)%conc(i,j,k) / (var(iiCell)%conc(i,j,k)+1.d-80)

                ! Calculation of the quota functions for population growth (cell division):
                if (Q1 >= 0.d0)   fQ1 = 1.d0 - Qmin1/(Q1+1.d-80) 

                fQ1 = max(0.d0,fQ1)

                !-- calculation of the minimum quota :
                fQ =  fQ1

                !-- Calculation of the effective cell division rate  (cell.l-1.s-1)
                f_growth_cell_droop_2024(i,j,k) =  mumax * fQ * var(iiCell)%conc(i,j,k)

                !-- Save the cell specific  division rate  (s-1)
                mu_growth(lab)%val(i,j,k) = mumax * fQ
            enddo
        enddo
    enddo
!-- Calcul of the limitation function  :
elseif (nbnut == 2) then
    do k=1,nz_max
        do j=ny_min,ny_max
            do i=nx_min,nx_max

                fQ1 = 0.d0
                fQ2 = 0.d0

                Q1 = var(iindX1)%conc(i,j,k) / (var(iiCell)%conc(i,j,k)+1.d-80)
                Q2 = var(iindX2)%conc(i,j,k) / (var(iiCell)%conc(i,j,k)+1.d-80)

                ! Calculation of the quota functions for population growth (cell division):
                if (Q1 >= 0.d0)   fQ1 = 1.d0 - Qmin1/(Q1+1.d-80) 
                if (Q2 >= 0.d0)   fQ2 = 1.d0 - Qmin2/(Q2+1.d-80) 

                fQ1 = max(0.d0,fQ1)
                fQ2 = max(0.d0,fQ2)

                !-- calculation of the minimum quota :
                fQ =  min(fQ1, fQ2)

                !-- Calculation of the effective cell division rate  (cell.l-1.s-1)
                f_growth_cell_droop_2024(i,j,k) =  mumax * fQ * var(iiCell)%conc(i,j,k)

                !-- Save the cell specific  division rate  (s-1)
                mu_growth(lab)%val(i,j,k) = mumax * fQ
            enddo
        enddo
    enddo
elseif (nbnut == 3) then
    do k=1,nz_max
        do j=ny_min,ny_max
            do i=nx_min,nx_max

                fQ1 = 0.d0
                fQ2 = 0.d0
                fQ3 = 0.d0

                Q1 = var(iindX1)%conc(i,j,k) / (var(iiCell)%conc(i,j,k)+1.d-80)
                Q2 = var(iindX2)%conc(i,j,k) / (var(iiCell)%conc(i,j,k)+1.d-80)
                Q3 = var(iindX3)%conc(i,j,k) / (var(iiCell)%conc(i,j,k)+1.d-80)

                ! Calculation of the quota functions for population growth (cell division):
                if (Q1 >= 0.d0)   fQ1 = 1.d0 - Qmin1/(Q1+1.d-80) 
                if (Q2 >= 0.d0)   fQ2 = 1.d0 - Qmin2/(Q2+1.d-80) 
                if (Q3 >= 0.d0)   fQ3 = 1.d0 - Qmin3/(Q3+1.d-80) 

                fQ1 = max(0.d0,fQ1)
                fQ2 = max(0.d0,fQ2)
                fQ3 = max(0.d0,fQ3)

                !-- calculation of the minimum quota :
                fQ =  min(fQ1, fQ2, fQ3)

                !-- Calculation of the effective  cell division rate  (cell.l-1.s-1)
                f_growth_cell_droop_2024(i,j,k) =  mumax * fQ * var(iiCell)%conc(i,j,k)

                !-- Save the cell specific division rate  (s-1)
                mu_growth(lab)%val(i,j,k) = mumax * fQ
            enddo
        enddo
    enddo
endif
end function f_growth_cell_droop_2024
!--------------------------------------------------------------------
