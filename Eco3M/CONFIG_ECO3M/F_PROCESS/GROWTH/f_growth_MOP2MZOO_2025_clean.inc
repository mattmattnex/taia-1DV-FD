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
!---------------------------------------------------------------------------
!
     function f_growth_MOP2MZOO_2025_clean(ili,jcol,indcell,indC,QC_MIN,QC_MAX,indDETC)
!
!---------------------------------------------------------------------------
! Function that calculates the DET consumption rate per MZOO
!
! Melika Baklouti 16/10/2024; last revision, June 2026
!-------------------------------------------------------------------------
    USE mod_eco3m_user
    Implicit None

    !-- function
    Real(8):: f_growth_MOP2MZOO_2025_clean(nx_min:nx_max,ny_min:ny_max,nz_max)

    !-- arguments
    Real(8) :: indC,indcell,indDETC,QC_MIN,QC_MAX
    Integer :: ili,jcol

    !--local variables 
    Real(8) :: QC,fQ
    Integer :: iorgX,iorgCell
    Integer :: iindC,iiCell,iindDETC
    Integer :: i,j,k,l,ll
    Character(L_VAR)    :: chain
    Real    :: nn

    nn = 2!0.06 

    !-- Index of population abundance
    iiCell = nint(indcell)
    iorgCell =var(iiCell)%idorg


    !-- MZOO carbon concentration index
    iindC = nint(indC)
    iorgX=var(iindC)%idorg

    !-- DET concentration index
    iindDETC = nint(indDETC)

    !-- Verification  
    if (iorgCell /= iorgX ) then
        write(*,*) 'pb avec f_growth_MOP2MZOO_2025_clean'
        STOP
    endif

    !-- Search for the mu_growth index that corresponds to iiCell
    do l = 1,nscp_cell
        if (mu_growth(l)%idorg == iiCell) then !! A VOIR iorg plutôt que iiCell ds le test ?
           ll = l
        endif
    enddo

    do k=1,nz_max
        do j=ny_min,ny_max
            do i=nx_min,nx_max

                QC = var(iindC)%conc(i,j,k) / (var(iiCell)%conc(i,j,k)+1.d-80)
                fQ = (QC_MAX - QC)/(QC_MAX - QC_MIN)
                fQ = min(1.d0,fQ)
                fQ = max(0.d0,fQ)

            !-- DET consumption is proportional to growth rate and carbon quota 
                f_growth_MOP2MZOO_2025_clean(i,j,k) =  mu_growth(ll)%val(i,j,k)* var(iindC)%conc(i,j,k) * fQ**nn
                !f_growth_MOP2MZOO_2025_clean(i,j,k) =  mu_growth(ll)%val(i,j,k)* var(iiCell)%conc(i,j,k) *QC_MIN* fQ**nn
            !-- ensuring that the concentration in the grid cell is sufficient for a time step:
                if (VAR(iindDETC)%conc(i,j,k) <  1.1 * f_growth_MOP2MZOO_2025_clean(i,j,k) * dt_bio ) then 
                        f_growth_MOP2MZOO_2025_clean(i,j,k) = 0.d0
                endif
            enddo
        enddo
    enddo
end function f_growth_MOP2MZOO_2025_clean
!--------------------------------------------------------------------
