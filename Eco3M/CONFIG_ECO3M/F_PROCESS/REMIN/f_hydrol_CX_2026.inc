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
   function f_hydrol_CX_2026(ili,jcol,indMOP,TT,indBACC,QminC,QmaxC,QminN,QmaxN,QminP,QmaxP)
!
!---------------------------------------------------------------------------
! This function calcules the rate of particulate detrital organic matter (DET) by bacteria
! based on a fixed turnover time
!
! Variables and parameters and their units:
!-----------------------------------------
! indMOP  : index of DET concentration
! indBACC: index of bacteria C biomass 
! mumax_B: max specific growth rate of bacteria 
!
!
! R. Mauriac  07/09/2009
! Last revision:  M. Baklouti 11/01/2026
!-------------------------------------------------------------------------

    Implicit None

    !-- function:
    Real(8):: f_hydrol_CX_2026(nx_min:nx_max,ny_min:ny_max,nz_max)

    !-- arguments
    Real(8) :: TT,indMOP,indBACC,QminC,QmaxC
    Real(8),optional :: QminN,QmaxN,QminP,QmaxP
    Integer :: ili,jcol

    !--local variables:
    integer :: i,j,k,l,iiMOP,igB
    integer :: iiC,iiCell,iielmt,iiorg
    real(8) :: nnC,nnX
    real(8) :: QC,fQC,QX,fQX,fQ,Qmin,Qmax
    character(L_VAR) :: chain,elMOP
    real(8):: v_hydrol,v_hydrol_ref

    !-- Index of  concentration of detrital particulate organic matter
    iiMOP = nint(indMOP)
    
    !element of MOP:
    elMOP = VAR(iiMOP)%elmt 
    elMOP = trim(adjustl(elMOP))
    
    !-- Index of  cabon concentration of bacteria
    iiC = nint(indBACC)

    ! index of bacterial concentration in cells and in the element of MOP
    iiorg = var(iiC)%idorg
    chain='cell'
    iiCell = f_idorg2id(iiorg,chain)
    iielmt = f_idorg2id(iiorg,elMOP)
    
    ! hydrolysis rate (s-1)
    v_hydrol_ref = 1.0d0/(24 * 3600 * TT)
    nnC=0.2

    !-- Find the specific cell division rate of BAC :
    do l = 1,nscp_cell
        if (mu_growth(l)%idorg == iiCell) then
             igB = l
             exit
        endif
    enddo


  f_hydrol_CX_2026 = 0.0d0
  
  if (elMOP == 'C') then
     Qmin = QminC
     Qmax = QmaxC
     nnX = nnC
  elseif (elMOP == 'N') then
     Qmin = QminN
     Qmax = QmaxN
     nnX = 6
  elseif (elMOP == 'P') then
     Qmin = QminP
     Qmax = QmaxP
     nnX = 6
  endif

 !-- Calculation of the hydrolysis function :
 do k=1,nz_max
   do j=ny_min,ny_max
     do i=nx_min,nx_max
        QC = var(iiC)%conc(i,j,k)/(var(iiCell)%conc(i,j,k)+1d-80)
   ! quota function for bacteria in carbon 
        fQC = ((QmaxC - QC)/(QmaxC-QminC))
   ! Elemental quota and quota function in element X for bacteria
        QX = var(iielmt)%conc(i,j,k)/(var(iiCell)%conc(i,j,k)+1d-80)
        fQX = ((Qmax - QX)/(Qmax-Qmin))

        fQC = max(0.0e0,fQC)
        fQC = min(1.0e0,fQC)
        fQX = max(0.0e0,fQX)
        fQX = min(1.0e0,fQX)

	fQ = max (fQC**nnC,fQX**nnX)
        
	v_hydrol = min(v_hydrol_ref,mu_growth(igB)%val(i,j,k))

!-- test on  the presence of bacteria and on DET availability                 
        if (var(iiCell)%conc(i,j,k) .gt. 100  .and. &
                var(iiMOP)%conc(i,j,k) > 1.1 * v_hydrol * var(iiMOP)%conc(i,j,k) *fQ  *dt_bio + 1.0d-50) then

            f_hydrol_CX_2026(i,j,k) =  v_hydrol * var(iiMOP)%conc(i,j,k) * fQ

        endif
     enddo
   enddo
 enddo

end function f_hydrol_CX_2026
!
