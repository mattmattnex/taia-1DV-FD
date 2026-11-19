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
function f_hQPPX_geid_2026(ili,jcol,indX,Qmin,Qmax,minCX_ratio,maxCX_ratio)
    !
    !--------------------------------------------------------------------
    ! Quota function governed by intracellular ratios and quotas.
    ! This function is based on the formulation given by Geider et al. (1998) 
    !
    ! Variables,   parameters and units:
    !-----------------------------------
    ! indX: index of the variable corresponding to the intracell. conc. in element X 
    ! Qmin: minimum intracellular quota in element X
    ! Qmax: maximum intracellular quota in element X
    !
    ! Ref: Geider et al.(1998) Limnol. & Oceanogr. 43(4):679-694
    !      Baklouti, M. et al. (2006) Progress in Oceanography 71:1-33
    !
    ! Romain Mauriac 4/09/2007; Melika Baklouti 2020;2024j
    !-------------------------------------------------------------------------

    Implicit None

    !-- function:
    Real(8):: f_hQPPX_geid_2026(nx_min:nx_max,ny_min:ny_max,nz_max)

    !-- arguments:
    Real(8) :: indX,Qmin,Qmax
    Real(8),optional :: minCX_ratio,maxCX_ratio
    Integer :: ili,jcol,iindX

    !-- local variables:
    Integer :: identite
    Integer :: iiC,i,j,k,iiCell,iiX
    Real(8) :: QXC,hQXC,fQXC,QX,hQX,fQX,nn
    Real(8) :: XC_max,XC_min
    Character(L_VAR)::chain,chain2

    ! Index of the biogenic element
    iindX = nint(indX)

    ! Organism number:
    identite = var(iindX)%idorg

    ! Index for carbon and cellular concentrations 
    chain = 'C'
    chain2 = 'cell'
    iiC = f_idorg2id(identite,chain)
    iiCell = f_idorg2id(identite,chain2)
    iiX = iindX

    ! Conversion 
    XC_max = 1.d0/minCX_ratio
    XC_min = 1.d0/maxCX_ratio

    nn = 0.06
    
  !-- Calculation of the limitation function:
    do k=1,nz_max
        do j=ny_min,ny_max
            do i=nx_min,nx_max ! good position
                QX = var(iiX)%conc(i,j,k) /  (var(iiCell)%conc(i,j,k)+1.d-80)

                fQX = ((Qmax - QX) / (Qmax - Qmin))      
                fQX = max(0.0d0,fQX)
                fQX = min(1.0d0,fQX)

                hQX = 1.d0 - fQX**nn

                if (QX <= Qmin) then
                     f_hQPPX_geid_2026(i,j,k) = 0.0d0
                elseif (QX >= Qmax) then
                     f_hQPPX_geid_2026(i,j,k) = hQX
                else
                     QXC = var(iiX)%conc(i,j,k) / (var(iiC)%conc(i,j,k)+1.d-80) ! X:C ratio
                     fQXC = ((XC_max - QXC) / (XC_max -XC_min))
                     fQXC = min(1.0d0,fQXC)
                     fQXC = max(0.0d0,fQXC)

                     hQXC = 1.d0 - fQXC**nn
                     f_hQPPX_geid_2026(i,j,k) = max(hQXC,hQX)
                endif
            enddo
        enddo
    enddo

    !WRITE(*,*) 'FUNCTION f_hQPPX_geid_2026 '
end function f_hQPPX_geid_2026
!--------------------------------------------------------------------
