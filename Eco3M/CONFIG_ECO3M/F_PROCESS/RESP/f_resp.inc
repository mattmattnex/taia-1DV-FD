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
function f_resp_phy_new(ili,jcol,costm,QminC,QmaxC,costg,mumax,QminN,QmaxN,QminP,QmaxP)
!
!--------------------------------------------------------------------
! Cette fonction calcule la respiration en fonction de la croissance brute
!
! Variables et parametres et leurs unites:
!----------------------------------------- 
! costg     = cout de la croissance (proportion carbone incorpore/carbone respire)
!  
! Ref: Thingstad,1987: Utilization of N;P and Organic C by heterotrophic bacteria. I. Outline of a chemostat 
!       theory with a consistent concept of 'maintenance metabolism'. MEPS Vol.35 99-109

!Authors: R. Mauriac, M. Baklouti
!--------------------------------------------------------------------


    Implicit None

    !-- fonction:
    Real(8) :: f_resp_phy_new(nx_min:nx_max,ny_min:ny_max,nz_max)

    !-- arguments
    Real(8)          :: costg,costm,QminC,mumax,QmaxC
    Real(8),optional :: QmaxN,QmaxP,QminN,QminP
    Integer          :: ili,jcol

    !--variables locales:
    Integer            :: i,j,k
    Integer            :: iorg,iiCell,iiC,iiN,iiP
    Real(8)            :: QC,quota,hQ,mu,QN,QP,fQ1,fQ2,fQ3,hQ1,hQ2,hQ3,nn
    Character(L_CHAIN) :: chain1,chain2,chain3,chain4
    !logical :: isnan

    !-- Calcul du quota intracelulaire en carbone
    chain1 = 'cell'
    chain2 = 'C'
    chain3 = 'N'
    chain4 = 'P'
    iorg = VAR(ili)%idorg
    iiCell = f_idorg2id(iorg,chain1)
    iiC    = f_idorg2id(iorg,chain2)
    iiN    = f_idorg2id(iorg,chain3)
    iiP    = f_idorg2id(iorg,chain4)
    nn = 1.00
    do k = 1,nz_max
        do j = ny_min,ny_max
            do i = nx_min,nx_max
                QC = VAR(iiC)%conc(i,j,k)/(VAR(iiCell)%conc(i,j,k)+1d-80)
                if (present (QmaxN)) QN = VAR(iiN)%conc(i,j,k)/(VAR(iiCell)%conc(i,j,k)+1d-80)
                if (present (QmaxP)) QP = VAR(iiP)%conc(i,j,k)/(VAR(iiCell)%conc(i,j,k)+1d-80)

                quota = QC - QminC
                ! Calcul des fQ de croissance
                if (QC .LE. QmaxC) then
                    fQ1 = ((QmaxC - QC)/(QmaxC - QminC))**(1/nn)
                else
                    fQ1 = -(abs((QmaxC - QC) / (QmaxC -QminC)))**(1/nn)  
                endif
                if (present (QmaxN))  then                        
                    if (QN .LE. QmaxN) then
                        fQ2 = ((QmaxN - QN)/(QmaxN - QminN))**(1/nn)
                    else
                        fQ2 = -(abs((QmaxN - QN) / (QmaxN -QminN)))**(1/nn)
                    endif
                else                                    
                    fQ2 = 0.e0
                endif
                if (present (QmaxP))  then                        
                    if (QP .LE. QmaxP) then
                        fQ3 = ((QmaxP - QP)/(QmaxP - QminP))**(1/nn)
                    else
                        fQ3 = -(abs((QmaxP - QP) / (QmaxP -QminP)))**(1/nn)
                    endif
                else
                    fQ3 = 0.e0
                endif                        
                ! Calcul de hQ de croissance

                hQ1 = 1.e0 - fQ1
                hQ2 = 1.e0 - fQ2
                hQ3 = 1.e0 - fQ3

                !-- calcul du quota minimum:
                hQ =  min(1.0e0, hQ1, hQ2, hQ3)
                mu = mumax * max(0.0e0,hQ)

                if (quota > 0.0e0 .and. quota < (QmaxC-QminC)) then
                    f_resp_phy_new(i,j,k) = (costm*quota)*VAR(iiCell)%conc(i,j,k) + mu*costg*QminC
                else
                    f_resp_phy_new(i,j,k) = mu*costg*QminC
                endif

                !IF (ISNAN(f_resp_phy_new(i,j,k))) THEN
                !WRITE(*,*) 'FUNCTION f_resp_phy_new contains NaN',i,j,k
                !STOP
                !ENDIF
            enddo
        enddo
    enddo

end function f_resp_phy_new
