!***************************************************************************
!***************************************************************************
!Copyright or © or Copr. CNRS/IRD/Université de la Méditerranée
!contributor(s) : Melika BAKLOUTI & Vincent FAURE (10/10/n006)
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
function f_fQPPX_geid_new(ili,jcol,indX,Qmin,Qmax,minCN_ratio,maxCN_ratio,minCP_ratio,maxCP_ratio)
    !
    !--------------------------------------------------------------------
    ! Cette fonction permet de contraidre les rapport C/N et C/P et N/P entre le valeurs 
    ! min_ratio et max_ratio
    !
    ! La fonction de quota est donnée par le modèle de Geider et al. (1968) 
    !
    ! Variables et parametres et leurs unites:
    !-----------------------------------------
    ! indX: indice de la conc. interne en nutriment X dont le rapport est a contraindre
    !
    !
    ! Ref: Geider et al.(1998) Limnol. & Oceanogr. 43(4):679-694
    !      Baklouti, M. et al. (2006) Progress in Oceanography 71:1-33
    !
    ! Romain MAuriac 4/09/2007
    !-------------------------------------------------------------------------

    Implicit None

    !-- fonction:
    Real(8):: f_fQPPX_geid_new(nx_min:nx_max,ny_min:ny_max,nz_max)

    !-- arguments
    Real(8) :: minCN_ratio,maxCN_ratio,indX,Qmin,Qmax
    Real(8) :: minCP_ratio,maxCP_ratio,NC_max,NC_min,PC_max,PC_min
    Integer :: ili,jcol,iindX

    !--variables locales:
    Integer :: identite
    Integer :: iiC,i,j,k,iiCell,iiX,iiN,iiP
    Real(8) :: Q1,hQ1,fQ1,Q2,hQ2,fQ2,Q3,hQ3,fQ3,nn
    Character(L_VAR)::chain,chain2,chain3,chain4

    ! Indice de l'element a contraindre
    iindX = nint(indX)
    ! Numero de l'Organismes
    identite = var(iindX)%idorg
    ! Indice de pour le carbone et les cellules
    chain = 'C'
    chain2 = 'N'
    chain3 = 'P'
    chain4 = 'cell'
    iiC = f_idorg2id(identite,chain)
    iiN =  f_idorg2id(identite,chain2)
    iiP =  f_idorg2id(identite,chain3)
    iiCell = f_idorg2id(identite,chain4)
    iiX = iindX

    ! Conversion de C/X min et max en X/C min et max
    NC_max = 1.d0/minCN_ratio
    NC_min = 1.d0/maxCN_ratio
    PC_max = 1.d0/minCP_ratio
    PC_min = 1.d0/maxCP_ratio

    nn = 0.06
    !-- Calcul de la fonction de limitation :
    do k=1,nz_max
        do j=ny_min,ny_max
            do i=nx_min,nx_max ! good position
                Q3 = var(iiX)%conc(i,j,k) /  (var(iiCell)%conc(i,j,k)+1.d-80)

                ! Calcul de fQ d'uptake
                if (Q3 <= Qmax) then
                    fQ3 = ((Qmax - Q3) / (Qmax - Qmin))**nn      
                else
                    fQ3 = -(abs((Qmax - Q3) / (Qmax -Qmin)))**(nn);
                endif

                hQ3 = 1.d0 - fQ3

                if (hQ3 < 0.0d0) hQ3 = 0.0d0

                !-- cas de l'azote
                if (iiX == iiN) then
                    if (Q3 <= Qmin) then
                        f_fQPPX_geid_new(i,j,k) = 0.0d0
                    elseif (Q3 >= Qmax) then
                        f_fQPPX_geid_new(i,j,k) = hQ3
                    else
                        Q1 = var(iiN)%conc(i,j,k) / (var(iiC)%conc(i,j,k)+1.d-80) !quota C:N
                        if (Q1 .LE. NC_max) then
                            fQ1 = ((NC_max - Q1) / (NC_max -NC_min))**nn
                        else
                            fQ1 =  -(abs((NC_max - Q1) / (NC_max - NC_min)))**(nn);
                        endif
                        hQ1 = 1.d0 - fQ1
                        if (hQ1 < 0.0d0) hQ1 = 0.0d0
                        f_fQPPX_geid_new(i,j,k) = hQ1
                    endif
                elseif (iiX == iiP) then
                    if (Q3 <= Qmin) then
                        f_fQPPX_geid_new(i,j,k) = 0.0d0
                    elseif (Q3 >= Qmax) then
                        f_fQPPX_geid_new(i,j,k) = hQ3
                    else
                        Q2 = var(iiP)%conc(i,j,k) / (var(iiC)%conc(i,j,k)+1.d-80) !quota C:P
                        if (Q2 .LE. PC_max) then
                            fQ2 = ((PC_max - Q2) / (PC_max - PC_min))**nn
                        else
                            fQ2 = -(abs((PC_max - Q2) / (PC_max - PC_min)))**(nn);
                        endif
                        hQ2 = 1.d0 - fQ2
                        if (hQ2 < 0.0d0) hQ2 = 0.0d0
                        f_fQPPX_geid_new(i,j,k) = hQ2 
                    endif
                else
                    STOP "pb dans l'utilisation de la fonction f_fQPPX_geid_new"
                endif
                if (f_fQPPX_geid_new(i,j,k)>1.0) THEN 
                    f_fQPPX_geid_new(i,j,k) = 1.d0
                endif
            enddo
        enddo
    enddo

    !WRITE(*,*) 'FUNCTION f_fQPPX_geid_new '
end function f_fQPPX_geid_new
!--------------------------------------------------------------------
