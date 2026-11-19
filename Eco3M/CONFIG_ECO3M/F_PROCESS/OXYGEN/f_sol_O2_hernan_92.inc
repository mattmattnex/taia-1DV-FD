!***************************************************************************
!***************************************************************************
!Copyright or © or Copr. CNRS/IRD/Université de la Méditerranée
!contributor(s) : Melika BAKLOUTI & Vincent FAURE (10/10/2006)
!
!m.baklouti@univmed.fr vincent.faure@univmed.fr
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
       function f_sol_O2_hernan_92()
!
!--------------------------------------------------------------------
!
!  This function provides the solubility (saturation) of oxygen (O2) as a function
!  of  salinity and temperature in sea water at 1-atm pressure of air including 
!  saturated water vapor
!=========================================================================
!
! Variables, parameters and their units:
!-----------------------------------------
!   S = salinity    [PSS]
!   T = temperature [degree C]
!   f_sol_O2_hernan_92 = solubility of O2  [umol/kg] 
! 
! REFERENCE:
!    Hernan E. Garcia and Louis I. Gordon, 1992.
!    "Oxygen solubility in seawater: Better fitting equations"
!    Limnology and Oceanography, 37, pp. 1307-1312.
!
! AUTHOR:  Roberta Hamme from Scripps Inst of Oceanography (rhamme@ucsd.edu)
! O2sol Version 1.1 4/4/2005
!=========================================================================
!
 USE COUPLEUR_PHY_BIO
! USE VAR_GLOBAL
! USE MOD_FCHAIN
  USE  comdynmod

 Implicit None

!-- function:
 Real(8):: f_sol_O2_hernan_92(nx_min:nx_max,ny_min:ny_max,nz_max)

!-- Arguments

!-- Local variables:
 real(8) :: temp_S(nx_min:nx_max,ny_min:ny_max,nz_max),RHO(nx_min:nx_max,ny_min:ny_max,nz_max)
 real(8):: A0_o2,A1_o2,A2_o2,A3_o2,A4_o2,A5_o2,B0_o2,B1_o2,B2_o2,B3_o2,C0_o2
 integer:: i,j,k

! convert T to scaled temperature
temp_S = log((298.15 - TEMP_BIO(:,:,:))/(273.15 + TEMP_BIO(:,:,:)))

! constants from Table 1 of Garcia & Gordon for the fit to Benson and Krause (1984)
A0_o2 = 5.80871 
A1_o2 = 3.20291
A2_o2 = 4.17887
A3_o2 = 5.10006
A4_o2 = -9.86643e-2
A5_o2 = 3.80369
B0_o2 = -7.01577e-3
B1_o2 = -7.70028e-3
B2_o2 = -1.13864e-2
B3_o2 = -9.51519e-3
C0_o2 = -2.75915e-7

! Corrected Eqn (8) of Garcia and Gordon 1992
f_sol_O2_hernan_92(:,:,:) = exp(A0_o2 + A1_o2*temp_S + A2_o2*temp_S**2 + A3_o2*temp_S**3 + A4_o2*temp_S**4 + A5_o2*temp_S**5 &
     + SAL_BIO(:,:,:)*(B0_o2 + B1_o2*temp_S + B2_o2*temp_S**2 + B3_o2*temp_S**3) + C0_o2*SAL_BIO(:,:,:)**2)

do k=1,nz_max
  do j=ny_min, ny_max
    do i= nx_min, nx_max
      RHO(i,j,k) = RAUT(k)
    enddo
 enddo
enddo

! conversion de µmol/kg en mol/l sachant que la masse vol est en kg/m3:
f_sol_O2_hernan_92(:,:,:) = f_sol_O2_hernan_92(:,:,:)  / (RHO(:,:,:)/1000.0) * 1d-6
!write(*,*) TEMP_BIO(1,1,nz_max),f_sol_O2_hernan_92(1,1,nz_max)

END function f_sol_O2_hernan_92
