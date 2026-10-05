subroutine ATMOCEOXY(idO2,deltaF1DV)
Use comrunmod
Use comdynmod
Use mod_eco3m_vartypes
Use mod_eco3m
Use eco3m_string, only : L_CHAIN
Use mod_eco3m_id_extract, only:f_scomp2id

Real(8),dimension(nx_min:nx_max,ny_min:ny_max)::Fd_o2,Fu_o2,deltaF
Real(8)::deltaF1DV
character(L_CHAIN)::scomp,elmt
INTEGER::idO2
REAL(8):: Tc,scho2,kgo2,rho,solo2,Patm,salprac,tempis

!-- chemical constants
REAL(8),parameter ::   oxyco  = 1. / 22.4144  ! converts from liters of an ideal gas to moles
REAL(8),parameter ::   atcox  = 0.20946         ! units atm
REAL(8),parameter ::   o2atm  = 1. / ( 1000. * 0.20946 )

rho= RAUT(1)  !!! densité potentielle au niveau de la surface
salprac = SN(1) * 35.0 / 35.16504 ! salinité pratique en surface
zpres = rho / 1000.

a1  = 0.04 * ( 1.0 + 0.185 * TN(1) + 0.035 * (salprac-35.0)) 

a2  = 0.0075 * ( 1.0 - TN(1)/30.0)

!!TN(1):SST , SN(1):SSS, TN:temperature now


tempis =  TN(1) - a1*zpres + a2 * zpres**2

!!!convert T to scaled temperature(Tg)
Tg=log((298.15-tempis)/(tempis+273.15))

zoxy  = 2.00856 + 3.22400 * Tg + 3.99063 * Tg**2 + 4.80299 * Tg**3   &
                        &       + 9.78188e-1 * Tg**4 + 1.71069 * Tg**5 + salprac * ( -6.24097e-3   &
                        &       - 6.93498e-3 * Tg - 6.90358e-3 * Tg**2 - 4.29155e-3 * Tg**3 )   &
                        &       - 3.11680e-7 * salprac**2


! solubilité de l'oxy en  mol/(L atm)
solo2 = ( EXP( zoxy ) * o2atm ) * oxyco * atcox     ! mol/(L atm)


Patm=1.013!*(1e5) !!!pression atmosphérique, constante pour l'instant
!! elle doit être exprimée en atmosphères vu que dans NEMOi,  on a:
! p4zflx.F90:      zpco2atm(ji,jj) = satmco2(ji,jj) * ( patm(ji,jj) - zvapsw )
!      or  zvapsw = EXP(24.4543 - 67.4509*(100.0/ztkel) - 4.8489*LOG(ztkel/100) - 0.000544*zsal)
! = loi de Weiss & Price 1980 donnant la pression de vapeur saturante en atm (j'ai vérifié ds l'article d'origine)
!!! VERIFIER DANS NEMO 3D QUE patm EST BIEN UNE PRESSION EN ATM !!!

Tc=min(35.0,TN(1))
scho2=     1920.4 - 135.6  *  Tc + 5.2122*(Tc**2)- 0.109390*(Tc**3) + 0.0009377*(Tc**4) !!!nb de schmidt

! wind=SQRT(SQRT(taux*taux+tauy*tauy)/(drag*rauair) !!!module de la vitesse du vent 
! coefficient de transfert à l'interface air-mer (vitesse piston convertie de cm/h en m/s)
kgo2 = 0.251*(wind**2)*(0.01/3600)*SQRT(660./scho2)

!!!!!!!!!!!!!!!!!!
 !!!  Flux d'oxy d'atm vers ocean
Fd_o2=Patm*solo2*kgo2

#ifdef MODTEST
!write(*,*) 'Fd_o2=',Fd_o2
#endif 

! Determination de l'indice correspondant à la variable O2
id02 = 0
scomp ='OXY'
elmt='O'
idO2 = f_scomp2id (scomp,elmt)
!!!concentration de l'oxy au niveau de z=0(surface) ???reste à trouver la concentration de o2

IF (idO2 /= 0 .and. id02 <=nbvar) THEN
   Fu_o2=VAR(idO2)%conc(:,:,1) * kgo2      
ENDIF

#ifdef MODTEST
!write(*,*) 'Fu_o2',Fu_o2
#endif 
deltaF=Fd_o2-Fu_o2

deltaF1DV = deltaF(1,1)
!!!!!!
return
end

