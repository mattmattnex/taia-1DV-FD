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
function f_extinc_morel88(Chl_tot)
!    
!-------------------------------------------------------------------------
!> This function computes the light extinction coefficients through diffusion according to the
!! \cite Mor88 model, which takes into account the absorption by Chlorophyll A and by Pheotinine A.
!! It also includes dispersion effects (scattering):
!!
!! \f[
!! K_{ext} =  K_{bio}
!! \f]
!! \f[
!! K_{bio} = 0.0518\times Chl_{tot}^{1-0.572}
!! \f]
!!
!! \f$K_{bio}\f$ combines the diffuse attenuation coefficient in pure marine water and the light attenuation due to Chl absorption
!! \ingroup lightex
!! \warning Valid only if \f$0.03 mg.m^{-3} < Chl_{tot} < 30 mg.m^{-3}\f$ (see \cite Mor88
!! , page 10,757)
!! Ref:  Morel, A. et al. (1988) J. Geophys. Res., 93(C9):10,749:10,768
!! \warning Chlorophyll concentration \f$Chl_{tot}\f$ must be in \f$g.m^{-3}\f$
!!
!! \author Melika Baklouti 
!! \date 2017-09-25
!-------------------------------------------------------------------------
!
    Implicit None

    ! function:
    Real(8):: f_extinc_morel88(nx_min:nx_max, ny_min:ny_max, 1:nz_max) 

    !-- arguments
    real(8):: Chl_tot(nx_min:nx_max, ny_min:ny_max, 1:nz_max)

    ! local variables
    integer :: i, j, k

    do k = 1, nz_max
        do j = ny_min, ny_max
            do i = nx_min, nx_max
              
              f_extinc_morel88(i,j,k) = 0.d0
              ! extinction coefficient through diffusion by biological particles (Chl exprime in mg m-3)
              if (Chl_tot(i,j,k) > 0.d0)  f_extinc_morel88(i,j,k) =  0.0518/(1-0.572) * (Chl_tot(i,j,k)*1.d6)**(1-0.572)  
           enddo
       enddo
   enddo

end function f_extinc_morel88
!--------------------------------------------------------------------
