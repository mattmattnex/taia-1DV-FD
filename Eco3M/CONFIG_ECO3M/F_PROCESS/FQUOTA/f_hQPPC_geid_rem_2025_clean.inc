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
function f_hQPPC_geid_rem_2025_clean(ili,jcol,indC,Qmin,Qmax)
    !
    !--------------------------------------------------------------------
    ! This function calculates the quota function for remineralization as a funtion
    ! of the carbon cell quota 
    !
    ! This function is based on the formulation provided in Geider et al. (1968) 
    !
    ! Variables,  parameters and their units :
    !-----------------------------------------
    ! indC   : index of carbon concentration 
    !
    ! Ref: Geider et al.(1998) Limnol. & Oceanogr. 43(4):679-694
    !      Baklouti, M. et al. (2006) Progress in Oceanography 71:1-33
    !
    ! Melika Baklouti 01/07/2011; Last version July 2025
    !-------------------------------------------------------------------------

    Implicit None

    !-- function:
    Real(8):: f_hQPPC_geid_rem_2025_clean(nx_min:nx_max,ny_min:ny_max,nz_max)

    !-- arguments
    Real(8) :: indC,Qmin,Qmax
    Integer :: ili,jcol


                f_hQPPC_geid_rem_2025_clean = 1.0d0 -  f_fQPPC_geid_rem_2025_clean(ili,jcol,indC,Qmin,Qmax)



end function f_hQPPC_geid_rem_2025_clean
!--------------------------------------------------------------------
