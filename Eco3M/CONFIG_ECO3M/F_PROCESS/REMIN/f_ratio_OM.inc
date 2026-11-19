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
         function f_ratio_OM(ili,jcol,iOM_X,iOM_C)
!
!---------------------------------------------------------------------------
! Calculates the X:C ratio in organic matter
!
! Variables et parametres et leurs unites:
!-----------------------------------------
! iOM_X  : index of conc. in element X in organic matter
! iOM_C  : index of conc. in element C in organic matter
!
!
! M Baklouti  17/02/2025
!-------------------------------------------------------------------------

    Implicit None

    !-- fonction:
    Real(8):: f_ratio_OM(nx_min:nx_max,ny_min:ny_max,nz_max)

    !-- arguments
    Real(8) :: iOM_X,iOM_C
    Integer :: ili,jcol

    !--local variables :
    integer:: iiOMX, iiOMC,i,j,k
    
    
    iiOMX = nint(iOM_X)
    iiOMC = nint(iOM_C)


    !-- Calcul de la fonction de remineralisation:
 do k=1,nz_max
   do j=ny_min,ny_max
     do i=nx_min,nx_max
        f_ratio_OM(i,j,k) = var(iiOMX)%conc(i,j,k) / (var(iiOMC)%conc(i,j,k) + 1.d-80) 
     enddo
   enddo
 enddo

end function f_ratio_OM
!
