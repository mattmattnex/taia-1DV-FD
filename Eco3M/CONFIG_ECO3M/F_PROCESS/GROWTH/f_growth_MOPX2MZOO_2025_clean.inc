!***************************************************************************
!***************************************************************************
!Copyright or © or Copr. CNRS/IRD/Université de la Méditerranée
!contributor(s) : Melika BAKLOUTI & Vincent FAURE (10/10/2006)
!
!melika.baklouti@univ-amu.fr; vincent.faure@univmed.fr
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
   function f_growth_MOPX2MZOO_2025_clean(ili,jcol,indcell,indC,QC_MIN,QC_MAX,indDETX)
!
!---------------------------------------------------------------------------
! Function that converts the DET consumption rate per MZOO from ingested carbon to X ingested
!
! Melika Baklouti 16/10/2024; last revision  june 2026
!-------------------------------------------------------------------------
    USE mod_eco3m_user
    Implicit None

    !-- function:
    Real(8):: f_growth_MOPX2MZOO_2025_clean(nx_min:nx_max,ny_min:ny_max,nz_max)

    !-- arguments
    Real(8) :: indC,indcell,indDETX,QC_MIN,QC_MAX,indDETC
    Integer :: ili,jcol

    !-- local variables
    Integer :: iorgX,iorgCell
    Integer :: iindX,iiCell,iindDETX,iindDETC
    Integer :: i,j,k,l,ll
    Character(L_VAR)    :: chain
    Character(L_CHAIN)  :: sc,elmt


    !-- DETX concentration index
    iindDETX = nint(indDETX)
    sc = var(iindDETX)%scomp
    elmt='C'
    iindDETC = f_scomp2id(sc,elmt)
    indDETC=iindDETC

    f_growth_MOPX2MZOO_2025_clean=f_growth_MOP2MZOO_2025_clean(ili,jcol,indcell,indC,QC_MIN,QC_MAX,indDETC) * &
                                                        VAR(iindDETX)%conc(:,:,:)  /(VAR(iindDETC)%conc(:,:,:) + 1d-50)

end function f_growth_MOPX2MZOO_2025_clean
!--------------------------------------------------------------------
