!***************************************************************************
!***************************************************************************
!Copyright or © or Copr. CNRS/IRD/Université de la Méditerranée
!contributor(s) : Melika BAKLOUTI & Vincent FAURE (10/10/2006)
!
!melika.baklouti@univ-amu.fr; 
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
!> Module that contains the functions related to the handling of state variables ID.  
!! \author Melika Baklouti, Vincent Faure, Nicolas Barrier
!----------------------------------------------------------------------------------
    MODULE mod_eco3m_id_extract
!----------------------------------------------------------------------------------
    use mod_eco3m , only: VAR, PROC_MOD  ! Module containing the VAR matrix
    use eco3m_string, only: l_chain, l_var  ! Module for string manipulation

    Implicit None

CONTAINS

!----------------------------------------------------------------------------------
    Integer  Function f_scomp2id (scomp_in,elmt_in)

    !> Finds the ID within the state variable array VAR corresponding to a given element within
    !! a given sub-compartment
    !! For instance, we seek the ID of the variable corresponding to the N concentration 
    !! of the diatom sub-compartment!
!----------------------------------------------------------------------------------
        ! arguments
        character(L_CHAIN), intent(in) :: elmt_in !< Element name
        character(L_CHAIN), intent(in) :: scomp_in  !< Sub-compartment name

        ! local variables
        integer              :: i

        f_scomp2id = 0
        i=1
        do while (i<=size(VAR))
            if(trim(adjustl(VAR(i)%scomp)) == trim(adjustl(scomp_in)) .AND. & 
                trim(adjustl(VAR(i)%elmt)) == trim(adjustl(elmt_in))) then
                f_scomp2id = i
                exit
            endif
            i=i+1
        End do

        Return 
    End Function f_scomp2id

!----------------------------------------------------------------------------------
    Integer  Function f_idorg2id(idorg_in,elmt_in)

    !> Finds the ID within the state variable array VAR corresponding to the variable
    !! which biomass is expressed in "elmt_in", for the sub-compartment of ID "idorg"
    !! For instance, we seek nitrogen N of the diatom sub-compartment!
!----------------------------------------------------------------------------------

        ! arguments
        character(L_VAR), intent(in) :: elmt_in !< Element name
        integer, intent(in) :: idorg_in !< Sub-compartment integer ID
        ! local variables
        integer              :: i, long, long_in

        f_idorg2id = 0
        i=1

        do while (i<=size(VAR))  
            long = len_trim(trim(adjustl(VAR(i)%elmt)))
            long_in = len_trim(trim(elmt_in))

            if (VAR(i)%idorg == idorg_in .AND. &
                (VAR(i)%elmt(1:long) == elmt_in(1:long_in))) then
                f_idorg2id = i
                exit
            endif
            i=i+1
        enddo

        return 
    End Function f_idorg2id

!----------------------------------------------------------------------------------
    Function f_idorg2id_vect (idorg_in, elmt_in, size_idorg)
    !> Finds the ID *array* within the state variable array VAR corresponding to the variables
    !! whose biomasses are expressed in the element "elmt_in", for the *array* of sub-compartment IDs "idorg_in"
    !! \note It is a vectorial form of the f_idorg2id function
!----------------------------------------------------------------------------------

        ! arguments
        character(L_VAR), intent(in) :: elmt_in !< Element name
        integer, intent(in) :: size_idorg !< Size of the idorg_in array
        integer, intent(in) :: idorg_in(size_idorg) !< Array integer of sub-compartments IDs

        ! output
        Integer              :: f_idorg2id_vect(size_idorg)

        ! local variables
        integer              :: i, long, long_in,ll

        f_idorg2id_vect = 0

        do ll=1,size_idorg
            i=1
            do while (i<=size(VAR))  
                long = len_trim(trim(adjustl(VAR(i)%elmt)))
                long_in = len_trim(trim(elmt_in))

                if (VAR(i)%idorg == idorg_in(ll) .AND. & 
                    (VAR(i)%elmt(1:long) == elmt_in(1:long_in))) then
                    f_idorg2id_vect(ll) = i
                    exit
                endif
                i=i+1
            enddo
        enddo
        return 
    End Function f_idorg2id_vect
    
!----------------------------------------------------------------------------------
    Integer function f_proc2id(nom_proc)
    !> This function returns the process ID number associated with a process name.
    !! If the input name is not valid, it stops the program
!----------------------------------------------------------------------------------

        Implicit none

        Character(L_VAR), intent(in) :: nom_proc !< Name of the processus

        ! local variables
        character(L_VAR) :: tempo
        integer          :: i

        f_proc2id=0
        do i=1,size(PROC_MOD)
            tempo = PROC_MOD(i)%nomproc
            ! write(*,*) 'tempo ', tempo
            if(tempo==trim(adjustl(nom_proc))) then
                f_proc2id = PROC_MOD(i)%idproc
                exit 
            end if
        end do

        if(f_proc2id==0) then
            write(*,*)"this process is not present in your", &
                &" model :",nom_proc
            write(*,*) "This  program will be stopped"
            read(*,*)
            stop
        endif

    End function f_proc2id

END MODULE mod_eco3m_id_extract
