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
!----------------------------------------------------------------------------------------------------------
    SUBROUTINE eco3m_mod_init

    !> Subroutine that defines the different processus models available, from the reading of the
    !! fileconfmod file (default name is "modele.def")
    !! The different processus model are then referenced into the PROC_MOD matrix, which associates
    !! processus names and ids (integer numbers)
    !!
    !! \warning All the processus names appearing in the fileconfmod file ("modele.def" by default) must be
    !! declared in the "mod_process.F90" file.
    !! \date 2007-07-06
    !! \author Melika Baklouti, Vincent Faure, Nicolas Barrier
!----------------------------------------------------------------------------------------------------------
    use eco3m_string, only: l_chain, l_var, f_chain ! Module for string manipulation
    use mod_eco3m_files
    use mod_eco3m
    
    implicit none

        ! Local variables
        integer            :: errlec,istat
        integer            :: i,j,itemp, k
        character(L_CHAIN) :: chaine,chaine_tempo
        logical            :: file_exist

  ! Initialisation of the number of processes defined in the modele.def file
        nb_proc =0
  ! Opening the modele.def file 
        filename_confmod = trim(adjustl(eco3m_root_dir))//trim(adjustl(filename_confmod))
        inquire(file=filename_confmod, exist=file_exist)
        if (.not.(file_exist)) then
            write(*,*) "****************************************"
            write(*,*) "The ", trim(filename_confmod), " file does not exist"
            write(*,*) "This program will stop"
            write(*,*) "****************************************"
            stop
        end if

        open(file_confmod_id, file=filename_confmod)
  ! Reads the modele.def file
        !  reads the file in order to count the number of lines
        !  corresponding to an effective process
        do
            read(file_confmod_id, *, iostat=errlec) chaine
            if (errlec /=0) then
                write(*,*) "Reading problem with the modele.def file.",& 
                &" The program will be stopped"
                stop
            endif
            if (chaine(1:1)=='#') then
                cycle
            elseif (chaine(1:4)=='!fin' .OR. chaine(1:4)=='!FIN'.OR. chaine(1:4)=='!Fin') then
                exit
            else
                nb_proc = nb_proc + 1  ! if good line, iteration of the nb_proc value
            endif
        end do


        ! Creation of the process matrix, which links the name, ID and number
        ! of parameters associated with a process
        if (nb_proc /=0)  Allocate(PROC_MOD(nb_proc))

        ! We go back to the beginning of the file, in order
        ! to now define the elements of the PROC_MOD array
        rewind(file_confmod_id)

        ! integer i is the index of the current process
        ! initialised at 0
        i=0

        ! we loop until we have the right number of processes
        do while (i < nb_proc)   
            read(file_confmod_id, *, iostat=errlec) chaine ! reading of a string corresping to a process
            if (errlec /=0) then
                write(*,*) "Reading problem with the modele.def file.",& 
                &" The program will be stopped"
                stop
            endif
            if (chaine(1:1)=='#') cycle

            i = i+1 
            PROC_MOD(i)%idproc = i  ! ID of the process

            ! Name of the process as used in the config.ini file
            chaine_tempo=f_chain(chaine,1,':')  
            read(chaine_tempo,*)PROC_MOD(i)%nomproc

            ! Name of the f_process function used to compute the sub-flux
            chaine_tempo=f_chain(chaine,2,':')
            read(chaine_tempo,*)PROC_MOD(i)%nomsub

            ! Number of arguments (of parameters)
            chaine_tempo=f_chain(chaine,3,':')
            read(chaine_tempo,*)PROC_MOD(i)%nbpar

            itemp = PROC_MOD(i)%nbpar
            ! Allocation of the array of parameter names
            if (associated (PROC_MOD(i)%nompar)) NULLIFY(PROC_MOD(i)%nompar)
            Allocate(PROC_MOD(i)%nompar(itemp),STAT=istat)
            if (istat /= 0) write(*,*) "Problem with the allocation of the",&
                &" PROC_MOD(",i,")%nompar array"

            ! Now we fill the values of the PROC_MOD(',i,')%nompar array
            do j=1,PROC_MOD(i)%nbpar
                chaine_tempo=f_chain(chaine,3+j,':')
                read(chaine_tempo,*)PROC_MOD(i)%nompar(j)
            enddo

        end do

        ! WRITING LOG
#ifdef INI
        write(file_CR_id,*) "================ Eco3M: Initialisation of the",&
            &" model configuration ================"
        write(file_CR_id,*) "Configuration file: ", trim(filename_confmod)
        write(file_CR_id,*)
        write(file_CR_id,*) "Number of model processes: ", nb_proc
        write(file_CR_id,*) "Description of the PROC_MOD array"
        write(file_CR_id,*) "Id : Name : Function name : Number of par. : Par."
        do i=1,size(PROC_MOD)
            write(file_CR_id,'(I3,$)') PROC_MOD(i)%idproc
            ! write(file_CR_id,*)' : ',trim(adjustl(PROC_MOD(i)%nomproc)),' : ',&
            write(file_CR_id,*)' : ',trim(adjustl(PROC_MOD(i)%nomproc)),' : ',&
                trim(adjustl(PROC_MOD(i)%nomsub)),       &
                PROC_MOD(i)%nbpar,' : ',                 &
                (trim(adjustl(PROC_MOD(i)%nompar(k))),'|', k=1,size(PROC_MOD(i)%nompar))
        enddo 
        write(file_CR_id,*)
#endif
        close(file_confmod_id)    
    
    End Subroutine eco3m_mod_init

