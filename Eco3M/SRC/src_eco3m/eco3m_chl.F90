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
!! \author Melika Baklouti
!------------------------------------------------------------------------------------------------------
    subroutine eco3m_chl_init
!------------------------------------------------------------------------------------------------------

    use mod_eco3m
    use mod_eco3m_files, only: file_CR_id
    use mod_eco3m_id_extract
    implicit none

    !> Initialises all the variables relative to Chlorophyll. 
    !!
    !! - Reading of the  Chl:C ratio settings in the Eco3M namelist (see namelist for details)
    !! - logical variable *CHL_C_BOOL* (defined in the eco3m_namelist): true if the Chl:C varie and false otherwise 
    !! - If *CHL_C_BOOL* and *CHL_C_FUNC* are true (i.e. we use a variable and function-defined 
    !!     Chl:C ratio),  the Chl:C subprocess name is read in the config.ini file.
    !!     In INI mode, we also write the "calc_ChlC.inc" file.
    !! - Initialisation (allocation and attribution) of the index of phytoplankton Chl elements (*iphy_Chl*, if CHL_C_BOOL is true) 
    !!     or C element (*iphy_C*, if CHL_C_BOOL is false)
    !! - Allocation of the total chlorophyll array (*Chl_tot*)
    !! - If CHL_C_BOOL is True, allocation of the Chl:C ratio (allocation of the *CHL_C* array and of the *CHL_C\%val*, *CHL_C\%idorg*
    !!  array for each phyto elements)
    !!
    !! \author Melika Baklouti

   character(l_chain) :: elmt
   integer :: k
   integer :: istat


    ! Read the chl-related text in the config.ini file
            call fscanf_chl
#ifdef INI
    ! In INI mode, if Chl_C_Func == True, the Chl:C ratio is calculated through a mathematical 
    ! function which name has to be provided in the config.ini file.
    ! The mathematical expression of Chl:C is written in the file calc_ChlC.inc which is included in the model if the boolean Chl_C_Func is true.
    ! Otherwise, an empty file calc_ChlC.inc  has to be created for the model for technical considerations
           call fwritef_chlc_inc
#endif
        
    ! If the model includes Phytoplankton groups, we allocate the table containing the indices of 
    ! the corresponding Chl concentrations (if Chl:C varies) or C conc. (if Chl:C is constant)
        if (nscp_phy /=0) then
            if (CHL_C_BOOL) then
                allocate(iphy_Chl(nscp_phy))
                iphy_Chl = 0
                elmt = 'Chl'
                iphy_Chl = f_idorg2id_vect(iscp_phy,elmt,nscp_phy)
            else
                allocate(iphy_C(nscp_phy))
                iphy_C = 0
                elmt = 'C'
                iphy_C = f_idorg2id_vect(iscp_phy,elmt,nscp_phy)
                allocate(iphy_Chl(nscp_phy))
                iphy_Chl = 0
            endif
        endif

#ifdef CALC
        
        ! If the Chl:C ratio  varies, allocation of the Chl_C array value for
        ! each phytoplankton element
        If (CHL_C_BOOL) then
            Allocate(CHL_C(nscp_phy))
            do k = 1, nscp_phy
             Allocate(Chl_C(k)%val(nx_min:nx_max,ny_min:ny_max,1:nz_max),stat = istat)
             if (istat /=0) write(*,*)"problem with the allocation of Chl_C(",k,")"
                CHL_C(k)%val(:,:,:)=0.d0
                Chl_C(k)%idorg=iscp_phy(k)
            enddo
        Endif
#endif
        
        
        ! WRITING LOG FILE
#ifdef INI
        write(file_CR_id,*) "================= Eco3M: Chl initialisation ======================="
        write(file_CR_id,*) "Variable Chl:C ratio: ", chl_c_bool
        if(chl_c_bool) then
            if(chl_c_func) then
                write(file_CR_id,*) "Chl:C ratio computed by the function ", &
                    trim(proc_mod(chl_c_par%idproc)%nomsub), " function"
                write(file_CR_id,*) "Parameters: "
                do k=1, size(proc_mod(chl_c_par%idproc)%nompar)
                    write(file_CR_id,"(A,A,F10.5)") proc_mod(chl_c_par%idproc)%nompar(k), &
                        " = ", chl_c_par%valpar(k)
                end do
            end if
        else
            write(file_CR_id,"(A, F10.5)") "Constant ratio: ", Chl_C0
        end if
        if (nscp_phy /= 0) then
            if(Chl_C_Bool) then
                write(file_CR_id,*) "Index of phyto Chl variables: "
                write(file_CR_id,*) (iphy_chl(k), k=1, nscp_phy)
            else
                write(file_CR_id,*) "Index of phyto C elements for calculation of Chl conc.: "
                write(file_CR_id,*) (iphy_c(k), k=1, nscp_phy)
            end if
        end if
        write(file_CR_id,*)
#endif

    end subroutine eco3m_chl_init


!-------------------------------------------------------------------------------------------------------------
#ifdef CALC

    subroutine compute_chl_tot(Chl_tot)

    !> **Only if there is phytoplankton in the model.** Calculation of total Cholorophyll. 
    !!
    !! - If *CHL_C_BOOL* is True, it sums the concentration of the Chlorophyll concentrations over all the 
    !! pythoplankton compartments. 
    !! - If *CHL_C_BOOL* is False, it sums the concentration of the Carbon elements 
    !! over all the phytoplankton concentrations, multiplied by the constant Chl:C ratio (*CHL_C0* variable, defined in the
    !! namelist)
    !!
    !! \author Melika Baklouti, N Barrier
!-------------------------------------------------------------------------------------------------------------
    use mod_eco3m
    implicit none
    integer :: k, i, j, nn
    Real(8):: Chl_tot(nx_min:nx_max,ny_min:ny_max,1:nz_max)


    Chl_tot = 0.d0
    do k = 1, nz_max
        do j = ny_min, ny_max
            do i = nx_min, nx_max
                if (nscp_phy /=0) then
                    do nn = 1, nscp_phy
                        if (CHL_C_BOOL .and. maxval(iphy_Chl) /= 0) then
                            ! If variable Chl:C ratio, we directly take Chl elements
                            Chl_tot(i, j, k) = Chl_tot(i, j, k) + &
                                VAR(iphy_Chl(nn))%conc(i, j, k) 
                        elseif(.NOT. (CHL_C_BOOL) .and. maxval(iphy_C) /= 0) then
                            ! If constant ratio, we convert Carbon into Chl
                            Chl_tot(i, j, k) = Chl_tot(i, j, k) + &
                                VAR(iphy_C(nn))%conc(i,j,k)*CHL_C0
                        endif
                    enddo
                endif
            enddo
        enddo
    enddo
    end subroutine compute_chl_tot
#endif

!-------------------------------------------------------------------------------------------------------------
#ifdef CALC
    
    SUBROUTINE compute_chl_c_ratio

    !> **Only if *CHL_C_BOOL* is True**. Computes the Chl:C ratio for each phytoplankton compartment  using the
    !! phytoplankton concentrations at previous time-step. 
    !!
    !!  - If *CHL_C_FUNC* is False, then the ratio is computed by using the concentrations of Chl and C 
    !!  - If *CHL_C_FUNC* is True, then the function defined in the "config.ini" is called through the "calc_ChlC.inc" file.
    !!
    !! \date 2007-06-27
    !! \author Melika Baklouti, Vincent Faure
!-------------------------------------------------------------------------------------------------------------
     use mod_eco3m
     use mod_eco3m_id_extract, only:f_idorg2id
     Implicit None
   !-- variables globales:
     Character(L_VAR) :: chaine,chaine2
     integer :: i, ichl, ic ,iorg

        ! If Chl is a state variable, we simply compute the ratio
        if ((CHL_C_BOOL .EQV. .TRUE.) .AND. (CHL_C_FUNC .EQV. .FALSE.)) then
            if (Allocated (CHL_C)) then  
                ic = 0
                ichl = 0
                do i = 1, nscp_phy 
                    iorg = CHL_C(i)%idorg
                    chaine = 'C'
                    ic = f_idorg2id(iorg,chaine)
                    chaine2 = 'Chl'
                    ichl = f_idorg2id(iorg,chaine2)
                    CHL_C(i)%val(:,:,:) = var(ichl)%conc(:,:,:)/(var(ic)%conc(:,:,:) + 1.d-15)
                enddo
            endif
        ! If Chl is not a state variable,  we call the Chl:C function through the calc_ChlC.inc file 
        elseif  ((CHL_C_BOOL .EQV. .TRUE.) .AND. (CHL_C_FUNC .EQV. .TRUE.)) then
            include "calc_ChlC.inc"
        endif
    End Subroutine compute_chl_c_ratio
#endif
!-------------------------------------------------------------------------------------------------------------

    subroutine fscanf_chl

    !> **Only if  both *CHL_C_BOOL* and *CHL_C_FUNC* are true**
    !! Reads in the *config.ini* file the mathematical function which is used for the Chl:C ratio calculation, and its parameters
    !! (initialisation of the *CHL_C_PAR* array). 
    !! \note The function must also be defined in the *model.def* file and
    !  declared in the eco3m_fprocess module.
    !
    !! \author Melika Baklouti
!-------------------------------------------------------------------------------------------------------------
    use mod_eco3m, only: CHL_C_BOOL,CHL_C_FUNC,CHL_C_PAR,PROC_MOD,CHL_C0
    use mod_eco3m_files, only:file_config_id
    use eco3m_string
    use mod_eco3m_id_extract
    implicit none

    Character(L_CHAIN) :: chaine, tempo, tempo2, tempopo
    integer :: errlec,istat
    integer :: idtemp, inbelmt
    integer :: k
    integer :: CHL_C_VAL 
errlec = 0
do
 read(file_config_id,*,iostat=errlec) chaine
#ifdef MODTEST
     write(*,*) chaine
#endif 
 if (errlec /=0) stop 'pb when reading the Chl:C ratio in the config file'
 if (chaine(1:1)=='#') cycle
 tempo=f_chain(chaine,1,':')
! -- Reads the CHL_C_VAL value in the config file
 Read(tempo,*) CHL_C_VAL
 if ( CHL_C_VAL == 0) then  ! Chl:C ratio is constant
       CHL_C_BOOL = .FALSE.
       tempo=f_chain(chaine,2,':')
       Read(tempo,*)CHL_C0
 elseif ( CHL_C_VAL == 1) then
      CHL_C_BOOL = .TRUE.
      tempo=f_chain(chaine,2,':')
      if (tempo == '') then
          CHL_C_FUNC = .FALSE.
      else
          CHL_C_FUNC = .TRUE.

!-- looks for the process index associated with the function calculating the
!   Chl:C ratio:
          tempo= f_chain(tempo,1,'(')
          tempo = trim(adjustl(tempo))

#ifdef MODTEST
          write(*,*) 'tempo =',tempo
#endif

         tempo2=f_chain(chaine,2,'(')
         tempo2 = tempo2(1:len_trim(tempo2)-1)

#ifdef MODTEST
     write(*,*) 'tempo2 =',tempo2
#endif
          idtemp = f_proc2id(tempo)
          CHL_C_PAR%idproc = idtemp
          inbelmt= PROC_MOD(idtemp)%nbpar
          Allocate(CHL_C_PAR%valpar(inbelmt),STAT=istat)
          if (istat /= 0) write(*,*) 'pb with the allocation of CHL_C_PAR%valpar'
#ifdef MODTEST
          write(*,*) 'alloc of IRR_PAR%valpar at ',inbelmt
#endif
          do k= 1, inbelmt
            tempopo=f_chain(tempo2,k,'>')
            Read(tempopo,*)CHL_C_PAR%valpar(k)
          enddo
      endif
 else
       stop 'Problem with the CHL:C ratio'
     endif
     exit
 enddo
#ifdef MODTEST
     write(*,*) 'CHL_C_FUNC', CHL_C_FUNC
#endif       

 end subroutine fscanf_chl
!-------------------------------------------------------------------------------------------------------------
#ifdef INI
    Subroutine fwritef_chlc_inc

    !> **Only if *CHL_C_FUNC* and *CHL_C_BOOL* are true and if the *INI* key is activated.**
    !! Writes in the *calc_ChlC.inc* file, which is dedicated to the
    !! calculation of Chlorophyll to Carbon ratio 
    !! \author Melika Baklouti
    !! \author Vincent Faure
!-------------------------------------------------------------------------------------------------------------
   use mod_eco3m, only :CHL_C_PAR,PROC_MOD,CHL_C_FUNC
   use eco3m_string
   implicit none

   ! local variables 
   Integer             :: idtemp
   integer             :: kk,dim_param
   character(L_VAR_LG)  :: nomsub
   integer             :: chlcinc_id

   chlcinc_id = 1008 
   open(chlcinc_id, file="calc_ChlC.inc", status="replace")

   !  case where a function is used to calculate the Chl:C ratio from the C
   !  concentration
   IF (CHL_C_FUNC) then
      write(chlcinc_id,*)'! This program is automatically created by the Eco3M model'
      write(chlcinc_id,*)'! when both CHL_C_BOOL  and CHL_C_FUNC logical variables are true.'
      write(chlcinc_id,*)'! It contains the expression of the function used for the calculation of the Chl:C ratio'
      write(chlcinc_id,*)'! when Chl is not a state variable'
      write(chlcinc_id,*)'! The name of the Chl:C function is given in the config.ini file (nearly end of file)'

        write(chlcinc_id,*) 
        idtemp = CHL_C_PAR%idproc
        write(nomsub,*) PROC_MOD(idtemp)%nomsub
        nomsub=trim(adjustl(nomsub))
        dim_param = PROC_MOD(idtemp)%nbpar

        write(chlcinc_id,'(A27)')'CHL_C(1)%val = &'
        write(chlcinc_id,*) nomsub(1:len_trim(nomsub)),'(& '
        do kk = 1, dim_param-1
            write(chlcinc_id,*)'     CHL_C_PAR%valpar(',  kk  ,'),  & '
        enddo
        write(chlcinc_id,*)'     CHL_C_PAR%valpar(',  dim_param  ,')  ) '
   !  case where Chl:C ratio is constant or is a pronostic variable 
    ELSE
        write(chlcinc_id,*) "! This file has been automatically generated but is"
        write(chlcinc_id,*) "! not used in this mode (though empty, this file is needed by the"
        write(chlcinc_id,*) "! program for technical reasons"
    ENDIF
    close( chlcinc_id)
    end subroutine fwritef_chlc_inc
#endif
!-------------------------------------------------------------------------------------------------------------

