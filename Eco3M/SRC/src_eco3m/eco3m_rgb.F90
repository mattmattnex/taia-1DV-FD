
!> Subroutines for the implementation of the RGB light extinction RBG scheme
!! \date 2015-12
!! \author Nicolas Barrier


#ifdef RGB
!-------------------------------------------------------------------------------------------------    
    subroutine eco3m_rgb_init

    !> Initialisation of the RGB module.
    !! 
    !! - Reading the number of wavelengths in the namelist (*rgb_nbands* variable)
    !! - Reads the number of lines in the kRGB.txt file
    !! - Allocates the *rgb_coef_table* 2D-array
    !! - Filling the values of the *rgb_coef_table* by reading the file 
    !! - Allocation of the *E_PARZ_RGB* and *kectinc_rgb* arrays
    !!
    !! \author Nicolas Barrier, Melika Baklouti
!------------------------------------------------------------------------------------------------
    use mod_eco3m_files
    use mod_eco3m
    implicit none
    integer :: rstatus !< reading status
    integer :: iline, iband  !< line and column indices
    real(8) :: dummyvar      !< Dummy variable that is used for the first column
                             ! of the kRGB file !<
    integer :: file_rgb_id      !< ID of the KRGB file
    logical :: file_exist
    integer :: istat

    namelist/nam_eco3m_rgb/filename_krgb,rgb_nbands 

    REWIND(namelist_eco3m_id)
    READ  (namelist_eco3m_id, nam_eco3m_rgb)
#if defined CALC
    write(file_CR_id,*) "=================Eco3M: Initialisation of the RGB module================" 
    write(file_CR_id,*) "Extinction coefficient file: ", trim(filename_krgb)
    write(file_CR_id,*) "Number of wavelengths: ", rgb_nbands
#endif
! Opening of the KRGB file
    filename_krgb = trim(adjustl(eco3m_root_dir))//trim(adjustl(filename_krgb))
    inquire(file=filename_krgb, exist=file_exist)
    if (.not.(file_exist)) then
        write(*,*) "****************************************"
        write(*,*) "The ", trim(filename_krgb), " file does not exist"
        write(*,*) "This program will stop"
        write(*,*) "****************************************"
        stop
    end if
        
    file_rgb_id = 1012   
    open(file_rgb_id, file=filename_krgb)

! initialisation of the number of lines
    rgb_ncoeff = 0
! reading status is set to 0. We loop until rstatus is not 0 anymore,
! i.e until the end of file is reached
    rstatus = 0
    do while(rstatus.eq.0)
        read(file_rgb_id, *, iostat=rstatus)
        rgb_ncoeff = rgb_ncoeff + 1
    end do 
! Correction of line number and allocation of the table
    rgb_ncoeff = rgb_ncoeff -1 
    
    Allocate (rgb_coef_table(1:rgb_nbands,1:rgb_ncoeff),stat = istat)
    if (istat /=0) write(*,*) "pb with the allocation of rgb_coef_table"
! reading the file back at the start, and then filling the array
    rewind(file_rgb_id) 
    do iline=1, rgb_ncoeff
       read(file_rgb_id, *) dummyvar, (rgb_coef_table(iband, iline), iband=1, rgb_nbands)
    end do 
     
! Closing of the krgb file
  close(file_rgb_id)
#ifdef CALC
  write(file_CR_id,*) "Number of lines in the file: ", rgb_ncoeff
  write(file_CR_id,*) 
#endif
  Allocate(e_parz_rgb(1:rgb_nbands,nx_min:nx_max,ny_min:ny_max,1:nz_max),stat=istat) 
  if (istat /=0) write(*,*) "pb with the allocation of e_parz_rgb"

  Allocate(kextinc_rgb(1:rgb_nbands,nx_min:nx_max,ny_min:ny_max,1:nz_max),stat=istat)
  if (istat /=0) write(*,*) "pb with the allocation of kextinc_rgb"
  kextinc_rgb(:,:,:,:) = 0.0d0
    end subroutine eco3m_rgb_init
!----------------------------------------------------------------------------------------------
    subroutine compute_rgb_kextinc(Chl_tot)

    !> Computes the light extinction coefficients by using the
    !! value of total Chl (*Chl_tot* variable). It looks for the line index, in the 
    !! *rgb_coef_table* array,  that corresponds to the current
    !! concentration. It then extracts the coefficients. for each wavelength (i.e column), and
    !! stores them in the *kextinc_rgb* variable.
    !! \author Nicolas Barrier
    !! \note The line index is corrected so that \f$1 \leq index \leq nb_{line}\f$
!----------------------------------------------------------------------------------------------
    use mod_eco3m
    implicit none
    integer :: i, j, k, irgb 
    integer :: tab_index
    Real(8):: Chl_tot(nx_min:nx_max,ny_min:ny_max,1:nz_max)

        do k = 1, nz_max-1
            do j = ny_min, ny_max
                do i = nx_min, nx_max
                    ! We look for the tab index that corresponds to the
                    ! Chltot concentration (expressed in µmol/l) 
                    !tab_index = NINT( 41 + 20.* LOG10(max(Chl_tot(i, j, k)*1e6,1d-10)))
                    tab_index = NINT( 41 + 20.* LOG10(max(Chl_tot(i, j, k)*1e6,1d-10)))
                    ! Correction of the index so that it remains within the file
                    ! limits 
                    tab_index = min(tab_index, rgb_ncoeff) 
                    tab_index = max(1, tab_index)

                    ! Looping over the different wavelength: picking up the
                    ! coefficient within the table
                    do irgb = 1, rgb_nbands   
                        kextinc_rgb(irgb, i, j, k) = rgb_coef_table(irgb, tab_index)
                    end do

                end do
            end do
        end do

    end subroutine compute_rgb_kextinc

!----------------------------------------------------------------------------------------------
    subroutine compute_eparz_rgb(Chl_tot)


    !> Computes the Photosynthetically Available Radiation for each wavelengths. 
    !!  
    !! - Calculation of the light extinction coefficients for each wavelength (compute_rgb_kextinc() subroutine)
    !! by using total Chlorophyl values (*Chl_tot* variable)
    !! - Calculation of the EPAR for each wavelength 
    !! - Calculation of the PAR that will be used as the mean computed over each wavelengths.
    !! 
    !! This subroutine is called by the eco3m_irrad::eco3m_compute_eparz() 
    !! (see the latter documentation for a detailed description of EPAR calculation)
    !! \author Melika Baklouti, Nicolas Barrier
    !! \note See the eco3m_irrad::eco3m_compute_eparz() documentation
    !! for a detailed description of the EPAR calculation.
!----------------------------------------------------------------------------------------------
    use mod_eco3m
    implicit none
    integer :: i, j, k, irgb  ! Loop indexes
    Real(8):: Chl_tot(nx_min:nx_max,ny_min:ny_max,1:nz_max)

        if (pos_nzmax == 'SURF') then 
        ! Calculation of EPARZ if nz_max is at the surface
            do k = nz_max, 1, -1
                do j = ny_min, ny_max
                    do i = nx_min, nx_max

                        ! We also loop over the wavelengths to compute the
                        ! E_PARZ for each wavelength, given the coefficient
                        ! provided in the RGB module. Then, we compute the
                        ! average of the EPAR over each wavelength
                        do irgb = 1, rgb_nbands
                            if (k == nz_max) then
                                E_PARZ_RGB(irgb,i,j,k)=E_PAR(i,j)*exp(-0.5*kextinc_rgb(irgb,i,j,k)*dz(i,j,k))
                            else
                                E_PARZ_RGB(irgb,i,j,k) = E_PARZ_RGB(irgb,i,j,k+1)*exp(-0.5*(kextinc_rgb(irgb,i,j,k+1)*dz(i,j,k+1)+&
                                    kextinc_rgb(irgb,i,j,k)*dz(i,j,k)))
                            endif
                        end do
                    end do
                end do
            end do

        elseif (pos_nzmax == 'BOTT') then
            ! Calculation of EPARZ if nz_max is at the bottom
            do k = 1, nz_max
                do j = ny_min, ny_max
                    do i = nx_min, nx_max

                        ! We also loop over the wavelengths to compute the
                        ! E_PARZ for each wavelength, given the coefficient
                        ! provided in the RGB module. Then, we compute the
                        ! average of the EPAR over each wavelength
                        do irgb = 1, rgb_nbands
                            if (k == 1) then
                                E_PARZ_RGB(irgb,i,j,k)=E_PAR(i,j)*exp(-0.5*kextinc_rgb(irgb,i,j,k)*dz(i,j,k))
                            else
                                E_PARZ_RGB(irgb,i,j,k) = E_PARZ_RGB(irgb,i,j,k-1)*exp(-0.5*(kextinc_rgb(irgb,i,j,k-1)*dz(i,j,k-1)+&
                                    kextinc_rgb(irgb,i,j,k)*dz(i,j,k)))
                            endif
                        end do
                    end do 
                end do
            end do
        end if

        ! Now we compute the E_PARZ value as the mean computed over all the
        ! wavelengths
        E_PARZ = sum(E_PARZ_RGB, dim=1)/(rgb_nbands)

        end subroutine compute_eparz_rgb
#endif
