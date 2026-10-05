!=======================================================================
! eco3m_capi.F90
!
! C-callable interface to Eco3M, for driving the biogeochemistry from an
! external model (here a Python 1DV physics model).
!
! Compiled only in the ECO3M_SUB + COUPL + CALC configuration
! (see Makefile.conf_PY). In that configuration:
!
!   - eco3m_main is a SUBROUTINE with no time loop. Its whole body is
!         eco3m_compute_eparz -> [compute_chl_c_ratio] -> eco3m_compute_sms
!   - TEND holds a RATE in per-second units (not an increment)
!   - Eco3M reads no initial conditions and owns no vertical grid;
!     the caller supplies both
!
! Existence of this file is why the caller cannot talk to Eco3M directly:
! VAR(:)%conc is an array of derived types with allocatable components,
! whose memory layout is a compiler-internal array descriptor.
!
! CALL ORDER (the first two must not be swapped):
!
!   eco3m_c_set_dims(...)      <-- sizes every array allocated below
!   eco3m_c_initialize()       <-- reads namelist + config.ini, allocates
!   eco3m_c_get_nbvar(...) / eco3m_c_get_varname(...)
!   eco3m_c_set_dz(...)
!   eco3m_c_set_conc(...)
!   loop:
!       eco3m_c_set_physics(...)
!       eco3m_c_set_conc(...)
!       eco3m_c_step()
!       eco3m_c_get_tend(...)
!
! Array convention at the boundary: (nz_max, nbvar), Fortran-ordered
! (first index varies fastest), matching TEND's native layout.
!
! 1DV only: assumes nx_min = nx_max = ny_min = ny_max = 1.
!=======================================================================
module eco3m_capi

    use iso_c_binding
    use mod_eco3m
    use mod_eco3m_irrad
    implicit none

contains

!-----------------------------------------------------------------------
! Internal helper: Fortran string -> null-terminated C buffer
!-----------------------------------------------------------------------
    subroutine put_cstring(src, dst, nbuf)
        integer,                intent(in)  :: nbuf
        character(len=*),       intent(in)  :: src
        character(kind=c_char), intent(out) :: dst(nbuf)
        integer :: i, n
        character(len=len_trim(adjustl(src))) :: t

        t = trim(adjustl(src))
        n = min(len(t), nbuf - 1)
        do i = 1, n
            dst(i) = t(i:i)
        end do
        dst(n+1) = c_null_char
    end subroutine put_cstring

!=======================================================================
! 1. GRID DIMENSIONS  --  must be called before eco3m_c_initialize
!=======================================================================
! Under ECO3M_SUB the namelist group nam_eco3m_run_main is NOT read
! (see eco3m_runpar_init), so nx_min..nz_max arrive only from here.
! eco3m_init_config allocates VAR%conc, TEND, TEMP_BIO, SAL_BIO, dz,
! irrad, E_PARZ and kextinc_rgb using these values.
!-----------------------------------------------------------------------
    subroutine eco3m_c_set_dims(nxmin, nxmax, nymin, nymax, nzmax) &
               bind(C, name="eco3m_c_set_dims")
        integer(c_int), value :: nxmin, nxmax, nymin, nymax, nzmax

        nx_min = nxmin
        nx_max = nxmax
        ny_min = nymin
        ny_max = nymax
        nz_max = nzmax
    end subroutine eco3m_c_set_dims

!=======================================================================
! 2. INITIALISATION
!=======================================================================
    subroutine eco3m_c_initialize() bind(C, name="eco3m_c_initialize")
        call eco3m_init_config
        tps = 0.d0
    end subroutine eco3m_c_initialize

!=======================================================================
! 3. MODEL DESCRIPTION  --  build the index map at runtime
!=======================================================================
    subroutine eco3m_c_get_nbvar(n) bind(C, name="eco3m_c_get_nbvar")
        integer(c_int), intent(out) :: n
        n = nbvar
    end subroutine eco3m_c_get_nbvar

    subroutine eco3m_c_get_nzmax(n) bind(C, name="eco3m_c_get_nzmax")
        integer(c_int), intent(out) :: n
        n = nz_max
    end subroutine eco3m_c_get_nzmax

    subroutine eco3m_c_get_dt_bio(v) bind(C, name="eco3m_c_get_dt_bio")
        real(c_double), intent(out) :: v
        v = dt_bio
    end subroutine eco3m_c_get_dt_bio

!-- "SURF" or "BOTT": which end of the k index is the surface
    subroutine eco3m_c_get_pos_nzmax(buf, nbuf) bind(C, name="eco3m_c_get_pos_nzmax")
        integer(c_int), value                :: nbuf
        character(kind=c_char), intent(out)  :: buf(nbuf)
        call put_cstring(pos_nzmax, buf, nbuf)
    end subroutine eco3m_c_get_pos_nzmax

!-- Name triplet of state variable ivar (1-based), e.g. "phy"/"PHYS"/"N"
    subroutine eco3m_c_get_varname(ivar, comp, scomp, elmt, nbuf) &
               bind(C, name="eco3m_c_get_varname")
        integer(c_int), value               :: ivar, nbuf
        character(kind=c_char), intent(out) :: comp(nbuf), scomp(nbuf), elmt(nbuf)

        if (ivar < 1 .or. ivar > nbvar) then
            call put_cstring("", comp,  nbuf)
            call put_cstring("", scomp, nbuf)
            call put_cstring("", elmt,  nbuf)
            return
        end if

        call put_cstring(VAR(ivar)%comp,  comp,  nbuf)
        call put_cstring(VAR(ivar)%scomp, scomp, nbuf)
        call put_cstring(VAR(ivar)%elmt,  elmt,  nbuf)
    end subroutine eco3m_c_get_varname

!=======================================================================
! 4. VERTICAL GRID
!=======================================================================
! dz MUST be positive. compute_eparz_rgb uses dz without abs(), unlike
! the non-RGB path, so negative thicknesses make light grow with depth.
!-----------------------------------------------------------------------
    subroutine eco3m_c_set_dz(arr) bind(C, name="eco3m_c_set_dz")
        real(c_double), intent(in) :: arr(nz_max)
        integer :: k
        do k = 1, nz_max
            dz(nx_min:nx_max, ny_min:ny_max, k) = arr(k)
        end do
    end subroutine eco3m_c_set_dz

!-- Not read by eco3m_compute_eparz, but allocated and possibly used by
!-- user-supplied process functions. Set it for safety.
    subroutine eco3m_c_set_depth(arr) bind(C, name="eco3m_c_set_depth")
        real(c_double), intent(in) :: arr(nz_max)
        integer :: k
        do k = 1, nz_max
            depth_z(nx_min:nx_max, ny_min:ny_max, k) = arr(k)
        end do
    end subroutine eco3m_c_set_depth

    subroutine eco3m_c_set_bathymetry(v) bind(C, name="eco3m_c_set_bathymetry")
        real(c_double), value :: v
        bathymetry(nx_min:nx_max, ny_min:ny_max) = v
    end subroutine eco3m_c_set_bathymetry

!=======================================================================
! 5. PER-STEP INPUTS
!=======================================================================
! irrad_surf must be NET (albedo-corrected) surface shortwave in W/m2.
! Under COUPL, eco3m_compute_eparz does not apply albedo -- it assumes
! the physical model already did. The namelist "albedo" entry is dead.
!-----------------------------------------------------------------------
    subroutine eco3m_c_set_physics(temp, sal, irrad_surf, time_s) &
               bind(C, name="eco3m_c_set_physics")
        real(c_double), intent(in) :: temp(nz_max), sal(nz_max)
        real(c_double), value      :: irrad_surf, time_s
        integer :: k

        do k = 1, nz_max
            TEMP_BIO(nx_min:nx_max, ny_min:ny_max, k) = temp(k)
            SAL_BIO (nx_min:nx_max, ny_min:ny_max, k) = sal(k)
        end do

        irrad(nx_min:nx_max, ny_min:ny_max) = irrad_surf
        tps = time_s
    end subroutine eco3m_c_set_physics

!-- conc(k, ivar), Fortran-ordered. Needs a loop: VAR is an array of
!-- derived types, so concentrations are not one contiguous block.
    subroutine eco3m_c_set_conc(conc) bind(C, name="eco3m_c_set_conc")
        real(c_double), intent(in) :: conc(nz_max, nbvar)
        integer :: ivar, k

        do ivar = 1, nbvar
            do k = 1, nz_max
                VAR(ivar)%conc(nx_min:nx_max, ny_min:ny_max, k) = conc(k, ivar)
            end do
        end do
    end subroutine eco3m_c_set_conc

!=======================================================================
! 6. THE STEP
!=======================================================================
! eco3m_main under ECO3M_SUB:
!     eco3m_compute_eparz -> [compute_chl_c_ratio] -> eco3m_compute_sms
!-----------------------------------------------------------------------
    subroutine eco3m_c_step() bind(C, name="eco3m_c_step")
        call eco3m_main
    end subroutine eco3m_c_step

!=======================================================================
! 7. OUTPUTS
!=======================================================================
! TEND is a RATE in per-second units.
!-----------------------------------------------------------------------
    subroutine eco3m_c_get_tend(tend_out) bind(C, name="eco3m_c_get_tend")
        real(c_double), intent(out) :: tend_out(nz_max, nbvar)
        integer :: ivar, k

        do ivar = 1, nbvar
            do k = 1, nz_max
                tend_out(k, ivar) = TEND(nx_min, ny_min, k, ivar)
            end do
        end do
    end subroutine eco3m_c_get_tend

    subroutine eco3m_c_get_conc(conc_out) bind(C, name="eco3m_c_get_conc")
        real(c_double), intent(out) :: conc_out(nz_max, nbvar)
        integer :: ivar, k

        do ivar = 1, nbvar
            do k = 1, nz_max
                conc_out(k, ivar) = VAR(ivar)%conc(nx_min, ny_min, k)
            end do
        end do
    end subroutine eco3m_c_get_conc

!-- PAR profile.
    subroutine eco3m_c_get_eparz(arr) bind(C, name="eco3m_c_get_eparz")
        real(c_double), intent(out) :: arr(nz_max)
        integer :: k
        do k = 1, nz_max
            arr(k) = E_PARZ(nx_min, ny_min, k)
        end do
    end subroutine eco3m_c_get_eparz

!-- Surface PAR just below the interface
    subroutine eco3m_c_get_epar(v) bind(C, name="eco3m_c_get_epar")
        real(c_double), intent(out) :: v
        v = E_PAR(nx_min, ny_min)
    end subroutine eco3m_c_get_epar

!-- A single element of the flux matrix, for debugging individual
    subroutine eco3m_c_get_flux(ivar, jvar, arr) bind(C, name="eco3m_c_get_flux")
        integer(c_int), value       :: ivar, jvar
        real(c_double), intent(out) :: arr(nz_max)
        integer :: k

        arr = 0.d0
        if (ivar < 1 .or. ivar > nbvar) return
        if (jvar < 1 .or. jvar > nbvar) return
        if (.not. allocated(FLUX_VAL(ivar, jvar)%idproc)) return

        do k = 1, nz_max
            arr(k) = FLUX_VAL(ivar, jvar)%val(nx_min, ny_min, k)
        end do
    end subroutine eco3m_c_get_flux

!=======================================================================
! 8. PARAMETERS
!=======================================================================
! Values come from config.ini at init, but nothing stops overwriting
! them afterwards. This is the hook for calibration or parameter
! estimation without re-running text2numb and the INI stage.
!-----------------------------------------------------------------------
    subroutine eco3m_c_get_nflux_par(n) bind(C, name="eco3m_c_get_nflux_par")
        integer(c_int), intent(out) :: n
        n = 0
        if (allocated(FLUX_PAR)) n = size(FLUX_PAR)
    end subroutine eco3m_c_get_nflux_par

    subroutine eco3m_c_get_nparam(iproc, n) bind(C, name="eco3m_c_get_nparam")
        integer(c_int), value       :: iproc
        integer(c_int), intent(out) :: n

        n = 0
        if (.not. allocated(FLUX_PAR)) return
        if (iproc < 1 .or. iproc > size(FLUX_PAR)) return
        if (.not. allocated(FLUX_PAR(iproc)%valpar)) return
        n = size(FLUX_PAR(iproc)%valpar)
    end subroutine eco3m_c_get_nparam

    subroutine eco3m_c_get_param(iproc, ipar, v, ierr) bind(C, name="eco3m_c_get_param")
        integer(c_int), value       :: iproc, ipar
        real(c_double), intent(out) :: v
        integer(c_int), intent(out) :: ierr

        v = 0.d0
        ierr = 1
        if (.not. allocated(FLUX_PAR)) return
        if (iproc < 1 .or. iproc > size(FLUX_PAR)) return
        if (.not. allocated(FLUX_PAR(iproc)%valpar)) return
        if (ipar < 1 .or. ipar > size(FLUX_PAR(iproc)%valpar)) return

        v = FLUX_PAR(iproc)%valpar(ipar)
        ierr = 0
    end subroutine eco3m_c_get_param

    subroutine eco3m_c_set_param(iproc, ipar, v, ierr) bind(C, name="eco3m_c_set_param")
        integer(c_int), value       :: iproc, ipar
        real(c_double), value       :: v
        integer(c_int), intent(out) :: ierr

        ierr = 1
        if (.not. allocated(FLUX_PAR)) return
        if (iproc < 1 .or. iproc > size(FLUX_PAR)) return
        if (.not. allocated(FLUX_PAR(iproc)%valpar)) return
        if (ipar < 1 .or. ipar > size(FLUX_PAR(iproc)%valpar)) return

        FLUX_PAR(iproc)%valpar(ipar) = v
        ierr = 0
    end subroutine eco3m_c_set_param

    subroutine eco3m_c_set_dt_bio(v) bind(C, name="eco3m_c_set_dt_bio")
        real(c_double), value :: v
        dt_bio = v
    end subroutine eco3m_c_set_dt_bio

end module eco3m_capi