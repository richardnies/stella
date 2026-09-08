!###############################################################################
!############### X-DEPENDENT KROOK SPONGE FOR TERTIARY MODES ###################
!###############################################################################
! Namelist: &dissipation
!
! Purpose
! -------
! The tertiary-mode theory comes in two flavours.  Either the zonal profile is
! PERIODIC in x, in which case the eigenvalue problem is a Bloch band problem on
! the lattice of zonal extrema and the periodic stella box is already the right
! domain; or the mode is localised about a SINGLE zonal-flow extremum and decays
! as |x| grows.  A periodic box cannot represent the second case on its own: the
! box always contains at least one maximum and one minimum of u_Z, and the mode
! tunnels between neighbouring (equivalent) extrema.
!
! This module supplies the missing ingredient for the localised case: a Krook
! drag
!
!     dg/dt = ... - nu(x) g ,
!
! with nu(x) = 0 in an interior window centred on the extremum of interest and
! rising smoothly to <nu_sponge> towards the edges of the radial box.  The
! damped ("sponge" or "buffer") region absorbs the mode before it can reach the
! neighbouring extremum, so what stella converges to is the eigenmode of the
! single extremum sitting in the middle of the box, decaying at large |x|.
!
! The seven input parameters are read as part of the &dissipation namelist, in
! dissipation.f90 alongside the collision and hyperdissipation switches; this
! module holds only the operator.
!
! The profile is
!
!     xi(x)  = 2 |x - x_c| / Lx                (0 at the centre, 1 at the edge,
!                                               distance measured periodically)
!     nu(x)  = 0                                        for xi <= 1 - w
!            = nu_sponge * [(xi - (1-w)) / w]**p        for xi >  1 - w
!
! with w = <sponge_width>, p = <sponge_exponent> and x_c set by
! <sponge_centre_frac>.
!
! Implementation notes
! --------------------
! * The operator is applied as  F^-1 nu(x) F  on the UNPADDED radial grid x_d
!   (nakx points, exactly nakx Fourier modes).  On that grid F is a full unitary
!   DFT, so the discrete operator is exactly Hermitian and positive
!   semi-definite for nu(x) >= 0: it can only damp, never drive, and there is no
!   aliasing ambiguity.  This is the same transform pair used by
!   calculations_kxky::multiply_by_rho.
! * By default only the NON-ZONAL modes (ky /= 0) are damped.  In a tertiary run
!   the zonal component is a prescribed, frozen background (see <freeze_zonal>,
!   <triangular_ZF>, <cos_ZF> in &parameters_physics) and must not be touched.
!   Set <sponge_zonal> = .true. to damp it too.
! * nu(x) is written to <run_name>.sponge (two columns: x, nu) at
!   initialisation so that post-processing can overlay the buffer region on the
!   eigenfunction.
!
! Convergence check
! -----------------
! The sponge is an approximation: it adds -i nu(x) to the local eigenvalue
! wherever nu /= 0.  For a properly localised mode the eigenfunction is
! exponentially small there, so omega must be insensitive to <nu_sponge> and to
! <sponge_width>.  Scan both before quoting a frequency.
!###############################################################################

module tertiary_sponge

   implicit none

   public :: init_tertiary_sponge
   public :: finish_tertiary_sponge
   public :: add_tertiary_sponge
   public :: nu_sponge_x

   private

   !> The input parameters live in the &dissipation namelist and are read and
   !> broadcast by dissipation::read_parameters:
   !>
   !>   include_tertiary_sponge  master switch
   !>   nu_sponge                peak Krook rate, in units of v_thermal/a
   !>                            (the units of omega)
   !>   sponge_width             fraction of each half-box over which nu ramps
   !>                            from 0 up to nu_sponge; 0 < sponge_width <= 1.
   !>                            0.5 damps the outer half of each side and
   !>                            leaves the middle half of the box clean
   !>   sponge_exponent          power of the ramp; 2.0 gives a C^1 profile,
   !>                            larger is gentler
   !>   sponge_centre_frac       centre of the UNDAMPED window, as a fraction of
   !>                            the radial box length measured from the first
   !>                            radial grid point x_d(1).  The default 0.5 puts
   !>                            it at the middle of the box, which is where the
   !>                            triangular/cosine zonal profiles initialised by
   !>                            dist_fn::init_gxyz put one of their two extrema
   !>                            (the other sits at x = x_d(1), on the box edge,
   !>                            and is the one the sponge eats).  Flip the sign
   !>                            of <triangular_ZF_g_exb> to swap which of the
   !>                            two extrema is the one being studied
   !>   sponge_zonal             damp the ky = 0 component as well; off by
   !>                            default, since the zonal background of a
   !>                            tertiary run is prescribed and frozen
   !>   write_sponge_profile     write <run_name>.sponge at initialisation

   !> nu(x) on the unpadded radial grid x_d(1:nakx).
   real, dimension(:), allocatable :: nu_sponge_x

   logical :: initialised = .false.

contains

   !======================================================================
   !======================== BUILD THE nu(x) PROFILE =====================
   !======================================================================
   subroutine init_tertiary_sponge

      use mp, only: proc0, mp_abort
      use file_utils, only: open_output_file, close_output_file
      use parameters_kxky_grids, only: nakx
      use parameters_physics, only: full_flux_surface, nonlinear
      use grids_kxky, only: x_d, box
      use dissipation, only: include_tertiary_sponge, nu_sponge, sponge_width
      use dissipation, only: sponge_exponent, sponge_centre_frac, write_sponge_profile

      implicit none

      integer :: ikx, unit
      real :: length_x, dx_local, x_centre, offset, xi, xi_flat

      if (initialised) return
      if (.not. include_tertiary_sponge) return
      initialised = .true.

      !> The sponge needs a real-space radial grid, which only exists for
      !> grid_option = 'box'.
      if (.not. box) call mp_abort('tertiary_sponge requires grid_option = "box". aborting.')
      !> the explicit RHS is still held in real space in y at the point where the
      !> sponge is added, so a ky-space operator cannot be applied there
      if (full_flux_surface) call mp_abort('tertiary_sponge is not implemented for full_flux_surface. aborting.')
      if (nakx < 2) call mp_abort('tertiary_sponge requires nakx > 1. aborting.')
      !> the kx <-> x transform pair used below only has its FFTW plans built when
      !> stella decides it needs transforms at all (see stella::check_transforms);
      !> nonlinear = .true. is required for a tertiary run anyway, since the zonal
      !> background can only reach the mode through the ExB term
      if (.not. nonlinear) call mp_abort('tertiary_sponge requires nonlinear = .true. aborting.')
      if (sponge_width <= 0.0 .or. sponge_width > 1.0) &
         call mp_abort('tertiary_sponge requires 0 < sponge_width <= 1. aborting.')

      dx_local = x_d(2) - x_d(1)
      length_x = real(nakx) * dx_local

      x_centre = x_d(1) + sponge_centre_frac * length_x

      !> fraction of the half-box that stays undamped
      xi_flat = 1.0 - sponge_width

      if (.not. allocated(nu_sponge_x)) allocate (nu_sponge_x(nakx))
      nu_sponge_x = 0.0

      do ikx = 1, nakx
         !> signed distance from the centre of the undamped window, wrapped into
         !> [-Lx/2, Lx/2) so that nu(x) respects the periodicity of the box
         offset = x_d(ikx) - x_centre
         offset = offset - length_x * nint(offset / length_x)
         !> normalised distance: 0 at the centre, 1 at the box edge
         xi = 2.0 * abs(offset) / length_x
         if (xi > xi_flat) then
            nu_sponge_x(ikx) = nu_sponge * ((xi - xi_flat) / (1.0 - xi_flat))**sponge_exponent
         end if
      end do

      if (proc0 .and. write_sponge_profile) then
         call open_output_file(unit, '.sponge')
         write (unit, '(a)') '# x   nu_sponge(x)'
         do ikx = 1, nakx
            write (unit, '(2es16.8)') x_d(ikx), nu_sponge_x(ikx)
         end do
         call close_output_file(unit)
      end if

   end subroutine init_tertiary_sponge

   !======================================================================
   !=================== ADD -nu(x) g TO THE RHS OF THE GKE ===============
   !======================================================================
   !> On entry <gke_rhs> already holds code_dt times the RHS assembled so far,
   !> exactly as for sources::add_krook_operator, so the contribution added here
   !> carries an explicit factor of code_dt.
   subroutine add_tertiary_sponge(g, gke_rhs)

      use zgrid, only: nzgrid, ntubes
      use parameters_kxky_grids, only: naky, nakx
      use grids_kxky, only: zonal_mode
      use stella_layouts, only: vmu_lo
      use stella_time, only: code_dt
      use stella_transforms, only: transform_kx2x_unpadded, transform_x2kx_unpadded
      use dissipation, only: include_tertiary_sponge, sponge_zonal

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gke_rhs

      complex, dimension(:, :), allocatable :: g0k, g0x
      integer :: iz, it, ivmu, iky

      if (.not. include_tertiary_sponge) return
      if (.not. allocated(nu_sponge_x)) return

      allocate (g0k(naky, nakx))
      allocate (g0x(naky, nakx))

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               g0k = g(:, :, iz, it, ivmu)
               !> leave the (frozen, prescribed) zonal background alone unless
               !> the user asks for it to be damped as well
               if (zonal_mode(1) .and. .not. sponge_zonal) g0k(1, :) = 0.0
               call transform_kx2x_unpadded(g0k, g0x)
               do iky = 1, naky
                  g0x(iky, :) = nu_sponge_x * g0x(iky, :)
               end do
               call transform_x2kx_unpadded(g0x, g0k)
               gke_rhs(:, :, iz, it, ivmu) = gke_rhs(:, :, iz, it, ivmu) - code_dt * g0k
            end do
         end do
      end do

      deallocate (g0k, g0x)

   end subroutine add_tertiary_sponge

   !======================================================================
   subroutine finish_tertiary_sponge

      implicit none

      if (allocated(nu_sponge_x)) deallocate (nu_sponge_x)
      initialised = .false.

   end subroutine finish_tertiary_sponge

end module tertiary_sponge
