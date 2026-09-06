!###############################################################################
!##################### ROSENBLUTH-HINTON RESPONSE FUNCTIONS ####################
!###############################################################################
!
! Calculation of the "Rosenbluth-Hinton" (RH) response functions that govern the
! evolution of long-wavelength stationary zonal flows:
!
!    RH_inertia         the RH polarisation/inertia,   I(kx,z,tube,s)
!    RH_integrand_even  the transit-averaged response, even in vpa
!    RH_integrand_odd   the transit-averaged response, odd  in vpa
!    RH_U_parallel_fac  the kx -> 0 limit used by the parallel-shear term
!
! together with the routines that evaluate the RH fluxes and the RH potential
! from a given distribution function.
!
! These are NOT diagnostics.  <RH_integrand_even/odd> and <RH_inertia> are read
! by dist_fn::init_gxyz when building a prescribed zonal profile, and
! <RH_U_parallel_fac> feeds prl_shear in flow_shear.  They must therefore be
! initialised before ginit, independently of whether any RH diagnostic was
! requested.  The diagnostics module <diagnostics_RH_inertia_fluxes> writes
! these quantities out; it does not compute them.
!
!###############################################################################

module rosenbluth_hinton

   implicit none

   public :: init_rosenbluth_hinton
   public :: finish_rosenbluth_hinton
   public :: get_RH_inertia_fluxtube
   public :: get_RH_fluxes_fluxtube
   public :: get_RH_phi_I_fluxtube
   public :: eval_transit_ints
   public :: eval_transit_int_integrand_RH
   public :: eval_Q_fac
   public :: eval_Q_profile_hat
   public :: RH_U_parallel_fac
   public :: RH_inertia
   public :: RH_integrand_even, RH_integrand_odd
   public :: RH_LW_even, RH_LW_odd
   public :: RH_SW_even, RH_SW_odd
   public :: RH_pmom_weight, RH_pmom_inertia
   public :: get_RH_pmom, get_RH_pmom_fluxes_fluxtube
   public :: RH_drift_bounce_avg
   public :: RH_drift_is_trapped
   public :: eval_bounce_averaged_drift

   real, dimension(:,:), allocatable :: RH_U_parallel_fac
   ! (-nzgrid:nzgrid, -vmu-layout-)
   ! No tube index: the construction below depends only on (iz, ivmu).

   complex, dimension(:,:,:,:), allocatable :: RH_inertia
   ! (nakx, -nzgrid:nzgrid, ntubes, nspec)

   complex, dimension(:,:,:,:), allocatable :: RH_integrand_even, RH_integrand_odd

   !> Long-wavelength approximations to the same weights, filled only when
   !> <write_RH_asymptotics> asks for them.
   complex, dimension(:,:,:,:), allocatable :: RH_LW_even, RH_LW_odd
   complex, dimension(:,:,:,:), allocatable :: RH_SW_even, RH_SW_odd

   !> The parallel-flow counterpart of the Rosenbluth-Hinton projection.
   !> <RH_pmom_weight> is V_sigma(z) = <vpa J0 exp(-Q)>_tau exp(Q(z)), the
   !> sigma-odd member of the same family of weights the linear streaming and
   !> radial drift annihilate; <RH_pmom_inertia> is the same projection applied
   !> to a unit-flow shifted Maxwellian, so that the ratio of their field-line
   !> averages is the parallel flow.  See DOCUMENTATION/RH_parallel_flow.
   complex, dimension(:,:,:,:), allocatable :: RH_pmom_weight
   complex, dimension(:,:,:,:), allocatable :: RH_pmom_inertia
   ! (nakx, -nzgrid:nzgrid, ntubes, -vmu-layout-)

   !> Transit-averaged radial magnetic drift: the bounce average over its own
   !> well for a trapped particle, the average along the whole field line for a
   !> passing one.  Zero only where no complete orbit could be identified, which
   !> is a well running off the end of the field line.  Needs no drift-orbit
   !> phase, so it is available in geometries where the rest of the
   !> Rosenbluth-Hinton machinery is not.
   real, dimension(:,:), allocatable :: RH_drift_bounce_avg
   ! (-nzgrid:nzgrid, -vmu-layout-)

   !> Which of those orbits are trapped.  The drift drive is reported separately
   !> for the two populations because they do not stand on the same footing.  A
   !> trapped particle's average is over its own well, which is wholly inside the
   !> simulated tube, so it is the orbit average whatever the tube.  A passing
   !> particle's is taken over the tube, and on an irrational surface the field
   !> line never closes: the true average is over the flux surface, and what the
   !> tube gives instead is an artefact of the flux-tube construction.  On a
   !> rational surface, with the tube spanning the closed line, the two coincide
   !> and the passing contribution is physical.
   logical, dimension(:,:), allocatable :: RH_drift_is_trapped
   ! (-nzgrid:nzgrid, -vmu-layout-)

   private

   ! Debugging
   logical :: debug = .false.

   !> Surrogate for kx -> 0 used when evaluating RH_U_parallel_fac.  The
   !> response is expanded about kx = 0 and divided by kx, so kxsmall must be
   !> small enough that the O(kxsmall^2) error is negligible, yet large enough
   !> that cancellation in (1 - <...>)/kxsmall stays well inside double
   !> precision.  1e-8 sits near the sqrt(epsilon) sweet spot for both.
   real, parameter :: kxsmall = 1.e-8

   ! Has this module been initialised?
   logical :: rosenbluth_hinton_initialized = .false.

   !> Whether Q comes from the closed form or from integration along the field
   !> line.  Resolved once in init_rosenbluth_hinton, because the default depends
   !> on what the active geometry can supply.
   logical :: use_analytic_drift_phase = .true.

contains


!###############################################################################
!############################ INITALIZE & FINALIZE #############################
!###############################################################################

   !============================================================================
   !======================== INITALIZE THE DIAGNOSTICS =========================
   !============================================================================
   subroutine init_rosenbluth_hinton()

      use mp, only: proc0, mp_abort
      use geometry, only: RH_drift_phase_defined
      use parameters_physics, only: RH_analytic_drift_phase, RH_analytic_drift_phase_specified
      use parameters_physics, only: full_flux_surface, radial_variation
      use parameters_diagnostics, only: write_RH_bounce_drift

      ! Dimensions
      use parameters_kxky_grids, only: nakx
      use grids_kxky, only: akx
      use zgrid, only: nzgrid, ntubes, nztot
      use species, only: nspec
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use vpamu_grids, only: vpa, vperp2, mu
      use geometry, only: bmag
      use species, only: spec
      use constants, only: zi

      use geometry, only: RH_drift_phase_fac
      use parameters_diagnostics, only: write_RH_asymptotics
      implicit none

      real :: energyval, muval, bmag_max, drift_average
      real, dimension(:), allocatable :: Q_hat_z
      complex :: integrand_tmp_pls, integrand_tmp_min, integrand_tmp_v
      complex :: LW_even_tmp, LW_odd_tmp
      complex :: SW_even_tmp, SW_odd_tmp
      logical :: trapped, well_found

      integer :: ivmu, iv, imu, is, ia, iz, it, ikx

      ia = 1

      !> The RH response functions are needed when they are diagnosed, when the
      !> parallel-shear term uses them (<omprimfac_RH>), or when a prescribed
      !> zonal profile is built from the RH closure.  Building them costs an
      !> O(nvmu * nz^2 * nakx) transit-average loop, so skip it otherwise.
      if (.not. rosenbluth_hinton_needed()) return

      !> The bounce-averaged radial drift comes first: it needs only the geometry,
      !> it is what the drift-orbit phase has to be generalised with, and it is
      !> the quantity that distinguishes a general stellarator from a tokamak.  So
      !> it is computed even where the rest of the machinery cannot run.
      allocate (RH_drift_bounce_avg(-nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      RH_drift_bounce_avg = 0.
      allocate (RH_drift_is_trapped(-nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      RH_drift_is_trapped = .false.

      !> Filled here only to be reported.  Where Q is integrated along the field
      !> line the main loop below overwrites this with the average that build
      !> actually subtracted, which is the one the drift flux has to use; and
      !> where Q is the closed form the drift flux is zero by construction, so
      !> nothing but the diagnostic wants these numbers.  The loop costs
      !> O(nvmu nz^2), so it is skipped when nobody asked for them.
      if (write_RH_bounce_drift) then
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         do iz = -nzgrid, nzgrid
            !> The transit average, not the bounce average: what survives the
            !> Rosenbluth-Hinton projection is i kx <v_Mx> averaged over the orbit
            !> the particle actually executes, which is its own well if it is
            !> trapped and the whole field line if it is passing.  Passing
            !> particles are most of velocity space and their transit-averaged
            !> radial drift does not vanish in a stellarator, so leaving them at
            !> zero -- as taking only the bounce average over a well does -- drops
            !> the bulk of the drive.
            call eval_drift_transit_average(vpa(iv)**2 + vperp2(ia, iz, imu), mu(imu), iz, &
                                            drift_average, well_found)
            if (well_found) RH_drift_bounce_avg(iz, ivmu) = drift_average
         end do
      end do
      end if

      !> Choose how Q is obtained.  Unless the input file asked for one, follow the
      !> geometry: the closed form where the geometry supplies it, which is
      !> Miller, and integration along the field line where it does not, which is
      !> VMEC.  Asking for the closed form where it does not exist is an error
      !> rather than something to silently substitute.
      if (RH_analytic_drift_phase_specified) then
         use_analytic_drift_phase = RH_analytic_drift_phase
      else
         use_analytic_drift_phase = RH_drift_phase_defined
      end if
      if (use_analytic_drift_phase .and. .not. RH_drift_phase_defined) call mp_abort &
         ('RH_analytic_drift_phase = .true. was requested, but the active geometry does &
          &not provide the closed-form drift-orbit phase.  Set it to .false. to integrate &
          &Q along the field line instead.  Aborting.')
      if (proc0 .and. debug) write (*, *) 'rosenbluth_hinton: analytic drift phase = ', use_analytic_drift_phase

      if (full_flux_surface) call mp_abort &
         ('Rosenbluth-Hinton diagnostics are not implemented for full_flux_surface.  Aborting.')
      if (radial_variation) call mp_abort &
         ('Rosenbluth-Hinton diagnostics are not implemented for radial_variation.  Aborting.')

      ! Only initialize once
      if (rosenbluth_hinton_initialized) return
      rosenbluth_hinton_initialized = .true.


      ! Only debug on the first processor
      debug = debug .and. proc0

      ! Allocate the arrays for the Rosenbluth-Hinton integrand term
      allocate (RH_integrand_even(nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); RH_integrand_even = 0.
      allocate (RH_integrand_odd( nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); RH_integrand_odd  = 0.
      allocate (RH_pmom_weight(  nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); RH_pmom_weight = 0.
      allocate (RH_pmom_inertia( nakx, -nzgrid:nzgrid, ntubes, nspec)); RH_pmom_inertia = 0.
      if (write_RH_asymptotics) then
         allocate (RH_LW_even(nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); RH_LW_even = 0.
         allocate (RH_LW_odd( nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); RH_LW_odd  = 0.
         allocate (RH_SW_even(nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); RH_SW_even = 0.
         allocate (RH_SW_odd( nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); RH_SW_odd  = 0.
      end if

      ! Allocate array for RH_U_parallel_fac
      allocate (RH_U_parallel_fac( -nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); RH_U_parallel_fac = 0.

      ! Allocate the array for the RH_inertia
      allocate (RH_inertia(nakx, -nzgrid:nzgrid, ntubes, nspec)); RH_inertia = 0

      !> Trapped/circulating label for the two drift-flux output channels.  The
      !> transit averages themselves do their own per-well decomposition (see
      !> find_well), so this global-maximum test is not a physics gate: both
      !> branches at the bottom of get_RH_bounce_drift_fluxes_fluxtube evaluate
      !> the same expression and only the channel differs.  In a stellarator a
      !> particle trapped in a shallow local well is therefore counted in the
      !> circulating channel, which shifts the split but not the sum that enters
      !> the budget.
      bmag_max = maxval(bmag(ia,:))

      allocate (Q_hat_z(-nzgrid:nzgrid)); Q_hat_z = 0.

      ! Evaluate the transit averages
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)

         do iz = -nzgrid, nzgrid

            ! Determine energy and magnetic moment
            energyval = vpa(iv)**2 + vperp2(ia,iz,imu)
            muval     = mu(imu)

            ! Is this particle trapped in the well?
            trapped = energyval <= 2*muval*bmag_max

            !> The drift-orbit phase is linear in kx and, at fixed pitch angle,
            !> scales as sqrt(energy); the species enters only as a scalar and
            !> the sign of v_par only as an overall sign.  What is left depends
            !> on the field line and on lambda = mu / energy alone, so it is
            !> built once here rather than nakx + 1 times inside the loops below.
            if (.not. use_analytic_drift_phase) then
               if (energyval > epsilon(0.)) then
                  call eval_Q_profile_hat(muval / energyval, iz, Q_hat_z, drift_average)
                  !> eval_Q_profile_hat works at unit energy and the drift is
                  !> linear in energy, so restore it here.  Taking the average
                  !> from the Q build rather than recomputing it is what keeps the
                  !> drift flux consistent with the phase, and so what lets the
                  !> budget close.
                  RH_drift_bounce_avg(iz, ivmu) = energyval * drift_average
                  RH_drift_is_trapped(iz, ivmu) = trapped
               else
                  Q_hat_z = 0.
                  RH_drift_bounce_avg(iz, ivmu) = 0.
               end if
            end if

            do it = 1, ntubes

               do ikx = 1, nakx

                  call get_RH_transit_integrands(energyval, muval, vpa(iv), akx(ikx), iz, is, trapped, &
                                                 integrand_tmp_pls, integrand_tmp_min, Q_hat_z, &
                                                 integrand_tmp_v)

                  ! Split into contributions that are even and odd in vpa
                  RH_integrand_even(ikx,iz,it,ivmu) = 0.5*(integrand_tmp_pls+integrand_tmp_min)
                  RH_integrand_odd( ikx,iz,it,ivmu) = 0.5*(integrand_tmp_pls-integrand_tmp_min)
                  !> Apply the geometric factor I; only the geometry knows it, and
                  !> for a quasisymmetric field it is (MG+NI)/(N-iota M).
                  RH_pmom_weight(   ikx,iz,it,ivmu) = integrand_tmp_v * RH_drift_phase_fac(iz)

                  if (write_RH_asymptotics) then
                     call get_RH_LW_weights(energyval, muval, vpa(iv), akx(ikx), iz, is, &
                                                    trapped, LW_even_tmp, LW_odd_tmp, Q_hat_z)
                     RH_LW_even(ikx,iz,it,ivmu) = LW_even_tmp
                     RH_LW_odd( ikx,iz,it,ivmu) = LW_odd_tmp
                     call get_RH_SW_weights(energyval, muval, vpa(iv), akx(ikx), iz, is, &
                                            trapped, SW_even_tmp, SW_odd_tmp, Q_hat_z)
                     RH_SW_even(ikx,iz,it,ivmu) = SW_even_tmp
                     RH_SW_odd( ikx,iz,it,ivmu) = SW_odd_tmp
                  end if

               end do !ikx

               !> RH_U_parallel_fac is the same construction evaluated at a tiny
               !> kx, so it does not depend on <ikx> and is evaluated once per
               !> (iz, it, ivmu) rather than nakx times.
               call get_RH_transit_integrands(energyval, muval, vpa(iv), kxsmall, iz, is, trapped, &
                                              integrand_tmp_pls, integrand_tmp_min, Q_hat_z)

               RH_U_parallel_fac(iz,ivmu) = real( (1 - 0.5*(integrand_tmp_pls-integrand_tmp_min))/(zi*kxsmall) &
                                               * spec(is)%z/spec(is)%mass )

            end do !it
         end do !iz
      end do !ivmu

      deallocate (Q_hat_z)

      !> Flux tube only.  init_rosenbluth_hinton aborts above for
      !> full_flux_surface and radial_variation, so this is reached only where
      !> the flux-tube form is the right one.
      call get_RH_inertia_fluxtube()
      call get_RH_pmom_inertia_fluxtube()

   end subroutine init_rosenbluth_hinton


   !============================================================================
   !======================== FINALIZE THE DIAGNOSTICS ==========================
   !============================================================================
   subroutine finish_rosenbluth_hinton()

      use mp, only: proc0

      implicit none

      ! Deallocate the arrays for the Rosenbluth-Hinton integrand term
      if (allocated(RH_integrand_even)) deallocate (RH_integrand_even)
      if (allocated(RH_integrand_odd )) deallocate (RH_integrand_odd)
      if (allocated(RH_pmom_weight )) deallocate (RH_pmom_weight)
      if (allocated(RH_pmom_inertia)) deallocate (RH_pmom_inertia)
      if (allocated(RH_LW_even)) deallocate (RH_LW_even)
      if (allocated(RH_LW_odd )) deallocate (RH_LW_odd)
      if (allocated(RH_SW_even)) deallocate (RH_SW_even)
      if (allocated(RH_SW_odd )) deallocate (RH_SW_odd)
      if (allocated(RH_U_parallel_fac)) deallocate (RH_U_parallel_fac)
      if (allocated(RH_inertia))        deallocate (RH_inertia)
      if (allocated(RH_drift_bounce_avg)) deallocate (RH_drift_bounce_avg)
      if (allocated(RH_drift_is_trapped)) deallocate (RH_drift_is_trapped)

      rosenbluth_hinton_initialized = .false.

   end subroutine finish_rosenbluth_hinton

   !> Whether anything has asked for a quantity built on the drift-orbit phase.
   !> The bounce-averaged radial drift is not one of them: it needs only the
   !> geometry, which is what lets it be diagnosed in a stellarator where the

   !============================================================================
   !=========== IS THE ROSENBLUTH-HINTON MACHINERY NEEDED AT ALL? ==============
   !============================================================================
   !> Single source of truth for the init/finish guard, so the two can never
   !> disagree and leak the (large) response arrays.
   logical function rosenbluth_hinton_needed()

      use parameters_diagnostics, only: write_RH_inertia_fluxes, write_RH_bounce_drift
      use parameters_physics, only: omprimfac_RH
      use parameters_physics, only: zonal_init_option_switch, zonal_init_none
      use parameters_physics, only: zonal_closure_option_switch, zonal_closure_rh

      implicit none

      rosenbluth_hinton_needed = write_RH_inertia_fluxes &
                                 .or. write_RH_bounce_drift &
                                 .or. abs(omprimfac_RH) > epsilon(0.) &
                                 .or. (zonal_init_option_switch /= zonal_init_none &
                                       .and. zonal_closure_option_switch == zonal_closure_rh)

   end function rosenbluth_hinton_needed


   !============================================================================
   !====================== GET RH_inertia FOR THE FLUX TUBE =====================
   !============================================================================
   subroutine get_RH_inertia_fluxtube()

      use zgrid, only: nzgrid, ntubes
      use species, only: spec, nspec
      use vpamu_grids, only: vpa, vperp2, integrate_vmu
      use vpamu_grids, only: maxwell_mu, ztmax, maxwell_fac, maxwell_vpa
      use parameters_kxky_grids, only: naky, nakx, nx
      use grids_kxky, only: aky
      use calculations_kxky, only: multiply_by_rho
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use gyro_averages, only: aj0x, gyro_average
      use arrays_fields, only: phi
      use parameters_numerical, only: fphi
      use parameters_numerical, only: maxwellian_normalization
      use constants, only: zi
      
      ! Import temp array g1 with dimension (nky, nkx, -nzgrid:nzgrid, ntubes, -vmu-layout-)
      use arrays_dist_fn, only: integrand_vpamu => g1

      implicit none

!      ! The RH inertia is returned with dimensions (kx, z, tube, spec)
!      complex, dimension(:, -nzgrid:, :, :), intent(out) :: RH_inertia

      ! Temp variable holding RH inertia with dimensions (ky, kx, z, tube, spec) (1st is dummy)
      complex, dimension(:, :, :, :, :), allocatable :: RH_inertia_tmp

      ! Local variables
      integer :: ivmu, iv, imu, is, ia, iz, it
      
      ! We only have one field line because <full_flux_surface> = .false.
      ia = 1

      allocate (RH_inertia_tmp(naky, nakx, -nzgrid:nzgrid, ntubes, nspec)); RH_inertia_tmp = 0.

      if (.not. allocated(integrand_vpamu)) &
         allocate (integrand_vpamu(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
 
      !=========================================================================
      !                     ROSENBLUTH-HINTON INERTIA                          !
      !=========================================================================
      ! The Rosenbluth-Hinton inertia is calculated as:
      !		<RH_inertia> = - sum_s Z_s^2 e/T_s * velocity_integral( F_Ms * (1 - <J_0s exp(-i*Q_s)>_tau * J_0s exp(i*Q_s)) )
      !=========================================================================
      
      integrand_vpamu = 0.

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)

         do it = 1, ntubes
            do iz = -nzgrid, nzgrid

                integrand_vpamu(1, :, iz, it, ivmu) = spec(is)%zt &
                     * (1 - aj0x(1, :, iz, ivmu) * (RH_integrand_even(:, iz, it, ivmu) &
                                                    + RH_integrand_odd(:, iz, it, ivmu)))
                !> integrate_vmu folds the Maxwellian into its own weights when
                !> the evolved pdf is normalised by one, so applying it here too
                !> would count it twice.
                if (.not. maxwellian_normalization) &
                   integrand_vpamu(1, :, iz, it, ivmu) = integrand_vpamu(1, :, iz, it, ivmu) &
                      * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is)

            end do
         end do

      end do
      
      ! Calculate RH_inertia
      call integrate_vmu(integrand_vpamu, spec%dens_psi0*spec%z, RH_inertia_tmp)

      RH_inertia(:,:,:,:) = RH_inertia_tmp(1,:,:,:,:)

      deallocate (RH_inertia_tmp)

   end subroutine get_RH_inertia_fluxtube

   !============================================================================
   !====================== GET RH_fluxes FOR THE FLUX TUBE =====================
   !============================================================================
   subroutine get_RH_fluxes_fluxtube(g, RH_fluxes_phi_even,  RH_fluxes_phi_odd, &
                                        RH_fluxes_apar_even, RH_fluxes_apar_odd, &
                                        RH_fluxes_bpar_even, RH_fluxes_bpar_odd, &
                                        RH_fluxes_coll, RH_fluxes_drift_trapped, RH_fluxes_drift_passing, &
                                        RH_fluxes_coll_even, RH_fluxes_coll_odd, &
                                        RH_fluxes_phi_even_LW, RH_fluxes_phi_odd_LW, &
                                        RH_fluxes_coll_even_LW, RH_fluxes_coll_odd_LW, &
                                        RH_fluxes_phi_even_SW, RH_fluxes_phi_odd_SW, &
                                        RH_fluxes_coll_even_SW, RH_fluxes_coll_odd_SW, &
                                        RH_fluxes_phi_even_rey, RH_fluxes_phi_odd_rey, &
                                        RH_fluxes_phi_even_dia, RH_fluxes_phi_odd_dia)

      use zgrid, only: nzgrid, ntubes
      use species, only: spec, nspec
      use vpamu_grids, only: vpa, mu, vperp2, integrate_vmu
      use vpamu_grids, only: maxwell_mu, ztmax, maxwell_fac, maxwell_vpa
      use parameters_kxky_grids, only: naky, nakx, nx
      use grids_kxky, only: aky, akx
      use calculations_kxky, only: multiply_by_rho
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use gyro_averages, only: gyro_average, gyro_average_j1, aj0x
      use arrays_fields, only: phi, apar, bpar
      use parameters_numerical, only: maxwellian_normalization
      use stella_transforms, only: transform_kx2x_xfirst, transform_x2kx_xfirst
      use constants, only: zi
      use parameters_physics, only: nonlinear, xdriftknob
      use geometry, only: exb_nonlin_fac, geo_surf, q_as_x
      use parameters_numerical, only: fphi
      use parameters_physics, only: include_apar, include_bpar
      use dissipation, only: include_collisions, collisions_implicit
      use dissipation, only: advance_collisions_explicit, advance_collisions_implicit
      use stella_time, only: code_dt

      ! Import temp arrays g1, g2 with dimensions (nky, nkx, -nzgrid:nzgrid, ntubes, -vmu-layout-)
      use arrays_dist_fn, only: gvmu
      use arrays_dist_fn, only: integrand_even   => g0
      use arrays_dist_fn, only: integrand_odd    => g1
      use parameters_diagnostics, only: write_RH_asymptotics, write_RH_stress_split

      implicit none

      ! Gyroaveraged ExB and NL term in k-space
      complex, dimension(naky, nakx) :: vchix_gyro, NL_term, coll_term

      ! Gyroaveraged ExB term, distribution function, and integrand in x and ky
      complex, dimension(naky, nx) :: vchix_gyro_ky_x, g_ky_x, NL_term_ky_x

      ! The distribution function enters with dimensions (ky, kx, z, tube, ivmus)
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g

      ! The RH fluxes are returned with dimensions (ky, kx, z, tube, s)
      complex, dimension(:, :, -nzgrid:, :, :), intent(out) :: RH_fluxes_phi_even,  RH_fluxes_phi_odd
      complex, dimension(:, :, -nzgrid:, :, :), intent(out) :: RH_fluxes_apar_even, RH_fluxes_apar_odd
      complex, dimension(:, :, -nzgrid:, :, :), intent(out) :: RH_fluxes_bpar_even, RH_fluxes_bpar_odd

      ! The RH collisional flux is returned with dimensions (kx, z, tube, s)
      complex, dimension(:, :, :, :, :), allocatable :: RH_fluxes_coll_tmp

      ! Copies shielding the simulation state from the collision time-advance
      complex, dimension(:, :, :, :), allocatable :: phi_copy, apar_copy, bpar_copy
      complex, dimension(:, :, :), allocatable :: gvmu_saved
      complex, dimension(   :, -nzgrid:, :, :), intent(out) :: RH_fluxes_coll
      complex, dimension(   :, -nzgrid:, :, :), intent(out), optional :: RH_fluxes_coll_even, RH_fluxes_coll_odd
      !> The same fluxes formed with the long-wavelength weights instead of the
      !> exact ones.  Comparing the two is a test of the expansion that does not
      !> care what the turbulence looks like: both are driven by the same fields,
      !> so any difference is the weight and nothing else.
      complex, dimension(:, :, -nzgrid:, :, :), intent(out), optional :: RH_fluxes_phi_even_LW, RH_fluxes_phi_odd_LW
      complex, dimension(   :, -nzgrid:, :, :), intent(out), optional :: RH_fluxes_coll_even_LW, RH_fluxes_coll_odd_LW
      !> The same again with the short-wavelength (stationary-phase) weights.
      complex, dimension(:, :, -nzgrid:, :, :), intent(out), optional :: RH_fluxes_phi_even_SW, RH_fluxes_phi_odd_SW
      complex, dimension(   :, -nzgrid:, :, :), intent(out), optional :: RH_fluxes_coll_even_SW, RH_fluxes_coll_odd_SW
      complex, dimension(:, :, :, :, :), allocatable :: work_coll
      complex, dimension(:, :, :, :, :), allocatable :: LW_int_even, LW_int_odd
      complex, dimension(:, :, :, :, :), allocatable :: SW_int_even, SW_int_odd

      !> The even nonlinear channel is a stress, and it is not the Reynolds
      !> stress alone.  The gyroaverage on the ExB velocity carries the FLR
      !> weight, and splitting it as J0 = 1 + (J0 - 1) splits the flux by
      !> velocity moment: the J0 -> 1 half is a density moment, which
      !> quasineutrality slaves to phi and which is therefore the Reynolds
      !> channel; the remainder carries (J0 - 1) ~ -kperp^2 vperp^2 / 4 Omega^2,
      !> so it is the vperp^2 moment -- the perpendicular pressure -- and is the
      !> diamagnetic channel.  Both are the same order in kperp*rho, so neither
      !> is a correction to the other, and the sign of their sum decides whether
      !> the zonal flow is driven or damped.
      !>
      !> Note this cannot instead be split into the adiabatic and non-adiabatic
      !> parts of g.  stella evolves g = h - (Z/T)<phi>F_M, so the adiabatic
      !> piece is -(Z/T)<phi>F_M, and advecting it gives identically zero: both
      !> factors of the beat then carry the same J0, the summand
      !> (k' x k'') J0(k') J0(k'') phi phi is antisymmetric, and it cancels.
      !> That is {<phi>,<phi>} = 0.  The whole flux comes from h, and the
      !> Reynolds/diamagnetic distinction lives in which moment of h is taken.
      complex, dimension(:, :, -nzgrid:, :, :), intent(out), optional :: RH_fluxes_phi_even_rey, RH_fluxes_phi_odd_rey
      complex, dimension(:, :, -nzgrid:, :, :), intent(out), optional :: RH_fluxes_phi_even_dia, RH_fluxes_phi_odd_dia
      complex, dimension(:, :, :, :, :), allocatable :: rey_int_even, rey_int_odd
      complex, dimension(naky, nakx) :: vchix_bare, NL_term_rey
      complex, dimension(naky, nx) :: vchix_bare_ky_x

      !> Drive from the transit-averaged radial magnetic drift, reported
      !> separately for the trapped and passing populations; their sum is the
      !> whole drive.  See RH_drift_is_trapped for why they are kept apart.
      complex, dimension(   :, -nzgrid:, :, :), intent(out) :: RH_fluxes_drift_trapped
      complex, dimension(   :, -nzgrid:, :, :), intent(out) :: RH_fluxes_drift_passing
      complex, dimension(:, :, :, :, :), allocatable :: RH_fluxes_drift_tmp
      real, dimension(:), allocatable :: drift_weight
      complex, dimension(:), allocatable :: boltzmann

      ! Local variables
      integer :: ivmu, iv, imu, is, ia, iz, it

      !> The caller passes these arrays whether or not it means to write them,
      !> so present() alone is always true and the extra work would be done on
      !> every step regardless.  Gate on the diagnostic flags as well.
      logical :: do_LW, do_stress

      ! We only have one field line because <full_flux_surface> = .false.
      ia = 1

      do_LW = write_RH_asymptotics .and. present(RH_fluxes_phi_even_LW)
      do_stress = write_RH_stress_split .and. present(RH_fluxes_phi_even_rey)

      !> Allocated here rather than lower down: the nonlinear block below is the
      !> first user, and allocating after it segfaults.
      if (do_stress) then
         RH_fluxes_phi_even_rey = 0.; RH_fluxes_phi_odd_rey = 0.
         RH_fluxes_phi_even_dia = 0.; RH_fluxes_phi_odd_dia = 0.
         allocate (rey_int_even(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         allocate (rey_int_odd( naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         rey_int_even = 0.; rey_int_odd = 0.
      end if
      if (do_LW) then
         RH_fluxes_phi_even_SW = 0.; RH_fluxes_phi_odd_SW = 0.
         allocate (SW_int_even(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         allocate (SW_int_odd( naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         SW_int_even = 0.; SW_int_odd = 0.
         RH_fluxes_phi_even_LW = 0.; RH_fluxes_phi_odd_LW = 0.
         allocate (LW_int_even(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         allocate (LW_int_odd( naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         LW_int_even = 0.; LW_int_odd = 0.
      end if

      ! Only compute RH fluxes for nonlinear run
      if (.not. nonlinear) then
         RH_fluxes_phi_even  = 0.
         RH_fluxes_phi_odd   = 0.
         RH_fluxes_apar_even = 0.
         RH_fluxes_apar_odd  = 0.
         RH_fluxes_bpar_even = 0.
         RH_fluxes_bpar_odd  = 0.
      else

         !=========================================================================
         !                     ROSENBLUTH-HINTON FLUXES                           !
         !=========================================================================
         ! The Rosenbluth-Hinton fluxes are calculated as:
         !		<RH_fluxes>(even/odd) = -Z_s * velocity_integral( <J_0s exp(-i*Q_s)>_tau * exp(i*Q_s) 
         !		                                    ( <vchi_ky>_R * nabla(x)  * conj(h_s(even/odd))_ky )_kx )
         ! We do this in the following steps
         ! 		ivmu, it, iz loop: <vchix_gyro> = vchix*J0 = i*ky*<chi> * <aj0x(iky, ikx, iz, ivmu)>
         ! 		FFT: g(ky,kx) -> g(ky, x)
         ! 		FFT: vchix_gyro(ky,kx) -> vchix_gyro(ky, x)
         ! 		<vchix_g_NL>(ky,x)  = <vchix_gyro>(ky,x)*conj(g(ky,x))
         ! 		IFFT: vchix_g_NL(ky,x) -> vchix_g_NL(ky,kx)
         ! 		<integrand_even> = -Zs*(<J_0s exp(-i*Q_s)>_tau * exp(i*Q_s))(even)*<vchix_g_NL>
         ! 		<integrand_odd>  = -Zs*(<J_0s exp(-i*Q_s)>_tau * exp(i*Q_s))(odd )*<vchix_g_NL>
         ! 		RH_fluxes  = integrate_vmu(integrand)
         !=========================================================================

         !!!!!!!!!!!!!!!!!!!!!!!!!
         !!! phi contribution  !!!
         !!!!!!!!!!!!!!!!!!!!!!!!!

         NL_term = 0.

         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            iv = iv_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            is = is_idx(vmu_lo, ivmu)
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                   call gyro_average(zi*fphi*spread(aky,2,nakx)*phi(:,:,iz,it), iz, ivmu, vchix_gyro)
                   call transform_kx2x_xfirst(vchix_gyro, vchix_gyro_ky_x)
                   call transform_kx2x_xfirst(g(:,:,iz,it,ivmu), g_ky_x)
                   NL_term_ky_x = 2*real(vchix_gyro_ky_x * conjg(g_ky_x)) *exb_nonlin_fac
                   call transform_x2kx_xfirst(NL_term_ky_x, NL_term)
                   integrand_even(:,:,iz,it,ivmu) = NL_term * spread(RH_integrand_even(:,iz,it,ivmu), 1, naky)
                   integrand_odd( :,:,iz,it,ivmu) = NL_term * spread(RH_integrand_odd( :,iz,it,ivmu), 1, naky)
                   if (do_LW) then
                      LW_int_even(:,:,iz,it,ivmu) = NL_term * spread(RH_LW_even(:,iz,it,ivmu), 1, naky)
                      LW_int_odd( :,:,iz,it,ivmu) = NL_term * spread(RH_LW_odd( :,iz,it,ivmu), 1, naky)
                      SW_int_even(:,:,iz,it,ivmu) = NL_term * spread(RH_SW_even(:,iz,it,ivmu), 1, naky)
                      SW_int_odd( :,:,iz,it,ivmu) = NL_term * spread(RH_SW_odd( :,iz,it,ivmu), 1, naky)
                   end if
                   if (do_stress) then
                      !> The same flux with the gyroaverage taken off the ExB
                      !> velocity.  This is the density-moment (Reynolds) half;
                      !> the difference from the full flux carries (J0 - 1) and
                      !> so is the vperp^2, or pressure, half.
                      vchix_bare = zi*fphi*spread(aky,2,nakx)*phi(:,:,iz,it)
                      call transform_kx2x_xfirst(vchix_bare, vchix_bare_ky_x)
                      NL_term_ky_x = 2*real(vchix_bare_ky_x * conjg(g_ky_x)) * exb_nonlin_fac
                      call transform_x2kx_xfirst(NL_term_ky_x, NL_term_rey)
                      rey_int_even(:,:,iz,it,ivmu) = NL_term_rey * spread(RH_integrand_even(:,iz,it,ivmu), 1, naky)
                      rey_int_odd( :,:,iz,it,ivmu) = NL_term_rey * spread(RH_integrand_odd( :,iz,it,ivmu), 1, naky)
                   end if
               end do
            end do
         end do

         ! Calculate <RH_fluxes>(even/odd)
         call integrate_vmu(integrand_even, spec%dens_psi0*spec%z, RH_fluxes_phi_even)
         call integrate_vmu(integrand_odd,  spec%dens_psi0*spec%z, RH_fluxes_phi_odd)
         if (do_LW) then
            call integrate_vmu(LW_int_even, spec%dens_psi0*spec%z, RH_fluxes_phi_even_LW)
            call integrate_vmu(LW_int_odd,  spec%dens_psi0*spec%z, RH_fluxes_phi_odd_LW)
            call integrate_vmu(SW_int_even, spec%dens_psi0*spec%z, RH_fluxes_phi_even_SW)
            call integrate_vmu(SW_int_odd,  spec%dens_psi0*spec%z, RH_fluxes_phi_odd_SW)
         end if
         if (do_stress) then
            call integrate_vmu(rey_int_even, spec%dens_psi0*spec%z, RH_fluxes_phi_even_rey)
            call integrate_vmu(rey_int_odd,  spec%dens_psi0*spec%z, RH_fluxes_phi_odd_rey)
            !> Every step is linear in the velocity, so the pressure half is
            !> exactly the remainder.  Taking it by difference rather than by a
            !> second convolution makes total = Reynolds + diamagnetic an
            !> identity rather than a check that could drift, and costs one
            !> transform instead of two.
            RH_fluxes_phi_even_dia = RH_fluxes_phi_even - RH_fluxes_phi_even_rey
            RH_fluxes_phi_odd_dia  = RH_fluxes_phi_odd  - RH_fluxes_phi_odd_rey
            deallocate (rey_int_even, rey_int_odd)
         end if

         !!!!!!!!!!!!!!!!!!!!!!!!!
         !!! apar contribution !!!
         !!!!!!!!!!!!!!!!!!!!!!!!!
         if (include_apar) then
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               iv = iv_idx(vmu_lo, ivmu)
               imu = imu_idx(vmu_lo, ivmu)
               is = is_idx(vmu_lo, ivmu)
               do it = 1, ntubes
                  do iz = -nzgrid, nzgrid
                         call gyro_average(-2.0 * vpa(iv)*spec(is)%stm_psi0 &
                                      * zi*spread(aky,2,nakx)*apar(:,:,iz,it), iz, ivmu, vchix_gyro)
                         call transform_kx2x_xfirst(vchix_gyro, vchix_gyro_ky_x)
                         call transform_kx2x_xfirst(g(:,:,iz,it,ivmu), g_ky_x)
                         NL_term_ky_x = 2*real(vchix_gyro_ky_x * conjg(g_ky_x)) *exb_nonlin_fac
                         call transform_x2kx_xfirst(NL_term_ky_x, NL_term)
                         ! Note odd/even swap because of v_parallel factor in vchi_x
                         integrand_odd( :,:,iz,it,ivmu) = NL_term *  spread(RH_integrand_even(:,iz,it,ivmu), 1, naky)
                         integrand_even(:,:,iz,it,ivmu) = NL_term *  spread(RH_integrand_odd( :,iz,it,ivmu), 1, naky)
                  end do
               end do
            end do

            ! Calculate <RH_fluxes>(even/odd)
            call integrate_vmu(integrand_even, spec%dens_psi0*spec%z, RH_fluxes_apar_even)
            call integrate_vmu(integrand_odd,  spec%dens_psi0*spec%z, RH_fluxes_apar_odd)

         endif

         !!!!!!!!!!!!!!!!!!!!!!!!!
         !!! bpar contribution !!!
         !!!!!!!!!!!!!!!!!!!!!!!!!
         if (include_bpar) then
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               iv = iv_idx(vmu_lo, ivmu)
               imu = imu_idx(vmu_lo, ivmu)
               is = is_idx(vmu_lo, ivmu)
               do it = 1, ntubes
                  do iz = -nzgrid, nzgrid
                         call gyro_average_j1( 4.0*mu(imu)*spec(is)%tz &
                                      * zi*spread(aky,2,nakx)*bpar(:,:,iz,it), iz, ivmu, vchix_gyro)

                         call transform_kx2x_xfirst(vchix_gyro, vchix_gyro_ky_x)
                         call transform_kx2x_xfirst(g(:,:,iz,it,ivmu), g_ky_x)
                         NL_term_ky_x = 2*real(vchix_gyro_ky_x * conjg(g_ky_x)) *exb_nonlin_fac
                         call transform_x2kx_xfirst(NL_term_ky_x, NL_term)
                         integrand_even(:,:,iz,it,ivmu) = NL_term *  spread(RH_integrand_even(:,iz,it,ivmu), 1, naky)
                         integrand_odd( :,:,iz,it,ivmu) = NL_term *  spread(RH_integrand_odd( :,iz,it,ivmu), 1, naky)
                  end do
               end do
            end do

            ! Calculate <RH_fluxes>(even/odd)
            call integrate_vmu(integrand_even, spec%dens_psi0*spec%z, RH_fluxes_bpar_even)
            call integrate_vmu(integrand_odd,  spec%dens_psi0*spec%z, RH_fluxes_bpar_odd)

         endif
      endif


      if (present(RH_fluxes_coll_even)) then
         RH_fluxes_coll_even = 0.; RH_fluxes_coll_odd = 0.
         allocate (work_coll(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      end if
      if (present(RH_fluxes_coll_even_LW)) then
         RH_fluxes_coll_even_LW = 0.; RH_fluxes_coll_odd_LW = 0.
      end if


      ! Only compute RH collisional flux when collisions are included
      if (.not. include_collisions) then
         RH_fluxes_coll = 0.
         if (allocated(work_coll)) deallocate (work_coll)
      else

         !!!!!!!!!!!!!!!!!!!!!!!!!
         !!! coll contribution !!!
         !!!!!!!!!!!!!!!!!!!!!!!!!
         allocate (RH_fluxes_coll_tmp(naky, nakx, -nzgrid:nzgrid, ntubes, nspec)); RH_fluxes_coll_tmp = 0.

         !> Evaluate dt * collision operator (stored in integrand_even).
         !>
         !> <advance_collisions_implicit> is a time-advance routine, not a
         !> side-effect-free evaluation of C[g]: it declares phi, apar and bpar
         !> intent(in out) and updates them as part of the implicit solve, and it
         !> overwrites the module-level <gvmu> through the scatter/gather to the
         !> kxkyz layout.  Calling it here with the live fields let a diagnostic
         !> corrupt the simulation state.  Electrostatically the damage was
         !> survivable, but with apar evolving it produced NaNs within ten steps.
         !> So hand it copies and put <gvmu> back afterwards.
         if (collisions_implicit) then
            allocate (phi_copy, source=phi)
            allocate (apar_copy, source=apar)
            allocate (bpar_copy, source=bpar)
            allocate (gvmu_saved, source=gvmu)

            integrand_even = g
            call advance_collisions_implicit(.false., phi_copy, apar_copy, bpar_copy, integrand_even)
            integrand_even = integrand_even - g

            gvmu = gvmu_saved
            deallocate (phi_copy, apar_copy, bpar_copy, gvmu_saved)
         else
            integrand_even = 0.
            call advance_collisions_explicit(g, phi, bpar, integrand_even)
         end if

         ! Remove nonzonal terms to save computational time
         integrand_even(2:,:,:,:,:) = 0.0

         !> Split the collisional flux the same way as the nonlinear one.  The
         !> two halves have opposite orderings at long wavelength: the operator
         !> conserves particles, so the even part loses its leading term exactly
         !> and is O(kx) while the odd part is O(1).  Keeping them apart is what
         !> lets that be checked.
         if (present(RH_fluxes_coll_even)) then
            work_coll = 1/code_dt * integrand_even * spread(RH_integrand_even, 1, naky)
            call integrate_vmu(work_coll, spec%dens_psi0*spec%z, RH_fluxes_coll_tmp)
            RH_fluxes_coll_even = RH_fluxes_coll_tmp(1,:,:,:,:)
            work_coll = 1/code_dt * integrand_even * spread(RH_integrand_odd, 1, naky)
            call integrate_vmu(work_coll, spec%dens_psi0*spec%z, RH_fluxes_coll_tmp)
            RH_fluxes_coll_odd = RH_fluxes_coll_tmp(1,:,:,:,:)
            if (present(RH_fluxes_coll_even_LW)) then
               work_coll = 1/code_dt * integrand_even * spread(RH_LW_even, 1, naky)
               call integrate_vmu(work_coll, spec%dens_psi0*spec%z, RH_fluxes_coll_tmp)
               RH_fluxes_coll_even_LW = RH_fluxes_coll_tmp(1,:,:,:,:)
               work_coll = 1/code_dt * integrand_even * spread(RH_LW_odd, 1, naky)
               call integrate_vmu(work_coll, spec%dens_psi0*spec%z, RH_fluxes_coll_tmp)
               RH_fluxes_coll_odd_LW = RH_fluxes_coll_tmp(1,:,:,:,:)
               work_coll = 1/code_dt * integrand_even * spread(RH_SW_even, 1, naky)
               call integrate_vmu(work_coll, spec%dens_psi0*spec%z, RH_fluxes_coll_tmp)
               RH_fluxes_coll_even_SW = RH_fluxes_coll_tmp(1,:,:,:,:)
               work_coll = 1/code_dt * integrand_even * spread(RH_SW_odd, 1, naky)
               call integrate_vmu(work_coll, spec%dens_psi0*spec%z, RH_fluxes_coll_tmp)
               RH_fluxes_coll_odd_SW = RH_fluxes_coll_tmp(1,:,:,:,:)
            end if
         end if

         ! Evaluate integrand in RH collisional flux
         integrand_odd = 1/code_dt * integrand_even * spread(RH_integrand_even+RH_integrand_odd, 1, naky)

         ! Integrate over velocity space
         call integrate_vmu(integrand_odd, spec%dens_psi0*spec%z, RH_fluxes_coll_tmp)
         RH_fluxes_coll = RH_fluxes_coll_tmp(1,:,:,:,:)

         ! Note : extra factor of -1/(1j*kx) to match definition of nonlinear fluxes
         if (abs(akx(1)) < epsilon(0.)) then
             RH_fluxes_coll(1, :,:,:) = 0.0
             RH_fluxes_coll(2:,:,:,:) = -RH_fluxes_coll(2:,:,:,:)/(zi*spread(spread(spread(akx(2:),2,2*nzgrid+1),3,ntubes),4,nspec))
         else
             RH_fluxes_coll(1:,:,:,:) = -RH_fluxes_coll(1:,:,:,:)/(zi*spread(spread(spread(akx(1:),2,2*nzgrid+1),3,ntubes),4,nspec))
         end if

         !> The asymptotic versions need the same -1/(i kx) as the exact ones;
         !> without it they are short by exactly a factor of kx.
         if (present(RH_fluxes_coll_even_LW)) then
            if (abs(akx(1)) < epsilon(0.)) then
               RH_fluxes_coll_even_LW(1,:,:,:) = 0.; RH_fluxes_coll_odd_LW(1,:,:,:) = 0.
               RH_fluxes_coll_even_LW(2:,:,:,:) = -RH_fluxes_coll_even_LW(2:,:,:,:) &
                  / (zi*spread(spread(spread(akx(2:),2,2*nzgrid+1),3,ntubes),4,nspec))
               RH_fluxes_coll_odd_LW(2:,:,:,:) = -RH_fluxes_coll_odd_LW(2:,:,:,:) &
                  / (zi*spread(spread(spread(akx(2:),2,2*nzgrid+1),3,ntubes),4,nspec))
            else
               RH_fluxes_coll_even_LW = -RH_fluxes_coll_even_LW &
                  / (zi*spread(spread(spread(akx(1:),2,2*nzgrid+1),3,ntubes),4,nspec))
               RH_fluxes_coll_odd_LW = -RH_fluxes_coll_odd_LW &
                  / (zi*spread(spread(spread(akx(1:),2,2*nzgrid+1),3,ntubes),4,nspec))
            end if
            if (abs(akx(1)) < epsilon(0.)) then
               RH_fluxes_coll_even_SW(1,:,:,:) = 0.; RH_fluxes_coll_odd_SW(1,:,:,:) = 0.
               RH_fluxes_coll_even_SW(2:,:,:,:) = -RH_fluxes_coll_even_SW(2:,:,:,:) &
                  / (zi*spread(spread(spread(akx(2:),2,2*nzgrid+1),3,ntubes),4,nspec))
               RH_fluxes_coll_odd_SW(2:,:,:,:) = -RH_fluxes_coll_odd_SW(2:,:,:,:) &
                  / (zi*spread(spread(spread(akx(2:),2,2*nzgrid+1),3,ntubes),4,nspec))
            else
               RH_fluxes_coll_even_SW = -RH_fluxes_coll_even_SW &
                  / (zi*spread(spread(spread(akx(1:),2,2*nzgrid+1),3,ntubes),4,nspec))
               RH_fluxes_coll_odd_SW = -RH_fluxes_coll_odd_SW &
                  / (zi*spread(spread(spread(akx(1:),2,2*nzgrid+1),3,ntubes),4,nspec))
            end if
         end if

         if (present(RH_fluxes_coll_even)) then
            if (abs(akx(1)) < epsilon(0.)) then
               RH_fluxes_coll_even(1,:,:,:) = 0.; RH_fluxes_coll_odd(1,:,:,:) = 0.
               RH_fluxes_coll_even(2:,:,:,:) = -RH_fluxes_coll_even(2:,:,:,:) &
                  / (zi*spread(spread(spread(akx(2:),2,2*nzgrid+1),3,ntubes),4,nspec))
               RH_fluxes_coll_odd(2:,:,:,:) = -RH_fluxes_coll_odd(2:,:,:,:) &
                  / (zi*spread(spread(spread(akx(2:),2,2*nzgrid+1),3,ntubes),4,nspec))
            else
               RH_fluxes_coll_even = -RH_fluxes_coll_even &
                  / (zi*spread(spread(spread(akx(1:),2,2*nzgrid+1),3,ntubes),4,nspec))
               RH_fluxes_coll_odd = -RH_fluxes_coll_odd &
                  / (zi*spread(spread(spread(akx(1:),2,2*nzgrid+1),3,ntubes),4,nspec))
            end if
         end if

         deallocate(RH_fluxes_coll_tmp)
         if (allocated(work_coll)) deallocate (work_coll)

      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! drift contribution !!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!
      !> The transit average is built so that parallel streaming together with
      !> the radial drift annihilates the Rosenbluth-Hinton projection: that is
      !> what v_par b.grad Q = i kx (v_Mx - <v_Mx>_tau) says.  What survives is
      !> exactly the bounce-averaged part, i kx <v_Mx>_b, which is why this term
      !> is absent from a tokamak -- where quasisymmetry makes <v_Mx>_b vanish --
      !> and present in a general stellarator.  It is not a nonlinear term and
      !> does not need collisions, so unlike the two above it is evaluated for
      !> every run.
      !>
      !> The other fluxes are defined with a factor -1/(i kx) relative to their
      !> source term.  Here the source is itself -i kx sum_s Z_s n_s
      !> integral(W <v_Mx>_b g), so the two factors of kx cancel and what is left
      !> is a plain velocity integral.
      !> Only where Q was integrated along the field line.  The closed form is
      !> derived on the assumption that the transit-averaged drift vanishes, so
      !> pairing it with a numerically non-zero one would leave the cancellation
      !> between streaming and the drift incomplete and put a spurious term in
      !> the budget.  Zero here is the consistent answer, and in the axisymmetric
      !> geometry the closed form applies to, the true value anyway.
      if (use_analytic_drift_phase) then
         RH_fluxes_drift_trapped = 0.
         RH_fluxes_drift_passing = 0.
         return
      end if

      allocate (RH_fluxes_drift_tmp(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))
      RH_fluxes_drift_tmp = 0.

      !> RH_drift_bounce_avg holds the bounce average of stella's geometric
      !> drift, cvdrift0 vpa^2 + gbdrift0 vperp^2 / 2.  time_advance turns that
      !> into a drift frequency with 0.5 * tz_psi0, divided by shat unless
      !> q_as_x; that normalisation is applied here rather than being folded into
      !> the stored array, which stays the pure geometric quantity the
      !> write_RH_bounce_drift diagnostic reports.
      allocate (drift_weight(nspec))
      allocate (boltzmann(nakx))
      drift_weight = 0.5 * xdriftknob * spec%dens_psi0 * spec%z * spec%tz_psi0
      if (.not. q_as_x) drift_weight = drift_weight / geo_surf%shat

      !> The drift acts on the full perturbed distribution, not on g alone.
      !> time_advance applies it twice: wdriftx_g against g, and wdriftx_phi
      !> against the gyroaveraged potential, the latter carrying an extra
      !> zt F_M.  Both are the same geometric drift, so what the projection
      !> leaves behind is i kx <v_Mx>_tau acting on
      !>
      !>    h = g + (Z/T) J_0 phi F_M,
      !>
      !> and keeping only the g piece of it accounts for a fraction of the drive.
      integrand_even = 0.
      integrand_odd = 0.
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               !> The Boltzmann part of the response, (Z/T) J_0 phi F_M.  Its
               !> Maxwellian is dropped when the evolved pdf already carries
               !> one, for the same reason as above.
               boltzmann = fphi * aj0x(1, :, iz, ivmu) * phi(1, :, iz, it) * spec(is)%zt
               if (.not. maxwellian_normalization) &
                  boltzmann = boltzmann * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is)

               if (RH_drift_is_trapped(iz, ivmu)) then
                  integrand_even(1, :, iz, it, ivmu) = &
                     (g(1, :, iz, it, ivmu) + boltzmann) &
                     * (RH_integrand_even(:, iz, it, ivmu) + RH_integrand_odd(:, iz, it, ivmu)) &
                     * RH_drift_bounce_avg(iz, ivmu)
               else
                  integrand_odd(1, :, iz, it, ivmu) = &
                     (g(1, :, iz, it, ivmu) + boltzmann) &
                     * (RH_integrand_even(:, iz, it, ivmu) + RH_integrand_odd(:, iz, it, ivmu)) &
                     * RH_drift_bounce_avg(iz, ivmu)
               end if
            end do
         end do
      end do

      call integrate_vmu(integrand_even, drift_weight, RH_fluxes_drift_tmp)
      RH_fluxes_drift_trapped = RH_fluxes_drift_tmp(1, :, :, :, :)
      call integrate_vmu(integrand_odd, drift_weight, RH_fluxes_drift_tmp)
      RH_fluxes_drift_passing = RH_fluxes_drift_tmp(1, :, :, :, :)

      deallocate (RH_fluxes_drift_tmp, drift_weight, boltzmann)

      if (allocated(LW_int_even)) deallocate (LW_int_even, LW_int_odd)
      if (allocated(SW_int_even)) deallocate (SW_int_even, SW_int_odd)

   end subroutine get_RH_fluxes_fluxtube
 

   !============================================================================
   !====================== GET RH_phi_I FOR THE FLUX TUBE ========================
   !============================================================================
   !> The toroidal-momentum inertia.  A rigid toroidal rotation of frequency
   !> omega has V_par = omega I/B, that is g = (m vpa/T)(omega I/B) F_M, so
   !> applying the projection to that state and dividing it out returns omega
   !> itself: <RH_pmom>/<RH_pmom_inertia> is the residual toroidal rotation
   !> frequency.  Built once from the geometry and the equilibrium Maxwellian.
   subroutine get_RH_pmom_inertia_fluxtube()

      use zgrid, only: nzgrid, ntubes
      use species, only: spec, nspec
      use vpamu_grids, only: vpa, integrate_vmu
      use vpamu_grids, only: maxwell_mu, maxwell_vpa, maxwell_fac
      use parameters_kxky_grids, only: naky, nakx
      use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      use parameters_numerical, only: maxwellian_normalization
      use geometry, only: RH_drift_phase_fac, bmag
      use arrays_dist_fn, only: integrand_vpamu => g1

      implicit none

      complex, dimension(:, :, :, :, :), allocatable :: inertia_tmp
      integer :: ivmu, iv, imu, is, ia, iz, it

      ia = 1
      allocate (inertia_tmp(naky, nakx, -nzgrid:nzgrid, ntubes, nspec)); inertia_tmp = 0.
      if (.not. allocated(integrand_vpamu)) &
         allocate (integrand_vpamu(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

      integrand_vpamu = 0.
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         do it = 1, ntubes
            do iz = -nzgrid, nzgrid
               !> The rotating state carries a further I/B, matching the one
               !> already inside the weight, so the inertia is quadratic in the
               !> geometric factor and the ratio returns a frequency.
               integrand_vpamu(1, :, iz, it, ivmu) = RH_pmom_weight(:, iz, it, ivmu) &
                    * vpa(iv) * spec(is)%mass / spec(is)%temp &
                    * RH_drift_phase_fac(iz) / bmag(ia, iz)
               !> integrate_vmu folds the Maxwellian into its own weights when
               !> the evolved pdf is normalised by one.
               if (.not. maxwellian_normalization) &
                  integrand_vpamu(1, :, iz, it, ivmu) = integrand_vpamu(1, :, iz, it, ivmu) &
                     * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is)
            end do
         end do
      end do

      call integrate_vmu(integrand_vpamu, spec%dens_psi0, inertia_tmp)
      RH_pmom_inertia(:,:,:,:) = inertia_tmp(1,:,:,:,:)

      deallocate (inertia_tmp)

   end subroutine get_RH_pmom_inertia_fluxtube


   !> The parallel-flow Rosenbluth-Hinton invariant,
   !>     RH_pmom = n_s * velocity_integral( V_sigma g_s ) ,
   !> which the linear parallel streaming and the non-secular part of the radial
   !> magnetic drift annihilate exactly.  Divided by RH_pmom_inertia it is the
   !> residual parallel flow.
   subroutine get_RH_pmom(g, RH_pmom, RH_pmom_g)

      use zgrid, only: nzgrid, ntubes
      use species, only: spec, nspec
      use vpamu_grids, only: integrate_vmu, vpa
      use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use parameters_kxky_grids, only: naky, nakx
      use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      use parameters_physics, only: include_apar
      use parameters_numerical, only: maxwellian_normalization
      use arrays_fields, only: apar
      use gyro_averages, only: gyro_average

      use arrays_dist_fn, only: RH_pmom_tmp => g1

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:, -nzgrid:, :, :), intent(out) :: RH_pmom
      complex, dimension(:, -nzgrid:, :, :), intent(out), optional :: RH_pmom_g

      complex, dimension(:, :, :, :, :), allocatable :: pmom_tmp, pdf
      complex, dimension(naky, nakx) :: gyro_apar, correction
      integer :: ivmu, iv, imu, is, iz, it

      allocate (pmom_tmp(naky, nakx, -nzgrid:nzgrid, ntubes, nspec)); pmom_tmp = 0.
      if (.not. allocated(RH_pmom_tmp)) &
         allocate (RH_pmom_tmp(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

      RH_pmom_tmp = 0.
      RH_pmom_tmp(1,:,:,:,:) = RH_pmom_weight

      !> Electromagnetically it is gbar, not g, that the parallel streaming
      !> advects: stella carries gbar = g + 2 (Z/T) stm vpa <Apar> F_M and
      !> converts between the two around the implicit solve.  So it is the
      !> projection of gbar that the annihilation leaves constant.  The
      !> difference is odd in vpa, which is why an even-dominated weight barely
      !> notices it while this one, being odd, sees it in full.
      allocate (pdf(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      pdf = g
      if (include_apar) then
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            iv = iv_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            is = is_idx(vmu_lo, ivmu)
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  call gyro_average(apar(:, :, iz, it), iz, ivmu, gyro_apar)
                  correction = 2.0 * spec(is)%zt * spec(is)%stm_psi0 * vpa(iv) * gyro_apar
                  if (.not. maxwellian_normalization) &
                     correction = correction * maxwell_vpa(iv, is) &
                                * maxwell_mu(1, iz, imu, is) * maxwell_fac(is)
                  pdf(:, :, iz, it, ivmu) = pdf(:, :, iz, it, ivmu) + correction
               end do
            end do
         end do
      end if

      call integrate_vmu(pdf * RH_pmom_tmp, spec%dens_psi0, pmom_tmp)
      RH_pmom = pmom_tmp(1,:,:,:,:)

      !> The same projection taken of g alone.  The difference is the Apar
      !> part of gbar, and its time derivative is a source of the gbar
      !> projection that no flux channel carries; writing both is what lets
      !> that be separated from the fluxes rather than lumped into a residual.
      if (present(RH_pmom_g)) then
         call integrate_vmu(g * RH_pmom_tmp, spec%dens_psi0, pmom_tmp)
         RH_pmom_g = pmom_tmp(1,:,:,:,:)
      end if
      deallocate (pdf)

      deallocate (pmom_tmp)

   end subroutine get_RH_pmom


   !> The fluxes driving the parallel-flow invariant.  Same three channels as the
   !> potential-like fluxes -- nonlinear, collisional and the transit-averaged
   !> radial drift -- and the same conventions, in particular the factor
   !> -1/(i kx) that turns a source term into a flux so that
   !>
   !>     d<RH_pmom>/dt|_channel = -i kx <F_channel>.
   !>
   !> The only differences from get_RH_fluxes_fluxtube are the weight, which is
   !> V_sigma rather than the sigma-even and sigma-odd pair, and the species
   !> factor, which is n_s alone because a flow carries no charge weighting.
   !> There is correspondingly no even/odd split: V_sigma is one array.
   subroutine get_RH_pmom_fluxes_fluxtube(g, RH_pmom_flux_nl, RH_pmom_flux_coll, RH_pmom_flux_drift)

      use zgrid, only: nzgrid, ntubes
      use species, only: spec, nspec
      use vpamu_grids, only: vpa, mu, vperp2, integrate_vmu
      use vpamu_grids, only: maxwell_mu, maxwell_fac, maxwell_vpa
      use parameters_kxky_grids, only: naky, nakx, nx
      use grids_kxky, only: aky, akx
      use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      use gyro_averages, only: gyro_average, gyro_average_j1, aj0x
      use arrays_fields, only: phi, apar, bpar
      use parameters_numerical, only: maxwellian_normalization, fphi
      use stella_transforms, only: transform_kx2x_xfirst, transform_x2kx_xfirst
      use constants, only: zi
      use parameters_physics, only: nonlinear, xdriftknob
      use parameters_physics, only: include_apar, include_bpar
      use geometry, only: exb_nonlin_fac, geo_surf, q_as_x
      use dissipation, only: include_collisions, collisions_implicit
      use dissipation, only: advance_collisions_explicit, advance_collisions_implicit
      use stella_time, only: code_dt

      use arrays_dist_fn, only: gvmu
      use arrays_dist_fn, only: integrand   => g0
      use arrays_dist_fn, only: work        => g1

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      !> Keeps its ky axis: the nonlinear drive of a zonal mode is carried by the
      !> ky beat, and vchix ~ i ky phi is identically zero on the ky = 0 row, so
      !> collapsing to that row would discard the whole term.  The consumer sums
      !> over ky, exactly as for the potential-like nonlinear fluxes.
      complex, dimension(:, :, -nzgrid:, :, :), intent(out) :: RH_pmom_flux_nl
      complex, dimension(:, -nzgrid:, :, :), intent(out) :: RH_pmom_flux_coll
      complex, dimension(:, -nzgrid:, :, :), intent(out) :: RH_pmom_flux_drift

      complex, dimension(naky, nakx) :: vchix_gyro, vchix_part, NL_term
      complex, dimension(naky, nx)   :: vchix_gyro_ky_x, g_ky_x, NL_term_ky_x
      complex, dimension(:, :, :, :, :), allocatable :: flux_tmp
      complex, dimension(:, :, :, :), allocatable :: phi_copy, apar_copy, bpar_copy
      complex, dimension(:, :, :), allocatable :: gvmu_saved
      real,    dimension(:), allocatable :: drift_weight
      complex, dimension(:), allocatable :: boltzmann

      integer :: ivmu, iv, imu, is, ia, iz, it

      ia = 1
      allocate (flux_tmp(naky, nakx, -nzgrid:nzgrid, ntubes, nspec))

      !----------------------------- nonlinear -------------------------------
      RH_pmom_flux_nl = 0.
      if (nonlinear) then
         integrand = 0.
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            iv = iv_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            is = is_idx(vmu_lo, ivmu)
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  !> The gyroaveraged radial ExB velocity from the whole of chi.
                  !> The three field contributions add; gyro_average overwrites
                  !> its output, so they are accumulated rather than chained.
                  call gyro_average(zi*fphi*spread(aky,2,nakx)*phi(:,:,iz,it), iz, ivmu, vchix_gyro)
                  if (include_apar) then
                     call gyro_average(-2.0 * vpa(iv)*spec(is)%stm_psi0 &
                                       * zi*spread(aky,2,nakx)*apar(:,:,iz,it), iz, ivmu, vchix_part)
                     vchix_gyro = vchix_gyro + vchix_part
                  end if
                  if (include_bpar) then
                     call gyro_average_j1(4.0*mu(imu)*spec(is)%tz &
                                          * zi*spread(aky,2,nakx)*bpar(:,:,iz,it), iz, ivmu, vchix_part)
                     vchix_gyro = vchix_gyro + vchix_part
                  end if
                  call transform_kx2x_xfirst(vchix_gyro, vchix_gyro_ky_x)
                  !> g, not gbar, here.  The nonlinearity advects g by the ExB
                  !> velocity; it is the conserved projection, not the flux, that
                  !> carries the Apar term.  Substituting gbar here was tried and
                  !> makes the electromagnetic budget worse (0.29 -> 0.44).
                  call transform_kx2x_xfirst(g(:,:,iz,it,ivmu), g_ky_x)
                  NL_term_ky_x = 2*real(vchix_gyro_ky_x * conjg(g_ky_x)) * exb_nonlin_fac
                  call transform_x2kx_xfirst(NL_term_ky_x, NL_term)
                  integrand(:,:,iz,it,ivmu) = NL_term * spread(RH_pmom_weight(:,iz,it,ivmu), 1, naky)
               end do
            end do
         end do
         call integrate_vmu(integrand, spec%dens_psi0, RH_pmom_flux_nl)
      end if

      !---------------------------- collisional ------------------------------
      RH_pmom_flux_coll = 0.
      if (include_collisions) then
         !> advance_collisions_implicit is a time advance, not a side-effect-free
         !> evaluation of C[g]; hand it copies and restore gvmu, exactly as the
         !> potential-like flux does.
         if (collisions_implicit) then
            allocate (phi_copy, source=phi)
            allocate (apar_copy, source=apar)
            allocate (bpar_copy, source=bpar)
            allocate (gvmu_saved, source=gvmu)
            integrand = g
            call advance_collisions_implicit(.false., phi_copy, apar_copy, bpar_copy, integrand)
            integrand = integrand - g
            gvmu = gvmu_saved
            deallocate (phi_copy, apar_copy, bpar_copy, gvmu_saved)
         else
            integrand = 0.
            call advance_collisions_explicit(g, phi, bpar, integrand)
         end if

         ! Zonal modes only
         integrand(2:,:,:,:,:) = 0.0

         work = 1/code_dt * integrand * spread(RH_pmom_weight, 1, naky)
         flux_tmp = 0.
         call integrate_vmu(work, spec%dens_psi0, flux_tmp)
         RH_pmom_flux_coll = flux_tmp(1,:,:,:,:)

         !> -1/(i kx), matching the nonlinear flux convention.
         if (abs(akx(1)) < epsilon(0.)) then
            RH_pmom_flux_coll(1, :,:,:) = 0.0
            RH_pmom_flux_coll(2:,:,:,:) = -RH_pmom_flux_coll(2:,:,:,:) &
               / (zi*spread(spread(spread(akx(2:),2,2*nzgrid+1),3,ntubes),4,nspec))
         else
            RH_pmom_flux_coll(1:,:,:,:) = -RH_pmom_flux_coll(1:,:,:,:) &
               / (zi*spread(spread(spread(akx(1:),2,2*nzgrid+1),3,ntubes),4,nspec))
         end if
      end if

      !------------------------------- drift ---------------------------------
      !> As for the potential-like drift flux, the two factors of kx cancel and
      !> what is left is a plain velocity integral; and it is only meaningful
      !> where Q was integrated along the field line, the closed form being
      !> derived on the assumption that the transit-averaged drift vanishes.
      RH_pmom_flux_drift = 0.
      if (.not. use_analytic_drift_phase) then
         allocate (drift_weight(nspec))
         allocate (boltzmann(nakx))
         drift_weight = 0.5 * xdriftknob * spec%dens_psi0 * spec%tz_psi0
         if (.not. q_as_x) drift_weight = drift_weight / geo_surf%shat

         integrand = 0.
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            iv = iv_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            is = is_idx(vmu_lo, ivmu)
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  !> The drift acts on h = g + (Z/T) J0 phi F_M, not on g alone.
                  boltzmann = fphi * aj0x(1, :, iz, ivmu) * phi(1, :, iz, it) * spec(is)%zt
                  if (.not. maxwellian_normalization) &
                     boltzmann = boltzmann * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is)
                  integrand(1, :, iz, it, ivmu) = &
                     (g(1, :, iz, it, ivmu) + boltzmann) &
                     * RH_pmom_weight(:, iz, it, ivmu) * RH_drift_bounce_avg(iz, ivmu)
               end do
            end do
         end do

         flux_tmp = 0.
         call integrate_vmu(integrand, drift_weight, flux_tmp)
         RH_pmom_flux_drift = flux_tmp(1,:,:,:,:)
         deallocate (drift_weight, boltzmann)
      end if

      deallocate (flux_tmp)

   end subroutine get_RH_pmom_fluxes_fluxtube


   subroutine get_RH_phi_I_fluxtube(g, RH_phi_I)

      use zgrid, only: nzgrid, ntubes
      use species, only: spec, nspec
      use vpamu_grids, only: vpa, vperp2, integrate_vmu
      use vpamu_grids, only: maxwell_mu, ztmax, maxwell_fac, maxwell_vpa
      use parameters_kxky_grids, only: naky, nakx, nx
      use grids_kxky, only: aky
      use calculations_kxky, only: multiply_by_rho
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use arrays_fields, only: phi
      use parameters_numerical, only: maxwellian_normalization
      use constants, only: zi

      ! Import temp array g1 with dimension (nky, nkx, -nzgrid:nzgrid, ntubes, -vmu-layout-)
      use arrays_dist_fn, only: RH_integrand_tmp => g1

      implicit none

      ! The distribution function enters with dimensions (ky, kx, z, tube, ivmus)
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g

      ! The RH phi is returned with dimensions (kx, z, tube, s)
      complex, dimension(:, -nzgrid:, :, :), intent(out) :: RH_phi_I

      ! Temp variable holding RH_phi_I with dimensions (ky, kx, z, tube, spec) (1st is dummy)
      complex, dimension(:, :, :, :, :), allocatable :: RH_phi_I_tmp

      ! Local variables
      integer :: ivmu, iv, imu, is, ia, iz, it

      ! We only have one field line because <full_flux_surface> = .false.
      ia = 1

      allocate (RH_phi_I_tmp(naky, nakx, -nzgrid:nzgrid, ntubes, nspec)); RH_phi_I_tmp = 0.

      if (.not. allocated(RH_integrand_tmp)) &
         allocate (RH_integrand_tmp(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))

      !=========================================================================
      !                     ROSENBLUTH-HINTON POTENTIAL                        !
      !=========================================================================
      ! The Rosenbluth-Hinton potential is calculated (for zonal g_s) as:
      !		RH_phi_I = Z_s * velocity_integral( <J_0s exp(-i*Q_s)>_tau * exp(i*Q_s) g_s )
      !=========================================================================
      RH_integrand_tmp = 0.
      RH_integrand_tmp(1,:,:,:,:) = RH_integrand_even+RH_integrand_odd

      !> Weighted by Z_s n_s, as the inertia and every flux are.  This carried
      !> only Z_s, so RH_phi_I was short a density factor relative to everything
      !> it is paired with -- invisible wherever dens = 1, wrong otherwise, and
      !> exactly the kind of mismatch that stops RH_phi_I being the inertia times
      !> phi.  It cancels from the budget slope, appearing in dE_RH/dt and in
      !> P_RH alike, which is why the benchmarks never saw it.
      call integrate_vmu(g * RH_integrand_tmp, spec%dens_psi0 * spec%z, RH_phi_I_tmp)
      RH_phi_I = RH_phi_I_tmp(1,:,:,:,:)

      deallocate (RH_phi_I_tmp)

   end subroutine get_RH_phi_I_fluxtube

   !==============================================
   !============== BOUNCE AVERAGES ===============
   !============================================================================
   !=========== BOUNCE INTEGRALS WITH THE TURNING POINTS RESOLVED ==============
   !============================================================================
   !> Locate the magnetic well containing <iz> for a particle whose turning
   !> points are where B = B_c, and return the grid indices spanning it.
   !>
   !> A trapped particle is confined to one connected stretch of field line where
   !> B < B_c, and its bounce average belongs to that stretch alone.  Summing
   !> over every accessible stretch at once, as a plain masked sum over the whole
   !> domain does, mixes in wells the particle can never reach.  In a tokamak
   !> flux tube there is only one well, so the distinction has never mattered;
   !> in a stellarator it does.
   subroutine find_well(B_c, iz, iz_lo, iz_hi, well_found)

      use geometry, only: bmag
      use zgrid, only: nzgrid

      implicit none

      real,    intent(in)  :: B_c
      integer, intent(in)  :: iz
      integer, intent(out) :: iz_lo, iz_hi
      logical, intent(out) :: well_found

      integer :: ia
      ia = 1

      !> Not accessible at all, or the well runs off the end of the simulated
      !> field line, in which case it is not a complete bounce and the caller
      !> should fall back rather than integrate a truncated well.
      well_found = .false.
      iz_lo = iz; iz_hi = iz
      if (bmag(ia, iz) >= B_c) return

      do while (iz_lo > -nzgrid)
         if (bmag(ia, iz_lo - 1) >= B_c) exit
         iz_lo = iz_lo - 1
      end do
      do while (iz_hi < nzgrid)
         if (bmag(ia, iz_hi + 1) >= B_c) exit
         iz_hi = iz_hi + 1
      end do

      !> Three interior points is the working floor.  Two is the hard limit --
      !> the endpoint values of the well profiles are extrapolated from the two
      !> interior points nearest each turning point, so a well holding only one
      !> of them cannot be built at all -- but admitting two-point wells was
      !> measured and is not worth it: over the five configurations they are
      !> 0.01-0.6% of trapped orbits, and the four-node spline they produce is a
      !> poor enough interpolant that the budget closes marginally *less* well
      !> than when they take the caller's fallback instead.
      !>
      !> The two rejections that do carry weight are both correct as they stand.
      !> A well running off the end of the tube (1.7% of orbits in ITER, 14% in
      !> W7-X and TJ-II, 33% in QA) is not a complete bounce and cannot honestly
      !> be integrated.  A well spanning a single cell -- 6.7% in TJ-II, under
      !> 0.3% elsewhere -- has the fallback as its exact limit, not as an
      !> approximation: the drift at the point the particle sits, and no phase
      !> accumulated across a vanishing well.
      well_found = (iz_lo > -nzgrid) .and. (iz_hi < nzgrid) .and. (iz_hi - iz_lo >= 2)

   end subroutine find_well

   !============================================================================
   !> dB/dz at <z>, from the quadratic through B at the three grid points <i0>,
   !> <i1>, <i2>.  Used at a turning point, where the three points bracketing it
   !> are the only ones the well is guaranteed to own.
   !============================================================================
   real function dbdz_local(z, i0, i1, i2)

      use geometry, only: bmag
      use zgrid, only: zed

      implicit none

      real,    intent(in) :: z
      integer, intent(in) :: i0, i1, i2

      real    :: x0, x1, x2
      integer :: ia
      ia = 1

      x0 = zed(i0); x1 = zed(i1); x2 = zed(i2)
      dbdz_local = bmag(ia, i0) * (2.*z - x1 - x2) / ((x0 - x1) * (x0 - x2)) &
                 + bmag(ia, i1) * (2.*z - x0 - x2) / ((x1 - x0) * (x1 - x2)) &
                 + bmag(ia, i2) * (2.*z - x0 - x1) / ((x2 - x0) * (x2 - x1))

   end function dbdz_local

   !> Turning point between <iz_out> (where B >= B_c) and <iz_in> (where B < B_c),
   !> found on the interpolated B rather than by joining the two grid values with
   !> a straight line.
   !>
   !> The turning points are where the quadrature weight is largest, so their
   !> position is what limits the bounce average.  Linear interpolation places
   !> them to O(dz^2), which is coarse enough to dominate the answer for a well
   !> spanning only a few grid points.  Here B is splined over a window around the
   !> bracket, sampled on a fine sub-grid in a single call, and the crossing taken
   !> from the sub-interval that contains it -- accurate to O((dz/n_refine)^2),
   !> which is far below every other error in the scheme.
   real function turning_point(B_c, iz_out, iz_in)

      use geometry, only: bmag
      use zgrid, only: nzgrid, zed
      use splines, only: geo_spline

      implicit none

      real,    intent(in) :: B_c
      integer, intent(in) :: iz_out, iz_in

      !> Points either side of the bracket used to build the spline, and the
      !> number of samples across the bracketing interval.
      integer, parameter :: n_halo = 3, n_refine = 256

      real, dimension(:), allocatable :: z_window, B_window
      real, dimension(n_refine) :: z_fine, B_fine
      real    :: z_out, z_in, B_out, B_in
      integer :: ia, iz, iz_first, iz_last, n_window, i
      ia = 1

      z_out = zed(iz_out)
      z_in = zed(iz_in)

      ! Fall back on the straight line if the window would run off the domain
      iz_first = max(-nzgrid, min(iz_out, iz_in) - n_halo)
      iz_last = min(nzgrid, max(iz_out, iz_in) + n_halo)
      n_window = iz_last - iz_first + 1
      if (n_window < 4) then
         B_out = bmag(ia, iz_out); B_in = bmag(ia, iz_in)
         turning_point = z_in + (z_out - z_in) * (B_c - B_in) / (B_out - B_in)
         return
      end if

      allocate (z_window(n_window), B_window(n_window))
      do iz = iz_first, iz_last
         i = iz - iz_first + 1
         z_window(i) = zed(iz)
         B_window(i) = bmag(ia, iz)
      end do

      do i = 1, n_refine
         z_fine(i) = z_in + (z_out - z_in) * real(i - 1) / real(n_refine - 1)
      end do
      call geo_spline(z_window, B_window, z_fine, B_fine)
      deallocate (z_window, B_window)

      !> Walk out from the interior point to the first crossing of B_c, then place
      !> it within that sub-interval.  Starting from the inside matters: if the
      !> spline wobbles near the outer point, the crossing nearest the well is the
      !> physical one.
      turning_point = z_out
      do i = 1, n_refine - 1
         if ((B_fine(i) - B_c) * (B_fine(i + 1) - B_c) <= 0.) then
            turning_point = z_fine(i) + (z_fine(i + 1) - z_fine(i)) &
                            * (B_c - B_fine(i)) / (B_fine(i + 1) - B_fine(i))
            return
         end if
      end do

   end function turning_point

   !============================================================================
   !============ DRIFT-ORBIT PHASE BY INTEGRATION ALONG THE FIELD LINE =========
   !============================================================================
   !> Q(z) for a general field, from its defining relation
   !>
   !>     v_par b.grad Q = i kx (vMx - <vMx>_tau)
   !>
   !> i.e.  Q(z) = i kx int [vMx - <vMx>_tau] / (v_par b.grad z) dz.
   !>
   !> Subtracting the transit average is what makes Q single valued: it removes
   !> the secular part of the radial drift, which is exactly the part that does
   !> not average away over an orbit.  In a quasisymmetric field that average is
   !> zero and the integral has the closed form eval_Q_fac uses; in general it
   !> does not, and the leftover is what drives the zonal flow through
   !> F_RH_drift.
   !>
   !> The integration constant is irrelevant: Q enters only as
   !> <J0 exp(-Q)>_tau exp(+Q), so a shift Q -> Q + c cancels between the two
   !> factors.  The integral is therefore started wherever is convenient.
   !> The geometric part of the drift-orbit phase Q.
   !>
   !> Q is linear in kx, and at fixed pitch angle lambda = mu / energy it scales
   !> as sqrt(energy): the drift goes as energy and v_par as sqrt(energy), so the
   !> radial excursion per unit parallel length goes as their ratio.  The species
   !> enters only through the scalar smz, and the sign of v_par only as an
   !> overall sign.  Factoring all four out leaves a profile that depends on the
   !> field line and on lambda alone, which is why this is built once per pitch
   !> angle rather than once per (energy, mu, sigma, kx, species).
   subroutine eval_Q_profile_hat(lambda, iz_ref, Q_hat, drift_average_out)

      use geometry, only: bmag, gradpar, cvdrift0, gbdrift0, geo_surf, q_as_x, dbdzed
      use zgrid, only: nzgrid, zed
      use constants, only: pi
      use splines, only: geo_spline
      use parameters_physics, only: xdriftknob

      implicit none

      real,    intent(in)  :: lambda
      integer, intent(in)  :: iz_ref
      real, dimension(-nzgrid:), intent(out) :: Q_hat

      !> The transit-averaged drift this routine subtracted, at unit energy.  The
      !> Rosenbluth-Hinton drift flux has to use this very number: Q is built so
      !> that streaming and the radial drift cancel against each other except for
      !> the part averaged away here, so a flux formed from a drift average
      !> computed by any other quadrature leaves that cancellation incomplete and
      !> the budget does not close.
      real, intent(out), optional :: drift_average_out

      !> Nodes in the angle variable.  What is integrated there is smooth, so the
      !> uniform trapezoidal rule on it is second order and this is well past
      !> converged for any well a stella grid resolves.
      integer, parameter :: n_theta = 512

      real, dimension(-nzgrid:nzgrid) :: integrand
      real    :: drift_average, vpa2, vperp2, drift, dz, drift_norm, B_c
      real    :: z_l, z_r, mid, half, dtheta, arg
      integer :: ia, iz, iz_lo, iz_hi, n_well, i, k, n_target
      logical :: well_found, trapped
      real, dimension(:), allocatable :: z_well, g_well, num_well, den_well
      real, dimension(:), allocatable :: theta_target, Q_target
      real, dimension(n_theta) :: theta, z_theta, g_theta, num_theta, den_theta, Psi
      real :: sum_num, sum_den

      ia = 1

      !> cvdrift0 and gbdrift0 carry the geometry of the radial drift but not its
      !> normalisation.  time_advance builds the drift coefficient as
      !> fac * (cvdrift0 vpa^2 + gbdrift0 vperp^2 / 2), with
      !> fac = -xdriftknob * 0.5 * code_dt * tz_psi0, divided by shat unless
      !> q_as_x.  The time step belongs to the time advance, not to the orbit,
      !> but the knob does: it scales the drift the equations are actually
      !> solving, so a run with xdriftknob = 0 has no drift and hence no
      !> drift-orbit phase.  tz does not belong here -- the phase needs
      !> v_drift / v_par,
      !> and tz / stm is precisely smz, which the caller applies -- carrying tz
      !> here as well would count the same factor twice.
      drift_norm = 0.5 * xdriftknob
      if (.not. q_as_x) drift_norm = drift_norm / geo_surf%shat

      Q_hat = 0.
      drift_average = 0.

      !> Is this particle trapped, and if so, in which well?
      trapped = .false.
      if (lambda > epsilon(0.)) then
         B_c = 1. / (2.*lambda)
         if (B_c < maxval(bmag(ia, :))) then
            call find_well(B_c, iz_ref, iz_lo, iz_hi, trapped)
         end if
      end if

      !> A trapped particle whose well find_well could not close, because it runs
      !> off the end of the simulated field line.  Q is left at zero: the particle
      !> must not fall through to the passing branch, which anchors Q at the end
      !> of the line rather than at a turning point, and that is not invariant
      !> under the two-sign average the caller applies to trapped particles.
      !>
      !> The drift average is a different matter.  Zero is not its limit -- these
      !> orbits are counted as trapped by the flux split, so reporting nothing for
      !> them dilutes the trapped drive by their share of it, which is 14% of
      !> trapped weight in W7-X and 27% in TJ-II.  What can be had is the average
      !> over the part of the orbit that does lie inside the domain, which
      !> find_well returns in iz_lo and iz_hi whether or not it closed the well.
      !> That is the honest estimate at both extremes: for a well spanning barely
      !> a cell it reduces to the local value of the drift, and for a wide well
      !> truncated by the end of the tube it is the average over what was
      !> simulated.
      if (lambda > epsilon(0.)) then
         if (1. / (2.*lambda) < maxval(bmag(ia, :)) .and. .not. trapped) then
            if (present(drift_average_out)) then
               vperp2 = 2.*lambda*bmag(ia, iz_ref)
               vpa2 = 1. - vperp2
               drift_average_out = cvdrift0(ia, iz_ref) * vpa2 + gbdrift0(ia, iz_ref) * 0.5 * vperp2
            end if
            return
         end if
      end if

      !> A passing particle keeps |v_par| away from zero, so the integrand is
      !> bounded and a plain cumulative trapezoid on the z grid is second order.
      if (.not. trapped) then
         !> The secular part of the radial drift, here over the whole field line,
         !> at unit energy so that it is the same pure function of lambda as the
         !> rest of the integrand.  A trapped particle instead needs the average
         !> over its own well, and the trapped branch builds that on its own
         !> quadrature nodes, so this is deliberately not done for both.
         call eval_drift_transit_average(1., lambda, iz_ref, drift_average, well_found)
         do iz = -nzgrid, nzgrid
            vperp2 = 2.*lambda*bmag(ia, iz)
            vpa2 = 1. - vperp2
            if (vpa2 <= epsilon(0.)) then
               integrand(iz) = 0.                ! forbidden region
            else
               drift = cvdrift0(ia, iz) * vpa2 + gbdrift0(ia, iz) * 0.5 * vperp2
               integrand(iz) = drift_norm * (drift - drift_average) / (sqrt(vpa2) * gradpar(iz))
            end if
         end do
         !> Only differences of Q along an orbit are physical -- the transit
         !> average multiplies exp(-Q) by exp(+Q) at the same pitch angle -- so
         !> the constant of integration is free, and is fixed here by starting
         !> from zero at the left end of the field line.
         do iz = -nzgrid + 1, nzgrid
            dz = zed(iz) - zed(iz - 1)
            Q_hat(iz) = Q_hat(iz - 1) + 0.5 * (integrand(iz) + integrand(iz - 1)) * dz
         end do
         if (present(drift_average_out)) drift_average_out = drift_average
         return
      end if

      !> Trapped.  Now 1/|v_par| has an inverse-square-root singularity at each
      !> turning point: Q stays finite, but a trapezoid on the z grid converges
      !> only as sqrt(dz), which at usable resolution is a several-per-cent error
      !> in Q and so in everything built from it.  Writing B_c - B as
      !> (z - z_l)(z_r - z) g(z) with g smooth and substituting
      !> z = mid + half cos(theta) gives
      !>
      !>    sqrt(v_par^2) = sqrt(2 lambda) half sin(theta) sqrt(g),
      !>    dz            = -half sin(theta) dtheta,
      !>
      !> so the sin(theta) cancels identically and what is left to integrate,
      !> num / (sqrt(2 lambda) sqrt(g)), is smooth in theta.  This is the same
      !> decomposition the bounce integrals use, applied cumulatively.
      z_l = turning_point(B_c, iz_lo - 1, iz_lo)
      z_r = turning_point(B_c, iz_hi + 1, iz_hi)
      mid = 0.5 * (z_l + z_r)
      half = 0.5 * (z_r - z_l)

      n_well = (iz_hi - iz_lo + 1) + 2
      allocate (z_well(n_well), g_well(n_well), num_well(n_well), den_well(n_well))

      z_well(1) = z_l
      z_well(n_well) = z_r
      do iz = iz_lo, iz_hi
         i = iz - iz_lo + 2
         z_well(i) = zed(iz)
         g_well(i) = (B_c - bmag(ia, iz)) / ((z_well(i) - z_l) * (z_r - z_well(i)))
         vperp2 = 2.*lambda*bmag(ia, iz)
         vpa2 = 1. - vperp2
         drift = cvdrift0(ia, iz) * vpa2 + gbdrift0(ia, iz) * 0.5 * vperp2
         !> Numerator and the plain bounce-time weight kept apart, so that the
         !> drift average subtracted below can be formed on these same nodes.
         num_well(i) = drift_norm * drift / gradpar(iz)
         den_well(i) = drift_norm / gradpar(iz)
      end do

      !> g at the turning points is the limit |dB/dz| / (z_r - z_l); the numerator
      !> is smooth there and is extrapolated from the two nearest interior points,
      !> since copying the neighbour would be a first-order error sitting exactly
      !> where the weight is largest.
      !> g at the turning points, where its own definition is 0/0 and its limit
      !> is |dB/dz| / (z_r - z_l).  The difference quotient across the bracketing
      !> cell returns dB/dz at that cell's midpoint rather than at the turning
      !> point, which is first order and sits exactly where 1/sqrt(g) weights the
      !> integrand most heavily -- it was the leading error of the whole scheme.
      !> A quadratic through B on the bracket is second order and, just as
      !> importantly, stays local: it never reads a value from outside the well,
      !> so it returns the same number however many poloidal turns the flux tube
      !> spans.  Neither splining g nor reading the grid's own dbdzed can promise
      !> that, since the turning points of the barely trapped sit at the maximum
      !> of B, which for a single-turn tube is the end of the domain.
      g_well(1) = abs(dbdz_local(z_l, iz_lo - 1, iz_lo, iz_lo + 1)) / (z_r - z_l)
      g_well(n_well) = abs(dbdz_local(z_r, iz_hi - 1, iz_hi, iz_hi + 1)) / (z_r - z_l)
      num_well(1) = extrapolate(z_well(1), z_well(2), z_well(3), num_well(2), num_well(3))
      num_well(n_well) = extrapolate(z_well(n_well), z_well(n_well - 1), z_well(n_well - 2), &
                                     num_well(n_well - 1), num_well(n_well - 2))
      den_well(1) = extrapolate(z_well(1), z_well(2), z_well(3), den_well(2), den_well(3))
      den_well(n_well) = extrapolate(z_well(n_well), z_well(n_well - 1), z_well(n_well - 2), &
                                     den_well(n_well - 1), den_well(n_well - 2))

      !> theta runs from 0 at the right turning point to pi at the left one.
      dtheta = pi / (n_theta - 1)
      do k = 1, n_theta
         theta(k) = (k - 1) * dtheta
         z_theta(k) = mid + half * cos(theta(k))
      end do
      call geo_spline(z_well, num_well, z_theta, num_theta)
      call geo_spline(z_well, den_well, z_theta, den_theta)
      !> g = (B_c - B)/((z-z_l)(z_r-z)) is positive throughout the well by
      !> construction, but its cubic spline is not.  Where the well contains
      !> interior maxima of B lying just below B_c -- the ordinary situation on a
      !> stellarator field line, and where g varies over orders of magnitude --
      !> the spline overshoots and returns negative values; clamping those at
      !> tiny(0.) then gives that node a weight of 1/sqrt(tiny), some 1e153,
      !> which swamps the numerator and the bounce time alike and collapses the
      !> transit average onto a single point.  Confining the interpolant to the
      !> range of the data it was built from keeps it positive and shape
      !> preserving.  Interpolating log g instead also enforces positivity, but
      !> is unstable here: g approaches zero at those interior barriers, so its
      !> logarithm spikes and the overshoot merely moves into the exponent.
      call geo_spline(z_well, g_well, z_theta, g_theta)
      g_theta = min(max(g_theta, minval(g_well)), maxval(g_well))
      num_theta = num_theta / sqrt(g_theta)
      den_theta = den_theta / sqrt(g_theta)

      !> The bounce-averaged drift, formed on these very nodes.  It has to be
      !> this quadrature and no other: the trapped orbit visits both signs of
      !> v_par, so the caller averages the two branches, and that average is not
      !> invariant under adding a constant to Q.  Q must therefore genuinely
      !> vanish at both turning points, as the closed form v_par / B does, and it
      !> does so only if the drift average subtracted here is exactly the one
      !> this rule integrates to zero.  Taking it from a separately quadratured
      !> bounce average leaves a residue that shifts Q bodily along the orbit.
      sum_num = 0.; sum_den = 0.
      do k = 1, n_theta - 1
         sum_num = sum_num + 0.5 * dtheta * (num_theta(k) + num_theta(k + 1))
         sum_den = sum_den + 0.5 * dtheta * (den_theta(k) + den_theta(k + 1))
      end do
      if (abs(sum_den) > tiny(0.)) then
         drift_average = sum_num / sum_den
      else
         drift_average = 0.
      end if
      num_theta = num_theta - drift_average * den_theta

      !> Psi(theta) = integral of the smooth integrand from theta out to pi,
      !> which is Q measured from the left turning point.  By the choice of
      !> drift average above, Psi(1) -- the right turning point -- is zero to
      !> roundoff.
      Psi(n_theta) = 0.
      do k = n_theta - 1, 1, -1
         Psi(k) = Psi(k + 1) + 0.5 * dtheta * (num_theta(k) + num_theta(k + 1))
      end do
      Psi = Psi / sqrt(2.*lambda)

      !> Back onto the z grid.  Points outside the well are in the forbidden
      !> region, where the transit integrand vanishes and Q is never used; they
      !> are held at the nearest turning-point value so nothing downstream sees a
      !> discontinuity.
      n_target = iz_hi - iz_lo + 1
      allocate (theta_target(n_target), Q_target(n_target))
      do iz = iz_lo, iz_hi
         arg = (zed(iz) - mid) / half
         theta_target(iz - iz_lo + 1) = acos(max(-1., min(1., arg)))
      end do
      call geo_spline(theta, Psi, theta_target, Q_target)

      do iz = iz_lo, iz_hi
         Q_hat(iz) = Q_target(iz - iz_lo + 1)
      end do
      if (iz_lo > -nzgrid) Q_hat(-nzgrid:iz_lo - 1) = 0.
      if (iz_hi < nzgrid) Q_hat(iz_hi + 1:nzgrid) = Psi(1)

      if (present(drift_average_out)) drift_average_out = drift_average

      deallocate (z_well, g_well, num_well, den_well, theta_target, Q_target)

   end subroutine eval_Q_profile_hat

   !> Transit average of the radial magnetic drift over the orbit: the bounce
   !> average within the particle's own well when it is trapped, and the average
   !> along the whole field line when it is passing.
   subroutine eval_drift_transit_average(energy, mu, iz_ref, drift_average, well_found)

      use geometry, only: bmag, gradpar, cvdrift0, gbdrift0
      use zgrid, only: nzgrid, zed

      implicit none

      real,    intent(in)  :: energy, mu
      integer, intent(in)  :: iz_ref
      real,    intent(out) :: drift_average
      logical, intent(out) :: well_found

      real    :: B_c, vpa2, vperp2, drift, weight, total_weight, dz
      integer :: ia, iz
      ia = 1

      drift_average = 0.
      well_found = .false.

      !> Trapped: the bounce average over its own well, with the turning-point
      !> singularity resolved.
      if (mu > epsilon(0.)) then
         B_c = energy / (2.*mu)
         if (B_c < maxval(bmag(ia, :))) then
            call eval_bounce_averaged_drift(energy, mu, iz_ref, drift_average, well_found)
            return
         end if
      end if

      !> Passing: |v_par| never vanishes, so a plain weighted sum along the field
      !> line is enough and there is no singularity to resolve.
      total_weight = 0.
      do iz = -nzgrid, nzgrid
         vperp2 = 2.*mu*bmag(ia, iz)
         vpa2 = energy - vperp2
         if (vpa2 <= epsilon(0.)) cycle
         !> Trapezoidal weight: half a cell at each end of the line, a full cell
         !> in between.  A one-sided difference here is only first order, and
         !> since the resulting drift average is subtracted from the integrand of
         !> Q it sets the order of Q itself.
         if (iz == -nzgrid) then
            dz = 0.5 * (zed(iz + 1) - zed(iz))
         else if (iz == nzgrid) then
            dz = 0.5 * (zed(iz) - zed(iz - 1))
         else
            dz = 0.5 * (zed(iz + 1) - zed(iz - 1))
         end if
         weight = dz / (abs(gradpar(iz)) * sqrt(vpa2))
         drift = cvdrift0(ia, iz) * vpa2 + gbdrift0(ia, iz) * 0.5 * vperp2
         drift_average = drift_average + drift * weight
         total_weight = total_weight + weight
      end do
      if (total_weight > 0.) then
         drift_average = drift_average / total_weight
         well_found = .true.
      end if

   end subroutine eval_drift_transit_average

   !============================================================================
   !=============== BOUNCE-AVERAGED RADIAL MAGNETIC DRIFT ======================
   !============================================================================
   !> <vMx>_b for a trapped particle, over the well it actually occupies.
   !>
   !> This is what separates a general stellarator from a tokamak.  Axisymmetry
   !> makes the bounce-averaged radial drift vanish exactly -- canonical toroidal
   !> momentum is conserved, so there is no secular radial motion -- and
   !> quasisymmetry does the same.  In a general field it does not vanish, the
   !> transit average no longer annihilates the radial drift on its own, and the
   !> leftover drives the zonal flow through F_RH_drift.
   !>
   !> It needs no drift-orbit phase, only the geometry, so it can be evaluated in
   !> a geometry where Q is not yet available -- which is the point, since this is
   !> the quantity Q has to be generalised with.
   subroutine eval_bounce_averaged_drift(energy, mu, iz_ref, drift_average, well_found)

      use geometry, only: bmag, gradpar, cvdrift0, gbdrift0, dbdzed
      use zgrid, only: nzgrid, zed
      use constants, only: pi
      use splines, only: geo_spline

      implicit none

      real,    intent(in)  :: energy, mu
      integer, intent(in)  :: iz_ref
      real,    intent(out) :: drift_average
      logical, intent(out) :: well_found

      integer, parameter :: n_nodes = 64

      real    :: B_c, z_l, z_r, mid, half, vpa2, vperp2
      integer :: ia, iz, iz_lo, iz_hi, n_well, i
      real, dimension(:), allocatable :: z_well, g_well, drift_well, weight_well
      real, dimension(n_nodes) :: t_node, z_node, g_node, drift_node, weight_node

      ia = 1
      drift_average = 0.0
      well_found = .false.

      if (mu <= epsilon(0.)) return
      B_c = energy / (2.*mu)
      if (B_c >= maxval(bmag(ia, :))) return          ! passing, not trapped

      call find_well(B_c, iz_ref, iz_lo, iz_hi, well_found)
      if (.not. well_found) return

      !> The well's own grid points, with the two turning points appended.  Only
      !> points inside the well are used: outside it vpa^2 = energy - 2 mu B is
      !> negative, so the drift there changes character and splining through it
      !> distorts the interpolant back inside the well.
      z_l = turning_point(B_c, iz_lo - 1, iz_lo)
      z_r = turning_point(B_c, iz_hi + 1, iz_hi)
      n_well = (iz_hi - iz_lo + 1) + 2
      allocate (z_well(n_well), g_well(n_well), drift_well(n_well), weight_well(n_well))

      z_well(1) = z_l
      z_well(n_well) = z_r
      do iz = iz_lo, iz_hi
         i = iz - iz_lo + 2
         z_well(i) = zed(iz)
         !> stella's decomposition of the radial drift: the curvature piece goes
         !> with vpa^2 and the grad-B piece with vperp^2 / 2.
         vperp2 = 2.*mu*bmag(ia, iz)
         vpa2 = energy - vperp2
         drift_well(i) = cvdrift0(ia, iz) * vpa2 + gbdrift0(ia, iz) * 0.5 * vperp2
         g_well(i) = (B_c - bmag(ia, iz)) / ((z_well(i) - z_l) * (z_r - z_well(i)))
         weight_well(i) = 1.0 / abs(gradpar(iz))
      end do

      !> g at the turning points, where its own definition is 0/0 and its limit
      !> is |dB/dz| / (z_r - z_l).  The difference quotient across the bracketing
      !> cell returns dB/dz at that cell's midpoint rather than at the turning
      !> point, which is first order and sits exactly where 1/sqrt(g) weights the
      !> integrand most heavily -- it was the leading error of the whole scheme.
      !> A quadratic through B on the bracket is second order and, just as
      !> importantly, stays local: it never reads a value from outside the well,
      !> so it returns the same number however many poloidal turns the flux tube
      !> spans.  Neither splining g nor reading the grid's own dbdzed can promise
      !> that, since the turning points of the barely trapped sit at the maximum
      !> of B, which for a single-turn tube is the end of the domain.
      g_well(1) = abs(dbdz_local(z_l, iz_lo - 1, iz_lo, iz_lo + 1)) / (z_r - z_l)
      g_well(n_well) = abs(dbdz_local(z_r, iz_hi - 1, iz_hi, iz_hi + 1)) / (z_r - z_l)

      weight_well(1) = extrapolate(z_well(1), z_well(2), z_well(3), weight_well(2), weight_well(3))
      weight_well(n_well) = extrapolate(z_well(n_well), z_well(n_well - 1), z_well(n_well - 2), &
                                        weight_well(n_well - 1), weight_well(n_well - 2))
      drift_well(1) = extrapolate(z_well(1), z_well(2), z_well(3), drift_well(2), drift_well(3))
      drift_well(n_well) = extrapolate(z_well(n_well), z_well(n_well - 1), z_well(n_well - 2), &
                                       drift_well(n_well - 1), drift_well(n_well - 2))

      mid = 0.5 * (z_l + z_r)
      half = 0.5 * (z_r - z_l)
      do i = 1, n_nodes
         t_node(i) = cos((2.*i - 1.) * pi / (2.*n_nodes))
         z_node(i) = mid + half * t_node(i)
      end do

      call geo_spline(z_well, drift_well, z_node, drift_node)
      call geo_spline(z_well, weight_well, z_node, weight_node)
      call geo_spline(z_well, g_well, z_node, g_node)

      !> g = (B_c - B)/((z-z_l)(z_r-z)) is positive throughout the well by
      !> construction, but its cubic spline is not.  Where the well contains
      !> interior maxima of B lying just below B_c -- the ordinary situation on a
      !> stellarator field line, and where g varies over orders of magnitude --
      !> the spline overshoots and returns negative values; clamping those at
      !> tiny(0.) then gives that node a weight of 1/sqrt(tiny), some 1e153,
      !> which swamps the numerator and the bounce time alike and collapses the
      !> transit average onto a single point.  Confining the interpolant to the
      !> range of the data it was built from keeps it positive and shape
      !> preserving.  Interpolating log g instead also enforces positivity, but
      !> is unstable here: g approaches zero at those interior barriers, so its
      !> logarithm spikes and the overshoot merely moves into the exponent.
      g_node = min(max(g_node, minval(g_well)), maxval(g_well))
      weight_node = weight_node / sqrt(g_node)

      drift_average = sum(drift_node * weight_node) / sum(weight_node)

      deallocate (z_well, g_well, drift_well, weight_well)

   end subroutine eval_bounce_averaged_drift

   !============================================================================
   !=============== TRANSIT-AVERAGED RESPONSE AT ONE (kx, z, v) ================
   !============================================================================
   !> Evaluate <J0 exp(-iQ)>_tau / <1>_tau * exp(+/-iQ), the transit-averaged
   !> response entering the RH integrands, for both signs of vpa.
   !>
   !> The bounce-time integrand is 1/|vpa|, which does not depend on the sign of
   !> vpa, so the two calls to eval_transit_ints return the same bounce time and
   !> only one is kept.
   !> The transit average the toroidal-momentum invariant needs.  The flow part
   !> of the canonical toroidal angular momentum projects onto the field line as
   !> R v_zeta = vpar R B_zeta / B = vpar I / B, with I = R B_zeta a flux
   !> function, so the weight is vpar/B rather than vpar:
   !>     T_v = sigma * contour_integral( dl J_0 exp(-Q) / B ) ,
   !> and Pi_sigma(z) = I (T_v / tau_b) exp(Q(z)), the factor I being applied by
   !> the caller since only the geometry knows it.
   !>
   !> Weighting by v_par cancels the 1/|v_par| of the orbit measure, so unlike
   !> the density transit average this integrand is bounded at the turning
   !> points and a plain sum on the z grid is enough -- no singularity to
   !> subtract, and no need for the well quadrature.  A trapped particle is
   !> restricted to its own well; a passing one runs over the whole line.
   subroutine eval_pmom_transit_int(energy, mu, sigma, akx, iz_ref, is, Q_profile, T_v)

      use geometry, only: bmag, dl_over_b, gds22, geo_surf, q_as_x
      use species, only: spec
      use spfunc, only: j0
      use zgrid, only: nzgrid

      implicit none

      real,    intent(in)  :: energy, mu, sigma, akx
      integer, intent(in)  :: iz_ref, is
      complex, dimension(-nzgrid:), intent(in) :: Q_profile
      complex, intent(out) :: T_v

      real    :: B_c, vpa2, vperp2, kperp2, aj0
      integer :: ia, iz, iz_lo, iz_hi
      logical :: trapped, well_found

      ia = 1
      T_v = 0.

      !> Which stretch of field line this orbit occupies.
      iz_lo = -nzgrid; iz_hi = nzgrid
      trapped = .false.
      if (mu > epsilon(0.)) then
         B_c = energy / (2.*mu)
         trapped = B_c < maxval(bmag(ia, :))
         if (trapped) then
            call find_well(B_c, iz_ref, iz_lo, iz_hi, well_found)
            !> An unresolvable well contributes nothing here, consistently with
            !> the phase being left at zero for it elsewhere.
            if (.not. well_found) return
         end if
      end if

      do iz = iz_lo, iz_hi
         vpa2 = energy - 2.*mu*bmag(ia, iz)
         if (vpa2 <= epsilon(0.)) cycle
         vperp2 = 2.*mu*bmag(ia, iz)
         if (q_as_x) then
            kperp2 = akx**2 * gds22(ia, iz)
         else
            kperp2 = akx**2 * gds22(ia, iz) / (geo_surf%shat**2)
         end if
         kperp2 = max(kperp2, 0.)
         aj0 = j0(sqrt(kperp2*vperp2) * spec(is)%bess_fac * spec(is)%smz_psi0 / bmag(ia, iz))
         !> dl = (dl/B) * B, matching the convention of eval_transit_ints.
         !> The 1/B of the toroidal projection cancels the B of the arc-length
         !> measure, dl = (dl/B) * B.
         T_v = T_v + sigma * aj0 * exp(-Q_profile(iz)) * dl_over_b(ia, iz)
      end do

   end subroutine eval_pmom_transit_int


   !> The long-wavelength (order kx^2) approximation to the projection weight.
   !>
   !> Expanding W = <J0 exp(-Q)>_tau exp(Q) with Q = i kx dx and splitting into
   !> parts even and odd in sgn(v_par), Q being odd,
   !>
   !>   even(W) - 1 = -<a^2>_tau/4 - (kx^2/2)[ <dx^2>_tau + dx^2 - 2 dxbar dx ]
   !>   odd(W)      =  i kx ( dx - dxbar )
   !>
   !> so the two halves are O(kx^2) and O(kx) respectively.  Which order a given
   !> quantity needs is not uniform, and getting it wrong is easy: the leading
   !> term survives for phi_RH, but cancels for the inertia (against the 1 in
   !> 1 - J0 W), for the nonlinear flux (by quasineutrality, the gyrocentre
   !> charge density being itself O(b)) and for the collisional flux (because
   !> the collision operator conserves particles).  Having both forms available
   !> in the same run is what lets such a claim be checked rather than asserted.
   subroutine get_RH_LW_weights(energyval, muval, vpaval, akxval, iz, is, trapped, &
                                        even_LW, odd_LW, Q_hat_in)

      use species, only: spec
      use zgrid, only: nzgrid
      use geometry, only: bmag, gds22, geo_surf, q_as_x
      use constants, only: zi

      implicit none

      real,    intent(in)  :: energyval, muval, vpaval, akxval
      integer, intent(in)  :: iz, is
      logical, intent(in)  :: trapped
      complex, intent(out) :: even_LW, odd_LW
      real, dimension(-nzgrid:), intent(in), optional :: Q_hat_in

      real,    dimension(-nzgrid:nzgrid) :: Q_hat
      complex, dimension(-nzgrid:nzgrid) :: Q_profile
      real    :: dx_here, dxbar, dx2bar, tau_b, a2, kperp2, vperp2
      complex :: tmp
      integer :: ia

      ia = 1
      even_LW = 0.; odd_LW = 0.
      if (abs(akxval) <= epsilon(0.)) return

      !> The excursion profile, from the same phase the exact weight uses.
      if (use_analytic_drift_phase) then
         do ia = -nzgrid, nzgrid
            call eval_Q_fac(sign(1., vpaval) * sqrt(max(energyval - 2.*muval*bmag(1, ia), 0.)), &
                            akxval, ia, is, Q_profile(ia))
         end do
         ia = 1
      else
         if (present(Q_hat_in)) then
            Q_hat = Q_hat_in
         else if (energyval > epsilon(0.)) then
            call eval_Q_profile_hat(muval/energyval, iz, Q_hat)
         else
            Q_hat = 0.
         end if
         Q_profile = zi * akxval * sqrt(max(energyval, 0.)) * spec(is)%smz_psi0 &
                   * sign(1., vpaval) * Q_hat
      end if

      dx_here = aimag(Q_profile(iz)) / akxval

      !> Transit averages of the excursion and its square, through the same
      !> quadrature as the exact weight.
      call eval_transit_ints(energyval, muval, sign(1., vpaval), akxval, iz, is, tmp, tau_b, &
                             Q_profile, moment=1)
      dxbar = real(tmp)
      call eval_transit_ints(energyval, muval, sign(1., vpaval), akxval, iz, is, tmp, tau_b, &
                             Q_profile, moment=2)
      dx2bar = real(tmp)
      if (tau_b > 0.) then
         dxbar = dxbar / tau_b
         dx2bar = dx2bar / tau_b
      else
         return
      end if
      !> A trapped orbit samples both signs of v_par, so its mean excursion
      !> vanishes; the code's symmetrisation of the exact weight does the same.
      if (trapped) dxbar = 0.

      !> The finite-Larmor-radius piece, -<a^2>/4, with a^2 evaluated locally --
      !> its transit and field-line averages coincide under the orbit measure.
      vperp2 = 2.*muval*bmag(ia, iz)
      if (q_as_x) then
         kperp2 = akxval**2 * gds22(ia, iz)
      else
         kperp2 = akxval**2 * gds22(ia, iz) / (geo_surf%shat**2)
      end if
      kperp2 = max(kperp2, 0.)
      a2 = kperp2 * vperp2 * (spec(is)%smz_psi0 / bmag(ia, iz))**2

      even_LW = 1. - 0.25 * a2 &
                - 0.5 * akxval**2 * (dx2bar + dx_here**2 - 2.*dxbar*dx_here)
      odd_LW  = zi * akxval * (dx_here - dxbar)

   end subroutine get_RH_LW_weights


   !**********************************************************************
   !            SHORT-WAVELENGTH WEIGHTS -- STATIONARY PHASE             !
   !**********************************************************************
   !> Once the drift-orbit phase Q = kx*dx is large, the transit average
   !>
   !>     W0 = <J0 exp(-i kx dx)>_tau
   !>        = (1/tau_b) Int dz  J0 exp(-i kx dx) / (|vpa| |gradpar|)
   !>
   !> is dominated by the neighbourhoods of the stationary points of dx(z),
   !> everywhere else cancelling against itself.  Ordinary stationary phase then
   !> gives, for each stationary point z_s,
   !>
   !>     W0 = (1/tau_b) Sum_s  A(z_s) sqrt( 2 pi / (|kx| |dx''(z_s)|) )
   !>                           exp[ -i kx dx(z_s) - i sgn(kx dx''_s) pi/4 ],
   !>     A(z) = J0(z) / |vpa(z)|  times the parallel measure,
   !>
   !> so each stationary point contributes at order kx^(-1/2), and the whole
   !> weight follows from the annihilation lemma as W(z) = W0 exp(i kx dx(z)).
   !>
   !> Three things are worth saying about this.
   !>
   !> The sum is kept coherent.  The modulus quoted in the write-up drops the
   !> interference between stationary points, which is what survives averaging
   !> over a Maxwellian -- but that averaging happens afterwards, in the
   !> velocity integral, so a single particle's weight must keep its phases.
   !>
   !> The stationary points are located from dx itself rather than from B.
   !> Where dx ~ vpa/B the two coincide, since
   !>     d(vpa/B)/dz = -B' (2E - mu B) / (vpa B^2)
   !> vanishes only where B' does (2E = mu B would need vpa^2 < 0), so they are
   !> the extrema of B.  That relation needs quasisymmetry, though, and the
   !> numerically integrated drift phase does not assume it, so working from dx
   !> keeps the two paths on the same footing.
   !>
   !> Turning points are not stationary points and are not summed.  There
   !> dx ~ vpa vanishes linearly in vpa while the measure dz/|vpa| becomes
   !> regular in vpa, so the endpoint behaves as a plain Fourier integral with a
   !> linear phase and contributes at O(1/kx) -- down by kx^(-1/2) on the terms
   !> kept here.
   !> Y0 and K0, needed for the uniform (Bessel) form of the transit average
   !> near the trapped-passing boundary.  stella's spfunc carries only J0 and J1,
   !> and these are the standard Abramowitz & Stegun rational approximations
   !> (9.4.1-9.4.4 and 9.8.1-9.8.6), good to about 1e-8 -- far finer than the
   !> asymptotics they serve.  Kept private here rather than added to spfunc so
   !> that nothing outside this module is disturbed.
   real function bessel_y0(x)

      use spfunc, only: j0
      use constants, only: pi

      implicit none
      real, intent(in) :: x
      real :: y, f0, th0, r

      if (x <= 0.) then
         bessel_y0 = 0.
      else if (x <= 3.) then
         y = (x / 3.)**2
         bessel_y0 = (2./pi) * log(0.5*x) * j0(x) &
            + 0.36746691 + y*(0.60559366 + y*(-0.74350384 + y*(0.25300117 &
            + y*(-0.04261214 + y*(0.00427916 + y*(-0.00024846))))))
      else
         r = 3. / x
         f0 = 0.79788456 + r*(-0.00000077 + r*(-0.00552740 + r*(-0.00009512 &
            + r*(0.00137237 + r*(-0.00072805 + r*0.00014476)))))
         th0 = x - 0.78539816 + r*(-0.04166397 + r*(-0.00003954 + r*(0.00262573 &
            + r*(-0.00054125 + r*(-0.00029333 + r*0.00013558)))))
         bessel_y0 = f0 * sin(th0) / sqrt(x)
      end if

   end function bessel_y0

   real function bessel_k0(x)

      implicit none
      real, intent(in) :: x
      real :: y, t, i0

      if (x <= 0.) then
         bessel_k0 = 0.
      else if (x <= 2.) then
         t = x / 3.75
         y = t*t
         i0 = 1. + y*(3.5156229 + y*(3.0899424 + y*(1.2067492 + y*(0.2659732 &
            + y*(0.0360768 + y*0.0045813)))))
         y = 0.25 * x * x
         bessel_k0 = -log(0.5*x) * i0 - 0.57721566 + y*(0.42278420 + y*(0.23069756 &
            + y*(0.03488590 + y*(0.00262698 + y*(0.00010750 + y*0.0000074)))))
      else
         y = 2. / x
         bessel_k0 = exp(-x) / sqrt(x) * (1.25331414 + y*(-0.07832358 + y*(0.02189568 &
            + y*(-0.01062446 + y*(0.00587872 + y*(-0.00251540 + y*0.00053208))))))
      end if

   end function bessel_k0

   subroutine get_RH_SW_weights(energyval, muval, vpaval, akxval, iz, is, trapped, &
                                even_SW, odd_SW, Q_hat_in)

      use species, only: spec
      use zgrid, only: nzgrid, zed
      use geometry, only: bmag, gradpar, dl_over_b, gds22, geo_surf, q_as_x
      use spfunc, only: j0
      use constants, only: zi, pi

      implicit none

      real,    intent(in)  :: energyval, muval, vpaval, akxval
      integer, intent(in)  :: iz, is
      logical, intent(in)  :: trapped
      complex, intent(out) :: even_SW, odd_SW
      real, dimension(-nzgrid:), intent(in), optional :: Q_hat_in

      real,    dimension(-nzgrid:nzgrid) :: Q_hat, dx_prof
      complex, dimension(-nzgrid:nzgrid) :: Q_profile
      complex :: W0, W, tmp, hank
      real    :: tau_b, B_c, dz, sigma
      real    :: dm, d0, dp, curv, shift, dx_s, frac
      real    :: bmag_s, vpa2_s, vperp2_s, kperp2_s, aj0_s, measure_s, amp
      real    :: X_barrier, amp_u, Bc, Bm, Bpp, kappa, Xt, vpa2_in
      integer :: izb, iztp, istep, izin
      integer :: ia, izz, izm, izp, iznb, izz_start, izz_end, iz_lo, iz_hi, nstat
      logical :: well_found, is_trapped

      ia = 1
      even_SW = 0.; odd_SW = 0.
      if (abs(akxval) <= epsilon(0.)) return
      sigma = sign(1., vpaval)
      dz = zed(1) - zed(0)
      if (abs(dz) <= epsilon(0.)) return

      !> The excursion profile, from the same phase the exact weight uses.
      if (use_analytic_drift_phase) then
         do izz = -nzgrid, nzgrid
            call eval_Q_fac(sigma * sqrt(max(energyval - 2.*muval*bmag(ia, izz), 0.)), &
                            akxval, izz, is, Q_profile(izz))
         end do
      else
         if (present(Q_hat_in)) then
            Q_hat = Q_hat_in
         else if (energyval > epsilon(0.)) then
            call eval_Q_profile_hat(muval / energyval, iz, Q_hat)
         else
            Q_hat = 0.
         end if
         Q_profile = zi * akxval * sqrt(max(energyval, 0.)) * spec(is)%smz_psi0 &
                   * sigma * Q_hat
      end if
      dx_prof = aimag(Q_profile) / akxval

      !> The stretch of field line the orbit actually covers, and the bounce time
      !> in whichever normalisation that branch uses -- the trapped branch
      !> weights by 1/|gradpar| and the circulating one by B*dl_over_b, which
      !> differ by a constant, so the measure below has to match branch for
      !> branch or the ratio W0 = numerator/tau_b comes out rescaled.
      is_trapped = .false.
      if (muval > epsilon(0.)) then
         B_c = energyval / (2.*muval)
         is_trapped = B_c < maxval(bmag(ia, :))
      end if
      iz_lo = -nzgrid; iz_hi = nzgrid
      if (is_trapped) then
         call find_well(B_c, iz, iz_lo, iz_hi, well_found)
         if (.not. well_found) return
      end if

      call eval_transit_ints(energyval, muval, sigma, akxval, iz, is, tmp, tau_b, Q_profile)
      if (tau_b <= 0.) return

      !> Walk the interior of the orbit looking for sign changes of dx', and fit
      !> a parabola through each bracket to place the stationary point and to
      !> read off dx'' there.
      !> A circulating orbit runs right round the line, so its stationary points
      !> have to be looked for periodically: for an axisymmetric tube the maximum
      !> of B sits at z = +-pi, which is the domain edge, and an interior-only
      !> scan finds the minimum at z = 0 and misses it -- half the sum, and the
      !> wrong phase.  The grid duplicates the endpoint, so the distinct points
      !> are -nzgrid .. nzgrid-1 and the neighbours wrap between them.
      W0 = 0.; nstat = 0
      if (is_trapped) then
         izz_start = max(iz_lo, -nzgrid) + 1; izz_end = min(iz_hi, nzgrid) - 1
      else
         izz_start = -nzgrid; izz_end = nzgrid - 1
      end if
      do izz = izz_start, izz_end
         izm = izz - 1; izp = izz + 1
         if (.not. is_trapped) then
            if (izm < -nzgrid) izm = nzgrid - 1
            if (izp > nzgrid - 1) izp = -nzgrid
         end if
         dm = dx_prof(izm); d0 = dx_prof(izz); dp = dx_prof(izp)
         if ((dp - d0) * (d0 - dm) > 0.) cycle          ! no turning point bracketed
         curv = (dm - 2.*d0 + dp)
         if (abs(curv) <= epsilon(0.)) cycle            ! degenerate: dx'' vanishes too
         shift = 0.5 * (dm - dp) / curv                 ! vertex offset in grid cells
         if (abs(shift) > 1.) cycle
         dx_s = d0 - 0.125 * (dm - dp)**2 / curv        ! dx at the stationary point
         curv = curv / dz**2                            ! dx'' in z units

         !> Amplitude and measure interpolated to the stationary point.  These
         !> vary on the equilibrium scale while the phase turns on the much
         !> shorter scale 1/(kx dx'), which is what makes the expansion work, so
         !> linear interpolation is enough for them.
         frac = abs(shift)
         iznb = izp
         if (shift < 0.) iznb = izm
         bmag_s = (1. - frac) * bmag(ia, izz) + frac * bmag(ia, iznb)
         measure_s = (1. - frac) * parallel_measure(izz, is_trapped) &
                   + frac * parallel_measure(iznb, is_trapped)
         kperp2_s = (1. - frac) * gds22(ia, izz) + frac * gds22(ia, iznb)

         vpa2_s = energyval - 2.*muval*bmag_s
         if (vpa2_s <= epsilon(0.)) cycle               ! stationary point outside the orbit
         vperp2_s = 2.*muval*bmag_s
         if (.not. q_as_x) kperp2_s = kperp2_s / (geo_surf%shat**2)
         kperp2_s = max(akxval**2 * kperp2_s, 0.)
         aj0_s = j0(sqrt(kperp2_s * vperp2_s) * spec(is)%bess_fac * spec(is)%smz_psi0 / bmag_s)

         amp = aj0_s * measure_s / sqrt(vpa2_s)

         !> A stationary point where |dx| is a minimum is a maximum of B, and
         !> there the turning point sits nearby -- at the separatrix it sits on
         !> top of it.  The Gaussian form then diverges as 1/sqrt(|vpa|), because
         !> the amplitude blows up as 1/|vpa| faster than sqrt(1/|dx''|) can
         !> shrink.  Use the uniform Bessel form there instead, which is finite.
         !> A stationary point where |dx| is a maximum is a minimum of B: |vpa| is
         !> largest there, no turning point is near, and the Gaussian is right.
         if (dx_s * curv > 0.) then
            X_barrier = abs(akxval * dx_s)
            !> a = mu |B''| of the local model, written through the relation
            !> |dx''| = a |dx| / vpa^2 so that only quantities already in hand
            !> are used.  amp_u = J0 * measure / sqrt(a).
            amp_u = aj0_s * measure_s / sqrt(abs(curv) * vpa2_s / abs(dx_s))
            !> The whole term conjugates when the phase runs the other way, the
            !> -i pi prefactor with it -- conjugating only H0 leaves the sign of
            !> the term wrong.
            hank = -zi * pi * amp_u * cmplx(j0(X_barrier), -bessel_y0(X_barrier))
            if (akxval * dx_s < 0.) hank = conjg(hank)
            W0 = W0 + hank
         else
            W0 = W0 + amp * sqrt(2.*pi / (abs(akxval) * abs(curv))) &
                    * exp(-zi * akxval * dx_s - zi * sign(1., akxval * curv) * 0.25 * pi)
         end if
         nstat = nstat + 1
      end do

      !> The trapped side of the same coalescence.  A trapped particle never
      !> reaches the barrier top, so that stationary point is off its orbit and
      !> the sum above misses it entirely; what the orbit has instead is a
      !> turning point sitting just below the barrier.  The same local model,
      !>     vpa^2 = Delta + a zeta^2,   Delta = 2(E - mu B_m) < 0,
      !> integrated from the turning point inwards gives K0 in place of the
      !> Hankel function, which is the analytic continuation of it through
      !> Delta = 0.  It is exponentially small away from the separatrix -- which
      !> is why an endpoint contribution could be dropped there -- and grows
      !> logarithmically as the separatrix is approached, exactly where the
      !> Gaussian form was diverging on the other side.
      !>
      !> Only the real part is needed: the imaginary part of this endpoint is
      !> odd in sigma and the trapped symmetrisation below discards it.
      if (is_trapped) then
         Bc = energyval / (2.*muval)
         do istep = 1, 2
            !> Walk outward from each end of the well to the barrier that closes
            !> it, i.e. the first maximum of B beyond the turning point.
            if (istep == 1) then
               iztp = iz_lo; izb = iz_lo
               do while (izb > -nzgrid + 1)
                  if (bmag(ia, izb) >= bmag(ia, izb - 1) .and. bmag(ia, izb) >= bmag(ia, izb + 1)) exit
                  izb = izb - 1
               end do
               izin = min(iztp + 1, nzgrid)
            else
               iztp = iz_hi; izb = iz_hi
               do while (izb < nzgrid - 1)
                  if (bmag(ia, izb) >= bmag(ia, izb - 1) .and. bmag(ia, izb) >= bmag(ia, izb + 1)) exit
                  izb = izb + 1
               end do
               izin = max(iztp - 1, -nzgrid)
            end if
            if (izb <= -nzgrid .or. izb >= nzgrid) cycle
            Bm = bmag(ia, izb)
            if (Bm <= Bc) cycle                       ! no barrier above the turning point
            Bpp = (bmag(ia, izb + 1) - 2.*bmag(ia, izb) + bmag(ia, izb - 1)) / dz**2
            if (abs(Bpp) <= epsilon(0.)) cycle
            !> kappa = |dx| / |vpa| of the local model, read just inside the
            !> turning point, where B is still within (B_m - B_c) of the barrier
            !> and the ratio is therefore the barrier's own.
            vpa2_in = energyval - 2.*muval*bmag(ia, izin)
            if (vpa2_in <= epsilon(0.)) cycle
            kappa = abs(dx_prof(izin)) / sqrt(vpa2_in)
            Xt = abs(akxval) * kappa * sqrt(2.*muval*(Bm - Bc))
            W0 = W0 + j0(sqrt(max(akxval**2 * gds22(ia, izin) / merge(1., geo_surf%shat**2, q_as_x), 0.) &
                              * 2.*muval*bmag(ia, izin)) * spec(is)%bess_fac * spec(is)%smz_psi0 / bmag(ia, izin)) &
                    * parallel_measure(izin, is_trapped) / sqrt(muval * abs(Bpp)) * bessel_k0(Xt)
            nstat = nstat + 1
         end do
      end if

      !> No stationary point on this orbit: the transit average is then smaller
      !> than anything this expansion describes, and zero is the right answer to
      !> the order kept.
      if (nstat == 0) return

      W0 = W0 / tau_b

      !> A trapped particle covers both signs of vpa within one bounce, so the
      !> exact weight averages the +vpa and -vpa transit integrals, leaving W0
      !> even in sigma.  Do the same here, or the trapped population carries a
      !> spurious odd part.
      if (trapped) W0 = cmplx(real(W0), 0.)
      W = W0 * exp(zi * akxval * dx_prof(iz))

      !> dx flips sign with sigma while the measure does not, so W(-sigma) is the
      !> conjugate of W(sigma): the even part is its real part and the odd part
      !> its imaginary one.
      even_SW = real(W)
      odd_SW = zi * aimag(W)

   contains

      !> The parallel measure per unit z, in the normalisation of whichever
      !> branch of eval_transit_ints will supply tau_b.
      real function parallel_measure(izloc, in_well)
         integer, intent(in) :: izloc
         logical, intent(in) :: in_well
         integer :: izuse
         if (in_well) then
            parallel_measure = 1.0 / max(abs(gradpar(izloc)), epsilon(0.))
         else
            !> dl_over_b drops the duplicated endpoint, so read the other end of
            !> the line for it -- the same physical point.
            izuse = izloc
            if (izloc == nzgrid) izuse = -nzgrid
            parallel_measure = bmag(ia, izuse) * dl_over_b(ia, izuse) / abs(dz)
         end if
      end function parallel_measure

   end subroutine get_RH_SW_weights

   subroutine get_RH_transit_integrands(energyval, muval, vpaval, akxval, iz, is, trapped, &
                                        integrand_pls, integrand_min, Q_hat_in, integrand_v)

      use geometry, only: bmag
      use species, only: spec
      use zgrid, only: nzgrid
      use constants, only: zi

      implicit none

      real,    intent(in)  :: energyval, muval, vpaval, akxval
      integer, intent(in)  :: iz, is
      logical, intent(in)  :: trapped
      complex, intent(out) :: integrand_pls, integrand_min
      real, dimension(-nzgrid:), intent(in), optional :: Q_hat_in
      !> The parallel-flow weight V_sigma on the branch matching sgn(vpa).
      complex, intent(out), optional :: integrand_v

      real    :: transit_int_tau_b
      complex :: transit_int_eiQJ0_pls, transit_int_eiQJ0_min
      complex :: Q_fac, tmp

      real,    dimension(-nzgrid:nzgrid) :: Q_hat
      complex, dimension(-nzgrid:nzgrid) :: Q_profile
      complex :: T_v_pls, T_v_min
      integer :: iz_q

      !> The drift-orbit phase, either from the axisymmetric closed form or
      !> integrated along the field line.  It has to be the same Q in the transit
      !> average and in the exp(+/-Q) that multiplies it below: mixing the two
      !> forms leaves a spurious net phase and changes the answer outright.
      if (use_analytic_drift_phase) then
         call eval_Q_fac(vpaval, akxval, iz, is, Q_fac)
         if (present(integrand_v)) then
            do iz_q = -nzgrid, nzgrid
               call eval_Q_fac(sign(1., vpaval) * sqrt(max(energyval - 2.*muval*bmag(1, iz_q), 0.)), &
                               akxval, iz_q, is, Q_profile(iz_q))
            end do
         end if
         call eval_transit_ints(energyval, muval, sign(1., vpaval), akxval, iz, is, transit_int_eiQJ0_pls, transit_int_tau_b)
         call eval_transit_ints(energyval, muval, sign(1.,-vpaval), akxval, iz, is, transit_int_eiQJ0_min, transit_int_tau_b)
      else
         if (present(Q_hat_in)) then
            Q_hat = Q_hat_in
         else if (energyval > epsilon(0.)) then
            call eval_Q_profile_hat(muval/energyval, iz, Q_hat)
         else
            Q_hat = 0.
         end if
         !> Restore the kx, energy, species and sign-of-v_par dependence that
         !> eval_Q_profile_hat factors out.
         Q_profile = zi * akxval * sqrt(max(energyval, 0.)) * spec(is)%smz_psi0 * sign(1., vpaval) * Q_hat
         Q_fac = Q_profile(iz)
         call eval_transit_ints(energyval, muval, sign(1., vpaval), akxval, iz, is, transit_int_eiQJ0_pls, transit_int_tau_b, &
                                Q_profile)
         !> Q is odd in the sign of v_par, so the reversed branch reuses the same
         !> profile negated rather than integrating it again.
         call eval_transit_ints(energyval, muval, sign(1.,-vpaval), akxval, iz, is, transit_int_eiQJ0_min, transit_int_tau_b, &
                                -Q_profile)
      end if

      !> A particle whose integrand vanishes at every z on the grid (it is in the
      !> forbidden region everywhere) has zero bounce time and contributes
      !> nothing.  Guard the division rather than producing a NaN; this is
      !> reachable when vpa = 0 at a maximum of B.
      if (transit_int_tau_b <= 0.) then
         integrand_pls = 0.
         integrand_min = 0.
         if (present(integrand_v)) integrand_v = 0.
         return
      end if

      !> The parallel-flow weight.  Same structure as the density one -- an
      !> orbit invariant times exp(Q(z)) -- with the transit average taken of
      !> vpa*J0*exp(-Q) rather than J0*exp(-Q).
      if (present(integrand_v)) then
         call eval_pmom_transit_int(energyval, muval, sign(1., vpaval), akxval, iz, is, &
                                     Q_profile, T_v_pls)
         call eval_pmom_transit_int(energyval, muval, sign(1., -vpaval), akxval, iz, is, &
                                     -Q_profile, T_v_min)
         !> A trapped particle traverses both legs within one bounce, so its
         !> bounce average is the mean of the two -- the same symmetrisation the
         !> density weight gets, and what makes the flow weight vanish for a
         !> trapped orbit at long wavelength.
         if (trapped) then
            integrand_v = 0.5 * (T_v_pls + T_v_min) / transit_int_tau_b * exp(Q_fac)
         else
            integrand_v = T_v_pls / transit_int_tau_b * exp(Q_fac)
         end if
      end if

      !> A trapped particle traverses both signs of vpa within one bounce, so its
      !> transit average is the mean of the +vpa and -vpa averages, and is
      !> therefore even in vpa.
      if (trapped) then
         tmp = 0.5*(transit_int_eiQJ0_pls + transit_int_eiQJ0_min)
         transit_int_eiQJ0_pls = tmp
         transit_int_eiQJ0_min = tmp
      end if

      ! Evaluate integrands in the vpa-mu integral
      integrand_pls = transit_int_eiQJ0_pls/transit_int_tau_b * exp( Q_fac)
      integrand_min = transit_int_eiQJ0_min/transit_int_tau_b * exp(-Q_fac)

   end subroutine get_RH_transit_integrands


   !==============================================

   ! Evaluate RH transit averages
   !> Transit (passing) or bounce (trapped) integrals of exp(-Q) J0 and of unity.
   !>
   !> A passing particle samples the whole domain and |vpa| never vanishes, so the
   !> plain sum over the z grid is fine.  A trapped particle is a different
   !> problem on two counts: it is confined to one well, and dl/|vpa| has an
   !> integrable inverse-square-root singularity at each of its turning points.
   !> Summing that on the z grid converges only as sqrt(dz), which is easily the
   !> largest error in the transit average.  For trapped particles this therefore
   !> hands off to the well-resolved quadrature below, falling back to the plain
   !> sum only if no complete well can be found -- which happens when the well
   !> runs off the end of the simulated field line.
   subroutine eval_transit_ints(energy, mu, sigma, akx, iz_ref, is, transit_int_eiQJ0, bounce_time, Q_profile, moment)

      use geometry, only: bmag, dl_over_b
      use zgrid, only: nzgrid

      implicit none

      real,    intent(in)  :: energy, mu, sigma, akx
      integer, intent(in)  :: iz_ref, is
      complex, intent(out) :: transit_int_eiQJ0
      real,    intent(out) :: bounce_time

      complex, dimension(-nzgrid:), intent(in), optional :: Q_profile
      integer, intent(in), optional :: moment

      complex, dimension(-nzgrid:nzgrid) :: integrand_eiQJ0
      complex, dimension(-nzgrid:nzgrid) :: integrand_tau_b
      real    :: B_c
      integer :: ia, iz, iz_lo, iz_hi
      logical :: trapped, well_found
      ia = 1

      !> B_c = energy / (2 mu) is where this particle turns; it is trapped if the
      !> field line reaches that value anywhere.
      trapped = .false.
      if (mu > epsilon(0.)) then
         B_c = energy / (2.*mu)
         trapped = B_c < maxval(bmag(ia, :))
      end if

      if (trapped) then
         call find_well(B_c, iz_ref, iz_lo, iz_hi, well_found)
         if (well_found) then
            if (present(Q_profile)) then
               call bounce_ints_in_well(energy, mu, sigma, akx, B_c, iz_lo, iz_hi, is, &
                                        transit_int_eiQJ0, bounce_time, Q_profile, moment)
            else
               call bounce_ints_in_well(energy, mu, sigma, akx, B_c, iz_lo, iz_hi, is, &
                                        transit_int_eiQJ0, bounce_time, moment=moment)
            end if
            return
         end if
      end if

      ! Evaluate integrands on z-grid
      do iz = -nzgrid, nzgrid

         if (present(Q_profile)) then
            call eval_transit_int_integrand_RH(energy, mu, sigma, akx, iz, is, .false., integrand_eiQJ0(iz), Q_profile(iz), moment)
            call eval_transit_int_integrand_RH(energy, mu, sigma, akx, iz, is, .true.,  integrand_tau_b(iz), Q_profile(iz))
         else
            call eval_transit_int_integrand_RH(energy, mu, sigma, akx, iz, is, .false., integrand_eiQJ0(iz))
            call eval_transit_int_integrand_RH(energy, mu, sigma, akx, iz, is, .true.,  integrand_tau_b(iz))
         end if

      end do

      ! Evaluate integrals (integrand has 1/vpa factor, need to integrate dl/vpa (...) = dl/B * B (...) )
      transit_int_eiQJ0 = sum(integrand_eiQJ0 * bmag(ia,:) * dl_over_b(ia, :))
      bounce_time       = sum(integrand_tau_b * bmag(ia,:) * dl_over_b(ia, :))

   end subroutine eval_transit_ints

   !> The smooth part of the transit-average integrand: exp(-Q) J0, without the
   !> 1/|vpa| that carries the turning-point singularity.  bounce_ints_in_well
   !> needs the two separated, because the singular factor is absorbed into the
   !> quadrature weight rather than evaluated.
   subroutine eval_transit_int_numerator(energy, mu, sigma, akx, iz, is, numerator, Q_at_z, moment)

      use geometry, only: bmag
      use species, only: spec
      use spfunc, only: j0
      use geometry, only: gds22, geo_surf, q_as_x

      implicit none

      real,    intent(in)  :: energy, mu, sigma, akx
      integer, intent(in)  :: iz, is
      complex, intent(out) :: numerator
      complex, intent(in), optional :: Q_at_z
      !> See eval_transit_int_integrand_RH: 0 or absent gives J0 exp(-Q), 1 and 2
      !> give the radial excursion and its square, which the long-wavelength
      !> expansion needs.  Without this the trapped branch silently returned the
      !> exact integrand for a moment request.
      integer, intent(in), optional :: moment

      real    :: vpa2, vpa, vperp2, kperp2, aj0x_local, dxloc
      integer :: mom
      complex :: Q_fac
      integer :: ia
      ia = 1
      mom = 0
      if (present(moment)) mom = moment

      vpa2 = energy - 2.*mu*bmag(ia, iz)
      if (vpa2 <= epsilon(0.)) then
         numerator = 0.
         return
      end if
      vpa = sigma * sqrt(vpa2)

      vperp2 = 2 * mu * bmag(ia, iz)
      if (q_as_x) then
         kperp2 = akx**2 * gds22(ia, iz)
      else
         kperp2 = akx**2 * gds22(ia, iz) / (geo_surf%shat**2)
      end if
      kperp2 = max(kperp2, 0.)
      aj0x_local = j0(sqrt(kperp2 * vperp2) * spec(is)%bess_fac * spec(is)%smz_psi0 / bmag(ia, iz))

      if (present(Q_at_z)) then
         Q_fac = Q_at_z
      else
         call eval_Q_fac(vpa, akx, iz, is, Q_fac)
      end if

      if (mom == 0) then
      numerator = exp(-Q_fac) * aj0x_local
      else
         dxloc = 0.
         if (abs(akx) > epsilon(0.)) dxloc = aimag(Q_fac) / akx
         if (mom == 1) then
            numerator = dxloc
         else
            numerator = dxloc**2
         end if
      end if

   end subroutine eval_transit_int_numerator

   !> Linear extrapolation of f to <z>, from its values at z1 and z2.
   pure real function extrapolate(z, z1, z2, f1, f2)

      implicit none

      real, intent(in) :: z, z1, z2, f1, f2

      extrapolate = f1 + (z - z1) * (f2 - f1) / (z2 - z1)

   end function extrapolate

   !> Bounce integrals over a single well, with the turning-point singularity
   !> removed analytically rather than integrated through.
   !>
   !> Writing B_c - B(z) = (z - z_l)(z_r - z) g(z), with z_l and z_r the turning
   !> points and g smooth and positive inside the well, and mapping the well onto
   !> t in [-1, 1], the singular factor becomes exactly the Gauss-Chebyshev
   !> weight:
   !>
   !>     int F(z) dz / sqrt(B_c - B)  =  int F(z(t)) / sqrt(g(z(t))) dt/sqrt(1-t^2)
   !>
   !> so Gauss-Chebyshev nodes integrate what is left, which is smooth, to
   !> spectral accuracy.  The smooth part of the integrand is splined from its
   !> grid values; g is taken to its analytic limit at the two turning points,
   !> where the definition above is 0/0.
   subroutine bounce_ints_in_well(energy, mu, sigma, akx, B_c, iz_lo, iz_hi, is, &
                                  transit_int_eiQJ0, bounce_time, Q_profile, moment)

      use geometry, only: bmag, gradpar, dbdzed
      use zgrid, only: nzgrid, zed
      use constants, only: pi
      use splines, only: geo_spline

      implicit none

      real,    intent(in)  :: energy, mu, sigma, akx, B_c
      integer, intent(in)  :: iz_lo, iz_hi, is
      complex, intent(out) :: transit_int_eiQJ0
      real,    intent(out) :: bounce_time
      complex, dimension(-nzgrid:), intent(in), optional :: Q_profile

      !> Nodes in the Chebyshev sum.  The integrand left after the substitution is
      !> smooth, so this converges quickly; 64 is far into the converged regime
      !> for the wells a stella grid resolves.
      integer, intent(in), optional :: moment

      integer, parameter :: n_nodes = 64

      real    :: z_l, z_r, mid, half
      integer :: ia, iz, n_well, i
      real,    dimension(:), allocatable :: z_well, g_well, weight_well
      complex, dimension(:), allocatable :: numerator_well
      real,    dimension(n_nodes) :: t_node, z_node, g_node, weight_node
      complex :: value
      real    :: real_part, imag_part
      real,    dimension(:), allocatable :: tmp_real, tmp_imag
      real,    dimension(n_nodes) :: node_real, node_imag, node_weight

      ia = 1

      ! Turning points, and the well's grid points with them appended at each end
      z_l = turning_point(B_c, iz_lo - 1, iz_lo)
      z_r = turning_point(B_c, iz_hi + 1, iz_hi)
      n_well = (iz_hi - iz_lo + 1) + 2
      allocate (z_well(n_well), g_well(n_well), weight_well(n_well), numerator_well(n_well))

      z_well(1) = z_l
      z_well(n_well) = z_r
      do iz = iz_lo, iz_hi
         z_well(iz - iz_lo + 2) = zed(iz)
      end do

      !> g = (B_c - B) / ((z - z_l)(z_r - z)) on the interior points, and its
      !> limit |dB/dz| / (z_r - z_l) at the turning points.
      do iz = iz_lo, iz_hi
         i = iz - iz_lo + 2
         g_well(i) = (B_c - bmag(ia, iz)) / ((z_well(i) - z_l) * (z_r - z_well(i)))
         weight_well(i) = 1.0 / abs(gradpar(iz))
         if (present(Q_profile)) then
            call eval_transit_int_numerator(energy, mu, sigma, akx, iz, is, numerator_well(i), Q_profile(iz), moment)
         else
            call eval_transit_int_numerator(energy, mu, sigma, akx, iz, is, numerator_well(i), moment=moment)
         end if
      end do

      !> g at the turning points, where its own definition is 0/0 and its limit
      !> is |dB/dz| / (z_r - z_l).  The difference quotient across the bracketing
      !> cell returns dB/dz at that cell's midpoint rather than at the turning
      !> point, which is first order and sits exactly where 1/sqrt(g) weights the
      !> integrand most heavily -- it was the leading error of the whole scheme.
      !> A quadratic through B on the bracket is second order and, just as
      !> importantly, stays local: it never reads a value from outside the well,
      !> so it returns the same number however many poloidal turns the flux tube
      !> spans.  Neither splining g nor reading the grid's own dbdzed can promise
      !> that, since the turning points of the barely trapped sit at the maximum
      !> of B, which for a single-turn tube is the end of the domain.
      g_well(1) = abs(dbdz_local(z_l, iz_lo - 1, iz_lo, iz_lo + 1)) / (z_r - z_l)
      g_well(n_well) = abs(dbdz_local(z_r, iz_hi - 1, iz_hi, iz_hi + 1)) / (z_r - z_l)

      !> The smooth quantities at the turning points, by linear extrapolation from
      !> the two nearest interior points.  Copying the neighbour instead is an
      !> O(dz) error, and since it sits right where the quadrature weight is
      !> largest it dominates everything else -- it held the whole scheme to
      !> first order.
      weight_well(1) = extrapolate(z_well(1), z_well(2), z_well(3), weight_well(2), weight_well(3))
      weight_well(n_well) = extrapolate(z_well(n_well), z_well(n_well - 1), z_well(n_well - 2), &
                                        weight_well(n_well - 1), weight_well(n_well - 2))
      numerator_well(1) = cmplx( &
         extrapolate(z_well(1), z_well(2), z_well(3), real(numerator_well(2)), real(numerator_well(3))), &
         extrapolate(z_well(1), z_well(2), z_well(3), aimag(numerator_well(2)), aimag(numerator_well(3))))
      numerator_well(n_well) = cmplx( &
         extrapolate(z_well(n_well), z_well(n_well - 1), z_well(n_well - 2), &
                     real(numerator_well(n_well - 1)), real(numerator_well(n_well - 2))), &
         extrapolate(z_well(n_well), z_well(n_well - 1), z_well(n_well - 2), &
                     aimag(numerator_well(n_well - 1)), aimag(numerator_well(n_well - 2))))

      ! Gauss-Chebyshev nodes mapped onto the well
      mid = 0.5 * (z_l + z_r)
      half = 0.5 * (z_r - z_l)
      do i = 1, n_nodes
         t_node(i) = cos((2.*i - 1.) * pi / (2.*n_nodes))
         z_node(i) = mid + half * t_node(i)
      end do

      allocate (tmp_real(n_well), tmp_imag(n_well))
      tmp_real = real(numerator_well) * weight_well
      tmp_imag = aimag(numerator_well) * weight_well
      call geo_spline(z_well, tmp_real, z_node, node_real)
      call geo_spline(z_well, tmp_imag, z_node, node_imag)
      call geo_spline(z_well, g_well, z_node, g_node)
      deallocate (tmp_real, tmp_imag)

      !> g = (B_c - B)/((z-z_l)(z_r-z)) is positive throughout the well by
      !> construction, but its cubic spline is not.  Where the well contains
      !> interior maxima of B lying just below B_c -- the ordinary situation on a
      !> stellarator field line, and where g varies over orders of magnitude --
      !> the spline overshoots and returns negative values; clamping those at
      !> tiny(0.) then gives that node a weight of 1/sqrt(tiny), some 1e153,
      !> which swamps the numerator and the bounce time alike and collapses the
      !> transit average onto a single point.  Confining the interpolant to the
      !> range of the data it was built from keeps it positive and shape
      !> preserving.  Interpolating log g instead also enforces positivity, but
      !> is unstable here: g approaches zero at those interior barriers, so its
      !> logarithm spikes and the overshoot merely moves into the exponent.
      g_node = min(max(g_node, minval(g_well)), maxval(g_well))
      node_weight = 1.0 / sqrt(g_node)

      !> The Chebyshev rule carries a common factor pi/n_nodes and a common
      !> 1/sqrt(2 mu) from |vpa|; both cancel in the ratio the caller forms, but
      !> are kept so the two returned integrals are individually the integrals
      !> they claim to be.
      real_part = sum(node_real * node_weight) * pi / n_nodes / sqrt(2.*mu)
      imag_part = sum(node_imag * node_weight) * pi / n_nodes / sqrt(2.*mu)
      transit_int_eiQJ0 = cmplx(real_part, imag_part)

      call geo_spline(z_well, weight_well, z_node, node_weight)
      node_weight = node_weight / sqrt(g_node)
      bounce_time = sum(node_weight) * pi / n_nodes / sqrt(2.*mu)

      deallocate (z_well, g_well, weight_well, numerator_well)

   end subroutine bounce_ints_in_well


   ! Evaluate integrand in RH transit average
   !> <moment> selects what is transit-averaged: absent or 0 gives the
   !> Rosenbluth-Hinton integrand J_0 exp(-Q), 1 gives the radial excursion
   !> delta_x = Q/(i kx), and 2 gives delta_x^2.  The last two are what the
   !> long-wavelength expansion needs, and taking them through the same
   !> quadrature is the point: the turning-point treatment that makes the
   !> exact weight accurate makes the asymptotic one accurate too.
   subroutine eval_transit_int_integrand_RH(energy, mu, sigma, akx, iz, is, bounce_time_bool, transit_avg_integrand, Q_at_z, moment)

      use geometry, only: bmag
      use species, only: spec
      use spfunc, only: j0
      use geometry, only: gds22, geo_surf, q_as_x

      implicit none

      real,    intent(in)  :: energy, mu, sigma, akx ! energy=vpa^2+vperp^2, mu=vperp^2/(2B)
      integer, intent(in)  :: iz, is
      logical, intent(in)  :: bounce_time_bool ! if true, evaluate integrand for bounce time
      complex, intent(out) :: transit_avg_integrand
      complex, intent(in), optional :: Q_at_z
      integer, intent(in), optional :: moment

      real    :: dxloc
      integer :: mom
      real    :: vpa2, vpa, vperp2, kperp2
      complex :: Q_fac, aj0x
      integer :: ia
      ia = 1
      mom = 0
      if (present(moment)) mom = moment

      !> Evaluate the integrand, zero in the forbidden region.  This is correct
      !> for a field line with many wells as it stands: which well the orbit
      !> belongs to is settled by find_well in eval_transit_ints, and all this
      !> has to do is refuse to contribute where v_par^2 would be negative.
      vpa2 = energy-2.*mu*bmag(ia,iz)
      if (vpa2 <= epsilon(0.)) then
         transit_avg_integrand = 0
      else
         ! Parallel velocity
         vpa = sigma*sqrt(vpa2)

         if (bounce_time_bool) then
            transit_avg_integrand = 1./abs(vpa)

         else
            ! Evaluate Bessel function
            vperp2 = 2*mu*bmag(ia,iz)
            !> Note this cannot simply reuse gyro_averages::aj0x: the transit
            !> average is also evaluated at <kxsmall>, which is not a grid kx.
            if (q_as_x) then
               kperp2 = akx**2 * gds22(ia,iz)
            else
               kperp2 = akx**2 * gds22(ia,iz) / (geo_surf%shat**2)
            end if
            ! gds22 can carry small negative interpolation noise; kperp2 >= 0
            kperp2 = max(kperp2, 0.)
            aj0x = j0( sqrt(kperp2*vperp2) * spec(is)%bess_fac * spec(is)%smz_psi0 / bmag(ia,iz) )

            ! Evaluate Q factor
            if (present(Q_at_z)) then
               Q_fac = Q_at_z
            else
               call eval_Q_fac(vpa, akx, iz, is, Q_fac)
            end if

            ! Integrand
            if (mom == 0) then
               transit_avg_integrand = exp(-Q_fac) * aj0x / abs(vpa)
            else
               !> Q is purely imaginary, so the excursion is its imaginary part
               !> divided by kx.  Guard the division: at kx = 0 the excursion is
               !> not resolved by the phase and the moments are not needed.
               dxloc = 0.
               if (abs(akx) > epsilon(0.)) dxloc = aimag(Q_fac) / akx
               if (mom == 1) then
                  transit_avg_integrand = dxloc / abs(vpa)
               else
                  transit_avg_integrand = dxloc**2 / abs(vpa)
               end if
            end if

         end if
      end if

   end subroutine eval_transit_int_integrand_RH

   ! Evaluate Q factor (i*kx*vmx = vpa*nabla_par(Q))
   !> Drift-orbit phase Q_s, defined so that transit-averaging annihilates the
   !> radial magnetic drift:
   !>
   !>     Q_s = i kx (v_par / Omega_s) * RH_drift_phase_fac
   !>
   !> The geometry-dependent half is <RH_drift_phase_fac>, which the geometry
   !> module builds -- in a quasisymmetric field it is (MG+NI)/(N-iota*M), and in
   !> a tokamak that reduces to the q R Btor form.  Keeping it there rather than
   !> here means this routine does not care which equilibrium it is looking at,
   !> and a geometry that learns to provide the factor needs no change to the
   !> Rosenbluth-Hinton code.  Compare diagnostics_fluxes_fluxtube, which consumes
   !> b_dot_grad_zeta_RR the same way.
   subroutine eval_Q_fac(vpa, akx, iz, is, Q_fac)

      use geometry, only: bmag, RH_drift_phase_fac
      use species, only: spec
      use constants, only: zi
      use parameters_physics, only: xdriftknob

      implicit none

      real,    intent(in)  :: vpa, akx
      integer, intent(in)  :: iz, is
      complex, intent(out) :: Q_fac

      integer :: ia
      ia = 1

      !> The normalisation is checked against the phase integrated along the
      !> field line, which is built from the geometry independently of this
      !> closed form: in a Miller tokamak the two agree on the coefficient of
      !> v_par / B to 0.06%, and the resulting RH inertia to 0.014% at nzed = 48
      !> and 0.0019% at nzed = 256.  RH_drift_phase_fac is q * btor * Rmajor /
      !> rhoc, which is q * bi / rhoc and so constant on the surface; smz = sqrt(m
      !> T) / Z is what converts v_drift / v_par into a radial excursion in units
      !> of rho_ref, and carrying tz here as well would count that factor twice.
      Q_fac = zi*akx * vpa/bmag(ia,iz) * spec(is)%smz_psi0 * RH_drift_phase_fac(iz) * xdriftknob

   end subroutine eval_Q_fac

end module rosenbluth_hinton
