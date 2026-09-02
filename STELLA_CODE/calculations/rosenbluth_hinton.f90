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
   public :: RH_U_parallel_fac
   public :: RH_inertia
   public :: RH_integrand_even, RH_integrand_odd

   real, dimension(:,:,:), allocatable :: RH_U_parallel_fac
   ! (-nzgrid:nzgrid, ntubes, -vmu-layout-)

   complex, dimension(:,:,:,:), allocatable :: RH_inertia
   ! (nakx, -nzgrid:nzgrid, ntubes, nspec)

   complex, dimension(:,:,:,:), allocatable :: RH_integrand_even, RH_integrand_odd
   ! (nakx, -nzgrid:nzgrid, ntubes, -vmu-layout-)

   private

   ! Debugging
   logical :: debug = .false.

   ! Has this module been initialised?
   logical :: rosenbluth_hinton_initialized = .false.

contains


!###############################################################################
!############################ INITALIZE & FINALIZE #############################
!###############################################################################

   !============================================================================
   !======================== INITALIZE THE DIAGNOSTICS =========================
   !============================================================================
   subroutine init_rosenbluth_hinton()

      use mp, only: proc0

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

      implicit none

      real :: energyval, muval
      complex :: Q_fac
      real :: transit_int_tau_b_pls, transit_int_tau_b_min
      complex :: transit_int_eiQJ0_pls, transit_int_eiQJ0_min
      complex :: integrand_tmp_pls, integrand_tmp_min
      complex :: tmp
      real :: kxsmall

      integer :: ivmu, iv, imu, is, ia, iz, it, ikx

      kxsmall = 1.d-8
      ia = 1

      !> The RH response functions are needed when they are diagnosed, when the
      !> parallel-shear term uses them (<omprimfac_RH>), or when a prescribed
      !> zonal profile is built from the RH closure.  Building them costs an
      !> O(nvmu * nz^2 * nakx) transit-average loop, so skip it otherwise.
      if (.not. rosenbluth_hinton_needed()) return

      ! Only initialize once
      if (rosenbluth_hinton_initialized) return
      rosenbluth_hinton_initialized = .true.


      ! Only debug on the first processor
      debug = debug .and. proc0

      ! Allocate the arrays for the Rosenbluth-Hinton integrand term
      allocate (RH_integrand_even(nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); RH_integrand_even = 0.
      allocate (RH_integrand_odd( nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); RH_integrand_odd  = 0.

      ! Allocate array for RH_U_parallel_fac
      allocate (RH_U_parallel_fac( -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); RH_U_parallel_fac = 0.

      ! Allocate the array for the RH_inertia
      allocate (RH_inertia(nakx, -nzgrid:nzgrid, ntubes, nspec)); RH_inertia = 0

      ! Evaluate the transit averages
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)

         do iz = -nzgrid, nzgrid

            ! Determine energy and magnetic moment
            energyval = vpa(iv)**2 + vperp2(ia,iz,imu)
            muval     = mu(imu)

            do it = 1, ntubes
               do ikx = 1, nakx

                  ! Evaluate transit averages for vpa and -vpa
                  call eval_transit_ints(energyval, muval, sign(1., vpa(iv)), akx(ikx), is, transit_int_eiQJ0_pls, transit_int_tau_b_pls)
                  call eval_transit_ints(energyval, muval, sign(1.,-vpa(iv)), akx(ikx), is, transit_int_eiQJ0_min, transit_int_tau_b_min)

                  ! For trapped particles, the result is the average of +vpa and -vpa transit averages
                  ! => transit_int_eiQJ0_{pls,min} is even in vpa for trapped particles
                  if (energyval <= 2*muval*maxval(bmag(ia,:))) then
                     tmp = 0.5*(transit_int_eiQJ0_pls + transit_int_eiQJ0_min)
                     transit_int_eiQJ0_pls = tmp
                     transit_int_eiQJ0_min = tmp
                  end if

                  ! Get Q factor
                  call eval_Q_fac(vpa(iv), akx(ikx), iz, is, Q_fac)

                  ! Evaluate integrands in vpa-mu integral
                  integrand_tmp_pls = transit_int_eiQJ0_pls/transit_int_tau_b_pls * exp( Q_fac)
                  integrand_tmp_min = transit_int_eiQJ0_min/transit_int_tau_b_pls * exp(-Q_fac)

                  ! Split into contributions that are even and odd in vpa
                  RH_integrand_even(ikx,iz,it,ivmu) = 0.5*(integrand_tmp_pls+integrand_tmp_min)
                  RH_integrand_odd( ikx,iz,it,ivmu) = 0.5*(integrand_tmp_pls-integrand_tmp_min)


               ! Evaluate RH_U_parallel_fac (same as above but tiny kx!)

               ! Evaluate transit averages for vpa and -vpa
               call eval_transit_ints(energyval, muval, sign(1., vpa(iv)), kxsmall, is, transit_int_eiQJ0_pls, transit_int_tau_b_pls)
               call eval_transit_ints(energyval, muval, sign(1.,-vpa(iv)), kxsmall, is, transit_int_eiQJ0_min, transit_int_tau_b_min)

               ! Get Q factor
               call eval_Q_fac(vpa(iv), kxsmall, iz, is, Q_fac)

               ! For trapped particles, the result is the average of +vpa and -vpa transit averages
               ! Note: it follows that transit_int_eiQJ0_{pls,min} should be even in vpa
               if (energyval <= 2*muval*maxval(bmag(ia,:))) then
                  tmp = 0.5*(transit_int_eiQJ0_pls + transit_int_eiQJ0_min)
                  transit_int_eiQJ0_pls = tmp
                  transit_int_eiQJ0_min = tmp
               end if

               ! Evaluate integrands in vpa-mu integral
               integrand_tmp_pls = transit_int_eiQJ0_pls/transit_int_tau_b_pls * exp( Q_fac)
               integrand_tmp_min = transit_int_eiQJ0_min/transit_int_tau_b_pls * exp(-Q_fac)

               ! RH_U_parallel_fac
               RH_U_parallel_fac(iz,it,ivmu) = real( (1 - 0.5*(integrand_tmp_pls-integrand_tmp_min))/(zi*kxsmall) &
                                               * spec(is)%z/spec(is)%mass )

               end do !ikx
            end do !it
         end do !iz
      end do !ivmu

      ! TODO-RN : implement for radial variation and full flux surface
      ! Calculate the RH_inertia for a flux tube simulation
      call get_RH_inertia_fluxtube()

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
      if (allocated(RH_U_parallel_fac)) deallocate (RH_U_parallel_fac)
      if (allocated(RH_inertia))        deallocate (RH_inertia)

      rosenbluth_hinton_initialized = .false.

   end subroutine finish_rosenbluth_hinton

   !============================================================================
   !=========== IS THE ROSENBLUTH-HINTON MACHINERY NEEDED AT ALL? ==============
   !============================================================================
   !> Single source of truth for the init/finish guard, so the two can never
   !> disagree and leak the (large) response arrays.
   logical function rosenbluth_hinton_needed()

      use parameters_diagnostics, only: write_RH_inertia_fluxes
      use parameters_physics, only: omprimfac_RH
      use parameters_physics, only: triangular_ZF, cos_ZF, triangular_ZF_RH

      implicit none

      rosenbluth_hinton_needed = write_RH_inertia_fluxes &
                                 .or. abs(omprimfac_RH) > epsilon(0.) &
                                 .or. ((triangular_ZF .or. cos_ZF) .and. triangular_ZF_RH)

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

                integrand_vpamu(1, :, iz, it, ivmu) = (1 - aj0x(1,:,iz,ivmu)*(RH_integrand_even(:,iz,it,ivmu)+RH_integrand_odd(:,iz,it,ivmu))) * &
                                       maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is)*maxwell_fac(is) * spec(is)%zt

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
                                        RH_fluxes_coll)

      use zgrid, only: nzgrid, ntubes
      use species, only: spec, nspec
      use vpamu_grids, only: vpa, mu, vperp2, integrate_vmu
      use vpamu_grids, only: maxwell_mu, ztmax, maxwell_fac, maxwell_vpa
      use parameters_kxky_grids, only: naky, nakx, nx
      use grids_kxky, only: aky, akx
      use calculations_kxky, only: multiply_by_rho
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use gyro_averages, only: gyro_average, gyro_average_j1
      use arrays_fields, only: phi, apar, bpar
      use parameters_numerical, only: maxwellian_normalization
      use stella_transforms, only: transform_kx2x_xfirst, transform_x2kx_xfirst
      use constants, only: zi
      use parameters_physics, only: nonlinear
      use geometry, only: exb_nonlin_fac
      use parameters_numerical, only: fphi
      use parameters_physics, only: include_apar, include_bpar
      use dissipation, only: include_collisions, collisions_implicit
      use dissipation, only: advance_collisions_explicit, advance_collisions_implicit
      use stella_time, only: code_dt

      ! Import temp arrays g1, g2 with dimensions (nky, nkx, -nzgrid:nzgrid, ntubes, -vmu-layout-)
      use arrays_dist_fn, only: integrand_even   => g0
      use arrays_dist_fn, only: integrand_odd    => g1

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
      complex, dimension(   :, -nzgrid:, :, :), intent(out) :: RH_fluxes_coll

      ! Local variables
      integer :: ivmu, iv, imu, is, ia, iz, it

      ! We only have one field line because <full_flux_surface> = .false.
      ia = 1

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
               end do
            end do
         end do

         ! Calculate <RH_fluxes>(even/odd)
         call integrate_vmu(integrand_even, spec%dens_psi0*spec%z, RH_fluxes_phi_even)
         call integrate_vmu(integrand_odd,  spec%dens_psi0*spec%z, RH_fluxes_phi_odd)

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


      ! Only compute RH collisional flux when collisions are included
      if (.not. include_collisions) then
         RH_fluxes_coll = 0.
      else

         !!!!!!!!!!!!!!!!!!!!!!!!!
         !!! coll contribution !!!
         !!!!!!!!!!!!!!!!!!!!!!!!!
         allocate (RH_fluxes_coll_tmp(naky, nakx, -nzgrid:nzgrid, ntubes, nspec)); RH_fluxes_coll_tmp = 0.

         ! Evaluate dt * collision operator (stored in integrand_even)
         if (collisions_implicit) then
            integrand_even = g
            call advance_collisions_implicit(.false., phi, apar, bpar, integrand_even)
            integrand_even = integrand_even - g
         else
            integrand_even = 0.
            call advance_collisions_explicit(g, phi, bpar, integrand_even)
         end if

         ! Remove nonzonal terms to save computational time
         integrand_even(2:,:,:,:,:) = 0.0

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

         deallocate(RH_fluxes_coll_tmp)

      endif


   end subroutine get_RH_fluxes_fluxtube
 

   !============================================================================
   !====================== GET RH_phi_I FOR THE FLUX TUBE ========================
   !============================================================================
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

      call integrate_vmu(g * RH_integrand_tmp, spec%z, RH_phi_I_tmp)
      RH_phi_I = RH_phi_I_tmp(1,:,:,:,:)

      deallocate (RH_phi_I_tmp)

   end subroutine get_RH_phi_I_fluxtube

   !==============================================
   !============== BOUNCE AVERAGES ===============
   !==============================================

   ! Evaluate RH transit averages
   subroutine eval_transit_ints(energy, mu, sigma, akx, is, transit_int_eiQJ0, bounce_time)

      use geometry, only: bmag, dl_over_b
      use zgrid, only: nzgrid
      use constants, only: zi

      implicit none

      real,    intent(in)  :: energy, mu, sigma, akx
      integer, intent(in)  :: is
      complex, intent(out) :: transit_int_eiQJ0
      real,    intent(out) :: bounce_time

      complex, dimension(-nzgrid:nzgrid) :: integrand_eiQJ0
      complex, dimension(-nzgrid:nzgrid) :: integrand_tau_b
      real    :: vpa2, vpa
      complex :: Q_fac
      integer :: ia, iz
      ia = 1

      ! Evaluate integrands on z-grid
      do iz = -nzgrid, nzgrid

         call eval_transit_int_integrand_RH(energy, mu, sigma, akx, iz, is, .false., integrand_eiQJ0(iz))
         call eval_transit_int_integrand_RH(energy, mu, sigma, akx, iz, is, .true.,  integrand_tau_b(iz))

      end do

      ! Evaluate integrals (integrand has 1/vpa factor, need to integrate dl/vpa (...) = dl/B * B (...) )
      transit_int_eiQJ0 = sum(integrand_eiQJ0 * bmag(ia,:) * dl_over_b(ia, :))
      bounce_time       = sum(integrand_tau_b * bmag(ia,:) * dl_over_b(ia, :))

   end subroutine eval_transit_ints


   ! Evaluate integrand in RH transit average
   subroutine eval_transit_int_integrand_RH(energy, mu, sigma, akx, iz, is, bounce_time_bool, transit_avg_integrand)

      use geometry, only: bmag
      use species, only: spec
      use spfunc, only: j0
      use geometry, only: gds22, geo_surf, q_as_x

      implicit none

      real,    intent(in)  :: energy, mu, sigma, akx ! energy=vpa^2+vperp^2, mu=vperp^2/(2B)
      integer, intent(in)  :: iz, is
      logical, intent(in)  :: bounce_time_bool ! if true, evaluate integrand for bounce time
      complex, intent(out) :: transit_avg_integrand

      real    :: vpa2, vpa, vperp2, kperp2
      complex :: Q_fac, aj0x
      integer :: ia
      ia = 1

      ! Evaluate integrand (=0 if in forbidden region)
      ! TODO-RN: implement for multiple wells
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
            if (q_as_x) then
               kperp2 = akx**2 * gds22(ia,iz)
            else
               kperp2 = akx**2 * gds22(ia,iz) / (geo_surf%shat**2)
            end if
            aj0x = j0( sqrt(kperp2*vperp2) * spec(is)%bess_fac * spec(is)%smz_psi0 / bmag(ia,iz) )

            ! Evaluate Q factor
            call eval_Q_fac(vpa, akx, iz, is, Q_fac)

            ! Integrand
            transit_avg_integrand = exp(-Q_fac) * aj0x / abs(vpa)

         end if
      end if

   end subroutine eval_transit_int_integrand_RH

   ! Evaluate Q factor (i*kx*vmx = vpa*nabla_par(Q))
   subroutine eval_Q_fac(vpa, akx, iz, is, Q_fac)
      ! TODO-RN : implement correctly for general geometry

      use geometry, only: bmag, geo_surf, btor, Rmajor
      use species, only: spec
      use constants, only: zi

      implicit none

      real,    intent(in)  :: vpa, akx
      integer, intent(in)  :: iz, is
      complex, intent(out) :: Q_fac

      integer :: ia
      ia = 1

      ! TODO-RN : Normalisation OK?
      Q_fac = zi*akx * vpa/bmag(ia,iz) * spec(is)%smz_psi0 &
              * geo_surf%qinp_psi0*btor(iz)*Rmajor(iz)/geo_surf%rhoc ! Note Btor*Rmajor should be constant along field-line

   end subroutine eval_Q_fac

end module rosenbluth_hinton
