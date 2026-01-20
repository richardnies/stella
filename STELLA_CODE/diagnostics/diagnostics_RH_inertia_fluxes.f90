
!###############################################################################
!############### DIAGNOSE ROSENBLUTH-HINTON INERTIA AND FLUXES #################
!###############################################################################
! 
! Routines for calculating and writing the "Rosenbluth-Hinton" inertia and fluxes,
! governing the evolution of long-wavelength stationary zonal flows.
! 
! The RH_inertia_vs_kxzts             is denoted by RH_inertia
! The RH_integrand_even_vs_kxztsvpamu is denoted by RH_integrand_even
! The RH_integrand_odd_vs_kxztsvpamu  is denoted by RH_integrand_odd
! The RH_fluxes_phi_even_vs_kykxzts       is denoted by RH_fluxes_phi_even
! The RH_fluxes_phi_odd_vs_kykxzts        is denoted by RH_fluxes_phi_odd
! The RH_fluxes_apar_even_vs_kykxzts      is denoted by RH_fluxes_apar_even
! The RH_fluxes_apar_odd_vs_kykxzts       is denoted by RH_fluxes_apar_odd
! The RH_fluxes_bpar_even_vs_kykxzts      is denoted by RH_fluxes_bpar_even
! The RH_fluxes_bpar_odd_vs_kykxzts       is denoted by RH_fluxes_bpar_odd
! 
!###############################################################################
 
module diagnostics_RH_inertia_fluxes

   implicit none
 
   public :: init_diagnostics_RH_inertia_fluxes
   public :: finish_diagnostics_RH_inertia_fluxes
   public :: eval_transit_ints
   public :: eval_transit_int_integrand
   public :: eval_Q_fac
   public :: write_RH_fluxes_to_netcdf_file
   public :: write_RH_inertia_to_netcdf_file
   public :: write_RH_phi_I_to_netcdf_file

   private 

   ! Debugging
   logical :: debug = .false.

   complex, dimension(:,:,:,:), allocatable :: RH_integrand_even, RH_integrand_odd
   ! (nakx, -nzgrid:nzgrid, ntubes, -vmu-layout-)

contains

!###############################################################################
!###################### WRITE RH_INERTIA_FLUXES ################################
!###############################################################################

   !============================================================================
   !========== CALCULATE AND WRITE RH_INTEGRANDS TO NETCDF FILE ================
   !============================================================================
   subroutine write_RH_integrands_to_netcdf_file()

      ! Dimensions
      use parameters_kxky_grids, only: naky, nakx
      use vpamu_grids, only: nvpa, nmu
      use zgrid, only: nztot, ntubes
      use species, only: nspec
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use mp, only: sum_reduce
      
      ! Flags 
      use parameters_physics, only: full_flux_surface

      ! Write to netcdf file 
      use stella_io, only: write_RH_integrands_nc
      
      ! Routines
      use mp, only: proc0
      
      ! Input file
      use parameters_diagnostics, only: write_RH_inertia_fluxes

      implicit none 

      integer :: ivmu, iv, imu, is

      ! Variables needed to write and calculate diagnostics 
      complex, dimension(:, :, :, :, :, :), allocatable :: RH_integrand_even_vs_kxztsvpamu
      complex, dimension(:, :, :, :, :, :), allocatable :: RH_integrand_odd_vs_kxztsvpamu

      !---------------------------------------------------------------------- 

      ! Only continue if the RH_inertia_fluxes have to be written
      if (.not. write_RH_inertia_fluxes) return

      ! Allocate the arrays for the RH_integrands
      allocate (RH_integrand_even_vs_kxztsvpamu(nakx, nztot, ntubes, nspec, nvpa, nmu))
      allocate (RH_integrand_odd_vs_kxztsvpamu( nakx, nztot, ntubes, nspec, nvpa, nmu))

      ! Calculate the RH inertia (kx,tube,s); RH fluxes(kx,tube,s)
      if (debug) write (*, *) 'diagnostics::diagnostics_stella::write_RH_integrands'

      ! TODO-RN : implement for radial variation and full flux surface

      ! Put the RH_integrands in form expected in stella_io
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)

         RH_integrand_even_vs_kxztsvpamu(:,:,:,is,iv,imu) = RH_integrand_even(:,:,:,ivmu)
         RH_integrand_odd_vs_kxztsvpamu( :,:,:,is,iv,imu) = RH_integrand_odd( :,:,:,ivmu)
      end do

      ! Make sure proc0 has full array 
      ! TODO-RN: might need too much memory? Average first in zed?
      call sum_reduce(RH_integrand_even_vs_kxztsvpamu, 0)
      call sum_reduce(RH_integrand_odd_vs_kxztsvpamu,  0)

      ! Write the RH_integrand to the netcdf file
      if (proc0 .and. write_RH_inertia_fluxes) then 
          call write_RH_integrands_nc(RH_integrand_even_vs_kxztsvpamu, &
                                      RH_integrand_odd_vs_kxztsvpamu)
      end if

      ! Deallocate the arrays for the RH_integrand
      deallocate (RH_integrand_even_vs_kxztsvpamu)
      deallocate (RH_integrand_odd_vs_kxztsvpamu)

   end subroutine write_RH_integrands_to_netcdf_file
 

   !============================================================================
   !========== CALCULATE AND WRITE RH_INERTIA TO NETCDF FILE ===================
   !============================================================================
   subroutine write_RH_inertia_to_netcdf_file()

      ! Dimensions
      use parameters_kxky_grids, only: naky, nakx
      use zgrid, only: nztot, ntubes
      use species, only: nspec
      
      ! Flags 
      use parameters_physics, only: full_flux_surface

      ! Write to netcdf file 
      use stella_io, only: write_RH_inertia_nc
      
      ! Routines
      use mp, only: proc0
      
      ! Input file
      use parameters_diagnostics, only: write_RH_inertia_fluxes

      implicit none 

      ! Variables needed to write and calculate diagnostics 
      complex, dimension(:, :, :, :), allocatable :: RH_inertia_vs_kxzts

      !---------------------------------------------------------------------- 

      ! Only continue if the RH_inertia_fluxes have to be written
      if (.not. write_RH_inertia_fluxes) return

      ! Allocate the arrays for the RH_inertia
      allocate (RH_inertia_vs_kxzts(nakx, nztot, ntubes, nspec))

      ! Calculate the RH inertia (kx,tube,s); RH fluxes(kx,tube,s)
      if (debug) write (*, *) 'diagnostics::diagnostics_stella::write_RH_inertia'

      ! TODO-RN : implement for radial variation and full flux surface

      ! Calculate the RH_inertia for a flux tube simulation
      if (write_RH_inertia_fluxes) then
         call get_RH_inertia_fluxtube(RH_inertia_vs_kxzts)
      end if

      ! Write the RH_inertia to the netcdf file
      if (proc0 .and. write_RH_inertia_fluxes) call write_RH_inertia_nc(RH_inertia_vs_kxzts)

      ! Deallocate the arrays for the RH_inertia
      deallocate (RH_inertia_vs_kxzts)

   end subroutine write_RH_inertia_to_netcdf_file
 
   !============================================================================
   !========== CALCULATE AND WRITE RH_FLUXES TO NETCDF FILE ====================
   !============================================================================
   subroutine write_RH_fluxes_to_netcdf_file(nout, timer)

      ! Data
      use arrays_dist_fn, only: gnew

      ! Dimensions
      use parameters_kxky_grids, only: naky, nakx
      use zgrid, only: nztot, ntubes
      use species, only: nspec
      
      ! Flags 
      use parameters_physics, only: full_flux_surface

      ! Write to netcdf file 
      use stella_io, only: write_RH_fluxes_phi_nc, write_RH_fluxes_apar_nc, write_RH_fluxes_bpar_nc
      
      ! Routines
      use job_manage, only: time_message
      use mp, only: proc0
      
      ! Input file
      use parameters_diagnostics, only: write_RH_inertia_fluxes
      
      ! Physics parameters
      use parameters_physics, only: include_apar, include_bpar

      implicit none 

      ! The pointer in the netcdf file and a timer
      real, dimension(:), intent(in out) :: timer   
      integer, intent(in) :: nout    

      ! Variables needed to write and calculate diagnostics 
      complex, dimension(:, :, :, :, :), allocatable :: RH_fluxes_phi_even_vs_kykxzts, RH_fluxes_phi_odd_vs_kykxzts
      complex, dimension(:, :, :, :, :), allocatable :: RH_fluxes_apar_even_vs_kykxzts, RH_fluxes_apar_odd_vs_kykxzts
      complex, dimension(:, :, :, :, :), allocatable :: RH_fluxes_bpar_even_vs_kykxzts, RH_fluxes_bpar_odd_vs_kykxzts

      !---------------------------------------------------------------------- 

      ! Only continue if the RH_inertia_fluxes have to be written
      if (.not. write_RH_inertia_fluxes) return  

      ! Start timer
      if (proc0) call time_message(.false., timer(:), 'Write RH_fluxes')
      
      ! Allocate the arrays for the RH_fluxes
      allocate (RH_fluxes_phi_even_vs_kykxzts( naky, nakx, nztot, ntubes, nspec))
      allocate (RH_fluxes_phi_odd_vs_kykxzts(  naky, nakx, nztot, ntubes, nspec))
      allocate (RH_fluxes_apar_even_vs_kykxzts(naky, nakx, nztot, ntubes, nspec))
      allocate (RH_fluxes_apar_odd_vs_kykxzts( naky, nakx, nztot, ntubes, nspec))
      allocate (RH_fluxes_bpar_even_vs_kykxzts(naky, nakx, nztot, ntubes, nspec))
      allocate (RH_fluxes_bpar_odd_vs_kykxzts( naky, nakx, nztot, ntubes, nspec))

      ! Calculate the RH inertia (kx,tube,s); RH fluxes(kx,tube,s)
      if (debug) write (*, *) 'diagnostics::diagnostics_stella::write_RH_fluxes'

      ! TODO-RN : implement for radial variation and full flux surface

      ! Calculate the RH_fluxes for a flux tube simulation
      if (write_RH_inertia_fluxes) then
         call get_RH_fluxes_fluxtube(gnew, &
                RH_fluxes_phi_even_vs_kykxzts,  RH_fluxes_phi_odd_vs_kykxzts, &
                RH_fluxes_apar_even_vs_kykxzts, RH_fluxes_apar_odd_vs_kykxzts,&
                RH_fluxes_bpar_even_vs_kykxzts, RH_fluxes_bpar_odd_vs_kykxzts)
      end if

      ! Write the RH_fluxes to the netcdf file
      if (proc0 .and. write_RH_inertia_fluxes) then
         call write_RH_fluxes_phi_nc(nout, RH_fluxes_phi_even_vs_kykxzts, RH_fluxes_phi_odd_vs_kykxzts)
         if (include_apar) call write_RH_fluxes_apar_nc(nout, RH_fluxes_apar_even_vs_kykxzts, RH_fluxes_apar_odd_vs_kykxzts)
         if (include_bpar) call write_RH_fluxes_bpar_nc(nout, RH_fluxes_bpar_even_vs_kykxzts, RH_fluxes_bpar_odd_vs_kykxzts)

      end if

      ! Deallocate the arrays for the RH_fluxes
      deallocate (RH_fluxes_phi_even_vs_kykxzts,  RH_fluxes_phi_odd_vs_kykxzts)
      deallocate (RH_fluxes_apar_even_vs_kykxzts, RH_fluxes_apar_odd_vs_kykxzts)
      deallocate (RH_fluxes_bpar_even_vs_kykxzts, RH_fluxes_bpar_odd_vs_kykxzts)

       ! End timer
       if (proc0) call time_message(.false., timer(:), 'Write RH_fluxes')
 
   end subroutine write_RH_fluxes_to_netcdf_file
   
 
   !============================================================================
   !========== CALCULATE AND WRITE RH_PHI TO NETCDF FILE =======================
   !============================================================================
   subroutine write_RH_phi_I_to_netcdf_file(nout, timer)

      ! Data
      use arrays_dist_fn, only: gnew

      ! Dimensions
      use parameters_kxky_grids, only: nakx
      use zgrid, only: nztot, ntubes
      use species, only: nspec
      
      ! Flags 
      use parameters_physics, only: full_flux_surface

      ! Write to netcdf file 
      use stella_io, only: write_RH_phi_I_nc
      
      ! Routines
      use job_manage, only: time_message
      use mp, only: proc0
      
      ! Input file
      use parameters_diagnostics, only: write_RH_inertia_fluxes

      implicit none 

      ! The pointer in the netcdf file and a timer
      real, dimension(:), intent(in out) :: timer   
      integer, intent(in) :: nout    

      ! Variables needed to write and calculate diagnostics 
      complex, dimension(:, :, :, :), allocatable :: RH_phi_I_vs_kxzts

      !---------------------------------------------------------------------- 

      ! Only continue if the RH_inertia_fluxes have to be written
      if (.not. write_RH_inertia_fluxes) return  

      ! Start timer
      if (proc0) call time_message(.false., timer(:), 'Write RH_phi_I')

      ! Allocate the array for RH_phi_I
      allocate (RH_phi_I_vs_kxzts(nakx, nztot, ntubes, nspec))

      ! Calculate the RH phi
      if (debug) write (*, *) 'diagnostics::diagnostics_stella::write_RH_phi_I'

      ! TODO-RN : implement for radial variation and full flux surface

      ! Calculate the RH_phi_I for a flux tube simulation
      if (write_RH_inertia_fluxes) then
         call get_RH_phi_I_fluxtube(gnew, RH_phi_I_vs_kxzts)
      end if

      ! Write the RH_phi_I to the netcdf file
      if (proc0 .and. write_RH_inertia_fluxes) call write_RH_phi_I_nc(nout, RH_phi_I_vs_kxzts)

      ! Deallocate the arrays for the RH_phi_I
      deallocate (RH_phi_I_vs_kxzts)

       ! End timer
       if (proc0) call time_message(.false., timer(:), 'Write RH_phi_I')

   end subroutine write_RH_phi_I_to_netcdf_file


   !============================================================================
   !====================== GET RH_inertia FOR THE FLUX TUBE =====================
   !============================================================================
   subroutine get_RH_inertia_fluxtube(RH_inertia)

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

      ! The RH inertia is returned with dimensions (kx, z, tube, spec)
      complex, dimension(:, -nzgrid:, :, :), intent(out) :: RH_inertia

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
                                        RH_fluxes_apar_even, RH_fluxes_apar_odd,&
                                        RH_fluxes_bpar_even, RH_fluxes_bpar_odd)

      use zgrid, only: nzgrid, ntubes
      use species, only: spec, nspec
      use vpamu_grids, only: vpa, mu, vperp2, integrate_vmu
      use vpamu_grids, only: maxwell_mu, ztmax, maxwell_fac, maxwell_vpa
      use parameters_kxky_grids, only: naky, nakx, nx
      use grids_kxky, only: aky
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

      ! Import temp arrays g1, g2 with dimensions (nky, nkx, -nzgrid:nzgrid, ntubes, -vmu-layout-)
      use arrays_dist_fn, only: integrand_even   => g0
      use arrays_dist_fn, only: integrand_odd    => g1

      implicit none

      ! Gyroaveraged ExB and NL term in k-space
      complex, dimension(naky, nakx) :: vchix_gyro, NL_term

      ! Gyroaveraged ExB term, distribution function, and integrand in x and ky
      complex, dimension(naky, nx) :: vchix_gyro_ky_x, g_ky_x, NL_term_ky_x

      ! The distribution function enters with dimensions (ky, kx, z, tube, ivmus)
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g

      ! The RH fluxes are returned with dimensions (ky, kx, z, tube, s)
      complex, dimension(:, :, -nzgrid:, :, :), intent(out) :: RH_fluxes_phi_even,  RH_fluxes_phi_odd
      complex, dimension(:, :, -nzgrid:, :, :), intent(out) :: RH_fluxes_apar_even, RH_fluxes_apar_odd
      complex, dimension(:, :, -nzgrid:, :, :), intent(out) :: RH_fluxes_bpar_even, RH_fluxes_bpar_odd

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
         return 
      end if

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
      !TODO: to implement
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


!###############################################################################
!############################ INITALIZE & FINALIZE #############################
!###############################################################################

   !============================================================================
   !======================== INITALIZE THE DIAGNOSTICS =========================
   !============================================================================
   subroutine init_diagnostics_RH_inertia_fluxes()
      !TODO-RN: call only when needed??

      use mp, only: proc0

      ! Dimensions
      use parameters_kxky_grids, only: nakx
      use grids_kxky, only: akx
      use zgrid, only: nzgrid, ntubes
      use stella_layouts, only: vmu_lo
      use stella_layouts, only: iv_idx, imu_idx, is_idx
      use vpamu_grids, only: vpa, vperp2, mu
      use geometry, only: bmag

      implicit none

      real :: energyval, muval
      complex :: Q_fac
      real :: transit_int_tau_b_pls, transit_int_tau_b_min
      complex :: transit_int_eiQJ0_pls, transit_int_eiQJ0_min
      complex :: integrand_tmp_pls, integrand_tmp_min
      complex :: tmp

      integer :: ivmu, iv, imu, is, ia, iz, it, ikx
      ia = 1

      ! Only debug on the first processor
      debug = debug .and. proc0

      ! Allocate the arrays for the Rosenbluth-Hinton integrand term
      allocate (RH_integrand_even(nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); RH_integrand_even = 0.
      allocate (RH_integrand_odd( nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); RH_integrand_odd  = 0.

      ! Evaluate the transit averages
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)

         do ikx = 1, nakx
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid

               ! Determine energy and magnetic moment
               energyval = vpa(iv)**2 + vperp2(ia,iz,imu)
               muval     = mu(imu)

               ! Evaluate transit averages for vpa and -vpa
               call eval_transit_ints(energyval, muval, sign(1., vpa(iv)), ikx, is, transit_int_eiQJ0_pls, transit_int_tau_b_pls)
               call eval_transit_ints(energyval, muval, sign(1.,-vpa(iv)), ikx, is, transit_int_eiQJ0_min, transit_int_tau_b_min)

               ! For trapped particles, the result is the average of +vpa and -vpa transit averages
               ! Note: it follows that transit_int_eiQJ0_pls is even in vpa
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
               !RH_integrand_odd( ikx,iz,it,ivmu) = -0.5*(integrand_tmp_pls-integrand_tmp_min)
               !TODO: Numerical results show we need a minus sign here. Why??

               !write(*,*) 'RH integrand even:'
               !write(*,*) RH_integrand_even(ikx,iz,it,ivmu)
               !write(*,*) 'RH integrand odd:'
               !write(*,*) RH_integrand_odd(ikx,iz,it,ivmu)

               end do !iz
            end do !it
         end do !ikx
      end do !ivmu

      ! Evaluate and write RH_inertia to netcdf file
      call write_RH_inertia_to_netcdf_file()

      ! Write RH_integrand_(even/odd) to netcdf file
      call write_RH_integrands_to_netcdf_file()

   end subroutine init_diagnostics_RH_inertia_fluxes


   !============================================================================
   !======================== FINALIZE THE DIAGNOSTICS ==========================
   !============================================================================
   subroutine finish_diagnostics_RH_inertia_fluxes()

      use mp, only: proc0

      implicit none

      ! Deallocate the arrays for the Rosenbluth-Hinton integrand term
      if (allocated(RH_integrand_even)) deallocate (RH_integrand_even)
      if (allocated(RH_integrand_odd )) deallocate (RH_integrand_odd)


   end subroutine finish_diagnostics_RH_inertia_fluxes

   !============================================================================
   !======================== USEFUL MISC FUNCTIONS =============================
   !============================================================================

   ! Evaluate RH transit averages
   subroutine eval_transit_ints(energy, mu, sigma, ikx, is, transit_int_eiQJ0, bounce_time)

      use geometry, only: bmag, dl_over_b
      use grids_kxky, only: akx
      use zgrid, only: nzgrid
      use constants, only: zi

      implicit none

      real,    intent(in)  :: energy, mu, sigma
      integer, intent(in)  :: ikx, is
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

         call eval_transit_int_integrand(energy, mu, sigma, ikx, iz, is, .false., integrand_eiQJ0(iz))
         call eval_transit_int_integrand(energy, mu, sigma, ikx, iz, is, .true.,  integrand_tau_b(iz))
         !write(*,*) integrand_eiQJ0(iz)

      end do

      ! Evaluate integrals (integrand has 1/vpa factor, need to integrate dl/vpa (...) = dl/B * B (...) )
      transit_int_eiQJ0 = sum(integrand_eiQJ0 * bmag(ia,:) * dl_over_b(ia, :))
      bounce_time       = sum(integrand_tau_b * bmag(ia,:) * dl_over_b(ia, :))

   end subroutine eval_transit_ints


   ! Evaluate integrand in RH transit average
   subroutine eval_transit_int_integrand(energy, mu, sigma, ikx, iz, is, bounce_time_bool, transit_avg_integrand)

      use geometry, only: bmag
      use species, only: spec
      use spfunc, only: j0
      use arrays_dist_fn, only: kperp2
      use grids_kxky, only: akx

      implicit none

      real,    intent(in)  :: energy, mu, sigma ! energy=vpa^2+vperp^2, mu=vperp^2/(2B)
      integer, intent(in)  :: ikx, iz, is
      logical, intent(in)  :: bounce_time_bool ! if true, evaluate integrand for bounce time
      complex, intent(out) :: transit_avg_integrand

      real    :: vpa2, vpa, vperp2
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
            aj0x = j0( sqrt(kperp2(1,ikx,ia,iz)*vperp2) * spec(is)%bess_fac * spec(is)%smz_psi0 / bmag(ia,iz) )

            ! Evaluate Q factor
            call eval_Q_fac(vpa, akx(ikx), iz, is, Q_fac)

            ! Integrand
            transit_avg_integrand = exp(-Q_fac) * aj0x / abs(vpa)

         end if
      end if

   end subroutine eval_transit_int_integrand

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
              * geo_surf%qinp*btor(iz)*Rmajor(iz)/geo_surf%rhoc ! Note Btor*Rmajor should be constant along field-line

!      write(*, *) btor(iz)*Rmajor(iz)
!      write(*, *) spec(is)%smz_psi0

   end subroutine eval_Q_fac



end module diagnostics_RH_inertia_fluxes

