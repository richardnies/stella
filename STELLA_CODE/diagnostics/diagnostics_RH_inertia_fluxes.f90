
!###############################################################################
!############### DIAGNOSE ROSENBLUTH-HINTON INERTIA AND FLUXES #################
!###############################################################################
! 
! Routines for calculating and writing the "Rosenbluth-Hinton" inertia and fluxes,
! governing the evolution of long-wavelength stationary zonal flows.
! 
! The RH_inertia             is denoted by RH_inertia
! The RH_integrand_even_vs_kxztsvpamu is denoted by RH_integrand_even
! The RH_integrand_odd_vs_kxztsvpamu  is denoted by RH_integrand_odd
! The RH_U_parallel_fac_vs_ztsvpamu    is denoted by RH_U_parallel_fac
! The RH_fluxes_phi_even_vs_kykxzts       is denoted by RH_fluxes_phi_even
! The RH_fluxes_phi_odd_vs_kykxzts        is denoted by RH_fluxes_phi_odd
! The RH_fluxes_apar_even_vs_kykxzts      is denoted by RH_fluxes_apar_even
! The RH_fluxes_apar_odd_vs_kykxzts       is denoted by RH_fluxes_apar_odd
! The RH_fluxes_bpar_even_vs_kykxzts      is denoted by RH_fluxes_bpar_even
! The RH_fluxes_bpar_odd_vs_kykxzts       is denoted by RH_fluxes_bpar_odd
! The RH_fluxes_coll_vs_kxzts             is denoted by RH_fluxes_coll
! 
!###############################################################################
 
module diagnostics_RH_inertia_fluxes

   implicit none
 
   public :: write_RH_integrands_to_netcdf_file
   public :: write_RH_bounce_drift_to_netcdf_file
   public :: write_RH_fluxes_to_netcdf_file
   public :: write_RH_inertia_to_netcdf_file
   public :: write_RH_phi_I_to_netcdf_file

   private

   ! Debugging
   logical :: debug = .false.

contains

!###############################################################################
!###################### WRITE RH_INERTIA_FLUXES ################################
!###############################################################################

   !============================================================================
   !========== CALCULATE AND WRITE RH_INTEGRANDS TO NETCDF FILE ================
   !============================================================================
   !============================================================================
   !=========== WRITE THE BOUNCE-AVERAGED RADIAL DRIFT TO NETCDF ===============
   !============================================================================
   !> Gathers RH_drift_bounce_avg from the vmu layout onto proc0 and writes it.
   !> Time-independent, so written once at initialisation like the integrands.
   subroutine write_RH_bounce_drift_to_netcdf_file()

      use parameters_diagnostics, only: write_RH_bounce_drift
      use rosenbluth_hinton, only: RH_drift_bounce_avg
      use vpamu_grids, only: nvpa, nmu
      use zgrid, only: nztot
      use species, only: nspec
      use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      use stella_io, only: write_RH_bounce_drift_nc
      use mp, only: sum_reduce, proc0

      implicit none

      integer :: ivmu, iv, imu, is
      real, dimension(:, :, :, :), allocatable :: drift_vs_ztsvpamu

      if (.not. write_RH_bounce_drift) return

      allocate (drift_vs_ztsvpamu(nztot, nspec, nvpa, nmu)); drift_vs_ztsvpamu = 0.

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         drift_vs_ztsvpamu(:, is, iv, imu) = RH_drift_bounce_avg(:, ivmu)
      end do

      call sum_reduce(drift_vs_ztsvpamu, 0)

      if (proc0) call write_RH_bounce_drift_nc(drift_vs_ztsvpamu)

      deallocate (drift_vs_ztsvpamu)

   end subroutine write_RH_bounce_drift_to_netcdf_file

   subroutine write_RH_integrands_to_netcdf_file()

      use rosenbluth_hinton, only: RH_integrand_even, RH_integrand_odd

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

      use rosenbluth_hinton, only: RH_inertia

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

      !---------------------------------------------------------------------- 

      ! Only continue if the RH_inertia_fluxes have to be written
      if (.not. write_RH_inertia_fluxes) return

      ! Calculate the RH inertia (kx,tube,s); RH fluxes(kx,tube,s)
      if (debug) write (*, *) 'diagnostics::diagnostics_stella::write_RH_inertia'

      ! Write the RH_inertia to the netcdf file
      if (proc0 .and. write_RH_inertia_fluxes) call write_RH_inertia_nc(RH_inertia)

   end subroutine write_RH_inertia_to_netcdf_file
 
   !============================================================================
   !========== CALCULATE AND WRITE RH_FLUXES TO NETCDF FILE ====================
   !============================================================================
   subroutine write_RH_fluxes_to_netcdf_file(nout, timer)

      use rosenbluth_hinton, only: get_RH_fluxes_fluxtube

      ! Data
      use arrays_dist_fn, only: gnew

      ! Dimensions
      use parameters_kxky_grids, only: naky, nakx
      use zgrid, only: nztot, ntubes
      use species, only: nspec
      
      ! Flags 
      use parameters_physics, only: full_flux_surface

      ! Write to netcdf file 
      use stella_io, only: write_RH_fluxes_phi_nc, write_RH_fluxes_apar_nc, write_RH_fluxes_bpar_nc, write_RH_fluxes_coll_nc
      use stella_io, only: write_RH_fluxes_drift_nc
      
      ! Routines
      use job_manage, only: time_message
      use mp, only: proc0
      
      ! Input file
      use parameters_diagnostics, only: write_RH_inertia_fluxes
      
      ! Physics parameters
      use parameters_physics, only: include_apar, include_bpar
      use dissipation, only: include_collisions

      implicit none 

      ! The pointer in the netcdf file and a timer
      real, dimension(:), intent(in out) :: timer   
      integer, intent(in) :: nout    

      ! Variables needed to write and calculate diagnostics 
      complex, dimension(:, :, :, :, :), allocatable :: RH_fluxes_phi_even_vs_kykxzts,  RH_fluxes_phi_odd_vs_kykxzts
      complex, dimension(:, :, :, :, :), allocatable :: RH_fluxes_apar_even_vs_kykxzts, RH_fluxes_apar_odd_vs_kykxzts
      complex, dimension(:, :, :, :, :), allocatable :: RH_fluxes_bpar_even_vs_kykxzts, RH_fluxes_bpar_odd_vs_kykxzts
      complex, dimension(:, :, :, :),    allocatable :: RH_fluxes_coll_vs_kxzts
      complex, dimension(:, :, :, :),    allocatable :: RH_fluxes_drift_vs_kxzts

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
      allocate (RH_fluxes_coll_vs_kxzts(                     nakx, nztot, ntubes, nspec))
      allocate (RH_fluxes_drift_vs_kxzts(                    nakx, nztot, ntubes, nspec))

      ! Calculate the RH inertia (kx,tube,s); RH fluxes(kx,tube,s)
      if (debug) write (*, *) 'diagnostics::diagnostics_stella::write_RH_fluxes'

      ! TODO-RN : implement for radial variation and full flux surface

      ! Calculate the RH_fluxes for a flux tube simulation
      if (write_RH_inertia_fluxes) then
         call get_RH_fluxes_fluxtube(gnew, &
                RH_fluxes_phi_even_vs_kykxzts,  RH_fluxes_phi_odd_vs_kykxzts, &
                RH_fluxes_apar_even_vs_kykxzts, RH_fluxes_apar_odd_vs_kykxzts, &
                RH_fluxes_bpar_even_vs_kykxzts, RH_fluxes_bpar_odd_vs_kykxzts, &
                RH_fluxes_coll_vs_kxzts, RH_fluxes_drift_vs_kxzts)
      end if

      ! Write the RH_fluxes to the netcdf file
      if (proc0 .and. write_RH_inertia_fluxes) then
         call write_RH_fluxes_phi_nc(nout, RH_fluxes_phi_even_vs_kykxzts, RH_fluxes_phi_odd_vs_kykxzts)
         if (include_apar) call write_RH_fluxes_apar_nc(nout, RH_fluxes_apar_even_vs_kykxzts, RH_fluxes_apar_odd_vs_kykxzts)
         if (include_bpar) call write_RH_fluxes_bpar_nc(nout, RH_fluxes_bpar_even_vs_kykxzts, RH_fluxes_bpar_odd_vs_kykxzts)
         if (include_collisions) call write_RH_fluxes_coll_nc(nout, RH_fluxes_coll_vs_kxzts)
         call write_RH_fluxes_drift_nc(nout, RH_fluxes_drift_vs_kxzts)

      end if

      ! Deallocate the arrays for the RH_fluxes
      deallocate (RH_fluxes_phi_even_vs_kykxzts,  RH_fluxes_phi_odd_vs_kykxzts)
      deallocate (RH_fluxes_apar_even_vs_kykxzts, RH_fluxes_apar_odd_vs_kykxzts)
      deallocate (RH_fluxes_bpar_even_vs_kykxzts, RH_fluxes_bpar_odd_vs_kykxzts)
      deallocate (RH_fluxes_coll_vs_kxzts)
      deallocate (RH_fluxes_drift_vs_kxzts)

      ! End timer
      if (proc0) call time_message(.false., timer(:), 'Write RH_fluxes')
 
   end subroutine write_RH_fluxes_to_netcdf_file
   
 
   !============================================================================
   !========== CALCULATE AND WRITE RH_PHI TO NETCDF FILE =======================
   !============================================================================
   subroutine write_RH_phi_I_to_netcdf_file(nout, timer)

      use rosenbluth_hinton, only: get_RH_phi_I_fluxtube

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

end module diagnostics_RH_inertia_fluxes
