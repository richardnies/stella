module implicit_solve

   implicit none

   public :: time_implicit_advance
   public :: sweep_zed_zonal
   public :: advance_implicit_terms
   public :: get_gke_rhs
   public :: sweep_g_zext

   private

   real, dimension(2, 3) :: time_implicit_advance = 0.

contains

   subroutine advance_implicit_terms(g, phi, apar)

      use mp, only: proc0
      use job_manage, only: time_message
      use stella_layouts, only: vmu_lo
      use physics_flags, only: include_apar
      use g_tofrom_h, only: gbar_to_g
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky, nakx
      use dist_fn_arrays, only: g1
      use run_parameters, only: stream_matrix_inversion
      use run_parameters, only: use_deltaphi_for_response_matrix
      use run_parameters, only: tupwnd_p => time_upwind_plus
      use run_parameters, only: tupwnd_m => time_upwind_minus
      use run_parameters, only: fphi
      use fields, only: advance_fields, fields_updated
      use extended_zgrid, only: map_to_extended_zgrid, map_from_extended_zgrid
      use extended_zgrid, only: nsegments, nzed_segment

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi, apar

      integer :: nz_ext
      complex, dimension(:, :, :, :), allocatable :: phi_old, apar_old
      complex, dimension(:, :, :, :), allocatable :: phi_source, apar_source
      character(5) :: dist_choice

      if (proc0) call time_message(.false., time_implicit_advance(:, 1), ' Implicit time advance')

      allocate (phi_source(naky, nakx, -nzgrid:nzgrid, ntubes))
      allocate (apar_source(naky, nakx, -nzgrid:nzgrid, ntubes))
      !> dist_choice indicates whether the non-Boltzmann part of the pdf (h) is evolved
      !> in parallel streaming or if the guiding centre distribution (g = <f>) is evolved
      dist_choice = 'g'
      !> if using delphi formulation for response matrix, then phi = phi^n replaces
      !> phi^{n+1} in the inhomogeneous GKE; else set phi_{n+1} to zero in inhomogeneous equation
      if (use_deltaphi_for_response_matrix) then
         phi_source = phi
      else
         phi_source = tupwnd_m * phi
      end if

      ! save the incoming pdf and phi, as they will be needed later
      g1 = g
      allocate (phi_old(naky, nakx, -nzgrid:nzgrid, ntubes))
      phi_old = phi

      if (include_apar) then
         ! save the incoming apar, as it will be needed later
         allocate (apar_old(naky, nakx, -nzgrid:nzgrid, ntubes))
         apar_old = apar
         ! when solving for 'inhomogeneous' piece of the pdf, use part of apar weighted towards
         ! previous time level
         apar_source = tupwnd_m * apar
         ! if including apar, the incoming pdf is gbar = g + (Ze/T)*<vpa*apar/c>*F0;
         ! convert from gbar to g
         call gbar_to_g(g1, apar, 1.0)
         ! set apar=0, as in update_pdf it expects apar to be apar^{n+1} or zero
         apar = 0.0
      end if

      ! solve for the 'inhomogeneous' piece of the pdf, stored in g_scratch
      call update_pdf

      fields_updated = .false.

      ! we now have g_{inh}^{n+1}
      ! calculate associated fields (phi_{inh}^{n+1})
      call advance_fields(g, phi, apar, dist=trim(dist_choice))

      ! solve response_matrix*(phi^{n+1}-phi^{n*}) = phi_{inh}^{n+1}-phi^{n*}
      ! phi = phi_{inh}^{n+1}-phi^{n*} is input and overwritten by phi = phi^{n+1}-phi^{n*}
      if (use_deltaphi_for_response_matrix) phi = phi - phi_old
      if (proc0) call time_message(.false., time_implicit_advance(:, 3), ' (back substitution)')
      call invert_parstream_response(phi, apar)
      if (proc0) call time_message(.false., time_implicit_advance(:, 3), ' (back substitution)')

      !> If using deltaphi formulation, must account for fact that phi = phi^{n+1}-phi^{n*}, but
      !> tupwnd_p should multiply phi^{n+1}
      if (use_deltaphi_for_response_matrix) phi = phi + phi_old
      ! get time-centered phi
      phi_source = tupwnd_m * phi_old + tupwnd_p * phi
      ! get time-centered apar
      if (include_apar) apar_source = tupwnd_m * apar_old + tupwnd_p * apar

      ! solve for the final, updated pdf now that we have phi^{n+1} and apar^{n+1}
      call update_pdf

      deallocate (phi_old)
      deallocate (phi_source, apar_source)
      if (allocated(apar_old)) deallocate (apar_old)

      if (proc0) call time_message(.false., time_implicit_advance(:, 1), ' Stream advance')

   contains

      subroutine update_pdf

         use extended_zgrid, only: neigen

         integer :: ie, it, iky, ivmu
         integer :: ulim
         complex, dimension(:), allocatable :: pdf1, pdf2
         complex, dimension(:), allocatable :: phiext
         complex, dimension(:), allocatable :: aparext, aparext_new

         ! start the timer for the pdf update
         if (proc0) call time_message(.false., time_implicit_advance(:, 2), ' (bidiagonal solve)')

         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            ! solve for the pdf, given the sources for phi and the pdf on the RHS of the GK equation
            ! we do this on set of connected zed segments at a time
            do iky = 1, naky
               do it = 1, ntubes
                  do ie = 1, neigen(iky)
                     ! nz_ext is the number of grid points in the extended zed domain
                     nz_ext = nsegments(ie, iky) * nzed_segment + 1
                     ! pdf1 and pdf2 will be scratch arrays needed to compute the pdf itself,
                     ! as well as contributions to the GK equation
                     allocate (pdf1(nz_ext), pdf2(nz_ext))
                     allocate (phiext(nz_ext))
                     allocate (aparext(nz_ext)); aparext = 0.0
                     allocate (aparext_new(nz_ext)); aparext_new = 0.0
                     ! map the incoming pdf 'g1' onto the extended zed domain and call it 'pdf1'
                     call map_to_extended_zgrid(it, ie, iky, g1(iky, :, :, :, ivmu), pdf1, ulim)
                     ! map the incoming potential 'phi_source' onto the extended zed domain and call it 'phiext'
                     call map_to_extended_zgrid(it, ie, iky, phi_source(iky, :, :, :), phiext, ulim)
                     ! map incoming parallel magnetic vector potetial 'apar_sosurce' onto
                     ! extended zed domain and call 'aparext'
                     if (include_apar) then
                        call map_to_extended_zgrid(it, ie, iky, apar_source(iky, :, :, :), aparext, ulim)
                        call map_to_extended_zgrid(it, ie, iky, apar(iky, :, :, :), aparext_new, ulim)
                     end if
                     ! calculate the RHS of the GK equation (using pdf1 and phi_source as the
                     ! pdf and potential, respectively) and store it in pdf2
                     call get_gke_rhs(ivmu, iky, ie, pdf1, phiext, aparext, aparext_new, pdf2)
                     ! given the RHS of the GK equation (pdf2), solve for the pdf at the
                     ! new time level by sweeping in zed on the extended domain;
                     ! the rhs is input as 'pdf2' and over-written with the updated solution for the pdf
                     call sweep_g_zext(iky, ie, it, ivmu, pdf2)
                     ! map the pdf 'pdf2' from the extended zed domain
                     ! to the standard zed domain; the mapped pdf is called 'g'
                     call map_from_extended_zgrid(it, ie, iky, pdf2, g(iky, :, :, :, ivmu))
                     deallocate (pdf1, pdf2, phiext, aparext, aparext_new)
                  end do
               end do
            end do
         end do

         ! stop the timer for the pdf update
         if (proc0) call time_message(.false., time_implicit_advance(:, 2), ' (bidiagonal solve)')

      end subroutine update_pdf

   end subroutine advance_implicit_terms

   subroutine get_gke_rhs(ivmu, iky, ie, pdf, phi, apar, aparnew, rhs)

      use kt_grids, only: naky, nakx
      use physics_flags, only: include_apar

      implicit none

      integer, intent(in) :: ivmu, iky, ie
      complex, dimension(:), intent(in) :: pdf
      complex, dimension(:), intent(in) :: phi, apar, aparnew
      complex, dimension(:), intent(out) :: rhs

      integer :: nz_ext
      complex, dimension(:), allocatable :: rhs_fields

      ! now have phi^{n+1} for non-negative kx
      ! obtain RHS of GK eqn

      ! nz_ext is the number of grid points in this extended zed domain
      nz_ext = size(pdf)
      allocate (rhs_fields(nz_ext))

      ! NB: rhs is used as a scratch array in get_contributions_from_fields
      ! so be careful not to move get_contributions_from_pdf before it, or rhs will be over-written
      call get_contributions_from_fields(phi, apar, aparnew, ivmu, iky, ie, rhs, rhs_fields)
      call get_contributions_from_pdf(pdf, ivmu, iky, ie, rhs)

      ! construct RHS of GK eqn
      rhs = rhs + rhs_fields

      deallocate (rhs_fields)

   end subroutine get_gke_rhs

   !> get_contributions_from_fields takes as input the appropriately averaged
   !> electrostatic potential phi and magnetic vector potential components apar
   !> and returns in rhs the sum of the source terms
   !> involving phi and apar that appear on the RHS of the GK equation when g is the pdf
   subroutine get_contributions_from_fields(phi, apar, aparnew, ivmu, iky, ie, scratch, rhs)

      use physics_flags, only: include_apar
      use extended_zgrid, only: map_to_iz_ikx_from_izext

      implicit none

      complex, dimension(:), intent(in) :: phi, apar, aparnew
      integer, intent(in) :: ivmu, iky, ie
      complex, dimension(:), intent(out) :: scratch, rhs

      integer, dimension(:), allocatable :: iz_from_izext, ikx_from_izext
      integer :: nz_ext

      ! nz_ext is the number of grid points in the extended zed domain
      nz_ext = size(phi)

      ! determine the mapping from the extended domain zed index (izext) to the
      ! zed and kx domain indices (iz, ikx)
      allocate (iz_from_izext(nz_ext))
      allocate (ikx_from_izext(nz_ext))
      call map_to_iz_ikx_from_izext(iky, ie, iz_from_izext, ikx_from_izext)

      ! calculate the contributions to the RHS of the GKE due to source terms proportional to phi
      call get_contributions_from_phi(phi, ivmu, iky, iz_from_izext, ikx_from_izext, scratch, rhs)
      ! if advancing apar, get its contribution to the RHS of the GKE and add to phi contribution
      if (include_apar) then
         call get_contributions_from_apar(apar, aparnew, ivmu, iky, iz_from_izext, ikx_from_izext, scratch, rhs)
      end if

      deallocate (iz_from_izext, ikx_from_izext)

   end subroutine get_contributions_from_fields

   !> get_contributions_from_phi takes as input the appropriately averaged
   !> electrostatic potential phi and returns in rhs the sum of the source terms
   !> involving phi that appear on the RHS of the GK equation when g is the pdf
   subroutine get_contributions_from_phi(phi, ivmu, iky, iz_from_izext, ikx_from_izext, scratch, rhs)

      use stella_time, only: code_dt
      use species, only: spec
      use zgrid, only: nzgrid, ntubes
      use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use vpamu_grids, only: vpa
      use kt_grids, only: naky, nakx
      use run_parameters, only: driftkinetic_implicit, maxwellian_normalization
      use run_parameters, only: maxwellian_inside_zed_derivative
      use run_parameters, only: drifts_implicit
      use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      use neoclassical_terms, only: include_neoclassical_terms
      use neoclassical_terms, only: dfneo_dvpa
      use extended_zgrid, only: map_to_iz_ikx_from_izext

      implicit none

      complex, dimension(:), intent(in) :: phi
      integer, intent(in) :: ivmu, iky
      integer, dimension(:), intent(in) :: iz_from_izext, ikx_from_izext
      complex, dimension(:), intent(out) :: scratch, rhs

      real, dimension(:), allocatable :: z_scratch
      integer :: ia, iz, iv, imu, is
      integer :: nz_ext

      ia = 1
      iv = iv_idx(vmu_lo, ivmu)
      imu = imu_idx(vmu_lo, ivmu)
      is = is_idx(vmu_lo, ivmu)

      ! nz_ext is the number of grid points in the extended zed domain
      nz_ext = size(phi)

      ! allocate a 1d array in zed for use as a scratch array
      allocate (z_scratch(-nzgrid:nzgrid))

      ! set scratc to be phi or <phi> depending on whether parallel streaming is
      ! implicit or only implicit in the kperp = 0 (drift kinetic) piece
      if (driftkinetic_implicit) then
         scratch = phi
      else
         call gyro_average_zext(iky, ivmu, ikx_from_izext, iz_from_izext, phi, scratch)
      end if

      call add_streaming_contribution_phi
      if (drifts_implicit) call add_drifts_contribution_phi

      deallocate (z_scratch)

   contains

      subroutine add_streaming_contribution_phi

         use extended_zgrid, only: fill_zext_ghost_zones
         use parallel_streaming, only: get_zed_derivative_extended_domain
         use parallel_streaming, only: center_zed
         use parallel_streaming, only: gradpar_c, stream_sign

         integer :: izext
         complex :: scratch_left, scratch_right

         ! fill ghost zones beyond ends of extended zed domain for <phi>
         ! and store values in scratch_left and scratch_right
         call fill_zext_ghost_zones(iky, scratch, scratch_left, scratch_right)

         ! obtain the zed derivative of <phi> (stored in scratch) and store in rhs
         call get_zed_derivative_extended_domain(iv, scratch, scratch_left, scratch_right, rhs)

         if (.not. maxwellian_normalization) then
            ! center Maxwellian factor in mu
            ! and store in dummy variable z_scratch
            z_scratch = maxwell_mu(ia, :, imu, is)
            call center_zed(iv, z_scratch, -nzgrid)
            ! multiply by Maxwellian factor
            do izext = 1, nz_ext
               rhs(izext) = rhs(izext) * z_scratch(iz_from_izext(izext))
            end do
         end if

         ! NB: could do this once at beginning of simulation to speed things up
         ! this is vpa*Z/T*exp(-vpa^2)
         z_scratch = vpa(iv) * spec(is)%zt
         if (.not. maxwellian_normalization) z_scratch = z_scratch * maxwell_vpa(iv, is) * maxwell_fac(is)
         ! if including neoclassical correction to equilibrium distribution function
         ! then must also account for -vpa*dF_neo/dvpa*Z/T
         ! CHECK TO ENSURE THAT DFNEO_DVPA EXCLUDES EXP(-MU*B/T) FACTOR !!
         if (include_neoclassical_terms) then
            do iz = -nzgrid, nzgrid
               z_scratch(iz) = z_scratch(iz) - 0.5 * dfneo_dvpa(ia, iz, ivmu) * spec(is)%zt
            end do
            call center_zed(iv, z_scratch, -nzgrid)
         end if

         if (stream_sign(iv) > 0) then
            z_scratch = z_scratch * gradpar_c(:, -1) * code_dt * spec(is)%stm_psi0
         else
            z_scratch = z_scratch * gradpar_c(:, 1) * code_dt * spec(is)%stm_psi0
         end if

         do izext = 1, nz_ext
            rhs(izext) = -z_scratch(iz_from_izext(izext)) * rhs(izext)
         end do

      end subroutine add_streaming_contribution_phi

      subroutine add_drifts_contribution_phi

         use constants, only: zi
         use kt_grids, only: nakx, naky
         use kt_grids, only: aky, akx
         use dist_fn_arrays, only: wstar, wdriftx_phi, wdrifty_phi
         use parallel_streaming, only: center_zed
         use extended_zgrid, only: periodic

         integer :: izext, iz, ikx

         ! 'scratch' starts out as the gyro-average of phi, evaluated at zed grid points
         do izext = 1, nz_ext
            ikx = ikx_from_izext(izext)
            iz = iz_from_izext(izext)
            scratch(izext) = zi * scratch(izext) * (akx(ikx) * wdriftx_phi(ia, iz, ivmu) &
                                                    + aky(iky) * (wdrifty_phi(ia, iz, ivmu) + wstar(ia, iz, ivmu)))
         end do
         call center_zed(iv, scratch, 1, periodic(iky))

         rhs = rhs + scratch

      end subroutine add_drifts_contribution_phi

   end subroutine get_contributions_from_phi

   !> get_contributions_from_apar takes as input the appropriately averaged
   !> parallel component of the vector potential, apar, and returns in rhs the sum of the source terms
   !> involving apar that appear on the RHS of the GK equation when g is the pdf
   subroutine get_contributions_from_apar(apar, aparnew, ivmu, iky, iz_from_izext, ikx_from_izext, scratch, rhs)

      use run_parameters, only: driftkinetic_implicit, drifts_implicit
      use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx

      implicit none

      complex, dimension(:), intent(in) :: apar, aparnew
      integer, intent(in) :: ivmu, iky
      integer, dimension(:), intent(in) :: iz_from_izext, ikx_from_izext
      complex, dimension(:), intent(out) :: scratch, rhs

      complex, dimension(:), allocatable :: scratch2
      integer :: ia, iv, imu, is
      integer :: nz_ext

      ia = 1
      iv = iv_idx(vmu_lo, ivmu)
      imu = imu_idx(vmu_lo, ivmu)
      is = is_idx(vmu_lo, ivmu)

      nz_ext = size(scratch)

      allocate (scratch2(nz_ext))

      ! set scratch to be apar or <apar> depending on whether parallel streaming is
      ! implicit or only implicit in the kperp = 0 (drift kinetic) piece
      if (driftkinetic_implicit) then
         scratch = apar
         scratch2 = aparnew
      else
         call gyro_average_zext(iky, ivmu, ikx_from_izext, iz_from_izext, apar, scratch)
         call gyro_average_zext(iky, ivmu, ikx_from_izext, iz_from_izext, aparnew, scratch2)
      end if

      call add_gbar_to_g_contribution_apar
      if (drifts_implicit) call add_drifts_contribution_apar

      deallocate (scratch2)

   contains

      subroutine add_gbar_to_g_contribution_apar

         use run_parameters, only: maxwellian_normalization
         use vpamu_grids, only: vpa, maxwell_vpa, maxwell_mu, maxwell_fac
         use parallel_streaming, only: center_zed
         use species, only: spec
         use extended_zgrid, only: periodic

         implicit none

         integer :: izext, iz
         real :: constant_factor

         constant_factor = -2.0 * spec(is)%zt_psi0 * spec(is)%stm_psi0 * vpa(iv)

         do izext = 1, nz_ext
            iz = iz_from_izext(izext)
            scratch2(izext) = constant_factor * scratch2(izext)
         end do

         ! if the pdf is not normalized by a Maxwellian then the source term contains a Maxwellian factor
         if (.not. maxwellian_normalization) then
            do izext = 1, nz_ext
               iz = iz_from_izext(izext)
               scratch2(izext) = scratch2(izext) * maxwell_vpa(iv, is) * maxwell_mu(ia, iz, imu, is) * maxwell_fac(is)
            end do
         end if
         call center_zed(iv, scratch2, 1, periodic(iky))
         rhs = rhs + scratch2

      end subroutine add_gbar_to_g_contribution_apar

      subroutine add_drifts_contribution_apar

         use constants, only: zi
         use species, only: spec
         use kt_grids, only: aky
         use dist_fn_arrays, only: wstar
         use parallel_streaming, only: center_zed
         use extended_zgrid, only: periodic
         use vpamu_grids, only: vpa

         implicit none

         integer :: izext, iz
         complex :: constant_factor

         constant_factor = -2.0 * zi * spec(is)%stm_psi0 * vpa(iv) * aky(iky)

         do izext = 1, nz_ext
            iz = iz_from_izext(izext)
            scratch(izext) = constant_factor * scratch(izext) * wstar(ia, iz, ivmu)
         end do
         call center_zed(iv, scratch, 1, periodic(iky))
         rhs = rhs + scratch

      end subroutine add_drifts_contribution_apar

   end subroutine get_contributions_from_apar

   subroutine gyro_average_zext(iky, ivmu, ikx_from_izext, iz_from_izext, phi, gyro_phi)

      use gyro_averages, only: gyro_average

      implicit none

      integer, intent(in) :: iky, ivmu
      integer, dimension(:), intent(in) :: ikx_from_izext, iz_from_izext
      complex, dimension(:), intent(in) :: phi
      complex, dimension(:), intent(out) :: gyro_phi

      integer :: izext, nz_ext

      nz_ext = size(phi)

      do izext = 1, nz_ext
         call gyro_average(phi(izext), iky, ikx_from_izext(izext), iz_from_izext(izext), ivmu, gyro_phi(izext))
      end do

   end subroutine gyro_average_zext

   !> get_contributions_from_pdf takes as an argument the evolved pdf
   !> (either guiding centre distribution g=<f> or maxwellian-normlized, non-Boltzmann distribution h/F0=f/F0+(Ze*phi/T))
   !> and the scratch array rhs, and returns the source terms that depend on the pdf in rhs
   subroutine get_contributions_from_pdf(pdf, ivmu, iky, ie, rhs)

      use constants, only: zi
      use stella_time, only: code_dt
      use species, only: spec
      use zgrid, only: nzgrid, ntubes
      use kt_grids, only: naky, nakx
      use kt_grids, only: aky, akx
      use vpamu_grids, only: vpa
      use stella_layouts, only: vmu_lo, iv_idx, is_idx
      use run_parameters, only: time_upwind_minus
      use run_parameters, only: drifts_implicit
      use parallel_streaming, only: get_zed_derivative_extended_domain, center_zed
      use parallel_streaming, only: gradpar_c, stream_sign
      use dist_fn_arrays, only: wdriftx_g, wdrifty_g
      use extended_zgrid, only: fill_zext_ghost_zones
      use extended_zgrid, only: map_to_iz_ikx_from_izext
      use extended_zgrid, only: periodic

      implicit none

      complex, dimension(:), intent(in) :: pdf
      integer, intent(in) :: ivmu, iky, ie
      complex, dimension(:), intent(out) :: rhs

      real, dimension(:), allocatable :: gradpar_fac
      complex, dimension(:), allocatable :: dpdf_dz
      real :: constant_factor
      integer :: iv, is, iz
      integer :: ia, ikx

      integer, dimension(:), allocatable :: iz_from_izext, ikx_from_izext
      integer :: nz_ext, izext
      complex :: pdf_left, pdf_right

      ia = 1
      iv = iv_idx(vmu_lo, ivmu)
      is = is_idx(vmu_lo, ivmu)

      ! nz_ext is the number of grid points in the extended zed domain
      nz_ext = size(pdf)

      ! determine the mapping from the extended domain zed index (izext) to the
      ! zed and kx domain indices (iz, ikx)
      allocate (iz_from_izext(nz_ext))
      allocate (ikx_from_izext(nz_ext))
      call map_to_iz_ikx_from_izext(iky, ie, iz_from_izext, ikx_from_izext)

      ! fill ghost zones beyond ends of extended zed domain for the pdf
      ! and store values in scratch_left and scratch_right
      call fill_zext_ghost_zones(iky, pdf, pdf_left, pdf_right)

      ! obtain the zed derivative of the pdf and store in dpdf_dz
      allocate (dpdf_dz(nz_ext))
      call get_zed_derivative_extended_domain(iv, pdf, pdf_left, pdf_right, dpdf_dz)

      ! compute the z-independent factor appearing in front of the d(pdf)/dz term on the RHS of the Gk equation
      constant_factor = -code_dt * spec(is)%stm_psi0 * vpa(iv) * time_upwind_minus
      ! use the correctly centred (b . grad z) pre-factor for this sign of vpa
      allocate (gradpar_fac(-nzgrid:nzgrid))
      if (stream_sign(iv) > 0) then
         gradpar_fac = gradpar_c(:, -1) * constant_factor
      else
         gradpar_fac = gradpar_c(:, 1) * constant_factor
      end if

      rhs = pdf
      if (drifts_implicit) then
         do izext = 1, nz_ext
            ikx = ikx_from_izext(izext)
            iz = iz_from_izext(izext)
            rhs(izext) = rhs(izext) * (1.0 + zi * time_upwind_minus &
                                       * (wdriftx_g(ia, iz, ivmu) * akx(ikx) + wdrifty_g(ia, iz, ivmu) * aky(iky)))
         end do
      end if

      ! cell-center the terms involving the pdf
      call center_zed(iv, rhs, 1, periodic(iky))

      ! construct the source term on the RHS of the GK equation coming from
      ! the pdf evaluated at the previous time level
      do izext = 1, nz_ext
         rhs(izext) = rhs(izext) + gradpar_fac(iz_from_izext(izext)) * dpdf_dz(izext)
      end do

      deallocate (dpdf_dz)
      deallocate (gradpar_fac)
      deallocate (iz_from_izext, ikx_from_izext)

   end subroutine get_contributions_from_pdf

   subroutine sweep_g_zext(iky, ie, it, ivmu, pdf)

      use constants, only: zi
      use zgrid, only: nzgrid, ntubes, delzed
      use run_parameters, only: drifts_implicit
      use run_parameters, only: zed_upwind_plus, zed_upwind_minus
      use run_parameters, only: time_upwind_plus
      use kt_grids, only: nakx, akx, aky
      use dist_fn_arrays, only: wdriftx_g, wdrifty_g
      use extended_zgrid, only: map_to_extended_zgrid
      use extended_zgrid, only: periodic, phase_shift
      use parallel_streaming, only: stream_sign, stream_c
      use parallel_streaming, only: center_zed
      use stella_layouts, only: vmu_lo, iv_idx, is_idx

      implicit none

      integer, intent(in) :: iky, ie, it, ivmu
      complex, dimension(:), intent(in out) :: pdf

      complex, dimension(:), allocatable :: wdrift_ext, pdf_cf
      complex, dimension(:, :), allocatable :: wdrift
      complex :: wd_factor, fac1, phase_factor
      real :: zupwnd_p, zupwnd_m, tupwnd_p
      real :: stream_term
      integer :: iz, ikx, ia
      integer :: iv, is
      integer :: ulim, sgn, iz1, iz2

      iv = iv_idx(vmu_lo, ivmu)
      is = is_idx(vmu_lo, ivmu)

      sgn = stream_sign(iv)
      ! avoid repeated calculation of constants
      zupwnd_p = 2.0 * zed_upwind_plus
      zupwnd_m = 2.0 * zed_upwind_minus
      tupwnd_p = 2.0 * time_upwind_plus

      ! if treating magentic drifts implicitly in time,
      ! get the drift frequency on the extended zed grid
      if (drifts_implicit) then
         allocate (wdrift(nakx, -nzgrid:nzgrid))
         ia = 1
         ! sum up the kx and ky contributions to the magnetic drift frequency
         do iz = -nzgrid, nzgrid
            do ikx = 1, nakx
               wdrift(ikx, iz) = -zi * (wdriftx_g(ia, iz, ivmu) * akx(ikx) + wdrifty_g(ia, iz, ivmu) * aky(iky))
            end do
         end do
         ! obtain the drift frequency on the extended zed domain
         allocate (wdrift_ext(size(pdf)))
         call map_to_extended_zgrid(it, ie, iky, spread(wdrift, 3, ntubes), wdrift_ext, ulim)
         ! NB: need to check if passing periodic(iky) is the right thing to do here
         call center_zed(iv, wdrift_ext, 1, periodic(iky))
      else
         ulim = size(pdf)
      end if
      ! determine the starting and ending indices for sweep over the extended zed grid.
      ! as we are using a zero-incoming BC, these indices depend on the sign of the advection velocity
      ! note that sgn < 0 actually corresponds to positive advection velocity
      if (sgn < 0) then
         iz1 = 1; iz2 = ulim
      else
         iz1 = ulim; iz2 = 1
      end if
      ! the case of periodic BC must be treated separately from the zero-incoming-BC case
      if (periodic(iky)) then
         ! to enforce periodicity, decompose the pdf into a particular integral
         ! and complementary function.

         ! calculate the particular integral, with zero BC, and store in pdf
         iz = sgn * nzgrid
         pdf(iz1) = 0.0
         call get_updated_pdf(iz, iv, is, sgn, iz1, iz2, wdrift_ext, pdf)
         ! calculate the complementary function, with unit BC, and store in pdf_cf
         allocate (pdf_cf(ulim))
         iz = sgn * nzgrid
         pdf_cf = 0.0; pdf_cf(iz1) = 1.0
         call get_updated_pdf(iz, iv, is, sgn, iz1, iz2, wdrift_ext, pdf_cf)
         ! construct pdf = pdf_PI + (pdf_PI(zend)/(1-pdf_CF(zend))) * pdf_CF
         phase_factor = phase_shift(iky)**(-sgn)
         pdf = pdf + (phase_factor * pdf(iz2) / (1.0 - phase_factor * pdf_cf(iz2))) * pdf_cf
         deallocate (pdf_cf)
      else
         ! specially treat the most upwind grid point
         iz = sgn * nzgrid
         wd_factor = 1.0
         if (drifts_implicit) wd_factor = 1.0 + 0.5 * tupwnd_p * wdrift_ext(iz1)
         stream_term = tupwnd_p * stream_c(iz, iv, is) / delzed(0)
         fac1 = zupwnd_p * wd_factor + sgn * stream_term
         pdf(iz1) = pdf(iz1) * 2.0 / fac1
         ! now that we have the pdf at the most upwind point, sweep over the
         ! rest of the extended zed domain to obtain the pdf(z)
         call get_updated_pdf(iz, iv, is, sgn, iz1, iz2, wdrift_ext, pdf)
      end if

      if (drifts_implicit) deallocate (wdrift, wdrift_ext)

   end subroutine sweep_g_zext

   subroutine get_updated_pdf(iz, iv, is, sgn, iz1, iz2, wdrift_ext, pdf)

      use zgrid, only: nzgrid, delzed
      use run_parameters, only: drifts_implicit
      use run_parameters, only: zed_upwind_plus, zed_upwind_minus
      use run_parameters, only: time_upwind_plus
      use parallel_streaming, only: stream_c

      implicit none

      integer, intent(in out) :: iz
      integer, intent(in) :: iv, is, sgn, iz1, iz2
      complex, dimension(:), intent(in) :: wdrift_ext
      complex, dimension(:), intent(in out) :: pdf

      integer :: izext
      real :: stream_term
      real :: tupwnd_p, zupwnd_p, zupwnd_m
      complex :: wd_factor, fac1, fac2

      tupwnd_p = 2.0 * time_upwind_plus
      zupwnd_p = 2.0 * zed_upwind_plus
      zupwnd_m = 2.0 * zed_upwind_minus

      ! wd_factor will be modified from below unity to account for magnetic drifts
      ! if the drifts are treated implicitly
      wd_factor = 1.0

      do izext = iz1 - sgn, iz2, -sgn
         if (iz == -sgn * nzgrid) then
            iz = sgn * nzgrid - sgn
         else
            iz = iz - sgn
         end if
         if (drifts_implicit) wd_factor = 1.0 + 0.5 * tupwnd_p * wdrift_ext(izext)
         stream_term = tupwnd_p * stream_c(iz, iv, is) / delzed(0)
         fac1 = zupwnd_p * wd_factor + sgn * stream_term
         fac2 = zupwnd_m * wd_factor - sgn * stream_term
         pdf(izext) = (-pdf(izext + sgn) * fac2 + 2.0 * pdf(izext)) / fac1
      end do

   end subroutine get_updated_pdf

   subroutine sweep_zed_zonal(iky, iv, is, sgn, g, llim)

      use zgrid, only: nzgrid, delzed
      use extended_zgrid, only: phase_shift
      use run_parameters, only: zed_upwind, time_upwind
      use parallel_streaming, only: stream_c

      implicit none

      integer, intent(in) :: iky, iv, is, sgn, llim
      complex, dimension(llim:), intent(in out) :: g

      integer :: iz, izext, iz1, iz2, npts, ulim
      real :: fac1, fac2
      complex :: pf
      complex, dimension(:), allocatable :: gcf, gpi

      npts = size(g)
      ulim = llim + npts - 1

      allocate (gpi(llim:ulim))
      allocate (gcf(llim:ulim))
      ! ky=0 is 2pi periodic (no extended zgrid)
      ! decompose into complementary function + particular integral
      ! zero BC for particular integral
      ! unit BC for complementary function (no source)
      if (sgn < 0) then
         iz1 = llim; iz2 = ulim
      else
         iz1 = ulim; iz2 = llim
      end if
      pf = phase_shift(iky)**(-sgn)
      gpi(iz1) = 0.; gcf(iz1) = 1.
      iz = sgn * nzgrid
      do izext = iz1 - sgn, iz2, -sgn
         if (iz == -sgn * nzgrid) then
            iz = sgn * nzgrid - sgn
         else
            iz = iz - sgn
         end if
         fac1 = 1.0 + zed_upwind + sgn * (1.0 + time_upwind) * stream_c(iz, iv, is) / delzed(0)
         fac2 = 1.0 - zed_upwind - sgn * (1.0 + time_upwind) * stream_c(iz, iv, is) / delzed(0)
         gpi(izext) = (-gpi(izext + sgn) * fac2 + 2.0 * g(izext)) / fac1
         gcf(izext) = -gcf(izext + sgn) * fac2 / fac1
      end do
      ! g = g_PI + (g_PI(pi)/(1-g_CF(pi))) * g_CF
      g = gpi + (pf * gpi(iz2) / (1.0 - pf * gcf(iz2))) * gcf
      deallocate (gpi, gcf)

   end subroutine sweep_zed_zonal

   subroutine invert_parstream_response(phi, apar)

      use linear_solve, only: lu_back_substitution
      use physics_flags, only: include_apar
      use zgrid, only: nzgrid, ntubes
      use extended_zgrid, only: neigen
      use extended_zgrid, only: nsegments
      use extended_zgrid, only: nzed_segment
      use extended_zgrid, only: map_to_extended_zgrid
      use extended_zgrid, only: map_from_extended_zgrid
      use extended_zgrid, only: ikxmod
      use extended_zgrid, only: periodic, phase_shift
      use kt_grids, only: naky
      use fields, only: nfields
      use fields_arrays, only: response_matrix

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in out) :: phi, apar

      integer :: iky, ie, it, ulim
      integer :: ikx
      integer :: nresponse_per_field, nresponse, nzp
      complex, dimension(:), allocatable :: fld_ext, phi_ext, apar_ext
      complex, dimension(:), allocatable :: fld

      ! need to put the fields into extended zed grid
      do iky = 1, naky
         ! avoid double counting of periodic endpoints for zonal (and any other periodic) modes
         if (periodic(iky)) then
            nzp = nzgrid - 1
            allocate (fld(nzp * nfields))
            do it = 1, ntubes
               do ie = 1, neigen(iky)
                  ikx = ikxmod(1, ie, iky)
                  ! construct the field vector, consisting of phi (and apar if evolving)
                  fld(:nzp) = phi(iky, ikx, :nzp, it)
                  if (include_apar) fld(nzp + 1:2 * nzp) = apar(iky, ikx, :nzp, it)
                  ! use LU back substitution to solve the linear response matrix system
                  call lu_back_substitution(response_matrix(iky)%eigen(ie)%zloc, &
                                            response_matrix(iky)%eigen(ie)%idx, fld)
                  ! unpack phi (and apar if evolving) from the field vector;
                  ! also apply phase shift at periodic point
                  phi(iky, ikx, :nzp, it) = fld(:nzp)
                  phi(iky, ikx, nzgrid, it) = phi(iky, ikx, -nzgrid, it) / phase_shift(iky)
                  if (include_apar) then
                     apar(iky, ikx, :nzp, it) = fld(nzp + 1:2 * nzp)
                     apar(iky, ikx, nzgrid, it) = apar(iky, ikx, -nzgrid, it) / phase_shift(iky)
                  end if
               end do
            end do
            deallocate (fld)
         else
            do it = 1, ntubes
               do ie = 1, neigen(iky)
                  ! solve response_matrix*(phi^{n+1}, apar^{n+1}) = (phi_{inh}^{n+1}, apar_{inh}^{n+1})
                  nresponse_per_field = nsegments(ie, iky) * nzed_segment + 1
                  nresponse = nresponse_per_field * nfields
                  ! fld_ext will contain phi or (phi, apar), depending on whether apar is evolved
                  allocate (fld_ext(nresponse))

                  ! phi_ext contains phi on the extended zed domain
                  allocate (phi_ext(nresponse_per_field))
                  call map_to_extended_zgrid(it, ie, iky, phi(iky, :, :, :), phi_ext, ulim)
                  ! include the phi_ext contribution to fld_ext
                  fld_ext(:nresponse_per_field) = phi_ext

                  ! if apar evolved, obtain apar on the extended zed domain (apar_ext)
                  ! and include its contribution to fld_ext
                  if (include_apar) then
                     allocate (apar_ext(nresponse_per_field))
                     call map_to_extended_zgrid(it, ie, iky, apar(iky, :, :, :), apar_ext, ulim)
                     fld_ext(nresponse_per_field + 1:2 * nresponse_per_field) = apar_ext
                  end if

                  ! use LU back substitution to solve linear response matrix system
                  call lu_back_substitution(response_matrix(iky)%eigen(ie)%zloc, &
                                            response_matrix(iky)%eigen(ie)%idx, fld_ext)

                  ! get phi_ext contribution from fld_ext and map back to normal (z, kx) grid
                  phi_ext = fld_ext(:nresponse_per_field)
                  call map_from_extended_zgrid(it, ie, iky, phi_ext, phi(iky, :, :, :))

                  ! if advancing apar, get apar_ext from fld_ext and map back to (z, kx) grid
                  if (include_apar) then
                     apar_ext = fld_ext(nresponse_per_field + 1:2 * nresponse_per_field)
                     call map_from_extended_zgrid(it, ie, iky, apar_ext, apar(iky, :, :, :))
                  end if

                  deallocate (fld_ext)
                  deallocate (phi_ext)
                  if (allocated(apar_ext)) deallocate (apar_ext)
               end do
            end do
         end if
      end do

   end subroutine invert_parstream_response

end module implicit_solve
