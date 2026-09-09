module dist_fn

  use debug_flags, only: debug => dist_fn_debug
  
   implicit none

   public :: init_gxyz
   public :: init_dist_fn, finish_dist_fn
   public :: checksum

   private

   interface checksum
      module procedure checksum_field
      module procedure checksum_dist
   end interface

   logical :: dist_fn_initialized = .false.
   logical :: gxyz_initialized = .false.
   logical :: kp2init = .false.
   logical :: dkp2drinit = .false.
   logical :: vp2init = .false.

contains

   subroutine init_gxyz(restarted)

      use arrays_dist_fn, only: gvmu, gold, gnew
      use redistribute, only: gather, scatter
      use mp, only: mp_abort
      use dist_redistribute, only: kxkyz2vmu
      use parameters_physics, only: radial_variation
      use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      use stella_transforms, only: transform_kx2x_xfirst, transform_x2kx_xfirst
      use parameters_kxky_grids, only: nalpha, nakx, naky, ikx_max
      use calculations_kxky, only: multiply_by_rho
      use vpamu_grids, only: mu, vpa, vperp2
      use vpamu_grids, only: maxwell_vpa, maxwell_mu, maxwell_fac
      use zgrid, only: nzgrid, ntubes, zed
      use species, only: spec, pfac, electron_species, nspec
      use geometry, only: dBdrho, gfac, geo_surf
      use geometry, only: PS_flow_fac, sym_flow_fac
      use geometry, only: PS_flow_defined, sym_flow_defined
      use gyro_averages, only: aj0x
      use arrays_dist_fn, only: kperp2
      use rosenbluth_hinton, only: RH_integrand_even, RH_integrand_odd, RH_inertia
      use parameters_physics, only: zonal_g_exb, zonal_nkx
      use parameters_physics, only: zonal_PS_fac, zonal_usym_fac
      use parameters_physics, only: zonal_rh_fac
      use parameters_physics, only: zonal_init_option_switch, zonal_init_none, zonal_init_triangular
      use parameters_physics, only: zonal_closure_option_switch, zonal_closure_rh, zonal_closure_flow
      use grids_kxky, only: akx
      use constants, only: zi

      implicit none

      real :: corr
      integer :: ivmu, is, imu, iv, it, iz, ia, ikx
      integer :: ikx_ZF, ikx_ZF_max, m_ZF
      real, dimension(:, :), allocatable :: energy
      complex, dimension(:, :), allocatable :: g0k
      complex, dimension(nakx) :: phi_ZF
      real, dimension(-nzgrid:nzgrid) :: u_parallel_ZF
      logical, intent(in) :: restarted

      !> The parallel flow the zonal Maxwellian carries, per unit dphi/dx.  Both
      !> profiles are available in any geometry; geometry warns if the symmetry
      !> one is being used where quasisymmetry does not hold.
      !>
      !> The Pfirsch-Schlueter weight belongs to the 'flow' closure alone, but
      !> the symmetry-direction one is available under 'rh' as well: a flow
      !> along the symmetry direction is divergence-free and carries no density,
      !> so it can be added to the relaxed Rosenbluth-Hinton state without
      !> disturbing the quasineutrality that state was built to satisfy.
      u_parallel_ZF = 0.

      !> Refuse rather than return zero: a geometry that cannot supply one of
      !> these profiles leaves its factor at zero, and a silently absent
      !> parallel flow looks exactly like a deliberate one.
      if (zonal_closure_option_switch == zonal_closure_flow) then
         if (abs(zonal_PS_fac) > epsilon(0.) .and. .not. PS_flow_defined) call mp_abort &
            ('the Pfirsch-Schlueter flow profile is not defined for this geometry; aborting')
         u_parallel_ZF = u_parallel_ZF + zonal_PS_fac * PS_flow_fac
      end if

      if (zonal_closure_option_switch == zonal_closure_flow .or. &
          zonal_closure_option_switch == zonal_closure_rh) then
         if (abs(zonal_usym_fac) > epsilon(0.) .and. .not. sym_flow_defined) call mp_abort &
            ('the symmetry-direction flow profile is not defined for this geometry; aborting')
         u_parallel_ZF = u_parallel_ZF + zonal_usym_fac * sym_flow_fac
      end if

      if (gxyz_initialized) return
      gxyz_initialized = .false.

      ! get version of g that has ky,kx,z local
      call gather(kxkyz2vmu, gvmu, gnew)

      ia = 1

      !calculate radial corrections to F0 for use in Krook operator, as well as g1 from initialization
      if (radial_variation) then
         !init_g uses maxwellians, so account for variation in temperature, density, and B

         allocate (energy(nalpha, -nzgrid:nzgrid))
         allocate (g0k(naky, nakx))

         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)
            energy = (vpa(iv)**2 + vperp2(:, :, imu)) * (spec(is)%temp_psi0 / spec(is)%temp)
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid

                  corr = -(pfac * (spec(is)%fprim + spec(is)%tprim * (energy(ia, iz) - 1.5)) &
                           + 2 * gfac * mu(imu) * dBdrho(iz))

                  if (.not. restarted) then
                     g0k = corr * gnew(:, :, iz, it, ivmu)
                     call multiply_by_rho(g0k)
                     gnew(:, :, iz, it, ivmu) = gnew(:, :, iz, it, ivmu) + g0k
                  end if
               end do
            end do
         end do
         deallocate (energy, g0k)

         if (.not. restarted) call scatter(kxkyz2vmu, gnew, gvmu)
      end if

      ! Treat zonal g differently for triangular ZF case
      ! We put this here because ginit uses layout in xyz
      if (zonal_init_option_switch /= zonal_init_none) then
         !> Prescribed zonal profile on kx = zonal_nkx * dkx, i.e. array index
         !> 1 + zonal_nkx.  Identical at every z and velocity-space point, so
         !> build it once here.  The harmonic needs a conjugate partner in the
         !> negative-kx half: past nakx-ikx_max+1 it is either the unpaired
         !> Nyquist mode, which cannot carry an imaginary amplitude, or already
         !> in that half, where the symmetry loop below overwrites it with zero.
         ikx_ZF_max = nakx - ikx_max + 1
         ikx_ZF = 1 + zonal_nkx
         if (zonal_nkx < 1 .or. ikx_ZF > ikx_ZF_max) call mp_abort &
            ('zonal_nkx must be at least 1 and no larger than nakx-ikx_max; aborting')
         phi_ZF(:) = 0
         phi_ZF(ikx_ZF) = zi*zonal_g_exb / akx(ikx_ZF)**2

         if (zonal_init_option_switch == zonal_init_triangular) then
            !> A triangular wave is the odd harmonics of the fundamental
            !> falling off as 1/m^3.
            do m_ZF = 3, nakx, 2
               ikx = 1 + m_ZF * zonal_nkx
               if (ikx > ikx_ZF_max) exit
               phi_ZF(ikx) = phi_ZF(ikx_ZF) / (akx(ikx)/akx(ikx_ZF))**3
            end do
         end if

         do ikx = 1, nakx - ikx_max
            phi_ZF(nakx - ikx + 1) = conjg(phi_ZF(ikx + 1))
         end do

         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            is = is_idx(vmu_lo, ivmu)
            imu = imu_idx(vmu_lo, ivmu)
            iv = iv_idx(vmu_lo, ivmu)
            do it = 1, ntubes
               do iz = -nzgrid, nzgrid
                  do ikx = 2, nakx
                     if (spec(is)%type == electron_species) then
                        gnew(1, ikx, iz, it, ivmu) = 0.

                     else

                        if (zonal_closure_option_switch == zonal_closure_rh) then
                           !> Rosenbluth-Hinton profile, with <H> = F_M * I_RH to satisfy
                           !> quasineutrality, plus the symmetry-direction parallel flow of
                           !> <zonal_usym_fac>, which is v_parallel-odd and so leaves that
                           !> quasineutrality untouched.  Zero by default.
                           gnew(1, ikx, iz, it, ivmu) = spec(is)%z * phi_ZF(ikx) &
                              * (2*(1-aj0x(1,ikx,iz,ivmu)) + zonal_rh_fac*(aj0x(1,ikx,iz,ivmu)-2 &
                                   + conjg(RH_integrand_even(ikx,iz,it,ivmu)+RH_integrand_odd(ikx,iz,it,ivmu)) &
                                   + real(RH_inertia(ikx,iz,it,is))) &
                                 - aj0x(1,ikx,iz,ivmu)*zi*akx(ikx)*vpa(iv)*u_parallel_ZF(iz) ) &
                              * maxwell_mu(ia, iz, imu, is) * maxwell_vpa(iv, is) * maxwell_fac(is)

                        else if (zonal_closure_option_switch == zonal_closure_flow) then
                           !> A Maxwellian carrying a parallel flow.  See
                           !> <zonal_PS_fac> and <zonal_usym_fac> in
                           !> parameters_physics for what the two weights select.
                           gnew(1, ikx, iz, it, ivmu) = spec(is)%z * phi_ZF(ikx) * ( 2*(1-aj0x(1,ikx,iz,ivmu) ) &
                               - aj0x(1,ikx,iz,ivmu)*zi*akx(ikx)*vpa(iv)*u_parallel_ZF(iz) ) &
                              * maxwell_mu(ia, iz, imu, is) * maxwell_vpa(iv, is) * maxwell_fac(is)

                        else
                           ! Choose g ~ phi*(1-J0) to satisfy quasineutrality
                           gnew(1, ikx, iz, it, ivmu) = spec(is)%z * phi_ZF(ikx) * 2*(1-aj0x(1,ikx,iz,ivmu) ) &
                              * maxwell_mu(ia, iz, imu, is) * maxwell_vpa(iv, is) * maxwell_fac(is)
                        end if

                     end if
                  end do


               end do
            end do
         end do
      end if
 
      gold = gnew

   end subroutine init_gxyz

   subroutine init_dist_fn

      use mp, only: proc0
      use parameters_physics, only: radial_variation
      use stella_layouts, only: init_dist_fn_layouts
      use gyro_averages, only: init_bessel

      implicit none

      if (dist_fn_initialized) return
      dist_fn_initialized = .true.

      debug = debug .and. proc0

      if (debug) write (*, *) 'dist_fn::init_dist_fn::allocate_arrays'
      call allocate_arrays

      !> allocate and initialise kperp2 and dkperp2dr
      if (debug) write (*, *) 'dist_fn::init_dist_fn::init_kperp2'
      call init_kperp2
      if (radial_variation) call init_dkperp2dr

      !> allocate and initialise vperp2
      if (debug) write (*, *) 'dist_fn::init_dist_fn::init_vperp2'
      call init_vperp2

      !> init_bessel sets up arrays needed for gyro-averaging;
      !> for a flux tube simulation, this is j0 and j1;
      !> for a flux annulus simulation, gyro-averaging is non-local in ky
      !> and so more effort is required
      if (debug) write (*, *) 'dist_fn::init_dist_fn::init_bessel'
      call init_bessel

   end subroutine init_dist_fn

   !> init_kperp2 allocates and initialises the kperp2 array
   subroutine init_kperp2

      use arrays_dist_fn, only: kperp2
      use geometry, only: gds2, gds21, gds22
      use geometry, only: geo_surf, q_as_x
      use zgrid, only: nzgrid
      use parameters_kxky_grids, only: naky, nakx, nalpha
      use grids_kxky, only: akx, aky, theta0
      use grids_kxky, only: zonal_mode

      implicit none

      integer :: iky, ikx

      if (kp2init) return
      kp2init = .true.

      !> allocate the kperp2 array to contain |k_perp|^2
      allocate (kperp2(naky, nakx, nalpha, -nzgrid:nzgrid))

      do iky = 1, naky
         if (zonal_mode(iky)) then
            do ikx = 1, nakx
               if (q_as_x) then
                  kperp2(iky, ikx, :, :) = akx(ikx) * akx(ikx) * gds22
               else
                  kperp2(iky, ikx, :, :) = akx(ikx) * akx(ikx) * gds22 / (geo_surf%shat**2)
               end if
            end do
         else
            do ikx = 1, nakx
               kperp2(iky, ikx, :, :) = aky(iky) * aky(iky) &
                                        * (gds2 + 2.0 * theta0(iky, ikx) * gds21 &
                                           + theta0(iky, ikx) * theta0(iky, ikx) * gds22)
            end do
         end if
      end do

      ! NB: should really avoid this by using higher resolution when reading in VMEC geometry and then
      ! NB: course-graining if necessary to map onto lower-resolution stella grid
      ! ensure kperp2 is positive everywhere (only might go negative if using full-flux-surface due to interpolation)
      where (kperp2 < 0.0)
         kperp2 = 0.0
      end where

      call enforce_single_valued_kperp2

   end subroutine init_kperp2

   !> init_dkperp2dr allocates and initialises the dkperp2dr array, needed for radial variation
   subroutine init_dkperp2dr

      use arrays_dist_fn, only: kperp2, dkperp2dr
      use geometry, only: dgds2dr, dgds21dr, dgds22dr
      use geometry, only: geo_surf, q_as_x
      use zgrid, only: nzgrid
      use parameters_kxky_grids, only: naky, nakx, nalpha
      use grids_kxky, only: akx, aky, theta0
      use grids_kxky, only: zonal_mode

      implicit none

      integer :: iky, ikx

      if (dkp2drinit) return
      dkp2drinit = .true.

      allocate (dkperp2dr(naky, nakx, nalpha, -nzgrid:nzgrid))
      do iky = 1, naky
         if (zonal_mode(iky)) then
            do ikx = 1, nakx
               if (q_as_x) then
                  where (kperp2(iky, ikx, :, :) > epsilon(0.0))
                     dkperp2dr(iky, ikx, :, :) = akx(ikx) * akx(ikx) * dgds22dr / kperp2(iky, ikx, :, :)
                  elsewhere
                     dkperp2dr(iky, ikx, :, :) = 0.0
                  end where
               else
                  where (kperp2(iky, ikx, :, :) > epsilon(0.0))
                     dkperp2dr(iky, ikx, :, :) = akx(ikx) * akx(ikx) * dgds22dr / (geo_surf%shat**2 * kperp2(iky, ikx, :, :))
                  elsewhere
                     dkperp2dr(iky, ikx, :, :) = 0.0
                  end where
               end if
            end do
         else
            do ikx = 1, nakx
               dkperp2dr(iky, ikx, :, :) = aky(iky) * aky(iky) &
                                           * (dgds2dr + 2.0 * theta0(iky, ikx) * dgds21dr &
                                              + theta0(iky, ikx) * theta0(iky, ikx) * dgds22dr)
               dkperp2dr(iky, ikx, :, :) = dkperp2dr(iky, ikx, :, :) / kperp2(iky, ikx, :, :)
               if (any(kperp2(iky, ikx, :, :) < epsilon(0.))) dkperp2dr(iky, ikx, :, :) = 0.
            end do
         end if
      end do

   end subroutine init_dkperp2dr

   subroutine enforce_single_valued_kperp2

      use arrays_dist_fn, only: kperp2
      use parameters_kxky_grids, only: naky, nalpha
      use zgrid, only: nzgrid
      use extended_zgrid, only: neigen, nsegments, ikxmod

      implicit none

      integer :: iky, ie, iseg
      real, dimension(:), allocatable :: tmp

      allocate (tmp(nalpha)); tmp = 0.0

      do iky = 1, naky
         do ie = 1, neigen(iky)
            if (nsegments(ie, iky) > 1) then
               do iseg = 2, nsegments(ie, iky)
                  tmp = 0.5 * (kperp2(iky, ikxmod(iseg - 1, ie, iky), :, nzgrid) + kperp2(iky, ikxmod(iseg, ie, iky), :, -nzgrid))
                  kperp2(iky, ikxmod(iseg, ie, iky), :, -nzgrid) = tmp
                  kperp2(iky, ikxmod(iseg - 1, ie, iky), :, nzgrid) = tmp
               end do
            end if
         end do
      end do

      deallocate (tmp)

   end subroutine enforce_single_valued_kperp2

   subroutine allocate_arrays

      use stella_layouts, only: kxkyz_lo, vmu_lo, kymus_lo
      use zgrid, only: nzgrid, ntubes
      use parameters_kxky_grids, only: naky, nakx
      use vpamu_grids, only: nvpa, nmu
      use arrays_dist_fn, only: gnew, gold, g_scratch
      use arrays_dist_fn, only: gvmu, g_kymus
      use parameters_numerical, only: split_parallel_dynamics
      
      implicit none

      if (.not. allocated(gnew)) &
         allocate (gnew(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      gnew = 0.
      if (.not. allocated(gold)) &
         allocate (gold(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      gold = 0.
      if (.not. allocated(g_scratch)) &
         allocate (g_scratch(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      g_scratch = 0.
      if (.not. allocated(gvmu)) &
         allocate (gvmu(nvpa, nmu, kxkyz_lo%llim_proc:kxkyz_lo%ulim_alloc))
      gvmu = 0.
      if (.not. allocated(g_kymus)) then
         if (.not. split_parallel_dynamics) then
            allocate (g_kymus(nakx, -nzgrid:nzgrid, ntubes, nvpa, kymus_lo%llim_proc:kymus_lo%ulim_alloc))
         else
            allocate (g_kymus(1, 1, 1, 1, 1))
         end if
         g_kymus = 0.
      end if

   end subroutine allocate_arrays

   subroutine init_vperp2

      use geometry, only: bmag
      use zgrid, only: nzgrid
      use vpamu_grids, only: vperp2
      use vpamu_grids, only: nmu, mu
      use parameters_kxky_grids, only: nalpha

      implicit none

      integer :: imu

      if (vp2init) return
      vp2init = .true.

      if (.not. allocated(vperp2)) allocate (vperp2(nalpha, -nzgrid:nzgrid, nmu)); vperp2 = 0.

      do imu = 1, nmu
         vperp2(:, :, imu) = 2.0 * mu(imu) * bmag
      end do

   end subroutine init_vperp2

   subroutine finish_dist_fn

      use gyro_averages, only: finish_bessel

      implicit none

      call finish_bessel
      call finish_kperp2
      call finish_vperp2
      call deallocate_arrays

      dist_fn_initialized = .false.
      gxyz_initialized = .false.

   end subroutine finish_dist_fn

   subroutine deallocate_arrays

      use arrays_dist_fn, only: gnew, gold, g_scratch, gvmu, g_kymus

      implicit none

      if (allocated(gnew)) deallocate (gnew)
      if (allocated(gold)) deallocate (gold)
      if (allocated(g_scratch)) deallocate (g_scratch)
      if (allocated(gvmu)) deallocate (gvmu)
      if (allocated(g_kymus)) deallocate (g_kymus)

   end subroutine deallocate_arrays

   subroutine finish_kperp2

      use arrays_dist_fn, only: kperp2, dkperp2dr

      implicit none

      if (allocated(kperp2)) deallocate (kperp2)
      if (allocated(dkperp2dr)) deallocate (dkperp2dr)

      kp2init = .false.
      dkp2drinit = .false.

   end subroutine finish_kperp2

   subroutine finish_vperp2

      use vpamu_grids, only: vperp2

      implicit none

      if (allocated(vperp2)) deallocate (vperp2)

      vp2init = .false.

   end subroutine finish_vperp2

   subroutine checksum_field(field, total)

      use zgrid, only: nzgrid, ntubes
      use parameters_kxky_grids, only: naky
      use extended_zgrid, only: neigen, nsegments, ikxmod
      use extended_zgrid, only: iz_low, iz_up

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: field
      real, intent(out) :: total

      integer :: it, iky, ie, iseg
      integer :: ikx

      total = 0.

      do iky = 1, naky
         do it = 1, ntubes
            do ie = 1, neigen(iky)
               iseg = 1
               ikx = ikxmod(iseg, ie, iky)
               total = total + sum(cabs(field(iky, ikx, iz_low(iseg):iz_up(iseg), it)))
               if (nsegments(ie, iky) > 1) then
                  do iseg = 2, nsegments(ie, iky)
                     ikx = ikxmod(iseg, ie, iky)
                     total = total + sum(cabs(field(iky, ikx, iz_low(iseg) + 1:iz_up(iseg), it)))
                  end do
               end if
            end do
         end do
      end do

   end subroutine checksum_field

   subroutine checksum_dist(dist, total, norm)

      use mp, only: sum_allreduce
      use zgrid, only: nzgrid, ntubes
      use stella_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      use parameters_kxky_grids, only: naky, nakx
      use vpamu_grids, only: maxwell_vpa, maxwell_mu

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: dist
      real, intent(out) :: total
      logical, intent(in), optional :: norm

      integer :: ivmu, iv, imu, is
      integer :: iky, ikx, it
      real :: subtotal

      complex, dimension(:, :, :, :), allocatable :: dist_single

      total = 0.

      allocate (dist_single(naky, nakx, -nzgrid:nzgrid, ntubes))
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         dist_single = dist(:, :, :, :, ivmu)
         if (present(norm)) then
            if (norm) then
               iv = iv_idx(vmu_lo, ivmu)
               imu = imu_idx(vmu_lo, ivmu)
               is = is_idx(vmu_lo, ivmu)
               do it = 1, ntubes
                  do ikx = 1, nakx
                     do iky = 1, naky
                        dist_single(iky, ikx, :, it) = dist_single(iky, ikx, :, it) * maxwell_vpa(iv, is) * maxwell_mu(1, :, imu, is)
                     end do
                  end do
               end do
            else
            end if
         end if
         call checksum(dist_single, subtotal)
         total = total + subtotal
      end do
      deallocate (dist_single)

      call sum_allreduce(total)

   end subroutine checksum_dist

end module dist_fn
