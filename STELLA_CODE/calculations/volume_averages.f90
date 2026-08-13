module volume_averages

   public :: init_volume_averages, finish_volume_averages
   public :: fieldline_average
   public :: eval_transit_ints
   public :: eval_transit_int_integrand_RH
   public :: eval_Q_fac
   public :: volume_average
   public :: flux_surface_average_ffs

   public :: mode_fac

   public :: jacobian_ky

   private

   interface fieldline_average
      module procedure fieldline_average_real_1d
      module procedure fieldline_average_real
      module procedure fieldline_average_complex
   end interface

   real, dimension(:), allocatable :: mode_fac
   !> Fourier coefficients in y of the Jacobian;
   !> needed for full flux surface simulations
   complex, dimension(:, :), allocatable :: jacobian_ky

contains

   subroutine init_volume_averages

      use zgrid, only: nzgrid, nztot, delzed
      use parameters_kxky_grids, only: nalpha, nakx, naky
      use grids_kxky, only: rho_d_clamped, aky
      use geometry, only: geo_surf, drhodpsip
      use geometry, only: geo_surf, jacob, djacdrho, q_as_x, dVolume
      use parameters_physics, only: full_flux_surface, radial_variation

      implicit none

      real :: dqdrho

      if (.not. allocated(mode_fac)) then
         allocate (mode_fac(naky)); mode_fac = 2.0
         if (aky(1) < epsilon(0.)) mode_fac(1) = 1.0
      end if
      if (full_flux_surface) then
         call init_flux_surface_average_ffs
      end if

      dqdrho = geo_surf%shat * geo_surf%qinp / geo_surf%rhoc
      if (.not. allocated(dVolume)) allocate (dVolume(nalpha, nakx, -nzgrid:nzgrid))

      !dVolume contains the volume element jacob, which may vary with x or alpha
      ! NB: dVolume does not contain the factor dx, as this should always be uniform
      dVolume = spread(jacob * spread(delzed, 1, nalpha), 2, nakx)
      if (q_as_x) then
         dVolume = dVolume / (dqdrho * drhodpsip)
      end if

      if (radial_variation) then
         if (q_as_x) then
            dVolume = dVolume * (1.+spread(spread(rho_d_clamped, 1, nalpha), 3, nztot) &
                                 * (spread(djacdrho / jacob, 2, nakx) - geo_surf%d2qdr2 / dqdrho &
                                    + geo_surf%d2psidr2 * drhodpsip))
         else
            dVolume = dVolume * (1.+spread(spread(rho_d_clamped, 1, nalpha), 3, nztot) &
                                 * spread(djacdrho / jacob, 2, nakx))
         end if
      end if

      if (full_flux_surface) then
         !something should go here
      end if

      !avoid the double counting at the zed boundaries
      dVolume(:, :, -nzgrid) = 0.5 * dVolume(:, :, -nzgrid)
      dVolume(:, :, nzgrid) = 0.5 * dVolume(:, :, nzgrid)

   end subroutine init_volume_averages

   subroutine finish_volume_averages

      use geometry, only: dVolume
      use parameters_physics, only: full_flux_surface

      implicit none

      if (allocated(mode_fac)) deallocate (mode_fac)
      if (allocated(dVolume)) deallocate (dVolume)
      if (full_flux_surface) then
         if (allocated(jacobian_ky)) deallocate (jacobian_ky)
      end if

   end subroutine finish_volume_averages

   !==============================================
   !============ FIELD LINE AVERAGE ==============
   !==============================================
   subroutine fieldline_average_real_1d(unavg, avg)

      use zgrid, only: nzgrid, ntubes
      use geometry, only: dl_over_b
      use species, only: nspec

      implicit none

      real, dimension(-nzgrid:, :), intent(in) :: unavg
      real, intent(out) :: avg

      integer :: it, ia

      ia = 1

      avg = 0.0
      do it = 1, ntubes
         avg = avg + sum( dl_over_b(ia, :) * unavg(:, it))
      end do
      avg = avg / real(ntubes)

   end subroutine fieldline_average_real_1d

   subroutine fieldline_average_real(unavg, avg)

      use zgrid, only: nzgrid, ntubes
      use parameters_kxky_grids, only: nakx, naky
      use geometry, only: dl_over_b

      implicit none

      real, dimension(:, :, -nzgrid:, :), intent(in) :: unavg
      real, dimension(:, :), intent(out) :: avg

      integer :: it, ia

      ia = 1

      avg = 0.0
      do it = 1, ntubes
         avg = avg + sum(spread(spread(dl_over_b(ia, :), 1, naky), 2, nakx) * unavg(:, :, :, it), dim=3)
      end do
      avg = avg / real(ntubes)

   end subroutine fieldline_average_real

   subroutine fieldline_average_complex(unavg, avg)

      use zgrid, only: nzgrid, ntubes
      use parameters_kxky_grids, only: nakx, naky
      use geometry, only: dl_over_b

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: unavg
      complex, dimension(:, :), intent(out) :: avg

      integer :: it, ia

      ia = 1

      avg = 0.0
      do it = 1, ntubes
         avg = avg + sum(spread(spread(dl_over_b(ia, :), 1, naky), 2, nakx) * unavg(:, :, :, it), dim=3)
      end do
      avg = avg / real(ntubes)

   end subroutine fieldline_average_complex

   !==============================================
   !============== VOLUME AVERAGE ================
   !==============================================
   subroutine volume_average(unavg, avg)

      use zgrid, only: nzgrid, ntubes
      use parameters_kxky_grids, only: naky, nakx
      use geometry, only: dl_over_b

      implicit none

      complex, dimension(:, :, -nzgrid:, :), intent(in) :: unavg
      real, intent(out) :: avg

      integer :: iky, ikx, iz, it, ia

      ia = 1

      avg = 0.
      do it = 1, ntubes
         do iz = -nzgrid, nzgrid
            do ikx = 1, nakx
               do iky = 1, naky
                  avg = avg + real(unavg(iky, ikx, iz, it) * conjg(unavg(iky, ikx, iz, it))) * mode_fac(iky) * dl_over_b(ia, iz)
               end do
            end do
         end do
      end do

      avg = avg / real(ntubes)

   end subroutine volume_average

   subroutine init_flux_surface_average_ffs

      use zgrid, only: nzgrid
      use parameters_kxky_grids, only: naky
      use extended_zgrid, only: periodic
      use geometry, only: jacob
      use stella_transforms, only: transform_alpha2kalpha

      implicit none

      integer :: iz

      if (.not. allocated(jacobian_ky)) allocate (jacobian_ky(naky, -nzgrid:nzgrid))

      !> calculate the Fourier coefficients in y of the Jacobian
      !> this is needed in the computation of the flux surface average of phi
      do iz = -nzgrid, nzgrid
         call transform_alpha2kalpha(jacob(:, iz), jacobian_ky(:, iz))
      end do
      where (periodic)
         jacobian_ky(:, nzgrid) = 0.0
      end where

   end subroutine init_flux_surface_average_ffs

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


   subroutine flux_surface_average_ffs(no_fsa, fsa)

      use zgrid, only: nzgrid, delzed 
      use parameters_kxky_grids, only: naky, naky_all 

      implicit none

      complex, dimension(:, -nzgrid:), intent(in) :: no_fsa
      complex, intent(out) :: fsa

      integer :: iky, ikymod
      complex :: area

      ! the normalising factor int dy dz Jacobian
      area = 0.0
      fsa = 0.0
      ! get contribution from negative ky values
      ! for no_fsa, iky=1 corresponds to -kymax, and iky=naky-1 to -dky
      do iky = 1, naky - 1
         ! jacobian_ky only defined for positive ky values
         ! use reality of the jacobian to fill in negative ky values
         ! i.e., jacobian_ky(-ky) = conjg(jacobian_ky(ky))
         ! ikymod runs from naky down to 2, which corresponds
         ! to ky values in jacobian_ky from kymax down to dky
         ikymod = naky - iky + 1
         ! for each ky, add the integral over zed
         fsa = fsa + sum(delzed * no_fsa(iky, :) * jacobian_ky(ikymod, :))
         area = area + sum(delzed * jacobian_ky(ikymod, :))
      end do
      ! get contribution from zero and positive ky values
      ! iky = naky correspond to ky=0 for no_fsa and iky=naky_all to ky=kymax
      do iky = naky, naky_all
         ! ikymod runs from 1 to naky
         ! ikymod = 1 corresponds to ky=0 for jacobian_ky and ikymod=naky to ky=kymax
         ikymod = iky - naky + 1
         ! for each ky, add the integral over zed
         fsa = fsa + sum(delzed * no_fsa(iky, :) * conjg(jacobian_ky(ikymod, :)))
         ! TODO-GA: WARNING: Possible change of value in conversion from COMPLEX(8) to REAL(8)
         area = area + sum(delzed * jacobian_ky(ikymod, :))
      end do
      ! normalise by the flux surface area
      fsa = fsa / area

   end subroutine flux_surface_average_ffs

end module volume_averages
