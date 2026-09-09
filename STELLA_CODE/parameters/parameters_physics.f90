!###############################################################################
!############################ READ PHYSICS PARAMETES ###########################
!###############################################################################
! Namelist: &parameters_physics
! These are different logicals and parameters that adjust the physics in the 
! problem. These will allow you to toggles whether you want to include different
! terms in the gyrokinetic equation, as well as allowing you to include different
! large scale effects such as whether the system allows for electromanetic 
! effect, full flux surface effect, or radially global effects.
!###############################################################################

module parameters_physics

   implicit none

   !> Public subroutines that are read by the main stella routine.
   public :: read_parameters_physics
   public :: finish_read_parameters_physics

   !> Available physics options: These are standard gyrokinetic terms that
   !> can be turned on/off with the following toggles.
   public :: include_parallel_streaming
   public :: include_mirror
   public :: nonlinear
   public :: xdriftknob, ydriftknob, wstarknob
   
   !> Adiabatic options: This is used when nspec = 1. The non-kinetic
   !> species (usually electrons) is set to have an adiabatic response.
   !> This can be either the classic adiabatic option, or the modified
   !> adiabatic option (i.e. modified Boltzmann electrons).
   public :: adiabatic_option_switch, adiabatic_option_fieldlineavg
   public :: zonal_init_option_switch, zonal_closure_option_switch
   public :: zonal_init_none, zonal_init_cosine, zonal_init_triangular
   public :: zonal_closure_density, zonal_closure_rh, zonal_closure_flow
   
   !> Additional physics effects
   public :: prp_shear_enabled
   public :: hammett_flow_shear
   public :: include_pressure_variation
   public :: include_geometric_variation
   public :: include_parallel_nonlinearity
   public :: suppress_zonal_interaction
   public :: only_zonal_interaction
   public :: freeze_nonzonal
   public :: freeze_zonal
   public :: freeze_zonal_factor
   public :: freeze_zonal_kmin
   public :: freeze_zonal_kmax
   public :: RH_analytic_drift_phase, RH_analytic_drift_phase_specified
   public :: zonal_PS_fac, zonal_usym_fac
   public :: zonal_rh_fac
   public :: zonal_g_exb
   public :: zonal_nkx
   
   !> Large scale physics options of the system - e.g. whether we have full flux effects, 
   !> electromagnetic effects, or radially global effects.
   public :: full_flux_surface
   public :: include_apar
   public :: include_bpar
   public :: radial_variation
   
   public :: beta, zeff, tite, nine, rhostar, vnew_ref
   public :: g_exb, g_exbfac, omprimfac, omprimfac_RH, omprimfac_PS
   
   private

   logical :: include_parallel_streaming
   logical :: include_mirror
   logical :: nonlinear
   real :: xdriftknob, ydriftknob, wstarknob
 
   !> How the zonal profile is launched, and what distribution is put under it.
   character(20) :: zonal_init_option, zonal_closure_option
   integer :: zonal_init_option_switch, zonal_closure_option_switch
   integer, parameter :: zonal_init_none = 1, &
                         zonal_init_cosine = 2, &
                         zonal_init_triangular = 3
   integer, parameter :: zonal_closure_density = 1, &
                         zonal_closure_rh = 2, &
                         zonal_closure_flow = 3

   integer :: adiabatic_option_switch
   integer, parameter :: adiabatic_option_periodic = 1, &
                       adiabatic_option_zero = 2, &
                       adiabatic_option_fieldlineavg = 3
 
 
   logical :: prp_shear_enabled
   logical :: hammett_flow_shear 
   logical :: include_pressure_variation 
   logical :: include_geometric_variation
   logical :: include_parallel_nonlinearity
   logical :: suppress_zonal_interaction
   logical :: only_zonal_interaction
   logical :: freeze_nonzonal
   logical :: freeze_zonal

   !> How the Rosenbluth-Hinton drift-orbit phase Q is obtained.  True takes the
   !> closed form the theory gives for a quasisymmetric field; false integrates
   !> Q along the field line from the magnetic drifts.  The two agree in a
   !> tokamak, where the closed form applies.
   !>
   !> Left unset it follows the geometry: analytic where the geometry supplies
   !> the closed form, which is Miller, and numerical where it does not, which is
   !> VMEC.  <RH_analytic_drift_phase_specified> records whether the input file
   !> asked for a particular one, so that default can be applied without
   !> overriding the user.
   logical :: RH_analytic_drift_phase
   logical :: RH_analytic_drift_phase_specified
   !> Initialise the zonal distribution as a Maxwellian carrying a parallel
   !> flow, u_par = zonal_PS_fac * PS_flow_fac
   !>             + zonal_usym_fac * sym_flow_fac,
   !> the two profiles being the Pfirsch-Schlueter return flow and the flow
   !> along the direction of symmetry (see geometry).  Together they span every
   !> divergence-free parallel flow that can accompany the ExB flow, so the two
   !> scalars reach any of them:
   !>
   !>     (0, 0)  a density perturbation with no flow
   !>     (1, 0)  pure Pfirsch-Schlueter flow          (the default)
   !>     (0, 1)  flow along the symmetry direction, toroidal in a tokamak
   !>
   !> Both profiles are built by the geometry module and hold at finite aspect
   !> ratio; 2 q cos(theta) and a constant are their large-aspect-ratio limits.
   real :: zonal_PS_fac, zonal_usym_fac
   
   logical :: full_flux_surface
   logical :: include_apar
   logical :: include_bpar
   logical :: radial_variation

   real :: beta, zeff, tite, nine, rhostar, irhostar, vnew_ref
   real :: g_exb, g_exbfac, omprimfac, omprimfac_RH, omprimfac_PS
   real :: zonal_g_exb, zonal_rh_fac 
   !> radial harmonic index carrying the prescribed zonal profile.  1 (the
   !> default) puts it on the lowest kx, so the zonal wavelength equals the
   !> box length and CANNOT be varied independently of Lx.  Setting it to n
   !> puts the profile on kx = n*dkx, so raising jtwist and n together holds
   !> the zonal wavelength -- hence u_Z(0) and the flow shear -- fixed while
   !> the box grows.  That is what a box-length convergence test needs.
   integer :: zonal_nkx
   real :: freeze_zonal_factor, freeze_zonal_kmin, freeze_zonal_kmax
   logical :: initialised = .false.

   !!> Need to fix for the warning messages
   logical :: debug = .false.

contains

  !======================================================================
  !====================== READ PHYSICS PARAMETERS =======================
  !======================================================================
  subroutine read_parameters_physics

   use mp, only: proc0
   use text_options, only: text_option, get_option_value
   use file_utils, only: input_unit, error_unit, input_unit_exist

   implicit none
   
   character(30) :: adiabatic_option

   if (initialised) return

   if (proc0) call set_default_parameters
   if (proc0) call read_input_file
   call broadcast_parameters

   initialised = .true.

 contains 
   
   !**********************************************************************
   !                        SET DEFAULT PARAMETERS                       !
   !**********************************************************************
   ! If not specified in the input file these are the default options that 
   ! will be set for all parameters under the namelist 
   ! &parameters_physics'.
   !**********************************************************************
   subroutine set_default_parameters

      implicit none 

      !> Standard gyrokinetic terms
      include_parallel_streaming = .true.
      include_mirror = .true.
      nonlinear = .false.
      xdriftknob = 1.0
      ydriftknob = 1.0
      wstarknob = 1.0

      !> If not chose we set adiabatic option to be adiabatic electrons (no modified Boltzmann response)
      adiabatic_option = 'field-line-average-term'

      !> Additional effects that can be included but are not by default
      prp_shear_enabled = .false.
      hammett_flow_shear = .true.
      include_pressure_variation = .false.
      include_geometric_variation = .true.
      include_parallel_nonlinearity = .false.
      suppress_zonal_interaction = .false.
      only_zonal_interaction = .false.
      freeze_nonzonal = .false.
      freeze_zonal = .false.
      freeze_zonal_factor = 1.0
      freeze_zonal_kmin = -1.0
      freeze_zonal_kmax = 1e10
      zonal_init_option = 'default'
      zonal_closure_option = 'default'
      RH_analytic_drift_phase = .true.
      zonal_PS_fac  = 1.0
      zonal_usym_fac = 0.0
      zonal_g_exb    = 0.0
      zonal_nkx      = 1
      zonal_rh_fac = 1.0
      
      full_flux_surface = .false.
      include_apar = .false.
      include_bpar = .false.
      radial_variation = .false.

      beta = 0.0 ! beta = 8 * pi * p_ref / B_ref^2
      zeff = 1.0
      tite = 1.0
      nine = 1.0
      rhostar = -1.0 ! = m_ref * vt_ref / (e * B_ref * a_ref), with refs in SI
      vnew_ref = -1.0 ! various input options will override this value if it is negative

      !> Zonal flow options -> TODO-HT: how to turn on/off
      g_exb = 0.0          ! ExB shear
      g_exbfac = 1.0       ! Scale factor for perp. flow shear
      omprimfac = 1.0      ! Scale factor for equilibrium parallel flow shear
      omprimfac_RH = 0.0   ! Scale factor for Rosenbluth-Hinton "parallel flow" (velocity-space dep.)
      omprimfac_PS = 0.0   ! Scale factor for Pfirsch-Schlueter parallel flow (~q*cos(theta))
      irhostar = -1.0 
      
   end subroutine set_default_parameters

   !**********************************************************************
   !                         READ INPUT OPTIONS                          !
   !**********************************************************************
   ! Overwrite any default options with those specified in the input file. 
   ! Then change the other parameters consistently.
   !**********************************************************************
   subroutine read_input_file

      use file_utils, only: input_unit_exist, error_unit

      implicit none

      type(text_option), dimension(6), parameter :: adiabaticopts = &
      (/text_option('default', adiabatic_option_fieldlineavg), &
      !> TODO-HT or TODO-GA: sed: adiabatic_option_default -> adiabatic_option_periodic
      text_option('no-field-line-average-term', adiabatic_option_periodic), &
      text_option('field-line-average-term', adiabatic_option_fieldlineavg), &
      text_option('iphi00=0', adiabatic_option_periodic), &
      text_option('iphi00=1', adiabatic_option_periodic), &
      text_option('iphi00=2', adiabatic_option_fieldlineavg)/)

      type(text_option), dimension(4), parameter :: zonalinitopts = &
      (/text_option('default', zonal_init_none), &
      text_option('none', zonal_init_none), &
      text_option('cosine', zonal_init_cosine), &
      text_option('triangular', zonal_init_triangular)/)

      type(text_option), dimension(4), parameter :: zonalclosureopts = &
      (/text_option('default', zonal_closure_rh), &
      text_option('density', zonal_closure_density), &
      text_option('rh', zonal_closure_rh), &
      text_option('flow', zonal_closure_flow)/)

      integer :: ierr, in_file
      logical :: nml_exist
      logical :: probe_analytic_drift_phase

      namelist /parameters_physics/ include_parallel_streaming, include_mirror, nonlinear, &
        xdriftknob, ydriftknob, wstarknob, adiabatic_option, prp_shear_enabled, &
        hammett_flow_shear, include_pressure_variation, include_geometric_variation, &
        include_parallel_nonlinearity, suppress_zonal_interaction, only_zonal_interaction, freeze_nonzonal, freeze_zonal, &
        freeze_zonal_factor, freeze_zonal_kmin, freeze_zonal_kmax, &
        zonal_init_option, zonal_closure_option, &
        zonal_g_exb, &
        RH_analytic_drift_phase, &
        zonal_PS_fac, zonal_usym_fac, zonal_rh_fac, zonal_nkx, &
        full_flux_surface, include_apar, include_bpar, radial_variation, &
        beta, zeff, tite, nine, rhostar, vnew_ref, &
        g_exb, g_exbfac, omprimfac, omprimfac_RH, omprimfac_PS, irhostar
        
     !> Overwrite the default options with any that are explicitly given in the input file
     !> under the heading '&parameters_physics'
     in_file = input_unit_exist("parameters_physics", nml_exist)
     if (nml_exist) read (unit=in_file, nml=parameters_physics)

     !> Read the namelist a second time with the opposite default, to find out
     !> whether the input file mentioned RH_analytic_drift_phase at all.  If it
     !> did, both reads return the file's value; if it did not, they return the
     !> two different defaults.  Everything else keeps the value it already has,
     !> so the second read changes nothing but this.
     RH_analytic_drift_phase_specified = .false.
     if (nml_exist) then
        probe_analytic_drift_phase = RH_analytic_drift_phase
        RH_analytic_drift_phase = .not. probe_analytic_drift_phase
        rewind (in_file)
        read (unit=in_file, nml=parameters_physics)
        RH_analytic_drift_phase_specified = (RH_analytic_drift_phase .eqv. probe_analytic_drift_phase)
        RH_analytic_drift_phase = probe_analytic_drift_phase
     end if

     call check_backwards_compatability

     if (irhostar > 0) rhostar = 1./irhostar
     !> Don't allow people to set rhostar when its not full flux                                                                                                                                        
     !> Otherwise phase_shift_angle will be changed in grids_kxky.f90
     if (.not. full_flux_surface) rhostar = 0

     ierr = error_unit()
     call get_option_value &
       (adiabatic_option, adiabaticopts, adiabatic_option_switch, &
         ierr, "adiabatic_option in parameters_physics")
     call get_option_value &
       (zonal_init_option, zonalinitopts, zonal_init_option_switch, &
         ierr, "zonal_init_option in parameters_physics")
     call get_option_value &
       (zonal_closure_option, zonalclosureopts, zonal_closure_option_switch, &
         ierr, "zonal_closure_option in parameters_physics")

   end subroutine


   !**********************************************************************
   !                    CHECK BACKWARDS COMPATIBILITY                    !
   !**********************************************************************
   ! Make sure stella either runs or aborts old names for variables or
   ! namelists are used
   !**********************************************************************
   subroutine check_backwards_compatability

      use mp, only: mp_abort, broadcast
      use debug_flags, only: const_alpha_geo
      implicit none

      logical :: old_nml_exist
      integer :: in_file
      logical :: probe_analytic_drift_phase

      ! These variables belonged to <time_advance_knobs> and are now read in <run_parameters>
      ! We define them here so we can read the namelist, but we will not use them.
      character(10) :: explicit_option
      logical :: flip_flop
      
      namelist /physics_flags/ full_flux_surface, radial_variation, &
         include_parallel_nonlinearity, include_parallel_streaming, &
         include_mirror, include_apar, include_bpar, nonlinear, &
         include_pressure_variation, include_geometric_variation, &
         adiabatic_option, const_alpha_geo, suppress_zonal_interaction, only_zonal_interaction, &
         freeze_nonzonal, freeze_zonal, freeze_zonal_factor, freeze_zonal_kmin, freeze_zonal_kmax, &
         zonal_init_option, zonal_closure_option, &
         zonal_g_exb, &
         RH_analytic_drift_phase, &
         zonal_PS_fac, zonal_usym_fac, zonal_rh_fac, zonal_nkx

      namelist /parameters/ beta, zeff, tite, nine, rhostar, vnew_ref, &
         g_exb, g_exbfac, omprimfac, omprimfac_RH, omprimfac_PS, irhostar

      namelist /time_advance_knobs/ xdriftknob, ydriftknob, wstarknob, explicit_option, flip_flop
      
      in_file = input_unit_exist("physics_flags", old_nml_exist)
      if (old_nml_exist) then
         read (unit=in_file, nml=physics_flags)
         !> The deprecated namelist is read after the probe in <read_parameters>,
         !> so repeat the probe here.  Without it a RH_analytic_drift_phase set
         !> under <physics_flags> would read as unspecified and be overridden by
         !> the geometry default, silently ignoring the user.
         probe_analytic_drift_phase = RH_analytic_drift_phase
         RH_analytic_drift_phase = .not. probe_analytic_drift_phase
         rewind (in_file)
         read (unit=in_file, nml=physics_flags)
         if (RH_analytic_drift_phase .eqv. probe_analytic_drift_phase) &
            RH_analytic_drift_phase_specified = .true.
         RH_analytic_drift_phase = probe_analytic_drift_phase
         if(debug) then 
            write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(*,*) 'Please change the namelist <phyiscs_flags> in the input file'
            write(*,*) 'to <parameters_physics>. You can inlclude the old flags under'
            write(*,*) 'this new namelist'
            write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         end if
!         call broadcast(const_alpha_geo) 
         !         write(*,*) "Aborting in parameters_physics.f90.&
         !              The namelist <physics_flags> does not exist. &
         !      Please replace this with the title <parameters_physics>"
         ! call mp_abort('Aborting in parameters_physics.f90.& 
         !      The namelist <physics_flags> does not exist. &
         !      Please replace this with the title <parameters_physics>')
      end if
      in_file = input_unit_exist("parameters", old_nml_exist)
      if (old_nml_exist) then
         read (unit=in_file, nml=parameters)
         if(debug) then
            write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(*,*) 'Please change the namelist <parameters> in the input file'
            write(*,*) 'to <parameters_physics>. You can inlclude the old flags under'
            write(*,*) 'this new namelist'
            write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         end if
         ! write(*,*) "Aborting in parameters_physics.f90.&
         !      The namelist <physics_parameters> does not exist. &
         !      Please replace this with the title <parameters_physics>"
         ! call mp_abort("Aborting in parameters_physics.f90.& 
         !      The namelist <physics_parameters> does not exist. &
         !      Please replace this with the title <parameters_physics>")
      end if
      
      in_file = input_unit_exist("time_advance_knobs", old_nml_exist)
      if (old_nml_exist) then
         read(unit=in_file, nml=time_advance_knobs) 
         if (debug) then 
            write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            write(*,*) 'Please replace the namelist <time_advance_knobs> in the input file.'
           write(*,*) 'Refer to the input paramters text file as to which namelist to use.'
            write(*,*) 'Some of these parameters have been moved to <run_parameters>'
            write(*,*) 'and others have been moves to <physics_parameters>.'
            write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         end if

           ! write(*,*) "Aborting in run_parameters.f90.&
         !      The namelist <time_advance_knobs> does not exist.&
         !      Please replace this with the title <numerical>"
         ! call mp_abort("Aborting in run_parameters.f90.&
         !      The namelist <time_advance_knobs> does not exist.&
         !      Please replace this with the title <run_parameters>")
      end if

   end subroutine check_backwards_compatability
    
   !**********************************************************************
   !                         BROADCAST OPTIONS                           !
   !**********************************************************************
   ! Broadcast these parameters to all the processors - necessary because
   ! the above was only done for the first processor (proc0).
   !**********************************************************************
   subroutine broadcast_parameters

     use mp, only: broadcast

     implicit none 

     call broadcast(include_parallel_streaming)
     call broadcast(include_mirror)
     call broadcast(nonlinear)
     call broadcast(xdriftknob)
     call broadcast(ydriftknob)
     call broadcast(wstarknob)

     call broadcast(adiabatic_option_switch)
     call broadcast(zonal_init_option_switch)
     call broadcast(zonal_closure_option_switch)

     call broadcast(prp_shear_enabled)
     call broadcast(hammett_flow_shear) 
     call broadcast(include_pressure_variation)
     call broadcast(include_geometric_variation)
     call broadcast(include_parallel_nonlinearity)
     call broadcast(suppress_zonal_interaction)
     call broadcast(only_zonal_interaction)
     call broadcast(freeze_nonzonal)
     call broadcast(freeze_zonal)
     call broadcast(freeze_zonal_factor)
     call broadcast(freeze_zonal_kmin)
     call broadcast(freeze_zonal_kmax)
     !> Both of these must be broadcast: only proc0 reads the input file, and
     !> <use_analytic_drift_phase> in rosenbluth_hinton gates a reduction, so a
     !> rank that disagrees about them deadlocks the run rather than getting a
     !> wrong answer.
     call broadcast(RH_analytic_drift_phase)
     call broadcast(RH_analytic_drift_phase_specified)
     call broadcast(zonal_PS_fac)
     call broadcast(zonal_usym_fac)
     call broadcast(zonal_g_exb)
     call broadcast(zonal_nkx)
     call broadcast(zonal_rh_fac)
     
     call broadcast(full_flux_surface)
     call broadcast(include_apar)
     call broadcast(include_bpar)
     call broadcast(radial_variation)

     call broadcast(beta)
     call broadcast(vnew_ref)
     call broadcast(zeff)
     call broadcast(rhostar)
     call broadcast(tite)
     call broadcast(nine)
     call broadcast(g_exb)
     call broadcast(g_exbfac)
     call broadcast(omprimfac)
     call broadcast(omprimfac_RH)
     call broadcast(omprimfac_PS)

   end subroutine broadcast_parameters

 end subroutine read_parameters_physics

 !**********************************************************************
 !                      FINISH READ PARAMETERS                         !
 !**********************************************************************
 ! Set the initialised flag to be false such that we do not initialise
 ! twice.
 !> TODO-HT or TODO-GA: sed: initialised -> initialized
 !**********************************************************************
 subroutine finish_read_parameters_physics
   implicit none
   initialised = .false.
 end subroutine finish_read_parameters_physics

end module parameters_physics
