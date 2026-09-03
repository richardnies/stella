################################################################################
#          ROSENBLUTH-HINTON ENERGY BUDGET:  dE_RH/dt  vs  sum P_RH            #
################################################################################
# The Rosenbluth-Hinton diagnostic writes the zonal-flow response RH_phi_I, the
# RH inertia, and the fluxes that drive them.  If the diagnostic is correct then
# the zonal-flow energy it implies must obey
#
#     d E_RH / dt  =  sum_kx P_RH
#
# with E_RH and P_RH built from those outputs (see ../rh_budget.py).  Both sides
# come from the same run, so this needs no reference data: it checks the
# diagnostic against the code's own time evolution, and it keeps working when
# the physics of a case legitimately changes.
#
# The budget is expected to close in tokamak geometry with hyperdissipation and
# the tertiary sponge switched off.
#
# Wavelength.  Every nonlinear deck now puts its zonal modes at kx rho_i <= 1,
# where the Rosenbluth-Hinton construction applies -- it targets kx rho_i << 1/q,
# about 0.71 at q = 1.4.  The adiabatic decks use jtwist = 5 so their single
# zonal mode sits at kx = 0.5 rather than 2.5; the two-species decks span
# kx = 0.5 upwards and are restricted with kx_max.
#
# Time windows.  The nonlinear cases start from noise, so before the ITG mode
# has grown E_RH sits at the level of numerical round-off and the budget there
# is meaningless.  They are also run in a 4x4 box, which has no cascade to
# saturate into, so past the end of the growth phase the fields run away and the
# run stops being resolved.  Each nonlinear case is therefore compared over the
# window in which the zonal flow is genuinely nonlinearly driven and the run is
# still well behaved.  Use plot_rh_budget.py to inspect a run and choose one.
#
# Channels.  P_RH splits into a nonlinear part (the phi/apar/bpar fluxes) and a
# collisional part.  The linear case is checked on the total, because there the
# collisional channel is the only source and the nonlinear one is identically
# zero.  The nonlinear cases are checked on the nonlinear channel alone, by
# subtracting the collisional channel from dE_RH/dt: the collisional channel is
# already verified in isolation by the linear benchmark, so re-testing it in a
# nonlinear run adds nothing, while the nonlinear channel is what those decks
# exist to exercise.  Isolating it also tightens the check -- the collisional
# channel is about a quarter of the nonlinear one here, and carrying it along
# imported its own error.
#
# Tolerances.  The nonlinear decks set rng_seed, so repeating one reproduces its
# residual exactly; without that the noise initial condition differs every run
# and, since these runs grow exponentially, the residual varied by an order of
# magnitude between invocations of the same deck.  What a fixed seed does not fix
# is a change of MPI decomposition, which reorders the reductions, so the
# tolerance still has to carry a factor of a few.  Observed values on the
# nonlinear channel run from 1.9e-3 to 3.2e-2, and 8% discriminates sharply
# against a real break, which gives O(1).
#
# Residual floor.  The budget does not close exactly.  Part of the mismatch is
# first order in delt by construction: RH_fluxes_collisional is evaluated as
# (g^{n+1} - g^n)/code_dt.  Refining delt on a collisions-only variant of case 1
# takes the residual from 1.3e-2 at delt=0.2 to 7.0e-3 at delt=0.0125, where it
# converges.  The remaining ~7e-3 is not understood: it is not upwinding (zeroing
# every upwind coefficient leaves it unchanged), not parallel or velocity
# resolution, and not a failure to annihilate the linear streaming and drift
# terms (switching both off leaves it at 7.1e-3).  Tolerances are set with that
# floor in mind.
################################################################################

# Python modules
import pytest
import os, sys
import pathlib
import numpy as np
import xarray as xr

# Package to run stella
module_path = str(pathlib.Path(__file__).parent.parent.parent / 'run_local_stella_simulation.py')
with open(module_path, 'r') as file: exec(file.read())

# Shared RH budget helper
budget_path = str(pathlib.Path(__file__).parent.parent / 'rh_budget.py')
with open(budget_path, 'r') as file: exec(file.read())

#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")


#-------------------------------------------------------------------------------
#                              SHARED CHECK                                    #
#-------------------------------------------------------------------------------
def check_rh_budget(input_filename, tmp_path, stella_version, tolerance,
                    time_min=None, time_max=None, require_decay=False,
                    channel='total', kx_max=None, error=False):
    '''Run <input_filename> and assert that the RH energy budget closes.

    channel='total'      compares dE_RH/dt with the whole of P_RH.  Right for the
                         linear case, where the collisional channel is the only
                         source.
    channel='nonlinear'  subtracts the collisional channel from dE_RH/dt and
                         compares the remainder with the nonlinear channel.  Right
                         for the nonlinear cases: the collisional channel is
                         already verified on its own by the linear benchmark, so
                         what is under test here is the nonlinear one.
    '''

    run_local_stella_simulation(input_filename, tmp_path, stella_version)
    local_netcdf_file = tmp_path / input_filename.replace('.in', '.out.nc')

    time, E_RH, dE_RH_dt, P_RH, P_nonlinear, P_collisional, P_drift = get_rh_budget(
        local_netcdf_file, time_min, time_max, kx_max)

    if channel == 'nonlinear':
        measured, expected, what = dE_RH_dt - P_collisional - P_drift, P_nonlinear, 'nonlinear channel'
    else:
        measured, expected, what = dE_RH_dt, P_RH, 'total budget'
    residual = np.linalg.norm(measured - expected) / np.linalg.norm(expected)

    #> Guard against a vacuous pass.  If the zonal flow barely moves over the
    #> window then both sides of the budget are near zero and it is satisfied
    #> without testing anything.  The meaningful statement is that the energy
    #> actually changed: integrate |dE_RH/dt| over the window and require it to
    #> be a decent fraction of the typical E_RH.  That works whether the flow is
    #> decaying (the linear case) or being driven (the nonlinear ones), and does
    #> not care how far E_RH happens to travel in a particular realisation --
    #> a fixed growth factor is not robust to that, since these runs are not
    #> reproducible between invocations.
    # Trapezoidal integral written out: np.trapz was removed in numpy 2.0 and
    # np.trapezoid does not exist in the numpy < 2 that requirements.txt pins.
    integrand = np.abs(dE_RH_dt)
    energy_turnover = np.sum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(time)) / E_RH.mean()
    if not (energy_turnover > 0.5):
        print('\nERROR: The zonal flow barely evolved, so the budget test is vacuous.'); error = True
        print(f'    integral |dE_RH/dt| dt / mean(E_RH) = {energy_turnover:.4f}   (need > 0.5)')
        print(f'    E_RH ranges over {E_RH.min():.6e} .. {E_RH.max():.6e}')

    # The linear case must additionally show the collisional decay it exists to test.
    if require_decay and not (E_RH[-1] < 0.5 * E_RH[0]):
        print('\nERROR: The zonal flow did not decay.'); error = True
        print(f'    E_RH(start) = {E_RH[0]:14.6e}')
        print(f'    E_RH(end)   = {E_RH[-1]:14.6e}')

    if not (residual < tolerance):
        print(f'\nERROR: The Rosenbluth-Hinton energy budget does not close for {input_filename}.'); error = True
        print(f'    {what}, relative L2 residual = {residual:14.6e}   (tolerance {tolerance:.1e})')
        print(f'    nonlinear channel peaks at {np.abs(P_nonlinear).max():.6e}, '
              f'collisional at {np.abs(P_collisional).max():.6e}, '
              f'drift at {np.abs(P_drift).max():.6e}')
        print(f'    {"time":>10} {"measured":>16} {"expected":>16} {"ratio":>10}')
        for i in range(0, len(time), max(1, len(time) // 12)):
            ratio = measured[i] / expected[i] if expected[i] != 0 else np.nan
            print(f'    {time[i]:10.3f} {measured[i]:16.6e} {expected[i]:16.6e} {ratio:10.4f}')

    assert (not error), f'The Rosenbluth-Hinton energy budget does not close for {input_filename}.'
    print(f'  -->  {input_filename}: {what} closes to {residual:.2e} (relative L2).')
    return residual


#-------------------------------------------------------------------------------
#                     LINEAR ZONAL FLOW WITH COLLISIONS                        #
#-------------------------------------------------------------------------------
def test_whether_rh_budget_closes_for_linear_collisional_zonal_flow(tmp_path, stella_version):
    '''The run is linear, so the nonlinear RH fluxes are identically zero and the
    only source in the budget is the collisional flux.  This is the cleanest
    test of the RH diagnostic and the tightest tolerance in this file.'''
    check_rh_budget('rh_linear_collisional.in', tmp_path, stella_version,
                    tolerance=0.05, require_decay=True, channel='total')
    return


#-------------------------------------------------------------------------------
#                MINIMAL NONLINEAR RUNS WITH ADIABATIC ELECTRONS               #
#-------------------------------------------------------------------------------
def test_whether_rh_budget_closes_for_nonlinear_modified_adiabatic_electrons(tmp_path, stella_version):
    '''Zonal flow driven nonlinearly by an ITG mode, with the flux-surface-average
    term retained in the adiabatic electron response.'''
    check_rh_budget('rh_nl_adiabatic_electrons.in', tmp_path, stella_version,
                    tolerance=0.08, time_min=20.0, time_max=28.0, channel='nonlinear')
    return


def test_whether_rh_budget_closes_for_nonlinear_unmodified_adiabatic_electrons(tmp_path, stella_version):
    '''As above, but with a plain Boltzmann electron response (no
    flux-surface-average term), which is the opposite adiabatic closure.'''
    check_rh_budget('rh_nl_adiabatic_ions.in', tmp_path, stella_version,
                    tolerance=0.08, time_min=20.0, time_max=28.0, channel='nonlinear')
    return


#-------------------------------------------------------------------------------
#                        NONLINEAR ELECTROMAGNETIC                             #
#-------------------------------------------------------------------------------
def test_whether_rh_budget_closes_for_nonlinear_kinetic_electrons(tmp_path, stella_version):
    '''Nonlinear, kinetic ions and kinetic electrons.

    Restricted to kx <= 1.1.  The Rosenbluth-Hinton construction targets
    kx rho_i << 1/q, and the budget degrades smoothly as kx grows past that, so
    a case with a broad kx spectrum is judged on the modes the theory addresses.
    '''
    check_rh_budget('rh_nl_kinetic.in', tmp_path, stella_version,
                    tolerance=0.08, time_min=10.0, time_max=20.0,
                    channel='nonlinear', kx_max=1.1)
    return


def test_whether_rh_budget_closes_for_nonlinear_electromagnetic(tmp_path, stella_version):
    '''Nonlinear electromagnetic, exercising the apar and bpar RH flux channels
    alongside the electrostatic one.'''
    #> The window matters more here than for the electrostatic decks.  Earlier
    #> than about t = 8 the zonal flow has barely moved -- the energy turnover
    #> over the window is only about 2 -- and the budget is then dominated by
    #> round-off rather than by the physics, giving a residual of order one.  By
    #> t = 10 the flow is genuinely driven, turnover is around 7, and it closes
    #> to about 2e-2.
    check_rh_budget('rh_nl_electromagnetic.in', tmp_path, stella_version,
                    tolerance=0.08, time_min=10.0, time_max=18.0,
                    channel='nonlinear', kx_max=1.1)
    return


#-------------------------------------------------------------------------------
#                    INDEPENDENCE OF THE PARALLEL DOMAIN LENGTH                 #
#-------------------------------------------------------------------------------
def test_whether_rh_diagnostics_are_independent_of_parallel_length(tmp_path, stella_version, error=False):
    '''The same physics over one poloidal turn and over three must give the same
    Rosenbluth-Hinton quantities.

    An axisymmetric equilibrium simply repeats, so extending the flux tube from
    z in [-pi, pi] to [-5pi, 5pi] changes nothing physical.  The diagnostic gets
    this right because dl_over_b is normalised to unit sum and therefore cancels
    in the transit-average ratio, and because both B and the drift phase Q are
    2pi-periodic in an axisymmetric field.  Nothing in the transit average
    assumes a single turn.

    This is the guard on that: it is what would break first if the parallel
    integration weight or its normalisation were changed, and it is the
    precondition for the stellarator work, where the domain is genuinely longer
    than one turn.

    The agreement is not exact at finite resolution, and cannot be.  For a
    single-turn tube the maximum of B sits at the very end of the domain, which
    is where the turning points of the barely trapped particles are; there
    stella has no grid points beyond the boundary, so the spline that locates a
    turning point works from a one-sided window, while for three turns the same
    physical point is interior and gets a symmetric one.  That difference falls
    away with resolution -- 1.0e-4, 2.2e-5, 1.2e-5 at nzed = 48, 96, 192 -- so
    the decks here run at nzed = 96, resolved enough for the tolerance below to
    be testing the physics rather than the boundary.
    '''

    budgets = {}
    inertias = {}
    for input_filename in ('rh_linear_collisional.in', 'rh_linear_collisional_nperiod3.in'):
        run_directory = tmp_path / input_filename.replace('.in', '')
        run_directory.mkdir()
        run_local_stella_simulation(input_filename, run_directory, stella_version)
        local_netcdf_file = run_directory / input_filename.replace('.in', '.out.nc')
        budgets[input_filename] = get_rh_budget(local_netcdf_file)
        inertias[input_filename] = field_line_averaged_rh_inertia(local_netcdf_file)

    one_turn, three_turns = 'rh_linear_collisional.in', 'rh_linear_collisional_nperiod3.in'

    # The RH inertia is time-independent and sets the scale of everything else
    inertia_difference = np.max(np.abs(inertias[one_turn] - inertias[three_turns])) \
                       / np.max(np.abs(inertias[one_turn]))
    if not (inertia_difference < 1e-4):
        print('\nERROR: The RH inertia depends on the length of the parallel domain.'); error = True
        print(f'    nperiod = 1: {inertias[one_turn]}')
        print(f'    nperiod = 3: {inertias[three_turns]}')
        print(f'    relative difference = {inertia_difference:.6e}   (tolerance 1.0e-04)')

    # And the whole E_RH trajectory must follow
    E_one, E_three = budgets[one_turn][1], budgets[three_turns][1]
    energy_difference = np.max(np.abs(E_one - E_three)) / np.max(np.abs(E_one))
    if not (energy_difference < 1e-4):
        print('\nERROR: E_RH(t) depends on the length of the parallel domain.'); error = True
        print(f'    relative difference = {energy_difference:.6e}   (tolerance 1.0e-04)')
        print(f'    {"time":>10} {"nperiod=1":>16} {"nperiod=3":>16}')
        time = budgets[one_turn][0]
        for i in range(0, len(time), max(1, len(time) // 10)):
            print(f'    {time[i]:10.3f} {E_one[i]:16.6e} {E_three[i]:16.6e}')

    assert (not error), 'The Rosenbluth-Hinton diagnostics depend on the parallel domain length.'
    print(f'  -->  RH diagnostics are independent of parallel length: inertia agrees to '
          f'{inertia_difference:.1e}, E_RH(t) to {energy_difference:.1e}.')
    return
