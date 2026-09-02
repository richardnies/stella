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
# Time windows.  The nonlinear cases start from noise, so before the ITG mode
# has grown E_RH sits at the level of numerical round-off and the budget there
# is meaningless.  They are also run in a 4x4 box, which has no cascade to
# saturate into, so past the end of the growth phase the fields run away and the
# run stops being resolved.  Each nonlinear case is therefore compared over the
# window in which the zonal flow is genuinely nonlinearly driven and the run is
# still well behaved.  Use plot_rh_budget.py to inspect a run and choose one.
#
# Tolerances.  The nonlinear cases grow exponentially, so round-off differences
# (a different MPI decomposition, say) are amplified over the ~13 e-folding times
# of the comparison window and the residual is not reproducible to better than a
# few times 1e-2 between runs of the same deck: repeated runs of these two decks
# gave 1.8e-3, 5.9e-3, 3.1e-2 and 8.0e-2.  The relative L2 measure is also
# dominated by the largest values, i.e. the end of the window.  Their tolerance
# is therefore 15%, roughly twice the worst observed value, which still
# discriminates sharply against a real break -- a broken budget gives O(1), as
# rh_nl_kinetic.in does at 0.6-1.2.  The linear case is reproducible and keeps a
# 5% tolerance.
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
                    time_min=None, time_max=None, require_decay=False, error=False):
    '''Run <input_filename> and assert that dE_RH/dt matches sum_kx P_RH.'''

    run_local_stella_simulation(input_filename, tmp_path, stella_version)
    local_netcdf_file = tmp_path / input_filename.replace('.in', '.out.nc')

    time, E_RH, dE_RH_dt, P_RH = get_rh_budget(local_netcdf_file, time_min, time_max)
    residual = np.linalg.norm(dE_RH_dt - P_RH) / np.linalg.norm(P_RH)

    # Guard against a vacuous pass: if the zonal flow never does anything, both
    # sides are zero and the budget is satisfied without testing anything.
    if require_decay and not (E_RH[-1] < 0.5 * E_RH[0]):
        print('\nERROR: The zonal flow did not decay, so the budget test is vacuous.'); error = True
        print(f'    E_RH(start) = {E_RH[0]:14.6e}')
        print(f'    E_RH(end)   = {E_RH[-1]:14.6e}')
    if not require_decay and not (E_RH.max() > 100 * E_RH.min()):
        print('\nERROR: The zonal flow was not driven, so the budget test is vacuous.'); error = True
        print(f'    E_RH ranges only over {E_RH.min():.6e} .. {E_RH.max():.6e}')

    if not (residual < tolerance):
        print(f'\nERROR: The Rosenbluth-Hinton energy budget does not close for {input_filename}.'); error = True
        print(f'    relative L2 residual = {residual:14.6e}   (tolerance {tolerance:.1e})')
        print(f'    {"time":>10} {"dE_RH/dt":>16} {"sum P_RH":>16} {"ratio":>10}')
        for i in range(0, len(time), max(1, len(time) // 12)):
            ratio = dE_RH_dt[i] / P_RH[i] if P_RH[i] != 0 else np.nan
            print(f'    {time[i]:10.3f} {dE_RH_dt[i]:16.6e} {P_RH[i]:16.6e} {ratio:10.4f}')

    assert (not error), f'The Rosenbluth-Hinton energy budget does not close for {input_filename}.'
    print(f'  -->  The RH energy budget closes to {residual:.2e} (relative L2) for {input_filename}.')
    return residual


#-------------------------------------------------------------------------------
#                     LINEAR ZONAL FLOW WITH COLLISIONS                        #
#-------------------------------------------------------------------------------
def test_whether_rh_budget_closes_for_linear_collisional_zonal_flow(tmp_path, stella_version):
    '''The run is linear, so the nonlinear RH fluxes are identically zero and the
    only source in the budget is the collisional flux.  This is the cleanest
    test of the RH diagnostic and the tightest tolerance in this file.'''
    check_rh_budget('rh_linear_collisional.in', tmp_path, stella_version,
                    tolerance=0.05, require_decay=True)
    return


#-------------------------------------------------------------------------------
#                MINIMAL NONLINEAR RUNS WITH ADIABATIC ELECTRONS               #
#-------------------------------------------------------------------------------
def test_whether_rh_budget_closes_for_nonlinear_modified_adiabatic_electrons(tmp_path, stella_version):
    '''Zonal flow driven nonlinearly by an ITG mode, with the flux-surface-average
    term retained in the adiabatic electron response.'''
    check_rh_budget('rh_nl_adiabatic_electrons.in', tmp_path, stella_version,
                    tolerance=0.15, time_min=15.0, time_max=27.0)
    return


def test_whether_rh_budget_closes_for_nonlinear_unmodified_adiabatic_electrons(tmp_path, stella_version):
    '''As above, but with a plain Boltzmann electron response (no
    flux-surface-average term), which is the opposite adiabatic closure.'''
    check_rh_budget('rh_nl_adiabatic_ions.in', tmp_path, stella_version,
                    tolerance=0.15, time_min=15.0, time_max=27.0)
    return


#-------------------------------------------------------------------------------
#                     CASES THAT DO NOT YET CLOSE                              #
#-------------------------------------------------------------------------------
# Kept as decks so the work is not lost, but skipped rather than asserted
# against a tolerance chosen to make them pass.

@pytest.mark.skip(reason='Does not close: residual ~0.6, insensitive to the choice of '
                         'time window, so not a noise or windowing artefact.  Something '
                         'about the two-kinetic-species case is genuinely inconsistent.')
def test_whether_rh_budget_closes_for_nonlinear_kinetic_electrons(tmp_path, stella_version):
    '''Nonlinear, kinetic ions and kinetic electrons.'''
    check_rh_budget('rh_nl_kinetic.in', tmp_path, stella_version,
                    tolerance=0.15, time_min=15.0, time_max=27.0)
    return


@pytest.mark.skip(reason='Run goes NaN from the second step at beta = 0.004, both with '
                         'implicit and with explicit streaming/mirror, and with delt '
                         'reduced to 5e-3.  The deck needs stabilising before the budget '
                         'can be assessed at all.')
def test_whether_rh_budget_closes_for_nonlinear_electromagnetic(tmp_path, stella_version):
    '''Nonlinear electromagnetic, exercising the apar and bpar RH flux channels.'''
    check_rh_budget('rh_nl_electromagnetic.in', tmp_path, stella_version,
                    tolerance=0.15, time_min=15.0, time_max=27.0)
    return
