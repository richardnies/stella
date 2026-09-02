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
# The budget is expected to close in tokamak geometry with no hyperdissipation
# and no sponge.  It does not close exactly: stella advances the zonal mode with
# a particular discretisation, while the RH projection that is supposed to
# annihilate the linear streaming and drift terms is a continuum construction.
# The leftover is a discretisation mismatch, not a missing term -- it is
# insensitive to upwinding and to velocity/parallel resolution, but it grows by
# an order of magnitude if the drifts are advanced explicitly rather than
# implicitly.
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
#                     LINEAR ZONAL FLOW WITH COLLISIONS                        #
#-------------------------------------------------------------------------------
def test_whether_rh_budget_closes_for_linear_collisional_zonal_flow(tmp_path, stella_version, error=False):
    '''The run is linear, so the nonlinear RH fluxes are identically zero and the
    only source in the budget is the collisional flux.  This is the cleanest
    test of the RH diagnostic and the tightest tolerance in this file.'''

    input_filename = 'rh_linear_collisional.in'
    run_local_stella_simulation(input_filename, tmp_path, stella_version)
    local_netcdf_file = tmp_path / input_filename.replace('.in', '.out.nc')

    time, E_RH, dE_RH_dt, P_RH = get_rh_budget(local_netcdf_file)
    residual = np.linalg.norm(dE_RH_dt - P_RH) / np.linalg.norm(P_RH)

    # The zonal flow must actually decay, otherwise the budget is trivially
    # satisfied by everything being zero and the test proves nothing.
    if not (E_RH[-1] < 0.5 * E_RH[0]):
        print('\nERROR: The zonal flow did not decay, so the budget test is vacuous.'); error = True
        print(f'    E_RH(start) = {E_RH[0]:14.6e}')
        print(f'    E_RH(end)   = {E_RH[-1]:14.6e}')

    if not (residual < 0.05):
        print('\nERROR: The Rosenbluth-Hinton energy budget does not close.'); error = True
        print(f'    relative L2 residual = {residual:14.6e}   (tolerance 5.0e-02)')
        print(f'    {"time":>10} {"dE_RH/dt":>16} {"sum P_RH":>16} {"ratio":>10}')
        for i in range(0, len(time), max(1, len(time) // 10)):
            ratio = dE_RH_dt[i] / P_RH[i] if P_RH[i] != 0 else np.nan
            print(f'    {time[i]:10.3f} {dE_RH_dt[i]:16.6e} {P_RH[i]:16.6e} {ratio:10.4f}')

    assert (not error), 'The Rosenbluth-Hinton energy budget does not close for the linear collisional case.'
    print(f'  -->  The RH energy budget closes to {residual:.2e} (relative L2) for the linear collisional case.')
    return


#-------------------------------------------------------------------------------
#                        MINIMAL NONLINEAR CASES                               #
#-------------------------------------------------------------------------------
# TODO: these three decks run and produce output, but they do not yet reach a
# state in which the zonal flow is nonlinearly driven: starting from noise at
# ny = nx = 4, the initial condition decays before the ITG mode grows, E_RH
# falls to the 1e-13 level and the budget is dominated by numerical noise
# (measured residuals 3.6, 120 and 217 respectively, with |dE/dt| exceeding
# |P_RH| by factors of 5-335).  They need a setup that actually sustains a
# zonal flow -- a longer run into saturation, or a seeded zonal mode -- before
# a tolerance can be set honestly.  Until then they are skipped rather than
# asserted against a tolerance chosen to make them pass.
NONLINEAR_CASES = [
    'rh_nl_adiabatic_electrons.in',
    'rh_nl_kinetic.in',
    'rh_nl_electromagnetic.in',
]

@pytest.mark.skip(reason='Deck does not yet sustain a nonlinearly driven zonal flow; '
                         'budget is dominated by numerical noise. See TODO above.')
@pytest.mark.parametrize('input_filename', NONLINEAR_CASES)
def test_whether_rh_budget_closes_for_nonlinear_zonal_flow(tmp_path, stella_version, input_filename, error=False):
    run_local_stella_simulation(input_filename, tmp_path, stella_version)
    local_netcdf_file = tmp_path / input_filename.replace('.in', '.out.nc')

    time, E_RH, dE_RH_dt, P_RH = get_rh_budget(local_netcdf_file)
    residual = np.linalg.norm(dE_RH_dt - P_RH) / np.linalg.norm(P_RH)

    if not (residual < 0.05):
        print(f'\nERROR: The RH energy budget does not close for {input_filename}.'); error = True
        print(f'    relative L2 residual = {residual:14.6e}   (tolerance 5.0e-02)')

    assert (not error), f'The Rosenbluth-Hinton energy budget does not close for {input_filename}.'
    print(f'  -->  The RH energy budget closes to {residual:.2e} for {input_filename}.')
    return
