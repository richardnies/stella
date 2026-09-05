################################################################################
#     ROSENBLUTH-HINTON ENERGY BUDGET IN STELLARATOR GEOMETRY                  #
################################################################################
# The companion of ../test_1_rh_budget/, run in five VMEC equilibria instead of
# a Miller tokamak.  The assertion is the same self-consistency statement,
#
#     d E_RH / dt  =  sum_kx P_RH ,
#
# so it needs no reference data.  What changes is which channels carry it.
#
# In an axisymmetric field the bounce-averaged radial drift vanishes, so
# P_RH_drift is zero and test_1 never exercises it.  Here it is a real source:
# over the linear decks it carries 1% of the budget in W7-X up to 28% in TJ-II.
# That channel is the whole point of the stellarator generalisation, and these
# are the tests that hold it.
#
# Configurations.  W7-X standard, a quasi-axisymmetric and a quasi-helically
# symmetric design, TJ-II, and ITER.  ITER is VMEC but axisymmetric, so it is
# the control: it runs the same code path with a drift channel that ought to
# stay negligible, and it would catch a generalisation that broke the tokamak
# limit it has to reduce to.  All are taken on the alpha0 = 0.7 field line,
# away from the symmetry at alpha0 = 0 where the bounce-averaged drift very
# nearly vanishes and both sides of the budget fall to the noise floor.
#
# Equilibria.  The VMEC files are not in the repository -- together they are
# some 24 MB -- so each test skips if its wout file is absent rather than
# failing.  Drop them beside the decks to enable the suite.
#
# Windows.  The nonlinear decks start from noise in an 8x8 box, so early on
# E_RH is round-off and late on the fields run away with no cascade to saturate
# into.  Each case is compared over the window where the zonal flow is
# genuinely nonlinearly driven and the run is still resolved.  These grow more
# slowly than the tokamak decks -- the ITG needs until about t = 25 to reach
# nonlinear amplitude -- so the windows sit later, around t = 28 to 44.  Past
# t = 45 the residual degrades exactly as the box runs away: for W7-X it goes
# 1.4e-3 over 30-40, 8.0e-2 over 45-55, and 3.9e-1 over 50-58.
#
# Channels.  As in test_1, the linear decks are checked on the total budget and
# the nonlinear ones on the nonlinear channel alone, subtracting the
# collisional and drift channels which the linear decks already verify.
#
# Tolerances.  Measured residuals are 2.0e-2 to 9.3e-2 on the linear decks and
# 1.4e-3 to 4.0e-2 on the nonlinear ones.  The tolerances below sit a factor of
# roughly two above those, which discriminates against a real break while
# absorbing the reordering of reductions that a different MPI decomposition
# causes.  The linear decks need the looser figure: they run at finite
# collisionality where RH_fluxes_collisional is evaluated as
# (g^{n+1} - g^n)/code_dt and so carries a first-order-in-delt error, the same
# residual floor documented in test_1.
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

#> Which VMEC equilibrium each configuration needs.  Keyed by the name that
#> appears in the deck filenames.
VMEC_FILE = {'ITER': 'wout_iter.nc',
             'W7X': 'wout_w7x_standard.nc',
             'QA': 'wout_QA.nc',
             'QH': 'wout_QH.nc',
             'TJII': 'wout_tjii.nc'}


#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")


def require_equilibrium(configuration):
    '''Skip rather than fail when the VMEC file for <configuration> is absent.

    The equilibria are too large to keep in the repository, so a checkout that
    does not have them should report these tests as skipped, not broken.
    '''
    wout = pathlib.Path(__file__).parent / VMEC_FILE[configuration]
    if not wout.exists():
        pytest.skip(f'VMEC equilibrium {VMEC_FILE[configuration]} is not present '
                    f'beside {pathlib.Path(__file__).parent.name}/; '
                    f'drop it there to run the {configuration} benchmarks.')


def check_stellarator_budget(configuration, input_filename, tmp_path, stella_version,
                             tolerance, time_min=None, time_max=None,
                             require_decay=False, channel='total', error=False):
    '''Run <input_filename> in <configuration> and assert the RH budget closes.

    Reports the share of the budget carried by the drift channel, which is what
    separates these cases from their tokamak counterparts.
    '''
    require_equilibrium(configuration)

    run_local_stella_simulation(input_filename, tmp_path, stella_version,
                                vmec_file=VMEC_FILE[configuration])
    local_netcdf_file = tmp_path / input_filename.replace('.in', '.out.nc')

    time, E_RH, dE_RH_dt, P_RH, P_nonlinear, P_collisional, P_drift, P_tr, P_pa = get_rh_budget(
        local_netcdf_file, time_min, time_max)

    if channel == 'nonlinear':
        measured, expected, what = dE_RH_dt - P_collisional - P_drift, P_nonlinear, 'nonlinear channel'
    else:
        measured, expected, what = dE_RH_dt, P_RH, 'total budget'
    residual = np.linalg.norm(measured - expected) / np.linalg.norm(expected)

    #> Same guard against a vacuous pass as in test_1: if the zonal flow barely
    #> moved then both sides are near zero and the budget closes without saying
    #> anything.
    integrand = np.abs(dE_RH_dt)
    energy_turnover = np.sum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(time)) / E_RH.mean()
    if not (energy_turnover > 0.5):
        print('\nERROR: The zonal flow barely evolved, so the budget test is vacuous.'); error = True
        print(f'    integral |dE_RH/dt| dt / mean(E_RH) = {energy_turnover:.4f}   (need > 0.5)')

    if require_decay and not (E_RH[-1] < 0.5 * E_RH[0]):
        print('\nERROR: The zonal flow did not decay.'); error = True
        print(f'    E_RH(start) = {E_RH[0]:14.6e}    E_RH(end) = {E_RH[-1]:14.6e}')

    if not (residual < tolerance):
        print(f'\nERROR: The RH energy budget does not close for {input_filename}.'); error = True
        print(f'    {what}, relative L2 residual = {residual:14.6e}   (tolerance {tolerance:.1e})')
        print(f'    nonlinear channel peaks at {np.abs(P_nonlinear).max():.6e}, '
              f'collisional at {np.abs(P_collisional).max():.6e}, '
              f'drift at {np.abs(P_drift).max():.6e} '
              f'(trapped {np.abs(P_tr).max():.3e}, passing {np.abs(P_pa).max():.3e})')
        print(f'    {"time":>10} {"measured":>16} {"expected":>16} {"ratio":>10}')
        for i in range(0, len(time), max(1, len(time) // 12)):
            ratio = measured[i] / expected[i] if expected[i] != 0 else np.nan
            print(f'    {time[i]:10.3f} {measured[i]:16.6e} {expected[i]:16.6e} {ratio:10.4f}')

    assert (not error), f'The RH energy budget does not close for {input_filename}.'

    drift_share = np.abs(P_drift).max() / max(np.abs(P_RH).max(), np.finfo(float).tiny)
    print(f'  -->  {input_filename}: {what} closes to {residual:.2e} (relative L2), '
          f'drift channel carries {100 * drift_share:.0f}% of the peak.')
    return residual


#-------------------------------------------------------------------------------
#            LINEAR ZONAL FLOW WITH COLLISIONS, IN EACH CONFIGURATION          #
#-------------------------------------------------------------------------------
# Checked on the total budget: with the run linear the nonlinear fluxes vanish
# and the two remaining sources are the collisional and drift channels, so this
# is the cleanest test of the drift channel that the stellarator work adds.
#-------------------------------------------------------------------------------
@pytest.mark.parametrize('configuration', ['ITER', 'W7X', 'QA', 'QH', 'TJII'])
def test_whether_rh_budget_closes_for_linear_collisional_zonal_flow(configuration, tmp_path, stella_version):
    '''A zonal flow relaxed by collisions and by the bounce-averaged magnetic
    drift.  ITER is the axisymmetric control, where the drift channel should
    stay small and the case reduces to the tokamak benchmark.'''
    check_stellarator_budget(configuration, f'{configuration}_linear_collisional.in',
                             tmp_path, stella_version,
                             tolerance=0.15, require_decay=True, channel='total')
    return


#-------------------------------------------------------------------------------
#           MINIMAL NONLINEAR RUNS, MODIFIED ADIABATIC ELECTRONS               #
#-------------------------------------------------------------------------------
#> Windows are per configuration because the ITG grows at a different rate in
#> each, so the interval over which the zonal flow is nonlinearly driven and the
#> box is still resolved does not sit at the same time in all of them.
NONLINEAR_WINDOW = {'ITER': (34.0, 44.0), 'W7X': (30.0, 40.0), 'QA': (28.0, 38.0),
                    'QH': (34.0, 44.0), 'TJII': (28.0, 38.0)}


@pytest.mark.parametrize('configuration', ['ITER', 'W7X', 'QA', 'QH', 'TJII'])
def test_whether_rh_budget_closes_for_nonlinear_modified_adiabatic_electrons(configuration, tmp_path, stella_version):
    '''Zonal flow driven nonlinearly by an ITG mode, with the flux-surface-average
    term retained in the adiabatic electron response.'''
    time_min, time_max = NONLINEAR_WINDOW[configuration]
    check_stellarator_budget(configuration, f'{configuration}_nl_adiabatic_electrons.in',
                             tmp_path, stella_version, tolerance=0.08,
                             time_min=time_min, time_max=time_max, channel='nonlinear')
    return


#-------------------------------------------------------------------------------
#          MINIMAL NONLINEAR RUNS, UNMODIFIED ADIABATIC ELECTRONS              #
#-------------------------------------------------------------------------------
NONLINEAR_WINDOW_IONS = {'ITER': (34.0, 44.0), 'W7X': (30.0, 40.0), 'QA': (28.0, 38.0),
                         'QH': (30.0, 40.0), 'TJII': (30.0, 40.0)}


@pytest.mark.parametrize('configuration', ['ITER', 'W7X', 'QA', 'QH', 'TJII'])
def test_whether_rh_budget_closes_for_nonlinear_unmodified_adiabatic_electrons(configuration, tmp_path, stella_version):
    '''As above, but with a plain Boltzmann electron response (no
    flux-surface-average term), which is the opposite adiabatic closure and
    weights the zonal part of the potential differently.'''
    time_min, time_max = NONLINEAR_WINDOW_IONS[configuration]
    check_stellarator_budget(configuration, f'{configuration}_nl_adiabatic_ions.in',
                             tmp_path, stella_version, tolerance=0.08,
                             time_min=time_min, time_max=time_max, channel='nonlinear')
    return
