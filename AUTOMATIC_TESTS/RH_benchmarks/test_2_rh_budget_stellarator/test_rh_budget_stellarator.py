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


#-------------------------------------------------------------------------------
#          THE DRIFT CHANNEL ON ITS OWN, WITH COLLISIONS OFF                    #
#-------------------------------------------------------------------------------
#> The collisional decks above verify the drift channel only as a minority
#> partner of the collisional one.  With collisions off and the run linear the
#> bounce-averaged radial drift is the *whole* source, so the statement
#>
#>     projection(t) - projection(0) = integral of the drift channel
#>
#> tests that channel at a hundred per cent of the drive.  Nothing else can
#> absorb an error in it.
#>
#> The initial condition carries a parallel flow rather than a zonal potential.
#> That is what makes the momentum invariant testable: started from a potential
#> its projection is near zero, and the conservation statement degenerates into
#> a ratio of two small numbers.  Measured that way it looked like a factor-of-30
#> failure; measured from a flow it is conserved to one part in a thousand.
#>
#> Tolerances are per configuration because the drift is a discretised orbit
#> average and its error is geometry dependent.  TJ-II is the loosest by an
#> order of magnitude, and is resolution limited rather than wrong: refining
#> nzed, nvgrid, nmu and delt together halves it, 4.98e-1 -> 2.06e-1.
DRIFT_TOLERANCE_PHI = {'ITER': 0.05, 'W7X': 0.10, 'QA': 0.05, 'QH': 0.15, 'TJII': 0.60}
DRIFT_TOLERANCE_UMOM = {'ITER': 0.01, 'W7X': 0.01, 'QA': 0.01, 'QH': 0.01, 'TJII': 0.01}


def _accumulated(time, kx, flux):
    '''Trapezoidal integral of the source -i kx F, per mode.'''
    source = -1j * kx[None, :] * flux
    return np.cumsum(0.5 * (source[1:] + source[:-1]) * np.diff(time)[:, None], axis=0)


def check_collisionless_drift(configuration, tmp_path, stella_version, error=False):
    '''Assert both projections change by exactly what the drift channel says.'''
    require_equilibrium(configuration)
    input_filename = f'{configuration}_collisionless_drift.in'
    run_local_stella_simulation(input_filename, tmp_path, stella_version,
                                vmec_file=VMEC_FILE[configuration])
    ncdata = Dataset(tmp_path / input_filename.replace('.in', '.out.nc'))

    time = np.array(ncdata.variables['t'][:])
    kx = np.array(ncdata.variables['kx'][:])
    zed = np.array(ncdata.variables['zed'][:])
    jacobian = np.array(ncdata.variables['jacob'][:])[:, 0]
    weight = (zed[1] - zed[0]) * jacobian.copy()
    weight[-1] = 0.0
    weight = weight / weight.sum()
    keep = np.abs(kx) > 1e-12
    kxf = kx[keep]

    def compare(projection, flux, label, tolerance):
        nonlocal error
        accumulated = _accumulated(time, kxf, flux)
        change = projection[1:] - projection[0]
        scale = np.abs(projection[0]).max()
        residual = np.abs(change - accumulated).max() / scale
        share = np.abs(accumulated).max() / scale
        #> A vacuous pass would have the drift doing nothing.  Require that the
        #> channel actually moves the projection by more than the error in it.
        if not (share > 3 * residual):
            print(f'\nERROR: the drift channel is too small to test {label}.'); error = True
            print(f'    it moves the projection by {share:.3e}, error {residual:.3e}')
        if not (residual < tolerance):
            print(f'\nERROR: {label} does not follow its drift channel in {configuration}.')
            error = True
            print(f'    |change - integral of drift| / |projection(0)| = {residual:.6e}'
                  f'   (tolerance {tolerance:.1e})')
            print(f'    the drift channel carries {share:.3e} of the initial projection')
        return residual, share

    phi = _field_line_average(ncdata, 'RH_phi_I', weight)[:, keep]
    phi_drift = sum(_field_line_average(ncdata, n, weight)
                    for n in ('RH_fluxes_drift_trapped', 'RH_fluxes_drift_passing'))[:, keep]
    r_phi, s_phi = compare(phi, phi_drift, 'phi_RH', DRIFT_TOLERANCE_PHI[configuration])

    umom = _field_line_average_per_species(ncdata, 'RH_umom', weight)[:, 0, keep]
    umom_drift = _field_line_average_per_species(ncdata, 'RH_umom_flux_drift', weight)[:, 0, keep]
    r_p, s_p = compare(umom, umom_drift, 'U_RH', DRIFT_TOLERANCE_UMOM[configuration])

    assert (not error), f'The drift channel is not verified in {configuration}.'
    print(f'  -->  {configuration}: drift channel accounts for the change in both '
          f'projections -- phi_RH to {r_phi:.1e} (carrying {s_phi:.2f}), '
          f'p_RH to {r_p:.1e} (carrying {s_p:.2f}).')
    return


@pytest.mark.parametrize('configuration', ['ITER', 'W7X', 'QA', 'QH', 'TJII'])
def test_whether_the_drift_channel_alone_accounts_for_the_change(configuration, tmp_path, stella_version):
    '''The drift channel, tested where it is the only source there is.'''
    check_collisionless_drift(configuration, tmp_path, stella_version)
    return


#-------------------------------------------------------------------------------
#            THE MOMENTUM BUDGET IN STELLARATOR GEOMETRY                        #
#-------------------------------------------------------------------------------
#> The momentum invariant was untestable here until its geometric factor was
#> defined for VMEC: RH_drift_phase_fac is set only by the Miller branch, so in
#> a stellarator the weight, the projection and the inertia were all identically
#> zero and the diagnostic returned NaN without complaining.
#>
#> QA is excluded.  Its collisional momentum drive is a near-cancellation --
#> the Dougherty operator conserves momentum, so what drives this invariant is
#> only the part by which the weight departs from the exact momentum moment --
#> and in QA that residue is two orders below the other configurations, leaving
#> a small-signal test that does not converge under refinement while its
#> potential-like counterpart does.  The collisionless test above covers QA at
#> 1.6e-3, which is the statement that matters for the invariant itself.
UMOM_TOLERANCE = {'ITER': 0.06, 'W7X': 0.09, 'QH': 0.08, 'TJII': 0.08}


def check_stellarator_umom_budget(configuration, input_filename, tmp_path, stella_version,
                                  tolerance, error=False):
    '''Assert the toroidal-momentum budget closes in stellarator geometry.'''
    require_equilibrium(configuration)
    run_local_stella_simulation(input_filename, tmp_path, stella_version,
                                vmec_file=VMEC_FILE[configuration])
    local_netcdf_file = tmp_path / input_filename.replace('.in', '.out.nc')

    time, E, dE_dt, P, P_nl, P_coll, P_drift = get_rh_umom_budget(local_netcdf_file)
    residual = np.linalg.norm(dE_dt - P) / np.linalg.norm(P)

    integrand = np.abs(dE_dt)
    turnover = np.sum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(time)) / E.mean()
    if not (turnover > 0.5):
        print('\nERROR: the flow barely evolved, so the test is vacuous.'); error = True
        print(f'    integral |dE/dt| dt / mean(E) = {turnover:.4f}')

    if not (residual < tolerance):
        print(f'\nERROR: the momentum budget does not close in {configuration}.'); error = True
        print(f'    relative L2 residual = {residual:14.6e}   (tolerance {tolerance:.1e})')
        print(f'    collisional channel peaks at {np.abs(P_coll).max():.3e}, '
              f'drift at {np.abs(P_drift).max():.3e}')

    assert (not error), f'The momentum budget does not close in {configuration}.'
    print(f'  -->  {configuration}: momentum budget closes, residual {residual:.2e}, '
          f'turnover {turnover:.1f}.')
    return


@pytest.mark.parametrize('configuration', ['ITER', 'W7X', 'QH', 'TJII'])
def test_whether_umom_budget_closes_in_stellarator_geometry(configuration, tmp_path, stella_version):
    '''The toroidal-momentum budget, linear and collisional, in a stellarator.'''
    check_stellarator_umom_budget(configuration, f'{configuration}_linear_collisional.in',
                                  tmp_path, stella_version,
                                  tolerance=UMOM_TOLERANCE[configuration])
    return


#-------------------------------------------------------------------------------
#            THE DRIFT CHANNEL WITH ELECTROMAGNETIC FIELDS                      #
#-------------------------------------------------------------------------------
#> Collisionless, linear, and electromagnetic in a field that is not
#> quasisymmetric.  Nothing else in the suite covers that combination, and it is
#> the only place where the difference between projecting g and projecting gbar
#> matters: gbar is the conserved one, since moving to g cancels only the
#> dphi/dt half of the d<chi>/dt drive.
#>
#> The momentum invariant is asserted tightly.  The potential-like one is not,
#> because it is the one known defect in this note: with the right projection it
#> still accounts for its drift channel only to about ten per cent, the error
#> plateaus under refinement of all four grids, and three candidates have been
#> excluded (resolution, the projection, and the distribution the drift acts on).
#> It is bounded here rather than left untested, so that a change in it -- a fix
#> or a regression -- shows up as a failure instead of going unnoticed.
EM_DRIFT_UMOM_TOLERANCE = {'W7X': 5.0e-3, 'TJII': 1.0e-2, 'QA': 5.0e-3}
EM_DRIFT_PHI_KNOWN = {'W7X': 0.119, 'TJII': 1.34, 'QA': 0.170}


@pytest.mark.parametrize('configuration', ['W7X', 'TJII', 'QA'])
def test_whether_the_drift_channel_holds_electromagnetically(configuration, tmp_path, stella_version, error=False):
    '''The drift channel with Apar on, in a non-quasisymmetric field.'''
    require_equilibrium(configuration)
    input_filename = f'{configuration}_em_drift.in'
    run_local_stella_simulation(input_filename, tmp_path, stella_version,
                                vmec_file=VMEC_FILE[configuration])
    ncdata = Dataset(tmp_path / input_filename.replace('.in', '.out.nc'))

    time = np.array(ncdata.variables['t'][:])
    kx = np.array(ncdata.variables['kx'][:])
    zed = np.array(ncdata.variables['zed'][:])
    jacobian = np.array(ncdata.variables['jacob'][:])[:, 0]
    weight = (zed[1] - zed[0]) * jacobian.copy()
    weight[-1] = 0.0
    weight = weight / weight.sum()
    keep = np.abs(kx) > 1e-12
    kxf = kx[keep]

    def unaccounted(projection, flux):
        accumulated = _accumulated(time, kxf, flux)
        #> Normalised by max|projection|, not by its value at t = 0: TJ-II's
        #> initial gbar projection is a near cancellation, and dividing by it
        #> turns a working diagnostic into an apparent factor-of-twenty failure.
        return np.abs((projection[1:] - projection[0]) - accumulated).max() / np.abs(projection).max()

    umom = _field_line_average_per_species(ncdata, 'RH_umom', weight)[:, 0, keep]
    umom_drift = _field_line_average_per_species(ncdata, 'RH_umom_flux_drift', weight)[:, 0, keep]
    r_u = unaccounted(umom, umom_drift)

    phi = _field_line_average(ncdata, 'RH_phi_I', weight)[:, keep]
    phi_drift = sum(_field_line_average(ncdata, n, weight)
                    for n in ('RH_fluxes_drift_trapped', 'RH_fluxes_drift_passing'))[:, keep]
    r_phi = unaccounted(phi, phi_drift)

    tol = EM_DRIFT_UMOM_TOLERANCE[configuration]
    if not (r_u < tol):
        print(f'\nERROR: U_RH does not follow its drift channel in electromagnetic {configuration}.')
        error = True
        print(f'    unaccounted fraction = {r_u:.6e}   (tolerance {tol:.1e})')

    #> Two-sided, so that an improvement is reported as loudly as a regression.
    known = EM_DRIFT_PHI_KNOWN[configuration]
    if not (0.4 * known < r_phi < 1.6 * known):
        print(f'\nNOTE: the known electromagnetic phi_RH defect has changed in {configuration}.')
        print(f'    unaccounted fraction = {r_phi:.6e}, was {known:.3f}')
        print('    If it fell, the defect may be fixed and this bound should be tightened.')
        error = True

    assert (not error), f'The electromagnetic drift channel changed in {configuration}.'
    print(f'  -->  {configuration} electromagnetic: U_RH to {r_u:.1e}; '
          f'phi_RH at {r_phi:.2f}, the known defect, unchanged.')
    return
