################################################################################
#     THE SEVEN PHYSICS CASES:  dE/dt = sum P, one ingredient at a time        #
################################################################################
# The budget identity is one statement, so testing it once proves little about
# where it would break.  These seven cases each add a single ingredient to the
# one before, in two configurations, so that a failure is attributable to the
# ingredient that introduced it:
#
#   1  linear, collisionless                     no drive at all in Miller
#   2  linear, collisional                       one source, no turbulence
#   3  nonlinear, modified-adiabatic electrons   the zonal-flow closure
#   4  nonlinear, adiabatic electrons            the opposite closure
#   5  nonlinear, kinetic ions and electrons     a second kinetic species
#   6  nonlinear electromagnetic, dApar          gbar, Ampere, vpar dApar
#   7  nonlinear electromagnetic, dApar + dBpar  the whole of chi
#
# Every case puts the zonal mode at kx rho = 0.2, where FLR corrections are
# small but not negligible, and drives turbulence at ky rho = 0.2.  Miller lands
# on 0.2 exactly (y0 = 5, jtwist = 5, shat = 0.796); W7-X's twist-and-shift
# boundary fixes dkx = 0.2244 geometrically.
#
# THE STATISTIC.  The residual is the MEDIAN pointwise relative error over the
# last 40% of the run, restricted to where |sum P| is within two decades of its
# peak.  A relative L2 norm is unusable on these runs: the electromagnetic cases
# grow through six decades before the CFL condition stops them, so an L2 norm is
# set by a handful of points.  On Miller case 7 it reads 0.99, of which a single
# time step contributes 86%, while the budget closes to 1.3e-02 at the median.
#
# TURNOVER.  A residual means nothing if the flow did not move.  Cases whose
# energy turnover is below 0.5 are asserted to be vacuous rather than passing,
# so that a run which starts driving the flow -- because the physics or the
# resolution changed -- fails here and gets looked at.
#
# INITIAL CONDITIONS.  The linear decks excite a density perturbation AND a
# parallel flow.  From a zonal potential alone the momentum projection starts
# near zero and dividing by it reported a residual of 40 where there was no
# defect; from a flow alone the same happens to the potential projection.
################################################################################

import pathlib
import sys

import numpy as np
import pytest
import xarray as xr

module_path = str(pathlib.Path(__file__).parent.parent.parent / 'run_local_stella_simulation.py')
with open(module_path, 'r') as file: exec(file.read())

budget_path = str(pathlib.Path(__file__).parent.parent / 'rh_budget.py')
with open(budget_path, 'r') as file: exec(file.read())


@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")


#> case -> (deck stem, phi tolerance, Omega tolerance or None if it is a known
#>          failure, set of configurations where the case is vacuous)
#>
#> A `None` momentum tolerance records a case that does NOT close and is not
#> expected to: case 5 in Miller and case 6 in both geometries.  They are
#> asserted to stay bad rather than silently ignored, so that a fix shows up
#> here as a failure and gets read.
CASES = {
    1: ('linear_collisionless', 0.30, None, {'miller', 'w7x'}),
    2: ('linear_collisional',   0.05, 0.05, set()),
    3: ('nl_modified_adiabatic', 0.08, 0.20, {'w7x'}),
    4: ('nl_adiabatic',          0.08, 0.08, {'w7x'}),
    5: ('nl_kinetic',            0.08, 0.08, set()),
    6: ('nl_em_apar',            0.12, None, set()),
    7: ('nl_em_apar_bpar',       0.08, 0.20, set()),
}

#> The momentum budget does not close in these, and the bound is two-sided so
#> that neither a regression nor a fix passes unnoticed.
KNOWN_MOMENTUM_FAILURES = {
    ('miller', 5): (0.10, 0.80),
    ('miller', 6): (0.40, 4.00),
    ('w7x', 6):    (1.50, 15.0),
}

CONFIGURATIONS = ('miller', 'w7x')
TURNOVER_FLOOR = 0.5


def _measure(netcdf_file, which, window=0.4):
    '''Median relative residual and energy turnover over the driven window.'''
    fn = get_rh_budget if which == 'phi' else get_rh_omega_budget
    t, E, dEdt, P = fn(netcdf_file)[:4]
    if len(t) > 4:
        keep = t >= t[0] + (1.0 - window) * (t[-1] - t[0])
        t, E, dEdt, P = t[keep], E[keep], dEdt[keep], P[keep]
    integrand = np.abs(dEdt)
    turnover = np.sum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(t)) / E.mean()
    magnitude = np.abs(P)
    if magnitude.max() <= 0.0:
        return float('nan'), turnover          # no drive at all
    big = magnitude > magnitude.max() / 100.0
    residual = float(np.median(np.abs(dEdt - P)[big] / magnitude[big]))
    return residual, turnover


VMEC_FILE = {'w7x': 'wout_w7x_standard.nc'}


def _run(configuration, case, tmp_path, stella_version):
    stem = CASES[case][0]
    deck = f'{configuration}_{case}_{stem}.in'
    vmec = VMEC_FILE.get(configuration)
    if vmec is not None and not (pathlib.Path(__file__).parent / vmec).exists():
        #> The VMEC equilibria are large and are not all in the repository, so a
        #> checkout without them reports these as skipped rather than broken.
        pytest.skip(f'VMEC equilibrium {vmec} is not present beside '
                    f'{pathlib.Path(__file__).parent.name}/')
    run_local_stella_simulation(deck, tmp_path, stella_version, vmec_file=vmec)
    return tmp_path / deck.replace('.in', '.out.nc')


@pytest.mark.parametrize('configuration', CONFIGURATIONS)
@pytest.mark.parametrize('case', sorted(CASES))
def test_whether_the_budget_closes_for_each_physics_case(configuration, case,
                                                         tmp_path, stella_version):
    stem, phi_tol, omega_tol, vacuous = CASES[case]
    netcdf_file = _run(configuration, case, tmp_path, stella_version)

    phi_residual, phi_turnover = _measure(netcdf_file, 'phi')
    omega_residual, omega_turnover = _measure(netcdf_file, 'omega')
    print(f'\n  -->  case {case} ({stem}), {configuration}: '
          f'phi {phi_residual:.2e} (turnover {phi_turnover:.2f}), '
          f'Omega {omega_residual:.2e} (turnover {omega_turnover:.2f})')

    if configuration in vacuous:
        #> Recorded as vacuous, not as passing.  If the flow starts moving here
        #> the case has become a real test and its tolerance needs setting.
        assert phi_turnover < TURNOVER_FLOOR or omega_turnover < TURNOVER_FLOOR, (
            f'case {case} in {configuration} was recorded as vacuous but now '
            f'drives the flow (turnovers {phi_turnover:.2f}, {omega_turnover:.2f}); '
            f'give it a tolerance instead of a free pass')
        return

    if phi_turnover >= TURNOVER_FLOOR:
        assert phi_residual < phi_tol, (
            f'phi_RH budget: {phi_residual:.3e} > {phi_tol} '
            f'(turnover {phi_turnover:.2f})')

    if omega_turnover < TURNOVER_FLOOR:
        return
    bounds = KNOWN_MOMENTUM_FAILURES.get((configuration, case))
    if bounds is not None:
        low, high = bounds
        assert low < omega_residual < high, (
            f'Omega_RH is a known failure here and is asserted between {low} and '
            f'{high} so that neither a regression nor a fix passes unnoticed; '
            f'measured {omega_residual:.3e}')
    elif omega_tol is not None:
        assert omega_residual < omega_tol, (
            f'Omega_RH budget: {omega_residual:.3e} > {omega_tol} '
            f'(turnover {omega_turnover:.2f})')
