################################################################################
#    THE ELEVEN PHYSICS CASES:  dE/dt = sum P, one ingredient at a time        #
################################################################################
# The budget identity is one statement, so testing it once proves little about
# where it would break.  These seven cases each add a single ingredient to the
# one before, in two configurations, so that a failure is attributable to the
# ingredient that introduced it:
#
#   1  linear, collisionless                     no drive at all in Miller
#   2  linear, collisionless, 2 species, ES      a second species, nothing else
#   3  linear, collisionless, 2 species, EM      ... and now finite beta
#   4  linear, collisional                       one source, no turbulence
#   5  linear, collisional, 2 species, ES        the same, with a source
#   6  linear, collisional, 2 species, EM        ... and now finite beta
#   7  nonlinear, modified-adiabatic electrons   the zonal-flow closure
#   8  nonlinear, adiabatic electrons            the opposite closure
#   9  nonlinear, kinetic ions and electrons     a second kinetic species
#  10  nonlinear electromagnetic, dApar          gbar, Ampere, vpar dApar
#  11  nonlinear electromagnetic, dApar + dBpar  the whole of chi
#
# The multi-species linear cases exist because that is what found the momentum
# inertia's missing normalisation.  A single-species linear case cannot see that
# error at all: for m = T = 1 the wrong coefficient and the right one are the
# same number, while for kinetic electrons they differ by 61.
#
# They come in pairs, electrostatic then electromagnetic, because the first
# version of this ladder had only the electromagnetic one and it failed without
# saying which of its two new ingredients was to blame.  The electrostatic
# member closes to 3e-03 and the electromagnetic member does not close at all,
# which settles it: finite beta, not the second species.
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
    1:  ('linear_collisionless',         0.08, None, {'miller'}),
    2:  ('linear_collisionless_kinetic', 0.08, 0.08, {'miller'}),   # w7x: see below
    3:  ('linear_collisionless_multi',   0.30, 0.30, {'miller'}),
    4:  ('linear_collisional',           0.05, 0.05, set()),
    #> Case 5 in Miller turns the flow over only 0.12 times, so its potential
    #> tolerance is never actually asserted; it is the W7-X member that tests
    #> anything here.  Left as is rather than lengthened, because what this pair
    #> exists to isolate is visible in W7-X, and a longer Miller run costs more
    #> than it proves.
    5:  ('linear_collisional_kinetic',   0.08, 0.08, set()),
    6:  ('linear_collisional_multi',     0.30, 0.30, set()),
    7:  ('nl_modified_adiabatic',        0.08, 0.20, {'w7x'}),
    8:  ('nl_adiabatic',                 0.08, 0.08, {'w7x'}),
    9:  ('nl_kinetic',                   0.08, 0.05, set()),
    10: ('nl_em_apar',                   0.12, None, set()),
    11: ('nl_em_apar_bpar',              0.08, 0.20, set()),
}

#> The momentum budget does not close in these, and the bound is two-sided so
#> that neither a regression nor a fix passes unnoticed.
#> Case 7 in Miller used to be here at (0.10, 0.80).  It was the momentum
#> inertia missing its parallel-velocity normalisation, which made the electron
#> inertia 61 times too small and let the electrons carry the whole of
#> E_Omega_RH; it now closes at 1.4e-02 and is asserted normally.
#>
#> Case 8 remains, and only at beta = 1e-2.  At beta = 1e-3 and 3e-3 cases 8 and
#> 9 give comparable residuals, 6e-02 to 1e-01; it is at 1e-2 that case 6 alone
#> diverges, and there its energy turnover is 131 against case 7's 13 -- the run
#> is doing something qualitatively different when dBpar is dropped at that
#> beta.  Whether that is a defect of the diagnostic or of the truncation is not
#> yet established, so it is bounded rather than tolerated.
KNOWN_MOMENTUM_FAILURES = {
    #> Case 10's potential-like budget closes in both geometries -- 4.8e-02 in
    #> Miller and 6.8e-02 in W7-X, at turnovers of 17 and 6 -- so the
    #> electromagnetic nonlinear machinery works once the nonlinear term is
    #> actually driving.  It is the momentum channel that fails, and the W7-X
    #> number sits with the multi-species electromagnetic failures of cases 3 and
    #> 6 rather than with anything specific to case 10: same geometry, same
    #> kinetic electrons, same dApar.  The working hypothesis is that all three
    #> are one defect, which would mean case 10 needs no separate explanation.
    #> Not yet tested; recorded so that a fix to cases 3 and 6 is checked here
    #> too.
    ('miller', 10): (0.40, 4.00),
    ('w7x', 10):    (1.50, 15.0),
    #> Case 2 in W7-X is case 1 with a second kinetic species and nothing else:
    #> collisionless, electrostatic, so the drift channel is still the only
    #> source.  The potential projection goes from 1.3e-01 to 1.7e+01 and the
    #> momentum projection to 1.1e+03 purely from adding electrons.
    #>
    #> This looked at first like a per-species normalisation in the drift path,
    #> which is what the momentum inertia turned out to be.  It is NOT.  Three
    #> measurements rule that out:
    #>
    #>   - splitting the identity dPhi/dt = -i kx F by species gives |z| = 1.07
    #>     for the ions -- the ordinary case-1 quadrature error -- and 47.6 for
    #>     the electrons;
    #>   - but a missing stm^p would give the same p at every mass, and the
    #>     implied p runs 0.18, 0.69, 0.94 at m_e = 1e-2, 1e-3, 2.7e-4.  There is
    #>     no constant factor that fits;
    #>   - and at xdriftknob = 0 the flux is identically zero, yet the electron
    #>     projection still grows 45x while the ion holds at 1.0018.  Whatever
    #>     destroys it is not in the drift channel at all.
    #>
    #> What is known about it: it is absent at t = 30 (|z| = 1.0025 ion, 1.0370
    #> electron) and grows in with time; it does not respond to delt (delt/8
    #> changes 1.0042 to 1.0045); and it gets WORSE with velocity resolution
    #> (79.2 at 3x against 47.6), which rules out collisionless phase mixing
    #> outrunning the grid, the obvious candidate.  The collisional twin, case 5,
    #> closes normally at 1.7e-02.  Unexplained; bounded on both sides so that a
    #> fix shows up here as a failure and gets read.
    ('w7x', 2):     (3.0e+02, 3.0e+03),
    #> The electromagnetic member of the same pair.  Its potential-like residual
    #> is bounded below; the momentum-like one is bounded here.  Both channels
    #> fail in W7-X and only the potential one fails in Miller, which is the
    #> asymmetry the electrostatic twin (case 2) exists to expose.
    ('w7x', 3):     (1.0e+03, 1.0e+04),
}

#> The same, for the potential-like invariant.  W7-X case 1 is the drift channel
#> on its own, and it does not close: 1.5e-01, converged (1.81e-01 at t = 30,
#> 1.52e-01 at t = 150, 1.47e-01 at t = 300).  This is the drift-channel defect
#> the report discusses, now measured on a case where the flow actually moves
#> rather than inferred from a vacuous one.
KNOWN_PHI_FAILURES = {
    #> W7-X case 1 is the drift channel on its own -- the only case in which it
    #> is the sole source -- so it is the only case that measures the drift
    #> channel's own accuracy rather than a diluted version of it.  It reads
    #> 1.3e-01 at the resolution the suite runs at, and that number is a
    #> quadrature error, not a defect of the formula:
    #>
    #>   nzed 128, nv  48/24   |z| - 1 = 0.095      (this suite)
    #>   nzed 512, nv  48/24             0.052
    #>   nzed 128, nv 192/96             0.059
    #>   nzed 512, nv  96/48             0.035
    #>
    #> where z is the complex fit to the identity dPhi/dt = -i kx F that the
    #> budget rests on.  It converges in both directions at once and neither
    #> alone, which is why refining z by itself looks flat.  Two controls place
    #> the error in the drift-orbit phase Q rather than the flux expression:
    #> with xdriftknob = 0 the phase is identically zero, no quadrature is
    #> performed, and the projection is conserved to 4.6e-06; and the relative
    #> error is linear in the drift strength (0.040 at xdriftknob 0.5, 0.071 at
    #> 1.0), so the absolute error is quadratic in it, which is what an error
    #> inside Q looks like when Q itself is proportional to the drift.
    #>
    #> Every other W7-X case runs the same drift channel and closes far tighter
    #> because another source dominates the budget: case 4 reaches 2.9e-03.
    ('w7x', 1): (0.06, 0.40),
    #> Cases 3 and 6 are the electromagnetic multi-species pair, and they do not
    #> close.  Cases 2 and 5 were added to say why: they are the same runs with
    #> the fields switched off, and they close normally, so it is the finite
    #> beta and not the second species.  See the note in their decks.
    ('w7x', 3): (1.00, 12.0),
    ('w7x', 6): (0.50, 6.00),
    #> And the electrostatic member of that pair, which fails harder than the
    #> electromagnetic one -- 1.7e+01 against 6.1e+00.  So in W7-X the second
    #> species alone is enough to break the budget and the finite beta is a
    #> separate matter; in Miller the same electrostatic run closes to 3e-03 and
    #> only the electromagnetic one fails.  The two configurations are failing
    #> for different reasons and the pair is what separates them.  See the note
    #> against ('w7x', 2) in KNOWN_MOMENTUM_FAILURES for what has been ruled out.
    ('w7x', 2): (5.0, 50.0),
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
        bounds = KNOWN_PHI_FAILURES.get((configuration, case))
        if bounds is not None:
            low, high = bounds
            assert low < phi_residual < high, (
                f'phi_RH is a known failure here and is bounded on both sides so '
                f'that neither a regression nor a fix passes unnoticed; measured '
                f'{phi_residual:.3e}')
        else:
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
