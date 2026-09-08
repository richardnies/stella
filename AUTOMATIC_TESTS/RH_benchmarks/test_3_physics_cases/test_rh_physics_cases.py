################################################################################
#    THE ELEVEN PHYSICS CASES:  dE/dt = sum P, one ingredient at a time        #
################################################################################
# The budget identity is one statement, so testing it once proves little about
# where it would break.  These eleven cases each add a single ingredient to the
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
# saying which of its two new ingredients was to blame.  In Miller the
# electrostatic member closes to 3e-03 and the electromagnetic member does not
# close at all, which settles it there: finite beta, not the second species.
# In W7-X the answer is less clean, because the electrostatic member found a
# third ingredient -- the flux tube -- that neither pair was designed to
# isolate; see the notes against cases 2 and 3 below.
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
# set by a handful of points.  On Miller case 11 it reads 0.99, of which a single
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
#> expected to: cases 1 and 10.  They are asserted to stay bad rather than
#> silently ignored, through the two-sided bounds in KNOWN_MOMENTUM_FAILURES,
#> so that a fix shows up here as a failure and gets read.
CASES = {
    1:  ('linear_collisionless',         0.08, None, {'miller'}),
    #> Case 2's potential tolerance is 0.12 rather than 0.08 because its W7-X
    #> deck is forced onto alpha0 = 0 (the zonal mode with kinetic electrons is
    #> numerically unstable on the alpha0 = 0.7 tube; see the deck), and on that
    #> line the parallel discretisation error case 1 shows is 0.13 for case 1
    #> and 9.0e-02 here.  Its momentum tolerance is never asserted: Miller by
    #> design, and W7-X because on alpha0 = 0 nothing drives the odd invariant.
    2:  ('linear_collisionless_kinetic', 0.12, 0.08, {'miller'}),
    3:  ('linear_collisionless_multi',   0.30, 0.30, {'miller'}),
    4:  ('linear_collisional',           0.05, 0.05, set()),
    #> Case 5 in Miller turns the flow over only 0.12 times, so its potential
    #> tolerance is never actually asserted; it is the W7-X member that tests
    #> anything here.  Left as is rather than lengthened, because what this pair
    #> exists to isolate is visible in W7-X, and a longer Miller run costs more
    #> than it proves.  Recorded as vacuous so that it cannot pass by default.
    5:  ('linear_collisional_kinetic',   0.08, 0.08, {'miller'}),
    #> Case 6 in Miller is the same: turnovers of 0.39 and 0.22, so neither
    #> channel is asserted there and it is the W7-X potential bound that holds
    #> the electromagnetic member.  Recorded as vacuous for the same reason.
    6:  ('linear_collisional_multi',     0.30, 0.30, {'miller'}),
    #> Cases 7 and 8 in W7-X start their noise at phiinit = 1.0 and run to
    #> t = 100, because from 0.01 the ITG mode on that tube has not saturated by
    #> the Miller decks' t = 50 and the window measured only the collisional
    #> decay of the noise (turnovers 0.05 to 0.47).  See the two decks.
    7:  ('nl_modified_adiabatic',        0.08, 0.20, set()),
    8:  ('nl_adiabatic',                 0.08, 0.08, set()),
    9:  ('nl_kinetic',                   0.08, 0.05, set()),
    10: ('nl_em_apar',                   0.12, None, set()),
    11: ('nl_em_apar_bpar',              0.08, 0.20, set()),
}

#> The momentum budget does not close in these, and the bound is two-sided so
#> that neither a regression nor a fix passes unnoticed.
#> Case 9 in Miller used to be here at (0.10, 0.80).  It was the momentum
#> inertia missing its parallel-velocity normalisation, which made the electron
#> inertia 61 times too small and let the electrons carry the whole of
#> E_Omega_RH; it now closes at 1.4e-02 and is asserted normally.
#>
#> Case 10 remains, and only at beta = 1e-2.  At beta = 1e-3 and 3e-3 cases 10
#> and 11 give comparable residuals, 6e-02 to 1e-01; it is at 1e-2 that case 10
#> alone diverges, and there its energy turnover is 131 against case 11's 13 --
#> the run is doing something qualitatively different when dBpar is dropped at
#> that beta.  Whether that is a defect of the diagnostic or of the truncation is not
#> yet established, so it is bounded rather than tolerated.
KNOWN_MOMENTUM_FAILURES = {
    #> W7-X case 1's momentum budget does not close, and it is the flux tube
    #> rather than the diagnostic.  The tube along alpha0 = 0.7 does not close
    #> on itself -- B jumps by 1% across the periodic join and kperp2 by a
    #> factor of four -- and the zonal mode is periodic in z regardless.  The
    #> odd invariant is carried entirely by passing particles, which cross the
    #> join once per transit, and E_Omega is not conserved by that dynamics: it
    #> swings between 0.36 and 1.0 of its initial value with a period of 39,
    #> undamped to t = 300, while the drift source integrates to a smooth
    #> -0.29.  The swing does not move with nzed or nvgrid, and it vanishes at
    #> alpha0 = 0, where the tube closes exactly and E_Omega is conserved to
    #> 2%.  At nfield_periods = 8 it was hidden inside a 3% drift of E_Omega
    #> over the run, which is what made the channel read as vacuous.  See the
    #> deck.  Measured 2.9 at t = 100 with drive 0.31; bounded, so that a
    #> change to the zonal boundary treatment shows up here.
    ('w7x', 1):     (1.0, 8.0),
    #> Case 10's potential-like budget closes in both geometries -- 4.8e-02 in
    #> Miller and 6.8e-02 in W7-X, at turnovers of 17 and 6 -- so the
    #> electromagnetic nonlinear machinery works once the nonlinear term is
    #> actually driving.  It is the momentum channel that fails, and the W7-X
    #> number sits with the multi-species electromagnetic failures of cases 3 and
    #> 6 rather than with anything specific to case 10: same geometry, same
    #> kinetic electrons, same dApar.  The working hypothesis is that all three
    #> are one defect, which would mean case 10 needs no separate explanation.
    #> Case 2 has since shown that the tube and the electrons are enough on
    #> their own: the zonal mode there is numerically unstable without any
    #> dApar (growth rate 0.025 at nfield_periods = 8; see its deck).  Whether
    #> that is what case 10 sees under its turbulence is not yet tested;
    #> recorded so that a fix to cases 3 and 6 is checked here too.
    ('miller', 10): (0.40, 4.00),
    ('w7x', 10):    (1.50, 15.0),
    #> Case 3 is the electromagnetic member of the pair whose electrostatic
    #> member is case 2.  Its potential-like residual is bounded below; the
    #> momentum-like one is bounded here.  Both channels fail in W7-X and only
    #> the potential one fails in Miller.  What case 2 turned out to show is
    #> that the alpha0 = 0.7 tube with kinetic electrons and no collisions is
    #> numerically unstable on its own, electrostatic or not (see its deck), so
    #> case 3 on that tube is measuring that growth as much as anything
    #> electromagnetic; moved to alpha0 = 0 it no longer grows but still fails,
    #> on its dApar transient (2.1e+01 at nfield_periods = 8, 8.0e-01 at 1), so
    #> the deck and the bound stay as they are until that is understood.
    ('w7x', 3):     (1.0e+03, 1.0e+04),
}

#> The same, for the potential-like invariant.
#>
#> W7-X case 1 used to be here at (0.06, 0.40): it is the drift channel on its
#> own -- the only case in which that is the sole source -- and at
#> nfield_periods = 8 it read 1.2e-01 to 1.5e-01, converged in time and flat
#> under nzed alone.  It is now asserted normally, at 2.5e-02, after the deck
#> was shortened to one field period.  The residual is a discretisation error
#> of the parallel dynamics, and it tracks how rugged the field line is per
#> grid point: at nzed = 128, nfp 1/2/4/8 give max |dB| per step of
#> 0.005/0.010/0.033/0.049 and residuals 0.025/0.058/0.086/0.150.  At nfp = 2
#> it converges with resolution, slowly and in nzed and nvgrid jointly (nzed
#> 64/128/256/512: 0.083/0.058/0.042/0.038; nvgrid 48 -> 96: 0.058 -> 0.045),
#> and it is independent of delt, of the mirror scheme and of implicit against
#> explicit drifts.  Restarting at t = 30 with one operator at a time locates
#> it: the drift projection agrees with the continuum to five digits, and the
#> whole of the error is the discrete streaming + mirror projection missing the
#> exact cancellation that the identity dPhi/dt = -i kx F rests on, by 10% of
#> |dPhi/dt| at nfp = 8.  Two earlier readings are superseded by this.  The
#> nfp = 5..12 scan's "resid ~ turnover^-1.44" was the drive varying with tube
#> length on top of the ruggedness, and the turnover alone accounts for 2x of
#> the 6x between nfp = 8 and nfp = 1.  And alpha0 = 0 is not where the drift
#> is significant: it is the stellarator-symmetric line, on which the
#> transit-averaged drift of every passing particle vanishes, so the drive is
#> trapped-only and the residual is worse (0.13 at nfp = 1, 0.25 at nfp = 2,
#> 0.39 at nfp = 8).  See the deck for the numbers behind each of these.
KNOWN_PHI_FAILURES = {
    #> Cases 3 and 6 are the electromagnetic multi-species pair, and they do not
    #> close.  Cases 2 and 5 are the same runs with the fields switched off.
    #> Case 5 closes on the same tube, 1.7e-02, so with collisions it is the
    #> finite beta.  Case 2 does not, but for a reason of its own: with kinetic
    #> electrons and no collisions the zonal mode on the alpha0 = 0.7 tube is
    #> numerically unstable, electrostatic or not (growth rate 0.025 at
    #> nfield_periods = 8, 0.145 at 1), and its old reading here, 1.7e+01, was
    #> that growth.  Its deck now runs on alpha0 = 0, where the mode is stable,
    #> and it is asserted normally at 9.0e-02.  Case 3 moved to alpha0 = 0 no
    #> longer grows but still fails on its dApar transient, so the finite beta
    #> is at least part of what these two bounds hold; how much of the rest is
    #> the instability is open.  See the notes in the decks of cases 2 and 3.
    ('w7x', 3): (1.00, 12.0),
    ('w7x', 6): (0.50, 6.00),
}

#> Only Miller and W7-X are asserted.  The other four equilibria have been run
#> on all eleven cases -- the decks are here and the measurements are below --
#> but no tolerances are set for them, because setting forty-four bounds from
#> one measurement each is bookkeeping dressed up as verification.  Add them
#> deliberately, case by case, as each is understood.
#>
#> The sweep was run to test a prediction, and refuted it.  The prediction was
#> that case 1, which is the drift channel on its own, would be worst in TJ-II
#> (the largest trapped fraction, so the most exposure to the trapped-channel
#> quadrature error) and weakest in QA and QH (quasi-symmetric, so small
#> bounce-averaged drift).  Measured, with * marking a turnover below 0.5 where
#> the case is not a real test:
#>
#>   case 1     miller  *vacuous   w7x  1.24e-01   iter *4.81e-01
#>              qa  5.40e-02       qh   1.39e+00   tjii  1.07e-01
#>              (w7x at nfield_periods = 8, which the other four still run; the
#>               w7x deck is now one field period and reads 2.5e-02)
#>
#> TJ-II is second best, not worst.  And QA and QH -- both quasi-symmetric, both
#> with small bounce-averaged drift -- differ by a factor of 26.  Two
#> configurations that share the property the explanation rests on cannot differ
#> by 26 because of that property, so the trapped-fraction account of case 1 is
#> wrong, whatever else is true.
#>
#> What the sweep does establish, across all six configurations:
#>
#>   - the nonlinear adiabatic and modified-adiabatic cases (7, 8) close
#>     everywhere: 1.2e-03 to 6.2e-02.  The core nonlinear machinery is sound in
#>     every geometry tried.
#>   - the linear collisional case (4) closes everywhere: 2.6e-03 to 1.8e-02.
#>   - the multi-species electromagnetic cases (3, 6) fail everywhere, 0.6 to
#>     3.2.  That is not a stellarator effect and not a W7-X peculiarity; it is
#>     the electromagnetic multi-species physics itself, which supports treating
#>     cases 3, 6 and the momentum channel of case 10 as one defect.
#>   - the multi-species electrostatic case (2) also failed in every geometry
#>     except QH, on the decks of the time (alpha0 = 0.7-style tubes, ginit
#>     'default').  What that was measuring turned out not to be the diagnostic
#>     -- see the case-2 deck -- and the four unasserted decks have not been
#>     revisited since.
CONFIGURATIONS = ('miller', 'w7x')

#> Below this, |P| is too small a rate of change of E_RH for a residual divided
#> by it to mean much.  Reported rather than asserted; see <_measure>.
DRIVE_FLOOR = 0.5
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
        return float('nan'), turnover, 0.0     # no drive at all
    big = magnitude > magnitude.max() / 100.0
    residual = float(np.median(np.abs(dEdt - P)[big] / magnitude[big]))

    #> Dividing by |P| is what makes this a relative error, and it is only
    #> meaningful while |P| is itself a significant rate of change of E_RH.  The
    #> mask above protects against individual small points -- it is relative to
    #> |P|'s own peak -- but not against |P| being globally weak, which is a real
    #> failure mode and not a hypothetical one.  Scanning nfield_periods in W7-X
    #> case 1 moves this residual over 1.5e-02 to 2.7e-01 while the diagnostic
    #> and the equilibrium are unchanged; most of that is the discretisation
    #> error growing with the ruggedness of the line (see the deck), but part
    #> of it is the drive varying sixfold with tube length under a discrepancy
    #> that does not.
    #>
    #> <drive> is the same construction as <turnover> but formed from P rather
    #> than from dE/dt: the fraction of E_RH that the accounted-for sources would
    #> move over the run.  Where it is small the residual is a ratio of two small
    #> numbers and should be read as such.
    drive = float(np.median(magnitude) * (t[-1] - t[0]) / E.mean())
    return residual, turnover, drive


VMEC_FILE = {'w7x':  'wout_w7x_standard.nc',
             'iter': 'wout_iter.nc',
             'qa':   'wout_QA.nc',
             'qh':   'wout_QH.nc',
             'tjii': 'wout_tjii.nc'}


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

    phi_residual, phi_turnover, phi_drive = _measure(netcdf_file, 'phi')
    omega_residual, omega_turnover, omega_drive = _measure(netcdf_file, 'omega')
    print(f'\n  -->  case {case} ({stem}), {configuration}: '
          f'phi {phi_residual:.2e} (turnover {phi_turnover:.2f}, drive {phi_drive:.2f}), '
          f'Omega {omega_residual:.2e} (turnover {omega_turnover:.2f}, drive {omega_drive:.2f})')
    #> Not an assertion: a weak drive does not make a case wrong, it makes its
    #> residual hard to read, and saying so beside the number is worth more than
    #> failing on it.
    for which, res, dr in (('phi', phi_residual, phi_drive),
                           ('Omega', omega_residual, omega_drive)):
        if dr < DRIVE_FLOOR and res == res:
            print(f'       NOTE: {which} drive is {dr:.2f}, below {DRIVE_FLOOR}; '
                  f'its residual divides by a |P| this small and is not a clean '
                  f'measure of the diagnostic')

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
