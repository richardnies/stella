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

#> Only Miller and W7-X are asserted.  The other four equilibria -- ITER, QA,
#> QH and TJ-II -- ship all eleven decks too, each the W7-X deck with the
#> equilibrium swapped (ITER's cases 1 and 2 also take nfield_periods = 2.505,
#> one poloidal turn, because one field period of an axisymmetric field is 0.4
#> of a turn and tests nothing).  The suite does not run them: forty-four runs
#> of up to ten minutes is more than a test suite should cost, and setting
#> forty-four bounds from one measurement each is bookkeeping dressed up as
#> verification.  make_case_figures.py measures and plots them with the same
#> statistic when their output is present; every deck header carries its last
#> measurement, and the report collects them (DOCUMENTATION/stella_RH_report,
#> "The four remaining configurations").
#>
#> What the four establish, read against the two asserted configurations:
#>
#>   - Nothing is special to Miller or W7-X.  The even invariant closes at
#>     1e-3 to 1e-1 wherever it is driven with adiabatic, collisional or
#>     turbulent electrons (cases 1, 4, 5, 7, 8, 9) and fails only with two
#>     kinetic species at finite beta (cases 3, 6, and phi_RH in 10/11 outside
#>     QA), in all six.  The dA_par transient is not a stellarator effect.
#>   - The odd invariant fails on every tube that does not close, in every
#>     equilibrium: case 1 reads 0.3 to 3 on the one-period alpha0 = 0.7 tubes
#>     (W7-X 2.9, QH 0.32, TJ-II 0.54, QA 1.0) and 1e1 to 1e3 with kinetic
#>     electrons (case 3).  ITER on a closed poloidal turn conserves it like
#>     Miller.  With kinetic electrons on an open tube E_Omega is pumped even
#>     with collisions (cases 9 and 10 in ITER, QH and W7-X, in the quiet phase
#>     before turbulence): the case-2 mechanism on the odd invariant.
#>   - The collisionless kinetic-electron growth of case 2 (see the w7x_2 deck)
#>     is strongest where the join is most mismatched: on the old eight-period
#>     tubes it blew up in QA and TJ-II (|grad x|^2 ratios 0.08 and 0.18
#>     across the join), read 7.9 in ITER (0.67) and 0.44 in QH (0.81), and it
#>     is absent on every closed tube.  All case-2 decks now sit on
#>     alpha0 = 0, where QA and QH have no source at all (quasi-symmetry: the
#>     drift channel cancels to 1e-19) and carry E(T)/E(0) instead of a
#>     residual; TJ-II is driven there and reads 0.53, halving at nzed = 256.
#>   - Three decks blow up: QA case 9 and TJ-II cases 8 and 9, under-resolved
#>     turbulence at a/L_T = 8 with no hyper-dissipation and twelve modes each
#>     way.  Their numbers are measured on the blow-up and say nothing about
#>     the diagnostic; the decks are kept identical to W7-X's on purpose.
CONFIGURATIONS = ('miller', 'w7x')

#> The potential budget is measured from the projection of g_s
#> ("RH_phi_I_g"), which is the one that equals the relaxed potential, and its
#> power sum carries the induction channel that the three fluxes -- which come
#> from the evolution equation for gbar_s -- do not account for.  See the header
#> of rh_budget.py.  Every electrostatic case is unaffected: there the two
#> projections are the same array and the induction channel is identically zero,
#> so cases 1, 2, 4, 5, 7, 8 and 9 read exactly what they read before.  Of the
#> four electromagnetic cases the residual is a wash -- Miller case 6 improves
#> from 5.37 to 9.30e-01 and W7-X case 3 from 6.35 to 5.67, W7-X cases 10 and 11
#> worsen slightly from 6.77e-02 to 7.50e-02 and 2.91e-02 to 3.40e-02 -- because
#> the induction channel is the size of the residual already present.  What the
#> change buys is that the quantity asserted is the relaxed potential.
#>
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
    #> No source at all: every channel is identically zero (the Miller
    #> collisionless cases) or cancels to round-off (the quasi-symmetric
    #> tubes of case 2, where |P| sits at 1e-19 against a dE/dt of 1e-5).
    #> A residual divided by that is 1e+15 of nothing; the run is a
    #> conservation test instead and is reported as such.
    if magnitude.max() <= 1e-10 * np.abs(dEdt).max():
        return float('nan'), turnover, 0.0
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
