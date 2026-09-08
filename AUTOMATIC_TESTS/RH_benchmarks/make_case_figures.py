#!/usr/bin/env python3
'''Generate one budget figure per physics case, for every configuration run,
and the summary figure that puts every case on one axis.

Usage:  python3 make_case_figures.py <run directory> <figure directory>

The run directory holds the output of test_3_physics_cases/*.in, one
<config>_<n>_<name>.out.nc per deck.  Each case becomes one figure,
fig_case<n>_<name>.pdf, with a row for every configuration that has output
-- the two the suite asserts (Miller and W7-X) and the four it only runs
(ITER, QA, QH, TJ-II) -- and the two invariants as columns.  The summary,
fig_case_summary.pdf, is measured here with the same statistic the test
asserts -- it imports the test module for the statistic, the tolerances and
the bounds -- so the figure and the assertions cannot drift apart; it
carries the asserted configurations only.  fig_case_configurations.pdf puts
all six on one axis per invariant, with the energy turnover beside each bar
and a dagger on the runs whose time step collapsed.  Every row is also
printed as a table.
'''
import importlib.util
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import plot_rh_report as P

TEST_FILE = pathlib.Path(__file__).resolve().parent / 'test_3_physics_cases' / 'test_rh_physics_cases.py'


def load_test_module():
    '''The statistic, the tolerances and the known-failure bounds live in the
    test file, and this is the only copy of them.'''
    spec = importlib.util.spec_from_file_location('test_rh_physics_cases', TEST_FILE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

#> Case number -> (deck stem, caption title, window applied to the budget).
#>
#> The windows keep the last 40% of each run.  Before that the zonal flow has
#> not been driven, so both sides of the budget are near zero and their ratio is
#> a ratio of noise -- measured over the first third instead, these same runs
#> report residuals of order 100, which says nothing about the diagnostic.  The
#> fraction is used rather than a fixed time because the electromagnetic cases
#> are cut short by the CFL condition, reaching t = 9.6 where the electrostatic
#> ones reach 28.  The linear cases need no window: they are driven from the
#> first step.
CASES = {
    1: ('linear_collisionless', 'Case 1: linear, collisionless', {}),
    2: ('linear_collisionless_kinetic',
        r'Case 2: linear, collisionless, two kinetic species, electrostatic', {}),
    3: ('linear_collisionless_multi',
        r'Case 3: linear, collisionless, two kinetic species, electromagnetic', {}),
    4: ('linear_collisional', 'Case 4: linear, collisional', {}),
    5: ('linear_collisional_kinetic',
        r'Case 5: linear, collisional, two kinetic species, electrostatic', {}),
    6: ('linear_collisional_multi',
        r'Case 6: linear, collisional, two kinetic species, electromagnetic', {}),
    7: ('nl_modified_adiabatic',
        'Case 7: nonlinear, modified-adiabatic electrons', dict(window=0.4)),
    8: ('nl_adiabatic',
        'Case 8: nonlinear, adiabatic electrons', dict(window=0.4)),
    9: ('nl_kinetic',
        'Case 9: nonlinear, kinetic ions and electrons', dict(window=0.4)),
    10: ('nl_em_apar',
        r'Case 10: nonlinear electromagnetic, $\delta A_\parallel$ only',
        dict(window=0.4)),
    11: ('nl_em_apar_bpar',
        r'Case 11: nonlinear electromagnetic, $\delta A_\parallel$ and $\delta B_\parallel$',
        dict(window=0.4)),
}

#> Every configuration with decks, asserted or not, in the order the figures
#> show them.  Which of them the suite asserts is the test module's business
#> (its CONFIGURATIONS); here it only decides which rows carry a tolerance.
CONFIGURATIONS = (('miller', 'Miller tokamak'), ('w7x', 'W7-X'), ('iter', 'ITER'),
                  ('qa', 'QA'), ('qh', 'QH'), ('tjii', 'TJ-II'))


SHORT = {
    1: 'linear collisionless', 2: 'linear collisionless, 2sp, ES',
    3: 'linear collisionless, 2sp, EM', 4: 'linear collisional',
    5: 'linear collisional, 2sp, ES', 6: 'linear collisional, 2sp, EM',
    7: 'NL modified-adiabatic', 8: 'NL adiabatic', 9: 'NL kinetic',
    10: 'NL EM, dApar', 11: 'NL EM, both',
}


def summary_rows(rundir):
    '''One row per (case, configuration) present in <rundir>, measured as the
    test measures it: (label, phi residual, Omega residual, phi tested, Omega
    tested, phi bound, Omega bound, the two turnovers, the two drives, the
    two conservation ratios E(T)/E(0) for a run with no source, else None,
    the configuration tag, the case number, and whether the time step
    collapsed during the run).  A bound is the
    tolerance as a float, a (low, high) pair for a known failure, or None
    where nothing is asserted -- which is every row of a configuration the
    suite does not assert.'''
    T = load_test_module()
    rows = []
    for case, (stem, phi_tol, omega_tol, vacuous) in sorted(T.CASES.items()):
        for tag, label in CONFIGURATIONS:
            nc = pathlib.Path(rundir) / f'{tag}_{case}_{stem}.out.nc'
            if not nc.exists():
                continue
            phi_res, phi_to, phi_dr = T._measure(nc, 'phi')
            om_res, om_to, om_dr = T._measure(nc, 'omega')
            phi_ok = phi_to >= T.TURNOVER_FLOOR
            om_ok = om_to >= T.TURNOVER_FLOOR
            if tag in vacuous or tag not in T.CONFIGURATIONS:
                phi_bound = om_bound = None
            else:
                phi_bound = T.KNOWN_PHI_FAILURES.get((tag, case), phi_tol)
                om_bound = T.KNOWN_MOMENTUM_FAILURES.get((tag, case), omega_tol)
            #> A run with no source at all (sum P identically zero: the linear
            #> collisionless cases in Miller) has no residual, only a
            #> conservation ratio E(T)/E(0).  Carried so the figure can say so
            #> instead of leaving the row blank -- Miller case 3 turns the flow
            #> over twenty times with nothing driving it, and that is the one
            #> result in the set the test's residual statistic cannot see.
            phi_cons = conservation_ratio(nc, 'phi') if phi_res != phi_res else None
            om_cons = conservation_ratio(nc, 'omega') if om_res != om_res else None
            rows.append((f'{case:>2}  {SHORT[case]} / {label}',
                         phi_res, om_res, phi_ok, om_ok, phi_bound, om_bound,
                         phi_to, om_to, phi_dr, om_dr, phi_cons, om_cons,
                         tag, case, P.time_step_collapsed(nc)[0]))
    return rows


def conservation_ratio(netcdf_file, which):
    '''E(T)/E(0) over the whole run, end points included, for a run that has
    no source.  The budget functions drop the end points by default, and on
    the two-species decks the first step moves the projection by a few per
    cent, so the ratio taken from t_1 would miss exactly what this reports.'''
    T = load_test_module()
    fn = T.get_rh_budget if which == 'phi' else T.get_rh_omega_budget
    t, E = fn(netcdf_file, interior=False)[:2]
    return float(E[-1] / E[0])


def print_summary(rows):
    print(f'{"case / configuration":42s} {"phi":>9s} {"turn":>5s} {"drive":>5s}   '
          f'{"Omega":>9s} {"turn":>5s} {"drive":>5s}')
    for (label, phi_res, om_res, phi_ok, om_ok, phi_b, om_b,
         phi_to, om_to, phi_dr, om_dr, phi_cons, om_cons, tag, case, blew) in rows:
        print(f'{label:42s} {phi_res:9.2e} {phi_to:5.2f} {phi_dr:5.2f}{" " if phi_ok else "*"}  '
              f'{om_res:9.2e} {om_to:5.2f} {om_dr:5.2f}{" " if om_ok else "*"}'
              + (f'   no source: E(T)/E(0) = {phi_cons:.3f} / {om_cons:.3f}'
                 if phi_cons is not None and om_cons is not None else '')
              + ('   + time step collapsed' if blew else ''))
    print('* turnover below the floor: shown hollow, not a result')


def main(rundir, figdir):
    rundir, figdir = pathlib.Path(rundir), pathlib.Path(figdir)
    figdir.mkdir(parents=True, exist_ok=True)
    asserted = load_test_module().CONFIGURATIONS
    for case, (stem, title, window) in sorted(CASES.items()):
        runs = []
        for tag, label in CONFIGURATIONS:
            nc = rundir / f'{tag}_{case}_{stem}.out.nc'
            if nc.exists():
                runs.append((label, nc, window))
            else:
                print(f'  missing: {nc.name}')
        if runs:
            #> Six rows have to fit a page with room for a caption.
            out = figdir / f'fig_case{case}_{stem}.pdf'
            P.figure_case_budget(runs, out, case_title=title,
                                 row_height=1.75 if len(runs) > 2 else 2.9)
            print('wrote', out.name)
    rows = summary_rows(rundir)
    if rows:
        print_summary(rows)
        out = figdir / 'fig_case_summary.pdf'
        P.figure_case_summary([r for r in rows if r[13] in asserted], out)
        print('wrote', out.name)
        if any(r[13] not in asserted for r in rows):
            out = figdir / 'fig_case_configurations.pdf'
            P.figure_case_configurations(rows, [c for c, _ in CONFIGURATIONS],
                                         dict(CONFIGURATIONS), out)
            print('wrote', out.name)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else '.')
