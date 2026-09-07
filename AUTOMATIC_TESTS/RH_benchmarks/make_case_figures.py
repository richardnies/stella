#!/usr/bin/env python3
'''Generate one budget figure per physics case, for every configuration run.

Usage:  python3 make_case_figures.py <run directory> <figure directory>

The run directory holds the output of test_3_physics_cases/*.in, one
<config>_<n>_<name>.out.nc per deck.  Each case becomes a single figure with
one row per configuration and one column per invariant, in the shape of the
budget figure the report already uses.
'''
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import plot_rh_report as P

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

CONFIGURATIONS = (('miller', 'Miller tokamak'), ('w7x', 'W7-X'))


def main(rundir, figdir):
    rundir, figdir = pathlib.Path(rundir), pathlib.Path(figdir)
    figdir.mkdir(parents=True, exist_ok=True)
    for case, (stem, title, window) in sorted(CASES.items()):
        runs = []
        for tag, label in CONFIGURATIONS:
            nc = rundir / f'{tag}_{case}_{stem}.out.nc'
            if nc.exists():
                runs.append((label, nc, window))
            else:
                print(f'  missing: {nc.name}')
        if not runs:
            continue
        out = figdir / f'fig_case{case}_{stem}.pdf'
        P.figure_case_budget(runs, out, case_title=title)
        print('wrote', out.name)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else '.')
