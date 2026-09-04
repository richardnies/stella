################################################################################
#              PLOT THE ROSENBLUTH-HINTON ENERGY BUDGET FOR A RUN              #
################################################################################
# Produces, for each run given on the command line, a two-panel figure:
#
#   top    E_RH(t), the zonal-flow energy carried by the RH response
#   bottom dE_RH/dt against sum_kx P_RH -- the budget the benchmarks assert --
#          with the nonlinear and collisional channels shown separately and the
#          running relative residual on a twin axis
#
# Usage:
#     python plot_rh_budget.py <run.out.nc> [<run.out.nc> ...] [--out DIR]
#                              [--time-min T] [--time-max T] [--kx-max K]
#
# This is a diagnostic aid for the benchmarks, not part of the pytest run.
################################################################################

import argparse
import pathlib
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from rh_budget import get_rh_budget


def plot_one(netcdf_file, out_dir, time_min=None, time_max=None, kx_max=None):
    time, E_RH, dE_RH_dt, P_RH, P_nonlinear, P_collisional, P_drift, P_drift_tr, P_drift_pa = get_rh_budget(
        netcdf_file, time_min, time_max, kx_max)

    residual = np.linalg.norm(dE_RH_dt - P_RH) / np.linalg.norm(P_RH)
    scale = np.maximum(np.abs(P_RH), np.abs(dE_RH_dt))
    local = np.abs(dE_RH_dt - P_RH) / np.where(scale > 0, scale, np.nan)

    fig, (ax_energy, ax_budget) = plt.subplots(
        2, 1, figsize=(9, 7), sharex=True, gridspec_kw={'height_ratios': [1, 1.4]})

    name = pathlib.Path(netcdf_file).name.replace('.out.nc', '')
    ax_energy.plot(time, E_RH, color='#1f4e79', lw=1.8)
    ax_energy.set_yscale('log')
    ax_energy.set_ylabel(r'$E_{\rm RH}$')
    ax_energy.set_title(f'{name}    (relative $L_2$ residual = {residual:.2e})')
    ax_energy.grid(alpha=0.3)

    ax_budget.plot(time, dE_RH_dt, color='#1f4e79', lw=2.2, label=r'$dE_{\rm RH}/dt$')
    ax_budget.plot(time, P_RH, color='#d1495b', lw=1.2, ls='--', label=r'$\sum_{k_x} P_{\rm RH}$')
    if np.any(P_nonlinear != 0):
        ax_budget.plot(time, P_nonlinear, color='#d1495b', lw=0.9, alpha=0.55, label=r'$P_{\rm RH}$ nonlinear')
    if np.any(P_drift != 0):
        ax_budget.plot(time, P_drift, color='#1f6d8a', lw=0.9, alpha=0.75, label=r'$P_{\rm RH}$ drift')
    if np.any(P_collisional != 0):
        ax_budget.plot(time, P_collisional, color='#8a6d1f', lw=0.9, alpha=0.75, label=r'$P_{\rm RH}$ collisional')
    # symlog: the linear case oscillates through zero, the nonlinear cases grow
    # exponentially over many decades, and this reads correctly for both.
    finite = np.abs(np.concatenate([dE_RH_dt, P_RH, P_nonlinear, P_collisional, P_drift]))
    finite = finite[finite > 0]
    if finite.size:
        ax_budget.set_yscale('symlog', linthresh=max(finite.min(), finite.max() * 1e-6))
    ax_budget.set_xlabel(r'$t \, v_{\rm th}/a$')
    ax_budget.set_ylabel('power into the zonal flow')
    ax_budget.legend(loc='upper left', frameon=False, fontsize=8)
    ax_budget.grid(alpha=0.3)

    ax_residual = ax_budget.twinx()
    ax_residual.plot(time, local, color='#8d8d8d', lw=0.9, alpha=0.8)
    ax_residual.set_yscale('log')
    ax_residual.set_ylabel('local relative residual', color='#8d8d8d')
    ax_residual.tick_params(axis='y', colors='#8d8d8d')

    fig.tight_layout()
    output = pathlib.Path(out_dir) / f'rh_budget_{name}.png'
    fig.savefig(output, dpi=150)
    plt.close(fig)
    print(f'{name}: residual = {residual:.3e}  ->  {output}')
    return residual


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('netcdf_files', nargs='+')
    parser.add_argument('--out', default='.')
    parser.add_argument('--time-min', type=float, default=None)
    parser.add_argument('--time-max', type=float, default=None)
    parser.add_argument('--kx-max', type=float, default=None)
    args = parser.parse_args()

    pathlib.Path(args.out).mkdir(parents=True, exist_ok=True)
    for netcdf_file in args.netcdf_files:
        plot_one(netcdf_file, args.out, args.time_min, args.time_max, args.kx_max)
