"""Per-configuration Rosenbluth-Hinton budget figures.

Top panel:    E_RH, the zonal-flow energy as reconstructed from the RH
              projection, against E_phi, the same energy formed directly from
              the potential.  RH_phi_I is by construction the RH inertia times
              phi, so the two curves coincide when the projection is faithful;
              anywhere they part company is a problem upstream of the budget.

Bottom panel: dE_RH/dt against the summed P_RH and its individual channels.  A
              closed budget puts dE_RH/dt on top of the total.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from netCDF4 import Dataset

sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent))
from rh_budget import get_rh_budget, _complex, _field_line_average


def energy_from_phi(netcdf_file):
    """The zonal-flow energy built from phi rather than from the projection.

    Same expression as E_RH -- the polarisation times |phi|^2 / 2, field-line
    averaged and summed over the finite kx -- so that any difference between the
    two is attributable to the projection and not to a change of definition.
    Returns None when the run did not write phi.
    """
    ncdata = Dataset(netcdf_file)
    if 'phi_vs_t' not in ncdata.variables:
        return None, None
    time = np.array(ncdata.variables['t'][:])
    kx = np.array(ncdata.variables['kx'][:])
    zed = np.array(ncdata.variables['zed'][:])
    jacobian = np.array(ncdata.variables['jacob'][:])[:, 0]
    bmag = np.array(ncdata.variables['bmag'][:])[:, 0]
    shat = float(np.array(ncdata.variables['shat'][...]))
    gds22 = np.array(ncdata.variables['gds22'][:])[:, 0] / shat**2

    weight = (zed[1] - zed[0]) * jacobian.copy()
    weight[-1] = 0.0
    weight = weight / weight.sum()

    phi, dims = _complex(ncdata, 'phi_vs_t')          # (t, tube, zed, kx, ky)
    phi = np.tensordot(phi, weight, axes=([dims.index('zed')], [0]))
    remaining = [d for d in dims if d != 'zed']
    #> Take the first tube and the zonal ky, leaving (t, kx).
    for axis_name in ('tube', 'ky'):
        if axis_name in remaining:
            i = remaining.index(axis_name)
            phi = np.take(phi, 0, axis=i)
            remaining.pop(i)
    charge = np.array(ncdata.variables['charge'][:])
    mass = np.array(ncdata.variables['mass'][:])
    temperature = np.array(ncdata.variables['temp'][:])
    density = np.array(ncdata.variables['dens'][:])

    finite = np.abs(kx) > 1e-12
    kx_f = kx[finite]
    b_reference = (kx_f[:, None] / bmag[None, :])**2 * gds22[None, :]
    polarisation = np.zeros_like(kx_f)
    for z_s, m_s, T_s, n_s in zip(charge, mass, temperature, density):
        b_s = b_reference * (m_s * T_s / z_s**2)
        Gamma0_s = (weight[None, :] * np.i0(b_s / 2) * np.exp(-b_s / 2)).sum(axis=1)
        polarisation += z_s**2 * n_s / T_s * (1 - Gamma0_s)

    phi = phi[..., finite]
    return time, (np.abs(phi)**2 * polarisation[None, :] / 2).sum(axis=-1)


def figure(netcdf_file, title, outfile, time_min=3.0):
    time, E, dEdt, P, P_nl, P_coll, P_drift = get_rh_budget(netcdf_file)
    t_phi, E_phi = energy_from_phi(netcdf_file)

    fig, (ax_e, ax_p) = plt.subplots(2, 1, figsize=(7.2, 6.4), sharex=True,
                                     gridspec_kw=dict(height_ratios=[1, 1.25], hspace=0.12))

    ax_e.plot(time, E, color='#1b3a5c', lw=1.7, label=r'$E_{\rm RH}$  (from the RH projection)')
    if E_phi is not None:
        ax_e.plot(t_phi, E_phi, color='#c2703a', lw=1.1, ls='--', label=r'$E_\phi$  (from $\phi$ directly)')
    ax_e.set_ylabel('zonal-flow energy')
    ax_e.legend(frameon=False, fontsize=8.5, loc='best')
    ax_e.grid(alpha=0.15, lw=0.6)
    ax_e.set_title(title, fontsize=11, loc='left')

    ax_p.plot(time, dEdt, color='#1b3a5c', lw=1.9, label=r'$dE_{\rm RH}/dt$  (measured)')
    ax_p.plot(time, P, color='#c2703a', lw=1.3, ls='--', label=r'$\sum P_{\rm RH}$  (predicted)')
    for series, colour, label in ((P_drift, '#2e7d6b', r'$P_{\rm RH}$ drift'),
                                  (P_coll, '#8a6d1f', r'$P_{\rm RH}$ collisional'),
                                  (P_nl, '#7d3c6b', r'$P_{\rm RH}$ nonlinear')):
        if np.any(series != 0):
            ax_p.plot(time, series, color=colour, lw=0.9, alpha=0.8, label=label)
    ax_p.axhline(0.0, color='k', lw=0.5, alpha=0.3)
    ax_p.set_xlabel(r'time  $[a/v_{\rm th}]$')
    ax_p.set_ylabel('power into the zonal flow')
    ax_p.legend(frameon=False, fontsize=8.5, loc='best', ncol=2)
    ax_p.grid(alpha=0.15, lw=0.6)

    window = time >= time_min
    if window.sum() > 2 and np.any(P[window] != 0):
        slope = (dEdt[window] * P[window]).sum() / (P[window] * P[window]).sum()
        corr = np.corrcoef(dEdt[window], P[window])[0, 1]
        ax_p.text(0.985, 0.04, f'$dE/dt$ on $\\sum P$:  slope {slope:.3f},  corr {corr:.4f}'
                              f'   (for $t>{time_min:g}$)',
                  transform=ax_p.transAxes, ha='right', va='bottom', fontsize=8.5, color='#444')

    fig.savefig(outfile, dpi=150, bbox_inches='tight')
    plt.close(fig)
    return outfile


if __name__ == '__main__':
    for spec in sys.argv[1:]:
        path, title, out = spec.split('::')
        print(figure(path, title, out))
