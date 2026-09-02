################################################################################
#            ROSENBLUTH-HINTON ENERGY BUDGET:  dE_RH/dt  vs  sum P_RH          #
################################################################################
# Shared helper for the RH budget benchmarks.
#
# The zonal-flow energy carried by the Rosenbluth-Hinton response and the power
# transferred into it are, per radial wavenumber and with <.> the dl/B
# field-line average,
#
#     E_RH(t,kx) = |<RH_phi_I>|^2 / (2 |<RH_inertia>|^2) * (1 - Gamma0)
#     P_RH(t,kx) = -Re[ i kx F <RH_phi_I>* ] / |<RH_inertia>|^2 * (1 - Gamma0)
#
# where F is the sum of every RH flux written to the netCDF file (phi/apar/bpar,
# even and odd in vpa, plus the collisional flux).  If the diagnostic is
# consistent then
#
#     d E_RH / dt  =  sum_kx P_RH
#
# and the benchmarks assert exactly that.  The test needs no reference data:
# both sides come from the same run, so it checks the diagnostic against the
# code's own time evolution.
#
# These definitions match stella_diagnostics_v2
# (stella_diagnostics/physics/rosenbluth_hinton.py: get_E_RH_t_kx / get_P_RH),
# so the two must be kept in step.
################################################################################

import numpy as np
from netCDF4 import Dataset

# Every RH flux channel written by stella.  Absent variables are skipped, so the
# same helper works for electrostatic and electromagnetic runs alike.
RH_FLUX_VARIABLES = [
    'RH_fluxes_phi_even', 'RH_fluxes_phi_odd',
    'RH_fluxes_apar_even', 'RH_fluxes_apar_odd',
    'RH_fluxes_bpar_even', 'RH_fluxes_bpar_odd',
    'RH_fluxes_collisional',
]


def _complex(ncdata, name):
    '''Read a stella complex variable, returning the array and its axis names.'''
    var = ncdata.variables[name]
    values = np.array(var[:])
    dimensions = list(var.dimensions)
    if dimensions[-1] != 'ri':
        raise ValueError(f'{name} is not a complex stella variable')
    return values[..., 0] + 1j * values[..., 1], dimensions[:-1]


def _field_line_average(ncdata, name, weight):
    '''dl/B average over zed, summed over species and ky, tube index dropped.'''
    if name not in ncdata.variables:
        return None
    array, dimensions = _complex(ncdata, name)
    array = np.tensordot(array, weight, axes=([dimensions.index('zed')], [0]))
    dimensions = [d for d in dimensions if d != 'zed']
    for axis in ('species', 'ky'):
        if axis in dimensions:
            array = array.sum(axis=dimensions.index(axis))
            dimensions.remove(axis)
    if 'tube' in dimensions:
        array = array.take(0, axis=dimensions.index('tube'))
        dimensions.remove('tube')
    return array


def get_rh_budget(netcdf_file):
    '''Return (time, E_RH, dE_RH/dt, sum P_RH) summed over kx.

    dE_RH/dt is a centred difference, so it is defined on the interior points;
    P_RH is returned on the same points.
    '''
    ncdata = Dataset(netcdf_file)

    time = np.array(ncdata.variables['t'][:])
    kx = np.array(ncdata.variables['kx'][:])
    zed = np.array(ncdata.variables['zed'][:])
    jacobian = np.array(ncdata.variables['jacob'][:])[:, 0]
    bmag = np.array(ncdata.variables['bmag'][:])[:, 0]
    shat = float(np.array(ncdata.variables['shat'][...]))
    gds22 = np.array(ncdata.variables['gds22'][:])[:, 0] / shat**2

    # stella's dl_over_b: delzed*jacob with the duplicated endpoint dropped,
    # normalised so that the field-line average of unity is unity.
    weight = (zed[1] - zed[0]) * jacobian.copy()
    weight[-1] = 0.0
    weight = weight / weight.sum()

    RH_phi_I = _field_line_average(ncdata, 'RH_phi_I', weight)
    RH_inertia = _field_line_average(ncdata, 'RH_inertia', weight)

    RH_fluxes = np.zeros_like(RH_phi_I)
    for name in RH_FLUX_VARIABLES:
        contribution = _field_line_average(ncdata, name, weight)
        if contribution is not None:
            RH_fluxes = RH_fluxes + contribution

    # kx = 0 carries no zonal-flow energy: 1-Gamma0 and the RH inertia both
    # vanish there, so the energy is 0/0.  Drop it.
    finite_kx = np.abs(kx) > 1e-12
    kx = kx[finite_kx]
    RH_phi_I = RH_phi_I[:, finite_kx]
    RH_inertia = RH_inertia[finite_kx]
    RH_fluxes = RH_fluxes[:, finite_kx]

    # Gamma0 = <I0(b) exp(-b)>, with b = kperp^2 rho^2 evaluated along the field line
    b = (kx[:, None] / bmag[None, :])**2 * gds22[None, :]
    Gamma0 = (weight[None, :] * np.i0(b / 2) * np.exp(-b / 2)).sum(axis=1)

    prefactor = (1 - Gamma0)[None, :] / np.abs(RH_inertia)[None, :]**2
    E_RH = np.abs(RH_phi_I)**2 / (2 * np.abs(RH_inertia)[None, :]**2) * (1 - Gamma0)[None, :]
    P_RH = -np.real(1j * kx[None, :] * RH_fluxes * np.conj(RH_phi_I)) * prefactor

    dt = time[1] - time[0]
    dE_RH_dt = (E_RH[2:] - E_RH[:-2]) / (2 * dt)

    return time[1:-1], E_RH.sum(axis=1), dE_RH_dt.sum(axis=1), P_RH[1:-1].sum(axis=1)


def budget_residual(netcdf_file):
    '''Relative L2 mismatch between dE_RH/dt and sum P_RH.'''
    _, _, dE_RH_dt, P_RH = get_rh_budget(netcdf_file)
    norm = np.linalg.norm(P_RH)
    if norm == 0.0:
        return np.inf
    return np.linalg.norm(dE_RH_dt - P_RH) / norm
