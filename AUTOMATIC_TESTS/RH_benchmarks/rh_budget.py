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

# The RH flux channels written by stella, split by the physics that drives them.
# Absent variables are skipped, so the same helper works for electrostatic and
# electromagnetic runs alike.
#
# The two are kept apart because they are verified differently.  The collisional
# channel is checked on its own by the linear benchmark, where it is the only
# source; in a nonlinear run it is a known, independently verified correction,
# and what needs testing is the nonlinear channel.
RH_FLUX_VARIABLES_NONLINEAR = [
    'RH_fluxes_phi_even', 'RH_fluxes_phi_odd',
    'RH_fluxes_apar_even', 'RH_fluxes_apar_odd',
    'RH_fluxes_bpar_even', 'RH_fluxes_bpar_odd',
]
RH_FLUX_VARIABLES_COLLISIONAL = ['RH_fluxes_collisional']

#> The bounce-averaged radial magnetic drift.  Absent from a tokamak, where
#> quasisymmetry makes <v_Mx>_b vanish, and the leading drive in a general
#> stellarator.  Older output files predate it, so it is read optionally.
RH_FLUX_VARIABLES_DRIFT = ['RH_fluxes_drift']


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


def get_rh_budget(netcdf_file, time_min=None, time_max=None, kx_max=None):
    '''Return (time, E_RH, dE_RH/dt, P_RH, P_RH_nonlinear, P_RH_collisional,
    P_RH_drift), all summed over kx.

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

    def summed_fluxes(names):
        total = np.zeros_like(RH_phi_I)
        for name in names:
            contribution = _field_line_average(ncdata, name, weight)
            if contribution is not None:
                total = total + contribution
        return total

    RH_fluxes_nonlinear = summed_fluxes(RH_FLUX_VARIABLES_NONLINEAR)
    RH_fluxes_collisional = summed_fluxes(RH_FLUX_VARIABLES_COLLISIONAL)
    RH_fluxes_drift = summed_fluxes(RH_FLUX_VARIABLES_DRIFT)

    # kx = 0 carries no zonal-flow energy: 1-Gamma0 and the RH inertia both
    # vanish there, so the energy is 0/0.  Drop it.
    finite_kx = np.abs(kx) > 1e-12
    if kx_max is not None:
        finite_kx &= np.abs(kx) <= kx_max
    kx = kx[finite_kx]
    RH_phi_I = RH_phi_I[:, finite_kx]
    RH_inertia = RH_inertia[finite_kx]
    RH_fluxes_nonlinear = RH_fluxes_nonlinear[:, finite_kx]
    RH_fluxes_collisional = RH_fluxes_collisional[:, finite_kx]
    RH_fluxes_drift = RH_fluxes_drift[:, finite_kx]

    #> The prefactor of eq (19), sum_s Z_s^2 e^2 n_s / T_s * <1 - Gamma_0s>_psi,
    #> summed over species.  Gamma_0s = I0(b_s) exp(-b_s) with b_s = kperp^2
    #> rho_s^2, and rho_s / rho_ref = sqrt(m_s T_s) / Z_s in stella's
    #> normalisation, so b_s is the reference b scaled by m_s T_s / Z_s^2.
    #>
    #> This cancels between E_RH and P_RH and so does not affect the residual
    #> the benchmarks assert, but it sets the absolute value of E_RH, which is
    #> what gets compared against the zonal-flow energy.
    charge = np.array(ncdata.variables['charge'][:])
    mass = np.array(ncdata.variables['mass'][:])
    temperature = np.array(ncdata.variables['temp'][:])
    density = np.array(ncdata.variables['dens'][:])

    b_reference = (kx[:, None] / bmag[None, :])**2 * gds22[None, :]
    polarisation = np.zeros_like(kx)
    for z_s, m_s, T_s, n_s in zip(charge, mass, temperature, density):
        b_s = b_reference * (m_s * T_s / z_s**2)
        Gamma0_s = (weight[None, :] * np.i0(b_s / 2) * np.exp(-b_s / 2)).sum(axis=1)
        polarisation += z_s**2 * n_s / T_s * (1 - Gamma0_s)

    prefactor = polarisation[None, :] / np.abs(RH_inertia)[None, :]**2
    E_RH = np.abs(RH_phi_I)**2 / (2 * np.abs(RH_inertia)[None, :]**2) * polarisation[None, :]

    def power(fluxes):
        return -np.real(1j * kx[None, :] * fluxes * np.conj(RH_phi_I)) * prefactor

    P_RH_nonlinear = power(RH_fluxes_nonlinear)
    P_RH_collisional = power(RH_fluxes_collisional)
    P_RH_drift = power(RH_fluxes_drift)

    # np.gradient rather than a fixed-step difference: a nonlinear run may adapt
    # delt, so the time axis is not guaranteed to be uniformly spaced.  Drop the
    # end points, where np.gradient falls back to a one-sided difference.
    E_RH_total = E_RH.sum(axis=1)
    dE_RH_dt = np.gradient(E_RH_total, time)

    interior = slice(1, -1)
    time, E_RH_total, dE_RH_dt = time[interior], E_RH_total[interior], dE_RH_dt[interior]
    P_nonlinear = P_RH_nonlinear[interior].sum(axis=1)
    P_collisional = P_RH_collisional[interior].sum(axis=1)
    P_drift = P_RH_drift[interior].sum(axis=1)

    window = np.ones_like(time, dtype=bool)
    if time_min is not None: window &= time >= time_min
    if time_max is not None: window &= time <= time_max

    return (time[window], E_RH_total[window], dE_RH_dt[window],
            (P_nonlinear + P_collisional + P_drift)[window],
            P_nonlinear[window], P_collisional[window], P_drift[window])


def budget_residual(netcdf_file, time_min=None, time_max=None, channel='total', kx_max=None):
    '''Relative L2 mismatch for the whole budget or for one channel.

    channel='total'      dE_RH/dt                                against P_RH
    channel='nonlinear'  dE_RH/dt - P_collisional - P_drift      against P_nonlinear
    channel='drift'      dE_RH/dt - P_collisional - P_nonlinear  against P_drift
    '''
    _, _, dE_RH_dt, P_RH, P_nonlinear, P_collisional, P_drift = get_rh_budget(
        netcdf_file, time_min, time_max, kx_max)

    if channel == 'nonlinear':
        measured, expected = dE_RH_dt - P_collisional - P_drift, P_nonlinear
    elif channel == 'drift':
        measured, expected = dE_RH_dt - P_collisional - P_nonlinear, P_drift
    else:
        measured, expected = dE_RH_dt, P_RH

    norm = np.linalg.norm(expected)
    if norm == 0.0:
        return np.inf
    return np.linalg.norm(measured - expected) / norm


def field_line_averaged_rh_inertia(netcdf_file):
    '''The dl/B-averaged RH inertia, summed over species, against kx.

    Time-independent, and the quantity every other RH result is scaled by, so it
    is the sharpest thing to compare between two runs of the same physics.
    '''
    ncdata = Dataset(netcdf_file)
    zed = np.array(ncdata.variables['zed'][:])
    jacobian = np.array(ncdata.variables['jacob'][:])[:, 0]
    weight = (zed[1] - zed[0]) * jacobian.copy()
    weight[-1] = 0.0
    weight = weight / weight.sum()
    return _field_line_average(ncdata, 'RH_inertia', weight)
