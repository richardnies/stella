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
#> Reported apart because the two populations do not stand on the same footing.
#> A trapped orbit lies wholly inside the simulated tube, so its average is the
#> orbit average whatever the tube is.  A passing orbit is averaged along the
#> tube, and on an irrational surface the field line never closes: the true
#> average is over the flux surface, and the tube gives an artefact of the
#> flux-tube construction instead.  On a rational surface, with the tube spanning
#> the closed line, the two coincide and the passing drive is physical.
#> 'RH_fluxes_drift' is the name written before the split; it is still read so
#> that earlier output files load.
RH_FLUX_VARIABLES_DRIFT = ['RH_fluxes_drift', 'RH_fluxes_drift_trapped', 'RH_fluxes_drift_passing']
RH_FLUX_VARIABLES_DRIFT_TRAPPED = ['RH_fluxes_drift_trapped']
RH_FLUX_VARIABLES_DRIFT_PASSING = ['RH_fluxes_drift_passing']


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
    RH_fluxes_drift_trapped = summed_fluxes(RH_FLUX_VARIABLES_DRIFT_TRAPPED)
    RH_fluxes_drift_passing = summed_fluxes(RH_FLUX_VARIABLES_DRIFT_PASSING)

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
    RH_fluxes_drift_trapped = RH_fluxes_drift_trapped[:, finite_kx]
    RH_fluxes_drift_passing = RH_fluxes_drift_passing[:, finite_kx]

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
    P_RH_drift_trapped = power(RH_fluxes_drift_trapped)
    P_RH_drift_passing = power(RH_fluxes_drift_passing)

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
    P_drift_trapped = P_RH_drift_trapped[interior].sum(axis=1)
    P_drift_passing = P_RH_drift_passing[interior].sum(axis=1)

    window = np.ones_like(time, dtype=bool)
    if time_min is not None: window &= time >= time_min
    if time_max is not None: window &= time <= time_max

    return (time[window], E_RH_total[window], dE_RH_dt[window],
            (P_nonlinear + P_collisional + P_drift)[window],
            P_nonlinear[window], P_collisional[window], P_drift[window],
            P_drift_trapped[window], P_drift_passing[window])


def budget_residual(netcdf_file, time_min=None, time_max=None, channel='total', kx_max=None):
    '''Relative L2 mismatch for the whole budget or for one channel.

    channel='total'      dE_RH/dt                                against P_RH
    channel='nonlinear'  dE_RH/dt - P_collisional - P_drift      against P_nonlinear
    channel='drift'      dE_RH/dt - P_collisional - P_nonlinear  against P_drift
    '''
    _, _, dE_RH_dt, P_RH, P_nonlinear, P_collisional, P_drift, _, _ = get_rh_budget(
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


################################################################################
#           THE PARALLEL-FLOW ROSENBLUTH-HINTON BUDGET                         #
################################################################################
# The sigma-odd member of the same family of projections: RH_upar is annihilated
# by parallel streaming and the non-secular radial drift exactly as RH_phi_I is,
# and RH_upar_inertia is the same projection applied to a unit-flow shifted
# Maxwellian, so their ratio is the residual parallel flow.  The energy is then
# the parallel kinetic energy that flow carries,
#
#     E_uRH = sum_s (1/2) m_s n_s |<RH_upar_s>/<RH_upar_inertia_s>|^2 ,
#
# and, the inertia being time independent,
#
#     dE_uRH/dt = sum_s m_s n_s Re[ conj(<RH_upar_s>) d<RH_upar_s>/dt ] / |<I_u,s>|^2
#               = sum_channels P_uRH .
#
# See DOCUMENTATION/RH_parallel_flow for the derivation.
################################################################################

RH_PMOM_FLUX_NONLINEAR   = 'RH_pmom_flux_nonlinear'
RH_PMOM_FLUX_COLLISIONAL = 'RH_pmom_flux_collisional'
RH_PMOM_FLUX_DRIFT       = 'RH_pmom_flux_drift'


def _field_line_average_per_species(ncdata, name, weight):
    '''As _field_line_average, but keeping the species axis.

    The flow energy weights each species by its own mass and normalises by its
    own inertia, so the species cannot be summed before those are applied.
    Returns an array with axes (time, species, kx), or (species, kx) for a
    variable with no time axis.
    '''
    if name not in ncdata.variables:
        return None
    array, dimensions = _complex(ncdata, name)
    array = np.tensordot(array, weight, axes=([dimensions.index('zed')], [0]))
    dimensions = [d for d in dimensions if d != 'zed']
    #> The nonlinear flux carries a ky axis, because the drive of a zonal mode
    #> is the ky beat rather than anything on the ky = 0 row; summing over it is
    #> what recovers the zonal drive.  The collisional and drift fluxes are
    #> already zonal and have no ky axis.
    if 'ky' in dimensions:
        array = array.sum(axis=dimensions.index('ky'))
        dimensions.remove('ky')
    if 'tube' in dimensions:
        array = array.take(0, axis=dimensions.index('tube'))
        dimensions.remove('tube')
    #> Order the remaining axes as (t, species, kx), or (species, kx).
    order = [d for d in ('t', 'species', 'kx') if d in dimensions]
    array = np.transpose(array, [dimensions.index(d) for d in order])
    return array


def get_rh_pmom_budget(netcdf_file, time_min=None, time_max=None, kx_max=None):
    '''The toroidal-momentum RH budget.

    Returns (time, E_uRH, dE_uRH_dt, P_total, P_nonlinear, P_collisional, P_drift).
    '''
    ncdata = Dataset(netcdf_file)

    time = np.array(ncdata.variables['t'][:])
    kx = np.array(ncdata.variables['kx'][:])
    zed = np.array(ncdata.variables['zed'][:])
    jacobian = np.array(ncdata.variables['jacob'][:])[:, 0]

    weight = (zed[1] - zed[0]) * jacobian.copy()
    weight[-1] = 0.0
    weight = weight / weight.sum()

    upar = _field_line_average_per_species(ncdata, 'RH_pmom', weight)
    inertia = _field_line_average_per_species(ncdata, 'RH_pmom_inertia', weight)
    if upar is None or inertia is None:
        raise KeyError('this run did not write the RH toroidal-momentum diagnostics')

    mass = np.array(ncdata.variables['mass'][:])
    density = np.array(ncdata.variables['dens'][:])

    #> Drop kx = 0, which carries no zonal flow, and anything above kx_max.
    keep = np.abs(kx) > 1e-12
    if kx_max is not None:
        keep &= np.abs(kx) <= kx_max
    kx = kx[keep]
    upar = upar[..., keep]
    inertia = inertia[..., keep]

    weight_s = (mass * density)[None, :, None]
    inertia2 = np.abs(inertia)[None, :, :]**2

    E = 0.5 * weight_s * np.abs(upar)**2 / inertia2          # (t, species, kx)

    def power(name):
        flux = _field_line_average_per_species(ncdata, name, weight)
        if flux is None:
            return np.zeros(E.shape)
        flux = flux[..., keep]
        return -weight_s * np.real(1j * kx[None, None, :] * flux * np.conj(upar)) / inertia2

    P_nl   = power(RH_PMOM_FLUX_NONLINEAR)
    P_coll = power(RH_PMOM_FLUX_COLLISIONAL)
    P_dr   = power(RH_PMOM_FLUX_DRIFT)

    E_total = E.sum(axis=(1, 2))
    dE_dt = np.gradient(E_total, time)

    interior = slice(1, -1)
    time, E_total, dE_dt = time[interior], E_total[interior], dE_dt[interior]
    P_nl, P_coll, P_dr = (P[interior].sum(axis=(1, 2)) for P in (P_nl, P_coll, P_dr))

    window = np.ones_like(time, dtype=bool)
    if time_min is not None: window &= time >= time_min
    if time_max is not None: window &= time <= time_max

    return (time[window], E_total[window], dE_dt[window],
            (P_nl + P_coll + P_dr)[window],
            P_nl[window], P_coll[window], P_dr[window])
