"""Standalone figures for the Rosenbluth-Hinton report.

Produces every numerical check made on the two RH invariants -- the
potential-like phi_RH and the parallel-flow upar_RH -- as separate PDFs that can
be included directly in the write-up.

    python plot_rh_report.py <output-directory>
"""
import sys
import pathlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from netCDF4 import Dataset

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from rh_budget import (get_rh_budget, get_rh_upar_budget, _complex,
                       _field_line_average, _field_line_average_per_species)

PHI = '#1b3a5c'      # the potential-like invariant
UPA = '#b0532a'      # the parallel-flow invariant
GREY = '#6b6b6b'

plt.rcParams.update({'font.size': 9, 'axes.grid': True, 'grid.alpha': 0.15,
                     'grid.linewidth': 0.6, 'axes.axisbelow': True,
                     'figure.dpi': 150, 'savefig.bbox': 'tight'})


def _weights(ncdata):
    zed = np.array(ncdata.variables['zed'][:])
    jacobian = np.array(ncdata.variables['jacob'][:])[:, 0]
    weight = (zed[1] - zed[0]) * jacobian.copy()
    weight[-1] = 0.0
    return weight / weight.sum()


def invariants(netcdf_file):
    """Field-line-averaged <phi_RH> and <upar_RH>, each divided by its inertia."""
    ncdata = Dataset(netcdf_file)
    weight = _weights(ncdata)
    time = np.array(ncdata.variables['t'][:])

    phi = _field_line_average(ncdata, 'RH_phi_I', weight)
    inertia = _field_line_average(ncdata, 'RH_inertia', weight)
    upar = _field_line_average_per_species(ncdata, 'RH_upar', weight)
    upar_inertia = _field_line_average_per_species(ncdata, 'RH_upar_inertia', weight)

    kx = np.array(ncdata.variables['kx'][:])
    keep = np.abs(kx) > 1e-12
    #> One zonal mode is enough to show conservation; take the smallest finite kx.
    index = np.argmin(np.where(keep, np.abs(kx), np.inf))
    phi_ratio = phi[:, index] / inertia[index]
    upar_ratio = upar[:, 0, index] / upar_inertia[0, index]
    return time, phi_ratio, upar_ratio


def figure_conservation(netcdf_file, outfile, convergence=None):
    """Both invariants against time in a linear collisionless run.

    Streaming and the non-secular radial drift annihilate both exactly, so each
    trace must be flat; what is plotted is the departure from the initial value,
    which is the discretisation error alone.
    """
    time, phi_ratio, upar_ratio = invariants(netcdf_file)
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(7.6, 2.9),
                                  gridspec_kw=dict(width_ratios=[1.45, 1]))

    for series, colour, label in ((phi_ratio, PHI, r'$\langle\varphi_{\rm RH}\rangle/\langle I\rangle$'),
                                  (upar_ratio, UPA, r'$\langle u_{\parallel\rm RH}\rangle/\langle I_u\rangle$')):
        drift = np.abs(series - series[0]) / np.abs(series[0])
        ax.plot(time, drift, color=colour, lw=1.6, label=label)
    ax.set_yscale('log')
    ax.set_xlabel(r'time  $[a/v_{\rm th}]$')
    ax.set_ylabel('departure from initial value')
    ax.legend(frameon=False, fontsize=8.5, loc='lower right')
    ax.set_title('Conservation, linear and collisionless', fontsize=9.5, loc='left')

    if convergence is not None:
        labels, phi_values, upar_values = convergence
        x = np.arange(len(labels))
        ax2.plot(x, phi_values, 'o-', color=PHI, lw=1.5, ms=4)
        ax2.plot(x, upar_values, 's-', color=UPA, lw=1.5, ms=4)
        ax2.set_xticks(x)
        ax2.set_xticklabels(labels, fontsize=7.5)
        ax2.set_yscale('log')
        ax2.set_ylabel('drift over the run')
        ax2.set_title('Convergence under refinement', fontsize=9.5, loc='left')
    else:
        ax2.axis('off')

    fig.tight_layout()
    fig.savefig(outfile)
    plt.close(fig)
    return outfile


def figure_budget(netcdf_file, outfile, which='phi', title='', time_min=None, time_max=None):
    """The budget itself: energy on top, dE/dt against the summed power below.

    `which` selects the invariant.  A closed budget puts the measured derivative
    on top of the predicted total; the individual channels show which term
    carries it.
    """
    if which == 'phi':
        t, E, dEdt, P, P_nl, P_coll, P_dr = get_rh_budget(netcdf_file, time_min, time_max)[:7]
        colour, energy_label = PHI, r'$E_{\rm RH}$'
    else:
        t, E, dEdt, P, P_nl, P_coll, P_dr = get_rh_upar_budget(netcdf_file, time_min, time_max)
        colour, energy_label = UPA, r'$E_{u\rm RH}$'

    fig, (ax_e, ax_p) = plt.subplots(2, 1, figsize=(6.6, 5.0), sharex=True,
                                     gridspec_kw=dict(height_ratios=[1, 1.3], hspace=0.12))

    ax_e.plot(t, E, color=colour, lw=1.8, label=energy_label)
    ax_e.set_ylabel('zonal energy')
    ax_e.legend(frameon=False, fontsize=8.5)
    if title:
        ax_e.set_title(title, fontsize=9.5, loc='left')

    ax_p.plot(t, dEdt, color=colour, lw=1.9, label=r'$dE/dt$  (measured)')
    ax_p.plot(t, P, color='#c2703a' if which == 'phi' else '#2e6f8e',
              lw=1.3, ls='--', label=r'$\sum P$  (predicted)')
    for series, style, label in ((P_nl, '#7d3c6b', r'$P$ nonlinear'),
                                 (P_coll, '#8a6d1f', r'$P$ collisional'),
                                 (P_dr, '#2e7d6b', r'$P$ drift')):
        if np.any(series != 0):
            ax_p.plot(t, series, color=style, lw=0.9, alpha=0.85, label=label)
    ax_p.axhline(0.0, color='k', lw=0.5, alpha=0.3)
    ax_p.set_xlabel(r'time  $[a/v_{\rm th}]$')
    ax_p.set_ylabel('power into the zonal flow')
    ax_p.legend(frameon=False, fontsize=8.5, ncol=2)

    nonlinear = np.any(P_nl != 0)
    measured, expected = (dEdt - P_coll - P_dr, P_nl) if nonlinear else (dEdt, P)
    residual = np.linalg.norm(measured - expected) / np.linalg.norm(expected)
    ax_p.text(0.985, 0.05, f'residual {residual:.2e}'
                           f'  ({"nonlinear channel" if nonlinear else "total budget"})',
              transform=ax_p.transAxes, ha='right', va='bottom', fontsize=8.5, color='#444')

    fig.savefig(outfile)
    plt.close(fig)
    return outfile


def per_kx_residual(netcdf_file, which, time_min, time_max):
    """Budget residual resolved by |kx|, rather than accumulated up to a cutoff.

    A cumulative cutoff hides which mode is responsible: it can only ever show
    the residual getting worse as modes are added.  This isolates each |kx|.
    """
    ncdata = Dataset(netcdf_file)
    weight = _weights(ncdata)
    time = np.array(ncdata.variables['t'][:])
    kx = np.array(ncdata.variables['kx'][:])

    if which == 'phi':
        signal = _field_line_average(ncdata, 'RH_phi_I', weight)
        inertia = _field_line_average(ncdata, 'RH_inertia', weight)
        fluxes = {n: _field_line_average(ncdata, n, weight) for n in
                  ('RH_fluxes_phi_even', 'RH_fluxes_phi_odd',
                   'RH_fluxes_collisional', 'RH_fluxes_drift_trapped',
                   'RH_fluxes_drift_passing')}
        fluxes = {k: v for k, v in fluxes.items() if v is not None}
        nonlinear = sum(v for k, v in fluxes.items() if 'phi_' in k or 'apar' in k or 'bpar' in k)
        other = sum(v for k, v in fluxes.items() if 'coll' in k or 'drift' in k)
        bmag = np.array(ncdata.variables['bmag'][:])[:, 0]
        gds22 = np.array(ncdata.variables['gds22'][:])[:, 0]
        shat = float(np.array(ncdata.variables['shat'][...]))
        charge = np.array(ncdata.variables['charge'][:]); mass = np.array(ncdata.variables['mass'][:])
        temp = np.array(ncdata.variables['temp'][:]); dens = np.array(ncdata.variables['dens'][:])
        b_ref = (kx[:, None] / bmag[None, :])**2 * (gds22 / shat**2)[None, :]
        norm = np.zeros_like(kx)
        for z_s, m_s, T_s, n_s in zip(charge, mass, temp, dens):
            b_s = b_ref * (m_s * T_s / z_s**2)
            G0 = (weight[None, :] * np.i0(b_s / 2) * np.exp(-b_s / 2)).sum(axis=1)
            norm += z_s**2 * n_s / T_s * (1 - G0)
        inertia2 = np.abs(inertia)**2
        E = np.abs(signal)**2 / (2 * inertia2)[None, :] * norm[None, :]
        pref = norm[None, :] / inertia2[None, :]
        P_nl = -np.real(1j * kx[None, :] * nonlinear * np.conj(signal)) * pref
        P_ot = -np.real(1j * kx[None, :] * other * np.conj(signal)) * pref
    else:
        signal = _field_line_average_per_species(ncdata, 'RH_upar', weight).sum(axis=1)
        inertia = _field_line_average_per_species(ncdata, 'RH_upar_inertia', weight)
        mass = np.array(ncdata.variables['mass'][:]); dens = np.array(ncdata.variables['dens'][:])
        u = _field_line_average_per_species(ncdata, 'RH_upar', weight)
        ws = (mass * dens)[None, :, None]
        i2 = np.abs(inertia)[None, :, :]**2
        E = (0.5 * ws * np.abs(u)**2 / i2).sum(axis=1)
        def pw(name):
            fl = _field_line_average_per_species(ncdata, name, weight)
            if fl is None: return np.zeros(E.shape)
            return (-ws * np.real(1j * kx[None, None, :] * fl * np.conj(u)) / i2).sum(axis=1)
        P_nl = pw('RH_upar_flux_nonlinear')
        P_ot = pw('RH_upar_flux_collisional') + pw('RH_upar_flux_drift')

    dEdt = np.gradient(E, time, axis=0)
    interior = slice(1, -1)
    time, E, dEdt, P_nl, P_ot = (a[interior] for a in (time, E, dEdt, P_nl, P_ot))
    window = (time >= time_min) & (time <= time_max)

    magnitudes = sorted({round(abs(k), 8) for k in kx if abs(k) > 1e-12})
    out = []
    for magnitude in magnitudes:
        columns = [i for i, k in enumerate(kx) if abs(abs(k) - magnitude) < 1e-8]
        measured = (dEdt[:, columns].sum(axis=1) - P_ot[:, columns].sum(axis=1))[window]
        expected = P_nl[:, columns].sum(axis=1)[window]
        if np.linalg.norm(expected) == 0:
            continue
        out.append((magnitude, np.linalg.norm(measured - expected) / np.linalg.norm(expected),
                    np.abs(expected).max()))
    return out


def figure_kx(netcdf_file, outfile, time_min, time_max, title=''):
    """Budget residual against |kx| for both invariants, mode by mode.

    The construction is not a long-wavelength one -- nothing in the annihilation
    it rests on expands in kx rho_i -- so a residual that grows with kx is a
    numerical statement about the implementation, not a property of the theory.
    """
    fig, (ax, ax_w) = plt.subplots(2, 1, figsize=(6.2, 4.4), sharex=True,
                                   gridspec_kw=dict(height_ratios=[2, 1], hspace=0.1))
    for which, colour, marker, label in (('phi', PHI, 'o', r'$\varphi_{\rm RH}$'),
                                         ('upar', UPA, 's', r'$u_{\parallel\rm RH}$')):
        rows = per_kx_residual(netcdf_file, which, time_min, time_max)
        if not rows:
            continue
        k = [r[0] for r in rows]; res = [r[1] for r in rows]; amp = [r[2] for r in rows]
        ax.plot(k, res, marker + '-', color=colour, lw=1.5, ms=5, label=label)
        ax_w.semilogy(k, amp, marker + '-', color=colour, lw=1.2, ms=4, alpha=0.8)
    ax.set_yscale('log')
    ax.set_ylabel('budget residual, per mode')
    ax.legend(frameon=False, fontsize=9)
    if title:
        ax.set_title(title, fontsize=9.5, loc='left')
    ax_w.set_ylabel('drive  $|P_{\\rm nl}|$')
    ax_w.set_xlabel(r'$k_x \rho_i$')
    fig.savefig(outfile)
    plt.close(fig)
    return outfile


def figure_summary(cases, outfile):
    """Residual for every benchmark case, both invariants side by side.

    `cases` is a list of (label, phi_residual, upar_residual); either may be
    None where that case was not measured.
    """
    labels = [c[0] for c in cases]
    y = np.arange(len(cases))
    fig, ax = plt.subplots(figsize=(7.0, 0.30 * len(cases) + 1.4))
    height = 0.38
    for offset, index, colour, label in ((+height/2, 1, PHI, r'$\varphi_{\rm RH}$'),
                                         (-height/2, 2, UPA, r'$u_{\parallel\rm RH}$')):
        values = [c[index] if c[index] is not None else np.nan for c in cases]
        ax.barh(y + offset, values, height=height, color=colour, label=label, alpha=0.9)
    ax.axvline(0.08, color=GREY, ls='--', lw=1.0)
    ax.text(0.084, len(cases) - 0.4, 'benchmark tolerance', fontsize=7.5, color=GREY, va='top')
    ax.set_xscale('log')
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel('relative budget residual')
    ax.legend(frameon=False, fontsize=9, loc='lower right')
    ax.grid(axis='y', alpha=0)
    fig.savefig(outfile)
    plt.close(fig)
    return outfile


def figure_quadrature(rows, outfile):
    """Conservation error against parallel resolution, at several kx.

    The construction is exact at every wavelength, so each curve must fall with
    resolution.  A curve that flattens -- or rises -- is a defect in the
    quadrature rather than a limit of the theory.

    `rows` maps kx to (nzed list, phi errors, upar errors).
    """
    kxs = sorted(rows)
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.1), sharey=True)
    shades = plt.cm.viridis(np.linspace(0.15, 0.85, len(kxs)))
    for ax, index, name, colour in ((axes[0], 1, r'$\varphi_{\rm RH}$', PHI),
                                    (axes[1], 2, r'$u_{\parallel\rm RH}$', UPA)):
        for kx, shade in zip(kxs, shades):
            nzed, *series = rows[kx]
            ax.loglog(nzed, series[index - 1], 'o-', color=shade, lw=1.4, ms=4,
                      label=rf'$k_x\rho_i={kx}$')
        ax.set_xlabel(r'$n_{\rm zed}$')
        ax.set_title(name, fontsize=10, loc='left', color=colour)
    axes[0].set_ylabel('conservation error')
    axes[1].legend(frameon=False, fontsize=7.5, loc='upper right', ncol=1)
    fig.tight_layout()
    fig.savefig(outfile)
    plt.close(fig)
    return outfile


def figure_drift_scan(xdrift, phi_err, upar_err, trapped, phi_tf, upar_tf, outfile):
    """What actually breaks the annihilation.

    Left: conservation error against the strength of the magnetic drift, at
    fixed kx.  It vanishes with the drift, which is the signature of an
    incomplete cancellation between the discrete streaming and drift operators
    rather than of a quadrature error.  Right: the same error against trapped
    fraction, where the mirror force -- the other half of discrete streaming --
    is strongest.
    """
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(7.6, 3.0))
    ax.loglog(xdrift, phi_err, 'o-', color=PHI, lw=1.6, ms=5, label=r'$\varphi_{\rm RH}$')
    ax.loglog(xdrift, upar_err, 's-', color=UPA, lw=1.6, ms=5, label=r'$u_{\parallel\rm RH}$')
    ax.set_xlabel('xdriftknob  (magnetic drift strength)')
    ax.set_ylabel('conservation error')
    ax.set_title(r'At fixed $k_x\rho_i=2$', fontsize=9.5, loc='left')
    ax.legend(frameon=False, fontsize=8.5)

    ax2.loglog(trapped, phi_tf, 'o-', color=PHI, lw=1.6, ms=5)
    ax2.loglog(trapped, upar_tf, 's-', color=UPA, lw=1.6, ms=5)
    ax2.set_xlabel('trapped fraction')
    ax2.set_title(r'At fixed $k_x\rho_i=2$, drift on', fontsize=9.5, loc='left')
    fig.tight_layout()
    fig.savefig(outfile)
    plt.close(fig)
    return outfile


def figure_asymptotics(datasets, outfile):
    """Inertia against the drift-orbit phase, for two aspect ratios.

    `datasets` maps a label to (Q, enhancement, I_upar, eps, q).  Plotting
    against Q rather than kx is the point: the neoclassical enhancement leaves
    its plateau when Q ~ 1, at values of kx rho that differ by q/eps between the
    two, so the controlling parameter is the orbit width against 1/kx and not
    the gyroradius against 1/kx.
    """
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(7.8, 3.1))
    colours = [PHI, UPA]
    for (label, (Q, enh, Iu, eps, q)), colour in zip(datasets.items(), colours):
        ax.loglog(Q, np.array(enh) - 1, 'o-', color=colour, lw=1.5, ms=4, label=label)
        plateau = 1.6 * q**2 / np.sqrt(eps)
        ax.axhline(plateau, color=colour, ls=':', lw=1.1)
        ax2.semilogx(Q, np.array(Iu), 's-', color=colour, lw=1.5, ms=4, label=label)
    ax.axvline(1.0, color=GREY, ls='--', lw=1.0)
    ax.text(1.15, ax.get_ylim()[0] * 2, r'$Q = 1$', fontsize=8, color=GREY)
    ax.set_xlabel(r'drift-orbit phase  $Q \sim k_x \rho\, q/\epsilon$')
    ax.set_ylabel(r'$\langle I\rangle/[(1-\Gamma_0)Z^2n/T] - 1$')
    ax.set_title('Neoclassical enhancement', fontsize=9.5, loc='left')
    ax.legend(frameon=False, fontsize=8)
    ax2.axhline(0.5, color=GREY, ls=':', lw=1.1)
    ax2.text(Q[0], 0.52, r'$1/2$', fontsize=8, color=GREY)
    ax2.axvline(1.0, color=GREY, ls='--', lw=1.0)
    ax2.set_xlabel(r'drift-orbit phase  $Q$')
    ax2.set_ylabel(r'$\langle I_u\rangle$')
    ax2.set_title('Parallel-flow inertia', fontsize=9.5, loc='left')
    fig.tight_layout()
    fig.savefig(outfile)
    plt.close(fig)
    return outfile


def figure_phase_mixing(Q, W_all, W_passing, outfile):
    """Coherence of the transit-average weight against the drift-orbit phase.

    Every flux channel carries exactly one factor of this weight, so the law it
    follows -- unity while Q < 1, then 1/Q -- sets the short-wavelength
    behaviour of all of them at once.
    """
    Q = np.asarray(Q)
    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    ax.loglog(Q, W_all, 'o-', color=PHI, lw=1.6, ms=5, label=r'all particles')
    ax.loglog(Q, W_passing, 's--', color=UPA, lw=1.4, ms=4, label=r'passing only')
    tail = Q[Q > 3]
    ax.loglog(tail, 4.0 / tail, ':', color=GREY, lw=1.4, label=r'$4/Q$')
    ax.axhline(1.0, color=GREY, ls=':', lw=1.0)
    ax.axvline(1.0, color=GREY, ls='--', lw=1.0)
    ax.set_xlabel(r'drift-orbit phase  $Q$')
    ax.set_ylabel(r'$\langle |W| \rangle$')
    ax.set_title('Phase mixing of the projection weight', fontsize=9.5, loc='left')
    ax.legend(frameon=False, fontsize=8.5)
    fig.tight_layout()
    fig.savefig(outfile)
    plt.close(fig)
    return outfile
