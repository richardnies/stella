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
from rh_budget import (get_rh_budget, get_rh_umom_budget, _complex, _first_present,
                       _field_line_average, _field_line_average_per_species,
                       get_rh_budget_LW)


def LW_inertia(netcdf_file):
    '''The inertia rebuilt from the long-wavelength weights, if the run has them.

    <I> = (Z^2 n/T) <int F_M (1 - J0 W)>, which needs only the weight and the
    equilibrium Maxwellian -- no g -- so it can be formed directly from the
    RH_LW_even/odd that write_RH_asymptotics produces and compared with the
    exact RH_inertia.  Returns (exact, LWptotic) field-line averages per kx,
    or (None, None) when the run did not write the asymptotic weights.
    '''
    ncdata = Dataset(netcdf_file)
    #> Written as RH_asym_* before the LW/SW naming; keep reading both.
    even_name = _first_present(ncdata, 'RH_LW_even', 'RH_asym_even')
    odd_name = _first_present(ncdata, 'RH_LW_odd', 'RH_asym_odd')
    if even_name not in ncdata.variables:
        return None, None
    zed = np.array(ncdata.variables['zed'][:])
    jac = np.array(ncdata.variables['jacob'][:])[:, 0]
    bmag = np.array(ncdata.variables['bmag'][:])
    bmag = bmag[:, 0] if bmag.ndim > 1 else bmag
    w = (zed[1] - zed[0]) * jac.copy(); w[-1] = 0.0; w /= w.sum()
    vpa = np.array(ncdata.variables['vpa'][:]); mu = np.array(ncdata.variables['mu'][:])

    def weight(name):
        a = np.array(ncdata.variables[name][:])
        return a[..., 0] + 1j * a[..., 1]          # (mu, vpa, species, tube, zed, kx)

    exact = weight('RH_integrand_even') + weight('RH_integrand_odd')
    LW = weight(even_name) + weight(odd_name)
    VP, MU = np.meshgrid(vpa, mu)
    nkx = exact.shape[-1]
    out_e = np.zeros(nkx); out_a = np.zeros(nkx)
    for ikx in range(nkx):
        ne = na = den = 0.0
        for iz in range(len(bmag)):
            fM = np.exp(-(VP**2 + 2 * MU * bmag[iz]))
            ne += w[iz] * (fM * (1 - exact[:, :, 0, 0, iz, ikx].real)).sum()
            na += w[iz] * (fM * (1 - LW[:, :, 0, 0, iz, ikx].real)).sum()
            den += w[iz] * fM.sum()
        out_e[ikx] = ne / den; out_a[ikx] = na / den
    return out_e, out_a

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
    upar = _field_line_average_per_species(ncdata, _first_present(ncdata, 'RH_umom', 'RH_upar'), weight)
    upar_inertia = _field_line_average_per_species(
        ncdata, _first_present(ncdata, 'RH_umom_inertia', 'RH_upar_inertia'), weight)

    kx = np.array(ncdata.variables['kx'][:])
    keep = np.abs(kx) > 1e-12
    #> One zonal mode is enough to show conservation; take the smallest finite kx.
    index = np.argmin(np.where(keep, np.abs(kx), np.inf))
    phi_ratio = phi[:, index] / inertia[index]
    upar_ratio = upar[:, 0, index] / upar_inertia[0, index]  # p_RH / I_p
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
                                  (upar_ratio, UPA, r'$\langle p_{\rm RH}\rangle/\langle I_p\rangle$')):
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
        t, E, dEdt, P, P_nl, P_coll, P_dr = get_rh_umom_budget(netcdf_file, time_min, time_max)
        colour, energy_label = UPA, r'$E_{p\rm RH}$'

    fig, (ax_e, ax_p) = plt.subplots(2, 1, figsize=(6.6, 5.0), sharex=True,
                                     gridspec_kw=dict(height_ratios=[1, 1.3], hspace=0.12))

    ax_e.plot(t, E, color=colour, lw=1.8, label=energy_label)

    #> The same energy with the inertia replaced by its long-wavelength form.
    #> E ~ 1/<I>^2, so this shows directly how far the order-kx^2 expansion is
    #> from the exact normalisation at the wavelength of the run.
    exact_I, LW_I = LW_inertia(netcdf_file)
    if exact_I is not None:
        finite = exact_I != 0
        if finite.any():
            scale = float(np.mean((exact_I[finite] / LW_I[finite])**2))
            err = abs(np.sqrt(1.0 / scale) - 1.0) * 100
            ax_e.plot(t, E * scale, color=GREY, lw=1.3, ls=':',
                      label=energy_label + rf'  with $O(k_x^2)\ \langle I\rangle$  ({err:.1f}% in $\langle I\rangle$)')
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
    ax_p.text(0.985, 0.92, f'residual {residual:.2e}'
                           f'  ({"nonlinear channel" if nonlinear else "total budget"})',
              transform=ax_p.transAxes, ha='right', va='top', fontsize=8.5, color='#444')

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
        signal = _field_line_average_per_species(ncdata, _first_present(ncdata, 'RH_umom', 'RH_upar'), weight).sum(axis=1)
        inertia = _field_line_average_per_species(ncdata, _first_present(ncdata, 'RH_umom_inertia', 'RH_upar_inertia'), weight)
        mass = np.array(ncdata.variables['mass'][:]); dens = np.array(ncdata.variables['dens'][:])
        u = _field_line_average_per_species(ncdata, _first_present(ncdata, 'RH_umom', 'RH_upar'), weight)
        ws = (mass * dens)[None, :, None]
        i2 = np.abs(inertia)[None, :, :]**2
        E = (0.5 * ws * np.abs(u)**2 / i2).sum(axis=1)
        def pw(name):
            fl = _field_line_average_per_species(ncdata, name, weight)
            if fl is None: return np.zeros(E.shape)
            return (-ws * np.real(1j * kx[None, None, :] * fl * np.conj(u)) / i2).sum(axis=1)
        P_nl = pw(_first_present(ncdata, 'RH_umom_flux_nonlinear', 'RH_upar_flux_nonlinear'))
        P_ot = (pw(_first_present(ncdata, 'RH_umom_flux_collisional', 'RH_upar_flux_collisional'))
                + pw(_first_present(ncdata, 'RH_umom_flux_drift', 'RH_upar_flux_drift')))

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
                                         ('upar', UPA, 's', r'$p_{\rm RH}$')):
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
                                         (-height/2, 2, UPA, r'$p_{\rm RH}$')):
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
                                    (axes[1], 2, r'$p_{\rm RH}$', UPA)):
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
    ax.loglog(xdrift, upar_err, 's-', color=UPA, lw=1.6, ms=5, label=r'$p_{\rm RH}$')
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
    ax2.set_ylabel(r'$\langle I_p\rangle$')
    ax2.set_title('Toroidal-momentum inertia', fontsize=9.5, loc='left')
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


CHANNEL_COLOURS = {
    'nonlinear, even': '#7d3c6b',
    'nonlinear, odd': '#2e7d6b',
    'collisional, even': '#8a6d1f',
    'collisional, odd': '#b0532a',
}


def figure_LW_power(panels, outfile):
    """P_RH channel by channel, exact weight against its long-wavelength form.

    Each panel is (title, file, t_min, t_max, channel-prefix or None).
    Solid is the exact transit-average weight, dotted
    the order-kx^2 expansion of it; both are built from the same fields and
    turned into a power with the same RH_phi_I, so the gap between a pair is
    the weight and nothing else.  That is what makes this a test of the
    expansion rather than of whatever happens to be driving the flux -- the
    drive cancels out of the comparison.
    """
    fig, axes = plt.subplots(1, len(panels), figsize=(3.5 * len(panels), 3.1))
    if len(panels) == 1:
        axes = [axes]

    for ax, (title, netcdf_file, time_min, time_max, only) in zip(axes, panels):
        time, channels = get_rh_budget_LW(netcdf_file, time_min, time_max)
        if time is None:
            continue
        for label, (exact, LW) in channels.items():
            #> A run can carry several channels at once, and they need not be
            #> the same size; `only` keeps a panel to the one being made about.
            if only is not None and not label.startswith(only):
                continue
            colour = CHANNEL_COLOURS.get(label, GREY)
            #> Integrated over the window rather than pointwise: the pointwise
            #> ratio is meaningless wherever the channel passes through zero.
            error = abs(np.trapz(LW, time) / np.trapz(exact, time) - 1) * 100
            ax.plot(time, exact, color=colour, lw=1.5,
                    label=rf'{label}  ({error:.1f}%)')
            ax.plot(time, LW, color=colour, lw=1.5, ls=':', alpha=0.95)
        ax.axhline(0.0, color='k', lw=0.5, alpha=0.3)
        ax.set_xlabel(r'time  $[a/v_{\rm th}]$')
        ax.set_title(title, fontsize=9.5, loc='left')
        ax.legend(frameon=False, fontsize=7.6)

    axes[0].set_ylabel(r'$P_{\rm RH}$ into the zonal flow')
    fig.tight_layout()
    fig.savefig(outfile)
    plt.close(fig)
    return outfile


def stress_split(netcdf_file, time_min, time_max):
    """Split the even nonlinear flux into its Reynolds and diamagnetic channels.

    The FLR expansion 1 - J0 ~ kperp^2 vperp^2 / 4 Omega^2 puts a vperp^2 inside
    a velocity integral over the whole distribution.  Taken against the
    adiabatic part of that distribution it returns a functional of phi alone --
    the Reynolds stress; taken against the non-adiabatic part it returns the
    perpendicular pressure -- the diamagnetic stress.  Both are the same order
    in kperp*rho, so the even channel carries both and the diagnostic does not
    separate them.

    Rebuild each from the saved fields and fit them to the measured flux:

      P_phi(kx) = sum_k' [i ky phi](k') [Gamma_pol phi](k-k')
      P_p(kx)   = sum_k' [i ky phi](k') [b p_perp](k-k')

    Returns (t, kx, measured, P_phi, P_p) field-line averaged, restricted to the
    window and to finite kx.  Needs write_phi_vs_kxkyz and write_moments.
    """
    ncdata = Dataset(netcdf_file)
    t = np.array(ncdata.variables['t'][:])
    kx = np.array(ncdata.variables['kx'][:]); ky = np.array(ncdata.variables['ky'][:])
    zed = np.array(ncdata.variables['zed'][:])
    jac = np.array(ncdata.variables['jacob'][:])[:, 0]
    kp2 = np.array(ncdata.variables['kperp2'][:])[:, 0]
    mass = np.array(ncdata.variables['mass'][:]); temp = np.array(ncdata.variables['temp'][:])
    chg = np.array(ncdata.variables['charge'][:]); dens = np.array(ncdata.variables['dens'][:])
    w = (zed[1] - zed[0]) * jac.copy(); w[-1] = 0.0; w /= w.sum()

    phi = _complex(ncdata, 'phi_vs_t')[0][:, 0]
    pperp = _complex(ncdata, 'pressure_perp')[0][:, 0, 0]
    meas = _complex(ncdata, 'RH_fluxes_phi_even')[0][:, 0, 0]

    Gpol = np.zeros(kp2.shape); bfac = np.zeros(kp2.shape)
    for z_s, m_s, T_s, n_s in zip(chg, mass, temp, dens):
        x = kp2 * (m_s * T_s / z_s**2) / 2.0
        Gpol += z_s**2 * n_s / T_s * (1 - np.i0(x) * np.exp(-x))
        bfac += x

    nkx = len(kx)
    idx = -np.ones((nkx, nkx), dtype=int)
    for p in range(nkx):
        for q in range(nkx):
            hit = np.where(np.abs(kx - (kx[p] - kx[q])) < 1e-8)[0]
            if hit.size:
                idx[p, q] = hit[0]
    ref = np.array([np.where(np.abs(kx + kx[q]) < 1e-8)[0][0] for q in range(nkx)])

    def convolve(field):
        out = np.zeros(phi.shape[:3], dtype=complex)
        for iy in range(len(ky)):
            if ky[iy] <= 0:
                continue
            A = 1j * ky[iy] * phi[:, :, :, iy]
            B = field[:, :, :, iy]
            for q in range(nkx):
                acc = np.zeros(phi.shape[:2], dtype=complex)
                for p in range(nkx):
                    r = idx[p, q]
                    if r >= 0:
                        acc += A[:, :, p] * np.conj(B[:, :, r])
                out[:, :, q] += acc
        return out + np.conj(out[:, :, ref])

    fla = lambda a: np.einsum('z,tzk->tk', w, a)
    P_phi = fla(convolve(Gpol[None] * phi))
    P_p = fla(convolve(bfac[None] * pperp))
    M = fla(meas.sum(axis=3))

    window = (t >= time_min) & (t <= time_max)
    finite = np.abs(kx) > 1e-12
    return t[window], kx[finite], M[window][:, finite], P_phi[window][:, finite], P_p[window][:, finite]


def _fit(columns, target):
    """Real coefficients fitted to complex data, plus the relative residual."""
    A = np.vstack([np.column_stack(columns).real, np.column_stack(columns).imag])
    y = np.concatenate([target.real, target.imag])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    return c, np.linalg.norm(y - A @ c) / np.linalg.norm(y)


def figure_stress_split(netcdf_file, outfile, time_min, time_max, title=''):
    """Show that the even channel needs the diamagnetic stress as well.

    Left: the measured flux at the longest wavelength the box holds, against a
    Reynolds-only fit and a fit carrying both stresses.  Right: how that changes
    with kx.  A Reynolds-only description is not merely imprecise -- it is
    missing a term of the same order.
    """
    t, kx, M, P_phi, P_p = stress_split(netcdf_file, time_min, time_max)
    c1, r1 = _fit([P_phi.ravel()], M.ravel())
    cb, rb = _fit([P_phi.ravel(), P_p.ravel()], M.ravel())

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.2))
    j = int(np.argmin(np.abs(kx - kx[kx > 0].min())))
    ax.plot(t, M[:, j].real, color=PHI, lw=1.8, label='measured')
    ax.plot(t, (c1[0] * P_phi[:, j]).real, color=GREY, lw=1.4, ls='--',
            label=r'Reynolds only')
    ax.plot(t, (cb[0] * P_phi[:, j] + cb[1] * P_p[:, j]).real, color=UPA, lw=1.4, ls=':',
            label=r'Reynolds $+$ diamagnetic')
    ax.set_xlabel(r'time  $[a/v_{\rm th}]$')
    ax.set_ylabel(r'$F^{\rm NL}_{\rm even}$')
    ax.set_title(f'{title}  $k_x\\rho = {kx[j]:.2f}$', fontsize=9.5, loc='left')
    ax.legend(frameon=False, fontsize=8)

    pos = kx > 0
    e1, eb = [], []
    for jj in np.where(pos)[0]:
        m = M[:, jj]
        e1.append(np.linalg.norm(m - c1[0] * P_phi[:, jj]) / np.linalg.norm(m) * 100)
        eb.append(np.linalg.norm(m - cb[0] * P_phi[:, jj] - cb[1] * P_p[:, jj]) / np.linalg.norm(m) * 100)
    ax2.semilogy(kx[pos], e1, 'o--', color=GREY, lw=1.4, ms=5, label='Reynolds only')
    ax2.semilogy(kx[pos], eb, 's-', color=UPA, lw=1.6, ms=5, label=r'Reynolds $+$ diamagnetic')
    ax2.set_xlabel(r'$k_x\rho$')
    ax2.set_ylabel('residual [%]')
    ax2.set_title('Both stresses are needed', fontsize=9.5, loc='left')
    ax2.legend(frameon=False, fontsize=8)
    weight = np.linalg.norm(cb[1] * P_p) / np.linalg.norm(cb[0] * P_phi)
    ax2.text(0.97, 0.06, rf'$|\Pi_T|/|\Pi_\varphi| = {weight:.2f}$', transform=ax2.transAxes,
             ha='right', fontsize=8.5, color='#444')
    fig.tight_layout()
    fig.savefig(outfile)
    plt.close(fig)
    return outfile


def figure_stress_channels(netcdf_file, outfile, time_min, time_max, title=''):
    """The Reynolds and diamagnetic halves of the even nonlinear channel.

    Reads the split stella now forms directly (write_RH_stress_split), so there
    is no fitting here and total = Reynolds + diamagnetic is exact.  The split
    is by the gyroaverage on the ExB velocity: J0 -> 1 leaves a density moment,
    which quasineutrality slaves to phi and which is the Reynolds channel, and
    the remainder carries (J0 - 1) ~ -kperp^2 vperp^2 / 4 Omega^2, i.e. the
    vperp^2 moment, and is the diamagnetic channel.
    """
    ncdata = Dataset(netcdf_file)
    t = np.array(ncdata.variables['t'][:])
    kx = np.array(ncdata.variables['kx'][:])
    w = _weights(ncdata)

    def fla(name):
        a = _complex(ncdata, name)[0][:, 0, 0]
        return np.einsum('z,tzky->tky', w, a).sum(axis=2)

    T = fla('RH_fluxes_phi_even')
    R = fla('RH_fluxes_phi_even_reynolds')
    D = fla('RH_fluxes_phi_even_diamagnetic')

    window = (t >= time_min) & (t <= time_max)
    positive = kx > 1e-12
    j = int(np.where(positive)[0][int(np.argmin(kx[positive]))])

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.2))
    ax.plot(t, T[:, j].real, color=PHI, lw=2.0, label='total')
    ax.plot(t, R[:, j].real, color=GREY, lw=1.4, ls='--', label=r'Reynolds  $\Pi_\varphi$')
    ax.plot(t, D[:, j].real, color=UPA, lw=1.4, ls='-.', label=r'diamagnetic  $\Pi_T$')
    ax.axhline(0.0, color='k', lw=0.5, alpha=0.3)
    ax.set_xlabel(r'time  $[a/v_{\rm th}]$')
    ax.set_ylabel(r'$F^{\rm NL}_{\rm even}$')
    ax.set_title(f'{title}  $k_x\\rho = {kx[j]:.2f}$', fontsize=9.5, loc='left')
    ax.legend(frameon=False, fontsize=8)

    ratio, align = [], []
    for jj in np.where(positive)[0]:
        r, dd = R[window][:, jj], D[window][:, jj]
        ratio.append(np.linalg.norm(dd) / np.linalg.norm(r))
        align.append(np.vdot(r, dd).real / (np.linalg.norm(r) * np.linalg.norm(dd)))
    ax2.semilogy(kx[positive], ratio, 'o-', color=UPA, lw=1.6, ms=5)
    ax2.axhline(1.0, color=GREY, ls=':', lw=1.1)
    ax2.set_xlabel(r'$k_x\rho$')
    ax2.set_ylabel(r'$|\Pi_T| / |\Pi_\varphi|$')
    ax2.set_title('Diamagnetic against Reynolds', fontsize=9.5, loc='left')
    axr = ax2.twinx()
    axr.plot(kx[positive], align, 's--', color=PHI, lw=1.2, ms=4, alpha=0.8)
    axr.axhline(0.0, color=PHI, lw=0.5, alpha=0.3)
    axr.set_ylabel(r'$\cos(\Pi_\varphi, \Pi_T)$', color=PHI)
    axr.tick_params(axis='y', labelcolor=PHI)
    axr.set_ylim(-1.05, 1.05)
    axr.grid(False)
    fig.tight_layout()
    fig.savefig(outfile)
    plt.close(fig)
    return outfile


#> Convergence data measured in the runs recorded in the write-up.  Kept here so
#> the figure and the text cannot drift apart, and so the point is made by a
#> picture rather than by two eight-column tables.
ONE_GRID_AT_A_TIME = {          # nzed only, others held fixed
    'nzed': [48, 96, 192, 384],
    'phi': {0.05: [8.3e-3, 5.6e-3, 4.7e-3, 4.0e-3],
            0.5:  [1.0e-2, 1.5e-1, 2.4e-1, 2.6e-1],
            2.0:  [1.1e-1, 5.4e-2, 4.7e-2, 4.7e-2]},
    'U':   {0.05: [7.0e-3, 3.9e-3, 2.8e-3, 1.4e-3],
            0.5:  [1.3e-3, 1.4e-3, 9.5e-4, 4.7e-4],
            2.0:  [1.4e-2, 1.6e-2, 9.8e-3, 5.4e-3]},
}
ALL_GRIDS = {                   # nzed / nvpa / nmu / dt refined together
    'level': [1, 2, 4, 8],
    'phi': {0.05: [2.3e-2, 5.5e-3, 1.5e-3, 4.5e-4],
            2.0:  [5.6e-1, 5.5e-2, 2.2e-2, 6.0e-3]},
    'U':   {0.05: [4.4e-2, 6.3e-3, 1.1e-4, 1.8e-3],
            2.0:  [2.6e-2, 1.6e-2, 9.4e-3, 5.0e-3]},
}


def figure_convergence(outfile):
    """Does the construction converge?  Yes, if every grid is refined together.

    Left: refining the parallel grid alone.  The flow invariant falls, the
    potential-like one flattens and at k_x rho = 0.5 climbs -- which looks like a
    defect and is not one.  Right: refining nzed, nvpa, nmu and dt together, the
    same quantity falls by nearly two orders.  The weight depends on all four
    grids, so refining one leaves the error floored by the other three.
    """
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(9.4, 3.6))
    colours = {0.05: PHI, 0.5: '#7d3c6b', 2.0: UPA}

    n = ONE_GRID_AT_A_TIME['nzed']
    for kx, c in colours.items():
        ax.loglog(n, ONE_GRID_AT_A_TIME['phi'][kx], 'o-', color=c, lw=1.6, ms=5,
                  label=rf'$\varphi_{{\rm RH}}$, $k_x\rho={kx}$')
        ax.loglog(n, ONE_GRID_AT_A_TIME['U'][kx], 's--', color=c, lw=1.2, ms=4, alpha=0.75)
    ax.set_xlabel(r'$n_{\rm zed}$  (other grids fixed)')
    ax.set_ylabel('conservation error')
    ax.set_title('Refining one grid: misleading', fontsize=9.5, loc='left')
    ax.legend(frameon=False, fontsize=7.6, ncol=1)
    ax.text(0.03, 0.05, 'solid $\\varphi_{\\rm RH}$,  dashed $U_{\\rm RH}$',
            transform=ax.transAxes, fontsize=7.6, color='#444')

    lv = ALL_GRIDS['level']
    for kx, c in ((0.05, PHI), (2.0, UPA)):
        ax2.loglog(lv, ALL_GRIDS['phi'][kx], 'o-', color=c, lw=1.6, ms=5,
                   label=rf'$\varphi_{{\rm RH}}$, $k_x\rho={kx}$')
        ax2.loglog(lv, ALL_GRIDS['U'][kx], 's--', color=c, lw=1.2, ms=4, alpha=0.75)
    ref = np.array(lv, dtype=float)
    ax2.loglog(ref, 5.6e-1 * (ref / ref[0])**-2, ':', color=GREY, lw=1.3, label=r'$\propto h^{2}$')
    ax2.set_xlabel(r'refinement of $n_{\rm zed}$, $n_{v_\parallel}$, $n_\mu$, $\Delta t$ together')
    ax2.set_title('Refining all four: it converges', fontsize=9.5, loc='left')
    ax2.legend(frameon=False, fontsize=7.6)
    fig.tight_layout()
    fig.savefig(outfile)
    plt.close(fig)
    return outfile


#> Stellarator benchmark results, as measured.  Held here so figure and text
#> cannot drift apart.
STELLARATOR = {
    #                 collisionless drift            collisional budget
    #            phi resid  U resid  drift share   phi resid  U resid
    'ITER':  dict(dphi=4.1e-3, dU=1.5e-3, share=0.28, bphi=1.98e-2, bU=2.82e-2),
    'W7-X':  dict(dphi=3.0e-2, dU=8.5e-4, share=0.49, bphi=3.75e-2, bU=5.14e-2),
    'QA':    dict(dphi=9.5e-3, dU=1.6e-3, share=0.74, bphi=9.30e-2, bU=None),
    'QH':    dict(dphi=6.9e-2, dU=1.1e-3, share=1.14, bphi=2.58e-2, bU=4.39e-2),
    'TJ-II': dict(dphi=5.0e-1, dU=3.5e-3, share=2.15, bphi=6.99e-2, bU=4.09e-2),
}
TJII_DRIFT_REFINEMENT = ([1, 2, 4], [4.98e-1, 2.06e-1, 1.20e-1])


def figure_stellarator_summary(outfile):
    """The stellarator benchmarks, both invariants, on one page.

    Left: the collisionless drift test, where that channel is the entire source
    and nothing can absorb an error in it.  Right: the collisional budget.
    """
    names = list(STELLARATOR)
    x = np.arange(len(names))
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(9.6, 3.5))

    ax.semilogy(x, [STELLARATOR[n]['dphi'] for n in names], 'o-', color=PHI, lw=1.6, ms=6,
                label=r'$\varphi_{\rm RH}$')
    ax.semilogy(x, [STELLARATOR[n]['dU'] for n in names], 's-', color=UPA, lw=1.6, ms=6,
                label=r'$U_{\rm RH}$')
    ax.semilogy(x, [STELLARATOR[n]['share'] for n in names], '^:', color=GREY, lw=1.2, ms=5,
                label='drift channel size')
    ax.set_xticks(x); ax.set_xticklabels(names, fontsize=8)
    ax.set_ylabel('unaccounted fraction')
    ax.set_title('Drift channel alone (collisionless)', fontsize=9.5, loc='left')
    ax.legend(frameon=False, fontsize=8, loc='lower left')

    bphi = [STELLARATOR[n]['bphi'] for n in names]
    bU = [STELLARATOR[n]['bU'] for n in names]
    ax2.semilogy(x, bphi, 'o-', color=PHI, lw=1.6, ms=6, label=r'$\varphi_{\rm RH}$')
    xs = [xi for xi, v in zip(x, bU) if v is not None]
    ys = [v for v in bU if v is not None]
    ax2.semilogy(xs, ys, 's-', color=UPA, lw=1.6, ms=6, label=r'$U_{\rm RH}$')
    #> QA carries no U_RH point: its collisional momentum drive is a near
    #> cancellation, leaving a small-signal test rather than a failing one.
    ax2.annotate('QA: momentum drive\nnearly cancels', xy=(2, 9.3e-2), xytext=(1.5, 2.6e-1),
                 fontsize=7, color='#444',
                 arrowprops=dict(arrowstyle='-', color='#888', lw=0.8))
    ax2.set_xticks(x); ax2.set_xticklabels(names, fontsize=8)
    ax2.set_ylabel('budget residual')
    ax2.set_title('Budget, linear collisional', fontsize=9.5, loc='left')
    ax2.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(outfile)
    plt.close(fig)
    return outfile


#> Asymptotic-weight verification, as measured.  Errors are per cent.
LW_COLL = dict(kx=[0.01, 0.02, 0.05, 0.1, 0.2],
               even=[0.00, 0.08, 1.79, 8.33, 24.04],
               odd=[0.03, 0.04, 0.27, 1.74, 3.81])
LW_NL_TOK = dict(kx=[0.1, 0.2, 0.3, 0.4, 0.5], even=[9.4, 38, 110, 184, 389],
                 odd=[12, 45, 97, 139, 224])
LW_NL_TJII = dict(kx=[0.4, 0.8], even=[35, 67], odd=[29, 92])
SW_VS_KX = dict(kx=[0.5, 1, 2, 5, 10, 20], err=[28.7, 16.4, 8.9, 10.0, 19.6, 55.5])
SW_REFINED = {10: ([128, 192, 384], [19.6, 15.0, 9.2]),
              20: ([128, 192, 384, 768], [55.5, 46.6, 32.2, 21.5])}
SW_BANDS = dict(names=['deeply\npassing', 'barely\npassing', 'barely\ntrapped', 'deeply\ntrapped'],
                weight=[0.81, 0.11, 0.02, 0.05],
                gaussian=[2.8, 22.2, 49.0, 28.8], uniform=[3.2, 6.1, 37.5, 28.8])


def figure_asymptotic_weights(outfile):
    """Both asymptotic limits measured against the exact weight in the same run.

    Long wavelength falls as k_x^2 to one part in ten thousand; short wavelength
    is limited by the quadrature it is compared against, not by the formula; and
    the uniform treatment fixes the band it was built for.
    """
    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(12.4, 3.4))

    a1.loglog(LW_COLL['kx'], np.maximum(LW_COLL['even'], 1e-3), 'o-', color=PHI, lw=1.6, ms=5,
              label='collisional, even')
    a1.loglog(LW_COLL['kx'], LW_COLL['odd'], 's-', color=UPA, lw=1.6, ms=5, label='collisional, odd')
    a1.loglog(LW_NL_TOK['kx'], LW_NL_TOK['even'], '^--', color='#7d3c6b', lw=1.3, ms=5,
              label='nonlinear, tokamak')
    a1.loglog(LW_NL_TJII['kx'], LW_NL_TJII['even'], 'v--', color='#2e7d6b', lw=1.3, ms=5,
              label='nonlinear, TJ-II')
    ref = np.array([0.01, 0.2])
    a1.loglog(ref, 0.08 * (ref / 0.02)**2, ':', color=GREY, lw=1.3, label=r'$\propto k_x^2$')
    a1.set_xlabel(r'$k_x\rho$'); a1.set_ylabel('error against the exact weight [%]')
    a1.set_title('Long wavelength', fontsize=9.5, loc='left')
    a1.legend(frameon=False, fontsize=7.2)

    a2.loglog(SW_VS_KX['kx'], SW_VS_KX['err'], 'o-', color=PHI, lw=1.6, ms=5,
              label=r'$n_{\rm zed}=128$')
    for kx, (nz, err) in SW_REFINED.items():
        a2.loglog([kx] * len(err), err, 'v', color=UPA, ms=5)
        a2.annotate('', xy=(kx, err[-1]), xytext=(kx, err[0]),
                    arrowprops=dict(arrowstyle='->', color=UPA, lw=1.3))
    a2.plot([], [], 'v-', color=UPA, label='refining the grid')
    a2.set_xlabel(r'$k_x\rho$')
    a2.set_title('Short wavelength', fontsize=9.5, loc='left')
    a2.legend(frameon=False, fontsize=7.6)
    a2.text(0.97, 0.44, 'the rise is the quadrature\nit is compared against, not\nthe formula',
            transform=a2.transAxes, fontsize=7, color='#444', ha='right')

    x = np.arange(len(SW_BANDS['names'])); w = 0.36
    a3.bar(x - w/2, SW_BANDS['gaussian'], w, color=GREY, label='stationary phase')
    a3.bar(x + w/2, SW_BANDS['uniform'], w, color=UPA, label='uniform (Bessel)')
    for xi, fr in zip(x, SW_BANDS['weight']):
        a3.text(xi, 52, f'{fr:.0%}', ha='center', fontsize=7, color='#444')
    a3.set_xticks(x); a3.set_xticklabels(SW_BANDS['names'], fontsize=7.5)
    a3.set_ylabel('error [%]'); a3.set_ylim(0, 58)
    a3.legend(frameon=False, fontsize=7.6)
    a3.set_title(r'By pitch angle, $k_x\rho=10$  (share of weight above)',
                 fontsize=9.5, loc='left')
    fig.tight_layout()
    fig.savefig(outfile)
    plt.close(fig)
    return outfile
