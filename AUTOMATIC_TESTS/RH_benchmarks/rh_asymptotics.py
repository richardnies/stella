"""Long-wavelength (order k_x^2) approximations to the Rosenbluth-Hinton weights.

Expanding the projection weight W = <J0 exp(-Q)>_tau exp(Q(z)) with
Q = i k_x delta_x(z), and splitting into parts even and odd in sgn(v_par) --
Q being odd -- gives

    even(W) - 1 = -<a^2>_tau/4 - (k_x^2/2)[ <dx^2>_tau + dx^2 - 2 dxbar dx ]
    odd(W)      =  i k_x ( dx - dxbar )

so the two halves carry different powers of k_x, O(k_x^2) and O(k_x), and can be
checked against the code separately.

Everything is evaluated on stella's own (v_par, mu, z) grid so that the averages
match the code's exactly and no velocity-space measure has to be reconstructed.

Which quantities need which order is the substantive point, and it is not
uniform.  A leading term survives only where nothing annihilates it:

    phi_RH        W -> 1 suffices; the leading term is the gyrocentre density.
    I_RH          leading term cancels against the 1 in (1 - J0 W): O(k_x^2).
    F_nonlinear   leading term cancels by quasineutrality, the gyrocentre charge
                  density being itself O(b): the even part is O(k_x^2), the odd
                  part O(k_x).
    F_collisional leading term cancels because the collision operator conserves
                  particles, int C[g] d^3v = 0: O(k_x^2) and O(k_x) likewise.
    F_drift       no conservation law removes it; W -> 1 suffices.
"""
import numpy as np
from netCDF4 import Dataset


def load(netcdf_file):
    """stella's grids and the field-line weight."""
    d = Dataset(netcdf_file)
    zed = np.array(d.variables['zed'][:])
    jac = np.array(d.variables['jacob'][:])[:, 0]
    bmag = np.array(d.variables['bmag'][:])
    bmag = bmag[:, 0] if bmag.ndim > 1 else bmag
    gds22 = np.array(d.variables['gds22'][:])
    gds22 = gds22[:, 0] if gds22.ndim > 1 else gds22
    w = (zed[1] - zed[0]) * jac.copy(); w[-1] = 0.0; w /= w.sum()
    return dict(zed=zed, jac=jac, bmag=bmag, gds22=gds22, weight=w,
                vpa=np.array(d.variables['vpa'][:]),
                mu=np.array(d.variables['mu'][:]),
                q=float(np.array(d.variables['q'][...])),
                shat=float(np.array(d.variables['shat'][...])))


def orbit_excursion(geo, rgeo, rhoc=0.5, smz=1.0, n_theta=256):
    """delta_x(z) and its transit average, on the (mu, vpa, z) grid.

    delta_x = (v_par/B) smz q I / rhoc from stella's closed form, with
    I = R B_tor a flux function.  The transit average uses the orbit measure
    dl/|v_par|, restricted to where the orbit can go; it vanishes for trapped
    particles, v_par/B being odd over a closed bounce.
    """
    b = geo['bmag']; zed = geo['zed']; dl = geo['jac'] * b
    fac = smz * geo['q'] * rgeo / rhoc
    vpa = geo['vpa']; mu = geo['mu']
    Bmax = b.max()

    nmu, nv, nz = len(mu), len(vpa), len(b)
    dx = np.zeros((nmu, nv, nz)); dxbar = np.zeros((nmu, nv))
    var = np.zeros((nmu, nv))
    access = np.zeros((nmu, nv, nz), dtype=bool)
    iz0 = nz // 2
    for i in range(nmu):
        for j in range(nv):
            E = vpa[j]**2 + 2 * mu[i] * b[iz0]
            v2 = E - 2 * mu[i] * b
            ok = v2 > 0
            access[i, j] = ok
            vp = np.sqrt(np.maximum(v2, 0.0))
            dx[i, j] = fac * np.where(ok, np.sign(vpa[j]) * vp / b, 0.0)
            trapped = mu[i] > 0 and (E / (2 * mu[i])) < Bmax
            if ok.sum() < 2:
                continue
            if trapped:
                #> The orbit measure dl/|v_par| is singular at the turning
                #> points, where v_par^2 vanishes linearly, so a plain sum over
                #> the grid converges as sqrt(dz).  Substituting
                #> z = mid + half cos(theta) puts the zero of v_par^2 at the
                #> endpoints of a cosine and removes it, which is the same
                #> device the Fortran bounce quadrature uses.
                B_c = E / (2 * mu[i])
                idx = np.where(ok)[0]
                zl, zr = _turning(zed, b, B_c, idx[0]), _turning(zed, b, B_c, idx[-1], right=True)
                th = (np.arange(n_theta) + 0.5) * np.pi / n_theta
                zt = 0.5 * (zl + zr) + 0.5 * (zr - zl) * np.cos(th)
                bt = np.interp(zt, zed, b); jt = np.interp(zt, zed, geo['jac'])
                v2t = np.maximum(E - 2 * mu[i] * bt, 0.0)
                #> dl/|v_par| with the Jacobian of the substitution folded in
                wt = jt * bt * 0.5 * (zr - zl) * np.sin(th) / np.sqrt(np.maximum(v2t, 1e-30))
                dxt = fac * np.sign(vpa[j]) * np.sqrt(v2t) / bt
                tau = wt.sum()
                if tau > 0:
                    var[i, j] = (wt * dxt**2).sum() / tau        # dxbar = 0 when trapped
            else:
                meas = np.where(ok, dl / np.maximum(vp, 1e-30), 0.0)
                tau = np.trapz(meas, zed)
                if tau > 0:
                    dxbar[i, j] = np.trapz(meas * dx[i, j], zed) / tau
                    var[i, j] = np.trapz(meas * (dx[i, j] - dxbar[i, j])**2, zed) / tau
    return dx, dxbar, var, access


def _turning(zed, b, B_c, i, right=False, n=200):
    '''Locate where B reaches B_c just outside the accessible index i.'''
    j = i + 1 if right else i - 1
    if j < 0 or j >= len(zed):
        return zed[i]
    zs = np.linspace(zed[min(i, j)], zed[max(i, j)], n)
    bs = np.interp(zs, zed, b)
    k = np.argmin(np.abs(bs - B_c))
    return zs[k]


def weight_asymptotics(geo, kx, rgeo, rhoc=0.5, smz=1.0):
    """Maxwellian- and field-line-averaged <even(W)>-1 and rms(odd(W))."""
    dx, dxbar, var, access = orbit_excursion(geo, rgeo, rhoc, smz)
    b = geo['bmag']; w = geo['weight']; vpa = geo['vpa']; mu = geo['mu']
    VP, MU = np.meshgrid(vpa, mu)
    #> a^2 = k_perp^2 vperp^2 / Omega^2, with k_perp^2 = kx^2 gds22/shat^2
    kperp2 = kx**2 * geo['gds22'] / geo['shat']**2

    num_e = num_o = den = 0.0
    for iz in range(len(b)):
        fM = np.exp(-(VP**2 + 2 * MU * b[iz])) * access[:, :, iz]
        a2 = kperp2[iz] * (2 * MU * b[iz]) / b[iz]**2
        d = dx[:, :, iz] - dxbar
        even = -0.25 * a2 - 0.5 * kx**2 * (d**2 + var)
        num_e += w[iz] * (fM * even).sum()
        num_o += w[iz] * (fM * (kx * d)**2).sum()
        den += w[iz] * fM.sum()
    return num_e / den, np.sqrt(num_o / den)
