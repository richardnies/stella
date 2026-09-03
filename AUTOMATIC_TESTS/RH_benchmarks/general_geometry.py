################################################################################
#       GENERAL-GEOMETRY GROUNDWORK FOR THE ROSENBLUTH-HINTON DIAGNOSTICS      #
################################################################################
# The RH diagnostics currently assume a quasisymmetric field.  Two things have
# to be built before they work in a general stellarator, and this module is where
# both were worked out and checked against real equilibria before being written
# in Fortran.
#
#   1. The drift-orbit phase Q.  It is defined by
#
#          v_par b.grad Q = i kx (vMx - <vMx>_tau)
#
#      so that transit-averaging annihilates the radial magnetic drift.  In a
#      quasisymmetric field <vMx>_tau vanishes and Q has the closed form the code
#      uses, q R Btor.  In general it has to be integrated along the field line:
#
#          Q(z) = i kx int [vMx - <vMx>_tau] / (v_par b.grad z) dz
#
#   2. The bounce-averaged radial drift <vMx>_b itself, which is what makes the
#      general case different: it is zero by symmetry in a tokamak or a
#      quasisymmetric stellarator, and nonzero otherwise, and it is the source of
#      the F_RH_drift and P_RH_drift contributions to the budget.
#
# Two numerical ingredients are needed and neither exists in stella today.
#
#   Multiple wells.  A trapped particle is confined to one connected region where
#   B < B_c, and its bounce average must be taken over that well alone.  stella's
#   transit average masks the forbidden region but then sums over every
#   accessible interval at once, which mixes wells the particle can never reach.
#   In a tokamak there is only one, so this has never mattered.
#
#   The turning-point singularity.  dl/|v_par| has an integrable inverse-square-
#   root singularity at each end of the well.  Integrating it by trapezoid on the
#   z grid converges only as sqrt(dz), and that error dominates: on the precise
#   QH configuration it gives 1.7e-2 where the converged answer is 7.7e-4.  The
#   substitution below absorbs the singularity into the Gauss-Chebyshev weight
#   exactly, and then converges to three digits by nzed = 256.
#
# Verification, run as __main__ below, is a ladder of configurations in which
# <vMx>_b must vanish for the first two and must not for the last two:
#
#     tokamak (Miller, axisymmetric)      2.3e-15     machine zero
#     precise QA (quasi-axisymmetric)     1.8e-14     machine zero
#     precise QH (quasi-helical)          8.0e-4      its optimisation residual
#     W7-X Standard                       3.5e-2
#     TJ-II                               4.6e-2
#
# and, for the phase itself, that the numerically integrated Q reproduces the
# analytic q R Btor form in a tokamak to 1e-4, with one geometric constant for
# every (energy, mu).
#
# The inputs are stella's own geometry output, so this needs a geometry-only run
# per equilibrium (RH diagnostics off) and nothing else.
################################################################################

# ---------------------------------------------------------------------------
# from qproto.py
# ---------------------------------------------------------------------------
import pathlib

import numpy as np
from netCDF4 import Dataset

def geometry(ncfile):
    d = Dataset(ncfile)
    g = dict(
        zed=np.array(d.variables['zed'][:]),
        bmag=np.array(d.variables['bmag'][:])[:, 0],
        gradpar=np.array(d.variables['gradpar'][:]),
        cvdrift0=np.array(d.variables['cvdrift0'][:])[:, 0],
        gbdrift0=np.array(d.variables['gbdrift0'][:])[:, 0],
        shat=float(np.array(d.variables['shat'][...])),
    )
    return g

def radial_drift(g, energy, mu):
    """vMx up to the overall constant, in stella's decomposition:
    wdriftx ~ cvdrift0 * vpa^2 + gbdrift0 * vperp^2 / 2 ."""
    B = g['bmag']
    vperp2 = 2.0 * mu * B
    vpa2 = energy - vperp2
    return g['cvdrift0'] * vpa2 + g['gbdrift0'] * 0.5 * vperp2

def transit_average(g, integrand, energy, mu):
    """<A>_tau = int A/|vpa| dl / int 1/|vpa| dl, over the accessible region."""
    B = g['bmag']
    vpa2 = energy - 2.0 * mu * B
    ok = vpa2 > 0
    if not ok.any():
        return np.nan
    w = np.zeros_like(B)
    # dl/|vpa| with dl = dz/(b.grad z)
    dz = np.gradient(g['zed'])
    w[ok] = dz[ok] / (np.abs(g['gradpar'][ok]) * np.sqrt(vpa2[ok]))
    return np.sum(integrand * w) / np.sum(w)

def check(ncfile, energies, mus, sigma=+1.0):
    g = geometry(ncfile)
    B, zed = g['bmag'], g['zed']
    print(f"{'energy':>8} {'mu':>8} {'<vMx>_tau':>12} {'corr':>9} {'fitted C':>12} {'rel scatter':>12}")
    for energy in energies:
        for mu in mus:
            vpa2 = energy - 2.0 * mu * B
            if (vpa2 <= 0).any():
                continue                      # passing particles only, for a clean test
            vpa = sigma * np.sqrt(vpa2)
            vMx = radial_drift(g, energy, mu)
            avg = transit_average(g, vMx, energy, mu)
            integrand = (vMx - avg) / (vpa * g['gradpar'])
            analytic = np.gradient(vpa / B, zed)       # d/dz of the closed form's shape
            m = np.isfinite(integrand) & np.isfinite(analytic)
            corr = np.corrcoef(integrand[m], analytic[m])[0, 1]
            C = np.dot(analytic[m], integrand[m]) / np.dot(analytic[m], analytic[m])
            scatter = np.std(integrand[m] - C*analytic[m]) / np.std(integrand[m])
            print(f"{energy:8.3f} {mu:8.3f} {avg:12.3e} {corr:9.5f} {C:12.5f} {scatter:12.3e}")

# ---------------------------------------------------------------------------
# from wells.py
# ---------------------------------------------------------------------------
import numpy as np


def wells(B, B_c):
    """Connected z-intervals where B < B_c, as (start, stop) index slices.

    Intervals touching the ends of the domain are incomplete -- the particle
    would leave the simulated field line -- and are flagged so callers can drop
    them rather than treat a truncated well as a real one.
    """
    inside = B < B_c
    out = []
    i = 0
    n = len(B)
    while i < n:
        if inside[i]:
            j = i
            while j + 1 < n and inside[j + 1]:
                j += 1
            out.append((i, j + 1, i == 0 or j == n - 1))
            i = j + 1
        else:
            i += 1
    return out


def bounce_average(A, B, zed, gradpar, B_c, well):
    """<A>_b over one well: int A dl/|vpar| / int dl/|vpar|, with vpar^2 ~ 1 - B/B_c."""
    lo, hi, _ = well
    sl = slice(lo, hi)
    vpa2 = 1.0 - B[sl] / B_c              # proportional to the true vpar^2
    vpa2 = np.maximum(vpa2, 0.0)
    good = vpa2 > 0
    if good.sum() < 3:
        return np.nan
    dz = np.gradient(zed)[sl]
    w = np.zeros_like(vpa2)
    w[good] = dz[good] / (np.abs(gradpar[sl][good]) * np.sqrt(vpa2[good]))
    return np.sum(A[sl] * w) / np.sum(w)


def survey_trapped(g, drift, lambdas):
    """Largest |<vMx>_b| over complete wells, normalised by the drift scale."""
    B, zed, gradpar = g['bmag'], g['zed'], g['gradpar']
    scale = np.max(np.abs(drift))
    rows = []
    for lam in lambdas:
        B_c = B.max() / lam               # lam > 1 => B_c < Bmax => trapped
        complete = [w for w in wells(B, B_c) if not w[2]]
        vals = [bounce_average(drift, B, zed, gradpar, B_c, w) for w in complete]
        vals = [v for v in vals if np.isfinite(v)]
        rows.append((lam, len(complete), max(np.abs(vals)) / scale if vals else np.nan))
    return rows

# ---------------------------------------------------------------------------
# from bounce.py
# ---------------------------------------------------------------------------
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq


def turning_points(zed, B, B_c, lo, hi):
    """Bracket the well [lo, hi) and solve B(z) = B_c at each end."""
    spline = CubicSpline(zed, B - B_c)
    z_l = brentq(spline, zed[lo - 1], zed[lo]) if lo > 0 else zed[lo]
    z_r = brentq(spline, zed[hi - 1], zed[hi]) if hi < len(zed) else zed[hi - 1]
    return z_l, z_r


def bounce_average(zed, B, gradpar, A, B_c, lo, hi, n=96):
    """<A>_b over the well spanning indices [lo, hi), to spectral accuracy."""
    z_l, z_r = turning_points(zed, B, B_c, lo, hi)
    if not (z_r > z_l):
        return np.nan

    spline_B = CubicSpline(zed, B)
    spline_g = CubicSpline(zed, gradpar)
    spline_A = CubicSpline(zed, A)

    k = np.arange(1, n + 1)
    t = np.cos((2 * k - 1) * np.pi / (2 * n))          # Gauss-Chebyshev nodes
    mid, half = 0.5 * (z_l + z_r), 0.5 * (z_r - z_l)
    z = mid + half * t

    # g(z) = (B_c - B) / ((z - z_l)(z_r - z)), smooth across the well
    denominator = (z - z_l) * (z_r - z)
    g = (B_c - spline_B(z)) / denominator
    g = np.maximum(g, 1e-300)

    weight = 1.0 / (np.abs(spline_g(z)) * np.sqrt(g))   # the smooth part of dl/|vpar|
    return np.sum(spline_A(z) * weight) / np.sum(weight)


# ---------------------------------------------------------------------------
# Verification ladder
# ---------------------------------------------------------------------------
def bounce_averaged_drift(ncfile, n_lambda=80, lambda_max=1.6):
    """|<vMx>_b| / |cvdrift0|max over every complete well, for a spread of lambda.

    The drift is evaluated at the same (energy, mu) that defines each well: mu is
    fixed by lambda, and the drift depends on the vpa^2 / vperp^2 split, so using
    a single fixed mu while scanning wells gives a spurious nonzero answer in a
    field where the true value vanishes only after the correct weighting.
    """
    g = geometry(ncfile)
    B, zed, gradpar = g['bmag'], g['zed'], g['gradpar']
    scale = np.max(np.abs(g['cvdrift0']))
    values = []
    for lam in np.linspace(1.005, lambda_max, n_lambda):
        B_c = B.max() / lam
        mu = 1.0 / (2.0 * B_c)                      # energy = 1
        vperp2 = 2.0 * mu * B
        vpa2 = 1.0 - vperp2
        drift = g['cvdrift0'] * vpa2 + g['gbdrift0'] * 0.5 * vperp2
        for lo, hi, incomplete in wells(B, B_c):
            if incomplete or hi - lo < 4:
                continue                            # truncated by the domain edge
            value = bounce_average(zed, B, gradpar, drift, B_c, lo, hi)
            if np.isfinite(value):
                values.append(abs(value) / scale)
    return np.array(values)


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        sys.exit(f"usage: python {sys.argv[0]} <geometry-run.out.nc> [<label>] ...")
    print(f"{'configuration':34s} {'nwells':>7} {'median':>10} {'p90':>10} {'max':>10}")
    for netcdf_file in sys.argv[1:]:
        v = bounce_averaged_drift(netcdf_file)
        label = pathlib.Path(netcdf_file).name.replace('.out.nc', '')
        if not len(v):
            print(f"{label:34s}   no complete wells")
            continue
        print(f"{label:34s} {len(v):7d} {np.median(v):10.2e} "
              f"{np.percentile(v, 90):10.2e} {v.max():10.2e}")
