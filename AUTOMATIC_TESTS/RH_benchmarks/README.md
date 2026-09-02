Rosenbluth-Hinton budget benchmarks
===================================

These benchmarks check the Rosenbluth-Hinton diagnostic against the code's own
time evolution, by asserting that the zonal-flow energy it implies obeys

    d E_RH / dt  =  sum_kx P_RH

with, `<.>` being the dl/B field-line average and `F` the sum of every RH flux
channel written to the netCDF file,

    E_RH(t,kx) = |<RH_phi_I>|^2 / (2 |<RH_inertia>|^2) * (1 - Gamma0)
    P_RH(t,kx) = -Re[ i kx F <RH_phi_I>* ] / |<RH_inertia>|^2 * (1 - Gamma0)

Both sides come from the same run, so no reference data is needed and the tests
keep working when the physics of a case legitimately changes.

P_RH splits into a nonlinear part (the phi/apar/bpar fluxes) and a collisional
part, and `get_rh_budget` returns them separately.  The linear case is checked on
the total, since there the collisional channel is the only source.  The nonlinear
cases are checked on the nonlinear channel alone, by subtracting the collisional
channel from dE_RH/dt -- the collisional channel is already verified in isolation
by the linear benchmark, and isolating the nonlinear one tightens the check.

`rh_budget.py` holds the shared evaluation.  These definitions mirror
`stella_diagnostics/physics/rosenbluth_hinton.py` in stella_diagnostics_v2
(`get_E_RH_t_kx` / `get_P_RH`) and the two must be kept in step.

Run the benchmarks with

    make rh-benchmarks

and plot any run with

    python plot_rh_budget.py <run>.out.nc --out <dir> [--time-min T --time-max T]

which draws E_RH(t) and dE_RH/dt against sum P_RH, with the running local
residual on a twin axis.

Cases
-----
| deck | what it covers | window | residual |
|------|----------------|--------|----------|
| `rh_linear_collisional.in`     | linear, collisional damping of a single zonal mode; the only source is `RH_fluxes_collisional` | whole run | 8.6e-3 |
| `rh_nl_adiabatic_electrons.in` | nonlinear, modified adiabatic electrons (flux-surface-average term kept) | t = 15..28 | 3.1e-2 |
| `rh_nl_adiabatic_ions.in`      | nonlinear, unmodified adiabatic electrons (plain Boltzmann) | t = 15..28 | 5.9e-3 |
| `rh_nl_kinetic.in`             | nonlinear, kinetic ions and electrons | -- | does not close, see below |
| `rh_nl_electromagnetic.in`     | nonlinear electromagnetic, apar and bpar channels | -- | see below |

The nonlinear decks are run in a 4x4 box, which has no cascade to saturate
into: the ITG mode grows and then the fields run away.  They are therefore
compared over the window in which the zonal flow is genuinely nonlinearly
driven and the run is still well resolved, rather than over the whole run.

Scope and the residual floor
----------------------------
The budget is expected to close in tokamak geometry with hyperdissipation and
the tertiary sponge switched off.  It does not close exactly, and the mismatch
has two parts.

*Explained.*  `RH_fluxes_collisional` is evaluated as `(g^{n+1} - g^n)/code_dt`,
a first-order difference, so it is only an O(delt)-accurate estimate of the
collisional rate.  Refining delt on a collisions-only variant of the linear case
takes the residual from 1.3e-2 at delt = 0.2 to 7.0e-3 at delt = 0.0125.

*Not explained.*  A floor near 7e-3 survives delt -> 0.  It is not numerical
upwinding (zeroing `zed_upwind`, `vpa_upwind` and `time_upwind` leaves it at
1.45e-2 -> 1.45e-2), not parallel or velocity resolution (no better at
nzed = 96, nmu = 24, nvgrid = 72), and not a failure of the transit-average
projection to annihilate the linear streaming and drift terms -- switching both
off, so that collisions are the only dynamics, leaves it at 7.1e-3.  Note that
switching off only one of streaming or drifts makes the residual much worse
(6.5), which is consistent with the two having to cancel against each other in
the RH construction.

Tolerances are 5%, set with that floor in mind.

Known gaps
----------
- `rh_nl_kinetic.in` does not close, and the failure is entirely in the electron
  species.  The budget reduces to the charge relation `dRH_phi_I/dt = -i kx F`,
  and checking that per species over the clean growth phase gives

      ions       |-i kx F| / |dRH_phi_I/dt| = 0.989    residual 1.1e-2
      electrons  |-i kx F| / |dRH_phi_I/dt| = 4 - 6    residual 0.78 - 0.84

  A linear collisionless two-species run isolates it further.  There the RH
  fluxes are identically zero, so any change in RH_phi_I is pure
  transit-average annihilation error:

      nzed        24        48        96
      ions      2.0e-2    4.1e-3    3.0e-4     converges
      electrons 8.3e-3    4.3e-3    3.3e-3     plateaus

  The electron error is also insensitive to velocity resolution (unchanged
  across nvgrid 24-96 and nmu 12-24), and only halves (3.3e-3 -> 1.9e-3) with
  the magnetic drifts switched off.  So it is a defect in the electron RH
  response rather than a resolution or setup problem.

  What makes it concrete is that the electron species is suppressed in the RH
  inertia but not in the quantities built from the same transit average.  Field-
  line averaged, at kx = 2.5:

      RH_inertia      electrons / ions = 0.0085     correctly negligible
      RH_phi_I        electrons / ions = 1.03
      nonlinear flux  electrons / ions = 1.0 - 2.7

  The inertia integrand carries a factor (1 - J0 <J0 exp(-iQ)>_tau), which tends
  to zero for electrons because J0 -> 1 and Q -> 0 as the electron gyroradius and
  drift vanish.  RH_phi_I and the fluxes are weighted by the bare
  <J0 exp(-iQ)>_tau exp(iQ), which tends to one instead, so the electron
  contribution to those is not suppressed at all.  Dropping the electron species
  from the sum makes the budget close: ion-only residual 1.1e-2, against 0.64 for
  the two summed and 0.76 for electrons alone.  Whether the electron weighting
  should carry the same suppression is a physics decision, not a coding one.

- `RH_phi_I` is integrated with weight `spec%z`, while `RH_inertia` and every RH
  flux use `spec%dens_psi0*spec%z` (rosenbluth_hinton.f90:581 against :301, :420
  and the rest).  Every deck here has `dens = 1.0`, so the two agree and the
  inconsistency is invisible, but it would break the budget for any run with a
  non-unit density.
- The unexplained ~7e-3 floor above.
