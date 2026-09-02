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
- `rh_nl_kinetic.in` does not close: the residual sits near 0.6 and is
  insensitive to the choice of time window, so it is not a noise or windowing
  artefact.  Something about the two-kinetic-species case is genuinely
  inconsistent and needs investigating before a tolerance is set.
- The unexplained ~7e-3 floor above.
