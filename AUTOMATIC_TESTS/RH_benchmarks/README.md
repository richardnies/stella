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
| `rh_linear_collisional.in`     | linear, collisional damping of one zonal mode; the only source is `RH_fluxes_collisional` | whole run | 8.6e-3 |
| `rh_nl_adiabatic_electrons.in` | nonlinear, modified adiabatic electrons (flux-surface-average term kept) | t = 15..27 | 6.5e-3 |
| `rh_nl_adiabatic_ions.in`      | nonlinear, unmodified adiabatic electrons (plain Boltzmann) | t = 15..27 | 3.9e-3 |
| `rh_nl_kinetic.in`             | nonlinear, kinetic ions and kinetic electrons | t = 10..20, kx <= 1.1 | 5.9e-3 |
| `rh_nl_electromagnetic.in`     | nonlinear electromagnetic, apar and bpar channels | -- | skipped, deck unstable |

Two things decide whether a nonlinear case is meaningful, and both were learned the
hard way:

*The box has to contain modes the theory addresses.*  The Rosenbluth-Hinton
construction targets kx rho_i << 1/q (about 0.71 at q = 1.4).  A 4x4 box with
jtwist = 1 holds a single zonal mode at kx rho_i = 2.5, and the budget there does
not close.  `rh_nl_kinetic.in` therefore uses jtwist = 5 (dkx = 0.5) and is
compared over kx <= 1.1.  Measured across a wide box, in the clean window:

    kx rho_i      0.5       1.0       1.5       2.5
    nzed = 24   4.6e-3    2.6e-2    3.5e-1    4.0e-1
    nzed = 48   1.0e-2    2.8e-2    2.6e-1    3.9e-1

Doubling the parallel resolution does not move it, so this is the range of
validity of the construction rather than a resolution artefact -- and it matches
the derivation, which scopes itself to kx rho_i << 1/q.  Both species follow the
same curve; this is not an electron effect.  A case with a broad kx spectrum must
therefore be restricted with `kx_max`.

*The run has to still be resolved.*  These boxes have no cascade to saturate
into, so past the growth phase the fields run away and the budget stops meaning
anything -- the same modes that close to 1.2e-2 during growth give 8.8e-1 once
phi2 reaches 1e3.  Each case is compared over its clean window.  Use
`plot_rh_budget.py` to inspect a run and choose one.

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
- `rh_nl_electromagnetic.in` goes NaN from the second step at beta = 0.004, under
  both implicit and explicit streaming/mirror and with delt cut to 5e-3.  The
  deck needs stabilising before the budget can be assessed at all.  Worth noting
  when it is: stella's own `advance_ExB_nonlinearity` converts g to h only when
  `include_apar .or. include_bpar` (time_advance.f90), because electrostatically
  the correction cancels out of the flux -- so the electromagnetic case is
  exactly where the g/h distinction in eq (23) starts to matter.

- The unexplained ~7e-3 floor described above.

- `RH_phi_I` is integrated with weight `spec%z`, while `RH_inertia` and every RH
  flux use `spec%dens_psi0*spec%z` (rosenbluth_hinton.f90:581 against :301, :420
  and the rest).  Eq (13) and (16) of the derivation both carry n_s, so :581
  looks like the odd one out.  Every deck here has `dens = 1.0`, so the two agree
  and nothing exercises it, and it has been left alone pending a decision.
