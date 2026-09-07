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
| `rh_nl_modified_adiabatic.in` | nonlinear, modified adiabatic electrons (flux-surface-average term kept) | t = 15..27 | 6.5e-3 |
| `rh_nl_adiabatic.in`      | nonlinear, unmodified adiabatic electrons (plain Boltzmann) | t = 15..27 | 3.9e-3 |
| `rh_nl_kinetic.in`             | nonlinear, kinetic ions and kinetic electrons | t = 10..20, kx <= 1.1 | 4.0e-3 |
| `rh_nl_electromagnetic.in`     | nonlinear electromagnetic, apar and bpar channels | t = 6..16, kx <= 1.1 | 5.9e-2 |

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
- stella goes NaN when `include_bpar` and `include_collisions` are both on in a
  nonlinear run.  The fields decay rather than blow up and then NaN abruptly.
  This reproduces with the RH diagnostics switched off and on the branch point,
  so it is not from this work.  `rh_nl_electromagnetic.in` therefore runs
  collisionless; nothing is lost, since the nonlinear cases assert the nonlinear
  channel and the collisional channel is verified on its own by the linear
  benchmark.

- The unexplained ~7e-3 floor described above.

Stellarator benchmarks
----------------------
`test_2_rh_budget_stellarator/` runs the same assertion in five VMEC equilibria
rather than a Miller tokamak, on the `alpha0 = 0.7` field line.  What it adds is
the magnetic-drift channel: in an axisymmetric field the bounce-averaged radial
drift vanishes and `P_RH_drift` with it, so `test_1` never exercises the channel
that the stellarator generalisation exists for.  Here it carries between 1% and
39% of the peak of the linear budget.

ITER is the control.  It is a VMEC equilibrium but axisymmetric, so it runs the
same code path with a drift channel that should stay negligible, and it would
catch a generalisation that failed to reduce to the tokamak limit.

| configuration | drift share | linear | nl, adiabatic electrons | nl, adiabatic ions |
|---------------|-------------|--------|-------------------------|--------------------|
| ITER (axisym) |  4%         | 2.0e-2 | 3.2e-3  (t = 34..44)    | 7.7e-3  (t = 34..44) |
| W7-X standard |  1%         | 3.8e-2 | 1.4e-3  (t = 30..40)    | 2.1e-3  (t = 30..40) |
| QA            | 32%         | 9.3e-2 | 7.4e-3  (t = 28..38)    | 7.4e-3  (t = 28..38) |
| QH            | 13%         | 2.6e-2 | 3.5e-2  (t = 34..44)    | 4.0e-2  (t = 30..40) |
| TJ-II         | 39%         | 7.0e-2 | 1.5e-2  (t = 28..38)    | 8.5e-3  (t = 30..40) |

Two differences from the tokamak decks are worth knowing.

*The windows sit later.*  These ITG modes grow more slowly, reaching nonlinear
amplitude only around t = 25, so the clean window is roughly t = 28..44 rather
than t = 10..27.  The degradation past it is the same runaway as in `test_1`:
for W7-X the nonlinear channel closes to 1.4e-3 over t = 30..40, 8.0e-2 over
t = 45..55, and 3.9e-1 over t = 50..58.

*The field line matters.*  At `alpha0 = 0` every one of these configurations
sits near a symmetry where the bounce-averaged drift very nearly vanishes.  Both
sides of the budget then fall to the time-integration noise floor -- for ITER
the drive `rms(P)/E` is 3e-4 and `rms(dE/dt)` exceeds `rms(P)` tenfold -- and the
comparison measures noise rather than physics.  The decks therefore use
`alpha0 = 0.7`.  `plot_rh_stellarator.py` reports the drive alongside the slope
and labels a panel rather than quoting a number when it is below a per-mille.

The VMEC equilibria are not in the repository, being some 24 MB together.  Each
test skips if its `wout` file is absent; drop them beside the decks to enable
the suite.
