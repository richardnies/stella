Rosenbluth-Hinton budget benchmarks
===================================

These benchmarks check the Rosenbluth-Hinton diagnostic against the code's own
time evolution, by asserting that the zonal-flow energy it implies obeys

    d E_RH / dt  =  sum_kx P_RH

with

    E_RH(t,kx) = |<RH_phi_I>|^2 / (2 |<RH_inertia>|^2) * (1 - Gamma0)
    P_RH(t,kx) = -Re[ i kx F <RH_phi_I>* ] / |<RH_inertia>|^2 * (1 - Gamma0)

where `<.>` is the dl/B field-line average and `F` is the sum of every RH flux
channel written to the netCDF file.  Both sides come from the same run, so no
reference data is needed and the tests keep working when the physics of a case
legitimately changes.

`rh_budget.py` holds the shared evaluation.  These definitions mirror
`stella_diagnostics/physics/rosenbluth_hinton.py` in stella_diagnostics_v2
(`get_E_RH_t_kx` / `get_P_RH`) and the two must be kept in step.

Run them with

    make rh-benchmarks

Scope
-----
The budget is expected to close in tokamak geometry with hyperdissipation and
the tertiary sponge switched off.  It does not close exactly: stella advances
the zonal mode with a particular discretisation while the RH projection that
annihilates the linear streaming and drift terms is a continuum construction.
Measurements on the linear collisional case put that mismatch near 1%: it is
insensitive to upwinding (1.45% -> 1.45% with all upwind coefficients zeroed)
and to resolution (no better at nzed = 96, nmu = 24, nvgrid = 72), but grows to
14% with `drifts_implicit = .false.`.  The tolerance is set at 5% accordingly.
