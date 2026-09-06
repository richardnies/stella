# namelist `parameters`

Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`beta`  | float | 0.0 | Plasma \\( \beta \\). Currently has no effect.
`zeff`  | float | 1.0 | Effective charge number for use with *effective* electron-ion and electron-impurity collisions in the Fokker-Planck collision operator (see `ecoll_zeff`).
`tite`  | float | 1.0 | Ratio of ion to electron temperature, \\( T_\mathrm{i}/T_\mathrm{e} \\). Used in quasineutrality when adiabatic species is used.
`nine`  | float |1.0 | Ratio of ion to electron density, \\( n_\mathrm{i}/n_\mathrm{e} \\). Used in quasineutrality when adiabatic species is used.
`rhostar`  | real | -1.0 | The gyrokinetic expansion parameter \\( \rho_\mathrm{th,ref}/a_\mathrm{ref} \\). For effects beyond the flux-tube limit (full-flux-surface, radially global, neoclassical terms, etc...). Overwritten if `irhostar` is positive.
`irhostar`  | real | -1.0 | Sets `rhostar = 1.0 / irhostar` if positive.
`vnew_ref`  | real | -1.0 | Reference collision frequency. Various input options will overwrite this if it is negative.
`g_exb`  | real | 1.0 | Equilibrium \\( \boldsymbol{E \times B} \\) shear rate. More specifically, \\( \gamma_\boldsymbol{ E \times B} = (r/q) (\textrm{d}\omega / \textrm{d}r) R_0/\sqrt{2}v_\mathrm{th,ref}\\). Uses the Hammett wavenumber shift method, with nonlinear corrections proposed by McMillan.
`g_exbfac`  | real | 1.0 | Prefactor for perpendicular component of equilibrium \\( \boldsymbol{E \times B} \\) flow shear. Setting to 0.0 turns this component off.
`omprimfac`  | real | 1.0 | Prefactor for parallel component of equilibrium \\( \boldsymbol{E \times B} \\) flow shear. Setting to 0.0 turns this component off.
`omprimfac_RH`  | real | 0.0 | Prefactor for a Rosenbluth--Hinton parallel flow shear drive, which unlike `omprimfac` carries the velocity-space dependence of the Rosenbluth--Hinton state. Setting to 0.0 turns it off.
`omprimfac_PS`  | real | 0.0 | Prefactor for a Pfirsch--Schl&uuml;ter parallel flow shear drive. Setting to 0.0 turns it off.
`RH_analytic_drift_phase`  | boolean | geometry-dependent | How the Rosenbluth--Hinton drift-orbit phase \\( Q \\) is obtained: `true` uses the closed form the theory gives for a quasisymmetric field, `false` integrates \\( Q \\) along the field line from the magnetic drifts. The two agree in a tokamak, where the closed form applies. Left unset it follows the geometry, i.e. analytic under Miller and numerical under VMEC.
