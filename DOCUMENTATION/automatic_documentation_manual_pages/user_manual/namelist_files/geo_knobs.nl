# namelist `geo_knobs`

Variable | Type | Default | Description
-------- | ---- | ------- | -----------
`geo_option` | string | `'local'` | Selects the geometry module. Should be one of <ul><li>`local` Miller equilibrium.</li><li>  `vmec` VMEC stellarator equilibrium. Requires VMEC netcdf input file. </li><li> `input.profiles` Reads in General Atomics `input.gacode` file. **This may need to be updated for newer `gacode` files.**</li><li>`miller`same as `local`</li><li>`default` same as `miller`</li></ul> 
`geo_file` | string | `'input.geometry'` | input file used to overwrite selected geometric coefficients below. File uses the same formatting as `.geometry` output.  
`q_as_x` | boolean | `radial_variation` | Uses the safety factor \\( q \\) as the radial coordinate, rather than \\( \psi \\). Used for radially global simulations.
`set_bmag_const` | boolean | `false` | Sets \\( B \\) uniformly to its value at the outboard midplane
`overwrite_bmag` | boolean | `false` | overwrite \\( B \\).
`overwrite_gradpar` | boolean | `false` | overwrite \\( \boldsymbol{ b \cdot \nabla z} \\).
`overwrite_gds2` | boolean | `false` | overwrite \\(  \lvert \nabla \alpha \rvert^2 (\textrm{d}\psi / \textrm{d}r)^2 \\).
`overwrite_gds21` | boolean | `false` | overwrite \\( \lvert \nabla q \cdot \nabla \alpha \rvert (\textrm{d}\psi / \textrm{d}r)^2 \\).
`overwrite_gds22` | boolean | `false` | overwrite \\( \lvert \nabla q \rvert^2 (\textrm{d}\psi / \textrm{d}r)^2 \\).
 `overwrite_gds23` | boolean | `false` | overwrite \\(  \nabla \theta \cdot [\nabla \alpha \times (\nabla r \times \nabla \alpha)] (\textrm{d}\psi / \textrm{d}r)^2 / B^2 \\).
 `overwrite_gds24` | boolean | `false` | overwrite \\( \nabla \theta \cdot [\nabla r \times (\nabla r \times \nabla \alpha)] (\textrm{d}\psi / \textrm{d}r)^2 (q/r) / B^2\\).
 `overwrite_gbdrift` | boolean | `false` | overwrite \\( 2 (\boldsymbol{b} \times \nabla B \cdot \nabla \alpha) (\textrm{d}\psi / \textrm{d}r) /B^2 \\).
  `overwrite_gbdrift0` | boolean | `false` | overwrite \\( 2 (\boldsymbol{b} \times \nabla B \cdot \nabla q) (\textrm{d}\psi / \textrm{d}r) /B^2  \\).
 `overwrite_cvdrift` | boolean | `false` | overwrite \\( 2(\textrm{d}\psi / \textrm{d}r)[\boldsymbol{b} \times (\boldsymbol{b \cdot \nabla b})]\cdot \nabla \alpha / B \\)
`one_sided_dbdz_at_ends`  | boolean | `false` | Evaluate \\( \partial B/\partial z \\) at the two ends of the flux tube with a second-order one-sided difference instead of a centred difference that reaches across the join, on a tube whose two ends do not carry the same \\( B \\). The centred difference presumes the two ends are the same physical point, which is true of a tube that closes on itself and false of one that does not: on a W7-X line with a one per cent mismatch in \\( B \\) it returns the same value at both ends, wrong by more than the whole scale of \\( \partial B/\partial z \\). This feeds the mirror term at every \\( z \\) including the ends. A tube that does close is bit-for-bit unchanged either way, since the periodicity of \\( B \\) is tested rather than assumed from the flag. **Note that correcting the coefficient does not stabilise the zonal mode on a non-closing tube -- it makes it grow faster** (a factor of five on the W7-X \\( \alpha_0 = 0.7 \\) line with kinetic electrons), because the erroneous derivative was incidentally damping it; the remedy for a zonal-flow calculation is a flux tube that closes.
