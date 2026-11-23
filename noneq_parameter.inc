* INCLUDE FILE for the coupled chemistry+mixing solve (noneq).
* Compile-time numerics: edit, then rebuild. Run time inputs are the
* deck (tMAX, DTMIN, DTMAX, CNVCRI, CNVRAT) and the data files.
* The values in use are printed at the start of every run.

* DLSODE error control: relative, and absolute [cm-3].
      real*8 noneq_rtol, noneq_atol
      parameter (noneq_rtol = 1.0d-4, noneq_atol = 1.0d-1)

* Convergence test (state at t vs t/2): species below noneq_conv_atol
* [cm-3] are ignored; the floored norm uses n_tot*noneq_conv_floor_frac;
* species below n_tot*noneq_conv_skip_frac are skipped.
      real*8 noneq_conv_atol, noneq_conv_floor_frac
      real*8 noneq_conv_skip_frac
      parameter (noneq_conv_atol = 1.0d-1)
      parameter (noneq_conv_floor_frac = 1.0d-15)
      parameter (noneq_conv_skip_frac = 1.0d-20)

* Time grid: chunk ends are tMAX/2**(max_chunks-1), doubling up to tMAX.
      integer noneq_max_chunks
      parameter (noneq_max_chunks = 25)

* DLSODE step budget per chunk and Jacobian: 24 = analytic banded,
* 25 = finite difference (debug).
      integer noneq_max_steps, noneq_jac_mode
      parameter (noneq_max_steps = 200000, noneq_jac_mode = 24)

* INITAB floor (4th column): relaxation time [s]; the push starts at
* floor*(1+band) and is smooth (C1) so the solver does not chatter.
      real*8 noneq_floor_tau, noneq_floor_band
      parameter (noneq_floor_tau = 1.0d2, noneq_floor_band = 5.0d-2)
