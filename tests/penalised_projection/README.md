# Penalised projection checks

`setPenalisedProjection(s, d_W, invEta)` turns on a Brinkman penalty `W / eta` in the pressure
projection. The corrector becomes `u = u** - (dtEff / rho) beta grad(phi)` with
`beta = c V / (c V + W invEta)`, `c = 1 / dtEff`, and the pressure operator takes the same nodal
coefficient (the Galerkin K with the element mean of beta).

**Inputs.** `d_W` has one finite, nonnegative value per local node (owned and ghosts, halo-complete),
in device memory. `invEta` is finite and nonnegative. `dt` is finite and positive. The setter checks
its arguments and `runNsStep` checks all of them (`validatePenalisedProjection`) before anything
changes. Invalid input raises `std::invalid_argument` on every rank, and the state is left as it was.

**Supported:** Hypre's assembled Galerkin operator on tets, with DDT (assembled) or K.
**Rejected:** matrix-free DDT CG, `MARS_HYPRE_USE_DDT`, an outlet beta, the FEM projection, and the
Rhie-Chow/VMS options (`useRhieChow`, `useVMSStab`, `useRcImplicit`, `useRcOnly`).

**The penalty in the velocity solve.** MARS does not add `W / eta` to the velocity matrix. The
caller does that, as ldp_fsi does.

```sh
cmake -S tests/penalised_projection -B build-penalised-projection \
    -DCMAKE_PREFIX_PATH=/path/to/candidate/install
cmake --build build-penalised-projection
ctest --test-dir build-penalised-projection --output-on-failure
```

Use a CUDA or HIP installation of this candidate with Hypre enabled. Configure the `MPIEXEC_*`
CMake variables for the system's launcher.

The twenty-one tests run seven modes at 1, 2 and 4 ranks on the public Tet4 box of
`tests/outlet_pressure`, with 6³ cubes: on 4³, BoomerAMG under GMRES fails the first pressure solve
with inflow (Hypre error 1), with or without the penalty.

| mode | checks |
|---|---|
| `baseline` | W = 0 leaves the setup operator (to 1e-14), and three converged inflow steps match a stepper without the penalty to the pressure solve's tolerance (checked at 1e-8) |
| `varying` | W piecewise constant in slabs, with jumps: every interior plane's rows sum to the layered operator's flux, to 1e-4 (MARS's CVFEM K meets the Galerkin identity to ~1e-5 even with a uniform beta); the same identity with the largest nodal beta per tet in place of the mean must miss by more than 1e-2 |
| `update` | W = w V gives a uniform beta, so K becomes beta K0; after changes of W, eta and dt, and of dt alone |
| `bdf2` | the first step at dtEff = dt, the following ones at 2 dt / 3 |
| `flow` | with inflow and beta = 1/2: the pressure solve converges, and continuity is no worse than without the penalty (two separate checks) |
| `failed` | a deliberately failed pressure solve (one GMRES iteration, an unreachable acceptance residual) returns -2, and the convergence check rejects it |
| `reject` | NaN or negative W, a bad invEta, dt zero, negative or NaN, W changed behind the setter, a wrong size, and unsupported configurations (including `useRhieChow`); each set on one rank only, so the refusal must be collective |

**Convergence** means a nonnegative return below the iteration cap: `solveOneComponent` returns -2
when a solve fails, which the cap alone would accept.
