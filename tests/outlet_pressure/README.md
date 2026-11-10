# Prescribed outlet pressure checks

`outletPressureLift` supplies physical pressure (Pa), including negative or zero
values, on the existing owned pressure-boundary DOFs. The caller sets a finite
target on each masked DOF before every `runNsStep`; unmasked entries are ignored.
The increment is `phi = target - p_old`, followed by `p_new = p_old + phi`.
All ranks must enable the flag together and keep the setup and masks fixed.
Disabling the flag restores the legacy zero-increment boundary behavior; it does
not reset an existing pressure to zero.

This version supports Pump/do-nothing outlets, one corrector, DDT selection with
Hypre's assembled K path, and GMRES. The K matrix retains boundary columns, so
its interior equations already include the nonzero boundary increment. Its
identity boundary rows make the full matrix nonsymmetric. The PCG environment
override cannot override GMRES for this solve. BoomerAMG uses l1-Jacobi on the
coarse level too, so Hypre builds the required row norms when no coarsening
occurs. This setting applies only to the prescribed-pressure solve.
The assembled AddT path has removed boundary columns without an interior RHS
lift, and matrix-free DDT uses CG on a row-modified operator; both are rejected.
FEM/PISO, rotational pressure correction, flux-Neumann, pressure-driven inlet,
pressure pin and average-outlet modes are also outside this contract. Ranks with
no assembled pressure entries are rejected.

Validation runs collectively before `runNsStep` changes pressure, velocity or
history. A malformed array, nonfinite masked target or unsupported mode on any
rank raises `std::invalid_argument` on all ranks. The low-level split-step API
requires callers to invoke `validateOutletPressureLift` before advancing any
sub-step. This guard is a configuration check, not a convergence/step-acceptance
API, and it does not establish continuity accuracy for the K approximation.

The host algebra test compiles the actual boundary kernels extracted from the
current header. It checks constant pressure, return to zero, unequal signed
outlet pressures, retained columns and the missing-lift negative control.

```sh
cmake -S tests/outlet_pressure -B build-outlet-pressure
cmake --build build-outlet-pressure
ctest --test-dir build-outlet-pressure --output-on-failure
```

With an installation of this candidate built with GPU and Hypre support, enable
`MARS_PRESSURE_GPU_CHECKS=ON` and point `CMAKE_PREFIX_PATH` at that installation.
The nine additional tests use a public procedural Tet4 box at 1/2/4 ranks. They
check default-off quiescence, whole-domain targets 2/-1/0, resistance to the PCG
override, and collective rejection without state changes. They use the existing
SFC coordinates; exact Tet4 geometry is a separate dependency. The original
125-node box exercises a one-level AMG hierarchy; the refined 729-node box
also checks pressure propagation above the coarse-size threshold. Configure
`MPIEXEC_EXECUTABLE`, `MPIEXEC_NUMPROC_FLAG`, `MPIEXEC_PREFLAGS` and
`MPIEXEC_POSTFLAGS` for the allocated system's launcher. GPU execution is required
before accepting this as a working fluid path; host algebra and syntax checks
alone do not validate it.
