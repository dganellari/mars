# Prescribed outlet pressure checks

`outletPressureLift` defaults off. Before each `runNsStep`, fill `d_pPhiTargetDof`
with `numOwnedDofs` physical-pressure (Pa) entries. Masked targets must be finite;
negative and zero values are valid, and unmasked entries are ignored. The update
uses `phi = target - p_old` and `p_new = p_old + phi`. All ranks must enable the
flag together and keep setup and masks fixed. Disabling it restores a zero
boundary increment; it does not reset existing pressure.

Supported: Pump/do-nothing outlets, one DDT corrector, Hypre's assembled K path
and GMRES. Retained boundary columns with identity boundary rows make the system
nonsymmetric, so the PCG override is ignored. The prescribed-pressure solve uses
coarse-level l1-Jacobi so BoomerAMG has row norms even without coarsening.
Rejected: matrix-free CG, assembled AddT, FEM/PISO, rotational correction,
flux-Neumann, pressure-driven inlet, pressure pin, average-outlet modes, and
ranks with no assembled pressure entries.

Invalid arrays, masked targets or configurations raise `std::invalid_argument`
on all ranks before state changes. Split-step callers must invoke
`validateOutletPressureLift` before advancing. This checks configuration; it does
not establish convergence or continuity accuracy for the K approximation.

```sh
cmake -S tests/outlet_pressure -B build-outlet-pressure \
    -DCMAKE_PREFIX_PATH=/path/to/candidate/install
cmake --build build-outlet-pressure
ctest --test-dir build-outlet-pressure --output-on-failure
```

Use a CUDA or HIP installation of this candidate with Hypre enabled. Configure
the `MPIEXEC_*` CMake variables for the allocated system's launcher.
The nine tests use real MARS/Hypre on public Tet4 boxes at 1/2/4 ranks: default-off
quiescence, whole-domain targets 2/-1/0 despite the PCG override, and collective
rejection without state changes. The 125-node and 729-node boxes exercise
one-level and multilevel AMG. They use SFC coordinates; exact Tet4 geometry and
full pump/FSI acceptance remain separate.
