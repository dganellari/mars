# Fixed graph Hypre lifecycle gate

Run from the repository root with an installed **sequential CPU Hypre**:

```sh
python3 tests/reference/openaccel/simple_performance/run_hypre_host_test.py \
  --hypre-prefix .local-worktrees/simple-hypre-audit/install
```

The script generates a header and executable under ignored
`.local-worktrees/simple-performance/hypre-host`. It retains the production GMRES
wrapper and shared packing kernels, adapting CUDA launches into CPU loops and
replacing device storage/MPI operations with single-process host equivalents.
Hypre matrix assembly, AMG setup, GMRES solve and destruction use the real library.
The compiler enables AddressSanitizer and UndefinedBehaviorSanitizer; the supplied
Hypre library itself may be uninstrumented. `--cxx` selects the compiler.
On macOS, match `MACOSX_DEPLOYMENT_TARGET` to the installed Hypre build if needed.

The gate checks changing nonsymmetric tridiagonal systems against manufactured
solutions, independently recomputed residuals, and fresh wrapper instances. It
covers nonzero-to-zero-to-nonzero entries, discarded local/map columns, stable
IJ/ParCSR/vector/GMRES/AMG/packing handles, relocated values, relocated graph
storage, explicit map invalidation, partition resizing, configuration invalidation,
mode transitions, the existing frozen-operator contract, Jacobi, zero RHS/guess,
rejected nonfinite/undersized inputs, empty owned partitions, and separate preconditioner matrices. Timing
scopes and removal of the two wrapper barriers are checked as well.

`enable_fixed_graph_updates()` requires unchanged CSR and local-to-global map
contents between solves. Call `invalidate_setup()` before changing either in
place. Pointer, dimensions and partition changes trigger collective rebuilds.
Every successful solve refreshes all matrix values and reruns GMRES/AMG setup;
there is no lagged AMG hierarchy. The host-vector overload and separate-K
preconditioner are deliberately rejected in this mode. Every rank must own at
least one row: the Hypre 3.1.0 device IJ reassembly path can skip its assembly
collective on an empty rank, so the wrapper rejects this case collectively.

Persistence covers the wrapper's graph packing and the outer IJ/ParCSR, vector,
GMRES and AMG handles. Hypre 3.1.0's device IJ SetValues/Assemble path still
rebuilds its inner diagonal/off-diagonal CSR arrays. Omitting Assemble would
leave queued values unapplied. This change does not claim to remove those
internal structural allocations or the necessary numerical AMG hierarchy setup.

`enable_timing()` reports per-call local wall API seconds through
`get_last_timing()`: `prepare_seconds`, `packing_seconds`, `setup_seconds`,
`solve_seconds`, `finish_seconds`. Packing is a subset of preparation. These
measurements add no device synchronization and are not kernel/device times.
`get_graph_build_count()`, `get_numeric_update_count()` and `get_setup_count()`
are cumulative. The fixed graph steady path currently has nine wrapper-level
Allreduces (including coordinated failure checks and the fused infinity norms),
plus Hypre's own communication; this is a cost to measure in MPI scaling tests.
When Hypre reports zero iterations, the fixed graph path recomputes the true
residual with a cached vector, public ParCSR matvec and global inner products.
This adds communication in that case and avoids Hypre 3.1.0's stale reported
residual on its exact-zero-residual early return. The gate covers both zero and
nonzero exact initial guesses after previous solves. For a zero RHS the reported
norm is absolute, matching the existing Hypre/wrapper convention.

This gate does not compile CUDA, exercise real MPI, validate Hypre's device
backend, or establish a speedup. Public native CUDA field/residual parity and
multi-rank lifecycle checks remain necessary on the user-run GPU environment.

## SIMPLE convergence with physical units

Both SIMPLE Hypre call sites enable `enable_true_residual_check(1e-13,1e-10)`.
Hypre still targets relative tolerance `1e-12`. Wrapper acceptance uses the
explicit residual `||b-Ax||_2 <= 1e-13 + 1e-10*||b||_2`, matching the existing
independent SIMPLE CSR check. A zero RHS uses the absolute tolerance. This
separates the tighter Krylov target from the application's acceptance limit:
Hypre can return successfully above its requested target. Nonfinite solutions
and API errors still fail. The no-argument overload retains strict acceptance
at the wrapper's configured relative tolerance (absolute when the RHS is zero).
The legacy `MARS_HYPRE_NULLX_RATIO` and `MARS_HYPRE_MAXX_RATIO` heuristics remain
unchanged for other callers; they do not apply in this mode. The caller must
still provide a pressure anchor. A small residual alone cannot detect a nullspace.
The existing SIMPLE check against its own CSR and exchanged solution also remains.
The wrapper forms `r = -Ax` with beta zero, then adds `b` with ParVectorAxpy.
Both operations use Hypre's compute stream. This avoids putting a runtime
device copy immediately before a matvec that reads and overwrites its destination.
It retains the same cached residual vector, inner products and coordinated error
check; no device-wide synchronization or host field copy is added. The acceptance
limits and independent MARS CSR check are unchanged.

The real-Hypre host test poisons that workspace with NaN and large finite values,
then compares repeated residuals with an independent original-CSR calculation.
It checks that neither input changes, including zero RHS/solution, changing
matrix values and RHS scaling from 1e-9 to 1e9, for GMRES and FlexGMRES. This
checks the algebra and storage contract, not CUDA stream ordering.

Hypre 2.x declares the exported `HYPRE_ParVectorAxpy` in its installed internal
`_hypre_parcsr_mv.h`; Hypre 3.x also declares it in the public header. The wrapper
includes the older header only for releases before 3.0 (or an unknown version).
The host adapter retains these includes and restores its instrumented MPI stubs
after Hypre's sequential aliases. The regression passes with sequential Hypre
2.32.0 and 3.1.0 under ASan/UBSan, including both Krylov backends. This checks
header compatibility and CPU behavior; the Daint CUDA rebuild remains separate.

A rejected solve in this mode prints the iteration count/limit, restart length,
Hypre's reported relative residual, the recomputed relative residual (absolute
when the RHS is zero), the Krylov target and the solve error code. Mixed mode
also prints the absolute residual and its acceptance limit. The report identifies
GMRES versus FlexGMRES, the RHS norm and the norm of the Krylov work residual.
The last quantity requires an extra global inner product on failure only; it is
unavailable after zero iterations and is not necessarily the final residual after
an early stop. No field is copied to the host. `MARS_HYPRE_VERBOSE=1` additionally
prints Hypre's iteration history. A residual mismatch does not, by itself,
identify whether the cause is roundoff, conditioning or an operator defect.

The distributed SIMPLE runner now unpacks and exchanges a dimensionally valid
rejected candidate before evaluating its original MARS CSR residual. A failure
prints `[simple-linear]` with that independent residual and both verdicts, then
stops. A passing CSR check never overrides a backend rejection. A missing
candidate is rejected collectively before it can be read. The existing verdict
collective carries both flags, so accepted solves gain no reduction or transfer.
The `marsSimpleLinearRejection{1,2,4}` host MPI tests cover valid candidates,
one-rank backend rejection, corrupted values, missing storage and nonfinite
values for both scalar pressure and three-component momentum systems.

The optional `MARS_HYPRE_FLEXGMRES=1` path uses FlexGMRES setters and getters.
Using GMRES APIs on that handle writes/reads a different Hypre data layout.
The real-Hypre regression exercises both backends, fresh and cached solves,
configuration getters, returned iterations and residual acceptance. This fixes
the optional API path; it does not establish the cause of the duct GPU mismatch.

In the reported duct-32 failure at momentum iteration 5 (Daint job 4884081),
the independent MARS residual and Hypre's Krylov work residual both equal
2.67277e-13, below the 4.05018e-11 acceptance limit. The former copy-and-matvec
wrapper check instead reports 1.88542e-7. This localizes the disagreement to
the wrapper's residual evaluation path; it does not prove a stream race.
The revised construction still needs the same cached GPU case to pass.

For the water/backflow fixture, the pseudo-time momentum diagonal is dominated
by `rho*V/(alpha_u*pseudo_dt)`. Thus `d=V/a` is about `6e-10` and pressure
matrix entries are about `1e-7`, while the mass RHS is about `125 kg/s`.
The independent host solve returns a pressure increment near `3.07003e9 Pa`
with relative residual `1.51e-14`. Comparing that increment directly to the
mass RHS with a fixed ratio ceiling rejects a valid linear solution. This
startup fixture is not a converged physical pressure prediction.

The wrapper regression scales both a nonsymmetric matrix and its RHS by
`1e-9`, `1` and `1e15`, preserving the manufactured solution and conditioning.
It exercises fresh and cached wrappers, zero RHS and a deliberately stalled
solve, using real CPU Hypre with ASan/UBSan. CUDA/MPI validation remains separate.

The mixed-tolerance regression forces Hypre to stop at supplied initial guesses
using an intentionally loose absolute stopping tolerance in the test only. It
checks actual returned fields against an independent CSR residual: candidates
above the Krylov target but inside SIMPLE's limit pass, candidates outside fail,
and the strict mode still rejects them. Fresh and cached wrappers cover small
and large RHS norms, zero RHS, and invalid tolerances. This is a controlled
early-stop test, not a reproduction of Hypre's GPU stagnation path.
