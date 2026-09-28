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

Both SIMPLE Hypre call sites enable `enable_true_residual_check()`. Acceptance
uses an explicitly recomputed `||b-Ax||_2 / ||b||_2` with the existing tolerance,
or the absolute residual for a zero RHS, and rejects nonfinite solutions.
The legacy `MARS_HYPRE_NULLX_RATIO` and `MARS_HYPRE_MAXX_RATIO` heuristics remain
unchanged for other callers; they do not apply in this mode. The caller must
still provide a pressure anchor. A small residual alone cannot detect a nullspace.
The existing SIMPLE check against its own CSR and exchanged solution also remains.
This adds a device ParCSR matvec, global inner products and a coordinated error
check per solve; scalar reductions and API control are host work, not field copies.

A rejected solve in this mode prints the iteration count/limit, restart length,
Hypre's reported relative residual, the recomputed relative residual (absolute
when the RHS is zero), the tolerance and the solve error code. This failure-only
report adds no reductions or field transfers. `MARS_HYPRE_VERBOSE=1` additionally
prints Hypre's iteration history. A residual mismatch does not, by itself,
identify whether the cause is roundoff, conditioning or an operator defect.

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
