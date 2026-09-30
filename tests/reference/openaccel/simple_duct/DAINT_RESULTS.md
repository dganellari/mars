# Daint duct validation, 2026-09-29

**PASS:** upwind Tet4 SIMPLE on the public rectangular duct, levels 8/16/32,
each on 1/2/4 GPU ranks. All nine runs converged and exited zero. Rechecking
the downloaded meshes, logs, metrics and complete distributed fields with
`duct_compare.py` at `0843226c` reproduced the cluster report, apart from its
run-directory heading. Final-field, nonlinear and refinement thresholds were unchanged.

## Evidence and provenance

Results are saved at
`/capstor/scratch/cscs/gandanie/simple-duct-native-G7Xqla`.
The recheck verifies mesh hashes, run identity, controls, exit status, final
metrics, field completeness, rank parity and refinement. The GPU solves ran
on Daint; the independent recheck reads their saved output and is not another
GPU execution.

- Levels 8/16 were preserved from `simple-duct-upwind-Ng67FX`. Their backend
  selection is not recorded in these logs; they must not be described as new
  native-SpMV runs.
- Level 32 records source `4e31ae15baed35a2462c6365fde93d801efd6869`, Hypre
  runtime 2.33.0 with CUDA/cuSPARSE/cuBLAS enabled, and explicit
  `MARS_HYPRE_SPMV_VENDOR=0`. It uses plain GMRES, cache and halo overlap enabled,
  distributed field output, and normal asynchronous launches.
- Recorded executable SHA-256:
  `55f9054880cba3fa7c58e67e3e17866c5c6ff5f1773b962eee81d9aae5694bd3`.
  The executable itself was not copied for the recheck.
- Verified level-32 Exodus SHA-256:
  `571c25a415c7e9dddbc468cdae41b345933261cd70234eda6cfee2e2e65d8f19`.
- Commit `0843226c` subsequently makes native GPU SpMV the default. These runs
  exercised the same selection explicitly in the earlier binary.

The case uses rho=1, mu=0.1, inlet speed=0.1, outlet pressure=0, reference
length=1, alpha_u=alpha_p=0.3, alpha_mass=0.75, beta=0.05 and pseudo_dt=0.01.
Nonlinear residual, mass and change tolerances are 1e-8; rank parity is 1e-6.
Hypre targets 1e-12 relative residual. Both linear checks must satisfy the
existing SIMPLE bound `||b-Ax|| <= 1e-13 + 1e-10*||b||`.

From a configured build directory, recheck without running a solver:

```bash
python3 ../tests/reference/openaccel/simple_duct/duct_compare.py study \
  /capstor/scratch/cscs/gandanie/simple-duct-native-G7Xqla \
  --levels 8,16,32 --ranks 1,2,4 --advection upwind
```

## Numerical results

Each row below applies to all three rank counts. Profile L2 error is relative
to the analytical profile's L2 norm; wall slip uses the analytical centreline
speed. G error is relative to the analytical gradient, 0.174915632 Pa/m.

| Cells across height | Iterations | Profile L2 error | G error | Wall slip |
|---|---:|---:|---:|---:|
| 8 | 3139 | 7.4461% | -13.985% | 10.857% |
| 16 | 4484 | 4.1965% | -8.0783% | 5.8314% |
| 32 | 10656 | 2.2504% | -4.3748% | 3.0418% |

The finest G is 0.16726345 Pa/m. Successive analytical-error orders are
0.83/0.90 for the profile and 0.79/0.88 for G. The three-grid G estimate has
order 0.67 and a 7.782% estimated band containing the 4.375% finest error;
its extrapolation is 1.851% above the analytical value. The profile estimate
also passes. These are finite-grid estimates, not a proof of asymptotic order.
The earlier host 4/8/16 study still fails its three-grid G-order threshold;
this finer 8/16/32 result does not change that historical verdict.

Across all six multi-rank comparisons, maximum velocity/U disagreement is
1.752e-15 and maximum absolute-pressure/(rho U^2) disagreement is 6.328e-13.
Each level has identical iteration counts across ranks. At level 32 the final
momentum residual is about 9.987e-9, continuity about 1.21e-14 and boundary
imbalance at most 2.78e-16. Small algebraic residuals do not remove the spatial
errors in the table.

## Backend workaround and remaining limits

Native Hypre SpMV computes the same sparse matrix-vector operation on the GPU.
It changes neither the CFD equations nor the residual acceptance policy.
The wrapper and independent MARS CSR checks remain active; failed checks are
not accepted after a retry. This study validates that backend choice for this
case, including nonlinear convergence and rank agreement.

The earlier vendor-enabled runs produced inconsistent residual evaluations.
The exact defect in the Hypre/cuSPARSE execution path remains unidentified;
workspace reuse is a hypothesis, not a demonstrated cause. A root-cause repair
requires a reproducible public matrix/vector case, identification of the
violated operation or API contract, and regression tests for the corrected
vendor path. See the [residual audit](../simple_performance/HYPRE_HOST_TEST.md#gpu-spmv-backend-probe).

This study covers upwind on one node. High-resolution duct refinement, level 64,
arbitrary meshes and pump accuracy are not established. The saved rank-maximum
level-32 iteration-loop times are 1625.83, 1938.42 and 2095.13 seconds on 1/2/4
ranks: these runs do not show a speedup. Node identities are absent from these
logs, and this is not a controlled performance or multi-node scaling study.
