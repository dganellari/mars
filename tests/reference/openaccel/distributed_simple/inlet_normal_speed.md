# Constant normal-speed inlet

The fixed-frame, steady Tet4 SIMPLE inlet follows OpenAccel `0d69041`:
`src/field/nodeField/vector/transport/flow/velocity.cpp`, constant `normal_speed`
branch, lines 4510–4592. It accumulates inward integration-point area vectors
at boundary nodes, normalizes them to the prescribed speed, then interpolates
the nodal values to face samples. The segregated momentum boundary assembler
uses the nodal values separately in its viscous term.

For a planar Tri3 face, every sample area vector is one third of the outward
face area vector. With inlet speed U [m/s], form at node i:

```
s_i = -sum(incident inlet sample area vectors)       [m²]
u_bc,i = U * s_i / |s_i|                            [m/s]
u_bc,f,s = sum_k N_s,k * u_bc,node(f,k)              [m/s]
mass_flux_f,s = rho * dot(u_bc,f,s, A_f,s)           [kg/s]
```

Only inlet faces contribute to the sum, including at inlet/wall intersections.
Standard face weights are 11/18 at the nearest node and 7/36 at the other two.
Shifted `linear-linear` sampling takes the nearest node's value. The velocity
override in the inlet viscous assembly uses `u_bc,i`, not the interpolated sample.
Zero or nonfinite accumulated inlet normals cause collective setup rejection.

The previous implementation used `-U*A/|A|` independently on every face. This
agrees on a planar inlet but generally differs on a nonplanar inlet. Interpolated
sample speed need not equal U when neighbouring nodal directions differ; neither
implementation promises total inflow `rho*U*geometric_area` on a bent inlet.
This change does not rescale flux to impose a different total flow rate.

The steady incompressible equations, physical pressure, density, viscosity,
quadrature geometry, relaxation, initialization and linear targets are unchanged.
There is no physical time discretization in this driver. The existing flux sign
and continuity assembly are retained: the outward inlet mass flux is negative.

## Device and MPI implementation

`SimpleInletArea` sums all held inlet-face contributions into existing device
scratch. Complete owned element stars supply each owned node's full inlet sum.
`SimpleInletVelocity` normalizes owned nodes, then one existing node-halo exchange
copies the three components to ghost nodes. Partial sums at ghost nodes are never
normalized for use. The one-rank reference runner uses the same construction.

A persistent 3-component device field is reused for every iteration. No mesh or
field data is copied to the host by this setup; its error flag uses the existing
collective check. No new per-iteration communication is added. Lifetime exchange
counts increase by one startup round. This scope covers the supported stationary
mesh with one inlet set and a constant prescribed speed, not moving boundaries.

## Validation

Local checks on 2026-09-30 passed with AppleClang 21, OpenMPI 5.0.5 and
AddressSanitizer/UndefinedBehaviorSanitizer:

- 2,104 independent SIMPLE algebra checks. An analytic two-face tetrahedron
  checks area-weighted nodal directions, both interpolation choices, the distinct
  nodal viscous values and sample fluxes, and invalid-normal rejection.
- 26 focused host tests, including bent-inlet three-iteration parity on 1/2/4
  ranks, configured planar and shifted regressions, and outlet reversal. The bent
  fixture shears x by `0.2*y² + 0.1*z²`, then rotates it; partitioning uses logical
  coordinates. The gate compares nodal inlet data, matrix rows, RHS and subsequent
  fields/fluxes with the one-rank reference. These host solves use test-only LU.
- 50 synthetic preparation/snapshot-comparison tests, including explicit rebuilt
  executable identity and continued rejection of mismatched deck/controls.

The formula is matched by source inspection and analytic tests. These checks are
not an independent OpenAccel run on a bent inlet. CUDA/Hypre execution for this
change and its effect on private field disagreement remain pending. Earlier
public channel and duct results do not establish private pump parity.

## Daint public GPU checks

In Bash, from the existing configured CUDA/Hypre build. All files are public
synthetic data. Generate new reference snapshots because the gate now records
the prescribed nodal inlet field as well. No OpenAccel rebuild is required.

```bash
(
set -euo pipefail
cd /capstor/scratch/cscs/gandanie/git/mars-v010-check/build-hypre
git pull --ff-only
cmake --build . --parallel 4 --target mars_segregated_simple \
  mars_segregated_simple_algebra_check mars_distributed_simple_cuda_gate
./examples/distributed/unstructured/mars_segregated_simple_algebra_check
run=$(mktemp -d /capstor/scratch/cscs/gandanie/simple-inlet-XXXXXX)
printf 'Results: %s\n' "$run"
git rev-parse HEAD > "$run/mars-revision.txt"
sha256sum ./mars_distributed_simple_cuda_gate \
  ./examples/distributed/unstructured/mars_segregated_simple > "$run/executables.sha256"
export MARS_HYPRE_SPMV_VENDOR=0
unset CUDA_LAUNCH_BLOCKING
for interpolation in trilinear linear-linear; do
  args=(--mesh 8x2x2 --iterations 3 --bent-inlet 1 --configured 1
        --high-resolution 1 --velocity-interpolation "$interpolation")
  srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=1 \
    --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
    "$HOME/affinity/bind_numa.sh" ./mars_distributed_simple_cuda_gate \
    "${args[@]}" --write-reference "$run/$interpolation.bin" \
    2>&1 | tee "$run/$interpolation-reference.log"
  for np in 1 2 4; do
    srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node="$np" \
      --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
      "$HOME/affinity/bind_numa.sh" ./mars_distributed_simple_cuda_gate \
      "${args[@]}" --reference "$run/$interpolation.bin" --builder 1 \
      2>&1 | tee "$run/$interpolation-$np.log"
  done
done
)
```

After these pass, rerun the [early saved-state comparison](snapshot_comparison.md)
with `--allow-executable-change` in the preparation command. This records both
binary hashes and preserves the original controls; never overwrite the baseline
hash to make a rebuilt binary appear identical. Only the public comparison JSON
may be shared. A changed error band tests relevance to the private case; it does
not alone establish full numerical parity or a converged solution.
