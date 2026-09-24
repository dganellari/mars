# Native geometry and nodal reconstruction

This is the first integration layer for the separate MARS SIMPLE solver.
Coordinates, connectivity and nodal field values are inputs; the device now
computes geometry, dual volumes and complete nodal gradients. Reference areas,
derivatives and gradients remain on the host as expected answers. The gate uses
existing public captures and needs no new OpenAccel build or simulation.

## Discrete contract

The selected profile is steady incompressible laminar flow, physical pressure
in Pa and dynamic viscosity in Pa s. This layer supplies geometry and gradients
to the already validated local momentum, continuity and update calculations.
For a positively oriented affine Tet4 of volume V and derivatives b_n=grad N_n:

```
local dual volume = V/4
interior A_LR     = (V/4) (b_R - b_L)
boundary A_sample = -V b_opposite
```

The second formula is the affine image of the median-dual quadrilateral.
The third gives one third of the outward triangular face area. Areas are m²,
volumes m³ and shape derivatives 1/m. Inverted/degenerate or nonfinite geometry
fails; nodes are never silently reordered.

Interior unshifted weights are 13/36 at edge endpoints and 5/36 elsewhere.
Shifted weights are 1/2 at each endpoint. Boundary unshifted weights are 11/18
at the nearest node and 7/36 elsewhere; shifted weights select that node.
These match `TetSCS` and `Tri3DSCS` at the pinned OpenAccel revision.
`Tet4CVFEM::scs_coords` is not used: its sampling table is different.

For each field component, the incremental reconstruction is

```
G(q)_i = [sum_interior sigma_is (I_s(q)-q_i) A_s
          + sum_boundary_at_i (I_b(q)-q_i) A_b] / V_i
```

Boundary samples interpolate the volume nodal field, not the independent
pressure trace. Pressure reconstruction uses shifted weights; velocity
reconstruction uses unshifted weights. Do not confuse velocity's interpolation
with its separately selected compact derivative interpolation. The constant
Tet4 derivative is shared by shifted/unshifted compact differentiation.
This matches `nodeField::updateGradientField` interior and ordinary boundary
terms, without symmetry, interfaces or limiting. In particular, shifted
boundary reconstruction need not reproduce affine fields exactly.

First reconstruction uses G(q) directly. Later calls may use
`(1-eta)*g_old+eta*G(q)`; eta=1 in this gate and for the pressure increment.
The first/unrelaxed paths avoid reading uninitialized gradient history.

## Device ownership and integration plan

`TetMeshView<KeyType,RealType>` accepts the four SoA arrays from
`ElementDomain::getElementToNodeConnectivity()` and the existing `getNodeX/Y/Z()`
arrays. It does not require a host numbering or geometry conversion. The gate
loads its public fixture on the host once; a future mesh-backed driver supplies
the view directly. No mesh-backed solver/controller is claimed by this gate.

Geometry and dual volumes persist for the mesh lifetime. Interior numerator
kernels scatter only the supplied owned-element range; boundary kernels consume
uniquely owned exterior facets. Caller preconditions: valid ranges, zeroed sum
buffers and successful geometry status before later kernels. For MPI, reverse-add
nodal contributions and publish completed sums before normalization. Distributed
ownership/halos are not exercised here. The caller checks launch/asynchronous
errors; only an error flag needs a host read in production. The gate deliberately
downloads computed values for validation.

The remaining single-iteration integration proceeds in dependency order:

1. Use these native geometry and reconstruction arrays with the existing
   `tet_interior`, `boundary_block` and `steady_momentum_node` functions.
   Construct sorted tetrahedral CSR on device and expand every 3x3 momentum
   block; retain the cross-component stress entries. Every expected CSR entry
   must exist. The old scalar projection sparsity is not the new block pattern.
2. Match boundary side/nodal-side preparation, wall coefficients and boundary
   RHS selection. These were inputs to the earlier gates and remain new work.
   Assemble, relax the momentum diagonal, then derive influence coefficients
   from the same matrix stage. Compare assembled matrix/RHS and influences
   before introducing linear-solver differences.
3. Solve the momentum increment, update velocity and its reconstruction, freeze
   the reference outlet state, and assemble/solve the pressure increment.
   Reuse `HypreGMRESSolver<double,int,cstone::GpuTag>` with device global maps;
   its full cross-component input uses scalar row `3*node+component`.
   Invalidate prepared setup whenever matrix values change. Verify the selected
   pressure reference policy against OpenAccel before inserting any gauge.
4. Apply the validated pressure/raw-increment, stored-flux/reversal and velocity
   updates in reference order. Reconstruct the raw increment separately from
   relaxed pressure. Compare each intermediate field and continuity using the
   same stored flux. Only then expose a complete SIMPLE iteration and outer loop.

Existing Chorin drivers and shared geometry helpers are unchanged. This plan
does not authorize another reference run; any additional assembled-state capture
should be specified once, then use the reusable instrumented build.

## Evidence and direct Daint gate

Local strict C++ builds pass 4,545 independent checks and 232,536 comparisons
against the actual saved public interior/boundary captures. Worst scaled errors:
element geometry 1.83e-17, boundary geometry 1.39e-17, velocity gradient 6.11e-16,
pressure gradient 4.02e-15. Actual CUDA compilation/execution is pending.
No full SIMPLE solve, convergence, MPI or performance result is implied.
The packer requires the recorded public deck/mesh hashes and complete capture
hashes, consistent shared nodes, closed boundary coverage and both field states.

After review/publication, from `mars/mlir` in the configured MARS CUDA environment:

```bash
(
set -euo pipefail
geometry_work=$(mktemp -d "$PWD/segregated-geometry-XXXXXX")
printf 'Geometry gate: %s\n' "$geometry_work"
python3 ../scripts/prepare_openaccel_geometry.py \
  $SCRATCH/git/OpenAccel-reference-updates-IxgJIp/run/exports \
  --boundary $SCRATCH/git/OpenAccel-reference-boundary-TljPP1/run/exports/boundary \
  --output "$geometry_work/inputs.txt"
cmake -S .. -B .
cmake --build . --target mars_segregated_geometry_check mars_segregated_geometry_algebra_check -j4
./examples/distributed/unstructured/mars_segregated_geometry_algebra_check
srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=1 \
  --export=ALL --kill-on-bad-exit=1 \
  ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_segregated_geometry_check \
  "$geometry_work/inputs.txt" 2>&1 | tee "$geometry_work/geometry.log"
git rev-parse HEAD > "$geometry_work/mars-revision.txt"
sha256sum ./examples/distributed/unstructured/mars_segregated_geometry_check \
  ../backend/distributed/unstructured/fem/segregated/mars_segregated_geometry.hpp \
  ../backend/distributed/unstructured/fem/segregated/mars_segregated_geometry_device.hpp \
  "$geometry_work/inputs.txt" > "$geometry_work/sha256.txt"
)
```
