# Native block assembly on the public channel

This gate combines native geometry and nodal reconstruction with the validated
interior, boundary and node terms. It assembles complete momentum and pressure
matrices/RHS for two captured iterations. Momentum influence coefficients are
computed from the assembled, relaxed matrix and then consumed by pressure
assembly. No linear solve or complete SIMPLE iteration runs yet.

The fixed public profile has 425 nodes and 1,536 Tet4 elements, steady laminar
incompressible flow, upwind advection, a prescribed normal inlet velocity,
no-slip walls and an average-pressure outlet. It uses physical pressure,
density and dynamic viscosity. It is not the private pump mesh.

## Discrete contract

For each stage, mapped local contributions add into the same node graph:

```
A = sum(element matrices) + sum(boundary matrices) + sum(node matrices)
b = sum(element RHS)      + sum(boundary RHS)      + sum(node RHS)
```

Momentum stores full 3x3 blocks, including the cross-component stress terms.
Its unknown is the velocity increment, not absolute velocity. Matrix units
are kg/s; RHS units are force (N). Only each component's scalar diagonal is
scaled by 1/alpha_u. There is no additional absolute-velocity relaxation RHS.
The influence coefficient is `d_i = V_i/(A_ii + SMALL)` after relaxation;
its units are m^3 s/kg. The fixed profile uses SIMPLE, not SIMPLEC.

Next, multiply momentum RHS values at the union of exterior nodes by 0.75,
once per node even at corners. This follows `navierStokesAssembler::postAssemble`:
base constraints/diagonal relaxation, influence coefficients, boundary RHS
relaxation, then symmetry (absent in this profile). The gate checks the selected
node union against the actual captured boundary-relaxation records.

Pressure RHS is minus the net outward mass flux (kg/s). Pressure matrix units
are m s, multiplying an increment in Pa. Its interior contributions cancel by
paired face signs. The outlet derivative retains the face-row/opposite-column
coupling; boundary column multipliers are local term masks, not identity rows.
No velocity-style diagonal relaxation or arbitrary pressure pin is added.
At the pinned reference, `pressureCorrectionEquation::solve` only pins a row
when `domain::pressureLevelRequired()` is true; the public inlet/outlet profile
makes that false. It is incompressible, so the compressibility diagonal is also
absent. Interface constraints, symmetry and inactive domains are absent here.
Neither assembled operator is assumed symmetric; solver integration must
respect that rather than selecting CG from the word "pressure".

## Device implementation and ownership

`DeviceBlockGraph` creates all 16 node pairs per tetrahedron plus nodal diagonals,
then sorts and deduplicates packed row/column keys on-device. It builds node CSR
with 4,921 blocks for this fixture; momentum stores nine values per block.
Only Thrust's unique count returns to the host for allocation. Geometry,
connectivity, indices, coefficients and assembled arrays remain on-device.
The graph is built once. The gate reuses persistent arrays across all four
states and zeros matrix/RHS storage before each assembly.

`scatter_block` verifies all required entries, even when their coefficient is
zero, before adding any part of that local block. A missing entry sets the
checked device error flag. Node/component mapping preserves the full block;
three independent scalar velocity matrices are not substituted.

Callers must supply valid connectivity/geometry, owned element/node ranges and
unique owned boundary records. This gate covers one rank only. Distributed
matrix/RHS reverse-add, ghost columns, collective failure and halo scheduling
still need integration. Error-flag and full result downloads in this executable
are validation boundaries, not a production solver's hot-loop policy.

## What remains captured

The stage velocity/pressure fields are prescribed reference states, not outputs
from new MARS solves. Boundary trace/velocity substitutions, wall coefficients,
reversal flags, stored mass-flux history, nodal mass divergence and pseudo-time
are still captured inputs. Their evolution is not certified by this gate.
Geometry, dual volumes and gradients are recomputed. The momentum matrix
produces the coefficients used in pressure assembly; captured `du` is used
only as an expected answer, never fed into the numerical path.

Pressure is unchanged between momentum and pressure assembly in a given
reference iteration. The packer joins pressure from that iteration's pressure
capture and verifies its reconstructed gradient against the momentum-node
capture. Velocity states remain distinct before and after the reference
momentum solve.

The complete-matrix oracle sums actual captured local LHS/RHS contributions by
global node ID in Python, independently of the C++ scatter implementation.
Captured global momentum diagonal, boundary RHS and influence records provide
additional cross-checks. These are not independently exported complete global
OpenAccel CSR matrices or solved fields. Pinned deck/mesh/export hashes and
stage joins are mandatory; mixing arbitrary runs is rejected.

## Validation and next integration

Local strict C++ host builds pass 176 independent scatter/relaxation checks and
104,370 comparisons across the four systems. Maximum scaled matrix difference
is 1.97e-15; maximum RHS difference 1.39e-16; influence difference 2.17e-18.
Wrong graph sizes, truncated inputs and perturbed expected matrix values fail.
CUDA compilation/execution is pending. No convergence, MPI or speedup claim.

Next: prepare boundary state and stored-flux history natively, expand the node
block graph for the selected linear-solver interface, solve increments and
connect the existing ordered updates. Compare the complete iteration's fields
and fluxes before enabling an outer SIMPLE loop. Reuse or extend the existing
instrumented reference build only when new expected state is actually needed.

From the configured MARS CUDA build (`mars/mlir`), after publication:

```bash
(
set -euo pipefail
git pull --ff-only
assembly_work=$(mktemp -d "$PWD/segregated-assembly-XXXXXX")
python3 ../scripts/prepare_openaccel_assembly.py \
  $SCRATCH/git/OpenAccel-reference-updates-IxgJIp/run/exports \
  --boundary $SCRATCH/git/OpenAccel-reference-boundary-TljPP1/run/exports/boundary \
  --nodes $SCRATCH/git/OpenAccel-reference-nodes-20260922-171816-KT5eld/run/exports/nodes \
  --output "$assembly_work/inputs.txt"
cmake -S .. -B .
cmake --build . --target mars_segregated_assembly_check mars_segregated_assembly_algebra_check -j4
./examples/distributed/unstructured/mars_segregated_assembly_algebra_check
srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=1 \
  --export=ALL --kill-on-bad-exit=1 \
  ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_segregated_assembly_check \
  "$assembly_work/inputs.txt" 2>&1 | tee "$assembly_work/assembly.log"
git rev-parse HEAD > "$assembly_work/mars-revision.txt"
sha256sum ./examples/distributed/unstructured/mars_segregated_assembly_check \
  "$assembly_work/inputs.txt" > "$assembly_work/sha256.txt"
printf 'Saved: %s\n' "$assembly_work"
)
```

No new OpenAccel build or simulation is needed for this gate.
