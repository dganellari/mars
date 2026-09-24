# Steady momentum node replay

This slice ports four local operations from the pinned OpenAccel reference:
steady momentum node terms, diagonal velocity relaxation, SIMPLE/SIMPLEC influence
coefficients, and the arithmetic of boundary RHS relaxation. The existing MARS
projection solvers are unchanged. This is not a complete SIMPLE solver.

For fixed-frame, incompressible steady flow, the nodal unknown is a velocity
increment. Pressure is physical pressure (Pa), `V` is dual volume (m³), `rho` is
kg/m³, `dt` is the steady pseudo-timescale (s), and `div_m` is integrated mass
flux divergence (kg/s). The local terms are

```
A_ii += max(-div_m, 0) + rho*V/dt
b_i  += div_m*u_i - V*grad(p)_i + V*(force_i + source_i)
```

The diagonal has units kg/s and the RHS has units N. Pseudo-time adds no steady
RHS. `force` and the non-redistributed `source` are force densities (N/m³).
Inputs are already gathered nodal fields; volume construction and pressure
reconstruction are not reproduced by this gate. Rotating and transient cases
are rejected.

After reference constraint assembly, diagonal velocity relaxation divides only
`A_ii` by `alpha_u`. In this increment formulation it adds no RHS term. Using the
resulting diagonal:

```
d_i       = V / (A_ii + SMALL)
d_tilde_i = V / (A_ii + sum(j != node, A_ij for component i) + SMALL)
```

`SMALL` is the reference's FP64 epsilon. No second `alpha_u` belongs in these
expressions. SIMPLEC excludes the whole diagonal node block from the neighbor
sum and ignores cross-component entries. Both coefficients have units m³ s/kg.
The reference subsequently multiplies RHS components at its selected boundary
nodes by 0.75. Its union selector applies this once to each owned node; this
slice ports the multiplication, not the mesh selector. Influence values are
captured before halo publication and symmetry handling.

The reference patch exports actual inputs and outputs at all four sites, by
global node ID and call. `openaccel_node_check.py` rejects incomplete stages,
missing rank files (including empty ranks), repeated owned nodes, invalid
scales, unsupported modes, and a coefficient row whose diagonal differs from
the captured relaxed diagonal. It does not calculate an oracle. Expected
outputs are never copied to device input buffers. The CUDA replay calls the
production header, with one thread per record. Host transfers are confined to
this test harness.

## Run after review and publication

In the existing working **OpenAccel Cray-MPICH environment**, from the MARS
build directory:

```bash
bash ../tests/reference/openaccel/node_replay/capture.sh
```

Defaults use the sibling `OpenAccel` checkout, its `prgenv` build cache, and the
already verified `mars-reference-inputs-20260920-v2` public channel bundle.
Override those three paths as positional arguments if necessary. The script
creates a fresh instrumented worktree, rebuilds only OpenAccel with its recorded
dependencies, and launches the same two-iteration public channel. It installs
no dependencies and does not alter the existing executable. It prints the exact
node replay command with the new export path. Output must say
`node_capture_completed`; two iterations establish neither convergence nor parity.

Then restore the existing **MARS CUDA build environment**, and from that build:

```bash
cmake -S .. -B .
cmake --build . --target mars_segregated_node_replay mars_segregated_node_algebra_check -j4
# Use the exact export path printed by capture.sh:
bash ../tests/reference/openaccel/node_replay/run.sh /absolute/reference/run/exports/nodes
```

Both scripts launch one rank for five minutes and save exit status, logs, and
provenance. The replay records binary/source/input hashes. No new Spack
installation, hardcoded CUDA compiler, or MPI-family switch is required.

## Evidence and remaining scope

Independent algebra checks and host writer/transport replay use hand-worked
fixtures, including off-component sentinels and nonzero diagonal block offsets.
They are not an executed OpenAccel reference comparison. Actual instrumented STK
compilation and CUDA node parity remain pending. The already passed interior
CUDA gate is separate.

Still outside this slice: boundary face operators, gauges and elimination,
geometry/reconstruction, global CSR accumulation, distributed halo semantics,
flux/pressure/velocity correction and relaxation, a full SIMPLE iteration,
convergence, performance, and pump validation.
