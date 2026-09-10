# Outlet CUDA/MPI validation

Author: GPT/Codex. Date: 2026-09-10.
Status: prepared; pull on Daint, build, and execution pending.

## Scope

One agent. Medium effort for source checks and preparing the existing gates.
The targeted source inspection of `assemble_outlet_continuity` found local-element
scattering, uniquely owned opening facets, and an unconditional reverse-add call.
`reverseExchangeNodeHaloAdd` packs ghost contributions, exchanges them with peers,
and adds them to owned nodes. This is source inspection, not distributed validation.

Distribution is through a commit pushed to `origin/cstone`; the user pulls it into
the Daint source checkout used by the current CUDA build. No rsync is needed.

## Build prerequisites

After pulling `cstone`, use the existing CUDA build and its existing dependency configuration. The three
targets are `mars_pump`, `mars_outlet_boundary_gate`, and `mars_outlet_kernel_gate`.
Reconfigure that build to pick up the two new targets. The example directory requires
`MARS_ENABLE_CUDA=ON`, `MARS_ENABLE_UNSTRUCTURED=ON`, and
`MARS_ENABLE_FEM_EXAMPLES=ON`. The outlet solver also requires Hypre support.
The user drives the build and launches below; no remote build was executed by GPT.

## Daint commands

Run from that build directory after all three targets compile. These gates use only
public synthetic fixtures. They do not load a mesh or write flow fields. Both gates
use the launcher's visible GPU; use the established NUMA/GPU binding wrapper.
Do not proceed past a failed command. Each launch has a five-minute scheduler limit.

```bash
srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=1 --export=ALL --kill-on-bad-exit=1 \
  ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_outlet_boundary_gate

srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=1 --export=ALL --kill-on-bad-exit=1 \
  ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_outlet_kernel_gate

srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=2 --export=ALL --kill-on-bad-exit=1 \
  ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_outlet_kernel_gate

srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=4 --export=ALL --kill-on-bad-exit=1 \
  ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_outlet_kernel_gate
```

## Acceptance and next boundary

Require zero exit status and PASS from every launch. The scalar gate runs its host
checks and then the production arithmetic on CUDA. The kernel gate checks actual
interior/boundary scatters, reporting, pressure derivatives, correction JVPs, and
anchor fault injection. Negative anchor fixtures must be detected without failing
the gate itself. On multiple ranks the injected failure originates only on rank 0;
empty-facet ranks must receive it through the reduction.

These tests use replicated tiny fixtures and real MPI reductions, not `ElementDomain`
communication. Passing them permits the next validation stage: a public synthetic
mesh through the real domain halo and full Hypre stepper, including BDF startup and
partition-independent continuity/boundary balance. Manufactured-flow accuracy and
refinement checks remain separate. No physical-case correctness claim follows from
these kernel gates alone.
