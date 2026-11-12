# Nonconforming AMR validation

These opt-in checks locate failures in the public cube Galerkin test. They are not
a performance benchmark or acceptance of nonconforming AMR on general meshes.

`mars_amr_galerkin --audit-mesh --element=tet` checks owned parent counts, each
parent's emitted child volumes and orientation, canonical four-node cell identities,
and consistently decoded volume across redistribution and halo extension. The normal
solver path does not run these checks. The count statistic now includes owned marks
only. Snapping's volume change is reported separately from raw refinement.

`--mesh-only` skips the hanging-node, DOF and solver paths and transfers no fields.
It keeps the same geometric marking and refinement loop. `--audit-algebra` checks
owned-corner coverage, independent affine bases, the P/P-transpose adjoint relation,
assembled versus composed matrix action, symmetry, positive energy on a probe and
both true solution residuals. Algebra comparisons use double precision with a
relative tolerance of 1e-10; the true residual gate is ten times the requested solve
tolerance. Exact-field residuals are reported without treating quantized geometry
as an exact affine patch test.

The mesh audit is enabled internally by `MARS_AMR_AUDIT_TET=1`. Use it only for the
tet diagnostic. Full cell records stay on device unless a check fails and
`MARS_AMR_AUDIT_DIR` names a dump directory. Hash sums screen for changes; exact
before/after CSV records allow independent comparison. Mesh failures abort all ranks
with code 77, algebra failures with 78; a launcher may report a different nonzero
status after terminating sibling ranks.

## Accepted native snapshots

| Check | NVIDIA GH200 jobs | AMD MI300A jobs | Scope |
| --- | --- | --- | --- |
| Cornerstone rebuild and constrained Poisson | 5337613 | 737287, 737304 | 15 cases per job; tet/hex, 1/4/8 ranks, audited and repeated tet cases |
| Immediate fields and empty local input | 5338560 build/control + 5339346 continuation | 737668 | 12 cases and an eight-GPU safety fixture per platform; 36 adaptation rounds |

The mesh campaign preserved counts and both cell fingerprints in all 39 audited
rebuilds. All 18 unaudited eight-rank tet repetitions passed. Worst maximum
solution error was `2.414034e-6` against `1e-4`; assembled/composed operator error
was `8.991060e-16` against `1e-10`. Its kit manifest SHA256 is
`87560d129a8f1c2a12cd8bec6e0dc106f83e433a6e53d100d710b4328a0f3eb2`.

The immediate-field campaign checks old values bit-for-bit, including an arbitrary
keyed field, on owners and needed ghost nodes. Constant error is zero; maximum
affine error is `7.085505526194247e-7`, below the accumulated coordinate-quantization
bound (largest bound `1.5735643724037992e-6`). Nonempty multi-rank cases must move
ownership. Seven ranks begin with no elements in the empty-input cases. The
staged cases stage AMR messages, not all Cornerstone/DOF communication.
Its kit manifest SHA256 is
`831a752b315a72ecae9efb6b72eefbb6d9670476c08f2077852104c5009603d0`.

Daint 5338560 built successfully and reproduced the intended old-code field
failure on all eight ranks. Its launcher misclassified Slurm's exit 143 after
`MPI_Abort(82)` as an unrelated failure. The failed receipt is retained. The
run-only continuation 5339346 pins the original binaries, all 565 source files,
42 original artifacts and runtime-library hashes; it runs the missing corrected
cases without rebuilding or repeating that control. Continuation manifest:
`3d47deec4388f965df507581a55de34193f0b69c3c9d37f54e7c8d550fe775a2`.

CPU references exercise the production sparse migration, direct/staged messages,
missing/duplicate source rejection and tree layout at 1/2/4/8 MPI ranks. Separate
CPU Cornerstone sync/addHalos/property exchange checks demonstrate changed owner
coverage. These are component checks, not execution of the CUDA/HIP manager on CPU.

The numerical snapshots above are accepted. Integration with newer source and
relocated installation/package checks have separate provenance; the old binaries
do not validate later packaging changes. No scaling or speedup claim follows.

## Maintained checks

Configure `tests/amr` with `MARS_AMR_MPI_CHECKS=ON` and the Cornerstone include
path for CPU checks. Labels include `field-transfer`, `tree-layout`, `tet-lepp`
and `diagnostic`.

For an installed CUDA/HIP MARS package, enable `MARS_AMR_GPU_CHECKS`, provide
`MARS_AMR_TRANSFER_MESH_DIR` containing the public `cube16.mesh` and
`cube16_tet.mesh` fixtures, and set `MARS_AMR_TEST_RANKS` to the allocated ranks.
The field-transfer label covers both marking policies, staged AMR messages and
empty local inputs; safety also exercises key-buffer initialization and empty
refiners. Follow [the AMR launch guidance](../README.md) for rank/GPU binding.

The scripts in this directory reproduce the earlier CUDA diagnostic campaign.
They are not a current release launcher: do not resubmit an old sealed kit to
validate a newer source revision. Retain its failed controls and pinned inputs.

For redistribution/halo failures, run `compare_cells.py CASE/cells RANKS` on the
retrieved results. It refuses incomplete rank inventories and reports exact missing
and added canonical cells. Do not interpret the parent-versus-child dump as an
identity-preserving comparison when the failing stage is raw refinement.
