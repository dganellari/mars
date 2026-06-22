# Handoff to the pump session — shared-code changes from the matrix-free / scaling session (2026-06-22)

Two changes this session touch code the pump path shares. Both need a pump-side check.
Ownership/DOF logic is **unchanged**; nothing here alters owned-DOF counts, BCs, or the pressure operator's owned rows.

## 1. `Tet4CVFEM::jacobian_and_dNdx` — I implemented a function your code CALLS (build was failing)

Your code calls `Tet4CVFEM::jacobian_and_dNdx<RealType>(coords, det, dNdx)` but `Tet4CVFEM`
had **no such member** -> the build failed (`class "mars::Tet4CVFEM" has no member
"jacobian_and_dNdx"` + a `variantName()` too-few-args error). I added the missing definition.

Callers — note these are on the **pump pressure-projection critical path**, not just tets:
- `mars_cvfem_tet_kernel_{graph,full,full_perip}.hpp`  (tet assembly)
- `mars_fem_projection.hpp`  (x4)   <- pressure projection
- `mars_vms_pressure_stab.hpp` (x2) <- VMS pressure stabilization

What I added (in `backend/distributed/unstructured/fem/mars_cvfem_kernel.hpp`): the standard
linear-tet form -- J columns = edge vectors from node 0 -> 3x3 adjugate inverse -> constant
`dNdx[n][i] = sum_j dNref[n][j] * (J^-1)[j][i]`, with `dNref = {(-1,-1,-1),(1,0,0),(0,1,0),(0,0,1)}`.
`det` is signed (= 6*volume); callers take `fabs(det)` for volume.

**CHECK (important):** the stiffness use (`dNdx . dNdx`) is sign-insensitive, but the
**projection / divergence / VMS uses are NOT**. Confirm the `dNdx` sign + index convention
matches what your projection and pressure-stab expect, then run your pump validation. It is the
textbook formula and I believe it is correct, but it is your physics on the critical path. If you
were about to write your own `jacobian_and_dNdx`, reconcile against mine (don't duplicate).

Also fixed: `examples/.../mars_cvfem_graph_tet.cu` `variantName()` -> `variantName(Assembler::Variant::GraphLump)`.

## 2. Shared node halo — receiver-driven send symmetrization (Option A)

`backend/distributed/unstructured/domain.cu` `buildFromCstoneHalos` (v2 path, `MARS_NODEHALO_V2=1`):
SEND lists are now **receiver-driven** -- each rank rebuilds its send list = exactly the global
keys peers request, so `A.send[B] == B.recv[A]` by construction. Fixes an `MPI_ERR_TRUNCATE`
over-claim that crashed the matrix-free solve at 64M DOF/GPU (verified fixed: gates pass at 1e-18,
64M clean, comm 19.2%). Build-time only, no per-matvec cost. **Ownership untouched.**

`domain.hpp`: added `waitallDiag` (`MARS_HALO_DEBUG=1`) -- turns an opaque `PMPI_Waitall` abort
into the exact per-request MPI error (truncation vs transport). Default-off. (A speculative
chunking experiment was added then removed -- net no change besides this instrument.)

The pump uses this exact `exchangeNodeHalo` / `reverseExchangeNodeHaloAdd`
(`mars_pump.cu`, `mars_ns_pump_solver.hpp`, `mars_ns_solver.hpp`).

**CHECK:** re-run pump baselines and confirm no regression -- values should be identical, or more
correct (the over-claim is removed). If the pump had latent multi-rank seam issues, the symmetric
topology may help.

## Regression oracle — NOT TGV
TGV multi-rank does not currently run, so it tells us nothing about regression. Use:
- the matrix-free **cube gates** (already PASS post-fix: A*1=0, A*linear=0 @ 1e-18, 64M/GPU clean), and
- **your pump baselines**: Rhie-Chow cavity Re=100 (max err), open-channel Poiseuille (dp), ACM-pump
  solve (residual / mass error) -- match against pre-change numbers.

## If a halo abort appears at scale
Rerun with `MARS_HALO_DEBUG=1` -> the `[halo-dbg] ... err="..."` line names the cause (truncation =
topology, CXI/GPU = transport) instead of an opaque crash.
