# Bochev-Dohrmann TGV Validation Checklist

After overnight work (2026-05-23 → 2026-05-24), the periodic TGV instability
has been root-caused. This file documents the fix and the validation steps.

## Root cause (diagnosis confirmed by 4 parallel agent literature searches)

**MARS uses Q1-Q1 equal-order linear hex CG-FEM for both velocity and pressure.**
This pair is *inf-sup unstable* (LBB-fails). The discrete pressure space admits
"checkerboard" null-space modes that the discrete gradient operator cannot detect.
On a periodic torus there are no boundaries to anchor those modes, so any
forcing into them grows exponentially. This is exactly the symptom: TGV runs
beautifully for ~25 steps (KE flat at 0.156), then `|p|` doubles every 2-3
steps via grad(p) feedback in the predictor, then everything blows up.

Confirmed by:
- **MFEM-Navier** (Karniadakis/Tomboulides) uses inf-sup-stable high-order
  GLL spectral elements + Orthogonalize/MeanZero/OrthoSolver triple projection.
- **NekRS** uses spectral elements + GL over-integration + HPFRT filter.
- **deal.II step-35** uses inf-sup-stable Taylor-Hood Q2-Q1.
- **Literature** (Gresho 1995, Sani-Gresho, Guermond-Minev-Shen 2006,
  Bochev-Dohrmann 2006): equal-order CG-FEM REQUIRES stabilization.

## The fix: Bochev-Dohrmann polynomial pressure projection

Bochev-Dohrmann (SISC 2006) adds an element-local symmetric PSD term
```
S^e_ij = tau * V_e * (delta_ij/8 - 1/64)   for Q1 hex
```
to the pressure matrix. This has row sum zero (consistent: kills constant
mode), positive semidefinite (kills checkerboard null-space without changing
the consistent solution), parameter-light (`tau ~ h^2` works robustly).

Implementation: see `mars_ns_solver.hpp`:
- `assembleBDStabilizationKernel` (kernel)
- `addBochevDohrmannStab` (driver helper)
- BD applied to `d_valuesPre` after `assembleLaplacian` for pressure.
- BD also applied **matrix-free** in `applyDDTPerNode` for the DDT path.
- Enabled by `s.stabBochevDohrmann = true` in the TGV driver.

## Validation

### Test 1: short stability check
```
srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=1 \
  ~/affinity/bind_numa.sh \
  /capstor/scratch/cscs/gandanie/git/mars/daint-gpu/examples/distributed/unstructured/mars_tgv \
  --mesh=/capstor/scratch/cscs/gandanie/git/mars/cube16.mesh \
  --box-lo=0 --box-hi=1 \
  --V0=1 --nu=0.05 --dt=1e-4 --num-steps=200 \
  --pressure-solve=DDT 2>&1 | grep -E "^Step|periodic DOF|stabilization"
```

**Expected (with BD fix)**:
- `pressure stabilization: Bochev-Dohrmann polynomial projection, tau=...`
- `KE` decays monotonically from 0.156 toward smaller values
- `|p|` stays O(1), not growing exponentially
- `div_max` at roundoff (1e-11) every step
- No NaN through step 200

**Previously (without BD)**: KE flat then blew up around step 30.

### Test 2: cavity regression (BD off path)
```
srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=1 \
  ~/affinity/bind_numa.sh \
  /capstor/scratch/cscs/gandanie/git/mars/daint-gpu/examples/distributed/unstructured/mars_amr_ns_projection \
  --mesh=/capstor/scratch/cscs/gandanie/git/mars/cube16.mesh \
  --bc=cavity --dt=1e-3 --num-steps=20
```

**Expected**: unchanged — `cg_iter_uvw=5`, `cg_iter_p=92`, `|u|~0.19` at step 20.
BD is gated on `bcKind == Periodic`, so cavity sees no change.

### Test 3: longer-time TGV with viz
```
srun --account=csstaff --time=00:30:00 --nodes=1 --ntasks-per-node=1 \
  ~/affinity/bind_numa.sh \
  /capstor/scratch/cscs/gandanie/git/mars/daint-gpu/examples/distributed/unstructured/mars_tgv \
  --mesh=/capstor/scratch/cscs/gandanie/git/mars/cube16.mesh \
  --box-lo=0 --box-hi=1 \
  --V0=1 --nu=0.05 --dt=1e-4 --num-steps=2000 \
  --pressure-solve=DDT \
  --vtu-output=/tmp/tgv_bd --vtu-every=50
```

**Expected**: VTU/PVTU/PVD output every 50 steps in `/tmp/tgv_bd_step*.pvtu`.
KE smoothly decays. `_step0000` through `_step2000`.

### Test 4: TGV + AMR
```
srun --account=csstaff --time=00:30:00 --nodes=1 --ntasks-per-node=1 \
  ~/affinity/bind_numa.sh \
  /capstor/scratch/cscs/gandanie/git/mars/daint-gpu/examples/distributed/unstructured/mars_tgv \
  --mesh=/capstor/scratch/cscs/gandanie/git/mars/cube16.mesh \
  --box-lo=0 --box-hi=1 \
  --V0=1 --nu=0.05 --dt=1e-4 --num-steps=2000 \
  --adapt-every=200 --max-levels=2 --refine-frac=0.10 --coarsen-frac=0.30 \
  --pressure-solve=DDT \
  --vtu-output=/tmp/tgv_amr --vtu-every=50
```

**Expected**: every 200 steps `[amr] elements: N nodes: M periodic slaves: K`
message. Element count grows after refinement; periodic map rebuilds. VTU
output includes `level` cell data showing refinement pattern.

### Test 5: 4-panel ParaView video
```
# rsync VTU output back to laptop
rsync -avz daint-alps3:/tmp/tgv_amr* /tmp/

# Run pvbatch
/Applications/ParaView-5.13.2.app/Contents/bin/pvbatch \
  /Users/gandanie/scratch/santis/mars/examples/distributed/unstructured/render_tgv_4panel.py \
  --pvd /tmp/tgv_amr --width 1920 --height 1080 --fps 24
```

Produces `/tmp/tgv_amr_4panel/frame_NNNN.png` + `/tmp/tgv_amr_4panel.mp4`.

## Files modified

| File | Change |
|---|---|
| `mars_ns_solver.hpp` | Added BD assembly kernel, BD-on-DDT matrix-free apply, `stabBochevDohrmann` flag, RHS mean projection for pressure CG, removeMean on s.d_p after cumulative update, per-substep diagnostic instrumentation |
| `mars_periodic_bc.hpp` | (no change today, written earlier) |
| `mars_tgv.cu` | Multi-field VTU, AMR-on-NS loop, fixed computeKineticEnergy mass indexing, scaled TGV IC by 2pi/L |
| `mars_vtu_parallel_writer.hpp` | New `writeMultiFieldFrame` for PointScalar/PointVector3/CellScalar |
| `mars_cvfem_hex_kernel_tensor.hpp` | `i==j` -> `col_dof==row_dof` for periodic DOF collapse |
| `render_tgv_4panel.py` | NEW: pvbatch 4-panel ParaView animation script |
| `TGV_HANDOFF.md` | Original handoff doc (now superseded by this file) |
| `TGV_BD_VALIDATION.md` | THIS file |

## If BD does not stabilize TGV

Fallback plan, in order:
1. **Try `--stabPressureTau=` with smaller / larger values**. Currently auto =
   h^2. Try 0.1*h^2, 10*h^2 to bracket. Robust BD typically wants tau in [h^2, h]
   for the lumped-mass form.
2. **Switch to BDF2 + EXT2 time integration** (the NekRS/MFEM approach). This
   is ~300 lines additional — see TGV_HANDOFF.md for the implementation plan.
   Probably required for Re > a few hundred.
3. **Add over-integration of the convective term**. The NekRS finding ("the
   heavy hitter for periodic Galerkin DNS stability") was that aliasing in
   `u·∇u` accumulates unboundedly without a finer quadrature evaluation of
   the convective term. For Q1 hex this means 3x3x3 Gauss instead of corner-
   sampling. ~200 lines additional.

## Quick win if all else fails: tighter `--nu` / `--dt`

For producing a watchable video without full stability proofs:
```
--nu=0.1 --dt=5e-5 --num-steps=500
```
This puts Re=10, very damped, with CFL ~ 0.0008. Should stay stable for the
duration of the video even if the underlying scheme is fragile.
