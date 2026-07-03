# mirage-mlir — the `mir` MLIR dialect (Stage 2)

The **real MLIR** path for Mirage: an out-of-tree MLIR dialect that lowers
high-level FEM operators down to the GPU / tensor-core level, i.e.

```
mir  →  linalg  →  gpu / nvgpu(tensor-core MMA)  →  nvvm  →  PTX
```

This is distinct from `mirage-compiler/` (the Python source-to-source emitter).
That one emits CUDA/HIP *text*; this one is a genuine compiler dialect with
lowering passes. Requires LLVM/MLIR (here: Homebrew `llvm` = MLIR 19).

## Build

```bash
cmake -G Ninja -B build \
  -DMLIR_DIR=$(brew --prefix llvm)/lib/cmake/mlir \
  -DLLVM_DIR=$(brew --prefix llvm)/lib/cmake/llvm \
  -DCMAKE_BUILD_TYPE=Release
ninja -C build mir-opt
```

## What exists

- **Dialect** `mir` (`include/mir/MirOps.td`) with the op
  `mir.contract %input, %op_matrix {axis} : (tensor, tensor) -> tensor` — the 1D
  reference-operator contraction that is the atom of tensor-product
  sum-factorization (Btil/Dtil/D/W applied along one axis).
- **Tool** `mir-opt` — like `mlir-opt`, with the `mir` dialect + its passes.
- **Lowering** `--convert-mir-to-linalg` (`lib/LowerMirToLinalg.cpp`) — rewrites
  `mir.contract` to a `linalg.generic` reduction. `linalg` is MLIR's on-ramp to
  GPU + tensor cores.

## Try it

```bash
OPT=./build/tools/mir-opt/mir-opt

# round-trip the dialect
$OPT test/contract.mlir

# the bridge: mir.contract -> linalg.generic
$OPT test/contract.mlir --convert-mir-to-linalg

# end-to-end to the GPU/CUDA level
$OPT test/contract.mlir --convert-mir-to-linalg \
  --one-shot-bufferize=bufferize-function-boundaries=true \
  --convert-linalg-to-parallel-loops --gpu-map-parallel-loops \
  --convert-parallel-loops-to-gpu --gpu-kernel-outlining
# -> gpu.launch_func / gpu.module / gpu.thread_id ...
```

## Tensor cores (FP64 DMMA) — working

`./run_dmma_pipeline.sh -q` runs the whole thing: `mir.contract` →
`--convert-mir-to-linalg` (emits the *unfolded matmul* form:
`tensor.collapse_shape` + `linalg.matmul`) → bufferize (identity layouts —
required, or the nvgpu patterns silently decline on dynamic strides) →
`test/dmma_schedule.mlir` (transform-dialect schedule: tile n×8, k×4 = the
m8n8k4 DMMA tile, then vectorize) → `convert-vector-to-gpu{use-nvgpu}` (per-lane
warp fragments + `nvgpu.mma.sync`) → `convert-nvgpu-to-nvvm` →
**`nvvm.mma.sync`**, the FP64 tensor-core instruction MARS hand-writes in
`mars_ho_laplacian_dmma.hpp`. Note: the *contractions* hit tensor cores;
`mir.flux` is elementwise and stays on the FMA pipes — that split is correct.

## Numerical execution gates — passing

`python3 test/exec_gate.py` lowers mir through the **CPU pipeline** (same
`--convert-mir-to-linalg`, incl. the matmul unfolding the DMMA path uses),
executes with `mlir-cpu-runner`, and checks against NumPy **in-IR** (a second
`mir.flux` computes `|y−expected|·10¹²`, immune to print precision):

- gate 1: `mir.contract` == einsum — machine precision
- gate 2: PA chain with the **real CVFEM flux** (`g2·deriv + g0·dt2 + g1·dt1`,
  3 directional contractions) — machine precision
- gate 3: **`specs/laplacian.op` → Python front-end (`--backend mlir`) → mir
  dialect → executed** == NumPy — machine precision. One `.op` feeds Stage 1
  (CUDA source) and Stage 2 (this dialect).

## PTX artifact

`./make_ptx.sh` serializes **real sm_90 PTX** via LLVM's NVPTX backend (no CUDA
toolkit needed) → `generated/sumfac_sm90.ptx`: 16× `mma.sync.aligned.m8n8k4.
row.col.f64.f64.f64.f64` with accumulator chaining — the FP64 DMMA instruction
MARS hand-writes, now compiler-generated.

## Remaining (honest)

- Run the PTX on an actual GPU (Alps) + measure performance — unmeasured.
- Ragged (tet, `p+q+r≤D`) contraction variant of `mir.contract`.
- FEM types: `!mir.block` as a real MLIR type.
- Building mir-opt on Alps needs an MLIR 19 install there (spack/uenv llvm+mlir).
