# marsir-mlir — the `mir` MLIR dialect (MARSIR Stage 2)

The **real MLIR** path for MARSIR: an out-of-tree MLIR dialect that lowers
high-level FEM operators down to the GPU / tensor-core level, i.e.

```
mir  →  linalg  →  gpu / nvgpu (tensor-core MMA)  →  nvvm  →  PTX
```

This is distinct from `marsir-compiler/` (the Python source-to-source emitter).
That one emits CUDA/HIP *text*; this one is a genuine compiler dialect with
lowering passes. Requires LLVM/MLIR (here: Homebrew `llvm` = MLIR 19.1.7).
The full walkthrough, with the why behind every mechanism, is
`internal-notes/marsir_full_tutorial.md`.

## Build

```bash
cmake -G Ninja -B build \
  -DMLIR_DIR=$(brew --prefix llvm)/lib/cmake/mlir \
  -DLLVM_DIR=$(brew --prefix llvm)/lib/cmake/llvm \
  -DCMAKE_BUILD_TYPE=Release
ninja -C build mir-opt
```

## The dialect (`include/mir/MirOps.td`)

- `mir.contract %u, %M {axis}` — the 1D sum-factorization contraction atom
  (Btil/Dtil/D/W along one axis).
- `mir.flux ins(...) { region }` + `mir.yield` — the authored pointwise flux
  (the D of B^T·D·B), region-carrying.
- `mir.simplex_contract` — the ragged collapsed-tet sweep
  (`r <= D - p - q` bounds, 4-D per-(p,q) factor table).
- `!mir.block<shape, p, deformed, batch>` type + `mir.bind` / `mir.unbind`
  with an `unbind(bind(x)) -> x` folder.

## Passes (all in `mir-opt`)

| pass | file | what it does |
|---|---|---|
| `--convert-mir-to-linalg` | `lib/LowerMirToLinalg.cpp` | contract -> `linalg.matmul` (axis-0 unfold) or `linalg.generic` reduction; flux -> all-parallel generic; simplex -> `scf.for` + iter_args |
| `--mir-hoist-transfer-pairs` | `lib/HoistTransferPairs.cpp` | loop-invariant transfer_read/write pair -> loop-carried iter_arg (the v4 accumulator registerization; upstream hoisting declines on aliasing); `assume-zero-init` option |
| `--mir-lower-copies` | `lib/LowerCopies.cpp` | `memref.copy` -> plain `scf.for` loops (strided copies otherwise lower to a memrefCopy host-runtime call, illegal in GPU kernels) |
| `--mir-warp-wrap` | `lib/WarpWrap.cpp` | wraps innermost `scf.parallel` bodies in `vector.warp_execute_on_lane_0` (lane = `gpu.lane_id`) + hoists uniform code out so boundary writes can distribute |
| `--mir-warp-distribute` | `lib/WarpDistribute.cpp` | lane-distributes warp regions via upstream VectorDistribution patterns (shipped in MLIR 19 but unexposed) + `WarpOpContractToMma`: bridges m8n8k4 f64 `vector.contract` to `nvgpu.mma.sync` on per-lane DMMA fragments |

## Pipelines / PTX scripts

All emit **real sm_90 PTX via LLVM's NVPTX backend** — no CUDA toolkit
anywhere; a Mac laptop produces deployable PTX. `test/extract_ptx.py` pulls
the PTX out of the `gpu.binary` attribute and asserts DMMA instructions are
present.

| script | what it builds |
|---|---|
| `run_dmma_pipeline.sh` | `mir.contract` -> `nvvm.mma.sync` end to end (`-q` counts the mmas) |
| `make_ptx.sh` | v1 sweep kernels -> `generated/sumfac_sm90.ptx` + `sumfac_batched_sm90.ptx` |
| `make_ptx_v4.sh` | the production sweep schedule: tile+vectorize -> `--mir-hoist-transfer-pairs` -> convert-vector-to-gpu *before* outlining -> `sumfac_batched_v4_sm90.ptx` |
| `make_ptx_v5.sh` | v5 experiment (zero-init + 8 warps/block; measured no-gain) |
| `make_ptx_full.sh [p]` | the fused full-operator kernel (whole Knaus apply, one kernel instead of 111 launches) -> `full_apply_p{p}_sm90.ptx` |
| `make_ptx_warp.sh` | the warp->tensor-core bridge gate -> `warp_mma_sm90.ptx` |

## Gates

**CPU execution** (`python3 test/exec_gate.py`, needs numpy): lowers mir
through the same `--convert-mir-to-linalg` the GPU path uses, executes with
`mlir-cpu-runner`, compares against NumPy *in-IR* (a second `mir.flux`
computes `|y−exp|·1e12`, immune to print precision). Six gates, all PASS:
contract==einsum, PA chain with the real CVFEM flux, `.op -> mir -> executed`,
ragged simplex (bit-exact), the **full HO operator vs a NumPy Knaus oracle**
(9.4e-16), and the **fused batched operator** at E=4 (1.8e-15).

**GPU** (`tools/`, driver-API-only harnesses: dlopen `libcuda.so.1`, `_v2`
symbols, build with `g++ -O2 <file>.cpp -ldl` on bare compute nodes):

- `run_ptx_gate.cpp` — single-warp sweep: **bit-exact on GH200**.
- `run_ptx_batched.cpp` — batched sweep + throughput
  (`<ptx> <E> [kernel] [warps/block]`; v4 kernel is `sweep_kernel`).
- `run_ptx_full.cpp` — fused full operator vs a built-in CPU Knaus oracle
  (`<ptx> <E> <p> [tpb]`).
- `run_warp_mma.cpp` — the warp->mma bridge: checks the per-lane fragment
  indexing against hardware (GH200 run pending).

## Measured (GH200)

Sweep tuning at E = 1M, all bit-exact — schedule changes only, zero dialect
changes: v1 23.9 ms (~539 GB/s, 13% of HBM peak) -> v3 elems/block packing
17.8 ms (~723 GB/s) -> **v4 accumulator hoist 4.404 ms (~2926 GB/s, ~73% of
peak)**. v5 showed no further gain (latency-bound; the toy sweep is tuned
out).

Fused full operator (scalar thread-per-element baseline): p=3 E=1M 34.8 ms
(~135 GB/s); p=7 E=131k 77.7 ms (~68 GB/s, spills). The ~30-40x gap to the
sweep ceiling is the quantified motivation for the warp/tensor-core work.

## Remaining (honest)

- GH200-run the warp->mma PTX gate (`run_warp_mma`); lowering to
  `nvvm.mma.sync {row,col,m8n8k4}` is verified, hardware is not.
- Integrate `--mir-warp-wrap` + `--mir-warp-distribute` into the full-operator
  payload; measure against the 2926 GB/s ceiling; smem-staging fallback for
  non-mma contract shapes.
- edof gather/scatter, then the head-to-head vs `mars_ho_dmma_bench`.
- Remaining simplex sweeps (q/p + transposes) and a ragged scheduling story.
- An MLIR 19 install on Alps for native iteration (PTX-only works today).
