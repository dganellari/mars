#!/bin/sh
# PTX for the warp-distribution -> tensor-core bridge gates:
#   test/warp_mma_kernel.mlir  standard  C = A*B + C       (B stored [k,n])
#   test/warp_mmt_kernel.mlir  transp-B  C = In @ W^T + C   (B stored [n,k], axis-1)
# Each: contract(s) inside vector.warp_execute_on_lane_0
#   -> --mir-warp-distribute (WarpOpContractToMma: per-lane fragments + nvgpu.mma.sync)
#   -> NVVM -> PTX (sm_90, LLVM NVPTX backend, no CUDA toolkit).
# Run with tools/run_warp_mma / tools/run_warp_mmt on a GPU node.
set -e
export PATH="/opt/homebrew/opt/llvm/bin:$PATH"
mkdir -p generated

emit() {
  # distribute -> gpu-lower-to-nvvm-pipeline DIRECTLY. The pipeline does the
  # nvgpu/vector + memory-space lowering itself; a manual gpu.module(convert-
  # nvgpu-to-nvvm,convert-vector-to-llvm) pre-stage ERRORS on workgroup (shared)
  # memory ("space conversion failed"). Clean path handles both global and
  # workgroup-memory kernels.
  build/tools/mir-opt/mir-opt "$1" \
      --mir-warp-distribute --canonicalize --cse --lower-affine \
  | mlir-opt --gpu-lower-to-nvvm-pipeline="cubin-chip=sm_90 cubin-format=isa" \
  | python3 test/extract_ptx.py "$2"
  grep -c "mma.sync.aligned.m8n8k4" "$2"
}

emit test/warp_mma_kernel.mlir generated/warp_mma_sm90.ptx
emit test/warp_mmt_kernel.mlir generated/warp_mmt_sm90.ptx
# one full direction of the HO operator (D + flux + W, both arms, staged): 6 mma
emit test/warp_one_dir.mlir    generated/warp_one_dir_sm90.ptx
# +/- plane scatter, 2 faces, loop-carried overlap (census's hardest blocker): 2 mma
emit test/warp_scatter_kernel.mlir generated/warp_scatter_sm90.ptx
# EMITTER-GENERATED (mlir_warp.py): std + transposed contract staged in .shared: 4 mma
python3 ../marsir-compiler/marsir/backends/mlir_warp.py > test/warp_emit_selftest.mlir
emit test/warp_emit_selftest.mlir generated/warp_emit_selftest_sm90.ptx
# EMITTER-GENERATED B-sweep: Btil(8x8) @ U(8x64) column-tiled: 16 mma
python3 ../marsir-compiler/marsir/backends/mlir_warp.py matmul > test/warp_bsweep.mlir
emit test/warp_bsweep.mlir generated/warp_bsweep_sm90.ptx
# EMITTER-GENERATED single face: 4 contracts + 3-term flux, staged in .shared: 8 mma
python3 ../marsir-compiler/marsir/backends/mlir_warp.py face > test/warp_face.mlir
emit test/warp_face.mlir generated/warp_face_sm90.ptx
# EMITTER-GENERATED full single direction (dir0): B-sweep + 7 faces + scatter: 88 mma
python3 ../marsir-compiler/marsir/backends/mlir_warp.py dir0 > test/warp_dir0.mlir
emit test/warp_dir0.mlir generated/warp_dir0_sm90.ptx
