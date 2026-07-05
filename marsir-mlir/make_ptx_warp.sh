#!/bin/sh
# PTX for the warp-distribution -> tensor-core bridge gate:
#   test/warp_mma_kernel.mlir (contract inside warp_execute_on_lane_0)
#   -> --mir-warp-distribute (WarpOpContractToMma: per-lane fragments + nvgpu.mma.sync)
#   -> NVVM -> PTX (sm_90, LLVM NVPTX backend, no CUDA toolkit).
# Output: generated/warp_mma_sm90.ptx. Run with tools/run_warp_mma on a GPU node.
set -e
export PATH="/opt/homebrew/opt/llvm/bin:$PATH"
mkdir -p generated

build/tools/mir-opt/mir-opt test/warp_mma_kernel.mlir \
    --mir-warp-distribute --canonicalize --cse --lower-affine \
| mlir-opt --pass-pipeline="builtin.module(gpu.module(convert-nvgpu-to-nvvm,convert-vector-to-llvm))" \
| mlir-opt --gpu-lower-to-nvvm-pipeline="cubin-chip=sm_90 cubin-format=isa" \
| python3 test/extract_ptx.py generated/warp_mma_sm90.ptx
grep -c "mma.sync.aligned.m8n8k4" generated/warp_mma_sm90.ptx
