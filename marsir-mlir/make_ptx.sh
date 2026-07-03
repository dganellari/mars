#!/bin/sh
# Serialize REAL PTX (sm_90, H100/GH200) for the tensor-core sum-factorization
# sweep, via LLVM's NVPTX backend (no CUDA toolkit needed):
#   test/ptx_kernel.mlir -> schedule -> nvgpu.mma.sync -> NVVM -> PTX
# Output: generated/sumfac_sm90.ptx (contains mma.sync.aligned.m8n8k4...f64).
set -e
export PATH="/opt/homebrew/opt/llvm/bin:$PATH"
mkdir -p generated

serialize() {
  mlir-opt "$1" --transform-interpreter \
  | mlir-opt --pass-pipeline="builtin.module(gpu.module(gpu.func(convert-vector-to-gpu{use-nvgpu=true})),canonicalize,cse)" \
  | mlir-opt --pass-pipeline="builtin.module(gpu.module(convert-nvgpu-to-nvvm,convert-vector-to-llvm))" \
  | mlir-opt --gpu-lower-to-nvvm-pipeline="cubin-chip=sm_90 cubin-format=isa" \
  | python3 test/extract_ptx.py "$2"
}
serialize test/ptx_kernel.mlir         generated/sumfac_sm90.ptx
serialize test/ptx_kernel_batched.mlir    generated/sumfac_batched_sm90.ptx
serialize test/ptx_kernel_batched_v3.mlir generated/sumfac_batched_v3_sm90.ptx
