#!/bin/sh
# v5: register accumulator (zero-init) + 8 elements/block via explicit GPU
# mapping (outer parallel -> block_x, inner -> thread_y). Launch block=(32,8).
set -e
export PATH="/opt/homebrew/opt/llvm/bin:$PATH"
mkdir -p generated
mlir-opt test/ptx_kernel_batched_v5.mlir --transform-interpreter \
| python3 test/strip_transform.py \
| ./build/tools/mir-opt/mir-opt --mir-hoist-transfer-pairs="assume-zero-init=true" \
| mlir-opt --convert-vector-to-gpu="use-nvgpu=true" --canonicalize --cse \
| mlir-opt --convert-parallel-loops-to-gpu --gpu-kernel-outlining \
| sed 's/ attributes {known_block_size = array<i32: [0-9]*, [0-9]*, [0-9]*>}//' \
| mlir-opt --pass-pipeline="builtin.module(gpu.module(convert-nvgpu-to-nvvm,convert-vector-to-llvm))" \
| mlir-opt --gpu-lower-to-nvvm-pipeline="cubin-chip=sm_90 cubin-format=isa" \
| python3 test/extract_ptx.py generated/sumfac_batched_v5_sm90.ptx
