#!/bin/sh
# v4 pipeline: the accumulator-registerized batched sweep.
#   tile+vectorize (schedule) -> OUR --mir-hoist-transfer-pairs (Y tile becomes
#   a k-loop iter_arg; upstream hoisting declines this) -> convert to nvgpu at
#   func level (where the iter_args pattern converts) -> forall -> gpu blocks ->
#   outline -> PTX. Y store traffic per element: 16 vs v1's 32.
# The known_block_size=1 attr from the parallel-loop mapping is stripped: the
# kernel is warp-cooperative (launch with 32 threads; lane distribution via
# gpu.lane_id), so a baked-in blockDim of 1 would be a lie to the optimizer.
set -e
export PATH="/opt/homebrew/opt/llvm/bin:$PATH"
mkdir -p generated

mlir-opt test/forall_hoist_experiment.mlir --transform-interpreter \
| python3 test/strip_transform.py \
| ./build/tools/mir-opt/mir-opt --mir-hoist-transfer-pairs \
| mlir-opt --convert-vector-to-gpu="use-nvgpu=true" --canonicalize --cse \
| mlir-opt --scf-forall-to-parallel --gpu-map-parallel-loops --convert-parallel-loops-to-gpu --gpu-kernel-outlining \
| sed 's/ attributes {known_block_size = array<i32: 1, 1, 1>}//' \
| mlir-opt --pass-pipeline="builtin.module(gpu.module(convert-nvgpu-to-nvvm,convert-vector-to-llvm))" \
| mlir-opt --gpu-lower-to-nvvm-pipeline="cubin-chip=sm_90 cubin-format=isa" \
| python3 test/extract_ptx.py generated/sumfac_batched_v4_sm90.ptx
