#!/bin/sh
# The mir -> FP64 tensor-core (DMMA) pipeline, end to end:
#   mir.contract -> linalg.matmul (unfolded) -> [dmma_schedule: tile m8n8k4 +
#   vectorize] -> vector.contract -> nvgpu.mma.sync -> nvvm.mma.sync
# Run from marsir-mlir/ after `ninja -C build mir-opt`. Prints the final IR;
# pass -q to just count the mma ops.
set -e
export PATH="/opt/homebrew/opt/llvm/bin:$PATH"
OPT=./build/tools/mir-opt/mir-opt

out=$($OPT test/contract.mlir --convert-mir-to-linalg \
      --one-shot-bufferize="bufferize-function-boundaries=true function-boundary-type-conversion=identity-layout-map" \
  | mlir-opt --transform-preload-library="transform-library-paths=test/dmma_schedule.mlir" \
             --transform-interpreter \
  | mlir-opt --convert-vector-to-gpu="use-nvgpu=true" --canonicalize --convert-nvgpu-to-nvvm)

if [ "$1" = "-q" ]; then
  echo "$out" | grep -oE "nvvm\.mma\.sync|vector\.contract" | sort | uniq -c
else
  echo "$out"
fi
