#!/bin/sh
# FUSED FULL-OPERATOR kernel: the whole Knaus Alg-2 apply, one thread per
# element, 128 threads/block, ONE gpu kernel (vs 111 naive launches).
# Usage: ./make_ptx_full.sh [p]   (default p=7; also builds p=3 gate variant)
set -e
export PATH="/opt/homebrew/opt/llvm/bin:$PATH"
mkdir -p generated
P=${1:-7}

python3 -c "
import sys; sys.path.insert(0, '../marsir-compiler')
from marsir import parse_spec_file, synthesize
from marsir.backends import mlir_ir
ea = synthesize(parse_spec_file('../marsir-compiler/specs/laplacian.op'))
sys.stdout.write(mlir_ir.emit_full_batched(ea, p=$P, tpb=128))
" \
| ./build/tools/mir-opt/mir-opt --convert-mir-to-linalg \
| mlir-opt --one-shot-bufferize \
| ./build/tools/mir-opt/mir-opt --mir-lower-copies \
| mlir-opt --promote-buffers-to-stack="max-alloc-size-in-bytes=65536" --convert-linalg-to-loops \
| mlir-opt --convert-parallel-loops-to-gpu --gpu-kernel-outlining \
| mlir-opt --gpu-lower-to-nvvm-pipeline="cubin-chip=sm_90 cubin-format=isa" \
| python3 test/extract_ptx.py generated/full_apply_p${P}_sm90.ptx 0
