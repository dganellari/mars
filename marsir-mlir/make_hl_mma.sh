#!/bin/sh
# HIGH-LEVEL mir -> FP64 TENSOR CORES, by passes only. No Python low-level emit.
#
#   laplacian.op -> mlir_ir.emit_full         the operator as high-level mir
#                                             (mir.contract + mir.flux only)
#     --convert-mir-to-linalg                 linalg
#     one-shot-bufferize                      memrefs
#     dmma_schedule.mlir (transform dialect)  tile m8n8k4 + vectorize
#     --mir-chain-contracts                   per-lane nvgpu.mma.sync
#     --canonicalize --cse                    drop the now-dead full-width reads
#     --loop-invariant-code-motion            so the fragment indices are invariant
#     --mir-hoist-transfer-pairs              C-fragment into a loop-carried register
#
# ORDER MATTERS: without the DCE the accumulator view still has the dead
# full-width read as a third user, and without LICM its column index is computed
# inside the loop -- either one makes the hoist decline.
#
# This is the route the tutorial's Chapter 0 calls the destination: the Python is
# a thin front-end emitting the dialect, everything below it is a pass.
#
# Expect 24 mma and 6 declined contracts. The 6 all carry P = 7 (the number of
# subcontrol-surface faces): 2 are the B-sweep at 7xN and 4 are rank-3 batched
# with a 7-extent. m8n8k4 fixes m = 8 in HARDWARE, so those need the front end to
# pad P to 8 -- which is exactly what mlir_warp.py does silently by using an 8x8
# Btil. Until then the pass DECLINES them rather than mangling them.
set -e
export PATH="/opt/homebrew/opt/llvm/bin:$PATH"
OPT=./build/tools/mir-opt/mir-opt
P=${1:-7}

python3 -c "
import sys; sys.path.insert(0, '../marsir-compiler')
from marsir import parse_spec_file, synthesize
from marsir.backends import mlir_ir
ea = synthesize(parse_spec_file('../marsir-compiler/specs/laplacian.op'))
sys.stdout.write(mlir_ir.emit_full(ea, p=$P))
" \
| $OPT --convert-mir-to-linalg \
       --one-shot-bufferize="bufferize-function-boundaries=true function-boundary-type-conversion=identity-layout-map" \
| mlir-opt --transform-preload-library="transform-library-paths=test/dmma_schedule.mlir" \
           --transform-interpreter \
| $OPT --mir-chain-contracts --canonicalize --cse \
| mlir-opt --loop-invariant-code-motion \
| $OPT --mir-hoist-transfer-pairs > generated/hl_mma_p$P.mlir

mma=$(grep -c "nvgpu.mma.sync" generated/hl_mma_p$P.mlir || true)
left=$(grep -c "vector.contract" generated/hl_mma_p$P.mlir || true)
$OPT generated/hl_mma_p$P.mlir -o /dev/null
acc=$(grep -c "iter_args(.*vector<1x2xf64>" generated/hl_mma_p$P.mlir || true)
echo "high-level mir -> mma (p=$P): $mma nvgpu.mma.sync, $left declined, $acc register-carried accumulators, verifies"
