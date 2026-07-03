// The FP64 tensor-core (DMMA) SCHEDULE for the sum-factorization sweep.
//
// This is a transform-dialect library: a reusable "how to compile it" script,
// separate from the "what to compute" (the mir ops). It matches the unfolded
// matmul that --convert-mir-to-linalg emits for a mir.contract and:
//   1. tiles n by 8 and k by 4  ->  m8n8k4 tiles (the FP64 DMMA shape), and
//   2. vectorizes each tile     ->  vector.contract,
// which --convert-vector-to-gpu{use-nvgpu} then maps onto nvgpu.mma.sync and
// --convert-nvgpu-to-nvvm onto nvvm.mma.sync (the DMMA instruction).
//
// Usage:
//   mlir-opt <bufferized payload> \
//     --transform-preload-library="transform-library-paths=test/dmma_schedule.mlir" \
//     --transform-interpreter
module attributes {transform.with_named_sequence} {
  transform.named_sequence @__transform_main(%root: !transform.any_op {transform.readonly}) {
    %mm = transform.structured.match ops{["linalg.matmul"]} in %root
          : (!transform.any_op) -> !transform.any_op
    %tiled, %loops:2 = transform.structured.tile_using_for %mm tile_sizes [0, 8, 4]
          : (!transform.any_op) -> (!transform.any_op, !transform.any_op, !transform.any_op)
    %f = transform.structured.match ops{["func.func"]} in %root
          : (!transform.any_op) -> !transform.any_op
    %fv = transform.structured.vectorize_children_and_apply_patterns %f
          : (!transform.any_op) -> !transform.any_op
    transform.yield
  }
}
