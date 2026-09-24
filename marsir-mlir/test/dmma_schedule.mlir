// The FP64 tensor-core (DMMA) SCHEDULE for the sum-factorization sweeps.
//
// A transform-dialect library: "how to compile it", separate from "what to
// compute" (the mir ops). It shapes every contraction into m8n8k4-compatible
// 2-D tiles and vectorizes them to vector.contract, which --mir-chain-contracts
// (or --convert-vector-to-gpu{use-nvgpu}) then maps onto the tensor cores.
//
//   1. The axis-0 sweep arrives as an unfolded linalg.matmul ([P x n] @ [n x n^2]):
//      tile n^2 by 8 and k by 4.
//   2. Sweeps along axis 1 or 2 arrive as a rank-3 linalg.generic tagged
//      mir.batch_contract: a batch of 2-D contractions, one per index of the
//      leading axis. Tile that axis by 1 and k by 4 -- the two sweep dimensions
//      (<= 8 each) stay whole -- then strip the unit batch dimension so each tile
//      is a plain 2-D contraction.
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
    %bc = transform.structured.match attributes {mir.batch_contract} in %root
          : (!transform.any_op) -> !transform.any_op
    %btiled, %bloops:2 = transform.structured.tile_using_for %bc tile_sizes [1, 0, 0, 4]
          : (!transform.any_op) -> (!transform.any_op, !transform.any_op, !transform.any_op)
    %f = transform.structured.match ops{["func.func"]} in %root
          : (!transform.any_op) -> !transform.any_op
    %fv = transform.structured.vectorize_children_and_apply_patterns %f
          : (!transform.any_op) -> !transform.any_op
    transform.apply_patterns to %fv {
      transform.apply_patterns.vector.cast_away_vector_leading_one_dim
    } : !transform.any_op
    transform.yield
  }
}
