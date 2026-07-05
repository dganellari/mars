// v5 = v4 (register accumulator via --mir-hoist-transfer-pairs) + 8 elements
// per block (explicit GPU mapping: outer parallel -> block_x, inner -> thread_y,
// so warp w = thread_id.y handles element e = b*8 + w) + zero-init accumulator
// (pass option assume-zero-init: the operator pre-zeroes Y, so semantics are
// Y = D*U and the 4KiB/elem zero-read disappears).
func.func @sweepv5(%D: memref<8x8xf64>, %U: memref<?x8x64xf64>, %Y: memref<?x8x64xf64>) {
  %c0 = arith.constant 0 : index
  %c1 = arith.constant 1 : index
  %c8 = arith.constant 8 : index
  %E = memref.dim %U, %c0 : memref<?x8x64xf64>
  %B = arith.divui %E, %c8 : index
  scf.parallel (%b) = (%c0) to (%B) step (%c1) {
    %be = arith.muli %b, %c8 : index
    scf.parallel (%w) = (%c0) to (%c8) step (%c1) {
      %e = arith.addi %be, %w : index
      %Ue = memref.subview %U[%e, 0, 0] [1, 8, 64] [1, 1, 1]
            : memref<?x8x64xf64> to memref<8x64xf64, strided<[64, 1], offset: ?>>
      %Ye = memref.subview %Y[%e, 0, 0] [1, 8, 64] [1, 1, 1]
            : memref<?x8x64xf64> to memref<8x64xf64, strided<[64, 1], offset: ?>>
      linalg.matmul ins(%D, %Ue : memref<8x8xf64>, memref<8x64xf64, strided<[64, 1], offset: ?>>)
                    outs(%Ye : memref<8x64xf64, strided<[64, 1], offset: ?>>)
    } {mapping = [#gpu.loop_dim_map<processor = thread_y, map = (d0) -> (d0), bound = (d0) -> (d0)>]}
  } {mapping = [#gpu.loop_dim_map<processor = block_x, map = (d0) -> (d0), bound = (d0) -> (d0)>]}
  return
}

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
