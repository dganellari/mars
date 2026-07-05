// Iteration C experiment: payload as func.func + scf.forall over elements, so
// hoist_redundant_vector_transfers (func.func-only) can lift the Y-tile
// accumulator out of the k-loop BEFORE mapping to GPU. Pipeline sketch:
//   tile matmul -> vectorize (func) -> HOIST -> map forall to gpu blocks
func.func @sweep(%D: memref<8x8xf64>, %U: memref<?x8x64xf64>, %Y: memref<?x8x64xf64>) {
  %c0 = arith.constant 0 : index
  %E = memref.dim %U, %c0 : memref<?x8x64xf64>
  scf.forall (%e) in (%E) {
    %Ue = memref.subview %U[%e, 0, 0] [1, 8, 64] [1, 1, 1]
          : memref<?x8x64xf64> to memref<8x64xf64, strided<[64, 1], offset: ?>>
    %Ye = memref.subview %Y[%e, 0, 0] [1, 8, 64] [1, 1, 1]
          : memref<?x8x64xf64> to memref<8x64xf64, strided<[64, 1], offset: ?>>
    linalg.matmul ins(%D, %Ue : memref<8x8xf64>, memref<8x64xf64, strided<[64, 1], offset: ?>>)
                  outs(%Ye : memref<8x64xf64, strided<[64, 1], offset: ?>>)
  }
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
    %fh = transform.structured.hoist_redundant_vector_transfers %fv
          : (!transform.any_op) -> !transform.any_op
    transform.yield
  }
}
