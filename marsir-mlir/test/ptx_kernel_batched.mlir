// BATCHED sum-factorization sweep: one warp per element (the hand-kernel shape).
//   Y[e] += D(8x8) * U[e](8x64)   for e = blockIdx.x in [0, E)
// The batch dim is dynamic (?), so one PTX serves any element count. The same
// m8n8k4 DMMA schedule as the single-element kernel applies to the per-element
// matmul on the subview.
gpu.module @mir_kernels [#nvvm.target<chip = "sm_90", O = 3>] {
  gpu.func @sumfac_sweep_batched(%D: memref<8x8xf64>,
                                 %U: memref<?x8x64xf64>,
                                 %Y: memref<?x8x64xf64>) kernel {
    %e = gpu.block_id x
    %Ue = memref.subview %U[%e, 0, 0] [1, 8, 64] [1, 1, 1]
          : memref<?x8x64xf64> to memref<8x64xf64, strided<[64, 1], offset: ?>>
    %Ye = memref.subview %Y[%e, 0, 0] [1, 8, 64] [1, 1, 1]
          : memref<?x8x64xf64> to memref<8x64xf64, strided<[64, 1], offset: ?>>
    linalg.matmul ins(%D, %Ue : memref<8x8xf64>, memref<8x64xf64, strided<[64, 1], offset: ?>>)
                  outs(%Ye : memref<8x64xf64, strided<[64, 1], offset: ?>>)
    gpu.return
  }
}

module attributes {transform.with_named_sequence} {
  transform.named_sequence @__transform_main(%root: !transform.any_op {transform.readonly}) {
    %mm = transform.structured.match ops{["linalg.matmul"]} in %root
          : (!transform.any_op) -> !transform.any_op
    %tiled, %loops:2 = transform.structured.tile_using_for %mm tile_sizes [0, 8, 4]
          : (!transform.any_op) -> (!transform.any_op, !transform.any_op, !transform.any_op)
    %f = transform.structured.match ops{["gpu.func"]} in %root
          : (!transform.any_op) -> !transform.any_op
    %fv = transform.structured.vectorize_children_and_apply_patterns %f
          : (!transform.any_op) -> !transform.any_op
    transform.yield
  }
}
