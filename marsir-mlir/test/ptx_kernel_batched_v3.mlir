// TUNED batched sweep, iteration A+B over ptx_kernel_batched.mlir:
//   A. hoist_redundant_vector_transfers: the Y accumulator tile is carried in
//      REGISTERS across the k-loop (read once before, write once after),
//      instead of a global read+write per k-tile.
//   B. 8 elements per block (256 threads, 8 warps): warp w = thread_id.y
//      handles element  e = block_id.x * 8 + w  -> occupancy like the hand
//      kernels' ElemsPerBlock, grid = E/8.
gpu.module @mir_kernels [#nvvm.target<chip = "sm_90", O = 3>] {
  gpu.func @sumfac_sweep_batched_v3(%D: memref<8x8xf64>,
                                    %U: memref<?x8x64xf64>,
                                    %Y: memref<?x8x64xf64>) kernel {
    %b = gpu.block_id x
    %w = gpu.thread_id y
    %c8 = arith.constant 8 : index
    %be = arith.muli %b, %c8 : index
    %e = arith.addi %be, %w : index
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
    // Unroll k BEFORE vectorizing (vectorize invalidates loop handles), so the
    // two mma tiles are straight-line and the accumulator can forward in regs.
    transform.loop.unroll %loops#1 { factor = 2 } : !transform.any_op
    %f = transform.structured.match ops{["gpu.func"]} in %root
          : (!transform.any_op) -> !transform.any_op
    %fv = transform.structured.vectorize_children_and_apply_patterns %f
          : (!transform.any_op) -> !transform.any_op
    transform.yield
  }
}
