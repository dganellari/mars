// PTX-emission payload: the sum-factorization sweep (unfolded matmul) as a GPU
// kernel, with an sm_90 target (H100/GH200). The embedded transform schedule
// tiles to the m8n8k4 DMMA shape and vectorizes; the driver script then converts
// vector -> nvgpu.mma.sync -> nvvm and serializes REAL PTX via LLVM's NVPTX
// backend (format=isa; no CUDA toolkit needed).
gpu.module @mir_kernels [#nvvm.target<chip = "sm_90", O = 3>] {
  gpu.func @sumfac_sweep(%D: memref<8x8xf64>, %U: memref<8x64xf64>, %Y: memref<8x64xf64>) kernel {
    linalg.matmul ins(%D, %U : memref<8x8xf64>, memref<8x64xf64>)
                  outs(%Y : memref<8x64xf64>)
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
