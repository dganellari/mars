// GPU-executable gate for the WarpOpContractToMma bridge: the m8n8k4 contract
// written INSIDE a warp region (as --mir-warp-wrap would leave it), lane taken
// from thread_id.x. --mir-warp-distribute must dissolve the region into
// per-lane fragment reads + nvgpu.mma.sync + a per-lane C write.
// Launch contract: grid 1, block (32,1,1). C := A*B + C.
#a = affine_map<(m, n, k) -> (m, k)>
#b = affine_map<(m, n, k) -> (k, n)>
#c = affine_map<(m, n, k) -> (m, n)>
gpu.module @mir_kernels [#nvvm.target<chip = "sm_90", O = 3>] {
  gpu.func @warp_mma(%A: memref<8x4xf64>, %B: memref<4x8xf64>, %C: memref<8x8xf64>) kernel {
    %c0 = arith.constant 0 : index
    %z = arith.constant 0.0 : f64
    %lane = gpu.thread_id x
    vector.warp_execute_on_lane_0(%lane)[32] {
      %va = vector.transfer_read %A[%c0, %c0], %z {in_bounds = [true, true]} : memref<8x4xf64>, vector<8x4xf64>
      %vb = vector.transfer_read %B[%c0, %c0], %z {in_bounds = [true, true]} : memref<4x8xf64>, vector<4x8xf64>
      %vc = vector.transfer_read %C[%c0, %c0], %z {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %r = vector.contract {indexing_maps = [#a, #b, #c],
            iterator_types = ["parallel", "parallel", "reduction"], kind = #vector.kind<add>}
            %va, %vb, %vc : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r, %C[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
    }
    gpu.return
  }
}
