// GPU gate for the TRANSPOSED-B arm of WarpOpContractToMma: an axis-1 rank-2
// contract out[m,n] = sum_k in[m,k]*W[n,k]  (= in @ W^T), the shape the sum-
// factorization emits for D/W applied on tensor axis 1. W is stored [n,k], so
// its mma B fragment is read at [L/4, L%4] (vs [L%4, L/4] for a standard [k,n]
// B). Two chained k=4 slabs cover the full 8-wide reduction; the bridge fires
// per slab and chains C in registers.
// Launch: grid 1, block (32,1,1). C := in @ W^T + C0.
#a  = affine_map<(m, n, k) -> (m, k)>
#bt = affine_map<(m, n, k) -> (n, k)>
#c  = affine_map<(m, n, k) -> (m, n)>
gpu.module @mir_kernels [#nvvm.target<chip = "sm_90", O = 3>] {
  gpu.func @warp_mmt(%In: memref<8x8xf64>, %W: memref<8x8xf64>, %C: memref<8x8xf64>) kernel {
    %c0 = arith.constant 0 : index
    %c4 = arith.constant 4 : index
    %z = arith.constant 0.0 : f64
    %lane = gpu.thread_id x
    vector.warp_execute_on_lane_0(%lane)[32] {
      %a0 = vector.transfer_read %In[%c0, %c0], %z {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x4xf64>
      %a1 = vector.transfer_read %In[%c0, %c4], %z {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x4xf64>
      %b0 = vector.transfer_read %W[%c0, %c0], %z {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1 = vector.transfer_read %W[%c0, %c4], %z {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x4xf64>
      %vc = vector.transfer_read %C[%c0, %c0], %z {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %r0 = vector.contract {indexing_maps = [#a, #bt, #c],
            iterator_types = ["parallel", "parallel", "reduction"], kind = #vector.kind<add>}
            %a0, %b0, %vc : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %r1 = vector.contract {indexing_maps = [#a, #bt, #c],
            iterator_types = ["parallel", "parallel", "reduction"], kind = #vector.kind<add>}
            %a1, %b1, %r0 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1, %C[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
    }
    gpu.return
  }
}
