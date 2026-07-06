#a = affine_map<(m, n, k) -> (m, k)>
#b = affine_map<(m, n, k) -> (k, n)>
#bt = affine_map<(m, n, k) -> (n, k)>
#c = affine_map<(m, n, k) -> (m, n)>
gpu.module @mir_kernels [#nvvm.target<chip = "sm_90", O = 3>] {
  gpu.func @two(%M: memref<8x8xf64>, %X: memref<8x8xf64>, %Out: memref<8x8xf64>)
      workgroup(%S1: memref<8x8xf64, #gpu.address_space<workgroup>>) kernel {
    %c0 = arith.constant 0 : index
    %c4 = arith.constant 4 : index
    %z = arith.constant 0.0 : f64
    %lane = gpu.thread_id x
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1 = arith.constant dense<0.0> : vector<8x8xf64>
      %a2 = vector.transfer_read %M[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b3 = vector.transfer_read %X[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<4x8xf64>
      %r4 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a2, %b3, %acc1 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a5 = vector.transfer_read %M[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b6 = vector.transfer_read %X[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<4x8xf64>
      %r7 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a5, %b6, %r4 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r7, %S1[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc8 = arith.constant dense<0.0> : vector<8x8xf64>
      %a9 = vector.transfer_read %S1[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b10 = vector.transfer_read %M[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r11 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a9, %b10, %acc8 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a12 = vector.transfer_read %S1[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b13 = vector.transfer_read %M[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r14 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a12, %b13, %r11 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r14, %Out[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64>
    }
    gpu.return
  }
}
