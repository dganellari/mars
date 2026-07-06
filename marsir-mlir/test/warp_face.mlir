#a = affine_map<(m, n, k) -> (m, k)>
#b = affine_map<(m, n, k) -> (k, n)>
#bt = affine_map<(m, n, k) -> (n, k)>
#c = affine_map<(m, n, k) -> (m, n)>
gpu.module @mir_kernels [#nvvm.target<chip = "sm_90", O = 3>] {
  gpu.func @face(%interp: memref<8x8xf64>, %deriv: memref<8x8xf64>, %Dm: memref<8x8xf64>, %W: memref<8x8xf64>, %g0: memref<8x8xf64>, %g1: memref<8x8xf64>, %g2: memref<8x8xf64>, %intf: memref<8x8xf64>) workgroup(%dt1g: memref<8x8xf64, #gpu.address_space<workgroup>>, %dt2g: memref<8x8xf64, #gpu.address_space<workgroup>>, %flux: memref<8x8xf64, #gpu.address_space<workgroup>>, %tmp: memref<8x8xf64, #gpu.address_space<workgroup>>) kernel {
    %c0 = arith.constant 0 : index
    %c4 = arith.constant 4 : index
    %z = arith.constant 0.0 : f64
    %lane = gpu.thread_id x
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1 = arith.constant dense<0.0> : vector<8x8xf64>
      %a2 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b3 = vector.transfer_read %interp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<4x8xf64>
      %r4 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a2, %b3, %acc1 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a5 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b6 = vector.transfer_read %interp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<4x8xf64>
      %r7 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a5, %b6, %r4 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g8 = vector.transfer_read %g1[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x8xf64>
      %sc9 = arith.mulf %r7, %g8 : vector<8x8xf64>
      vector.transfer_write %sc9, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc10 = arith.constant dense<0.0> : vector<8x8xf64>
      %a11 = vector.transfer_read %interp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b12 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r13 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a11, %b12, %acc10 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a14 = vector.transfer_read %interp[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b15 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r16 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a14, %b15, %r13 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g17 = vector.transfer_read %g0[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x8xf64>
      %sc18 = arith.mulf %r16, %g17 : vector<8x8xf64>
      vector.transfer_write %sc18, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d19 = vector.transfer_read %deriv[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x8xf64>
      %g220 = vector.transfer_read %g2[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x8xf64>
      %t21 = arith.mulf %d19, %g220 : vector<8x8xf64>
      %a22 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t223 = arith.addf %t21, %a22 : vector<8x8xf64>
      %b24 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl25 = arith.addf %t223, %b24 : vector<8x8xf64>
      vector.transfer_write %fl25, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc26 = arith.constant dense<0.0> : vector<8x8xf64>
      %a27 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b28 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r29 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a27, %b28, %acc26 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a30 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b31 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r32 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a30, %b31, %r29 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r32, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc33 = arith.constant dense<0.0> : vector<8x8xf64>
      %a34 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b35 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r36 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a34, %b35, %acc33 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a37 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b38 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r39 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a37, %b38, %r36 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r39, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64>
    }
    gpu.return
  }
}
