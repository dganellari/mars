#a = affine_map<(m, n, k) -> (m, k)>
#b = affine_map<(m, n, k) -> (k, n)>
#bt = affine_map<(m, n, k) -> (n, k)>
#c = affine_map<(m, n, k) -> (m, n)>
gpu.module @mir_kernels [#nvvm.target<chip = "sm_90", O = 3>] {
  gpu.func @bsweep(%Bt: memref<8x8xf64>, %U: memref<8x64xf64>, %Interp: memref<8x64xf64>) kernel {
    %z = arith.constant 0.0 : f64
    %lane = gpu.thread_id x
    %c0 = arith.constant 0 : index
    %c4 = arith.constant 4 : index
    %c8 = arith.constant 8 : index
    %c12 = arith.constant 12 : index
    %c16 = arith.constant 16 : index
    %c20 = arith.constant 20 : index
    %c24 = arith.constant 24 : index
    %c28 = arith.constant 28 : index
    %c32 = arith.constant 32 : index
    %c36 = arith.constant 36 : index
    %c40 = arith.constant 40 : index
    %c44 = arith.constant 44 : index
    %c48 = arith.constant 48 : index
    %c52 = arith.constant 52 : index
    %c56 = arith.constant 56 : index
    %c60 = arith.constant 60 : index
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1 = arith.constant dense<0.0> : vector<8x8xf64>
      %a2 = vector.transfer_read %Bt[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b3 = vector.transfer_read %U[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r4 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a2, %b3, %acc1 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a5 = vector.transfer_read %Bt[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b6 = vector.transfer_read %U[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r7 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a5, %b6, %r4 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r7, %Interp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc8 = arith.constant dense<0.0> : vector<8x8xf64>
      %a9 = vector.transfer_read %Bt[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b10 = vector.transfer_read %U[%c0, %c8], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r11 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a9, %b10, %acc8 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a12 = vector.transfer_read %Bt[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b13 = vector.transfer_read %U[%c4, %c8], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r14 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a12, %b13, %r11 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r14, %Interp[%c0, %c8] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc15 = arith.constant dense<0.0> : vector<8x8xf64>
      %a16 = vector.transfer_read %Bt[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b17 = vector.transfer_read %U[%c0, %c16], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r18 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a16, %b17, %acc15 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a19 = vector.transfer_read %Bt[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b20 = vector.transfer_read %U[%c4, %c16], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r21 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a19, %b20, %r18 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r21, %Interp[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc22 = arith.constant dense<0.0> : vector<8x8xf64>
      %a23 = vector.transfer_read %Bt[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b24 = vector.transfer_read %U[%c0, %c24], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r25 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a23, %b24, %acc22 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a26 = vector.transfer_read %Bt[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b27 = vector.transfer_read %U[%c4, %c24], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r28 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a26, %b27, %r25 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r28, %Interp[%c0, %c24] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc29 = arith.constant dense<0.0> : vector<8x8xf64>
      %a30 = vector.transfer_read %Bt[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b31 = vector.transfer_read %U[%c0, %c32], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r32 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a30, %b31, %acc29 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a33 = vector.transfer_read %Bt[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b34 = vector.transfer_read %U[%c4, %c32], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r35 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a33, %b34, %r32 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r35, %Interp[%c0, %c32] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc36 = arith.constant dense<0.0> : vector<8x8xf64>
      %a37 = vector.transfer_read %Bt[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b38 = vector.transfer_read %U[%c0, %c40], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r39 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a37, %b38, %acc36 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a40 = vector.transfer_read %Bt[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b41 = vector.transfer_read %U[%c4, %c40], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r42 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a40, %b41, %r39 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r42, %Interp[%c0, %c40] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc43 = arith.constant dense<0.0> : vector<8x8xf64>
      %a44 = vector.transfer_read %Bt[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b45 = vector.transfer_read %U[%c0, %c48], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r46 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a44, %b45, %acc43 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a47 = vector.transfer_read %Bt[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b48 = vector.transfer_read %U[%c4, %c48], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r49 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a47, %b48, %r46 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r49, %Interp[%c0, %c48] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc50 = arith.constant dense<0.0> : vector<8x8xf64>
      %a51 = vector.transfer_read %Bt[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b52 = vector.transfer_read %U[%c0, %c56], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r53 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a51, %b52, %acc50 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a54 = vector.transfer_read %Bt[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b55 = vector.transfer_read %U[%c4, %c56], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r56 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a54, %b55, %r53 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r56, %Interp[%c0, %c56] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64>
    }
    gpu.return
  }
}
