#a = affine_map<(m, n, k) -> (m, k)>
#b = affine_map<(m, n, k) -> (k, n)>
#bt = affine_map<(m, n, k) -> (n, k)>
#c = affine_map<(m, n, k) -> (m, n)>
gpu.module @mir_kernels [#nvvm.target<chip = "sm_90", O = 3>] {
  gpu.func @dir0(%u2: memref<8x64xf64>, %Btil: memref<8x8xf64>, %Dtil: memref<8x8xf64>, %Dm: memref<8x8xf64>, %W: memref<8x8xf64>, %g0all: memref<8x8x8xf64>, %g1all: memref<8x8x8xf64>, %g2all: memref<8x8x8xf64>, %y3: memref<8x8x8xf64>) workgroup(%interp_all: memref<8x64xf64, #gpu.address_space<workgroup>>, %deriv_all: memref<8x64xf64, #gpu.address_space<workgroup>>, %flux: memref<8x8xf64, #gpu.address_space<workgroup>>, %tmp: memref<8x8xf64, #gpu.address_space<workgroup>>, %intf: memref<8x8xf64, #gpu.address_space<workgroup>>) kernel {
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
      %a2 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b3 = vector.transfer_read %u2[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r4 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a2, %b3, %acc1 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a5 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b6 = vector.transfer_read %u2[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r7 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a5, %b6, %r4 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r7, %interp_all[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc8 = arith.constant dense<0.0> : vector<8x8xf64>
      %a9 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b10 = vector.transfer_read %u2[%c0, %c8], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r11 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a9, %b10, %acc8 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a12 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b13 = vector.transfer_read %u2[%c4, %c8], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r14 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a12, %b13, %r11 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r14, %interp_all[%c0, %c8] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc15 = arith.constant dense<0.0> : vector<8x8xf64>
      %a16 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b17 = vector.transfer_read %u2[%c0, %c16], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r18 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a16, %b17, %acc15 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a19 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b20 = vector.transfer_read %u2[%c4, %c16], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r21 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a19, %b20, %r18 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r21, %interp_all[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc22 = arith.constant dense<0.0> : vector<8x8xf64>
      %a23 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b24 = vector.transfer_read %u2[%c0, %c24], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r25 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a23, %b24, %acc22 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a26 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b27 = vector.transfer_read %u2[%c4, %c24], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r28 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a26, %b27, %r25 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r28, %interp_all[%c0, %c24] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc29 = arith.constant dense<0.0> : vector<8x8xf64>
      %a30 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b31 = vector.transfer_read %u2[%c0, %c32], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r32 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a30, %b31, %acc29 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a33 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b34 = vector.transfer_read %u2[%c4, %c32], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r35 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a33, %b34, %r32 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r35, %interp_all[%c0, %c32] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc36 = arith.constant dense<0.0> : vector<8x8xf64>
      %a37 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b38 = vector.transfer_read %u2[%c0, %c40], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r39 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a37, %b38, %acc36 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a40 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b41 = vector.transfer_read %u2[%c4, %c40], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r42 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a40, %b41, %r39 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r42, %interp_all[%c0, %c40] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc43 = arith.constant dense<0.0> : vector<8x8xf64>
      %a44 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b45 = vector.transfer_read %u2[%c0, %c48], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r46 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a44, %b45, %acc43 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a47 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b48 = vector.transfer_read %u2[%c4, %c48], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r49 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a47, %b48, %r46 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r49, %interp_all[%c0, %c48] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc50 = arith.constant dense<0.0> : vector<8x8xf64>
      %a51 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b52 = vector.transfer_read %u2[%c0, %c56], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r53 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a51, %b52, %acc50 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a54 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b55 = vector.transfer_read %u2[%c4, %c56], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r56 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a54, %b55, %r53 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r56, %interp_all[%c0, %c56] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc57 = arith.constant dense<0.0> : vector<8x8xf64>
      %a58 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b59 = vector.transfer_read %u2[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r60 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a58, %b59, %acc57 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a61 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b62 = vector.transfer_read %u2[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r63 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a61, %b62, %r60 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r63, %deriv_all[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc64 = arith.constant dense<0.0> : vector<8x8xf64>
      %a65 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b66 = vector.transfer_read %u2[%c0, %c8], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r67 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a65, %b66, %acc64 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a68 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b69 = vector.transfer_read %u2[%c4, %c8], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r70 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a68, %b69, %r67 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r70, %deriv_all[%c0, %c8] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc71 = arith.constant dense<0.0> : vector<8x8xf64>
      %a72 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b73 = vector.transfer_read %u2[%c0, %c16], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r74 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a72, %b73, %acc71 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a75 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b76 = vector.transfer_read %u2[%c4, %c16], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r77 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a75, %b76, %r74 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r77, %deriv_all[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc78 = arith.constant dense<0.0> : vector<8x8xf64>
      %a79 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b80 = vector.transfer_read %u2[%c0, %c24], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r81 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a79, %b80, %acc78 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a82 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b83 = vector.transfer_read %u2[%c4, %c24], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r84 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a82, %b83, %r81 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r84, %deriv_all[%c0, %c24] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc85 = arith.constant dense<0.0> : vector<8x8xf64>
      %a86 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b87 = vector.transfer_read %u2[%c0, %c32], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r88 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a86, %b87, %acc85 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a89 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b90 = vector.transfer_read %u2[%c4, %c32], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r91 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a89, %b90, %r88 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r91, %deriv_all[%c0, %c32] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc92 = arith.constant dense<0.0> : vector<8x8xf64>
      %a93 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b94 = vector.transfer_read %u2[%c0, %c40], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r95 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a93, %b94, %acc92 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a96 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b97 = vector.transfer_read %u2[%c4, %c40], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r98 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a96, %b97, %r95 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r98, %deriv_all[%c0, %c40] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc99 = arith.constant dense<0.0> : vector<8x8xf64>
      %a100 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b101 = vector.transfer_read %u2[%c0, %c48], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r102 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a100, %b101, %acc99 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a103 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b104 = vector.transfer_read %u2[%c4, %c48], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r105 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a103, %b104, %r102 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r105, %deriv_all[%c0, %c48] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc106 = arith.constant dense<0.0> : vector<8x8xf64>
      %a107 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b108 = vector.transfer_read %u2[%c0, %c56], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r109 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a107, %b108, %acc106 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a110 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b111 = vector.transfer_read %u2[%c4, %c56], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r112 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a110, %b111, %r109 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r112, %deriv_all[%c0, %c56] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %v3113 = memref.reinterpret_cast %interp_all to offset: [0], sizes: [8, 8, 8], strides: [64, 8, 1] : memref<8x64xf64, #gpu.address_space<workgroup>> to memref<8x8x8xf64, #gpu.address_space<workgroup>>
    %v3114 = memref.reinterpret_cast %deriv_all to offset: [0], sizes: [8, 8, 8], strides: [64, 8, 1] : memref<8x64xf64, #gpu.address_space<workgroup>> to memref<8x8x8xf64, #gpu.address_space<workgroup>>
    %sv115 = memref.subview %v3113[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>
    %sv116 = memref.subview %v3114[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>
    %sv117 = memref.subview %g0all[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 0>>
    %sv118 = memref.subview %g1all[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 0>>
    %sv119 = memref.subview %g2all[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 0>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc120 = arith.constant dense<0.0> : vector<8x8xf64>
      %a121 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b122 = vector.transfer_read %sv115[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r123 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a121, %b122, %acc120 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a124 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b125 = vector.transfer_read %sv115[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r126 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a124, %b125, %r123 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g127 = vector.transfer_read %sv118[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>>, vector<8x8xf64>
      %dt128 = arith.mulf %r126, %g127 : vector<8x8xf64>
      %acc129 = arith.constant dense<0.0> : vector<8x8xf64>
      %a130 = vector.transfer_read %sv115[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b131 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r132 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a130, %b131, %acc129 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a133 = vector.transfer_read %sv115[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b134 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r135 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a133, %b134, %r132 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g136 = vector.transfer_read %sv117[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>>, vector<8x8xf64>
      %dt137 = arith.mulf %r135, %g136 : vector<8x8xf64>
      %dv138 = vector.transfer_read %sv116[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g139 = vector.transfer_read %sv119[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>>, vector<8x8xf64>
      %t140 = arith.mulf %dv138, %g139 : vector<8x8xf64>
      %t2141 = arith.addf %t140, %dt137 : vector<8x8xf64>
      %fl142 = arith.addf %t2141, %dt128 : vector<8x8xf64>
      vector.transfer_write %fl142, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc143 = arith.constant dense<0.0> : vector<8x8xf64>
      %a144 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b145 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r146 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a144, %b145, %acc143 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a147 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b148 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r149 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a147, %b148, %r146 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r149, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc150 = arith.constant dense<0.0> : vector<8x8xf64>
      %a151 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b152 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r153 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a151, %b152, %acc150 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a154 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b155 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r156 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a154, %b155, %r153 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r156, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv157 = memref.subview %y3[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 0>>
    %sv158 = memref.subview %y3[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 64>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa159 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r160 = vector.transfer_read %sv157[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>>, vector<8x8xf64>
      %m161 = arith.subf %r160, %fa159 : vector<8x8xf64>
      vector.transfer_write %m161, %sv157[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 0>>
      %r162 = vector.transfer_read %sv158[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<8x8xf64>
      %p163 = arith.addf %r162, %fa159 : vector<8x8xf64>
      vector.transfer_write %p163, %sv158[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 64>>
    }
    gpu.barrier
    %sv164 = memref.subview %v3113[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv165 = memref.subview %v3114[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv166 = memref.subview %g0all[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 64>>
    %sv167 = memref.subview %g1all[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 64>>
    %sv168 = memref.subview %g2all[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 64>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc169 = arith.constant dense<0.0> : vector<8x8xf64>
      %a170 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b171 = vector.transfer_read %sv164[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r172 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a170, %b171, %acc169 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a173 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b174 = vector.transfer_read %sv164[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r175 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a173, %b174, %r172 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g176 = vector.transfer_read %sv167[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<8x8xf64>
      %dt177 = arith.mulf %r175, %g176 : vector<8x8xf64>
      %acc178 = arith.constant dense<0.0> : vector<8x8xf64>
      %a179 = vector.transfer_read %sv164[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b180 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r181 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a179, %b180, %acc178 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a182 = vector.transfer_read %sv164[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b183 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r184 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a182, %b183, %r181 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g185 = vector.transfer_read %sv166[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<8x8xf64>
      %dt186 = arith.mulf %r184, %g185 : vector<8x8xf64>
      %dv187 = vector.transfer_read %sv165[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g188 = vector.transfer_read %sv168[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<8x8xf64>
      %t189 = arith.mulf %dv187, %g188 : vector<8x8xf64>
      %t2190 = arith.addf %t189, %dt186 : vector<8x8xf64>
      %fl191 = arith.addf %t2190, %dt177 : vector<8x8xf64>
      vector.transfer_write %fl191, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc192 = arith.constant dense<0.0> : vector<8x8xf64>
      %a193 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b194 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r195 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a193, %b194, %acc192 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a196 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b197 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r198 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a196, %b197, %r195 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r198, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc199 = arith.constant dense<0.0> : vector<8x8xf64>
      %a200 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b201 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r202 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a200, %b201, %acc199 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a203 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b204 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r205 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a203, %b204, %r202 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r205, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv206 = memref.subview %y3[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 64>>
    %sv207 = memref.subview %y3[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 128>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa208 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r209 = vector.transfer_read %sv206[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<8x8xf64>
      %m210 = arith.subf %r209, %fa208 : vector<8x8xf64>
      vector.transfer_write %m210, %sv206[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 64>>
      %r211 = vector.transfer_read %sv207[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<8x8xf64>
      %p212 = arith.addf %r211, %fa208 : vector<8x8xf64>
      vector.transfer_write %p212, %sv207[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 128>>
    }
    gpu.barrier
    %sv213 = memref.subview %v3113[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv214 = memref.subview %v3114[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv215 = memref.subview %g0all[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 128>>
    %sv216 = memref.subview %g1all[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 128>>
    %sv217 = memref.subview %g2all[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 128>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc218 = arith.constant dense<0.0> : vector<8x8xf64>
      %a219 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b220 = vector.transfer_read %sv213[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r221 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a219, %b220, %acc218 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a222 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b223 = vector.transfer_read %sv213[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r224 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a222, %b223, %r221 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g225 = vector.transfer_read %sv216[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<8x8xf64>
      %dt226 = arith.mulf %r224, %g225 : vector<8x8xf64>
      %acc227 = arith.constant dense<0.0> : vector<8x8xf64>
      %a228 = vector.transfer_read %sv213[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b229 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r230 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a228, %b229, %acc227 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a231 = vector.transfer_read %sv213[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b232 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r233 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a231, %b232, %r230 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g234 = vector.transfer_read %sv215[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<8x8xf64>
      %dt235 = arith.mulf %r233, %g234 : vector<8x8xf64>
      %dv236 = vector.transfer_read %sv214[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g237 = vector.transfer_read %sv217[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<8x8xf64>
      %t238 = arith.mulf %dv236, %g237 : vector<8x8xf64>
      %t2239 = arith.addf %t238, %dt235 : vector<8x8xf64>
      %fl240 = arith.addf %t2239, %dt226 : vector<8x8xf64>
      vector.transfer_write %fl240, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc241 = arith.constant dense<0.0> : vector<8x8xf64>
      %a242 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b243 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r244 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a242, %b243, %acc241 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a245 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b246 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r247 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a245, %b246, %r244 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r247, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc248 = arith.constant dense<0.0> : vector<8x8xf64>
      %a249 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b250 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r251 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a249, %b250, %acc248 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a252 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b253 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r254 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a252, %b253, %r251 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r254, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv255 = memref.subview %y3[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 128>>
    %sv256 = memref.subview %y3[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 192>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa257 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r258 = vector.transfer_read %sv255[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<8x8xf64>
      %m259 = arith.subf %r258, %fa257 : vector<8x8xf64>
      vector.transfer_write %m259, %sv255[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 128>>
      %r260 = vector.transfer_read %sv256[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<8x8xf64>
      %p261 = arith.addf %r260, %fa257 : vector<8x8xf64>
      vector.transfer_write %p261, %sv256[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 192>>
    }
    gpu.barrier
    %sv262 = memref.subview %v3113[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv263 = memref.subview %v3114[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv264 = memref.subview %g0all[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 192>>
    %sv265 = memref.subview %g1all[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 192>>
    %sv266 = memref.subview %g2all[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 192>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc267 = arith.constant dense<0.0> : vector<8x8xf64>
      %a268 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b269 = vector.transfer_read %sv262[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r270 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a268, %b269, %acc267 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a271 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b272 = vector.transfer_read %sv262[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r273 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a271, %b272, %r270 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g274 = vector.transfer_read %sv265[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<8x8xf64>
      %dt275 = arith.mulf %r273, %g274 : vector<8x8xf64>
      %acc276 = arith.constant dense<0.0> : vector<8x8xf64>
      %a277 = vector.transfer_read %sv262[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b278 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r279 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a277, %b278, %acc276 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a280 = vector.transfer_read %sv262[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b281 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r282 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a280, %b281, %r279 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g283 = vector.transfer_read %sv264[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<8x8xf64>
      %dt284 = arith.mulf %r282, %g283 : vector<8x8xf64>
      %dv285 = vector.transfer_read %sv263[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g286 = vector.transfer_read %sv266[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<8x8xf64>
      %t287 = arith.mulf %dv285, %g286 : vector<8x8xf64>
      %t2288 = arith.addf %t287, %dt284 : vector<8x8xf64>
      %fl289 = arith.addf %t2288, %dt275 : vector<8x8xf64>
      vector.transfer_write %fl289, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc290 = arith.constant dense<0.0> : vector<8x8xf64>
      %a291 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b292 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r293 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a291, %b292, %acc290 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a294 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b295 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r296 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a294, %b295, %r293 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r296, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc297 = arith.constant dense<0.0> : vector<8x8xf64>
      %a298 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b299 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r300 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a298, %b299, %acc297 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a301 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b302 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r303 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a301, %b302, %r300 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r303, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv304 = memref.subview %y3[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 192>>
    %sv305 = memref.subview %y3[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 256>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa306 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r307 = vector.transfer_read %sv304[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<8x8xf64>
      %m308 = arith.subf %r307, %fa306 : vector<8x8xf64>
      vector.transfer_write %m308, %sv304[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 192>>
      %r309 = vector.transfer_read %sv305[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<8x8xf64>
      %p310 = arith.addf %r309, %fa306 : vector<8x8xf64>
      vector.transfer_write %p310, %sv305[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 256>>
    }
    gpu.barrier
    %sv311 = memref.subview %v3113[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv312 = memref.subview %v3114[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv313 = memref.subview %g0all[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 256>>
    %sv314 = memref.subview %g1all[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 256>>
    %sv315 = memref.subview %g2all[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 256>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc316 = arith.constant dense<0.0> : vector<8x8xf64>
      %a317 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b318 = vector.transfer_read %sv311[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r319 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a317, %b318, %acc316 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a320 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b321 = vector.transfer_read %sv311[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r322 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a320, %b321, %r319 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g323 = vector.transfer_read %sv314[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<8x8xf64>
      %dt324 = arith.mulf %r322, %g323 : vector<8x8xf64>
      %acc325 = arith.constant dense<0.0> : vector<8x8xf64>
      %a326 = vector.transfer_read %sv311[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b327 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r328 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a326, %b327, %acc325 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a329 = vector.transfer_read %sv311[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b330 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r331 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a329, %b330, %r328 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g332 = vector.transfer_read %sv313[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<8x8xf64>
      %dt333 = arith.mulf %r331, %g332 : vector<8x8xf64>
      %dv334 = vector.transfer_read %sv312[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g335 = vector.transfer_read %sv315[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<8x8xf64>
      %t336 = arith.mulf %dv334, %g335 : vector<8x8xf64>
      %t2337 = arith.addf %t336, %dt333 : vector<8x8xf64>
      %fl338 = arith.addf %t2337, %dt324 : vector<8x8xf64>
      vector.transfer_write %fl338, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc339 = arith.constant dense<0.0> : vector<8x8xf64>
      %a340 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b341 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r342 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a340, %b341, %acc339 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a343 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b344 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r345 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a343, %b344, %r342 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r345, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc346 = arith.constant dense<0.0> : vector<8x8xf64>
      %a347 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b348 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r349 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a347, %b348, %acc346 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a350 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b351 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r352 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a350, %b351, %r349 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r352, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv353 = memref.subview %y3[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 256>>
    %sv354 = memref.subview %y3[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 320>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa355 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r356 = vector.transfer_read %sv353[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<8x8xf64>
      %m357 = arith.subf %r356, %fa355 : vector<8x8xf64>
      vector.transfer_write %m357, %sv353[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 256>>
      %r358 = vector.transfer_read %sv354[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<8x8xf64>
      %p359 = arith.addf %r358, %fa355 : vector<8x8xf64>
      vector.transfer_write %p359, %sv354[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 320>>
    }
    gpu.barrier
    %sv360 = memref.subview %v3113[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv361 = memref.subview %v3114[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv362 = memref.subview %g0all[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 320>>
    %sv363 = memref.subview %g1all[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 320>>
    %sv364 = memref.subview %g2all[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 320>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc365 = arith.constant dense<0.0> : vector<8x8xf64>
      %a366 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b367 = vector.transfer_read %sv360[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r368 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a366, %b367, %acc365 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a369 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b370 = vector.transfer_read %sv360[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r371 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a369, %b370, %r368 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g372 = vector.transfer_read %sv363[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<8x8xf64>
      %dt373 = arith.mulf %r371, %g372 : vector<8x8xf64>
      %acc374 = arith.constant dense<0.0> : vector<8x8xf64>
      %a375 = vector.transfer_read %sv360[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b376 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r377 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a375, %b376, %acc374 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a378 = vector.transfer_read %sv360[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b379 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r380 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a378, %b379, %r377 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g381 = vector.transfer_read %sv362[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<8x8xf64>
      %dt382 = arith.mulf %r380, %g381 : vector<8x8xf64>
      %dv383 = vector.transfer_read %sv361[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g384 = vector.transfer_read %sv364[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<8x8xf64>
      %t385 = arith.mulf %dv383, %g384 : vector<8x8xf64>
      %t2386 = arith.addf %t385, %dt382 : vector<8x8xf64>
      %fl387 = arith.addf %t2386, %dt373 : vector<8x8xf64>
      vector.transfer_write %fl387, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc388 = arith.constant dense<0.0> : vector<8x8xf64>
      %a389 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b390 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r391 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a389, %b390, %acc388 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a392 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b393 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r394 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a392, %b393, %r391 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r394, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc395 = arith.constant dense<0.0> : vector<8x8xf64>
      %a396 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b397 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r398 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a396, %b397, %acc395 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a399 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b400 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r401 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a399, %b400, %r398 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r401, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv402 = memref.subview %y3[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 320>>
    %sv403 = memref.subview %y3[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 384>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa404 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r405 = vector.transfer_read %sv402[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<8x8xf64>
      %m406 = arith.subf %r405, %fa404 : vector<8x8xf64>
      vector.transfer_write %m406, %sv402[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 320>>
      %r407 = vector.transfer_read %sv403[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<8x8xf64>
      %p408 = arith.addf %r407, %fa404 : vector<8x8xf64>
      vector.transfer_write %p408, %sv403[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 384>>
    }
    gpu.barrier
    %sv409 = memref.subview %v3113[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv410 = memref.subview %v3114[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv411 = memref.subview %g0all[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 384>>
    %sv412 = memref.subview %g1all[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 384>>
    %sv413 = memref.subview %g2all[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 384>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc414 = arith.constant dense<0.0> : vector<8x8xf64>
      %a415 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b416 = vector.transfer_read %sv409[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r417 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a415, %b416, %acc414 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a418 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b419 = vector.transfer_read %sv409[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r420 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a418, %b419, %r417 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g421 = vector.transfer_read %sv412[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<8x8xf64>
      %dt422 = arith.mulf %r420, %g421 : vector<8x8xf64>
      %acc423 = arith.constant dense<0.0> : vector<8x8xf64>
      %a424 = vector.transfer_read %sv409[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b425 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r426 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a424, %b425, %acc423 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a427 = vector.transfer_read %sv409[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b428 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r429 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a427, %b428, %r426 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g430 = vector.transfer_read %sv411[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<8x8xf64>
      %dt431 = arith.mulf %r429, %g430 : vector<8x8xf64>
      %dv432 = vector.transfer_read %sv410[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g433 = vector.transfer_read %sv413[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<8x8xf64>
      %t434 = arith.mulf %dv432, %g433 : vector<8x8xf64>
      %t2435 = arith.addf %t434, %dt431 : vector<8x8xf64>
      %fl436 = arith.addf %t2435, %dt422 : vector<8x8xf64>
      vector.transfer_write %fl436, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc437 = arith.constant dense<0.0> : vector<8x8xf64>
      %a438 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b439 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r440 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a438, %b439, %acc437 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a441 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b442 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r443 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a441, %b442, %r440 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r443, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc444 = arith.constant dense<0.0> : vector<8x8xf64>
      %a445 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b446 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r447 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a445, %b446, %acc444 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a448 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b449 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r450 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a448, %b449, %r447 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r450, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv451 = memref.subview %y3[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 384>>
    %sv452 = memref.subview %y3[7, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 448>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa453 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r454 = vector.transfer_read %sv451[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<8x8xf64>
      %m455 = arith.subf %r454, %fa453 : vector<8x8xf64>
      vector.transfer_write %m455, %sv451[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 384>>
      %r456 = vector.transfer_read %sv452[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 448>>, vector<8x8xf64>
      %p457 = arith.addf %r456, %fa453 : vector<8x8xf64>
      vector.transfer_write %p457, %sv452[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 448>>
    }
    gpu.barrier
    gpu.return
  }
}
