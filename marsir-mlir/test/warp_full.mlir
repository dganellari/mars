#a = affine_map<(m, n, k) -> (m, k)>
#b = affine_map<(m, n, k) -> (k, n)>
#bt = affine_map<(m, n, k) -> (n, k)>
#c = affine_map<(m, n, k) -> (m, n)>
gpu.module @mir_kernels [#nvvm.target<chip = "sm_90", O = 3>] {
  gpu.func @full(%u0: memref<8x64xf64>, %u1: memref<8x64xf64>, %u2: memref<8x64xf64>, %Btil: memref<8x8xf64>, %Dtil: memref<8x8xf64>, %Dm: memref<8x8xf64>, %W: memref<8x8xf64>, %g00all: memref<8x8x8xf64>, %g01all: memref<8x8x8xf64>, %g02all: memref<8x8x8xf64>, %g10all: memref<8x8x8xf64>, %g11all: memref<8x8x8xf64>, %g12all: memref<8x8x8xf64>, %g20all: memref<8x8x8xf64>, %g21all: memref<8x8x8xf64>, %g22all: memref<8x8x8xf64>, %y3: memref<8x8x8xf64>) workgroup(%interp_all: memref<8x64xf64, #gpu.address_space<workgroup>>, %deriv_all: memref<8x64xf64, #gpu.address_space<workgroup>>, %dt1g: memref<8x8xf64, #gpu.address_space<workgroup>>, %dt2g: memref<8x8xf64, #gpu.address_space<workgroup>>, %flux: memref<8x8xf64, #gpu.address_space<workgroup>>, %tmp: memref<8x8xf64, #gpu.address_space<workgroup>>, %intf: memref<8x8xf64, #gpu.address_space<workgroup>>) kernel {
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
      %b3 = vector.transfer_read %u0[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r4 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a2, %b3, %acc1 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a5 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b6 = vector.transfer_read %u0[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r7 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a5, %b6, %r4 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r7, %interp_all[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc8 = arith.constant dense<0.0> : vector<8x8xf64>
      %a9 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b10 = vector.transfer_read %u0[%c0, %c8], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r11 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a9, %b10, %acc8 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a12 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b13 = vector.transfer_read %u0[%c4, %c8], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r14 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a12, %b13, %r11 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r14, %interp_all[%c0, %c8] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc15 = arith.constant dense<0.0> : vector<8x8xf64>
      %a16 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b17 = vector.transfer_read %u0[%c0, %c16], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r18 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a16, %b17, %acc15 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a19 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b20 = vector.transfer_read %u0[%c4, %c16], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r21 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a19, %b20, %r18 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r21, %interp_all[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc22 = arith.constant dense<0.0> : vector<8x8xf64>
      %a23 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b24 = vector.transfer_read %u0[%c0, %c24], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r25 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a23, %b24, %acc22 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a26 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b27 = vector.transfer_read %u0[%c4, %c24], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r28 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a26, %b27, %r25 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r28, %interp_all[%c0, %c24] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc29 = arith.constant dense<0.0> : vector<8x8xf64>
      %a30 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b31 = vector.transfer_read %u0[%c0, %c32], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r32 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a30, %b31, %acc29 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a33 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b34 = vector.transfer_read %u0[%c4, %c32], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r35 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a33, %b34, %r32 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r35, %interp_all[%c0, %c32] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc36 = arith.constant dense<0.0> : vector<8x8xf64>
      %a37 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b38 = vector.transfer_read %u0[%c0, %c40], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r39 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a37, %b38, %acc36 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a40 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b41 = vector.transfer_read %u0[%c4, %c40], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r42 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a40, %b41, %r39 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r42, %interp_all[%c0, %c40] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc43 = arith.constant dense<0.0> : vector<8x8xf64>
      %a44 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b45 = vector.transfer_read %u0[%c0, %c48], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r46 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a44, %b45, %acc43 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a47 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b48 = vector.transfer_read %u0[%c4, %c48], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r49 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a47, %b48, %r46 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r49, %interp_all[%c0, %c48] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc50 = arith.constant dense<0.0> : vector<8x8xf64>
      %a51 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b52 = vector.transfer_read %u0[%c0, %c56], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r53 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a51, %b52, %acc50 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a54 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b55 = vector.transfer_read %u0[%c4, %c56], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r56 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a54, %b55, %r53 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r56, %interp_all[%c0, %c56] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc57 = arith.constant dense<0.0> : vector<8x8xf64>
      %a58 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b59 = vector.transfer_read %u0[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r60 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a58, %b59, %acc57 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a61 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b62 = vector.transfer_read %u0[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r63 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a61, %b62, %r60 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r63, %deriv_all[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc64 = arith.constant dense<0.0> : vector<8x8xf64>
      %a65 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b66 = vector.transfer_read %u0[%c0, %c8], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r67 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a65, %b66, %acc64 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a68 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b69 = vector.transfer_read %u0[%c4, %c8], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r70 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a68, %b69, %r67 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r70, %deriv_all[%c0, %c8] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc71 = arith.constant dense<0.0> : vector<8x8xf64>
      %a72 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b73 = vector.transfer_read %u0[%c0, %c16], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r74 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a72, %b73, %acc71 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a75 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b76 = vector.transfer_read %u0[%c4, %c16], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r77 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a75, %b76, %r74 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r77, %deriv_all[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc78 = arith.constant dense<0.0> : vector<8x8xf64>
      %a79 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b80 = vector.transfer_read %u0[%c0, %c24], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r81 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a79, %b80, %acc78 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a82 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b83 = vector.transfer_read %u0[%c4, %c24], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r84 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a82, %b83, %r81 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r84, %deriv_all[%c0, %c24] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc85 = arith.constant dense<0.0> : vector<8x8xf64>
      %a86 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b87 = vector.transfer_read %u0[%c0, %c32], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r88 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a86, %b87, %acc85 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a89 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b90 = vector.transfer_read %u0[%c4, %c32], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r91 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a89, %b90, %r88 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r91, %deriv_all[%c0, %c32] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc92 = arith.constant dense<0.0> : vector<8x8xf64>
      %a93 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b94 = vector.transfer_read %u0[%c0, %c40], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r95 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a93, %b94, %acc92 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a96 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b97 = vector.transfer_read %u0[%c4, %c40], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r98 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a96, %b97, %r95 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r98, %deriv_all[%c0, %c40] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc99 = arith.constant dense<0.0> : vector<8x8xf64>
      %a100 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b101 = vector.transfer_read %u0[%c0, %c48], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r102 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a100, %b101, %acc99 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a103 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b104 = vector.transfer_read %u0[%c4, %c48], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r105 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a103, %b104, %r102 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r105, %deriv_all[%c0, %c48] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc106 = arith.constant dense<0.0> : vector<8x8xf64>
      %a107 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b108 = vector.transfer_read %u0[%c0, %c56], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r109 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a107, %b108, %acc106 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a110 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b111 = vector.transfer_read %u0[%c4, %c56], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r112 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a110, %b111, %r109 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r112, %deriv_all[%c0, %c56] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %v3113 = memref.reinterpret_cast %interp_all to offset: [0], sizes: [8, 8, 8], strides: [64, 8, 1] : memref<8x64xf64, #gpu.address_space<workgroup>> to memref<8x8x8xf64, #gpu.address_space<workgroup>>
    %v3114 = memref.reinterpret_cast %deriv_all to offset: [0], sizes: [8, 8, 8], strides: [64, 8, 1] : memref<8x64xf64, #gpu.address_space<workgroup>> to memref<8x8x8xf64, #gpu.address_space<workgroup>>
    %sv115 = memref.subview %v3113[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>
    %sv116 = memref.subview %v3114[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>
    %sv117 = memref.subview %g00all[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 0>>
    %sv118 = memref.subview %g01all[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 0>>
    %sv119 = memref.subview %g02all[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 0>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc120 = arith.constant dense<0.0> : vector<8x8xf64>
      %a121 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b122 = vector.transfer_read %sv115[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r123 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a121, %b122, %acc120 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a124 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b125 = vector.transfer_read %sv115[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r126 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a124, %b125, %r123 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g127 = vector.transfer_read %sv118[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>>, vector<8x8xf64>
      %sc128 = arith.mulf %r126, %g127 : vector<8x8xf64>
      vector.transfer_write %sc128, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc129 = arith.constant dense<0.0> : vector<8x8xf64>
      %a130 = vector.transfer_read %sv115[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b131 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r132 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a130, %b131, %acc129 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a133 = vector.transfer_read %sv115[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b134 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r135 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a133, %b134, %r132 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g136 = vector.transfer_read %sv117[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>>, vector<8x8xf64>
      %sc137 = arith.mulf %r135, %g136 : vector<8x8xf64>
      vector.transfer_write %sc137, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d138 = vector.transfer_read %sv116[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2139 = vector.transfer_read %sv119[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>>, vector<8x8xf64>
      %t140 = arith.mulf %d138, %g2139 : vector<8x8xf64>
      %a141 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2142 = arith.addf %t140, %a141 : vector<8x8xf64>
      %b143 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl144 = arith.addf %t2142, %b143 : vector<8x8xf64>
      vector.transfer_write %fl144, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc145 = arith.constant dense<0.0> : vector<8x8xf64>
      %a146 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b147 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r148 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a146, %b147, %acc145 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a149 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b150 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r151 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a149, %b150, %r148 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r151, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc152 = arith.constant dense<0.0> : vector<8x8xf64>
      %a153 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b154 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r155 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a153, %b154, %acc152 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a156 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b157 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r158 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a156, %b157, %r155 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r158, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv159 = memref.subview %y3[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 0>>
    %sv160 = memref.subview %y3[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 64>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa161 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r162 = vector.transfer_read %sv159[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>>, vector<8x8xf64>
      %m163 = arith.subf %r162, %fa161 : vector<8x8xf64>
      vector.transfer_write %m163, %sv159[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 0>>
      %fa164 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r165 = vector.transfer_read %sv160[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<8x8xf64>
      %p166 = arith.addf %r165, %fa164 : vector<8x8xf64>
      vector.transfer_write %p166, %sv160[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 64>>
    }
    gpu.barrier
    %sv167 = memref.subview %v3113[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv168 = memref.subview %v3114[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv169 = memref.subview %g00all[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 64>>
    %sv170 = memref.subview %g01all[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 64>>
    %sv171 = memref.subview %g02all[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 64>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc172 = arith.constant dense<0.0> : vector<8x8xf64>
      %a173 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b174 = vector.transfer_read %sv167[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r175 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a173, %b174, %acc172 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a176 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b177 = vector.transfer_read %sv167[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r178 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a176, %b177, %r175 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g179 = vector.transfer_read %sv170[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<8x8xf64>
      %sc180 = arith.mulf %r178, %g179 : vector<8x8xf64>
      vector.transfer_write %sc180, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc181 = arith.constant dense<0.0> : vector<8x8xf64>
      %a182 = vector.transfer_read %sv167[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b183 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r184 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a182, %b183, %acc181 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a185 = vector.transfer_read %sv167[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b186 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r187 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a185, %b186, %r184 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g188 = vector.transfer_read %sv169[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<8x8xf64>
      %sc189 = arith.mulf %r187, %g188 : vector<8x8xf64>
      vector.transfer_write %sc189, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d190 = vector.transfer_read %sv168[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2191 = vector.transfer_read %sv171[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<8x8xf64>
      %t192 = arith.mulf %d190, %g2191 : vector<8x8xf64>
      %a193 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2194 = arith.addf %t192, %a193 : vector<8x8xf64>
      %b195 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl196 = arith.addf %t2194, %b195 : vector<8x8xf64>
      vector.transfer_write %fl196, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc197 = arith.constant dense<0.0> : vector<8x8xf64>
      %a198 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b199 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r200 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a198, %b199, %acc197 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a201 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b202 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r203 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a201, %b202, %r200 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r203, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc204 = arith.constant dense<0.0> : vector<8x8xf64>
      %a205 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b206 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r207 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a205, %b206, %acc204 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a208 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b209 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r210 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a208, %b209, %r207 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r210, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv211 = memref.subview %y3[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 64>>
    %sv212 = memref.subview %y3[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 128>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa213 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r214 = vector.transfer_read %sv211[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<8x8xf64>
      %m215 = arith.subf %r214, %fa213 : vector<8x8xf64>
      vector.transfer_write %m215, %sv211[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 64>>
      %fa216 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r217 = vector.transfer_read %sv212[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<8x8xf64>
      %p218 = arith.addf %r217, %fa216 : vector<8x8xf64>
      vector.transfer_write %p218, %sv212[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 128>>
    }
    gpu.barrier
    %sv219 = memref.subview %v3113[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv220 = memref.subview %v3114[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv221 = memref.subview %g00all[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 128>>
    %sv222 = memref.subview %g01all[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 128>>
    %sv223 = memref.subview %g02all[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 128>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc224 = arith.constant dense<0.0> : vector<8x8xf64>
      %a225 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b226 = vector.transfer_read %sv219[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r227 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a225, %b226, %acc224 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a228 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b229 = vector.transfer_read %sv219[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r230 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a228, %b229, %r227 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g231 = vector.transfer_read %sv222[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<8x8xf64>
      %sc232 = arith.mulf %r230, %g231 : vector<8x8xf64>
      vector.transfer_write %sc232, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc233 = arith.constant dense<0.0> : vector<8x8xf64>
      %a234 = vector.transfer_read %sv219[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b235 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r236 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a234, %b235, %acc233 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a237 = vector.transfer_read %sv219[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b238 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r239 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a237, %b238, %r236 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g240 = vector.transfer_read %sv221[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<8x8xf64>
      %sc241 = arith.mulf %r239, %g240 : vector<8x8xf64>
      vector.transfer_write %sc241, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d242 = vector.transfer_read %sv220[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2243 = vector.transfer_read %sv223[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<8x8xf64>
      %t244 = arith.mulf %d242, %g2243 : vector<8x8xf64>
      %a245 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2246 = arith.addf %t244, %a245 : vector<8x8xf64>
      %b247 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl248 = arith.addf %t2246, %b247 : vector<8x8xf64>
      vector.transfer_write %fl248, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc249 = arith.constant dense<0.0> : vector<8x8xf64>
      %a250 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b251 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r252 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a250, %b251, %acc249 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a253 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b254 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r255 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a253, %b254, %r252 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r255, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc256 = arith.constant dense<0.0> : vector<8x8xf64>
      %a257 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b258 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r259 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a257, %b258, %acc256 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a260 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b261 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r262 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a260, %b261, %r259 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r262, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv263 = memref.subview %y3[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 128>>
    %sv264 = memref.subview %y3[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 192>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa265 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r266 = vector.transfer_read %sv263[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<8x8xf64>
      %m267 = arith.subf %r266, %fa265 : vector<8x8xf64>
      vector.transfer_write %m267, %sv263[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 128>>
      %fa268 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r269 = vector.transfer_read %sv264[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<8x8xf64>
      %p270 = arith.addf %r269, %fa268 : vector<8x8xf64>
      vector.transfer_write %p270, %sv264[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 192>>
    }
    gpu.barrier
    %sv271 = memref.subview %v3113[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv272 = memref.subview %v3114[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv273 = memref.subview %g00all[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 192>>
    %sv274 = memref.subview %g01all[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 192>>
    %sv275 = memref.subview %g02all[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 192>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc276 = arith.constant dense<0.0> : vector<8x8xf64>
      %a277 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b278 = vector.transfer_read %sv271[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r279 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a277, %b278, %acc276 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a280 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b281 = vector.transfer_read %sv271[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r282 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a280, %b281, %r279 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g283 = vector.transfer_read %sv274[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<8x8xf64>
      %sc284 = arith.mulf %r282, %g283 : vector<8x8xf64>
      vector.transfer_write %sc284, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc285 = arith.constant dense<0.0> : vector<8x8xf64>
      %a286 = vector.transfer_read %sv271[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b287 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r288 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a286, %b287, %acc285 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a289 = vector.transfer_read %sv271[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b290 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r291 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a289, %b290, %r288 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g292 = vector.transfer_read %sv273[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<8x8xf64>
      %sc293 = arith.mulf %r291, %g292 : vector<8x8xf64>
      vector.transfer_write %sc293, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d294 = vector.transfer_read %sv272[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2295 = vector.transfer_read %sv275[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<8x8xf64>
      %t296 = arith.mulf %d294, %g2295 : vector<8x8xf64>
      %a297 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2298 = arith.addf %t296, %a297 : vector<8x8xf64>
      %b299 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl300 = arith.addf %t2298, %b299 : vector<8x8xf64>
      vector.transfer_write %fl300, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc301 = arith.constant dense<0.0> : vector<8x8xf64>
      %a302 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b303 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r304 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a302, %b303, %acc301 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a305 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b306 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r307 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a305, %b306, %r304 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r307, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc308 = arith.constant dense<0.0> : vector<8x8xf64>
      %a309 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b310 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r311 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a309, %b310, %acc308 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a312 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b313 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r314 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a312, %b313, %r311 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r314, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv315 = memref.subview %y3[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 192>>
    %sv316 = memref.subview %y3[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 256>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa317 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r318 = vector.transfer_read %sv315[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<8x8xf64>
      %m319 = arith.subf %r318, %fa317 : vector<8x8xf64>
      vector.transfer_write %m319, %sv315[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 192>>
      %fa320 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r321 = vector.transfer_read %sv316[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<8x8xf64>
      %p322 = arith.addf %r321, %fa320 : vector<8x8xf64>
      vector.transfer_write %p322, %sv316[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 256>>
    }
    gpu.barrier
    %sv323 = memref.subview %v3113[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv324 = memref.subview %v3114[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv325 = memref.subview %g00all[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 256>>
    %sv326 = memref.subview %g01all[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 256>>
    %sv327 = memref.subview %g02all[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 256>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc328 = arith.constant dense<0.0> : vector<8x8xf64>
      %a329 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b330 = vector.transfer_read %sv323[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r331 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a329, %b330, %acc328 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a332 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b333 = vector.transfer_read %sv323[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r334 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a332, %b333, %r331 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g335 = vector.transfer_read %sv326[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<8x8xf64>
      %sc336 = arith.mulf %r334, %g335 : vector<8x8xf64>
      vector.transfer_write %sc336, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc337 = arith.constant dense<0.0> : vector<8x8xf64>
      %a338 = vector.transfer_read %sv323[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b339 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r340 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a338, %b339, %acc337 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a341 = vector.transfer_read %sv323[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b342 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r343 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a341, %b342, %r340 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g344 = vector.transfer_read %sv325[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<8x8xf64>
      %sc345 = arith.mulf %r343, %g344 : vector<8x8xf64>
      vector.transfer_write %sc345, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d346 = vector.transfer_read %sv324[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2347 = vector.transfer_read %sv327[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<8x8xf64>
      %t348 = arith.mulf %d346, %g2347 : vector<8x8xf64>
      %a349 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2350 = arith.addf %t348, %a349 : vector<8x8xf64>
      %b351 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl352 = arith.addf %t2350, %b351 : vector<8x8xf64>
      vector.transfer_write %fl352, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc353 = arith.constant dense<0.0> : vector<8x8xf64>
      %a354 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b355 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r356 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a354, %b355, %acc353 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a357 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b358 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r359 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a357, %b358, %r356 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r359, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc360 = arith.constant dense<0.0> : vector<8x8xf64>
      %a361 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b362 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r363 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a361, %b362, %acc360 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a364 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b365 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r366 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a364, %b365, %r363 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r366, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv367 = memref.subview %y3[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 256>>
    %sv368 = memref.subview %y3[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 320>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa369 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r370 = vector.transfer_read %sv367[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<8x8xf64>
      %m371 = arith.subf %r370, %fa369 : vector<8x8xf64>
      vector.transfer_write %m371, %sv367[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 256>>
      %fa372 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r373 = vector.transfer_read %sv368[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<8x8xf64>
      %p374 = arith.addf %r373, %fa372 : vector<8x8xf64>
      vector.transfer_write %p374, %sv368[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 320>>
    }
    gpu.barrier
    %sv375 = memref.subview %v3113[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv376 = memref.subview %v3114[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv377 = memref.subview %g00all[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 320>>
    %sv378 = memref.subview %g01all[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 320>>
    %sv379 = memref.subview %g02all[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 320>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc380 = arith.constant dense<0.0> : vector<8x8xf64>
      %a381 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b382 = vector.transfer_read %sv375[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r383 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a381, %b382, %acc380 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a384 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b385 = vector.transfer_read %sv375[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r386 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a384, %b385, %r383 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g387 = vector.transfer_read %sv378[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<8x8xf64>
      %sc388 = arith.mulf %r386, %g387 : vector<8x8xf64>
      vector.transfer_write %sc388, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc389 = arith.constant dense<0.0> : vector<8x8xf64>
      %a390 = vector.transfer_read %sv375[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b391 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r392 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a390, %b391, %acc389 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a393 = vector.transfer_read %sv375[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b394 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r395 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a393, %b394, %r392 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g396 = vector.transfer_read %sv377[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<8x8xf64>
      %sc397 = arith.mulf %r395, %g396 : vector<8x8xf64>
      vector.transfer_write %sc397, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d398 = vector.transfer_read %sv376[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2399 = vector.transfer_read %sv379[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<8x8xf64>
      %t400 = arith.mulf %d398, %g2399 : vector<8x8xf64>
      %a401 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2402 = arith.addf %t400, %a401 : vector<8x8xf64>
      %b403 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl404 = arith.addf %t2402, %b403 : vector<8x8xf64>
      vector.transfer_write %fl404, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc405 = arith.constant dense<0.0> : vector<8x8xf64>
      %a406 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b407 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r408 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a406, %b407, %acc405 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a409 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b410 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r411 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a409, %b410, %r408 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r411, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc412 = arith.constant dense<0.0> : vector<8x8xf64>
      %a413 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b414 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r415 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a413, %b414, %acc412 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a416 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b417 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r418 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a416, %b417, %r415 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r418, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv419 = memref.subview %y3[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 320>>
    %sv420 = memref.subview %y3[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 384>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa421 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r422 = vector.transfer_read %sv419[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<8x8xf64>
      %m423 = arith.subf %r422, %fa421 : vector<8x8xf64>
      vector.transfer_write %m423, %sv419[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 320>>
      %fa424 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r425 = vector.transfer_read %sv420[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<8x8xf64>
      %p426 = arith.addf %r425, %fa424 : vector<8x8xf64>
      vector.transfer_write %p426, %sv420[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 384>>
    }
    gpu.barrier
    %sv427 = memref.subview %v3113[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv428 = memref.subview %v3114[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv429 = memref.subview %g00all[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 384>>
    %sv430 = memref.subview %g01all[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 384>>
    %sv431 = memref.subview %g02all[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 384>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc432 = arith.constant dense<0.0> : vector<8x8xf64>
      %a433 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b434 = vector.transfer_read %sv427[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r435 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a433, %b434, %acc432 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a436 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b437 = vector.transfer_read %sv427[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r438 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a436, %b437, %r435 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g439 = vector.transfer_read %sv430[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<8x8xf64>
      %sc440 = arith.mulf %r438, %g439 : vector<8x8xf64>
      vector.transfer_write %sc440, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc441 = arith.constant dense<0.0> : vector<8x8xf64>
      %a442 = vector.transfer_read %sv427[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b443 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r444 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a442, %b443, %acc441 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a445 = vector.transfer_read %sv427[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b446 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r447 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a445, %b446, %r444 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g448 = vector.transfer_read %sv429[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<8x8xf64>
      %sc449 = arith.mulf %r447, %g448 : vector<8x8xf64>
      vector.transfer_write %sc449, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d450 = vector.transfer_read %sv428[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2451 = vector.transfer_read %sv431[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<8x8xf64>
      %t452 = arith.mulf %d450, %g2451 : vector<8x8xf64>
      %a453 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2454 = arith.addf %t452, %a453 : vector<8x8xf64>
      %b455 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl456 = arith.addf %t2454, %b455 : vector<8x8xf64>
      vector.transfer_write %fl456, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc457 = arith.constant dense<0.0> : vector<8x8xf64>
      %a458 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b459 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r460 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a458, %b459, %acc457 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a461 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b462 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r463 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a461, %b462, %r460 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r463, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc464 = arith.constant dense<0.0> : vector<8x8xf64>
      %a465 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b466 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r467 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a465, %b466, %acc464 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a468 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b469 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r470 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a468, %b469, %r467 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r470, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv471 = memref.subview %y3[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 384>>
    %sv472 = memref.subview %y3[7, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 448>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa473 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r474 = vector.transfer_read %sv471[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<8x8xf64>
      %m475 = arith.subf %r474, %fa473 : vector<8x8xf64>
      vector.transfer_write %m475, %sv471[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 384>>
      %fa476 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r477 = vector.transfer_read %sv472[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 448>>, vector<8x8xf64>
      %p478 = arith.addf %r477, %fa476 : vector<8x8xf64>
      vector.transfer_write %p478, %sv472[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 448>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc479 = arith.constant dense<0.0> : vector<8x8xf64>
      %a480 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b481 = vector.transfer_read %u1[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r482 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a480, %b481, %acc479 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a483 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b484 = vector.transfer_read %u1[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r485 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a483, %b484, %r482 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r485, %interp_all[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc486 = arith.constant dense<0.0> : vector<8x8xf64>
      %a487 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b488 = vector.transfer_read %u1[%c0, %c8], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r489 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a487, %b488, %acc486 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a490 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b491 = vector.transfer_read %u1[%c4, %c8], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r492 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a490, %b491, %r489 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r492, %interp_all[%c0, %c8] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc493 = arith.constant dense<0.0> : vector<8x8xf64>
      %a494 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b495 = vector.transfer_read %u1[%c0, %c16], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r496 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a494, %b495, %acc493 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a497 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b498 = vector.transfer_read %u1[%c4, %c16], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r499 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a497, %b498, %r496 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r499, %interp_all[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc500 = arith.constant dense<0.0> : vector<8x8xf64>
      %a501 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b502 = vector.transfer_read %u1[%c0, %c24], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r503 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a501, %b502, %acc500 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a504 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b505 = vector.transfer_read %u1[%c4, %c24], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r506 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a504, %b505, %r503 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r506, %interp_all[%c0, %c24] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc507 = arith.constant dense<0.0> : vector<8x8xf64>
      %a508 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b509 = vector.transfer_read %u1[%c0, %c32], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r510 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a508, %b509, %acc507 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a511 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b512 = vector.transfer_read %u1[%c4, %c32], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r513 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a511, %b512, %r510 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r513, %interp_all[%c0, %c32] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc514 = arith.constant dense<0.0> : vector<8x8xf64>
      %a515 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b516 = vector.transfer_read %u1[%c0, %c40], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r517 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a515, %b516, %acc514 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a518 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b519 = vector.transfer_read %u1[%c4, %c40], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r520 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a518, %b519, %r517 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r520, %interp_all[%c0, %c40] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc521 = arith.constant dense<0.0> : vector<8x8xf64>
      %a522 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b523 = vector.transfer_read %u1[%c0, %c48], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r524 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a522, %b523, %acc521 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a525 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b526 = vector.transfer_read %u1[%c4, %c48], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r527 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a525, %b526, %r524 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r527, %interp_all[%c0, %c48] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc528 = arith.constant dense<0.0> : vector<8x8xf64>
      %a529 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b530 = vector.transfer_read %u1[%c0, %c56], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r531 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a529, %b530, %acc528 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a532 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b533 = vector.transfer_read %u1[%c4, %c56], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r534 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a532, %b533, %r531 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r534, %interp_all[%c0, %c56] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc535 = arith.constant dense<0.0> : vector<8x8xf64>
      %a536 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b537 = vector.transfer_read %u1[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r538 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a536, %b537, %acc535 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a539 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b540 = vector.transfer_read %u1[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r541 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a539, %b540, %r538 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r541, %deriv_all[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc542 = arith.constant dense<0.0> : vector<8x8xf64>
      %a543 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b544 = vector.transfer_read %u1[%c0, %c8], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r545 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a543, %b544, %acc542 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a546 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b547 = vector.transfer_read %u1[%c4, %c8], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r548 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a546, %b547, %r545 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r548, %deriv_all[%c0, %c8] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc549 = arith.constant dense<0.0> : vector<8x8xf64>
      %a550 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b551 = vector.transfer_read %u1[%c0, %c16], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r552 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a550, %b551, %acc549 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a553 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b554 = vector.transfer_read %u1[%c4, %c16], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r555 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a553, %b554, %r552 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r555, %deriv_all[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc556 = arith.constant dense<0.0> : vector<8x8xf64>
      %a557 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b558 = vector.transfer_read %u1[%c0, %c24], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r559 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a557, %b558, %acc556 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a560 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b561 = vector.transfer_read %u1[%c4, %c24], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r562 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a560, %b561, %r559 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r562, %deriv_all[%c0, %c24] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc563 = arith.constant dense<0.0> : vector<8x8xf64>
      %a564 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b565 = vector.transfer_read %u1[%c0, %c32], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r566 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a564, %b565, %acc563 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a567 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b568 = vector.transfer_read %u1[%c4, %c32], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r569 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a567, %b568, %r566 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r569, %deriv_all[%c0, %c32] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc570 = arith.constant dense<0.0> : vector<8x8xf64>
      %a571 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b572 = vector.transfer_read %u1[%c0, %c40], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r573 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a571, %b572, %acc570 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a574 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b575 = vector.transfer_read %u1[%c4, %c40], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r576 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a574, %b575, %r573 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r576, %deriv_all[%c0, %c40] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc577 = arith.constant dense<0.0> : vector<8x8xf64>
      %a578 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b579 = vector.transfer_read %u1[%c0, %c48], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r580 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a578, %b579, %acc577 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a581 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b582 = vector.transfer_read %u1[%c4, %c48], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r583 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a581, %b582, %r580 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r583, %deriv_all[%c0, %c48] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc584 = arith.constant dense<0.0> : vector<8x8xf64>
      %a585 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b586 = vector.transfer_read %u1[%c0, %c56], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r587 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a585, %b586, %acc584 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a588 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b589 = vector.transfer_read %u1[%c4, %c56], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r590 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a588, %b589, %r587 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r590, %deriv_all[%c0, %c56] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %v3591 = memref.reinterpret_cast %interp_all to offset: [0], sizes: [8, 8, 8], strides: [64, 8, 1] : memref<8x64xf64, #gpu.address_space<workgroup>> to memref<8x8x8xf64, #gpu.address_space<workgroup>>
    %v3592 = memref.reinterpret_cast %deriv_all to offset: [0], sizes: [8, 8, 8], strides: [64, 8, 1] : memref<8x64xf64, #gpu.address_space<workgroup>> to memref<8x8x8xf64, #gpu.address_space<workgroup>>
    %sv593 = memref.subview %v3591[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>
    %sv594 = memref.subview %v3592[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>
    %sv595 = memref.subview %g10all[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 0>>
    %sv596 = memref.subview %g11all[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 0>>
    %sv597 = memref.subview %g12all[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 0>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc598 = arith.constant dense<0.0> : vector<8x8xf64>
      %a599 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b600 = vector.transfer_read %sv593[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r601 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a599, %b600, %acc598 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a602 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b603 = vector.transfer_read %sv593[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r604 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a602, %b603, %r601 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g605 = vector.transfer_read %sv596[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>>, vector<8x8xf64>
      %sc606 = arith.mulf %r604, %g605 : vector<8x8xf64>
      vector.transfer_write %sc606, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc607 = arith.constant dense<0.0> : vector<8x8xf64>
      %a608 = vector.transfer_read %sv593[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b609 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r610 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a608, %b609, %acc607 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a611 = vector.transfer_read %sv593[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b612 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r613 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a611, %b612, %r610 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g614 = vector.transfer_read %sv595[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>>, vector<8x8xf64>
      %sc615 = arith.mulf %r613, %g614 : vector<8x8xf64>
      vector.transfer_write %sc615, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d616 = vector.transfer_read %sv594[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2617 = vector.transfer_read %sv597[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>>, vector<8x8xf64>
      %t618 = arith.mulf %d616, %g2617 : vector<8x8xf64>
      %a619 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2620 = arith.addf %t618, %a619 : vector<8x8xf64>
      %b621 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl622 = arith.addf %t2620, %b621 : vector<8x8xf64>
      vector.transfer_write %fl622, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc623 = arith.constant dense<0.0> : vector<8x8xf64>
      %a624 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b625 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r626 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a624, %b625, %acc623 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a627 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b628 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r629 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a627, %b628, %r626 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r629, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc630 = arith.constant dense<0.0> : vector<8x8xf64>
      %a631 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b632 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r633 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a631, %b632, %acc630 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a634 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b635 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r636 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a634, %b635, %r633 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r636, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv637 = memref.subview %y3[0, 0, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 0>>
    %sv638 = memref.subview %y3[0, 1, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 8>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa639 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r640 = vector.transfer_read %sv637[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 0>>, vector<8x8xf64>
      %m641 = arith.subf %r640, %fa639 : vector<8x8xf64>
      vector.transfer_write %m641, %sv637[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: 0>>
      %fa642 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r643 = vector.transfer_read %sv638[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 8>>, vector<8x8xf64>
      %p644 = arith.addf %r643, %fa642 : vector<8x8xf64>
      vector.transfer_write %p644, %sv638[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: 8>>
    }
    gpu.barrier
    %sv645 = memref.subview %v3591[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv646 = memref.subview %v3592[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv647 = memref.subview %g10all[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 64>>
    %sv648 = memref.subview %g11all[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 64>>
    %sv649 = memref.subview %g12all[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 64>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc650 = arith.constant dense<0.0> : vector<8x8xf64>
      %a651 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b652 = vector.transfer_read %sv645[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r653 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a651, %b652, %acc650 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a654 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b655 = vector.transfer_read %sv645[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r656 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a654, %b655, %r653 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g657 = vector.transfer_read %sv648[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<8x8xf64>
      %sc658 = arith.mulf %r656, %g657 : vector<8x8xf64>
      vector.transfer_write %sc658, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc659 = arith.constant dense<0.0> : vector<8x8xf64>
      %a660 = vector.transfer_read %sv645[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b661 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r662 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a660, %b661, %acc659 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a663 = vector.transfer_read %sv645[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b664 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r665 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a663, %b664, %r662 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g666 = vector.transfer_read %sv647[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<8x8xf64>
      %sc667 = arith.mulf %r665, %g666 : vector<8x8xf64>
      vector.transfer_write %sc667, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d668 = vector.transfer_read %sv646[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2669 = vector.transfer_read %sv649[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<8x8xf64>
      %t670 = arith.mulf %d668, %g2669 : vector<8x8xf64>
      %a671 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2672 = arith.addf %t670, %a671 : vector<8x8xf64>
      %b673 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl674 = arith.addf %t2672, %b673 : vector<8x8xf64>
      vector.transfer_write %fl674, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc675 = arith.constant dense<0.0> : vector<8x8xf64>
      %a676 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b677 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r678 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a676, %b677, %acc675 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a679 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b680 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r681 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a679, %b680, %r678 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r681, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc682 = arith.constant dense<0.0> : vector<8x8xf64>
      %a683 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b684 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r685 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a683, %b684, %acc682 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a686 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b687 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r688 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a686, %b687, %r685 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r688, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv689 = memref.subview %y3[0, 1, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 8>>
    %sv690 = memref.subview %y3[0, 2, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 16>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa691 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r692 = vector.transfer_read %sv689[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 8>>, vector<8x8xf64>
      %m693 = arith.subf %r692, %fa691 : vector<8x8xf64>
      vector.transfer_write %m693, %sv689[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: 8>>
      %fa694 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r695 = vector.transfer_read %sv690[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 16>>, vector<8x8xf64>
      %p696 = arith.addf %r695, %fa694 : vector<8x8xf64>
      vector.transfer_write %p696, %sv690[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: 16>>
    }
    gpu.barrier
    %sv697 = memref.subview %v3591[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv698 = memref.subview %v3592[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv699 = memref.subview %g10all[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 128>>
    %sv700 = memref.subview %g11all[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 128>>
    %sv701 = memref.subview %g12all[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 128>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc702 = arith.constant dense<0.0> : vector<8x8xf64>
      %a703 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b704 = vector.transfer_read %sv697[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r705 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a703, %b704, %acc702 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a706 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b707 = vector.transfer_read %sv697[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r708 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a706, %b707, %r705 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g709 = vector.transfer_read %sv700[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<8x8xf64>
      %sc710 = arith.mulf %r708, %g709 : vector<8x8xf64>
      vector.transfer_write %sc710, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc711 = arith.constant dense<0.0> : vector<8x8xf64>
      %a712 = vector.transfer_read %sv697[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b713 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r714 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a712, %b713, %acc711 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a715 = vector.transfer_read %sv697[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b716 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r717 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a715, %b716, %r714 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g718 = vector.transfer_read %sv699[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<8x8xf64>
      %sc719 = arith.mulf %r717, %g718 : vector<8x8xf64>
      vector.transfer_write %sc719, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d720 = vector.transfer_read %sv698[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2721 = vector.transfer_read %sv701[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<8x8xf64>
      %t722 = arith.mulf %d720, %g2721 : vector<8x8xf64>
      %a723 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2724 = arith.addf %t722, %a723 : vector<8x8xf64>
      %b725 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl726 = arith.addf %t2724, %b725 : vector<8x8xf64>
      vector.transfer_write %fl726, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc727 = arith.constant dense<0.0> : vector<8x8xf64>
      %a728 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b729 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r730 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a728, %b729, %acc727 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a731 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b732 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r733 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a731, %b732, %r730 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r733, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc734 = arith.constant dense<0.0> : vector<8x8xf64>
      %a735 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b736 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r737 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a735, %b736, %acc734 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a738 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b739 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r740 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a738, %b739, %r737 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r740, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv741 = memref.subview %y3[0, 2, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 16>>
    %sv742 = memref.subview %y3[0, 3, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 24>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa743 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r744 = vector.transfer_read %sv741[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 16>>, vector<8x8xf64>
      %m745 = arith.subf %r744, %fa743 : vector<8x8xf64>
      vector.transfer_write %m745, %sv741[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: 16>>
      %fa746 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r747 = vector.transfer_read %sv742[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 24>>, vector<8x8xf64>
      %p748 = arith.addf %r747, %fa746 : vector<8x8xf64>
      vector.transfer_write %p748, %sv742[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: 24>>
    }
    gpu.barrier
    %sv749 = memref.subview %v3591[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv750 = memref.subview %v3592[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv751 = memref.subview %g10all[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 192>>
    %sv752 = memref.subview %g11all[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 192>>
    %sv753 = memref.subview %g12all[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 192>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc754 = arith.constant dense<0.0> : vector<8x8xf64>
      %a755 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b756 = vector.transfer_read %sv749[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r757 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a755, %b756, %acc754 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a758 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b759 = vector.transfer_read %sv749[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r760 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a758, %b759, %r757 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g761 = vector.transfer_read %sv752[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<8x8xf64>
      %sc762 = arith.mulf %r760, %g761 : vector<8x8xf64>
      vector.transfer_write %sc762, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc763 = arith.constant dense<0.0> : vector<8x8xf64>
      %a764 = vector.transfer_read %sv749[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b765 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r766 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a764, %b765, %acc763 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a767 = vector.transfer_read %sv749[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b768 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r769 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a767, %b768, %r766 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g770 = vector.transfer_read %sv751[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<8x8xf64>
      %sc771 = arith.mulf %r769, %g770 : vector<8x8xf64>
      vector.transfer_write %sc771, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d772 = vector.transfer_read %sv750[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2773 = vector.transfer_read %sv753[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<8x8xf64>
      %t774 = arith.mulf %d772, %g2773 : vector<8x8xf64>
      %a775 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2776 = arith.addf %t774, %a775 : vector<8x8xf64>
      %b777 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl778 = arith.addf %t2776, %b777 : vector<8x8xf64>
      vector.transfer_write %fl778, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc779 = arith.constant dense<0.0> : vector<8x8xf64>
      %a780 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b781 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r782 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a780, %b781, %acc779 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a783 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b784 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r785 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a783, %b784, %r782 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r785, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc786 = arith.constant dense<0.0> : vector<8x8xf64>
      %a787 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b788 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r789 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a787, %b788, %acc786 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a790 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b791 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r792 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a790, %b791, %r789 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r792, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv793 = memref.subview %y3[0, 3, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 24>>
    %sv794 = memref.subview %y3[0, 4, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 32>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa795 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r796 = vector.transfer_read %sv793[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 24>>, vector<8x8xf64>
      %m797 = arith.subf %r796, %fa795 : vector<8x8xf64>
      vector.transfer_write %m797, %sv793[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: 24>>
      %fa798 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r799 = vector.transfer_read %sv794[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 32>>, vector<8x8xf64>
      %p800 = arith.addf %r799, %fa798 : vector<8x8xf64>
      vector.transfer_write %p800, %sv794[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: 32>>
    }
    gpu.barrier
    %sv801 = memref.subview %v3591[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv802 = memref.subview %v3592[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv803 = memref.subview %g10all[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 256>>
    %sv804 = memref.subview %g11all[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 256>>
    %sv805 = memref.subview %g12all[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 256>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc806 = arith.constant dense<0.0> : vector<8x8xf64>
      %a807 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b808 = vector.transfer_read %sv801[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r809 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a807, %b808, %acc806 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a810 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b811 = vector.transfer_read %sv801[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r812 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a810, %b811, %r809 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g813 = vector.transfer_read %sv804[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<8x8xf64>
      %sc814 = arith.mulf %r812, %g813 : vector<8x8xf64>
      vector.transfer_write %sc814, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc815 = arith.constant dense<0.0> : vector<8x8xf64>
      %a816 = vector.transfer_read %sv801[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b817 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r818 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a816, %b817, %acc815 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a819 = vector.transfer_read %sv801[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b820 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r821 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a819, %b820, %r818 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g822 = vector.transfer_read %sv803[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<8x8xf64>
      %sc823 = arith.mulf %r821, %g822 : vector<8x8xf64>
      vector.transfer_write %sc823, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d824 = vector.transfer_read %sv802[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2825 = vector.transfer_read %sv805[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<8x8xf64>
      %t826 = arith.mulf %d824, %g2825 : vector<8x8xf64>
      %a827 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2828 = arith.addf %t826, %a827 : vector<8x8xf64>
      %b829 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl830 = arith.addf %t2828, %b829 : vector<8x8xf64>
      vector.transfer_write %fl830, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc831 = arith.constant dense<0.0> : vector<8x8xf64>
      %a832 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b833 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r834 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a832, %b833, %acc831 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a835 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b836 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r837 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a835, %b836, %r834 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r837, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc838 = arith.constant dense<0.0> : vector<8x8xf64>
      %a839 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b840 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r841 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a839, %b840, %acc838 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a842 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b843 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r844 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a842, %b843, %r841 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r844, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv845 = memref.subview %y3[0, 4, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 32>>
    %sv846 = memref.subview %y3[0, 5, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 40>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa847 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r848 = vector.transfer_read %sv845[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 32>>, vector<8x8xf64>
      %m849 = arith.subf %r848, %fa847 : vector<8x8xf64>
      vector.transfer_write %m849, %sv845[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: 32>>
      %fa850 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r851 = vector.transfer_read %sv846[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 40>>, vector<8x8xf64>
      %p852 = arith.addf %r851, %fa850 : vector<8x8xf64>
      vector.transfer_write %p852, %sv846[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: 40>>
    }
    gpu.barrier
    %sv853 = memref.subview %v3591[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv854 = memref.subview %v3592[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv855 = memref.subview %g10all[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 320>>
    %sv856 = memref.subview %g11all[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 320>>
    %sv857 = memref.subview %g12all[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 320>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc858 = arith.constant dense<0.0> : vector<8x8xf64>
      %a859 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b860 = vector.transfer_read %sv853[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r861 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a859, %b860, %acc858 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a862 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b863 = vector.transfer_read %sv853[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r864 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a862, %b863, %r861 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g865 = vector.transfer_read %sv856[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<8x8xf64>
      %sc866 = arith.mulf %r864, %g865 : vector<8x8xf64>
      vector.transfer_write %sc866, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc867 = arith.constant dense<0.0> : vector<8x8xf64>
      %a868 = vector.transfer_read %sv853[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b869 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r870 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a868, %b869, %acc867 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a871 = vector.transfer_read %sv853[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b872 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r873 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a871, %b872, %r870 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g874 = vector.transfer_read %sv855[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<8x8xf64>
      %sc875 = arith.mulf %r873, %g874 : vector<8x8xf64>
      vector.transfer_write %sc875, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d876 = vector.transfer_read %sv854[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2877 = vector.transfer_read %sv857[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<8x8xf64>
      %t878 = arith.mulf %d876, %g2877 : vector<8x8xf64>
      %a879 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2880 = arith.addf %t878, %a879 : vector<8x8xf64>
      %b881 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl882 = arith.addf %t2880, %b881 : vector<8x8xf64>
      vector.transfer_write %fl882, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc883 = arith.constant dense<0.0> : vector<8x8xf64>
      %a884 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b885 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r886 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a884, %b885, %acc883 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a887 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b888 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r889 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a887, %b888, %r886 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r889, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc890 = arith.constant dense<0.0> : vector<8x8xf64>
      %a891 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b892 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r893 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a891, %b892, %acc890 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a894 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b895 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r896 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a894, %b895, %r893 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r896, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv897 = memref.subview %y3[0, 5, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 40>>
    %sv898 = memref.subview %y3[0, 6, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 48>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa899 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r900 = vector.transfer_read %sv897[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 40>>, vector<8x8xf64>
      %m901 = arith.subf %r900, %fa899 : vector<8x8xf64>
      vector.transfer_write %m901, %sv897[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: 40>>
      %fa902 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r903 = vector.transfer_read %sv898[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 48>>, vector<8x8xf64>
      %p904 = arith.addf %r903, %fa902 : vector<8x8xf64>
      vector.transfer_write %p904, %sv898[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: 48>>
    }
    gpu.barrier
    %sv905 = memref.subview %v3591[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv906 = memref.subview %v3592[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv907 = memref.subview %g10all[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 384>>
    %sv908 = memref.subview %g11all[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 384>>
    %sv909 = memref.subview %g12all[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 384>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc910 = arith.constant dense<0.0> : vector<8x8xf64>
      %a911 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b912 = vector.transfer_read %sv905[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r913 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a911, %b912, %acc910 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a914 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b915 = vector.transfer_read %sv905[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r916 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a914, %b915, %r913 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g917 = vector.transfer_read %sv908[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<8x8xf64>
      %sc918 = arith.mulf %r916, %g917 : vector<8x8xf64>
      vector.transfer_write %sc918, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc919 = arith.constant dense<0.0> : vector<8x8xf64>
      %a920 = vector.transfer_read %sv905[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b921 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r922 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a920, %b921, %acc919 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a923 = vector.transfer_read %sv905[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b924 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r925 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a923, %b924, %r922 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g926 = vector.transfer_read %sv907[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<8x8xf64>
      %sc927 = arith.mulf %r925, %g926 : vector<8x8xf64>
      vector.transfer_write %sc927, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d928 = vector.transfer_read %sv906[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2929 = vector.transfer_read %sv909[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<8x8xf64>
      %t930 = arith.mulf %d928, %g2929 : vector<8x8xf64>
      %a931 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2932 = arith.addf %t930, %a931 : vector<8x8xf64>
      %b933 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl934 = arith.addf %t2932, %b933 : vector<8x8xf64>
      vector.transfer_write %fl934, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc935 = arith.constant dense<0.0> : vector<8x8xf64>
      %a936 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b937 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r938 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a936, %b937, %acc935 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a939 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b940 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r941 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a939, %b940, %r938 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r941, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc942 = arith.constant dense<0.0> : vector<8x8xf64>
      %a943 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b944 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r945 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a943, %b944, %acc942 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a946 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b947 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r948 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a946, %b947, %r945 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r948, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv949 = memref.subview %y3[0, 6, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 48>>
    %sv950 = memref.subview %y3[0, 7, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 56>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa951 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r952 = vector.transfer_read %sv949[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 48>>, vector<8x8xf64>
      %m953 = arith.subf %r952, %fa951 : vector<8x8xf64>
      vector.transfer_write %m953, %sv949[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: 48>>
      %fa954 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r955 = vector.transfer_read %sv950[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 56>>, vector<8x8xf64>
      %p956 = arith.addf %r955, %fa954 : vector<8x8xf64>
      vector.transfer_write %p956, %sv950[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: 56>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc957 = arith.constant dense<0.0> : vector<8x8xf64>
      %a958 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b959 = vector.transfer_read %u2[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r960 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a958, %b959, %acc957 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a961 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b962 = vector.transfer_read %u2[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r963 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a961, %b962, %r960 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r963, %interp_all[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc964 = arith.constant dense<0.0> : vector<8x8xf64>
      %a965 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b966 = vector.transfer_read %u2[%c0, %c8], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r967 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a965, %b966, %acc964 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a968 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b969 = vector.transfer_read %u2[%c4, %c8], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r970 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a968, %b969, %r967 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r970, %interp_all[%c0, %c8] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc971 = arith.constant dense<0.0> : vector<8x8xf64>
      %a972 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b973 = vector.transfer_read %u2[%c0, %c16], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r974 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a972, %b973, %acc971 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a975 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b976 = vector.transfer_read %u2[%c4, %c16], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r977 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a975, %b976, %r974 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r977, %interp_all[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc978 = arith.constant dense<0.0> : vector<8x8xf64>
      %a979 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b980 = vector.transfer_read %u2[%c0, %c24], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r981 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a979, %b980, %acc978 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a982 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b983 = vector.transfer_read %u2[%c4, %c24], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r984 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a982, %b983, %r981 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r984, %interp_all[%c0, %c24] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc985 = arith.constant dense<0.0> : vector<8x8xf64>
      %a986 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b987 = vector.transfer_read %u2[%c0, %c32], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r988 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a986, %b987, %acc985 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a989 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b990 = vector.transfer_read %u2[%c4, %c32], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r991 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a989, %b990, %r988 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r991, %interp_all[%c0, %c32] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc992 = arith.constant dense<0.0> : vector<8x8xf64>
      %a993 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b994 = vector.transfer_read %u2[%c0, %c40], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r995 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a993, %b994, %acc992 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a996 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b997 = vector.transfer_read %u2[%c4, %c40], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r998 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a996, %b997, %r995 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r998, %interp_all[%c0, %c40] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc999 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1000 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1001 = vector.transfer_read %u2[%c0, %c48], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r1002 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1000, %b1001, %acc999 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1003 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1004 = vector.transfer_read %u2[%c4, %c48], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r1005 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1003, %b1004, %r1002 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1005, %interp_all[%c0, %c48] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1006 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1007 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1008 = vector.transfer_read %u2[%c0, %c56], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r1009 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1007, %b1008, %acc1006 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1010 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1011 = vector.transfer_read %u2[%c4, %c56], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r1012 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1010, %b1011, %r1009 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1012, %interp_all[%c0, %c56] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1013 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1014 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1015 = vector.transfer_read %u2[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r1016 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1014, %b1015, %acc1013 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1017 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1018 = vector.transfer_read %u2[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r1019 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1017, %b1018, %r1016 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1019, %deriv_all[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1020 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1021 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1022 = vector.transfer_read %u2[%c0, %c8], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r1023 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1021, %b1022, %acc1020 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1024 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1025 = vector.transfer_read %u2[%c4, %c8], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r1026 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1024, %b1025, %r1023 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1026, %deriv_all[%c0, %c8] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1027 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1028 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1029 = vector.transfer_read %u2[%c0, %c16], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r1030 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1028, %b1029, %acc1027 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1031 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1032 = vector.transfer_read %u2[%c4, %c16], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r1033 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1031, %b1032, %r1030 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1033, %deriv_all[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1034 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1035 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1036 = vector.transfer_read %u2[%c0, %c24], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r1037 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1035, %b1036, %acc1034 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1038 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1039 = vector.transfer_read %u2[%c4, %c24], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r1040 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1038, %b1039, %r1037 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1040, %deriv_all[%c0, %c24] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1041 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1042 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1043 = vector.transfer_read %u2[%c0, %c32], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r1044 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1042, %b1043, %acc1041 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1045 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1046 = vector.transfer_read %u2[%c4, %c32], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r1047 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1045, %b1046, %r1044 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1047, %deriv_all[%c0, %c32] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1048 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1049 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1050 = vector.transfer_read %u2[%c0, %c40], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r1051 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1049, %b1050, %acc1048 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1052 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1053 = vector.transfer_read %u2[%c4, %c40], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r1054 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1052, %b1053, %r1051 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1054, %deriv_all[%c0, %c40] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1055 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1056 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1057 = vector.transfer_read %u2[%c0, %c48], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r1058 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1056, %b1057, %acc1055 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1059 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1060 = vector.transfer_read %u2[%c4, %c48], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r1061 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1059, %b1060, %r1058 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1061, %deriv_all[%c0, %c48] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1062 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1063 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1064 = vector.transfer_read %u2[%c0, %c56], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r1065 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1063, %b1064, %acc1062 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1066 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1067 = vector.transfer_read %u2[%c4, %c56], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<4x8xf64>
      %r1068 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1066, %b1067, %r1065 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1068, %deriv_all[%c0, %c56] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %v31069 = memref.reinterpret_cast %interp_all to offset: [0], sizes: [8, 8, 8], strides: [64, 8, 1] : memref<8x64xf64, #gpu.address_space<workgroup>> to memref<8x8x8xf64, #gpu.address_space<workgroup>>
    %v31070 = memref.reinterpret_cast %deriv_all to offset: [0], sizes: [8, 8, 8], strides: [64, 8, 1] : memref<8x64xf64, #gpu.address_space<workgroup>> to memref<8x8x8xf64, #gpu.address_space<workgroup>>
    %sv1071 = memref.subview %v31069[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>
    %sv1072 = memref.subview %v31070[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>
    %sv1073 = memref.subview %g20all[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 0>>
    %sv1074 = memref.subview %g21all[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 0>>
    %sv1075 = memref.subview %g22all[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 0>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1076 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1077 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1078 = vector.transfer_read %sv1071[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1079 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1077, %b1078, %acc1076 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1080 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1081 = vector.transfer_read %sv1071[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1082 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1080, %b1081, %r1079 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1083 = vector.transfer_read %sv1074[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>>, vector<8x8xf64>
      %sc1084 = arith.mulf %r1082, %g1083 : vector<8x8xf64>
      vector.transfer_write %sc1084, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1085 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1086 = vector.transfer_read %sv1071[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1087 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1088 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1086, %b1087, %acc1085 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1089 = vector.transfer_read %sv1071[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1090 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1091 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1089, %b1090, %r1088 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1092 = vector.transfer_read %sv1073[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>>, vector<8x8xf64>
      %sc1093 = arith.mulf %r1091, %g1092 : vector<8x8xf64>
      vector.transfer_write %sc1093, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d1094 = vector.transfer_read %sv1072[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g21095 = vector.transfer_read %sv1075[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>>, vector<8x8xf64>
      %t1096 = arith.mulf %d1094, %g21095 : vector<8x8xf64>
      %a1097 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t21098 = arith.addf %t1096, %a1097 : vector<8x8xf64>
      %b1099 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl1100 = arith.addf %t21098, %b1099 : vector<8x8xf64>
      vector.transfer_write %fl1100, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1101 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1102 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1103 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1104 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1102, %b1103, %acc1101 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1105 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1106 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1107 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1105, %b1106, %r1104 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1107, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1108 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1109 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1110 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1111 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1109, %b1110, %acc1108 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1112 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1113 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1114 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1112, %b1113, %r1111 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1114, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1115 = memref.subview %y3[0, 0, 0] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: 0>>
    %sv1116 = memref.subview %y3[0, 0, 1] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: 1>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1117 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1118 = vector.transfer_read %sv1115[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: 0>>, vector<8x8xf64>
      %m1119 = arith.subf %r1118, %fa1117 : vector<8x8xf64>
      vector.transfer_write %m1119, %sv1115[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: 0>>
      %fa1120 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1121 = vector.transfer_read %sv1116[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: 1>>, vector<8x8xf64>
      %p1122 = arith.addf %r1121, %fa1120 : vector<8x8xf64>
      vector.transfer_write %p1122, %sv1116[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: 1>>
    }
    gpu.barrier
    %sv1123 = memref.subview %v31069[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv1124 = memref.subview %v31070[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv1125 = memref.subview %g20all[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 64>>
    %sv1126 = memref.subview %g21all[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 64>>
    %sv1127 = memref.subview %g22all[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 64>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1128 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1129 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1130 = vector.transfer_read %sv1123[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1131 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1129, %b1130, %acc1128 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1132 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1133 = vector.transfer_read %sv1123[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1134 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1132, %b1133, %r1131 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1135 = vector.transfer_read %sv1126[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<8x8xf64>
      %sc1136 = arith.mulf %r1134, %g1135 : vector<8x8xf64>
      vector.transfer_write %sc1136, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1137 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1138 = vector.transfer_read %sv1123[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1139 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1140 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1138, %b1139, %acc1137 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1141 = vector.transfer_read %sv1123[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1142 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1143 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1141, %b1142, %r1140 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1144 = vector.transfer_read %sv1125[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<8x8xf64>
      %sc1145 = arith.mulf %r1143, %g1144 : vector<8x8xf64>
      vector.transfer_write %sc1145, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d1146 = vector.transfer_read %sv1124[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g21147 = vector.transfer_read %sv1127[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<8x8xf64>
      %t1148 = arith.mulf %d1146, %g21147 : vector<8x8xf64>
      %a1149 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t21150 = arith.addf %t1148, %a1149 : vector<8x8xf64>
      %b1151 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl1152 = arith.addf %t21150, %b1151 : vector<8x8xf64>
      vector.transfer_write %fl1152, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1153 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1154 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1155 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1156 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1154, %b1155, %acc1153 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1157 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1158 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1159 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1157, %b1158, %r1156 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1159, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1160 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1161 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1162 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1163 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1161, %b1162, %acc1160 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1164 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1165 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1166 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1164, %b1165, %r1163 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1166, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1167 = memref.subview %y3[0, 0, 1] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: 1>>
    %sv1168 = memref.subview %y3[0, 0, 2] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: 2>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1169 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1170 = vector.transfer_read %sv1167[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: 1>>, vector<8x8xf64>
      %m1171 = arith.subf %r1170, %fa1169 : vector<8x8xf64>
      vector.transfer_write %m1171, %sv1167[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: 1>>
      %fa1172 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1173 = vector.transfer_read %sv1168[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: 2>>, vector<8x8xf64>
      %p1174 = arith.addf %r1173, %fa1172 : vector<8x8xf64>
      vector.transfer_write %p1174, %sv1168[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: 2>>
    }
    gpu.barrier
    %sv1175 = memref.subview %v31069[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv1176 = memref.subview %v31070[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv1177 = memref.subview %g20all[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 128>>
    %sv1178 = memref.subview %g21all[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 128>>
    %sv1179 = memref.subview %g22all[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 128>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1180 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1181 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1182 = vector.transfer_read %sv1175[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1183 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1181, %b1182, %acc1180 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1184 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1185 = vector.transfer_read %sv1175[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1186 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1184, %b1185, %r1183 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1187 = vector.transfer_read %sv1178[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<8x8xf64>
      %sc1188 = arith.mulf %r1186, %g1187 : vector<8x8xf64>
      vector.transfer_write %sc1188, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1189 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1190 = vector.transfer_read %sv1175[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1191 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1192 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1190, %b1191, %acc1189 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1193 = vector.transfer_read %sv1175[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1194 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1195 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1193, %b1194, %r1192 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1196 = vector.transfer_read %sv1177[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<8x8xf64>
      %sc1197 = arith.mulf %r1195, %g1196 : vector<8x8xf64>
      vector.transfer_write %sc1197, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d1198 = vector.transfer_read %sv1176[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g21199 = vector.transfer_read %sv1179[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<8x8xf64>
      %t1200 = arith.mulf %d1198, %g21199 : vector<8x8xf64>
      %a1201 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t21202 = arith.addf %t1200, %a1201 : vector<8x8xf64>
      %b1203 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl1204 = arith.addf %t21202, %b1203 : vector<8x8xf64>
      vector.transfer_write %fl1204, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1205 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1206 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1207 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1208 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1206, %b1207, %acc1205 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1209 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1210 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1211 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1209, %b1210, %r1208 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1211, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1212 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1213 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1214 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1215 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1213, %b1214, %acc1212 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1216 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1217 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1218 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1216, %b1217, %r1215 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1218, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1219 = memref.subview %y3[0, 0, 2] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: 2>>
    %sv1220 = memref.subview %y3[0, 0, 3] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: 3>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1221 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1222 = vector.transfer_read %sv1219[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: 2>>, vector<8x8xf64>
      %m1223 = arith.subf %r1222, %fa1221 : vector<8x8xf64>
      vector.transfer_write %m1223, %sv1219[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: 2>>
      %fa1224 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1225 = vector.transfer_read %sv1220[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: 3>>, vector<8x8xf64>
      %p1226 = arith.addf %r1225, %fa1224 : vector<8x8xf64>
      vector.transfer_write %p1226, %sv1220[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: 3>>
    }
    gpu.barrier
    %sv1227 = memref.subview %v31069[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv1228 = memref.subview %v31070[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv1229 = memref.subview %g20all[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 192>>
    %sv1230 = memref.subview %g21all[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 192>>
    %sv1231 = memref.subview %g22all[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 192>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1232 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1233 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1234 = vector.transfer_read %sv1227[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1235 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1233, %b1234, %acc1232 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1236 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1237 = vector.transfer_read %sv1227[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1238 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1236, %b1237, %r1235 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1239 = vector.transfer_read %sv1230[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<8x8xf64>
      %sc1240 = arith.mulf %r1238, %g1239 : vector<8x8xf64>
      vector.transfer_write %sc1240, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1241 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1242 = vector.transfer_read %sv1227[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1243 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1244 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1242, %b1243, %acc1241 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1245 = vector.transfer_read %sv1227[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1246 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1247 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1245, %b1246, %r1244 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1248 = vector.transfer_read %sv1229[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<8x8xf64>
      %sc1249 = arith.mulf %r1247, %g1248 : vector<8x8xf64>
      vector.transfer_write %sc1249, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d1250 = vector.transfer_read %sv1228[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g21251 = vector.transfer_read %sv1231[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<8x8xf64>
      %t1252 = arith.mulf %d1250, %g21251 : vector<8x8xf64>
      %a1253 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t21254 = arith.addf %t1252, %a1253 : vector<8x8xf64>
      %b1255 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl1256 = arith.addf %t21254, %b1255 : vector<8x8xf64>
      vector.transfer_write %fl1256, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1257 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1258 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1259 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1260 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1258, %b1259, %acc1257 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1261 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1262 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1263 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1261, %b1262, %r1260 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1263, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1264 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1265 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1266 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1267 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1265, %b1266, %acc1264 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1268 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1269 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1270 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1268, %b1269, %r1267 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1270, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1271 = memref.subview %y3[0, 0, 3] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: 3>>
    %sv1272 = memref.subview %y3[0, 0, 4] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: 4>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1273 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1274 = vector.transfer_read %sv1271[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: 3>>, vector<8x8xf64>
      %m1275 = arith.subf %r1274, %fa1273 : vector<8x8xf64>
      vector.transfer_write %m1275, %sv1271[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: 3>>
      %fa1276 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1277 = vector.transfer_read %sv1272[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: 4>>, vector<8x8xf64>
      %p1278 = arith.addf %r1277, %fa1276 : vector<8x8xf64>
      vector.transfer_write %p1278, %sv1272[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: 4>>
    }
    gpu.barrier
    %sv1279 = memref.subview %v31069[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv1280 = memref.subview %v31070[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv1281 = memref.subview %g20all[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 256>>
    %sv1282 = memref.subview %g21all[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 256>>
    %sv1283 = memref.subview %g22all[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 256>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1284 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1285 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1286 = vector.transfer_read %sv1279[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1287 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1285, %b1286, %acc1284 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1288 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1289 = vector.transfer_read %sv1279[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1290 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1288, %b1289, %r1287 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1291 = vector.transfer_read %sv1282[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<8x8xf64>
      %sc1292 = arith.mulf %r1290, %g1291 : vector<8x8xf64>
      vector.transfer_write %sc1292, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1293 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1294 = vector.transfer_read %sv1279[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1295 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1296 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1294, %b1295, %acc1293 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1297 = vector.transfer_read %sv1279[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1298 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1299 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1297, %b1298, %r1296 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1300 = vector.transfer_read %sv1281[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<8x8xf64>
      %sc1301 = arith.mulf %r1299, %g1300 : vector<8x8xf64>
      vector.transfer_write %sc1301, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d1302 = vector.transfer_read %sv1280[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g21303 = vector.transfer_read %sv1283[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<8x8xf64>
      %t1304 = arith.mulf %d1302, %g21303 : vector<8x8xf64>
      %a1305 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t21306 = arith.addf %t1304, %a1305 : vector<8x8xf64>
      %b1307 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl1308 = arith.addf %t21306, %b1307 : vector<8x8xf64>
      vector.transfer_write %fl1308, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1309 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1310 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1311 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1312 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1310, %b1311, %acc1309 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1313 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1314 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1315 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1313, %b1314, %r1312 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1315, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1316 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1317 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1318 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1319 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1317, %b1318, %acc1316 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1320 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1321 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1322 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1320, %b1321, %r1319 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1322, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1323 = memref.subview %y3[0, 0, 4] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: 4>>
    %sv1324 = memref.subview %y3[0, 0, 5] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: 5>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1325 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1326 = vector.transfer_read %sv1323[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: 4>>, vector<8x8xf64>
      %m1327 = arith.subf %r1326, %fa1325 : vector<8x8xf64>
      vector.transfer_write %m1327, %sv1323[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: 4>>
      %fa1328 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1329 = vector.transfer_read %sv1324[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: 5>>, vector<8x8xf64>
      %p1330 = arith.addf %r1329, %fa1328 : vector<8x8xf64>
      vector.transfer_write %p1330, %sv1324[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: 5>>
    }
    gpu.barrier
    %sv1331 = memref.subview %v31069[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv1332 = memref.subview %v31070[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv1333 = memref.subview %g20all[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 320>>
    %sv1334 = memref.subview %g21all[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 320>>
    %sv1335 = memref.subview %g22all[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 320>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1336 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1337 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1338 = vector.transfer_read %sv1331[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1339 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1337, %b1338, %acc1336 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1340 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1341 = vector.transfer_read %sv1331[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1342 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1340, %b1341, %r1339 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1343 = vector.transfer_read %sv1334[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<8x8xf64>
      %sc1344 = arith.mulf %r1342, %g1343 : vector<8x8xf64>
      vector.transfer_write %sc1344, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1345 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1346 = vector.transfer_read %sv1331[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1347 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1348 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1346, %b1347, %acc1345 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1349 = vector.transfer_read %sv1331[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1350 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1351 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1349, %b1350, %r1348 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1352 = vector.transfer_read %sv1333[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<8x8xf64>
      %sc1353 = arith.mulf %r1351, %g1352 : vector<8x8xf64>
      vector.transfer_write %sc1353, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d1354 = vector.transfer_read %sv1332[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g21355 = vector.transfer_read %sv1335[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<8x8xf64>
      %t1356 = arith.mulf %d1354, %g21355 : vector<8x8xf64>
      %a1357 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t21358 = arith.addf %t1356, %a1357 : vector<8x8xf64>
      %b1359 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl1360 = arith.addf %t21358, %b1359 : vector<8x8xf64>
      vector.transfer_write %fl1360, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1361 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1362 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1363 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1364 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1362, %b1363, %acc1361 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1365 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1366 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1367 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1365, %b1366, %r1364 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1367, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1368 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1369 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1370 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1371 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1369, %b1370, %acc1368 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1372 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1373 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1374 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1372, %b1373, %r1371 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1374, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1375 = memref.subview %y3[0, 0, 5] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: 5>>
    %sv1376 = memref.subview %y3[0, 0, 6] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: 6>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1377 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1378 = vector.transfer_read %sv1375[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: 5>>, vector<8x8xf64>
      %m1379 = arith.subf %r1378, %fa1377 : vector<8x8xf64>
      vector.transfer_write %m1379, %sv1375[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: 5>>
      %fa1380 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1381 = vector.transfer_read %sv1376[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: 6>>, vector<8x8xf64>
      %p1382 = arith.addf %r1381, %fa1380 : vector<8x8xf64>
      vector.transfer_write %p1382, %sv1376[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: 6>>
    }
    gpu.barrier
    %sv1383 = memref.subview %v31069[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv1384 = memref.subview %v31070[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv1385 = memref.subview %g20all[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 384>>
    %sv1386 = memref.subview %g21all[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 384>>
    %sv1387 = memref.subview %g22all[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 384>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1388 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1389 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1390 = vector.transfer_read %sv1383[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1391 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1389, %b1390, %acc1388 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1392 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1393 = vector.transfer_read %sv1383[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1394 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1392, %b1393, %r1391 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1395 = vector.transfer_read %sv1386[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<8x8xf64>
      %sc1396 = arith.mulf %r1394, %g1395 : vector<8x8xf64>
      vector.transfer_write %sc1396, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1397 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1398 = vector.transfer_read %sv1383[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1399 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1400 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1398, %b1399, %acc1397 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1401 = vector.transfer_read %sv1383[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1402 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1403 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1401, %b1402, %r1400 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1404 = vector.transfer_read %sv1385[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<8x8xf64>
      %sc1405 = arith.mulf %r1403, %g1404 : vector<8x8xf64>
      vector.transfer_write %sc1405, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d1406 = vector.transfer_read %sv1384[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g21407 = vector.transfer_read %sv1387[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<8x8xf64>
      %t1408 = arith.mulf %d1406, %g21407 : vector<8x8xf64>
      %a1409 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t21410 = arith.addf %t1408, %a1409 : vector<8x8xf64>
      %b1411 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl1412 = arith.addf %t21410, %b1411 : vector<8x8xf64>
      vector.transfer_write %fl1412, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1413 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1414 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1415 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1416 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1414, %b1415, %acc1413 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1417 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1418 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1419 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1417, %b1418, %r1416 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1419, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1420 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1421 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1422 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1423 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1421, %b1422, %acc1420 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1424 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1425 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1426 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1424, %b1425, %r1423 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1426, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1427 = memref.subview %y3[0, 0, 6] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: 6>>
    %sv1428 = memref.subview %y3[0, 0, 7] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: 7>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1429 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1430 = vector.transfer_read %sv1427[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: 6>>, vector<8x8xf64>
      %m1431 = arith.subf %r1430, %fa1429 : vector<8x8xf64>
      vector.transfer_write %m1431, %sv1427[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: 6>>
      %fa1432 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1433 = vector.transfer_read %sv1428[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: 7>>, vector<8x8xf64>
      %p1434 = arith.addf %r1433, %fa1432 : vector<8x8xf64>
      vector.transfer_write %p1434, %sv1428[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: 7>>
    }
    gpu.barrier
    gpu.return
  }
}
