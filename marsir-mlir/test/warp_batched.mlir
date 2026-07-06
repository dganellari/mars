#a = affine_map<(m, n, k) -> (m, k)>
#b = affine_map<(m, n, k) -> (k, n)>
#bt = affine_map<(m, n, k) -> (n, k)>
#c = affine_map<(m, n, k) -> (m, n)>
gpu.module @mir_kernels [#nvvm.target<chip = "sm_90", O = 3>] {
  gpu.func @full_batched(%U: memref<?x8x8x8xf64>, %Btil: memref<8x8xf64>, %Dtil: memref<8x8xf64>, %Dm: memref<8x8xf64>, %W: memref<8x8xf64>, %g00all: memref<?x8x8x8xf64>, %g01all: memref<?x8x8x8xf64>, %g02all: memref<?x8x8x8xf64>, %g10all: memref<?x8x8x8xf64>, %g11all: memref<?x8x8x8xf64>, %g12all: memref<?x8x8x8xf64>, %g20all: memref<?x8x8x8xf64>, %g21all: memref<?x8x8x8xf64>, %g22all: memref<?x8x8x8xf64>, %Y: memref<?x8x8x8xf64>) workgroup(%interp_all: memref<8x64xf64, #gpu.address_space<workgroup>>, %deriv_all: memref<8x64xf64, #gpu.address_space<workgroup>>, %flux: memref<8x8xf64, #gpu.address_space<workgroup>>, %tmp: memref<8x8xf64, #gpu.address_space<workgroup>>, %intf: memref<8x8xf64, #gpu.address_space<workgroup>>) kernel {
    %z = arith.constant 0.0 : f64
    %lane = gpu.thread_id x
    %e = gpu.block_id x
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
    %el1 = memref.subview %U[%e, 0, 0, 0] [1, 8, 8, 8] [1, 1, 1, 1] : memref<?x8x8x8xf64> to memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>>
    %el2 = memref.subview %Y[%e, 0, 0, 0] [1, 8, 8, 8] [1, 1, 1, 1] : memref<?x8x8x8xf64> to memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>>
    %el3 = memref.subview %g00all[%e, 0, 0, 0] [1, 8, 8, 8] [1, 1, 1, 1] : memref<?x8x8x8xf64> to memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>>
    %el4 = memref.subview %g01all[%e, 0, 0, 0] [1, 8, 8, 8] [1, 1, 1, 1] : memref<?x8x8x8xf64> to memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>>
    %el5 = memref.subview %g02all[%e, 0, 0, 0] [1, 8, 8, 8] [1, 1, 1, 1] : memref<?x8x8x8xf64> to memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>>
    %el6 = memref.subview %g10all[%e, 0, 0, 0] [1, 8, 8, 8] [1, 1, 1, 1] : memref<?x8x8x8xf64> to memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>>
    %el7 = memref.subview %g11all[%e, 0, 0, 0] [1, 8, 8, 8] [1, 1, 1, 1] : memref<?x8x8x8xf64> to memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>>
    %el8 = memref.subview %g12all[%e, 0, 0, 0] [1, 8, 8, 8] [1, 1, 1, 1] : memref<?x8x8x8xf64> to memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>>
    %el9 = memref.subview %g20all[%e, 0, 0, 0] [1, 8, 8, 8] [1, 1, 1, 1] : memref<?x8x8x8xf64> to memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>>
    %el10 = memref.subview %g21all[%e, 0, 0, 0] [1, 8, 8, 8] [1, 1, 1, 1] : memref<?x8x8x8xf64> to memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>>
    %el11 = memref.subview %g22all[%e, 0, 0, 0] [1, 8, 8, 8] [1, 1, 1, 1] : memref<?x8x8x8xf64> to memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>>
    %sv12 = memref.subview %el1[0, 0, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc13 = arith.constant dense<0.0> : vector<8x8xf64>
      %a14 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b15 = vector.transfer_read %sv12[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
      %r16 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a14, %b15, %acc13 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a17 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b18 = vector.transfer_read %sv12[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
      %r19 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a17, %b18, %r16 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r19, %interp_all[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv20 = memref.subview %el1[0, 1, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc21 = arith.constant dense<0.0> : vector<8x8xf64>
      %a22 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b23 = vector.transfer_read %sv20[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
      %r24 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a22, %b23, %acc21 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a25 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b26 = vector.transfer_read %sv20[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
      %r27 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a25, %b26, %r24 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r27, %interp_all[%c0, %c8] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv28 = memref.subview %el1[0, 2, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc29 = arith.constant dense<0.0> : vector<8x8xf64>
      %a30 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b31 = vector.transfer_read %sv28[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
      %r32 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a30, %b31, %acc29 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a33 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b34 = vector.transfer_read %sv28[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
      %r35 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a33, %b34, %r32 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r35, %interp_all[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv36 = memref.subview %el1[0, 3, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc37 = arith.constant dense<0.0> : vector<8x8xf64>
      %a38 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b39 = vector.transfer_read %sv36[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
      %r40 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a38, %b39, %acc37 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a41 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b42 = vector.transfer_read %sv36[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
      %r43 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a41, %b42, %r40 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r43, %interp_all[%c0, %c24] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv44 = memref.subview %el1[0, 4, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc45 = arith.constant dense<0.0> : vector<8x8xf64>
      %a46 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b47 = vector.transfer_read %sv44[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
      %r48 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a46, %b47, %acc45 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a49 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b50 = vector.transfer_read %sv44[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
      %r51 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a49, %b50, %r48 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r51, %interp_all[%c0, %c32] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv52 = memref.subview %el1[0, 5, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc53 = arith.constant dense<0.0> : vector<8x8xf64>
      %a54 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b55 = vector.transfer_read %sv52[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
      %r56 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a54, %b55, %acc53 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a57 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b58 = vector.transfer_read %sv52[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
      %r59 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a57, %b58, %r56 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r59, %interp_all[%c0, %c40] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv60 = memref.subview %el1[0, 6, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc61 = arith.constant dense<0.0> : vector<8x8xf64>
      %a62 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b63 = vector.transfer_read %sv60[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
      %r64 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a62, %b63, %acc61 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a65 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b66 = vector.transfer_read %sv60[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
      %r67 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a65, %b66, %r64 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r67, %interp_all[%c0, %c48] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv68 = memref.subview %el1[0, 7, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc69 = arith.constant dense<0.0> : vector<8x8xf64>
      %a70 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b71 = vector.transfer_read %sv68[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
      %r72 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a70, %b71, %acc69 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a73 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b74 = vector.transfer_read %sv68[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
      %r75 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a73, %b74, %r72 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r75, %interp_all[%c0, %c56] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv76 = memref.subview %el1[0, 0, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc77 = arith.constant dense<0.0> : vector<8x8xf64>
      %a78 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b79 = vector.transfer_read %sv76[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
      %r80 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a78, %b79, %acc77 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a81 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b82 = vector.transfer_read %sv76[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
      %r83 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a81, %b82, %r80 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r83, %deriv_all[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv84 = memref.subview %el1[0, 1, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc85 = arith.constant dense<0.0> : vector<8x8xf64>
      %a86 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b87 = vector.transfer_read %sv84[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
      %r88 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a86, %b87, %acc85 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a89 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b90 = vector.transfer_read %sv84[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
      %r91 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a89, %b90, %r88 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r91, %deriv_all[%c0, %c8] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv92 = memref.subview %el1[0, 2, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc93 = arith.constant dense<0.0> : vector<8x8xf64>
      %a94 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b95 = vector.transfer_read %sv92[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
      %r96 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a94, %b95, %acc93 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a97 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b98 = vector.transfer_read %sv92[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
      %r99 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a97, %b98, %r96 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r99, %deriv_all[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv100 = memref.subview %el1[0, 3, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc101 = arith.constant dense<0.0> : vector<8x8xf64>
      %a102 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b103 = vector.transfer_read %sv100[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
      %r104 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a102, %b103, %acc101 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a105 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b106 = vector.transfer_read %sv100[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
      %r107 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a105, %b106, %r104 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r107, %deriv_all[%c0, %c24] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv108 = memref.subview %el1[0, 4, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc109 = arith.constant dense<0.0> : vector<8x8xf64>
      %a110 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b111 = vector.transfer_read %sv108[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
      %r112 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a110, %b111, %acc109 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a113 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b114 = vector.transfer_read %sv108[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
      %r115 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a113, %b114, %r112 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r115, %deriv_all[%c0, %c32] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv116 = memref.subview %el1[0, 5, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc117 = arith.constant dense<0.0> : vector<8x8xf64>
      %a118 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b119 = vector.transfer_read %sv116[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
      %r120 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a118, %b119, %acc117 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a121 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b122 = vector.transfer_read %sv116[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
      %r123 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a121, %b122, %r120 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r123, %deriv_all[%c0, %c40] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv124 = memref.subview %el1[0, 6, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc125 = arith.constant dense<0.0> : vector<8x8xf64>
      %a126 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b127 = vector.transfer_read %sv124[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
      %r128 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a126, %b127, %acc125 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a129 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b130 = vector.transfer_read %sv124[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
      %r131 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a129, %b130, %r128 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r131, %deriv_all[%c0, %c48] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv132 = memref.subview %el1[0, 7, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc133 = arith.constant dense<0.0> : vector<8x8xf64>
      %a134 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b135 = vector.transfer_read %sv132[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
      %r136 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a134, %b135, %acc133 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a137 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b138 = vector.transfer_read %sv132[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
      %r139 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a137, %b138, %r136 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r139, %deriv_all[%c0, %c56] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %v3140 = memref.reinterpret_cast %interp_all to offset: [0], sizes: [8, 8, 8], strides: [64, 8, 1] : memref<8x64xf64, #gpu.address_space<workgroup>> to memref<8x8x8xf64, #gpu.address_space<workgroup>>
    %v3141 = memref.reinterpret_cast %deriv_all to offset: [0], sizes: [8, 8, 8], strides: [64, 8, 1] : memref<8x64xf64, #gpu.address_space<workgroup>> to memref<8x8x8xf64, #gpu.address_space<workgroup>>
    %sv142 = memref.subview %v3140[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>
    %sv143 = memref.subview %v3141[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>
    %sv144 = memref.subview %el3[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv145 = memref.subview %el4[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv146 = memref.subview %el5[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc147 = arith.constant dense<0.0> : vector<8x8xf64>
      %a148 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b149 = vector.transfer_read %sv142[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r150 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a148, %b149, %acc147 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a151 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b152 = vector.transfer_read %sv142[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r153 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a151, %b152, %r150 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g154 = vector.transfer_read %sv145[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt155 = arith.mulf %r153, %g154 : vector<8x8xf64>
      %acc156 = arith.constant dense<0.0> : vector<8x8xf64>
      %a157 = vector.transfer_read %sv142[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b158 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r159 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a157, %b158, %acc156 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a160 = vector.transfer_read %sv142[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b161 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r162 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a160, %b161, %r159 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g163 = vector.transfer_read %sv144[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt164 = arith.mulf %r162, %g163 : vector<8x8xf64>
      %dv165 = vector.transfer_read %sv143[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g166 = vector.transfer_read %sv146[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t167 = arith.mulf %dv165, %g166 : vector<8x8xf64>
      %t2168 = arith.addf %t167, %dt164 : vector<8x8xf64>
      %fl169 = arith.addf %t2168, %dt155 : vector<8x8xf64>
      vector.transfer_write %fl169, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc170 = arith.constant dense<0.0> : vector<8x8xf64>
      %a171 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b172 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r173 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a171, %b172, %acc170 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a174 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b175 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r176 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a174, %b175, %r173 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r176, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc177 = arith.constant dense<0.0> : vector<8x8xf64>
      %a178 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b179 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r180 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a178, %b179, %acc177 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a181 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b182 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r183 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a181, %b182, %r180 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r183, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv184 = memref.subview %el2[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv185 = memref.subview %el2[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa186 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r187 = vector.transfer_read %sv184[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %m188 = arith.subf %r187, %fa186 : vector<8x8xf64>
      vector.transfer_write %m188, %sv184[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      %r189 = vector.transfer_read %sv185[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %p190 = arith.addf %r189, %fa186 : vector<8x8xf64>
      vector.transfer_write %p190, %sv185[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    }
    gpu.barrier
    %sv191 = memref.subview %v3140[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv192 = memref.subview %v3141[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv193 = memref.subview %el3[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv194 = memref.subview %el4[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv195 = memref.subview %el5[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc196 = arith.constant dense<0.0> : vector<8x8xf64>
      %a197 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b198 = vector.transfer_read %sv191[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r199 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a197, %b198, %acc196 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a200 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b201 = vector.transfer_read %sv191[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r202 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a200, %b201, %r199 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g203 = vector.transfer_read %sv194[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt204 = arith.mulf %r202, %g203 : vector<8x8xf64>
      %acc205 = arith.constant dense<0.0> : vector<8x8xf64>
      %a206 = vector.transfer_read %sv191[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b207 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r208 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a206, %b207, %acc205 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a209 = vector.transfer_read %sv191[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b210 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r211 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a209, %b210, %r208 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g212 = vector.transfer_read %sv193[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt213 = arith.mulf %r211, %g212 : vector<8x8xf64>
      %dv214 = vector.transfer_read %sv192[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g215 = vector.transfer_read %sv195[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t216 = arith.mulf %dv214, %g215 : vector<8x8xf64>
      %t2217 = arith.addf %t216, %dt213 : vector<8x8xf64>
      %fl218 = arith.addf %t2217, %dt204 : vector<8x8xf64>
      vector.transfer_write %fl218, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc219 = arith.constant dense<0.0> : vector<8x8xf64>
      %a220 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b221 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r222 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a220, %b221, %acc219 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a223 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b224 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r225 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a223, %b224, %r222 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r225, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc226 = arith.constant dense<0.0> : vector<8x8xf64>
      %a227 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b228 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r229 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a227, %b228, %acc226 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a230 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b231 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r232 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a230, %b231, %r229 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r232, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv233 = memref.subview %el2[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv234 = memref.subview %el2[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa235 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r236 = vector.transfer_read %sv233[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %m237 = arith.subf %r236, %fa235 : vector<8x8xf64>
      vector.transfer_write %m237, %sv233[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      %r238 = vector.transfer_read %sv234[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %p239 = arith.addf %r238, %fa235 : vector<8x8xf64>
      vector.transfer_write %p239, %sv234[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    }
    gpu.barrier
    %sv240 = memref.subview %v3140[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv241 = memref.subview %v3141[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv242 = memref.subview %el3[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv243 = memref.subview %el4[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv244 = memref.subview %el5[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc245 = arith.constant dense<0.0> : vector<8x8xf64>
      %a246 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b247 = vector.transfer_read %sv240[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r248 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a246, %b247, %acc245 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a249 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b250 = vector.transfer_read %sv240[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r251 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a249, %b250, %r248 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g252 = vector.transfer_read %sv243[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt253 = arith.mulf %r251, %g252 : vector<8x8xf64>
      %acc254 = arith.constant dense<0.0> : vector<8x8xf64>
      %a255 = vector.transfer_read %sv240[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b256 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r257 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a255, %b256, %acc254 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a258 = vector.transfer_read %sv240[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b259 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r260 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a258, %b259, %r257 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g261 = vector.transfer_read %sv242[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt262 = arith.mulf %r260, %g261 : vector<8x8xf64>
      %dv263 = vector.transfer_read %sv241[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g264 = vector.transfer_read %sv244[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t265 = arith.mulf %dv263, %g264 : vector<8x8xf64>
      %t2266 = arith.addf %t265, %dt262 : vector<8x8xf64>
      %fl267 = arith.addf %t2266, %dt253 : vector<8x8xf64>
      vector.transfer_write %fl267, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc268 = arith.constant dense<0.0> : vector<8x8xf64>
      %a269 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b270 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r271 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a269, %b270, %acc268 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a272 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b273 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r274 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a272, %b273, %r271 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r274, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc275 = arith.constant dense<0.0> : vector<8x8xf64>
      %a276 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b277 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r278 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a276, %b277, %acc275 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a279 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b280 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r281 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a279, %b280, %r278 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r281, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv282 = memref.subview %el2[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv283 = memref.subview %el2[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa284 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r285 = vector.transfer_read %sv282[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %m286 = arith.subf %r285, %fa284 : vector<8x8xf64>
      vector.transfer_write %m286, %sv282[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      %r287 = vector.transfer_read %sv283[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %p288 = arith.addf %r287, %fa284 : vector<8x8xf64>
      vector.transfer_write %p288, %sv283[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    }
    gpu.barrier
    %sv289 = memref.subview %v3140[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv290 = memref.subview %v3141[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv291 = memref.subview %el3[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv292 = memref.subview %el4[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv293 = memref.subview %el5[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc294 = arith.constant dense<0.0> : vector<8x8xf64>
      %a295 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b296 = vector.transfer_read %sv289[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r297 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a295, %b296, %acc294 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a298 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b299 = vector.transfer_read %sv289[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r300 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a298, %b299, %r297 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g301 = vector.transfer_read %sv292[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt302 = arith.mulf %r300, %g301 : vector<8x8xf64>
      %acc303 = arith.constant dense<0.0> : vector<8x8xf64>
      %a304 = vector.transfer_read %sv289[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b305 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r306 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a304, %b305, %acc303 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a307 = vector.transfer_read %sv289[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b308 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r309 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a307, %b308, %r306 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g310 = vector.transfer_read %sv291[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt311 = arith.mulf %r309, %g310 : vector<8x8xf64>
      %dv312 = vector.transfer_read %sv290[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g313 = vector.transfer_read %sv293[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t314 = arith.mulf %dv312, %g313 : vector<8x8xf64>
      %t2315 = arith.addf %t314, %dt311 : vector<8x8xf64>
      %fl316 = arith.addf %t2315, %dt302 : vector<8x8xf64>
      vector.transfer_write %fl316, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc317 = arith.constant dense<0.0> : vector<8x8xf64>
      %a318 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b319 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r320 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a318, %b319, %acc317 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a321 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b322 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r323 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a321, %b322, %r320 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r323, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc324 = arith.constant dense<0.0> : vector<8x8xf64>
      %a325 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b326 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r327 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a325, %b326, %acc324 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a328 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b329 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r330 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a328, %b329, %r327 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r330, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv331 = memref.subview %el2[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv332 = memref.subview %el2[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa333 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r334 = vector.transfer_read %sv331[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %m335 = arith.subf %r334, %fa333 : vector<8x8xf64>
      vector.transfer_write %m335, %sv331[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      %r336 = vector.transfer_read %sv332[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %p337 = arith.addf %r336, %fa333 : vector<8x8xf64>
      vector.transfer_write %p337, %sv332[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    }
    gpu.barrier
    %sv338 = memref.subview %v3140[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv339 = memref.subview %v3141[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv340 = memref.subview %el3[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv341 = memref.subview %el4[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv342 = memref.subview %el5[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc343 = arith.constant dense<0.0> : vector<8x8xf64>
      %a344 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b345 = vector.transfer_read %sv338[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r346 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a344, %b345, %acc343 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a347 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b348 = vector.transfer_read %sv338[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r349 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a347, %b348, %r346 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g350 = vector.transfer_read %sv341[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt351 = arith.mulf %r349, %g350 : vector<8x8xf64>
      %acc352 = arith.constant dense<0.0> : vector<8x8xf64>
      %a353 = vector.transfer_read %sv338[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b354 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r355 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a353, %b354, %acc352 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a356 = vector.transfer_read %sv338[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b357 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r358 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a356, %b357, %r355 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g359 = vector.transfer_read %sv340[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt360 = arith.mulf %r358, %g359 : vector<8x8xf64>
      %dv361 = vector.transfer_read %sv339[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g362 = vector.transfer_read %sv342[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t363 = arith.mulf %dv361, %g362 : vector<8x8xf64>
      %t2364 = arith.addf %t363, %dt360 : vector<8x8xf64>
      %fl365 = arith.addf %t2364, %dt351 : vector<8x8xf64>
      vector.transfer_write %fl365, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc366 = arith.constant dense<0.0> : vector<8x8xf64>
      %a367 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b368 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r369 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a367, %b368, %acc366 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a370 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b371 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r372 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a370, %b371, %r369 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r372, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc373 = arith.constant dense<0.0> : vector<8x8xf64>
      %a374 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b375 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r376 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a374, %b375, %acc373 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a377 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b378 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r379 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a377, %b378, %r376 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r379, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv380 = memref.subview %el2[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv381 = memref.subview %el2[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa382 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r383 = vector.transfer_read %sv380[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %m384 = arith.subf %r383, %fa382 : vector<8x8xf64>
      vector.transfer_write %m384, %sv380[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      %r385 = vector.transfer_read %sv381[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %p386 = arith.addf %r385, %fa382 : vector<8x8xf64>
      vector.transfer_write %p386, %sv381[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    }
    gpu.barrier
    %sv387 = memref.subview %v3140[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv388 = memref.subview %v3141[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv389 = memref.subview %el3[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv390 = memref.subview %el4[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv391 = memref.subview %el5[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc392 = arith.constant dense<0.0> : vector<8x8xf64>
      %a393 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b394 = vector.transfer_read %sv387[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r395 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a393, %b394, %acc392 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a396 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b397 = vector.transfer_read %sv387[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r398 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a396, %b397, %r395 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g399 = vector.transfer_read %sv390[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt400 = arith.mulf %r398, %g399 : vector<8x8xf64>
      %acc401 = arith.constant dense<0.0> : vector<8x8xf64>
      %a402 = vector.transfer_read %sv387[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b403 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r404 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a402, %b403, %acc401 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a405 = vector.transfer_read %sv387[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b406 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r407 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a405, %b406, %r404 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g408 = vector.transfer_read %sv389[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt409 = arith.mulf %r407, %g408 : vector<8x8xf64>
      %dv410 = vector.transfer_read %sv388[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g411 = vector.transfer_read %sv391[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t412 = arith.mulf %dv410, %g411 : vector<8x8xf64>
      %t2413 = arith.addf %t412, %dt409 : vector<8x8xf64>
      %fl414 = arith.addf %t2413, %dt400 : vector<8x8xf64>
      vector.transfer_write %fl414, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc415 = arith.constant dense<0.0> : vector<8x8xf64>
      %a416 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b417 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r418 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a416, %b417, %acc415 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a419 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b420 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r421 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a419, %b420, %r418 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r421, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc422 = arith.constant dense<0.0> : vector<8x8xf64>
      %a423 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b424 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r425 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a423, %b424, %acc422 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a426 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b427 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r428 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a426, %b427, %r425 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r428, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv429 = memref.subview %el2[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv430 = memref.subview %el2[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa431 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r432 = vector.transfer_read %sv429[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %m433 = arith.subf %r432, %fa431 : vector<8x8xf64>
      vector.transfer_write %m433, %sv429[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      %r434 = vector.transfer_read %sv430[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %p435 = arith.addf %r434, %fa431 : vector<8x8xf64>
      vector.transfer_write %p435, %sv430[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    }
    gpu.barrier
    %sv436 = memref.subview %v3140[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv437 = memref.subview %v3141[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv438 = memref.subview %el3[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv439 = memref.subview %el4[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv440 = memref.subview %el5[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc441 = arith.constant dense<0.0> : vector<8x8xf64>
      %a442 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b443 = vector.transfer_read %sv436[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r444 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a442, %b443, %acc441 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a445 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b446 = vector.transfer_read %sv436[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r447 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a445, %b446, %r444 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g448 = vector.transfer_read %sv439[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt449 = arith.mulf %r447, %g448 : vector<8x8xf64>
      %acc450 = arith.constant dense<0.0> : vector<8x8xf64>
      %a451 = vector.transfer_read %sv436[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b452 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r453 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a451, %b452, %acc450 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a454 = vector.transfer_read %sv436[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b455 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r456 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a454, %b455, %r453 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g457 = vector.transfer_read %sv438[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt458 = arith.mulf %r456, %g457 : vector<8x8xf64>
      %dv459 = vector.transfer_read %sv437[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g460 = vector.transfer_read %sv440[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t461 = arith.mulf %dv459, %g460 : vector<8x8xf64>
      %t2462 = arith.addf %t461, %dt458 : vector<8x8xf64>
      %fl463 = arith.addf %t2462, %dt449 : vector<8x8xf64>
      vector.transfer_write %fl463, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc464 = arith.constant dense<0.0> : vector<8x8xf64>
      %a465 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b466 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r467 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a465, %b466, %acc464 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a468 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b469 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r470 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a468, %b469, %r467 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r470, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc471 = arith.constant dense<0.0> : vector<8x8xf64>
      %a472 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b473 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r474 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a472, %b473, %acc471 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a475 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b476 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r477 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a475, %b476, %r474 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r477, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv478 = memref.subview %el2[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv479 = memref.subview %el2[7, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa480 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r481 = vector.transfer_read %sv478[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %m482 = arith.subf %r481, %fa480 : vector<8x8xf64>
      vector.transfer_write %m482, %sv478[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      %r483 = vector.transfer_read %sv479[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %p484 = arith.addf %r483, %fa480 : vector<8x8xf64>
      vector.transfer_write %p484, %sv479[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    }
    gpu.barrier
    %sv485 = memref.subview %el1[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc486 = arith.constant dense<0.0> : vector<8x8xf64>
      %a487 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b488 = vector.transfer_read %sv485[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r489 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a487, %b488, %acc486 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a490 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b491 = vector.transfer_read %sv485[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r492 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a490, %b491, %r489 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r492, %interp_all[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv493 = memref.subview %el1[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc494 = arith.constant dense<0.0> : vector<8x8xf64>
      %a495 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b496 = vector.transfer_read %sv493[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r497 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a495, %b496, %acc494 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a498 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b499 = vector.transfer_read %sv493[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r500 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a498, %b499, %r497 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r500, %interp_all[%c0, %c8] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv501 = memref.subview %el1[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc502 = arith.constant dense<0.0> : vector<8x8xf64>
      %a503 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b504 = vector.transfer_read %sv501[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r505 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a503, %b504, %acc502 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a506 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b507 = vector.transfer_read %sv501[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r508 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a506, %b507, %r505 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r508, %interp_all[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv509 = memref.subview %el1[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc510 = arith.constant dense<0.0> : vector<8x8xf64>
      %a511 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b512 = vector.transfer_read %sv509[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r513 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a511, %b512, %acc510 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a514 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b515 = vector.transfer_read %sv509[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r516 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a514, %b515, %r513 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r516, %interp_all[%c0, %c24] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv517 = memref.subview %el1[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc518 = arith.constant dense<0.0> : vector<8x8xf64>
      %a519 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b520 = vector.transfer_read %sv517[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r521 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a519, %b520, %acc518 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a522 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b523 = vector.transfer_read %sv517[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r524 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a522, %b523, %r521 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r524, %interp_all[%c0, %c32] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv525 = memref.subview %el1[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc526 = arith.constant dense<0.0> : vector<8x8xf64>
      %a527 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b528 = vector.transfer_read %sv525[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r529 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a527, %b528, %acc526 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a530 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b531 = vector.transfer_read %sv525[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r532 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a530, %b531, %r529 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r532, %interp_all[%c0, %c40] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv533 = memref.subview %el1[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc534 = arith.constant dense<0.0> : vector<8x8xf64>
      %a535 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b536 = vector.transfer_read %sv533[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r537 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a535, %b536, %acc534 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a538 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b539 = vector.transfer_read %sv533[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r540 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a538, %b539, %r537 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r540, %interp_all[%c0, %c48] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv541 = memref.subview %el1[7, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc542 = arith.constant dense<0.0> : vector<8x8xf64>
      %a543 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b544 = vector.transfer_read %sv541[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r545 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a543, %b544, %acc542 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a546 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b547 = vector.transfer_read %sv541[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r548 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a546, %b547, %r545 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r548, %interp_all[%c0, %c56] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv549 = memref.subview %el1[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc550 = arith.constant dense<0.0> : vector<8x8xf64>
      %a551 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b552 = vector.transfer_read %sv549[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r553 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a551, %b552, %acc550 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a554 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b555 = vector.transfer_read %sv549[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r556 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a554, %b555, %r553 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r556, %deriv_all[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv557 = memref.subview %el1[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc558 = arith.constant dense<0.0> : vector<8x8xf64>
      %a559 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b560 = vector.transfer_read %sv557[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r561 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a559, %b560, %acc558 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a562 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b563 = vector.transfer_read %sv557[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r564 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a562, %b563, %r561 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r564, %deriv_all[%c0, %c8] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv565 = memref.subview %el1[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc566 = arith.constant dense<0.0> : vector<8x8xf64>
      %a567 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b568 = vector.transfer_read %sv565[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r569 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a567, %b568, %acc566 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a570 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b571 = vector.transfer_read %sv565[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r572 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a570, %b571, %r569 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r572, %deriv_all[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv573 = memref.subview %el1[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc574 = arith.constant dense<0.0> : vector<8x8xf64>
      %a575 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b576 = vector.transfer_read %sv573[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r577 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a575, %b576, %acc574 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a578 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b579 = vector.transfer_read %sv573[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r580 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a578, %b579, %r577 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r580, %deriv_all[%c0, %c24] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv581 = memref.subview %el1[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc582 = arith.constant dense<0.0> : vector<8x8xf64>
      %a583 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b584 = vector.transfer_read %sv581[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r585 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a583, %b584, %acc582 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a586 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b587 = vector.transfer_read %sv581[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r588 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a586, %b587, %r585 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r588, %deriv_all[%c0, %c32] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv589 = memref.subview %el1[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc590 = arith.constant dense<0.0> : vector<8x8xf64>
      %a591 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b592 = vector.transfer_read %sv589[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r593 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a591, %b592, %acc590 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a594 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b595 = vector.transfer_read %sv589[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r596 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a594, %b595, %r593 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r596, %deriv_all[%c0, %c40] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv597 = memref.subview %el1[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc598 = arith.constant dense<0.0> : vector<8x8xf64>
      %a599 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b600 = vector.transfer_read %sv597[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r601 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a599, %b600, %acc598 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a602 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b603 = vector.transfer_read %sv597[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r604 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a602, %b603, %r601 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r604, %deriv_all[%c0, %c48] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv605 = memref.subview %el1[7, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc606 = arith.constant dense<0.0> : vector<8x8xf64>
      %a607 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b608 = vector.transfer_read %sv605[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r609 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a607, %b608, %acc606 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a610 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b611 = vector.transfer_read %sv605[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r612 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a610, %b611, %r609 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r612, %deriv_all[%c0, %c56] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %v3613 = memref.reinterpret_cast %interp_all to offset: [0], sizes: [8, 8, 8], strides: [64, 8, 1] : memref<8x64xf64, #gpu.address_space<workgroup>> to memref<8x8x8xf64, #gpu.address_space<workgroup>>
    %v3614 = memref.reinterpret_cast %deriv_all to offset: [0], sizes: [8, 8, 8], strides: [64, 8, 1] : memref<8x64xf64, #gpu.address_space<workgroup>> to memref<8x8x8xf64, #gpu.address_space<workgroup>>
    %sv615 = memref.subview %v3613[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>
    %sv616 = memref.subview %v3614[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>
    %sv617 = memref.subview %el6[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv618 = memref.subview %el7[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv619 = memref.subview %el8[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc620 = arith.constant dense<0.0> : vector<8x8xf64>
      %a621 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b622 = vector.transfer_read %sv615[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r623 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a621, %b622, %acc620 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a624 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b625 = vector.transfer_read %sv615[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r626 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a624, %b625, %r623 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g627 = vector.transfer_read %sv618[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt628 = arith.mulf %r626, %g627 : vector<8x8xf64>
      %acc629 = arith.constant dense<0.0> : vector<8x8xf64>
      %a630 = vector.transfer_read %sv615[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b631 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r632 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a630, %b631, %acc629 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a633 = vector.transfer_read %sv615[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b634 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r635 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a633, %b634, %r632 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g636 = vector.transfer_read %sv617[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt637 = arith.mulf %r635, %g636 : vector<8x8xf64>
      %dv638 = vector.transfer_read %sv616[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g639 = vector.transfer_read %sv619[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t640 = arith.mulf %dv638, %g639 : vector<8x8xf64>
      %t2641 = arith.addf %t640, %dt637 : vector<8x8xf64>
      %fl642 = arith.addf %t2641, %dt628 : vector<8x8xf64>
      vector.transfer_write %fl642, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc643 = arith.constant dense<0.0> : vector<8x8xf64>
      %a644 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b645 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r646 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a644, %b645, %acc643 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a647 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b648 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r649 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a647, %b648, %r646 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r649, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc650 = arith.constant dense<0.0> : vector<8x8xf64>
      %a651 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b652 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r653 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a651, %b652, %acc650 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a654 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b655 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r656 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a654, %b655, %r653 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r656, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv657 = memref.subview %el2[0, 0, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv658 = memref.subview %el2[0, 1, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa659 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r660 = vector.transfer_read %sv657[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %m661 = arith.subf %r660, %fa659 : vector<8x8xf64>
      vector.transfer_write %m661, %sv657[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
      %r662 = vector.transfer_read %sv658[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %p663 = arith.addf %r662, %fa659 : vector<8x8xf64>
      vector.transfer_write %p663, %sv658[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    }
    gpu.barrier
    %sv664 = memref.subview %v3613[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv665 = memref.subview %v3614[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv666 = memref.subview %el6[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv667 = memref.subview %el7[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv668 = memref.subview %el8[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc669 = arith.constant dense<0.0> : vector<8x8xf64>
      %a670 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b671 = vector.transfer_read %sv664[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r672 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a670, %b671, %acc669 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a673 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b674 = vector.transfer_read %sv664[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r675 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a673, %b674, %r672 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g676 = vector.transfer_read %sv667[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt677 = arith.mulf %r675, %g676 : vector<8x8xf64>
      %acc678 = arith.constant dense<0.0> : vector<8x8xf64>
      %a679 = vector.transfer_read %sv664[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b680 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r681 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a679, %b680, %acc678 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a682 = vector.transfer_read %sv664[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b683 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r684 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a682, %b683, %r681 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g685 = vector.transfer_read %sv666[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt686 = arith.mulf %r684, %g685 : vector<8x8xf64>
      %dv687 = vector.transfer_read %sv665[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g688 = vector.transfer_read %sv668[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t689 = arith.mulf %dv687, %g688 : vector<8x8xf64>
      %t2690 = arith.addf %t689, %dt686 : vector<8x8xf64>
      %fl691 = arith.addf %t2690, %dt677 : vector<8x8xf64>
      vector.transfer_write %fl691, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc692 = arith.constant dense<0.0> : vector<8x8xf64>
      %a693 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b694 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r695 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a693, %b694, %acc692 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a696 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b697 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r698 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a696, %b697, %r695 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r698, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc699 = arith.constant dense<0.0> : vector<8x8xf64>
      %a700 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b701 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r702 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a700, %b701, %acc699 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a703 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b704 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r705 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a703, %b704, %r702 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r705, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv706 = memref.subview %el2[0, 1, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv707 = memref.subview %el2[0, 2, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa708 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r709 = vector.transfer_read %sv706[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %m710 = arith.subf %r709, %fa708 : vector<8x8xf64>
      vector.transfer_write %m710, %sv706[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
      %r711 = vector.transfer_read %sv707[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %p712 = arith.addf %r711, %fa708 : vector<8x8xf64>
      vector.transfer_write %p712, %sv707[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    }
    gpu.barrier
    %sv713 = memref.subview %v3613[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv714 = memref.subview %v3614[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv715 = memref.subview %el6[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv716 = memref.subview %el7[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv717 = memref.subview %el8[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc718 = arith.constant dense<0.0> : vector<8x8xf64>
      %a719 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b720 = vector.transfer_read %sv713[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r721 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a719, %b720, %acc718 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a722 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b723 = vector.transfer_read %sv713[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r724 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a722, %b723, %r721 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g725 = vector.transfer_read %sv716[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt726 = arith.mulf %r724, %g725 : vector<8x8xf64>
      %acc727 = arith.constant dense<0.0> : vector<8x8xf64>
      %a728 = vector.transfer_read %sv713[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b729 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r730 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a728, %b729, %acc727 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a731 = vector.transfer_read %sv713[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b732 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r733 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a731, %b732, %r730 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g734 = vector.transfer_read %sv715[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt735 = arith.mulf %r733, %g734 : vector<8x8xf64>
      %dv736 = vector.transfer_read %sv714[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g737 = vector.transfer_read %sv717[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t738 = arith.mulf %dv736, %g737 : vector<8x8xf64>
      %t2739 = arith.addf %t738, %dt735 : vector<8x8xf64>
      %fl740 = arith.addf %t2739, %dt726 : vector<8x8xf64>
      vector.transfer_write %fl740, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc741 = arith.constant dense<0.0> : vector<8x8xf64>
      %a742 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b743 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r744 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a742, %b743, %acc741 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a745 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b746 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r747 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a745, %b746, %r744 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r747, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc748 = arith.constant dense<0.0> : vector<8x8xf64>
      %a749 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b750 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r751 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a749, %b750, %acc748 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a752 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b753 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r754 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a752, %b753, %r751 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r754, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv755 = memref.subview %el2[0, 2, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv756 = memref.subview %el2[0, 3, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa757 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r758 = vector.transfer_read %sv755[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %m759 = arith.subf %r758, %fa757 : vector<8x8xf64>
      vector.transfer_write %m759, %sv755[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
      %r760 = vector.transfer_read %sv756[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %p761 = arith.addf %r760, %fa757 : vector<8x8xf64>
      vector.transfer_write %p761, %sv756[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    }
    gpu.barrier
    %sv762 = memref.subview %v3613[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv763 = memref.subview %v3614[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv764 = memref.subview %el6[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv765 = memref.subview %el7[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv766 = memref.subview %el8[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc767 = arith.constant dense<0.0> : vector<8x8xf64>
      %a768 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b769 = vector.transfer_read %sv762[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r770 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a768, %b769, %acc767 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a771 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b772 = vector.transfer_read %sv762[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r773 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a771, %b772, %r770 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g774 = vector.transfer_read %sv765[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt775 = arith.mulf %r773, %g774 : vector<8x8xf64>
      %acc776 = arith.constant dense<0.0> : vector<8x8xf64>
      %a777 = vector.transfer_read %sv762[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b778 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r779 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a777, %b778, %acc776 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a780 = vector.transfer_read %sv762[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b781 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r782 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a780, %b781, %r779 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g783 = vector.transfer_read %sv764[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt784 = arith.mulf %r782, %g783 : vector<8x8xf64>
      %dv785 = vector.transfer_read %sv763[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g786 = vector.transfer_read %sv766[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t787 = arith.mulf %dv785, %g786 : vector<8x8xf64>
      %t2788 = arith.addf %t787, %dt784 : vector<8x8xf64>
      %fl789 = arith.addf %t2788, %dt775 : vector<8x8xf64>
      vector.transfer_write %fl789, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc790 = arith.constant dense<0.0> : vector<8x8xf64>
      %a791 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b792 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r793 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a791, %b792, %acc790 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a794 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b795 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r796 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a794, %b795, %r793 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r796, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc797 = arith.constant dense<0.0> : vector<8x8xf64>
      %a798 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b799 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r800 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a798, %b799, %acc797 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a801 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b802 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r803 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a801, %b802, %r800 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r803, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv804 = memref.subview %el2[0, 3, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv805 = memref.subview %el2[0, 4, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa806 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r807 = vector.transfer_read %sv804[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %m808 = arith.subf %r807, %fa806 : vector<8x8xf64>
      vector.transfer_write %m808, %sv804[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
      %r809 = vector.transfer_read %sv805[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %p810 = arith.addf %r809, %fa806 : vector<8x8xf64>
      vector.transfer_write %p810, %sv805[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    }
    gpu.barrier
    %sv811 = memref.subview %v3613[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv812 = memref.subview %v3614[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv813 = memref.subview %el6[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv814 = memref.subview %el7[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv815 = memref.subview %el8[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc816 = arith.constant dense<0.0> : vector<8x8xf64>
      %a817 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b818 = vector.transfer_read %sv811[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r819 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a817, %b818, %acc816 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a820 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b821 = vector.transfer_read %sv811[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r822 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a820, %b821, %r819 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g823 = vector.transfer_read %sv814[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt824 = arith.mulf %r822, %g823 : vector<8x8xf64>
      %acc825 = arith.constant dense<0.0> : vector<8x8xf64>
      %a826 = vector.transfer_read %sv811[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b827 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r828 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a826, %b827, %acc825 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a829 = vector.transfer_read %sv811[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b830 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r831 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a829, %b830, %r828 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g832 = vector.transfer_read %sv813[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt833 = arith.mulf %r831, %g832 : vector<8x8xf64>
      %dv834 = vector.transfer_read %sv812[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g835 = vector.transfer_read %sv815[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t836 = arith.mulf %dv834, %g835 : vector<8x8xf64>
      %t2837 = arith.addf %t836, %dt833 : vector<8x8xf64>
      %fl838 = arith.addf %t2837, %dt824 : vector<8x8xf64>
      vector.transfer_write %fl838, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc839 = arith.constant dense<0.0> : vector<8x8xf64>
      %a840 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b841 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r842 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a840, %b841, %acc839 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a843 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b844 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r845 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a843, %b844, %r842 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r845, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc846 = arith.constant dense<0.0> : vector<8x8xf64>
      %a847 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b848 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r849 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a847, %b848, %acc846 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a850 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b851 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r852 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a850, %b851, %r849 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r852, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv853 = memref.subview %el2[0, 4, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv854 = memref.subview %el2[0, 5, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa855 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r856 = vector.transfer_read %sv853[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %m857 = arith.subf %r856, %fa855 : vector<8x8xf64>
      vector.transfer_write %m857, %sv853[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
      %r858 = vector.transfer_read %sv854[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %p859 = arith.addf %r858, %fa855 : vector<8x8xf64>
      vector.transfer_write %p859, %sv854[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    }
    gpu.barrier
    %sv860 = memref.subview %v3613[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv861 = memref.subview %v3614[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv862 = memref.subview %el6[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv863 = memref.subview %el7[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv864 = memref.subview %el8[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc865 = arith.constant dense<0.0> : vector<8x8xf64>
      %a866 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b867 = vector.transfer_read %sv860[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r868 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a866, %b867, %acc865 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a869 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b870 = vector.transfer_read %sv860[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r871 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a869, %b870, %r868 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g872 = vector.transfer_read %sv863[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt873 = arith.mulf %r871, %g872 : vector<8x8xf64>
      %acc874 = arith.constant dense<0.0> : vector<8x8xf64>
      %a875 = vector.transfer_read %sv860[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b876 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r877 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a875, %b876, %acc874 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a878 = vector.transfer_read %sv860[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b879 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r880 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a878, %b879, %r877 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g881 = vector.transfer_read %sv862[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt882 = arith.mulf %r880, %g881 : vector<8x8xf64>
      %dv883 = vector.transfer_read %sv861[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g884 = vector.transfer_read %sv864[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t885 = arith.mulf %dv883, %g884 : vector<8x8xf64>
      %t2886 = arith.addf %t885, %dt882 : vector<8x8xf64>
      %fl887 = arith.addf %t2886, %dt873 : vector<8x8xf64>
      vector.transfer_write %fl887, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc888 = arith.constant dense<0.0> : vector<8x8xf64>
      %a889 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b890 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r891 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a889, %b890, %acc888 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a892 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b893 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r894 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a892, %b893, %r891 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r894, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc895 = arith.constant dense<0.0> : vector<8x8xf64>
      %a896 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b897 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r898 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a896, %b897, %acc895 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a899 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b900 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r901 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a899, %b900, %r898 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r901, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv902 = memref.subview %el2[0, 5, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv903 = memref.subview %el2[0, 6, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa904 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r905 = vector.transfer_read %sv902[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %m906 = arith.subf %r905, %fa904 : vector<8x8xf64>
      vector.transfer_write %m906, %sv902[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
      %r907 = vector.transfer_read %sv903[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %p908 = arith.addf %r907, %fa904 : vector<8x8xf64>
      vector.transfer_write %p908, %sv903[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    }
    gpu.barrier
    %sv909 = memref.subview %v3613[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv910 = memref.subview %v3614[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv911 = memref.subview %el6[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv912 = memref.subview %el7[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv913 = memref.subview %el8[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc914 = arith.constant dense<0.0> : vector<8x8xf64>
      %a915 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b916 = vector.transfer_read %sv909[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r917 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a915, %b916, %acc914 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a918 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b919 = vector.transfer_read %sv909[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r920 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a918, %b919, %r917 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g921 = vector.transfer_read %sv912[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt922 = arith.mulf %r920, %g921 : vector<8x8xf64>
      %acc923 = arith.constant dense<0.0> : vector<8x8xf64>
      %a924 = vector.transfer_read %sv909[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b925 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r926 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a924, %b925, %acc923 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a927 = vector.transfer_read %sv909[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b928 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r929 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a927, %b928, %r926 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g930 = vector.transfer_read %sv911[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt931 = arith.mulf %r929, %g930 : vector<8x8xf64>
      %dv932 = vector.transfer_read %sv910[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g933 = vector.transfer_read %sv913[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t934 = arith.mulf %dv932, %g933 : vector<8x8xf64>
      %t2935 = arith.addf %t934, %dt931 : vector<8x8xf64>
      %fl936 = arith.addf %t2935, %dt922 : vector<8x8xf64>
      vector.transfer_write %fl936, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc937 = arith.constant dense<0.0> : vector<8x8xf64>
      %a938 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b939 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r940 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a938, %b939, %acc937 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a941 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b942 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r943 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a941, %b942, %r940 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r943, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc944 = arith.constant dense<0.0> : vector<8x8xf64>
      %a945 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b946 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r947 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a945, %b946, %acc944 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a948 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b949 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r950 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a948, %b949, %r947 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r950, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv951 = memref.subview %el2[0, 6, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv952 = memref.subview %el2[0, 7, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa953 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r954 = vector.transfer_read %sv951[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %m955 = arith.subf %r954, %fa953 : vector<8x8xf64>
      vector.transfer_write %m955, %sv951[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
      %r956 = vector.transfer_read %sv952[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %p957 = arith.addf %r956, %fa953 : vector<8x8xf64>
      vector.transfer_write %p957, %sv952[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    }
    gpu.barrier
    %sv958 = memref.subview %el1[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc959 = arith.constant dense<0.0> : vector<8x8xf64>
      %a960 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b961 = vector.transfer_read %sv958[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r962 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a960, %b961, %acc959 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a963 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b964 = vector.transfer_read %sv958[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r965 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a963, %b964, %r962 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r965, %interp_all[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv966 = memref.subview %el1[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc967 = arith.constant dense<0.0> : vector<8x8xf64>
      %a968 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b969 = vector.transfer_read %sv966[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r970 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a968, %b969, %acc967 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a971 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b972 = vector.transfer_read %sv966[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r973 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a971, %b972, %r970 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r973, %interp_all[%c0, %c8] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv974 = memref.subview %el1[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc975 = arith.constant dense<0.0> : vector<8x8xf64>
      %a976 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b977 = vector.transfer_read %sv974[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r978 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a976, %b977, %acc975 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a979 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b980 = vector.transfer_read %sv974[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r981 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a979, %b980, %r978 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r981, %interp_all[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv982 = memref.subview %el1[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc983 = arith.constant dense<0.0> : vector<8x8xf64>
      %a984 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b985 = vector.transfer_read %sv982[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r986 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a984, %b985, %acc983 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a987 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b988 = vector.transfer_read %sv982[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r989 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a987, %b988, %r986 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r989, %interp_all[%c0, %c24] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv990 = memref.subview %el1[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc991 = arith.constant dense<0.0> : vector<8x8xf64>
      %a992 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b993 = vector.transfer_read %sv990[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r994 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a992, %b993, %acc991 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a995 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b996 = vector.transfer_read %sv990[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r997 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a995, %b996, %r994 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r997, %interp_all[%c0, %c32] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv998 = memref.subview %el1[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc999 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1000 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1001 = vector.transfer_read %sv998[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1002 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1000, %b1001, %acc999 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1003 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1004 = vector.transfer_read %sv998[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1005 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1003, %b1004, %r1002 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1005, %interp_all[%c0, %c40] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1006 = memref.subview %el1[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1007 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1008 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1009 = vector.transfer_read %sv1006[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1010 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1008, %b1009, %acc1007 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1011 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1012 = vector.transfer_read %sv1006[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1013 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1011, %b1012, %r1010 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1013, %interp_all[%c0, %c48] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1014 = memref.subview %el1[7, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1015 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1016 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1017 = vector.transfer_read %sv1014[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1018 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1016, %b1017, %acc1015 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1019 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1020 = vector.transfer_read %sv1014[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1021 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1019, %b1020, %r1018 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1021, %interp_all[%c0, %c56] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1022 = memref.subview %el1[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1023 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1024 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1025 = vector.transfer_read %sv1022[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1026 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1024, %b1025, %acc1023 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1027 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1028 = vector.transfer_read %sv1022[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1029 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1027, %b1028, %r1026 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1029, %deriv_all[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1030 = memref.subview %el1[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1031 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1032 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1033 = vector.transfer_read %sv1030[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1034 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1032, %b1033, %acc1031 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1035 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1036 = vector.transfer_read %sv1030[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1037 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1035, %b1036, %r1034 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1037, %deriv_all[%c0, %c8] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1038 = memref.subview %el1[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1039 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1040 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1041 = vector.transfer_read %sv1038[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1042 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1040, %b1041, %acc1039 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1043 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1044 = vector.transfer_read %sv1038[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1045 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1043, %b1044, %r1042 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1045, %deriv_all[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1046 = memref.subview %el1[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1047 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1048 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1049 = vector.transfer_read %sv1046[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1050 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1048, %b1049, %acc1047 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1051 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1052 = vector.transfer_read %sv1046[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1053 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1051, %b1052, %r1050 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1053, %deriv_all[%c0, %c24] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1054 = memref.subview %el1[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1055 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1056 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1057 = vector.transfer_read %sv1054[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1058 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1056, %b1057, %acc1055 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1059 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1060 = vector.transfer_read %sv1054[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1061 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1059, %b1060, %r1058 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1061, %deriv_all[%c0, %c32] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1062 = memref.subview %el1[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1063 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1064 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1065 = vector.transfer_read %sv1062[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1066 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1064, %b1065, %acc1063 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1067 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1068 = vector.transfer_read %sv1062[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1069 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1067, %b1068, %r1066 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1069, %deriv_all[%c0, %c40] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1070 = memref.subview %el1[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1071 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1072 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1073 = vector.transfer_read %sv1070[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1074 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1072, %b1073, %acc1071 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1075 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1076 = vector.transfer_read %sv1070[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1077 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1075, %b1076, %r1074 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1077, %deriv_all[%c0, %c48] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1078 = memref.subview %el1[7, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1079 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1080 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1081 = vector.transfer_read %sv1078[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1082 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1080, %b1081, %acc1079 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1083 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1084 = vector.transfer_read %sv1078[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1085 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1083, %b1084, %r1082 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1085, %deriv_all[%c0, %c56] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %v31086 = memref.reinterpret_cast %interp_all to offset: [0], sizes: [8, 8, 8], strides: [64, 8, 1] : memref<8x64xf64, #gpu.address_space<workgroup>> to memref<8x8x8xf64, #gpu.address_space<workgroup>>
    %v31087 = memref.reinterpret_cast %deriv_all to offset: [0], sizes: [8, 8, 8], strides: [64, 8, 1] : memref<8x64xf64, #gpu.address_space<workgroup>> to memref<8x8x8xf64, #gpu.address_space<workgroup>>
    %sv1088 = memref.subview %v31086[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>
    %sv1089 = memref.subview %v31087[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>
    %sv1090 = memref.subview %el9[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv1091 = memref.subview %el10[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv1092 = memref.subview %el11[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1093 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1094 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1095 = vector.transfer_read %sv1088[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1096 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1094, %b1095, %acc1093 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1097 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1098 = vector.transfer_read %sv1088[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1099 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1097, %b1098, %r1096 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1100 = vector.transfer_read %sv1091[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt1101 = arith.mulf %r1099, %g1100 : vector<8x8xf64>
      %acc1102 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1103 = vector.transfer_read %sv1088[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1104 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1105 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1103, %b1104, %acc1102 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1106 = vector.transfer_read %sv1088[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1107 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1108 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1106, %b1107, %r1105 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1109 = vector.transfer_read %sv1090[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt1110 = arith.mulf %r1108, %g1109 : vector<8x8xf64>
      %dv1111 = vector.transfer_read %sv1089[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g1112 = vector.transfer_read %sv1092[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t1113 = arith.mulf %dv1111, %g1112 : vector<8x8xf64>
      %t21114 = arith.addf %t1113, %dt1110 : vector<8x8xf64>
      %fl1115 = arith.addf %t21114, %dt1101 : vector<8x8xf64>
      vector.transfer_write %fl1115, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1116 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1117 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1118 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1119 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1117, %b1118, %acc1116 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1120 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1121 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1122 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1120, %b1121, %r1119 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1122, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1123 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1124 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1125 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1126 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1124, %b1125, %acc1123 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1127 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1128 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1129 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1127, %b1128, %r1126 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1129, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1130 = memref.subview %el2[0, 0, 0] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv1131 = memref.subview %el2[0, 0, 1] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1132 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1133 = vector.transfer_read %sv1130[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %m1134 = arith.subf %r1133, %fa1132 : vector<8x8xf64>
      vector.transfer_write %m1134, %sv1130[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
      %r1135 = vector.transfer_read %sv1131[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %p1136 = arith.addf %r1135, %fa1132 : vector<8x8xf64>
      vector.transfer_write %p1136, %sv1131[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    }
    gpu.barrier
    %sv1137 = memref.subview %v31086[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv1138 = memref.subview %v31087[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv1139 = memref.subview %el9[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv1140 = memref.subview %el10[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv1141 = memref.subview %el11[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1142 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1143 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1144 = vector.transfer_read %sv1137[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1145 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1143, %b1144, %acc1142 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1146 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1147 = vector.transfer_read %sv1137[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1148 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1146, %b1147, %r1145 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1149 = vector.transfer_read %sv1140[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt1150 = arith.mulf %r1148, %g1149 : vector<8x8xf64>
      %acc1151 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1152 = vector.transfer_read %sv1137[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1153 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1154 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1152, %b1153, %acc1151 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1155 = vector.transfer_read %sv1137[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1156 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1157 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1155, %b1156, %r1154 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1158 = vector.transfer_read %sv1139[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt1159 = arith.mulf %r1157, %g1158 : vector<8x8xf64>
      %dv1160 = vector.transfer_read %sv1138[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g1161 = vector.transfer_read %sv1141[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t1162 = arith.mulf %dv1160, %g1161 : vector<8x8xf64>
      %t21163 = arith.addf %t1162, %dt1159 : vector<8x8xf64>
      %fl1164 = arith.addf %t21163, %dt1150 : vector<8x8xf64>
      vector.transfer_write %fl1164, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1165 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1166 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1167 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1168 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1166, %b1167, %acc1165 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1169 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1170 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1171 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1169, %b1170, %r1168 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1171, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1172 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1173 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1174 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1175 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1173, %b1174, %acc1172 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1176 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1177 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1178 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1176, %b1177, %r1175 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1178, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1179 = memref.subview %el2[0, 0, 1] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv1180 = memref.subview %el2[0, 0, 2] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1181 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1182 = vector.transfer_read %sv1179[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %m1183 = arith.subf %r1182, %fa1181 : vector<8x8xf64>
      vector.transfer_write %m1183, %sv1179[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
      %r1184 = vector.transfer_read %sv1180[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %p1185 = arith.addf %r1184, %fa1181 : vector<8x8xf64>
      vector.transfer_write %p1185, %sv1180[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    }
    gpu.barrier
    %sv1186 = memref.subview %v31086[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv1187 = memref.subview %v31087[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv1188 = memref.subview %el9[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv1189 = memref.subview %el10[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv1190 = memref.subview %el11[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1191 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1192 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1193 = vector.transfer_read %sv1186[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1194 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1192, %b1193, %acc1191 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1195 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1196 = vector.transfer_read %sv1186[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1197 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1195, %b1196, %r1194 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1198 = vector.transfer_read %sv1189[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt1199 = arith.mulf %r1197, %g1198 : vector<8x8xf64>
      %acc1200 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1201 = vector.transfer_read %sv1186[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1202 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1203 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1201, %b1202, %acc1200 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1204 = vector.transfer_read %sv1186[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1205 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1206 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1204, %b1205, %r1203 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1207 = vector.transfer_read %sv1188[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt1208 = arith.mulf %r1206, %g1207 : vector<8x8xf64>
      %dv1209 = vector.transfer_read %sv1187[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g1210 = vector.transfer_read %sv1190[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t1211 = arith.mulf %dv1209, %g1210 : vector<8x8xf64>
      %t21212 = arith.addf %t1211, %dt1208 : vector<8x8xf64>
      %fl1213 = arith.addf %t21212, %dt1199 : vector<8x8xf64>
      vector.transfer_write %fl1213, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1214 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1215 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1216 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1217 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1215, %b1216, %acc1214 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1218 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1219 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1220 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1218, %b1219, %r1217 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1220, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1221 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1222 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1223 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1224 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1222, %b1223, %acc1221 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1225 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1226 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1227 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1225, %b1226, %r1224 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1227, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1228 = memref.subview %el2[0, 0, 2] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv1229 = memref.subview %el2[0, 0, 3] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1230 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1231 = vector.transfer_read %sv1228[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %m1232 = arith.subf %r1231, %fa1230 : vector<8x8xf64>
      vector.transfer_write %m1232, %sv1228[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
      %r1233 = vector.transfer_read %sv1229[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %p1234 = arith.addf %r1233, %fa1230 : vector<8x8xf64>
      vector.transfer_write %p1234, %sv1229[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    }
    gpu.barrier
    %sv1235 = memref.subview %v31086[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv1236 = memref.subview %v31087[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv1237 = memref.subview %el9[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv1238 = memref.subview %el10[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv1239 = memref.subview %el11[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1240 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1241 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1242 = vector.transfer_read %sv1235[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1243 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1241, %b1242, %acc1240 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1244 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1245 = vector.transfer_read %sv1235[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1246 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1244, %b1245, %r1243 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1247 = vector.transfer_read %sv1238[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt1248 = arith.mulf %r1246, %g1247 : vector<8x8xf64>
      %acc1249 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1250 = vector.transfer_read %sv1235[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1251 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1252 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1250, %b1251, %acc1249 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1253 = vector.transfer_read %sv1235[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1254 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1255 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1253, %b1254, %r1252 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1256 = vector.transfer_read %sv1237[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt1257 = arith.mulf %r1255, %g1256 : vector<8x8xf64>
      %dv1258 = vector.transfer_read %sv1236[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g1259 = vector.transfer_read %sv1239[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t1260 = arith.mulf %dv1258, %g1259 : vector<8x8xf64>
      %t21261 = arith.addf %t1260, %dt1257 : vector<8x8xf64>
      %fl1262 = arith.addf %t21261, %dt1248 : vector<8x8xf64>
      vector.transfer_write %fl1262, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1263 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1264 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1265 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1266 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1264, %b1265, %acc1263 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1267 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1268 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1269 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1267, %b1268, %r1266 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1269, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1270 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1271 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1272 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1273 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1271, %b1272, %acc1270 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1274 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1275 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1276 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1274, %b1275, %r1273 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1276, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1277 = memref.subview %el2[0, 0, 3] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv1278 = memref.subview %el2[0, 0, 4] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1279 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1280 = vector.transfer_read %sv1277[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %m1281 = arith.subf %r1280, %fa1279 : vector<8x8xf64>
      vector.transfer_write %m1281, %sv1277[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
      %r1282 = vector.transfer_read %sv1278[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %p1283 = arith.addf %r1282, %fa1279 : vector<8x8xf64>
      vector.transfer_write %p1283, %sv1278[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    }
    gpu.barrier
    %sv1284 = memref.subview %v31086[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv1285 = memref.subview %v31087[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv1286 = memref.subview %el9[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv1287 = memref.subview %el10[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv1288 = memref.subview %el11[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1289 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1290 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1291 = vector.transfer_read %sv1284[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1292 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1290, %b1291, %acc1289 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1293 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1294 = vector.transfer_read %sv1284[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1295 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1293, %b1294, %r1292 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1296 = vector.transfer_read %sv1287[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt1297 = arith.mulf %r1295, %g1296 : vector<8x8xf64>
      %acc1298 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1299 = vector.transfer_read %sv1284[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1300 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1301 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1299, %b1300, %acc1298 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1302 = vector.transfer_read %sv1284[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1303 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1304 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1302, %b1303, %r1301 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1305 = vector.transfer_read %sv1286[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt1306 = arith.mulf %r1304, %g1305 : vector<8x8xf64>
      %dv1307 = vector.transfer_read %sv1285[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g1308 = vector.transfer_read %sv1288[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t1309 = arith.mulf %dv1307, %g1308 : vector<8x8xf64>
      %t21310 = arith.addf %t1309, %dt1306 : vector<8x8xf64>
      %fl1311 = arith.addf %t21310, %dt1297 : vector<8x8xf64>
      vector.transfer_write %fl1311, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1312 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1313 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1314 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1315 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1313, %b1314, %acc1312 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1316 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1317 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1318 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1316, %b1317, %r1315 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1318, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1319 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1320 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1321 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1322 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1320, %b1321, %acc1319 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1323 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1324 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1325 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1323, %b1324, %r1322 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1325, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1326 = memref.subview %el2[0, 0, 4] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv1327 = memref.subview %el2[0, 0, 5] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1328 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1329 = vector.transfer_read %sv1326[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %m1330 = arith.subf %r1329, %fa1328 : vector<8x8xf64>
      vector.transfer_write %m1330, %sv1326[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
      %r1331 = vector.transfer_read %sv1327[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %p1332 = arith.addf %r1331, %fa1328 : vector<8x8xf64>
      vector.transfer_write %p1332, %sv1327[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    }
    gpu.barrier
    %sv1333 = memref.subview %v31086[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv1334 = memref.subview %v31087[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv1335 = memref.subview %el9[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv1336 = memref.subview %el10[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv1337 = memref.subview %el11[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1338 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1339 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1340 = vector.transfer_read %sv1333[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1341 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1339, %b1340, %acc1338 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1342 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1343 = vector.transfer_read %sv1333[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1344 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1342, %b1343, %r1341 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1345 = vector.transfer_read %sv1336[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt1346 = arith.mulf %r1344, %g1345 : vector<8x8xf64>
      %acc1347 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1348 = vector.transfer_read %sv1333[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1349 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1350 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1348, %b1349, %acc1347 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1351 = vector.transfer_read %sv1333[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1352 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1353 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1351, %b1352, %r1350 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1354 = vector.transfer_read %sv1335[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt1355 = arith.mulf %r1353, %g1354 : vector<8x8xf64>
      %dv1356 = vector.transfer_read %sv1334[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g1357 = vector.transfer_read %sv1337[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t1358 = arith.mulf %dv1356, %g1357 : vector<8x8xf64>
      %t21359 = arith.addf %t1358, %dt1355 : vector<8x8xf64>
      %fl1360 = arith.addf %t21359, %dt1346 : vector<8x8xf64>
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
    %sv1375 = memref.subview %el2[0, 0, 5] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv1376 = memref.subview %el2[0, 0, 6] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1377 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1378 = vector.transfer_read %sv1375[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %m1379 = arith.subf %r1378, %fa1377 : vector<8x8xf64>
      vector.transfer_write %m1379, %sv1375[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
      %r1380 = vector.transfer_read %sv1376[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %p1381 = arith.addf %r1380, %fa1377 : vector<8x8xf64>
      vector.transfer_write %p1381, %sv1376[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    }
    gpu.barrier
    %sv1382 = memref.subview %v31086[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv1383 = memref.subview %v31087[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv1384 = memref.subview %el9[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv1385 = memref.subview %el10[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv1386 = memref.subview %el11[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1387 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1388 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1389 = vector.transfer_read %sv1382[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1390 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1388, %b1389, %acc1387 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1391 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1392 = vector.transfer_read %sv1382[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1393 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1391, %b1392, %r1390 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1394 = vector.transfer_read %sv1385[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt1395 = arith.mulf %r1393, %g1394 : vector<8x8xf64>
      %acc1396 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1397 = vector.transfer_read %sv1382[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1398 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1399 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1397, %b1398, %acc1396 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1400 = vector.transfer_read %sv1382[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1401 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1402 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1400, %b1401, %r1399 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1403 = vector.transfer_read %sv1384[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt1404 = arith.mulf %r1402, %g1403 : vector<8x8xf64>
      %dv1405 = vector.transfer_read %sv1383[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g1406 = vector.transfer_read %sv1386[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t1407 = arith.mulf %dv1405, %g1406 : vector<8x8xf64>
      %t21408 = arith.addf %t1407, %dt1404 : vector<8x8xf64>
      %fl1409 = arith.addf %t21408, %dt1395 : vector<8x8xf64>
      vector.transfer_write %fl1409, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1410 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1411 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1412 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1413 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1411, %b1412, %acc1410 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1414 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1415 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1416 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1414, %b1415, %r1413 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1416, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1417 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1418 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1419 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1420 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1418, %b1419, %acc1417 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1421 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1422 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1423 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1421, %b1422, %r1420 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1423, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1424 = memref.subview %el2[0, 0, 6] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv1425 = memref.subview %el2[0, 0, 7] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1426 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1427 = vector.transfer_read %sv1424[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %m1428 = arith.subf %r1427, %fa1426 : vector<8x8xf64>
      vector.transfer_write %m1428, %sv1424[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
      %r1429 = vector.transfer_read %sv1425[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %p1430 = arith.addf %r1429, %fa1426 : vector<8x8xf64>
      vector.transfer_write %p1430, %sv1425[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    }
    gpu.barrier
    gpu.return
  }
}
