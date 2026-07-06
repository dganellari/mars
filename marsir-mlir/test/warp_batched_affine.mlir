#a = affine_map<(m, n, k) -> (m, k)>
#b = affine_map<(m, n, k) -> (k, n)>
#bt = affine_map<(m, n, k) -> (n, k)>
#c = affine_map<(m, n, k) -> (m, n)>
gpu.module @mir_kernels [#nvvm.target<chip = "sm_90", O = 3>] {
  gpu.func @full_batched_affine(%U: memref<?x8x8x8xf64>, %Btil: memref<8x8xf64>, %Dtil: memref<8x8xf64>, %Dm: memref<8x8xf64>, %W: memref<8x8xf64>, %g00all: memref<?x8x8xf64>, %g01all: memref<?x8x8xf64>, %g02all: memref<?x8x8xf64>, %g10all: memref<?x8x8xf64>, %g11all: memref<?x8x8xf64>, %g12all: memref<?x8x8xf64>, %g20all: memref<?x8x8xf64>, %g21all: memref<?x8x8xf64>, %g22all: memref<?x8x8xf64>, %Y: memref<?x8x8x8xf64>) workgroup(%interp_all: memref<8x64xf64, #gpu.address_space<workgroup>>, %deriv_all: memref<8x64xf64, #gpu.address_space<workgroup>>, %dt1g: memref<8x8xf64, #gpu.address_space<workgroup>>, %dt2g: memref<8x8xf64, #gpu.address_space<workgroup>>, %flux: memref<8x8xf64, #gpu.address_space<workgroup>>, %tmp: memref<8x8xf64, #gpu.address_space<workgroup>>, %intf: memref<8x8xf64, #gpu.address_space<workgroup>>) kernel {
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
    %el3 = memref.subview %g00all[%e, 0, 0] [1, 8, 8] [1, 1, 1] : memref<?x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %el4 = memref.subview %g01all[%e, 0, 0] [1, 8, 8] [1, 1, 1] : memref<?x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %el5 = memref.subview %g02all[%e, 0, 0] [1, 8, 8] [1, 1, 1] : memref<?x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %el6 = memref.subview %g10all[%e, 0, 0] [1, 8, 8] [1, 1, 1] : memref<?x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %el7 = memref.subview %g11all[%e, 0, 0] [1, 8, 8] [1, 1, 1] : memref<?x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %el8 = memref.subview %g12all[%e, 0, 0] [1, 8, 8] [1, 1, 1] : memref<?x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %el9 = memref.subview %g20all[%e, 0, 0] [1, 8, 8] [1, 1, 1] : memref<?x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %el10 = memref.subview %g21all[%e, 0, 0] [1, 8, 8] [1, 1, 1] : memref<?x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %el11 = memref.subview %g22all[%e, 0, 0] [1, 8, 8] [1, 1, 1] : memref<?x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: ?>>
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
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc144 = arith.constant dense<0.0> : vector<8x8xf64>
      %a145 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b146 = vector.transfer_read %sv142[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r147 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a145, %b146, %acc144 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a148 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b149 = vector.transfer_read %sv142[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r150 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a148, %b149, %r147 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g151 = vector.transfer_read %el4[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc152 = arith.mulf %r150, %g151 : vector<8x8xf64>
      vector.transfer_write %sc152, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc153 = arith.constant dense<0.0> : vector<8x8xf64>
      %a154 = vector.transfer_read %sv142[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b155 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r156 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a154, %b155, %acc153 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a157 = vector.transfer_read %sv142[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b158 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r159 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a157, %b158, %r156 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g160 = vector.transfer_read %el3[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc161 = arith.mulf %r159, %g160 : vector<8x8xf64>
      vector.transfer_write %sc161, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d162 = vector.transfer_read %sv143[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2163 = vector.transfer_read %el5[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t164 = arith.mulf %d162, %g2163 : vector<8x8xf64>
      %a165 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2166 = arith.addf %t164, %a165 : vector<8x8xf64>
      %b167 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl168 = arith.addf %t2166, %b167 : vector<8x8xf64>
      vector.transfer_write %fl168, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc169 = arith.constant dense<0.0> : vector<8x8xf64>
      %a170 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b171 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r172 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a170, %b171, %acc169 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a173 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b174 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r175 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a173, %b174, %r172 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r175, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc176 = arith.constant dense<0.0> : vector<8x8xf64>
      %a177 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b178 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r179 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a177, %b178, %acc176 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a180 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b181 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r182 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a180, %b181, %r179 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r182, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv183 = memref.subview %el2[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv184 = memref.subview %el2[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa185 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r186 = vector.transfer_read %sv183[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %m187 = arith.subf %r186, %fa185 : vector<8x8xf64>
      vector.transfer_write %m187, %sv183[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      %fa188 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r189 = vector.transfer_read %sv184[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %p190 = arith.addf %r189, %fa188 : vector<8x8xf64>
      vector.transfer_write %p190, %sv184[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    }
    gpu.barrier
    %sv191 = memref.subview %v3140[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv192 = memref.subview %v3141[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc193 = arith.constant dense<0.0> : vector<8x8xf64>
      %a194 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b195 = vector.transfer_read %sv191[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r196 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a194, %b195, %acc193 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a197 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b198 = vector.transfer_read %sv191[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r199 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a197, %b198, %r196 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g200 = vector.transfer_read %el4[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc201 = arith.mulf %r199, %g200 : vector<8x8xf64>
      vector.transfer_write %sc201, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc202 = arith.constant dense<0.0> : vector<8x8xf64>
      %a203 = vector.transfer_read %sv191[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b204 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r205 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a203, %b204, %acc202 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a206 = vector.transfer_read %sv191[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b207 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r208 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a206, %b207, %r205 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g209 = vector.transfer_read %el3[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc210 = arith.mulf %r208, %g209 : vector<8x8xf64>
      vector.transfer_write %sc210, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d211 = vector.transfer_read %sv192[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2212 = vector.transfer_read %el5[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t213 = arith.mulf %d211, %g2212 : vector<8x8xf64>
      %a214 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2215 = arith.addf %t213, %a214 : vector<8x8xf64>
      %b216 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl217 = arith.addf %t2215, %b216 : vector<8x8xf64>
      vector.transfer_write %fl217, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc218 = arith.constant dense<0.0> : vector<8x8xf64>
      %a219 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b220 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r221 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a219, %b220, %acc218 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a222 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b223 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r224 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a222, %b223, %r221 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r224, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc225 = arith.constant dense<0.0> : vector<8x8xf64>
      %a226 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b227 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r228 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a226, %b227, %acc225 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a229 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b230 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r231 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a229, %b230, %r228 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r231, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv232 = memref.subview %el2[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv233 = memref.subview %el2[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa234 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r235 = vector.transfer_read %sv232[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %m236 = arith.subf %r235, %fa234 : vector<8x8xf64>
      vector.transfer_write %m236, %sv232[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      %fa237 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r238 = vector.transfer_read %sv233[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %p239 = arith.addf %r238, %fa237 : vector<8x8xf64>
      vector.transfer_write %p239, %sv233[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    }
    gpu.barrier
    %sv240 = memref.subview %v3140[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv241 = memref.subview %v3141[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc242 = arith.constant dense<0.0> : vector<8x8xf64>
      %a243 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b244 = vector.transfer_read %sv240[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r245 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a243, %b244, %acc242 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a246 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b247 = vector.transfer_read %sv240[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r248 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a246, %b247, %r245 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g249 = vector.transfer_read %el4[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc250 = arith.mulf %r248, %g249 : vector<8x8xf64>
      vector.transfer_write %sc250, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc251 = arith.constant dense<0.0> : vector<8x8xf64>
      %a252 = vector.transfer_read %sv240[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b253 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r254 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a252, %b253, %acc251 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a255 = vector.transfer_read %sv240[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b256 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r257 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a255, %b256, %r254 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g258 = vector.transfer_read %el3[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc259 = arith.mulf %r257, %g258 : vector<8x8xf64>
      vector.transfer_write %sc259, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d260 = vector.transfer_read %sv241[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2261 = vector.transfer_read %el5[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t262 = arith.mulf %d260, %g2261 : vector<8x8xf64>
      %a263 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2264 = arith.addf %t262, %a263 : vector<8x8xf64>
      %b265 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl266 = arith.addf %t2264, %b265 : vector<8x8xf64>
      vector.transfer_write %fl266, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc267 = arith.constant dense<0.0> : vector<8x8xf64>
      %a268 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b269 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r270 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a268, %b269, %acc267 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a271 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b272 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r273 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a271, %b272, %r270 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r273, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc274 = arith.constant dense<0.0> : vector<8x8xf64>
      %a275 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b276 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r277 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a275, %b276, %acc274 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a278 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b279 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r280 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a278, %b279, %r277 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r280, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv281 = memref.subview %el2[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv282 = memref.subview %el2[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa283 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r284 = vector.transfer_read %sv281[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %m285 = arith.subf %r284, %fa283 : vector<8x8xf64>
      vector.transfer_write %m285, %sv281[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      %fa286 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r287 = vector.transfer_read %sv282[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %p288 = arith.addf %r287, %fa286 : vector<8x8xf64>
      vector.transfer_write %p288, %sv282[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    }
    gpu.barrier
    %sv289 = memref.subview %v3140[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv290 = memref.subview %v3141[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc291 = arith.constant dense<0.0> : vector<8x8xf64>
      %a292 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b293 = vector.transfer_read %sv289[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r294 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a292, %b293, %acc291 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a295 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b296 = vector.transfer_read %sv289[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r297 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a295, %b296, %r294 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g298 = vector.transfer_read %el4[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc299 = arith.mulf %r297, %g298 : vector<8x8xf64>
      vector.transfer_write %sc299, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc300 = arith.constant dense<0.0> : vector<8x8xf64>
      %a301 = vector.transfer_read %sv289[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b302 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r303 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a301, %b302, %acc300 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a304 = vector.transfer_read %sv289[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b305 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r306 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a304, %b305, %r303 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g307 = vector.transfer_read %el3[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc308 = arith.mulf %r306, %g307 : vector<8x8xf64>
      vector.transfer_write %sc308, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d309 = vector.transfer_read %sv290[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2310 = vector.transfer_read %el5[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t311 = arith.mulf %d309, %g2310 : vector<8x8xf64>
      %a312 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2313 = arith.addf %t311, %a312 : vector<8x8xf64>
      %b314 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl315 = arith.addf %t2313, %b314 : vector<8x8xf64>
      vector.transfer_write %fl315, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc316 = arith.constant dense<0.0> : vector<8x8xf64>
      %a317 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b318 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r319 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a317, %b318, %acc316 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a320 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b321 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r322 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a320, %b321, %r319 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r322, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc323 = arith.constant dense<0.0> : vector<8x8xf64>
      %a324 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b325 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r326 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a324, %b325, %acc323 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a327 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b328 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r329 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a327, %b328, %r326 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r329, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv330 = memref.subview %el2[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv331 = memref.subview %el2[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa332 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r333 = vector.transfer_read %sv330[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %m334 = arith.subf %r333, %fa332 : vector<8x8xf64>
      vector.transfer_write %m334, %sv330[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      %fa335 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r336 = vector.transfer_read %sv331[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %p337 = arith.addf %r336, %fa335 : vector<8x8xf64>
      vector.transfer_write %p337, %sv331[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    }
    gpu.barrier
    %sv338 = memref.subview %v3140[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv339 = memref.subview %v3141[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc340 = arith.constant dense<0.0> : vector<8x8xf64>
      %a341 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b342 = vector.transfer_read %sv338[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r343 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a341, %b342, %acc340 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a344 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b345 = vector.transfer_read %sv338[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r346 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a344, %b345, %r343 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g347 = vector.transfer_read %el4[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc348 = arith.mulf %r346, %g347 : vector<8x8xf64>
      vector.transfer_write %sc348, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc349 = arith.constant dense<0.0> : vector<8x8xf64>
      %a350 = vector.transfer_read %sv338[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b351 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r352 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a350, %b351, %acc349 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a353 = vector.transfer_read %sv338[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b354 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r355 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a353, %b354, %r352 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g356 = vector.transfer_read %el3[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc357 = arith.mulf %r355, %g356 : vector<8x8xf64>
      vector.transfer_write %sc357, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d358 = vector.transfer_read %sv339[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2359 = vector.transfer_read %el5[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t360 = arith.mulf %d358, %g2359 : vector<8x8xf64>
      %a361 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2362 = arith.addf %t360, %a361 : vector<8x8xf64>
      %b363 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl364 = arith.addf %t2362, %b363 : vector<8x8xf64>
      vector.transfer_write %fl364, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc365 = arith.constant dense<0.0> : vector<8x8xf64>
      %a366 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b367 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r368 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a366, %b367, %acc365 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a369 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b370 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r371 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a369, %b370, %r368 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r371, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc372 = arith.constant dense<0.0> : vector<8x8xf64>
      %a373 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b374 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r375 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a373, %b374, %acc372 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a376 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b377 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r378 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a376, %b377, %r375 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r378, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv379 = memref.subview %el2[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv380 = memref.subview %el2[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa381 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r382 = vector.transfer_read %sv379[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %m383 = arith.subf %r382, %fa381 : vector<8x8xf64>
      vector.transfer_write %m383, %sv379[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      %fa384 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r385 = vector.transfer_read %sv380[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %p386 = arith.addf %r385, %fa384 : vector<8x8xf64>
      vector.transfer_write %p386, %sv380[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    }
    gpu.barrier
    %sv387 = memref.subview %v3140[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv388 = memref.subview %v3141[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc389 = arith.constant dense<0.0> : vector<8x8xf64>
      %a390 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b391 = vector.transfer_read %sv387[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r392 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a390, %b391, %acc389 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a393 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b394 = vector.transfer_read %sv387[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r395 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a393, %b394, %r392 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g396 = vector.transfer_read %el4[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc397 = arith.mulf %r395, %g396 : vector<8x8xf64>
      vector.transfer_write %sc397, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc398 = arith.constant dense<0.0> : vector<8x8xf64>
      %a399 = vector.transfer_read %sv387[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b400 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r401 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a399, %b400, %acc398 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a402 = vector.transfer_read %sv387[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b403 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r404 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a402, %b403, %r401 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g405 = vector.transfer_read %el3[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc406 = arith.mulf %r404, %g405 : vector<8x8xf64>
      vector.transfer_write %sc406, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d407 = vector.transfer_read %sv388[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2408 = vector.transfer_read %el5[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t409 = arith.mulf %d407, %g2408 : vector<8x8xf64>
      %a410 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2411 = arith.addf %t409, %a410 : vector<8x8xf64>
      %b412 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl413 = arith.addf %t2411, %b412 : vector<8x8xf64>
      vector.transfer_write %fl413, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc414 = arith.constant dense<0.0> : vector<8x8xf64>
      %a415 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b416 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r417 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a415, %b416, %acc414 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a418 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b419 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r420 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a418, %b419, %r417 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r420, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc421 = arith.constant dense<0.0> : vector<8x8xf64>
      %a422 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b423 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r424 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a422, %b423, %acc421 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a425 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b426 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r427 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a425, %b426, %r424 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r427, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv428 = memref.subview %el2[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv429 = memref.subview %el2[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa430 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r431 = vector.transfer_read %sv428[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %m432 = arith.subf %r431, %fa430 : vector<8x8xf64>
      vector.transfer_write %m432, %sv428[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      %fa433 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r434 = vector.transfer_read %sv429[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %p435 = arith.addf %r434, %fa433 : vector<8x8xf64>
      vector.transfer_write %p435, %sv429[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    }
    gpu.barrier
    %sv436 = memref.subview %v3140[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv437 = memref.subview %v3141[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc438 = arith.constant dense<0.0> : vector<8x8xf64>
      %a439 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b440 = vector.transfer_read %sv436[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r441 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a439, %b440, %acc438 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a442 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b443 = vector.transfer_read %sv436[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r444 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a442, %b443, %r441 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g445 = vector.transfer_read %el4[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc446 = arith.mulf %r444, %g445 : vector<8x8xf64>
      vector.transfer_write %sc446, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc447 = arith.constant dense<0.0> : vector<8x8xf64>
      %a448 = vector.transfer_read %sv436[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b449 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r450 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a448, %b449, %acc447 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a451 = vector.transfer_read %sv436[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b452 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r453 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a451, %b452, %r450 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g454 = vector.transfer_read %el3[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc455 = arith.mulf %r453, %g454 : vector<8x8xf64>
      vector.transfer_write %sc455, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d456 = vector.transfer_read %sv437[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2457 = vector.transfer_read %el5[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t458 = arith.mulf %d456, %g2457 : vector<8x8xf64>
      %a459 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2460 = arith.addf %t458, %a459 : vector<8x8xf64>
      %b461 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl462 = arith.addf %t2460, %b461 : vector<8x8xf64>
      vector.transfer_write %fl462, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc463 = arith.constant dense<0.0> : vector<8x8xf64>
      %a464 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b465 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r466 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a464, %b465, %acc463 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a467 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b468 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r469 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a467, %b468, %r466 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r469, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc470 = arith.constant dense<0.0> : vector<8x8xf64>
      %a471 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b472 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r473 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a471, %b472, %acc470 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a474 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b475 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r476 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a474, %b475, %r473 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r476, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv477 = memref.subview %el2[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv478 = memref.subview %el2[7, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa479 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r480 = vector.transfer_read %sv477[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %m481 = arith.subf %r480, %fa479 : vector<8x8xf64>
      vector.transfer_write %m481, %sv477[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      %fa482 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r483 = vector.transfer_read %sv478[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %p484 = arith.addf %r483, %fa482 : vector<8x8xf64>
      vector.transfer_write %p484, %sv478[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
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
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc617 = arith.constant dense<0.0> : vector<8x8xf64>
      %a618 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b619 = vector.transfer_read %sv615[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r620 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a618, %b619, %acc617 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a621 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b622 = vector.transfer_read %sv615[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r623 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a621, %b622, %r620 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g624 = vector.transfer_read %el7[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc625 = arith.mulf %r623, %g624 : vector<8x8xf64>
      vector.transfer_write %sc625, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc626 = arith.constant dense<0.0> : vector<8x8xf64>
      %a627 = vector.transfer_read %sv615[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b628 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r629 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a627, %b628, %acc626 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a630 = vector.transfer_read %sv615[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b631 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r632 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a630, %b631, %r629 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g633 = vector.transfer_read %el6[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc634 = arith.mulf %r632, %g633 : vector<8x8xf64>
      vector.transfer_write %sc634, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d635 = vector.transfer_read %sv616[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2636 = vector.transfer_read %el8[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t637 = arith.mulf %d635, %g2636 : vector<8x8xf64>
      %a638 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2639 = arith.addf %t637, %a638 : vector<8x8xf64>
      %b640 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl641 = arith.addf %t2639, %b640 : vector<8x8xf64>
      vector.transfer_write %fl641, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc642 = arith.constant dense<0.0> : vector<8x8xf64>
      %a643 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b644 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r645 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a643, %b644, %acc642 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a646 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b647 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r648 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a646, %b647, %r645 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r648, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc649 = arith.constant dense<0.0> : vector<8x8xf64>
      %a650 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b651 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r652 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a650, %b651, %acc649 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a653 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b654 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r655 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a653, %b654, %r652 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r655, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv656 = memref.subview %el2[0, 0, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv657 = memref.subview %el2[0, 1, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa658 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r659 = vector.transfer_read %sv656[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %m660 = arith.subf %r659, %fa658 : vector<8x8xf64>
      vector.transfer_write %m660, %sv656[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
      %fa661 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r662 = vector.transfer_read %sv657[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %p663 = arith.addf %r662, %fa661 : vector<8x8xf64>
      vector.transfer_write %p663, %sv657[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    }
    gpu.barrier
    %sv664 = memref.subview %v3613[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv665 = memref.subview %v3614[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc666 = arith.constant dense<0.0> : vector<8x8xf64>
      %a667 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b668 = vector.transfer_read %sv664[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r669 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a667, %b668, %acc666 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a670 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b671 = vector.transfer_read %sv664[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r672 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a670, %b671, %r669 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g673 = vector.transfer_read %el7[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc674 = arith.mulf %r672, %g673 : vector<8x8xf64>
      vector.transfer_write %sc674, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc675 = arith.constant dense<0.0> : vector<8x8xf64>
      %a676 = vector.transfer_read %sv664[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b677 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r678 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a676, %b677, %acc675 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a679 = vector.transfer_read %sv664[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b680 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r681 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a679, %b680, %r678 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g682 = vector.transfer_read %el6[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc683 = arith.mulf %r681, %g682 : vector<8x8xf64>
      vector.transfer_write %sc683, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d684 = vector.transfer_read %sv665[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2685 = vector.transfer_read %el8[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t686 = arith.mulf %d684, %g2685 : vector<8x8xf64>
      %a687 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2688 = arith.addf %t686, %a687 : vector<8x8xf64>
      %b689 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl690 = arith.addf %t2688, %b689 : vector<8x8xf64>
      vector.transfer_write %fl690, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc691 = arith.constant dense<0.0> : vector<8x8xf64>
      %a692 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b693 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r694 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a692, %b693, %acc691 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a695 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b696 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r697 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a695, %b696, %r694 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r697, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc698 = arith.constant dense<0.0> : vector<8x8xf64>
      %a699 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b700 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r701 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a699, %b700, %acc698 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a702 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b703 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r704 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a702, %b703, %r701 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r704, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv705 = memref.subview %el2[0, 1, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv706 = memref.subview %el2[0, 2, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa707 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r708 = vector.transfer_read %sv705[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %m709 = arith.subf %r708, %fa707 : vector<8x8xf64>
      vector.transfer_write %m709, %sv705[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
      %fa710 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r711 = vector.transfer_read %sv706[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %p712 = arith.addf %r711, %fa710 : vector<8x8xf64>
      vector.transfer_write %p712, %sv706[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    }
    gpu.barrier
    %sv713 = memref.subview %v3613[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv714 = memref.subview %v3614[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc715 = arith.constant dense<0.0> : vector<8x8xf64>
      %a716 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b717 = vector.transfer_read %sv713[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r718 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a716, %b717, %acc715 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a719 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b720 = vector.transfer_read %sv713[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r721 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a719, %b720, %r718 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g722 = vector.transfer_read %el7[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc723 = arith.mulf %r721, %g722 : vector<8x8xf64>
      vector.transfer_write %sc723, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc724 = arith.constant dense<0.0> : vector<8x8xf64>
      %a725 = vector.transfer_read %sv713[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b726 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r727 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a725, %b726, %acc724 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a728 = vector.transfer_read %sv713[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b729 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r730 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a728, %b729, %r727 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g731 = vector.transfer_read %el6[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc732 = arith.mulf %r730, %g731 : vector<8x8xf64>
      vector.transfer_write %sc732, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d733 = vector.transfer_read %sv714[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2734 = vector.transfer_read %el8[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t735 = arith.mulf %d733, %g2734 : vector<8x8xf64>
      %a736 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2737 = arith.addf %t735, %a736 : vector<8x8xf64>
      %b738 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl739 = arith.addf %t2737, %b738 : vector<8x8xf64>
      vector.transfer_write %fl739, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc740 = arith.constant dense<0.0> : vector<8x8xf64>
      %a741 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b742 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r743 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a741, %b742, %acc740 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a744 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b745 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r746 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a744, %b745, %r743 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r746, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc747 = arith.constant dense<0.0> : vector<8x8xf64>
      %a748 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b749 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r750 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a748, %b749, %acc747 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a751 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b752 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r753 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a751, %b752, %r750 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r753, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv754 = memref.subview %el2[0, 2, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv755 = memref.subview %el2[0, 3, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa756 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r757 = vector.transfer_read %sv754[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %m758 = arith.subf %r757, %fa756 : vector<8x8xf64>
      vector.transfer_write %m758, %sv754[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
      %fa759 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r760 = vector.transfer_read %sv755[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %p761 = arith.addf %r760, %fa759 : vector<8x8xf64>
      vector.transfer_write %p761, %sv755[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    }
    gpu.barrier
    %sv762 = memref.subview %v3613[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv763 = memref.subview %v3614[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc764 = arith.constant dense<0.0> : vector<8x8xf64>
      %a765 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b766 = vector.transfer_read %sv762[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r767 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a765, %b766, %acc764 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a768 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b769 = vector.transfer_read %sv762[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r770 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a768, %b769, %r767 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g771 = vector.transfer_read %el7[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc772 = arith.mulf %r770, %g771 : vector<8x8xf64>
      vector.transfer_write %sc772, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc773 = arith.constant dense<0.0> : vector<8x8xf64>
      %a774 = vector.transfer_read %sv762[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b775 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r776 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a774, %b775, %acc773 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a777 = vector.transfer_read %sv762[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b778 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r779 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a777, %b778, %r776 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g780 = vector.transfer_read %el6[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc781 = arith.mulf %r779, %g780 : vector<8x8xf64>
      vector.transfer_write %sc781, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d782 = vector.transfer_read %sv763[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2783 = vector.transfer_read %el8[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t784 = arith.mulf %d782, %g2783 : vector<8x8xf64>
      %a785 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2786 = arith.addf %t784, %a785 : vector<8x8xf64>
      %b787 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl788 = arith.addf %t2786, %b787 : vector<8x8xf64>
      vector.transfer_write %fl788, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc789 = arith.constant dense<0.0> : vector<8x8xf64>
      %a790 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b791 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r792 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a790, %b791, %acc789 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a793 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b794 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r795 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a793, %b794, %r792 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r795, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc796 = arith.constant dense<0.0> : vector<8x8xf64>
      %a797 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b798 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r799 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a797, %b798, %acc796 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a800 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b801 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r802 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a800, %b801, %r799 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r802, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv803 = memref.subview %el2[0, 3, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv804 = memref.subview %el2[0, 4, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa805 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r806 = vector.transfer_read %sv803[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %m807 = arith.subf %r806, %fa805 : vector<8x8xf64>
      vector.transfer_write %m807, %sv803[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
      %fa808 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r809 = vector.transfer_read %sv804[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %p810 = arith.addf %r809, %fa808 : vector<8x8xf64>
      vector.transfer_write %p810, %sv804[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    }
    gpu.barrier
    %sv811 = memref.subview %v3613[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv812 = memref.subview %v3614[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc813 = arith.constant dense<0.0> : vector<8x8xf64>
      %a814 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b815 = vector.transfer_read %sv811[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r816 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a814, %b815, %acc813 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a817 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b818 = vector.transfer_read %sv811[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r819 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a817, %b818, %r816 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g820 = vector.transfer_read %el7[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc821 = arith.mulf %r819, %g820 : vector<8x8xf64>
      vector.transfer_write %sc821, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc822 = arith.constant dense<0.0> : vector<8x8xf64>
      %a823 = vector.transfer_read %sv811[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b824 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r825 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a823, %b824, %acc822 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a826 = vector.transfer_read %sv811[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b827 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r828 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a826, %b827, %r825 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g829 = vector.transfer_read %el6[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc830 = arith.mulf %r828, %g829 : vector<8x8xf64>
      vector.transfer_write %sc830, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d831 = vector.transfer_read %sv812[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2832 = vector.transfer_read %el8[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t833 = arith.mulf %d831, %g2832 : vector<8x8xf64>
      %a834 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2835 = arith.addf %t833, %a834 : vector<8x8xf64>
      %b836 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl837 = arith.addf %t2835, %b836 : vector<8x8xf64>
      vector.transfer_write %fl837, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc838 = arith.constant dense<0.0> : vector<8x8xf64>
      %a839 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b840 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r841 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a839, %b840, %acc838 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a842 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b843 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r844 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a842, %b843, %r841 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r844, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc845 = arith.constant dense<0.0> : vector<8x8xf64>
      %a846 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b847 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r848 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a846, %b847, %acc845 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a849 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b850 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r851 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a849, %b850, %r848 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r851, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv852 = memref.subview %el2[0, 4, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv853 = memref.subview %el2[0, 5, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa854 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r855 = vector.transfer_read %sv852[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %m856 = arith.subf %r855, %fa854 : vector<8x8xf64>
      vector.transfer_write %m856, %sv852[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
      %fa857 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r858 = vector.transfer_read %sv853[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %p859 = arith.addf %r858, %fa857 : vector<8x8xf64>
      vector.transfer_write %p859, %sv853[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    }
    gpu.barrier
    %sv860 = memref.subview %v3613[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv861 = memref.subview %v3614[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc862 = arith.constant dense<0.0> : vector<8x8xf64>
      %a863 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b864 = vector.transfer_read %sv860[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r865 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a863, %b864, %acc862 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a866 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b867 = vector.transfer_read %sv860[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r868 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a866, %b867, %r865 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g869 = vector.transfer_read %el7[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc870 = arith.mulf %r868, %g869 : vector<8x8xf64>
      vector.transfer_write %sc870, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc871 = arith.constant dense<0.0> : vector<8x8xf64>
      %a872 = vector.transfer_read %sv860[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b873 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r874 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a872, %b873, %acc871 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a875 = vector.transfer_read %sv860[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b876 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r877 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a875, %b876, %r874 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g878 = vector.transfer_read %el6[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc879 = arith.mulf %r877, %g878 : vector<8x8xf64>
      vector.transfer_write %sc879, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d880 = vector.transfer_read %sv861[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2881 = vector.transfer_read %el8[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t882 = arith.mulf %d880, %g2881 : vector<8x8xf64>
      %a883 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2884 = arith.addf %t882, %a883 : vector<8x8xf64>
      %b885 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl886 = arith.addf %t2884, %b885 : vector<8x8xf64>
      vector.transfer_write %fl886, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc887 = arith.constant dense<0.0> : vector<8x8xf64>
      %a888 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b889 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r890 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a888, %b889, %acc887 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a891 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b892 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r893 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a891, %b892, %r890 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r893, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc894 = arith.constant dense<0.0> : vector<8x8xf64>
      %a895 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b896 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r897 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a895, %b896, %acc894 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a898 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b899 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r900 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a898, %b899, %r897 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r900, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv901 = memref.subview %el2[0, 5, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv902 = memref.subview %el2[0, 6, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa903 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r904 = vector.transfer_read %sv901[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %m905 = arith.subf %r904, %fa903 : vector<8x8xf64>
      vector.transfer_write %m905, %sv901[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
      %fa906 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r907 = vector.transfer_read %sv902[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %p908 = arith.addf %r907, %fa906 : vector<8x8xf64>
      vector.transfer_write %p908, %sv902[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    }
    gpu.barrier
    %sv909 = memref.subview %v3613[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv910 = memref.subview %v3614[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc911 = arith.constant dense<0.0> : vector<8x8xf64>
      %a912 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b913 = vector.transfer_read %sv909[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r914 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a912, %b913, %acc911 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a915 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b916 = vector.transfer_read %sv909[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r917 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a915, %b916, %r914 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g918 = vector.transfer_read %el7[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc919 = arith.mulf %r917, %g918 : vector<8x8xf64>
      vector.transfer_write %sc919, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc920 = arith.constant dense<0.0> : vector<8x8xf64>
      %a921 = vector.transfer_read %sv909[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b922 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r923 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a921, %b922, %acc920 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a924 = vector.transfer_read %sv909[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b925 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r926 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a924, %b925, %r923 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g927 = vector.transfer_read %el6[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc928 = arith.mulf %r926, %g927 : vector<8x8xf64>
      vector.transfer_write %sc928, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d929 = vector.transfer_read %sv910[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2930 = vector.transfer_read %el8[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t931 = arith.mulf %d929, %g2930 : vector<8x8xf64>
      %a932 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2933 = arith.addf %t931, %a932 : vector<8x8xf64>
      %b934 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl935 = arith.addf %t2933, %b934 : vector<8x8xf64>
      vector.transfer_write %fl935, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc936 = arith.constant dense<0.0> : vector<8x8xf64>
      %a937 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b938 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r939 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a937, %b938, %acc936 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a940 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b941 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r942 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a940, %b941, %r939 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r942, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc943 = arith.constant dense<0.0> : vector<8x8xf64>
      %a944 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b945 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r946 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a944, %b945, %acc943 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a947 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b948 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r949 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a947, %b948, %r946 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r949, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv950 = memref.subview %el2[0, 6, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv951 = memref.subview %el2[0, 7, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa952 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r953 = vector.transfer_read %sv950[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %m954 = arith.subf %r953, %fa952 : vector<8x8xf64>
      vector.transfer_write %m954, %sv950[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
      %fa955 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r956 = vector.transfer_read %sv951[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %p957 = arith.addf %r956, %fa955 : vector<8x8xf64>
      vector.transfer_write %p957, %sv951[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
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
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1090 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1091 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1092 = vector.transfer_read %sv1088[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1093 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1091, %b1092, %acc1090 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1094 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1095 = vector.transfer_read %sv1088[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1096 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1094, %b1095, %r1093 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1097 = vector.transfer_read %el10[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc1098 = arith.mulf %r1096, %g1097 : vector<8x8xf64>
      vector.transfer_write %sc1098, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1099 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1100 = vector.transfer_read %sv1088[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1101 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1102 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1100, %b1101, %acc1099 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1103 = vector.transfer_read %sv1088[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1104 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1105 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1103, %b1104, %r1102 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1106 = vector.transfer_read %el9[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc1107 = arith.mulf %r1105, %g1106 : vector<8x8xf64>
      vector.transfer_write %sc1107, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d1108 = vector.transfer_read %sv1089[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g21109 = vector.transfer_read %el11[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t1110 = arith.mulf %d1108, %g21109 : vector<8x8xf64>
      %a1111 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t21112 = arith.addf %t1110, %a1111 : vector<8x8xf64>
      %b1113 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl1114 = arith.addf %t21112, %b1113 : vector<8x8xf64>
      vector.transfer_write %fl1114, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1115 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1116 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1117 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1118 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1116, %b1117, %acc1115 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1119 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1120 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1121 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1119, %b1120, %r1118 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1121, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1122 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1123 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1124 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1125 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1123, %b1124, %acc1122 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1126 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1127 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1128 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1126, %b1127, %r1125 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1128, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1129 = memref.subview %el2[0, 0, 0] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv1130 = memref.subview %el2[0, 0, 1] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1131 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1132 = vector.transfer_read %sv1129[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %m1133 = arith.subf %r1132, %fa1131 : vector<8x8xf64>
      vector.transfer_write %m1133, %sv1129[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
      %fa1134 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1135 = vector.transfer_read %sv1130[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %p1136 = arith.addf %r1135, %fa1134 : vector<8x8xf64>
      vector.transfer_write %p1136, %sv1130[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    }
    gpu.barrier
    %sv1137 = memref.subview %v31086[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv1138 = memref.subview %v31087[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1139 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1140 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1141 = vector.transfer_read %sv1137[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1142 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1140, %b1141, %acc1139 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1143 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1144 = vector.transfer_read %sv1137[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1145 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1143, %b1144, %r1142 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1146 = vector.transfer_read %el10[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc1147 = arith.mulf %r1145, %g1146 : vector<8x8xf64>
      vector.transfer_write %sc1147, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1148 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1149 = vector.transfer_read %sv1137[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1150 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1151 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1149, %b1150, %acc1148 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1152 = vector.transfer_read %sv1137[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1153 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1154 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1152, %b1153, %r1151 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1155 = vector.transfer_read %el9[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc1156 = arith.mulf %r1154, %g1155 : vector<8x8xf64>
      vector.transfer_write %sc1156, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d1157 = vector.transfer_read %sv1138[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g21158 = vector.transfer_read %el11[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t1159 = arith.mulf %d1157, %g21158 : vector<8x8xf64>
      %a1160 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t21161 = arith.addf %t1159, %a1160 : vector<8x8xf64>
      %b1162 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl1163 = arith.addf %t21161, %b1162 : vector<8x8xf64>
      vector.transfer_write %fl1163, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1164 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1165 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1166 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1167 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1165, %b1166, %acc1164 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1168 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1169 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1170 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1168, %b1169, %r1167 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1170, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1171 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1172 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1173 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1174 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1172, %b1173, %acc1171 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1175 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1176 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1177 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1175, %b1176, %r1174 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1177, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1178 = memref.subview %el2[0, 0, 1] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv1179 = memref.subview %el2[0, 0, 2] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1180 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1181 = vector.transfer_read %sv1178[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %m1182 = arith.subf %r1181, %fa1180 : vector<8x8xf64>
      vector.transfer_write %m1182, %sv1178[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
      %fa1183 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1184 = vector.transfer_read %sv1179[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %p1185 = arith.addf %r1184, %fa1183 : vector<8x8xf64>
      vector.transfer_write %p1185, %sv1179[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    }
    gpu.barrier
    %sv1186 = memref.subview %v31086[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv1187 = memref.subview %v31087[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1188 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1189 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1190 = vector.transfer_read %sv1186[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1191 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1189, %b1190, %acc1188 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1192 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1193 = vector.transfer_read %sv1186[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1194 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1192, %b1193, %r1191 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1195 = vector.transfer_read %el10[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc1196 = arith.mulf %r1194, %g1195 : vector<8x8xf64>
      vector.transfer_write %sc1196, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1197 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1198 = vector.transfer_read %sv1186[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1199 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1200 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1198, %b1199, %acc1197 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1201 = vector.transfer_read %sv1186[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1202 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1203 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1201, %b1202, %r1200 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1204 = vector.transfer_read %el9[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc1205 = arith.mulf %r1203, %g1204 : vector<8x8xf64>
      vector.transfer_write %sc1205, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d1206 = vector.transfer_read %sv1187[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g21207 = vector.transfer_read %el11[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t1208 = arith.mulf %d1206, %g21207 : vector<8x8xf64>
      %a1209 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t21210 = arith.addf %t1208, %a1209 : vector<8x8xf64>
      %b1211 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl1212 = arith.addf %t21210, %b1211 : vector<8x8xf64>
      vector.transfer_write %fl1212, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1213 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1214 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1215 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1216 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1214, %b1215, %acc1213 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1217 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1218 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1219 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1217, %b1218, %r1216 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1219, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1220 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1221 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1222 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1223 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1221, %b1222, %acc1220 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1224 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1225 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1226 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1224, %b1225, %r1223 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1226, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1227 = memref.subview %el2[0, 0, 2] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv1228 = memref.subview %el2[0, 0, 3] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1229 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1230 = vector.transfer_read %sv1227[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %m1231 = arith.subf %r1230, %fa1229 : vector<8x8xf64>
      vector.transfer_write %m1231, %sv1227[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
      %fa1232 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1233 = vector.transfer_read %sv1228[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %p1234 = arith.addf %r1233, %fa1232 : vector<8x8xf64>
      vector.transfer_write %p1234, %sv1228[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    }
    gpu.barrier
    %sv1235 = memref.subview %v31086[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv1236 = memref.subview %v31087[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1237 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1238 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1239 = vector.transfer_read %sv1235[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1240 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1238, %b1239, %acc1237 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1241 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1242 = vector.transfer_read %sv1235[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1243 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1241, %b1242, %r1240 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1244 = vector.transfer_read %el10[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc1245 = arith.mulf %r1243, %g1244 : vector<8x8xf64>
      vector.transfer_write %sc1245, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1246 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1247 = vector.transfer_read %sv1235[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1248 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1249 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1247, %b1248, %acc1246 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1250 = vector.transfer_read %sv1235[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1251 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1252 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1250, %b1251, %r1249 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1253 = vector.transfer_read %el9[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc1254 = arith.mulf %r1252, %g1253 : vector<8x8xf64>
      vector.transfer_write %sc1254, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d1255 = vector.transfer_read %sv1236[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g21256 = vector.transfer_read %el11[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t1257 = arith.mulf %d1255, %g21256 : vector<8x8xf64>
      %a1258 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t21259 = arith.addf %t1257, %a1258 : vector<8x8xf64>
      %b1260 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl1261 = arith.addf %t21259, %b1260 : vector<8x8xf64>
      vector.transfer_write %fl1261, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1262 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1263 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1264 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1265 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1263, %b1264, %acc1262 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1266 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1267 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1268 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1266, %b1267, %r1265 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1268, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1269 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1270 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1271 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1272 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1270, %b1271, %acc1269 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1273 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1274 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1275 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1273, %b1274, %r1272 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1275, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1276 = memref.subview %el2[0, 0, 3] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv1277 = memref.subview %el2[0, 0, 4] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1278 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1279 = vector.transfer_read %sv1276[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %m1280 = arith.subf %r1279, %fa1278 : vector<8x8xf64>
      vector.transfer_write %m1280, %sv1276[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
      %fa1281 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1282 = vector.transfer_read %sv1277[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %p1283 = arith.addf %r1282, %fa1281 : vector<8x8xf64>
      vector.transfer_write %p1283, %sv1277[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    }
    gpu.barrier
    %sv1284 = memref.subview %v31086[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv1285 = memref.subview %v31087[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1286 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1287 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1288 = vector.transfer_read %sv1284[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1289 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1287, %b1288, %acc1286 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1290 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1291 = vector.transfer_read %sv1284[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1292 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1290, %b1291, %r1289 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1293 = vector.transfer_read %el10[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc1294 = arith.mulf %r1292, %g1293 : vector<8x8xf64>
      vector.transfer_write %sc1294, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1295 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1296 = vector.transfer_read %sv1284[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1297 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1298 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1296, %b1297, %acc1295 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1299 = vector.transfer_read %sv1284[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1300 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1301 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1299, %b1300, %r1298 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1302 = vector.transfer_read %el9[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc1303 = arith.mulf %r1301, %g1302 : vector<8x8xf64>
      vector.transfer_write %sc1303, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d1304 = vector.transfer_read %sv1285[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g21305 = vector.transfer_read %el11[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t1306 = arith.mulf %d1304, %g21305 : vector<8x8xf64>
      %a1307 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t21308 = arith.addf %t1306, %a1307 : vector<8x8xf64>
      %b1309 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl1310 = arith.addf %t21308, %b1309 : vector<8x8xf64>
      vector.transfer_write %fl1310, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1311 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1312 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1313 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1314 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1312, %b1313, %acc1311 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1315 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1316 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1317 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1315, %b1316, %r1314 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1317, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1318 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1319 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1320 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1321 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1319, %b1320, %acc1318 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1322 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1323 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1324 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1322, %b1323, %r1321 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1324, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1325 = memref.subview %el2[0, 0, 4] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv1326 = memref.subview %el2[0, 0, 5] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1327 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1328 = vector.transfer_read %sv1325[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %m1329 = arith.subf %r1328, %fa1327 : vector<8x8xf64>
      vector.transfer_write %m1329, %sv1325[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
      %fa1330 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1331 = vector.transfer_read %sv1326[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %p1332 = arith.addf %r1331, %fa1330 : vector<8x8xf64>
      vector.transfer_write %p1332, %sv1326[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    }
    gpu.barrier
    %sv1333 = memref.subview %v31086[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv1334 = memref.subview %v31087[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1335 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1336 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1337 = vector.transfer_read %sv1333[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1338 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1336, %b1337, %acc1335 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1339 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1340 = vector.transfer_read %sv1333[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1341 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1339, %b1340, %r1338 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1342 = vector.transfer_read %el10[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc1343 = arith.mulf %r1341, %g1342 : vector<8x8xf64>
      vector.transfer_write %sc1343, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1344 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1345 = vector.transfer_read %sv1333[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1346 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1347 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1345, %b1346, %acc1344 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1348 = vector.transfer_read %sv1333[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1349 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1350 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1348, %b1349, %r1347 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1351 = vector.transfer_read %el9[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc1352 = arith.mulf %r1350, %g1351 : vector<8x8xf64>
      vector.transfer_write %sc1352, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d1353 = vector.transfer_read %sv1334[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g21354 = vector.transfer_read %el11[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t1355 = arith.mulf %d1353, %g21354 : vector<8x8xf64>
      %a1356 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t21357 = arith.addf %t1355, %a1356 : vector<8x8xf64>
      %b1358 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl1359 = arith.addf %t21357, %b1358 : vector<8x8xf64>
      vector.transfer_write %fl1359, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1360 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1361 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1362 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1363 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1361, %b1362, %acc1360 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1364 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1365 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1366 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1364, %b1365, %r1363 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1366, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1367 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1368 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1369 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1370 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1368, %b1369, %acc1367 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1371 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1372 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1373 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1371, %b1372, %r1370 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1373, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1374 = memref.subview %el2[0, 0, 5] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv1375 = memref.subview %el2[0, 0, 6] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1376 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1377 = vector.transfer_read %sv1374[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %m1378 = arith.subf %r1377, %fa1376 : vector<8x8xf64>
      vector.transfer_write %m1378, %sv1374[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
      %fa1379 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1380 = vector.transfer_read %sv1375[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %p1381 = arith.addf %r1380, %fa1379 : vector<8x8xf64>
      vector.transfer_write %p1381, %sv1375[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    }
    gpu.barrier
    %sv1382 = memref.subview %v31086[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv1383 = memref.subview %v31087[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1384 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1385 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1386 = vector.transfer_read %sv1382[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1387 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1385, %b1386, %acc1384 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1388 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1389 = vector.transfer_read %sv1382[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1390 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1388, %b1389, %r1387 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1391 = vector.transfer_read %el10[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc1392 = arith.mulf %r1390, %g1391 : vector<8x8xf64>
      vector.transfer_write %sc1392, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1393 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1394 = vector.transfer_read %sv1382[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1395 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1396 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1394, %b1395, %acc1393 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1397 = vector.transfer_read %sv1382[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1398 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1399 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1397, %b1398, %r1396 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1400 = vector.transfer_read %el9[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc1401 = arith.mulf %r1399, %g1400 : vector<8x8xf64>
      vector.transfer_write %sc1401, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d1402 = vector.transfer_read %sv1383[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g21403 = vector.transfer_read %el11[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t1404 = arith.mulf %d1402, %g21403 : vector<8x8xf64>
      %a1405 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t21406 = arith.addf %t1404, %a1405 : vector<8x8xf64>
      %b1407 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl1408 = arith.addf %t21406, %b1407 : vector<8x8xf64>
      vector.transfer_write %fl1408, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1409 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1410 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1411 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1412 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1410, %b1411, %acc1409 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1413 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1414 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1415 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1413, %b1414, %r1412 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1415, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1416 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1417 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1418 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1419 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1417, %b1418, %acc1416 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1420 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1421 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1422 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1420, %b1421, %r1419 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1422, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1423 = memref.subview %el2[0, 0, 6] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv1424 = memref.subview %el2[0, 0, 7] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1425 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1426 = vector.transfer_read %sv1423[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %m1427 = arith.subf %r1426, %fa1425 : vector<8x8xf64>
      vector.transfer_write %m1427, %sv1423[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
      %fa1428 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1429 = vector.transfer_read %sv1424[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %p1430 = arith.addf %r1429, %fa1428 : vector<8x8xf64>
      vector.transfer_write %p1430, %sv1424[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    }
    gpu.barrier
    gpu.return
  }
}
