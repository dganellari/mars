#a = affine_map<(m, n, k) -> (m, k)>
#b = affine_map<(m, n, k) -> (k, n)>
#bt = affine_map<(m, n, k) -> (n, k)>
#c = affine_map<(m, n, k) -> (m, n)>
gpu.module @mir_kernels [#nvvm.target<chip = "sm_90", O = 3>] {
  gpu.func @full_hybrid(%U: memref<?x8x8x8xf64>, %Btil: memref<8x8xf64>, %Dtil: memref<8x8xf64>, %Dm: memref<8x8xf64>, %W: memref<8x8xf64>, %g00all: memref<?x8x8xf64>, %g01all: memref<?x8x8xf64>, %g02all: memref<?x8x8xf64>, %g10all: memref<?x8x8xf64>, %g11all: memref<?x8x8xf64>, %g12all: memref<?x8x8xf64>, %g20all: memref<?x8x8xf64>, %g21all: memref<?x8x8xf64>, %g22all: memref<?x8x8xf64>, %Y: memref<?x8x8x8xf64>) workgroup(%interp_all: memref<8x64xf64, #gpu.address_space<workgroup>>, %deriv_all: memref<8x64xf64, #gpu.address_space<workgroup>>) kernel {
    %z = arith.constant 0.0 : f64
    %zc1 = arith.constant dense<0.0> : vector<1x1xf64>
    %zc2 = arith.constant dense<0.0> : vector<1x2xf64>
    %lane = gpu.thread_id x
    %e = gpu.block_id x
    %c2 = arith.constant 2 : index
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
    %i = arith.divui %lane, %c4 : index
    %k = arith.remui %lane, %c4 : index
    %k4 = arith.addi %k, %c4 : index
    %tc = arith.muli %k, %c2 : index
    %li = arith.index_cast %lane : index to i32
    %j1 = arith.constant 1 : i32
    %j2 = arith.constant 2 : i32
    %j3 = arith.constant 3 : i32
    %j4 = arith.constant 4 : i32
    %j16 = arith.constant 16 : i32
    %j28 = arith.constant 28 : i32
    %j32 = arith.constant 32 : i32
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
    %f144 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f145 = vector.transfer_read %sv142[%k, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m146 = nvgpu.mma.sync(%f144, %f145, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f147 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f148 = vector.transfer_read %sv142[%k4, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m149 = nvgpu.mma.sync(%f147, %f148, %m146) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f150 = vector.transfer_read %sv142[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f151 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m152 = nvgpu.mma.sync(%f150, %f151, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f153 = vector.transfer_read %sv142[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f154 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m155 = nvgpu.mma.sync(%f153, %f154, %m152) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %cf156 = vector.transfer_read %sv143[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<1x2xf64>
    %cf157 = vector.transfer_read %el5[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf158 = vector.transfer_read %el3[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf159 = vector.transfer_read %el4[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %x160 = arith.mulf %cf157, %cf156 : vector<1x2xf64>
    %x161 = arith.mulf %cf158, %m155 : vector<1x2xf64>
    %x162 = arith.mulf %cf159, %m149 : vector<1x2xf64>
    %x163 = arith.addf %x160, %x161 : vector<1x2xf64>
    %fl164 = arith.addf %x163, %x162 : vector<1x2xf64>
    %lo165 = vector.extract %fl164[0, 0] : f64 from vector<1x2xf64>
    %hi166 = vector.extract %fl164[0, 1] : f64 from vector<1x2xf64>
    %t168 = arith.andi %li, %j28 : i32
    %t169 = arith.andi %li, %j3 : i32
    %t170 = arith.shrui %t169, %j1 : i32
    %t171 = arith.ori %t168, %t170 : i32
    %s172, %unused172 = gpu.shuffle idx %lo165, %t171, %j32 : f64
    %s174, %unused174 = gpu.shuffle idx %hi166, %t171, %j32 : f64
    %t176 = arith.andi %li, %j1 : i32
    %t177 = arith.cmpi eq, %t176, %j1 : i32
    %sv178 = arith.select %t177, %s174, %s172 : f64
    %a179 = vector.insert %sv178, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f180 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m181 = nvgpu.mma.sync(%a179, %f180, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %lo182 = vector.extract %fl164[0, 0] : f64 from vector<1x2xf64>
    %hi183 = vector.extract %fl164[0, 1] : f64 from vector<1x2xf64>
    %t185 = arith.andi %li, %j28 : i32
    %t186 = arith.andi %li, %j3 : i32
    %t187 = arith.shrui %t186, %j1 : i32
    %t188 = arith.ori %t185, %t187 : i32
    %t189 = arith.ori %t188, %j2 : i32
    %s190, %unused190 = gpu.shuffle idx %lo182, %t189, %j32 : f64
    %s192, %unused192 = gpu.shuffle idx %hi183, %t189, %j32 : f64
    %t194 = arith.andi %li, %j1 : i32
    %t195 = arith.cmpi eq, %t194, %j1 : i32
    %sv196 = arith.select %t195, %s192, %s190 : f64
    %a197 = vector.insert %sv196, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f198 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m199 = nvgpu.mma.sync(%a197, %f198, %m181) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f200 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo201 = vector.extract %m199[0, 0] : f64 from vector<1x2xf64>
    %hi202 = vector.extract %m199[0, 1] : f64 from vector<1x2xf64>
    %t203 = arith.andi %li, %j3 : i32
    %t204 = arith.muli %t203, %j4 : i32
    %t205 = arith.shrui %li, %j3 : i32
    %t206 = arith.addi %t204, %t205 : i32
    %s207, %unused207 = gpu.shuffle idx %lo201, %t206, %j32 : f64
    %s209, %unused209 = gpu.shuffle idx %hi202, %t206, %j32 : f64
    %t211 = arith.shrui %li, %j2 : i32
    %t212 = arith.andi %t211, %j1 : i32
    %t213 = arith.cmpi eq, %t212, %j1 : i32
    %sv214 = arith.select %t213, %s209, %s207 : f64
    %b215 = vector.insert %sv214, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m216 = nvgpu.mma.sync(%f200, %b215, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f217 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo218 = vector.extract %m199[0, 0] : f64 from vector<1x2xf64>
    %hi219 = vector.extract %m199[0, 1] : f64 from vector<1x2xf64>
    %t220 = arith.andi %li, %j3 : i32
    %t221 = arith.muli %t220, %j4 : i32
    %t222 = arith.shrui %li, %j3 : i32
    %t223 = arith.addi %t221, %t222 : i32
    %t224 = arith.addi %t223, %j16 : i32
    %s225, %unused225 = gpu.shuffle idx %lo218, %t224, %j32 : f64
    %s227, %unused227 = gpu.shuffle idx %hi219, %t224, %j32 : f64
    %t229 = arith.shrui %li, %j2 : i32
    %t230 = arith.andi %t229, %j1 : i32
    %t231 = arith.cmpi eq, %t230, %j1 : i32
    %sv232 = arith.select %t231, %s227, %s225 : f64
    %b233 = vector.insert %sv232, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m234 = nvgpu.mma.sync(%f217, %b233, %m216) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %sv235 = memref.subview %el2[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv236 = memref.subview %el2[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %y237 = vector.transfer_read %sv235[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %y238 = arith.subf %y237, %m234 : vector<1x2xf64>
    vector.transfer_write %y238, %sv235[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    %y239 = vector.transfer_read %sv236[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %y240 = arith.addf %y239, %m234 : vector<1x2xf64>
    vector.transfer_write %y240, %sv236[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv241 = memref.subview %v3140[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv242 = memref.subview %v3141[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %f243 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f244 = vector.transfer_read %sv241[%k, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m245 = nvgpu.mma.sync(%f243, %f244, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f246 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f247 = vector.transfer_read %sv241[%k4, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m248 = nvgpu.mma.sync(%f246, %f247, %m245) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f249 = vector.transfer_read %sv241[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f250 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m251 = nvgpu.mma.sync(%f249, %f250, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f252 = vector.transfer_read %sv241[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f253 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m254 = nvgpu.mma.sync(%f252, %f253, %m251) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %cf255 = vector.transfer_read %sv242[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<1x2xf64>
    %cf256 = vector.transfer_read %el5[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf257 = vector.transfer_read %el3[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf258 = vector.transfer_read %el4[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %x259 = arith.mulf %cf256, %cf255 : vector<1x2xf64>
    %x260 = arith.mulf %cf257, %m254 : vector<1x2xf64>
    %x261 = arith.mulf %cf258, %m248 : vector<1x2xf64>
    %x262 = arith.addf %x259, %x260 : vector<1x2xf64>
    %fl263 = arith.addf %x262, %x261 : vector<1x2xf64>
    %lo264 = vector.extract %fl263[0, 0] : f64 from vector<1x2xf64>
    %hi265 = vector.extract %fl263[0, 1] : f64 from vector<1x2xf64>
    %t267 = arith.andi %li, %j28 : i32
    %t268 = arith.andi %li, %j3 : i32
    %t269 = arith.shrui %t268, %j1 : i32
    %t270 = arith.ori %t267, %t269 : i32
    %s271, %unused271 = gpu.shuffle idx %lo264, %t270, %j32 : f64
    %s273, %unused273 = gpu.shuffle idx %hi265, %t270, %j32 : f64
    %t275 = arith.andi %li, %j1 : i32
    %t276 = arith.cmpi eq, %t275, %j1 : i32
    %sv277 = arith.select %t276, %s273, %s271 : f64
    %a278 = vector.insert %sv277, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f279 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m280 = nvgpu.mma.sync(%a278, %f279, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %lo281 = vector.extract %fl263[0, 0] : f64 from vector<1x2xf64>
    %hi282 = vector.extract %fl263[0, 1] : f64 from vector<1x2xf64>
    %t284 = arith.andi %li, %j28 : i32
    %t285 = arith.andi %li, %j3 : i32
    %t286 = arith.shrui %t285, %j1 : i32
    %t287 = arith.ori %t284, %t286 : i32
    %t288 = arith.ori %t287, %j2 : i32
    %s289, %unused289 = gpu.shuffle idx %lo281, %t288, %j32 : f64
    %s291, %unused291 = gpu.shuffle idx %hi282, %t288, %j32 : f64
    %t293 = arith.andi %li, %j1 : i32
    %t294 = arith.cmpi eq, %t293, %j1 : i32
    %sv295 = arith.select %t294, %s291, %s289 : f64
    %a296 = vector.insert %sv295, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f297 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m298 = nvgpu.mma.sync(%a296, %f297, %m280) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f299 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo300 = vector.extract %m298[0, 0] : f64 from vector<1x2xf64>
    %hi301 = vector.extract %m298[0, 1] : f64 from vector<1x2xf64>
    %t302 = arith.andi %li, %j3 : i32
    %t303 = arith.muli %t302, %j4 : i32
    %t304 = arith.shrui %li, %j3 : i32
    %t305 = arith.addi %t303, %t304 : i32
    %s306, %unused306 = gpu.shuffle idx %lo300, %t305, %j32 : f64
    %s308, %unused308 = gpu.shuffle idx %hi301, %t305, %j32 : f64
    %t310 = arith.shrui %li, %j2 : i32
    %t311 = arith.andi %t310, %j1 : i32
    %t312 = arith.cmpi eq, %t311, %j1 : i32
    %sv313 = arith.select %t312, %s308, %s306 : f64
    %b314 = vector.insert %sv313, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m315 = nvgpu.mma.sync(%f299, %b314, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f316 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo317 = vector.extract %m298[0, 0] : f64 from vector<1x2xf64>
    %hi318 = vector.extract %m298[0, 1] : f64 from vector<1x2xf64>
    %t319 = arith.andi %li, %j3 : i32
    %t320 = arith.muli %t319, %j4 : i32
    %t321 = arith.shrui %li, %j3 : i32
    %t322 = arith.addi %t320, %t321 : i32
    %t323 = arith.addi %t322, %j16 : i32
    %s324, %unused324 = gpu.shuffle idx %lo317, %t323, %j32 : f64
    %s326, %unused326 = gpu.shuffle idx %hi318, %t323, %j32 : f64
    %t328 = arith.shrui %li, %j2 : i32
    %t329 = arith.andi %t328, %j1 : i32
    %t330 = arith.cmpi eq, %t329, %j1 : i32
    %sv331 = arith.select %t330, %s326, %s324 : f64
    %b332 = vector.insert %sv331, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m333 = nvgpu.mma.sync(%f316, %b332, %m315) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %sv334 = memref.subview %el2[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv335 = memref.subview %el2[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %y336 = vector.transfer_read %sv334[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %y337 = arith.subf %y336, %m333 : vector<1x2xf64>
    vector.transfer_write %y337, %sv334[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    %y338 = vector.transfer_read %sv335[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %y339 = arith.addf %y338, %m333 : vector<1x2xf64>
    vector.transfer_write %y339, %sv335[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv340 = memref.subview %v3140[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv341 = memref.subview %v3141[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %f342 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f343 = vector.transfer_read %sv340[%k, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m344 = nvgpu.mma.sync(%f342, %f343, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f345 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f346 = vector.transfer_read %sv340[%k4, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m347 = nvgpu.mma.sync(%f345, %f346, %m344) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f348 = vector.transfer_read %sv340[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f349 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m350 = nvgpu.mma.sync(%f348, %f349, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f351 = vector.transfer_read %sv340[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f352 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m353 = nvgpu.mma.sync(%f351, %f352, %m350) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %cf354 = vector.transfer_read %sv341[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<1x2xf64>
    %cf355 = vector.transfer_read %el5[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf356 = vector.transfer_read %el3[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf357 = vector.transfer_read %el4[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %x358 = arith.mulf %cf355, %cf354 : vector<1x2xf64>
    %x359 = arith.mulf %cf356, %m353 : vector<1x2xf64>
    %x360 = arith.mulf %cf357, %m347 : vector<1x2xf64>
    %x361 = arith.addf %x358, %x359 : vector<1x2xf64>
    %fl362 = arith.addf %x361, %x360 : vector<1x2xf64>
    %lo363 = vector.extract %fl362[0, 0] : f64 from vector<1x2xf64>
    %hi364 = vector.extract %fl362[0, 1] : f64 from vector<1x2xf64>
    %t366 = arith.andi %li, %j28 : i32
    %t367 = arith.andi %li, %j3 : i32
    %t368 = arith.shrui %t367, %j1 : i32
    %t369 = arith.ori %t366, %t368 : i32
    %s370, %unused370 = gpu.shuffle idx %lo363, %t369, %j32 : f64
    %s372, %unused372 = gpu.shuffle idx %hi364, %t369, %j32 : f64
    %t374 = arith.andi %li, %j1 : i32
    %t375 = arith.cmpi eq, %t374, %j1 : i32
    %sv376 = arith.select %t375, %s372, %s370 : f64
    %a377 = vector.insert %sv376, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f378 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m379 = nvgpu.mma.sync(%a377, %f378, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %lo380 = vector.extract %fl362[0, 0] : f64 from vector<1x2xf64>
    %hi381 = vector.extract %fl362[0, 1] : f64 from vector<1x2xf64>
    %t383 = arith.andi %li, %j28 : i32
    %t384 = arith.andi %li, %j3 : i32
    %t385 = arith.shrui %t384, %j1 : i32
    %t386 = arith.ori %t383, %t385 : i32
    %t387 = arith.ori %t386, %j2 : i32
    %s388, %unused388 = gpu.shuffle idx %lo380, %t387, %j32 : f64
    %s390, %unused390 = gpu.shuffle idx %hi381, %t387, %j32 : f64
    %t392 = arith.andi %li, %j1 : i32
    %t393 = arith.cmpi eq, %t392, %j1 : i32
    %sv394 = arith.select %t393, %s390, %s388 : f64
    %a395 = vector.insert %sv394, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f396 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m397 = nvgpu.mma.sync(%a395, %f396, %m379) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f398 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo399 = vector.extract %m397[0, 0] : f64 from vector<1x2xf64>
    %hi400 = vector.extract %m397[0, 1] : f64 from vector<1x2xf64>
    %t401 = arith.andi %li, %j3 : i32
    %t402 = arith.muli %t401, %j4 : i32
    %t403 = arith.shrui %li, %j3 : i32
    %t404 = arith.addi %t402, %t403 : i32
    %s405, %unused405 = gpu.shuffle idx %lo399, %t404, %j32 : f64
    %s407, %unused407 = gpu.shuffle idx %hi400, %t404, %j32 : f64
    %t409 = arith.shrui %li, %j2 : i32
    %t410 = arith.andi %t409, %j1 : i32
    %t411 = arith.cmpi eq, %t410, %j1 : i32
    %sv412 = arith.select %t411, %s407, %s405 : f64
    %b413 = vector.insert %sv412, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m414 = nvgpu.mma.sync(%f398, %b413, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f415 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo416 = vector.extract %m397[0, 0] : f64 from vector<1x2xf64>
    %hi417 = vector.extract %m397[0, 1] : f64 from vector<1x2xf64>
    %t418 = arith.andi %li, %j3 : i32
    %t419 = arith.muli %t418, %j4 : i32
    %t420 = arith.shrui %li, %j3 : i32
    %t421 = arith.addi %t419, %t420 : i32
    %t422 = arith.addi %t421, %j16 : i32
    %s423, %unused423 = gpu.shuffle idx %lo416, %t422, %j32 : f64
    %s425, %unused425 = gpu.shuffle idx %hi417, %t422, %j32 : f64
    %t427 = arith.shrui %li, %j2 : i32
    %t428 = arith.andi %t427, %j1 : i32
    %t429 = arith.cmpi eq, %t428, %j1 : i32
    %sv430 = arith.select %t429, %s425, %s423 : f64
    %b431 = vector.insert %sv430, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m432 = nvgpu.mma.sync(%f415, %b431, %m414) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %sv433 = memref.subview %el2[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv434 = memref.subview %el2[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %y435 = vector.transfer_read %sv433[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %y436 = arith.subf %y435, %m432 : vector<1x2xf64>
    vector.transfer_write %y436, %sv433[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    %y437 = vector.transfer_read %sv434[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %y438 = arith.addf %y437, %m432 : vector<1x2xf64>
    vector.transfer_write %y438, %sv434[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv439 = memref.subview %v3140[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv440 = memref.subview %v3141[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %f441 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f442 = vector.transfer_read %sv439[%k, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m443 = nvgpu.mma.sync(%f441, %f442, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f444 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f445 = vector.transfer_read %sv439[%k4, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m446 = nvgpu.mma.sync(%f444, %f445, %m443) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f447 = vector.transfer_read %sv439[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f448 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m449 = nvgpu.mma.sync(%f447, %f448, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f450 = vector.transfer_read %sv439[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f451 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m452 = nvgpu.mma.sync(%f450, %f451, %m449) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %cf453 = vector.transfer_read %sv440[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<1x2xf64>
    %cf454 = vector.transfer_read %el5[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf455 = vector.transfer_read %el3[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf456 = vector.transfer_read %el4[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %x457 = arith.mulf %cf454, %cf453 : vector<1x2xf64>
    %x458 = arith.mulf %cf455, %m452 : vector<1x2xf64>
    %x459 = arith.mulf %cf456, %m446 : vector<1x2xf64>
    %x460 = arith.addf %x457, %x458 : vector<1x2xf64>
    %fl461 = arith.addf %x460, %x459 : vector<1x2xf64>
    %lo462 = vector.extract %fl461[0, 0] : f64 from vector<1x2xf64>
    %hi463 = vector.extract %fl461[0, 1] : f64 from vector<1x2xf64>
    %t465 = arith.andi %li, %j28 : i32
    %t466 = arith.andi %li, %j3 : i32
    %t467 = arith.shrui %t466, %j1 : i32
    %t468 = arith.ori %t465, %t467 : i32
    %s469, %unused469 = gpu.shuffle idx %lo462, %t468, %j32 : f64
    %s471, %unused471 = gpu.shuffle idx %hi463, %t468, %j32 : f64
    %t473 = arith.andi %li, %j1 : i32
    %t474 = arith.cmpi eq, %t473, %j1 : i32
    %sv475 = arith.select %t474, %s471, %s469 : f64
    %a476 = vector.insert %sv475, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f477 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m478 = nvgpu.mma.sync(%a476, %f477, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %lo479 = vector.extract %fl461[0, 0] : f64 from vector<1x2xf64>
    %hi480 = vector.extract %fl461[0, 1] : f64 from vector<1x2xf64>
    %t482 = arith.andi %li, %j28 : i32
    %t483 = arith.andi %li, %j3 : i32
    %t484 = arith.shrui %t483, %j1 : i32
    %t485 = arith.ori %t482, %t484 : i32
    %t486 = arith.ori %t485, %j2 : i32
    %s487, %unused487 = gpu.shuffle idx %lo479, %t486, %j32 : f64
    %s489, %unused489 = gpu.shuffle idx %hi480, %t486, %j32 : f64
    %t491 = arith.andi %li, %j1 : i32
    %t492 = arith.cmpi eq, %t491, %j1 : i32
    %sv493 = arith.select %t492, %s489, %s487 : f64
    %a494 = vector.insert %sv493, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f495 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m496 = nvgpu.mma.sync(%a494, %f495, %m478) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f497 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo498 = vector.extract %m496[0, 0] : f64 from vector<1x2xf64>
    %hi499 = vector.extract %m496[0, 1] : f64 from vector<1x2xf64>
    %t500 = arith.andi %li, %j3 : i32
    %t501 = arith.muli %t500, %j4 : i32
    %t502 = arith.shrui %li, %j3 : i32
    %t503 = arith.addi %t501, %t502 : i32
    %s504, %unused504 = gpu.shuffle idx %lo498, %t503, %j32 : f64
    %s506, %unused506 = gpu.shuffle idx %hi499, %t503, %j32 : f64
    %t508 = arith.shrui %li, %j2 : i32
    %t509 = arith.andi %t508, %j1 : i32
    %t510 = arith.cmpi eq, %t509, %j1 : i32
    %sv511 = arith.select %t510, %s506, %s504 : f64
    %b512 = vector.insert %sv511, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m513 = nvgpu.mma.sync(%f497, %b512, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f514 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo515 = vector.extract %m496[0, 0] : f64 from vector<1x2xf64>
    %hi516 = vector.extract %m496[0, 1] : f64 from vector<1x2xf64>
    %t517 = arith.andi %li, %j3 : i32
    %t518 = arith.muli %t517, %j4 : i32
    %t519 = arith.shrui %li, %j3 : i32
    %t520 = arith.addi %t518, %t519 : i32
    %t521 = arith.addi %t520, %j16 : i32
    %s522, %unused522 = gpu.shuffle idx %lo515, %t521, %j32 : f64
    %s524, %unused524 = gpu.shuffle idx %hi516, %t521, %j32 : f64
    %t526 = arith.shrui %li, %j2 : i32
    %t527 = arith.andi %t526, %j1 : i32
    %t528 = arith.cmpi eq, %t527, %j1 : i32
    %sv529 = arith.select %t528, %s524, %s522 : f64
    %b530 = vector.insert %sv529, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m531 = nvgpu.mma.sync(%f514, %b530, %m513) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %sv532 = memref.subview %el2[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv533 = memref.subview %el2[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %y534 = vector.transfer_read %sv532[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %y535 = arith.subf %y534, %m531 : vector<1x2xf64>
    vector.transfer_write %y535, %sv532[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    %y536 = vector.transfer_read %sv533[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %y537 = arith.addf %y536, %m531 : vector<1x2xf64>
    vector.transfer_write %y537, %sv533[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv538 = memref.subview %v3140[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv539 = memref.subview %v3141[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %f540 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f541 = vector.transfer_read %sv538[%k, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m542 = nvgpu.mma.sync(%f540, %f541, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f543 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f544 = vector.transfer_read %sv538[%k4, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m545 = nvgpu.mma.sync(%f543, %f544, %m542) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f546 = vector.transfer_read %sv538[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f547 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m548 = nvgpu.mma.sync(%f546, %f547, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f549 = vector.transfer_read %sv538[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f550 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m551 = nvgpu.mma.sync(%f549, %f550, %m548) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %cf552 = vector.transfer_read %sv539[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<1x2xf64>
    %cf553 = vector.transfer_read %el5[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf554 = vector.transfer_read %el3[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf555 = vector.transfer_read %el4[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %x556 = arith.mulf %cf553, %cf552 : vector<1x2xf64>
    %x557 = arith.mulf %cf554, %m551 : vector<1x2xf64>
    %x558 = arith.mulf %cf555, %m545 : vector<1x2xf64>
    %x559 = arith.addf %x556, %x557 : vector<1x2xf64>
    %fl560 = arith.addf %x559, %x558 : vector<1x2xf64>
    %lo561 = vector.extract %fl560[0, 0] : f64 from vector<1x2xf64>
    %hi562 = vector.extract %fl560[0, 1] : f64 from vector<1x2xf64>
    %t564 = arith.andi %li, %j28 : i32
    %t565 = arith.andi %li, %j3 : i32
    %t566 = arith.shrui %t565, %j1 : i32
    %t567 = arith.ori %t564, %t566 : i32
    %s568, %unused568 = gpu.shuffle idx %lo561, %t567, %j32 : f64
    %s570, %unused570 = gpu.shuffle idx %hi562, %t567, %j32 : f64
    %t572 = arith.andi %li, %j1 : i32
    %t573 = arith.cmpi eq, %t572, %j1 : i32
    %sv574 = arith.select %t573, %s570, %s568 : f64
    %a575 = vector.insert %sv574, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f576 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m577 = nvgpu.mma.sync(%a575, %f576, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %lo578 = vector.extract %fl560[0, 0] : f64 from vector<1x2xf64>
    %hi579 = vector.extract %fl560[0, 1] : f64 from vector<1x2xf64>
    %t581 = arith.andi %li, %j28 : i32
    %t582 = arith.andi %li, %j3 : i32
    %t583 = arith.shrui %t582, %j1 : i32
    %t584 = arith.ori %t581, %t583 : i32
    %t585 = arith.ori %t584, %j2 : i32
    %s586, %unused586 = gpu.shuffle idx %lo578, %t585, %j32 : f64
    %s588, %unused588 = gpu.shuffle idx %hi579, %t585, %j32 : f64
    %t590 = arith.andi %li, %j1 : i32
    %t591 = arith.cmpi eq, %t590, %j1 : i32
    %sv592 = arith.select %t591, %s588, %s586 : f64
    %a593 = vector.insert %sv592, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f594 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m595 = nvgpu.mma.sync(%a593, %f594, %m577) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f596 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo597 = vector.extract %m595[0, 0] : f64 from vector<1x2xf64>
    %hi598 = vector.extract %m595[0, 1] : f64 from vector<1x2xf64>
    %t599 = arith.andi %li, %j3 : i32
    %t600 = arith.muli %t599, %j4 : i32
    %t601 = arith.shrui %li, %j3 : i32
    %t602 = arith.addi %t600, %t601 : i32
    %s603, %unused603 = gpu.shuffle idx %lo597, %t602, %j32 : f64
    %s605, %unused605 = gpu.shuffle idx %hi598, %t602, %j32 : f64
    %t607 = arith.shrui %li, %j2 : i32
    %t608 = arith.andi %t607, %j1 : i32
    %t609 = arith.cmpi eq, %t608, %j1 : i32
    %sv610 = arith.select %t609, %s605, %s603 : f64
    %b611 = vector.insert %sv610, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m612 = nvgpu.mma.sync(%f596, %b611, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f613 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo614 = vector.extract %m595[0, 0] : f64 from vector<1x2xf64>
    %hi615 = vector.extract %m595[0, 1] : f64 from vector<1x2xf64>
    %t616 = arith.andi %li, %j3 : i32
    %t617 = arith.muli %t616, %j4 : i32
    %t618 = arith.shrui %li, %j3 : i32
    %t619 = arith.addi %t617, %t618 : i32
    %t620 = arith.addi %t619, %j16 : i32
    %s621, %unused621 = gpu.shuffle idx %lo614, %t620, %j32 : f64
    %s623, %unused623 = gpu.shuffle idx %hi615, %t620, %j32 : f64
    %t625 = arith.shrui %li, %j2 : i32
    %t626 = arith.andi %t625, %j1 : i32
    %t627 = arith.cmpi eq, %t626, %j1 : i32
    %sv628 = arith.select %t627, %s623, %s621 : f64
    %b629 = vector.insert %sv628, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m630 = nvgpu.mma.sync(%f613, %b629, %m612) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %sv631 = memref.subview %el2[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv632 = memref.subview %el2[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %y633 = vector.transfer_read %sv631[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %y634 = arith.subf %y633, %m630 : vector<1x2xf64>
    vector.transfer_write %y634, %sv631[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    %y635 = vector.transfer_read %sv632[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %y636 = arith.addf %y635, %m630 : vector<1x2xf64>
    vector.transfer_write %y636, %sv632[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv637 = memref.subview %v3140[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv638 = memref.subview %v3141[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %f639 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f640 = vector.transfer_read %sv637[%k, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m641 = nvgpu.mma.sync(%f639, %f640, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f642 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f643 = vector.transfer_read %sv637[%k4, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m644 = nvgpu.mma.sync(%f642, %f643, %m641) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f645 = vector.transfer_read %sv637[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f646 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m647 = nvgpu.mma.sync(%f645, %f646, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f648 = vector.transfer_read %sv637[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f649 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m650 = nvgpu.mma.sync(%f648, %f649, %m647) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %cf651 = vector.transfer_read %sv638[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<1x2xf64>
    %cf652 = vector.transfer_read %el5[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf653 = vector.transfer_read %el3[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf654 = vector.transfer_read %el4[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %x655 = arith.mulf %cf652, %cf651 : vector<1x2xf64>
    %x656 = arith.mulf %cf653, %m650 : vector<1x2xf64>
    %x657 = arith.mulf %cf654, %m644 : vector<1x2xf64>
    %x658 = arith.addf %x655, %x656 : vector<1x2xf64>
    %fl659 = arith.addf %x658, %x657 : vector<1x2xf64>
    %lo660 = vector.extract %fl659[0, 0] : f64 from vector<1x2xf64>
    %hi661 = vector.extract %fl659[0, 1] : f64 from vector<1x2xf64>
    %t663 = arith.andi %li, %j28 : i32
    %t664 = arith.andi %li, %j3 : i32
    %t665 = arith.shrui %t664, %j1 : i32
    %t666 = arith.ori %t663, %t665 : i32
    %s667, %unused667 = gpu.shuffle idx %lo660, %t666, %j32 : f64
    %s669, %unused669 = gpu.shuffle idx %hi661, %t666, %j32 : f64
    %t671 = arith.andi %li, %j1 : i32
    %t672 = arith.cmpi eq, %t671, %j1 : i32
    %sv673 = arith.select %t672, %s669, %s667 : f64
    %a674 = vector.insert %sv673, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f675 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m676 = nvgpu.mma.sync(%a674, %f675, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %lo677 = vector.extract %fl659[0, 0] : f64 from vector<1x2xf64>
    %hi678 = vector.extract %fl659[0, 1] : f64 from vector<1x2xf64>
    %t680 = arith.andi %li, %j28 : i32
    %t681 = arith.andi %li, %j3 : i32
    %t682 = arith.shrui %t681, %j1 : i32
    %t683 = arith.ori %t680, %t682 : i32
    %t684 = arith.ori %t683, %j2 : i32
    %s685, %unused685 = gpu.shuffle idx %lo677, %t684, %j32 : f64
    %s687, %unused687 = gpu.shuffle idx %hi678, %t684, %j32 : f64
    %t689 = arith.andi %li, %j1 : i32
    %t690 = arith.cmpi eq, %t689, %j1 : i32
    %sv691 = arith.select %t690, %s687, %s685 : f64
    %a692 = vector.insert %sv691, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f693 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m694 = nvgpu.mma.sync(%a692, %f693, %m676) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f695 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo696 = vector.extract %m694[0, 0] : f64 from vector<1x2xf64>
    %hi697 = vector.extract %m694[0, 1] : f64 from vector<1x2xf64>
    %t698 = arith.andi %li, %j3 : i32
    %t699 = arith.muli %t698, %j4 : i32
    %t700 = arith.shrui %li, %j3 : i32
    %t701 = arith.addi %t699, %t700 : i32
    %s702, %unused702 = gpu.shuffle idx %lo696, %t701, %j32 : f64
    %s704, %unused704 = gpu.shuffle idx %hi697, %t701, %j32 : f64
    %t706 = arith.shrui %li, %j2 : i32
    %t707 = arith.andi %t706, %j1 : i32
    %t708 = arith.cmpi eq, %t707, %j1 : i32
    %sv709 = arith.select %t708, %s704, %s702 : f64
    %b710 = vector.insert %sv709, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m711 = nvgpu.mma.sync(%f695, %b710, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f712 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo713 = vector.extract %m694[0, 0] : f64 from vector<1x2xf64>
    %hi714 = vector.extract %m694[0, 1] : f64 from vector<1x2xf64>
    %t715 = arith.andi %li, %j3 : i32
    %t716 = arith.muli %t715, %j4 : i32
    %t717 = arith.shrui %li, %j3 : i32
    %t718 = arith.addi %t716, %t717 : i32
    %t719 = arith.addi %t718, %j16 : i32
    %s720, %unused720 = gpu.shuffle idx %lo713, %t719, %j32 : f64
    %s722, %unused722 = gpu.shuffle idx %hi714, %t719, %j32 : f64
    %t724 = arith.shrui %li, %j2 : i32
    %t725 = arith.andi %t724, %j1 : i32
    %t726 = arith.cmpi eq, %t725, %j1 : i32
    %sv727 = arith.select %t726, %s722, %s720 : f64
    %b728 = vector.insert %sv727, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m729 = nvgpu.mma.sync(%f712, %b728, %m711) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %sv730 = memref.subview %el2[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv731 = memref.subview %el2[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %y732 = vector.transfer_read %sv730[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %y733 = arith.subf %y732, %m729 : vector<1x2xf64>
    vector.transfer_write %y733, %sv730[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    %y734 = vector.transfer_read %sv731[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %y735 = arith.addf %y734, %m729 : vector<1x2xf64>
    vector.transfer_write %y735, %sv731[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv736 = memref.subview %v3140[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv737 = memref.subview %v3141[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %f738 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f739 = vector.transfer_read %sv736[%k, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m740 = nvgpu.mma.sync(%f738, %f739, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f741 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f742 = vector.transfer_read %sv736[%k4, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m743 = nvgpu.mma.sync(%f741, %f742, %m740) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f744 = vector.transfer_read %sv736[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f745 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m746 = nvgpu.mma.sync(%f744, %f745, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f747 = vector.transfer_read %sv736[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f748 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m749 = nvgpu.mma.sync(%f747, %f748, %m746) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %cf750 = vector.transfer_read %sv737[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<1x2xf64>
    %cf751 = vector.transfer_read %el5[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf752 = vector.transfer_read %el3[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf753 = vector.transfer_read %el4[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %x754 = arith.mulf %cf751, %cf750 : vector<1x2xf64>
    %x755 = arith.mulf %cf752, %m749 : vector<1x2xf64>
    %x756 = arith.mulf %cf753, %m743 : vector<1x2xf64>
    %x757 = arith.addf %x754, %x755 : vector<1x2xf64>
    %fl758 = arith.addf %x757, %x756 : vector<1x2xf64>
    %lo759 = vector.extract %fl758[0, 0] : f64 from vector<1x2xf64>
    %hi760 = vector.extract %fl758[0, 1] : f64 from vector<1x2xf64>
    %t762 = arith.andi %li, %j28 : i32
    %t763 = arith.andi %li, %j3 : i32
    %t764 = arith.shrui %t763, %j1 : i32
    %t765 = arith.ori %t762, %t764 : i32
    %s766, %unused766 = gpu.shuffle idx %lo759, %t765, %j32 : f64
    %s768, %unused768 = gpu.shuffle idx %hi760, %t765, %j32 : f64
    %t770 = arith.andi %li, %j1 : i32
    %t771 = arith.cmpi eq, %t770, %j1 : i32
    %sv772 = arith.select %t771, %s768, %s766 : f64
    %a773 = vector.insert %sv772, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f774 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m775 = nvgpu.mma.sync(%a773, %f774, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %lo776 = vector.extract %fl758[0, 0] : f64 from vector<1x2xf64>
    %hi777 = vector.extract %fl758[0, 1] : f64 from vector<1x2xf64>
    %t779 = arith.andi %li, %j28 : i32
    %t780 = arith.andi %li, %j3 : i32
    %t781 = arith.shrui %t780, %j1 : i32
    %t782 = arith.ori %t779, %t781 : i32
    %t783 = arith.ori %t782, %j2 : i32
    %s784, %unused784 = gpu.shuffle idx %lo776, %t783, %j32 : f64
    %s786, %unused786 = gpu.shuffle idx %hi777, %t783, %j32 : f64
    %t788 = arith.andi %li, %j1 : i32
    %t789 = arith.cmpi eq, %t788, %j1 : i32
    %sv790 = arith.select %t789, %s786, %s784 : f64
    %a791 = vector.insert %sv790, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f792 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m793 = nvgpu.mma.sync(%a791, %f792, %m775) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f794 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo795 = vector.extract %m793[0, 0] : f64 from vector<1x2xf64>
    %hi796 = vector.extract %m793[0, 1] : f64 from vector<1x2xf64>
    %t797 = arith.andi %li, %j3 : i32
    %t798 = arith.muli %t797, %j4 : i32
    %t799 = arith.shrui %li, %j3 : i32
    %t800 = arith.addi %t798, %t799 : i32
    %s801, %unused801 = gpu.shuffle idx %lo795, %t800, %j32 : f64
    %s803, %unused803 = gpu.shuffle idx %hi796, %t800, %j32 : f64
    %t805 = arith.shrui %li, %j2 : i32
    %t806 = arith.andi %t805, %j1 : i32
    %t807 = arith.cmpi eq, %t806, %j1 : i32
    %sv808 = arith.select %t807, %s803, %s801 : f64
    %b809 = vector.insert %sv808, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m810 = nvgpu.mma.sync(%f794, %b809, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f811 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo812 = vector.extract %m793[0, 0] : f64 from vector<1x2xf64>
    %hi813 = vector.extract %m793[0, 1] : f64 from vector<1x2xf64>
    %t814 = arith.andi %li, %j3 : i32
    %t815 = arith.muli %t814, %j4 : i32
    %t816 = arith.shrui %li, %j3 : i32
    %t817 = arith.addi %t815, %t816 : i32
    %t818 = arith.addi %t817, %j16 : i32
    %s819, %unused819 = gpu.shuffle idx %lo812, %t818, %j32 : f64
    %s821, %unused821 = gpu.shuffle idx %hi813, %t818, %j32 : f64
    %t823 = arith.shrui %li, %j2 : i32
    %t824 = arith.andi %t823, %j1 : i32
    %t825 = arith.cmpi eq, %t824, %j1 : i32
    %sv826 = arith.select %t825, %s821, %s819 : f64
    %b827 = vector.insert %sv826, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m828 = nvgpu.mma.sync(%f811, %b827, %m810) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %sv829 = memref.subview %el2[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv830 = memref.subview %el2[7, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %y831 = vector.transfer_read %sv829[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %y832 = arith.subf %y831, %m828 : vector<1x2xf64>
    vector.transfer_write %y832, %sv829[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    %y833 = vector.transfer_read %sv830[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %y834 = arith.addf %y833, %m828 : vector<1x2xf64>
    vector.transfer_write %y834, %sv830[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    gpu.barrier
    %sv835 = memref.subview %el1[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc836 = arith.constant dense<0.0> : vector<8x8xf64>
      %a837 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b838 = vector.transfer_read %sv835[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r839 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a837, %b838, %acc836 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a840 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b841 = vector.transfer_read %sv835[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r842 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a840, %b841, %r839 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r842, %interp_all[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv843 = memref.subview %el1[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc844 = arith.constant dense<0.0> : vector<8x8xf64>
      %a845 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b846 = vector.transfer_read %sv843[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r847 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a845, %b846, %acc844 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a848 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b849 = vector.transfer_read %sv843[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r850 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a848, %b849, %r847 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r850, %interp_all[%c0, %c8] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv851 = memref.subview %el1[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc852 = arith.constant dense<0.0> : vector<8x8xf64>
      %a853 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b854 = vector.transfer_read %sv851[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r855 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a853, %b854, %acc852 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a856 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b857 = vector.transfer_read %sv851[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r858 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a856, %b857, %r855 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r858, %interp_all[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv859 = memref.subview %el1[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc860 = arith.constant dense<0.0> : vector<8x8xf64>
      %a861 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b862 = vector.transfer_read %sv859[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r863 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a861, %b862, %acc860 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a864 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b865 = vector.transfer_read %sv859[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r866 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a864, %b865, %r863 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r866, %interp_all[%c0, %c24] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv867 = memref.subview %el1[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc868 = arith.constant dense<0.0> : vector<8x8xf64>
      %a869 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b870 = vector.transfer_read %sv867[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r871 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a869, %b870, %acc868 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a872 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b873 = vector.transfer_read %sv867[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r874 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a872, %b873, %r871 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r874, %interp_all[%c0, %c32] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv875 = memref.subview %el1[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc876 = arith.constant dense<0.0> : vector<8x8xf64>
      %a877 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b878 = vector.transfer_read %sv875[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r879 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a877, %b878, %acc876 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a880 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b881 = vector.transfer_read %sv875[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r882 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a880, %b881, %r879 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r882, %interp_all[%c0, %c40] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv883 = memref.subview %el1[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc884 = arith.constant dense<0.0> : vector<8x8xf64>
      %a885 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b886 = vector.transfer_read %sv883[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r887 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a885, %b886, %acc884 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a888 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b889 = vector.transfer_read %sv883[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r890 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a888, %b889, %r887 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r890, %interp_all[%c0, %c48] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv891 = memref.subview %el1[7, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc892 = arith.constant dense<0.0> : vector<8x8xf64>
      %a893 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b894 = vector.transfer_read %sv891[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r895 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a893, %b894, %acc892 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a896 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b897 = vector.transfer_read %sv891[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r898 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a896, %b897, %r895 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r898, %interp_all[%c0, %c56] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv899 = memref.subview %el1[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc900 = arith.constant dense<0.0> : vector<8x8xf64>
      %a901 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b902 = vector.transfer_read %sv899[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r903 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a901, %b902, %acc900 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a904 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b905 = vector.transfer_read %sv899[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r906 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a904, %b905, %r903 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r906, %deriv_all[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv907 = memref.subview %el1[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc908 = arith.constant dense<0.0> : vector<8x8xf64>
      %a909 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b910 = vector.transfer_read %sv907[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r911 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a909, %b910, %acc908 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a912 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b913 = vector.transfer_read %sv907[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r914 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a912, %b913, %r911 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r914, %deriv_all[%c0, %c8] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv915 = memref.subview %el1[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc916 = arith.constant dense<0.0> : vector<8x8xf64>
      %a917 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b918 = vector.transfer_read %sv915[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r919 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a917, %b918, %acc916 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a920 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b921 = vector.transfer_read %sv915[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r922 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a920, %b921, %r919 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r922, %deriv_all[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv923 = memref.subview %el1[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc924 = arith.constant dense<0.0> : vector<8x8xf64>
      %a925 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b926 = vector.transfer_read %sv923[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r927 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a925, %b926, %acc924 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a928 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b929 = vector.transfer_read %sv923[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r930 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a928, %b929, %r927 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r930, %deriv_all[%c0, %c24] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv931 = memref.subview %el1[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc932 = arith.constant dense<0.0> : vector<8x8xf64>
      %a933 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b934 = vector.transfer_read %sv931[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r935 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a933, %b934, %acc932 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a936 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b937 = vector.transfer_read %sv931[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r938 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a936, %b937, %r935 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r938, %deriv_all[%c0, %c32] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv939 = memref.subview %el1[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc940 = arith.constant dense<0.0> : vector<8x8xf64>
      %a941 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b942 = vector.transfer_read %sv939[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r943 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a941, %b942, %acc940 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a944 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b945 = vector.transfer_read %sv939[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r946 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a944, %b945, %r943 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r946, %deriv_all[%c0, %c40] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv947 = memref.subview %el1[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc948 = arith.constant dense<0.0> : vector<8x8xf64>
      %a949 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b950 = vector.transfer_read %sv947[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r951 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a949, %b950, %acc948 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a952 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b953 = vector.transfer_read %sv947[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r954 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a952, %b953, %r951 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r954, %deriv_all[%c0, %c48] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv955 = memref.subview %el1[7, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc956 = arith.constant dense<0.0> : vector<8x8xf64>
      %a957 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b958 = vector.transfer_read %sv955[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r959 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a957, %b958, %acc956 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a960 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b961 = vector.transfer_read %sv955[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r962 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a960, %b961, %r959 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r962, %deriv_all[%c0, %c56] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %v3963 = memref.reinterpret_cast %interp_all to offset: [0], sizes: [8, 8, 8], strides: [64, 8, 1] : memref<8x64xf64, #gpu.address_space<workgroup>> to memref<8x8x8xf64, #gpu.address_space<workgroup>>
    %v3964 = memref.reinterpret_cast %deriv_all to offset: [0], sizes: [8, 8, 8], strides: [64, 8, 1] : memref<8x64xf64, #gpu.address_space<workgroup>> to memref<8x8x8xf64, #gpu.address_space<workgroup>>
    %sv965 = memref.subview %v3963[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>
    %sv966 = memref.subview %v3964[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>
    %f967 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f968 = vector.transfer_read %sv965[%k, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m969 = nvgpu.mma.sync(%f967, %f968, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f970 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f971 = vector.transfer_read %sv965[%k4, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m972 = nvgpu.mma.sync(%f970, %f971, %m969) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f973 = vector.transfer_read %sv965[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f974 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m975 = nvgpu.mma.sync(%f973, %f974, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f976 = vector.transfer_read %sv965[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f977 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m978 = nvgpu.mma.sync(%f976, %f977, %m975) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %cf979 = vector.transfer_read %sv966[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<1x2xf64>
    %cf980 = vector.transfer_read %el8[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf981 = vector.transfer_read %el6[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf982 = vector.transfer_read %el7[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %x983 = arith.mulf %cf980, %cf979 : vector<1x2xf64>
    %x984 = arith.mulf %cf981, %m978 : vector<1x2xf64>
    %x985 = arith.mulf %cf982, %m972 : vector<1x2xf64>
    %x986 = arith.addf %x983, %x984 : vector<1x2xf64>
    %fl987 = arith.addf %x986, %x985 : vector<1x2xf64>
    %lo988 = vector.extract %fl987[0, 0] : f64 from vector<1x2xf64>
    %hi989 = vector.extract %fl987[0, 1] : f64 from vector<1x2xf64>
    %t991 = arith.andi %li, %j28 : i32
    %t992 = arith.andi %li, %j3 : i32
    %t993 = arith.shrui %t992, %j1 : i32
    %t994 = arith.ori %t991, %t993 : i32
    %s995, %unused995 = gpu.shuffle idx %lo988, %t994, %j32 : f64
    %s997, %unused997 = gpu.shuffle idx %hi989, %t994, %j32 : f64
    %t999 = arith.andi %li, %j1 : i32
    %t1000 = arith.cmpi eq, %t999, %j1 : i32
    %sv1001 = arith.select %t1000, %s997, %s995 : f64
    %a1002 = vector.insert %sv1001, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f1003 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1004 = nvgpu.mma.sync(%a1002, %f1003, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %lo1005 = vector.extract %fl987[0, 0] : f64 from vector<1x2xf64>
    %hi1006 = vector.extract %fl987[0, 1] : f64 from vector<1x2xf64>
    %t1008 = arith.andi %li, %j28 : i32
    %t1009 = arith.andi %li, %j3 : i32
    %t1010 = arith.shrui %t1009, %j1 : i32
    %t1011 = arith.ori %t1008, %t1010 : i32
    %t1012 = arith.ori %t1011, %j2 : i32
    %s1013, %unused1013 = gpu.shuffle idx %lo1005, %t1012, %j32 : f64
    %s1015, %unused1015 = gpu.shuffle idx %hi1006, %t1012, %j32 : f64
    %t1017 = arith.andi %li, %j1 : i32
    %t1018 = arith.cmpi eq, %t1017, %j1 : i32
    %sv1019 = arith.select %t1018, %s1015, %s1013 : f64
    %a1020 = vector.insert %sv1019, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f1021 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1022 = nvgpu.mma.sync(%a1020, %f1021, %m1004) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1023 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo1024 = vector.extract %m1022[0, 0] : f64 from vector<1x2xf64>
    %hi1025 = vector.extract %m1022[0, 1] : f64 from vector<1x2xf64>
    %t1026 = arith.andi %li, %j3 : i32
    %t1027 = arith.muli %t1026, %j4 : i32
    %t1028 = arith.shrui %li, %j3 : i32
    %t1029 = arith.addi %t1027, %t1028 : i32
    %s1030, %unused1030 = gpu.shuffle idx %lo1024, %t1029, %j32 : f64
    %s1032, %unused1032 = gpu.shuffle idx %hi1025, %t1029, %j32 : f64
    %t1034 = arith.shrui %li, %j2 : i32
    %t1035 = arith.andi %t1034, %j1 : i32
    %t1036 = arith.cmpi eq, %t1035, %j1 : i32
    %sv1037 = arith.select %t1036, %s1032, %s1030 : f64
    %b1038 = vector.insert %sv1037, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m1039 = nvgpu.mma.sync(%f1023, %b1038, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1040 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo1041 = vector.extract %m1022[0, 0] : f64 from vector<1x2xf64>
    %hi1042 = vector.extract %m1022[0, 1] : f64 from vector<1x2xf64>
    %t1043 = arith.andi %li, %j3 : i32
    %t1044 = arith.muli %t1043, %j4 : i32
    %t1045 = arith.shrui %li, %j3 : i32
    %t1046 = arith.addi %t1044, %t1045 : i32
    %t1047 = arith.addi %t1046, %j16 : i32
    %s1048, %unused1048 = gpu.shuffle idx %lo1041, %t1047, %j32 : f64
    %s1050, %unused1050 = gpu.shuffle idx %hi1042, %t1047, %j32 : f64
    %t1052 = arith.shrui %li, %j2 : i32
    %t1053 = arith.andi %t1052, %j1 : i32
    %t1054 = arith.cmpi eq, %t1053, %j1 : i32
    %sv1055 = arith.select %t1054, %s1050, %s1048 : f64
    %b1056 = vector.insert %sv1055, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m1057 = nvgpu.mma.sync(%f1040, %b1056, %m1039) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %sv1058 = memref.subview %el2[0, 0, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv1059 = memref.subview %el2[0, 1, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %y1060 = vector.transfer_read %sv1058[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<1x2xf64>
    %y1061 = arith.subf %y1060, %m1057 : vector<1x2xf64>
    vector.transfer_write %y1061, %sv1058[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    %y1062 = vector.transfer_read %sv1059[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<1x2xf64>
    %y1063 = arith.addf %y1062, %m1057 : vector<1x2xf64>
    vector.transfer_write %y1063, %sv1059[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv1064 = memref.subview %v3963[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv1065 = memref.subview %v3964[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %f1066 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f1067 = vector.transfer_read %sv1064[%k, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m1068 = nvgpu.mma.sync(%f1066, %f1067, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1069 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f1070 = vector.transfer_read %sv1064[%k4, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m1071 = nvgpu.mma.sync(%f1069, %f1070, %m1068) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1072 = vector.transfer_read %sv1064[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f1073 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1074 = nvgpu.mma.sync(%f1072, %f1073, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1075 = vector.transfer_read %sv1064[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f1076 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1077 = nvgpu.mma.sync(%f1075, %f1076, %m1074) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %cf1078 = vector.transfer_read %sv1065[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<1x2xf64>
    %cf1079 = vector.transfer_read %el8[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf1080 = vector.transfer_read %el6[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf1081 = vector.transfer_read %el7[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %x1082 = arith.mulf %cf1079, %cf1078 : vector<1x2xf64>
    %x1083 = arith.mulf %cf1080, %m1077 : vector<1x2xf64>
    %x1084 = arith.mulf %cf1081, %m1071 : vector<1x2xf64>
    %x1085 = arith.addf %x1082, %x1083 : vector<1x2xf64>
    %fl1086 = arith.addf %x1085, %x1084 : vector<1x2xf64>
    %lo1087 = vector.extract %fl1086[0, 0] : f64 from vector<1x2xf64>
    %hi1088 = vector.extract %fl1086[0, 1] : f64 from vector<1x2xf64>
    %t1090 = arith.andi %li, %j28 : i32
    %t1091 = arith.andi %li, %j3 : i32
    %t1092 = arith.shrui %t1091, %j1 : i32
    %t1093 = arith.ori %t1090, %t1092 : i32
    %s1094, %unused1094 = gpu.shuffle idx %lo1087, %t1093, %j32 : f64
    %s1096, %unused1096 = gpu.shuffle idx %hi1088, %t1093, %j32 : f64
    %t1098 = arith.andi %li, %j1 : i32
    %t1099 = arith.cmpi eq, %t1098, %j1 : i32
    %sv1100 = arith.select %t1099, %s1096, %s1094 : f64
    %a1101 = vector.insert %sv1100, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f1102 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1103 = nvgpu.mma.sync(%a1101, %f1102, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %lo1104 = vector.extract %fl1086[0, 0] : f64 from vector<1x2xf64>
    %hi1105 = vector.extract %fl1086[0, 1] : f64 from vector<1x2xf64>
    %t1107 = arith.andi %li, %j28 : i32
    %t1108 = arith.andi %li, %j3 : i32
    %t1109 = arith.shrui %t1108, %j1 : i32
    %t1110 = arith.ori %t1107, %t1109 : i32
    %t1111 = arith.ori %t1110, %j2 : i32
    %s1112, %unused1112 = gpu.shuffle idx %lo1104, %t1111, %j32 : f64
    %s1114, %unused1114 = gpu.shuffle idx %hi1105, %t1111, %j32 : f64
    %t1116 = arith.andi %li, %j1 : i32
    %t1117 = arith.cmpi eq, %t1116, %j1 : i32
    %sv1118 = arith.select %t1117, %s1114, %s1112 : f64
    %a1119 = vector.insert %sv1118, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f1120 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1121 = nvgpu.mma.sync(%a1119, %f1120, %m1103) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1122 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo1123 = vector.extract %m1121[0, 0] : f64 from vector<1x2xf64>
    %hi1124 = vector.extract %m1121[0, 1] : f64 from vector<1x2xf64>
    %t1125 = arith.andi %li, %j3 : i32
    %t1126 = arith.muli %t1125, %j4 : i32
    %t1127 = arith.shrui %li, %j3 : i32
    %t1128 = arith.addi %t1126, %t1127 : i32
    %s1129, %unused1129 = gpu.shuffle idx %lo1123, %t1128, %j32 : f64
    %s1131, %unused1131 = gpu.shuffle idx %hi1124, %t1128, %j32 : f64
    %t1133 = arith.shrui %li, %j2 : i32
    %t1134 = arith.andi %t1133, %j1 : i32
    %t1135 = arith.cmpi eq, %t1134, %j1 : i32
    %sv1136 = arith.select %t1135, %s1131, %s1129 : f64
    %b1137 = vector.insert %sv1136, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m1138 = nvgpu.mma.sync(%f1122, %b1137, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1139 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo1140 = vector.extract %m1121[0, 0] : f64 from vector<1x2xf64>
    %hi1141 = vector.extract %m1121[0, 1] : f64 from vector<1x2xf64>
    %t1142 = arith.andi %li, %j3 : i32
    %t1143 = arith.muli %t1142, %j4 : i32
    %t1144 = arith.shrui %li, %j3 : i32
    %t1145 = arith.addi %t1143, %t1144 : i32
    %t1146 = arith.addi %t1145, %j16 : i32
    %s1147, %unused1147 = gpu.shuffle idx %lo1140, %t1146, %j32 : f64
    %s1149, %unused1149 = gpu.shuffle idx %hi1141, %t1146, %j32 : f64
    %t1151 = arith.shrui %li, %j2 : i32
    %t1152 = arith.andi %t1151, %j1 : i32
    %t1153 = arith.cmpi eq, %t1152, %j1 : i32
    %sv1154 = arith.select %t1153, %s1149, %s1147 : f64
    %b1155 = vector.insert %sv1154, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m1156 = nvgpu.mma.sync(%f1139, %b1155, %m1138) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %sv1157 = memref.subview %el2[0, 1, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv1158 = memref.subview %el2[0, 2, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %y1159 = vector.transfer_read %sv1157[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<1x2xf64>
    %y1160 = arith.subf %y1159, %m1156 : vector<1x2xf64>
    vector.transfer_write %y1160, %sv1157[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    %y1161 = vector.transfer_read %sv1158[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<1x2xf64>
    %y1162 = arith.addf %y1161, %m1156 : vector<1x2xf64>
    vector.transfer_write %y1162, %sv1158[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv1163 = memref.subview %v3963[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv1164 = memref.subview %v3964[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %f1165 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f1166 = vector.transfer_read %sv1163[%k, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m1167 = nvgpu.mma.sync(%f1165, %f1166, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1168 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f1169 = vector.transfer_read %sv1163[%k4, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m1170 = nvgpu.mma.sync(%f1168, %f1169, %m1167) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1171 = vector.transfer_read %sv1163[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f1172 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1173 = nvgpu.mma.sync(%f1171, %f1172, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1174 = vector.transfer_read %sv1163[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f1175 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1176 = nvgpu.mma.sync(%f1174, %f1175, %m1173) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %cf1177 = vector.transfer_read %sv1164[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<1x2xf64>
    %cf1178 = vector.transfer_read %el8[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf1179 = vector.transfer_read %el6[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf1180 = vector.transfer_read %el7[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %x1181 = arith.mulf %cf1178, %cf1177 : vector<1x2xf64>
    %x1182 = arith.mulf %cf1179, %m1176 : vector<1x2xf64>
    %x1183 = arith.mulf %cf1180, %m1170 : vector<1x2xf64>
    %x1184 = arith.addf %x1181, %x1182 : vector<1x2xf64>
    %fl1185 = arith.addf %x1184, %x1183 : vector<1x2xf64>
    %lo1186 = vector.extract %fl1185[0, 0] : f64 from vector<1x2xf64>
    %hi1187 = vector.extract %fl1185[0, 1] : f64 from vector<1x2xf64>
    %t1189 = arith.andi %li, %j28 : i32
    %t1190 = arith.andi %li, %j3 : i32
    %t1191 = arith.shrui %t1190, %j1 : i32
    %t1192 = arith.ori %t1189, %t1191 : i32
    %s1193, %unused1193 = gpu.shuffle idx %lo1186, %t1192, %j32 : f64
    %s1195, %unused1195 = gpu.shuffle idx %hi1187, %t1192, %j32 : f64
    %t1197 = arith.andi %li, %j1 : i32
    %t1198 = arith.cmpi eq, %t1197, %j1 : i32
    %sv1199 = arith.select %t1198, %s1195, %s1193 : f64
    %a1200 = vector.insert %sv1199, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f1201 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1202 = nvgpu.mma.sync(%a1200, %f1201, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %lo1203 = vector.extract %fl1185[0, 0] : f64 from vector<1x2xf64>
    %hi1204 = vector.extract %fl1185[0, 1] : f64 from vector<1x2xf64>
    %t1206 = arith.andi %li, %j28 : i32
    %t1207 = arith.andi %li, %j3 : i32
    %t1208 = arith.shrui %t1207, %j1 : i32
    %t1209 = arith.ori %t1206, %t1208 : i32
    %t1210 = arith.ori %t1209, %j2 : i32
    %s1211, %unused1211 = gpu.shuffle idx %lo1203, %t1210, %j32 : f64
    %s1213, %unused1213 = gpu.shuffle idx %hi1204, %t1210, %j32 : f64
    %t1215 = arith.andi %li, %j1 : i32
    %t1216 = arith.cmpi eq, %t1215, %j1 : i32
    %sv1217 = arith.select %t1216, %s1213, %s1211 : f64
    %a1218 = vector.insert %sv1217, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f1219 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1220 = nvgpu.mma.sync(%a1218, %f1219, %m1202) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1221 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo1222 = vector.extract %m1220[0, 0] : f64 from vector<1x2xf64>
    %hi1223 = vector.extract %m1220[0, 1] : f64 from vector<1x2xf64>
    %t1224 = arith.andi %li, %j3 : i32
    %t1225 = arith.muli %t1224, %j4 : i32
    %t1226 = arith.shrui %li, %j3 : i32
    %t1227 = arith.addi %t1225, %t1226 : i32
    %s1228, %unused1228 = gpu.shuffle idx %lo1222, %t1227, %j32 : f64
    %s1230, %unused1230 = gpu.shuffle idx %hi1223, %t1227, %j32 : f64
    %t1232 = arith.shrui %li, %j2 : i32
    %t1233 = arith.andi %t1232, %j1 : i32
    %t1234 = arith.cmpi eq, %t1233, %j1 : i32
    %sv1235 = arith.select %t1234, %s1230, %s1228 : f64
    %b1236 = vector.insert %sv1235, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m1237 = nvgpu.mma.sync(%f1221, %b1236, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1238 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo1239 = vector.extract %m1220[0, 0] : f64 from vector<1x2xf64>
    %hi1240 = vector.extract %m1220[0, 1] : f64 from vector<1x2xf64>
    %t1241 = arith.andi %li, %j3 : i32
    %t1242 = arith.muli %t1241, %j4 : i32
    %t1243 = arith.shrui %li, %j3 : i32
    %t1244 = arith.addi %t1242, %t1243 : i32
    %t1245 = arith.addi %t1244, %j16 : i32
    %s1246, %unused1246 = gpu.shuffle idx %lo1239, %t1245, %j32 : f64
    %s1248, %unused1248 = gpu.shuffle idx %hi1240, %t1245, %j32 : f64
    %t1250 = arith.shrui %li, %j2 : i32
    %t1251 = arith.andi %t1250, %j1 : i32
    %t1252 = arith.cmpi eq, %t1251, %j1 : i32
    %sv1253 = arith.select %t1252, %s1248, %s1246 : f64
    %b1254 = vector.insert %sv1253, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m1255 = nvgpu.mma.sync(%f1238, %b1254, %m1237) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %sv1256 = memref.subview %el2[0, 2, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv1257 = memref.subview %el2[0, 3, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %y1258 = vector.transfer_read %sv1256[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<1x2xf64>
    %y1259 = arith.subf %y1258, %m1255 : vector<1x2xf64>
    vector.transfer_write %y1259, %sv1256[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    %y1260 = vector.transfer_read %sv1257[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<1x2xf64>
    %y1261 = arith.addf %y1260, %m1255 : vector<1x2xf64>
    vector.transfer_write %y1261, %sv1257[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv1262 = memref.subview %v3963[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv1263 = memref.subview %v3964[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %f1264 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f1265 = vector.transfer_read %sv1262[%k, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m1266 = nvgpu.mma.sync(%f1264, %f1265, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1267 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f1268 = vector.transfer_read %sv1262[%k4, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m1269 = nvgpu.mma.sync(%f1267, %f1268, %m1266) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1270 = vector.transfer_read %sv1262[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f1271 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1272 = nvgpu.mma.sync(%f1270, %f1271, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1273 = vector.transfer_read %sv1262[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f1274 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1275 = nvgpu.mma.sync(%f1273, %f1274, %m1272) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %cf1276 = vector.transfer_read %sv1263[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<1x2xf64>
    %cf1277 = vector.transfer_read %el8[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf1278 = vector.transfer_read %el6[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf1279 = vector.transfer_read %el7[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %x1280 = arith.mulf %cf1277, %cf1276 : vector<1x2xf64>
    %x1281 = arith.mulf %cf1278, %m1275 : vector<1x2xf64>
    %x1282 = arith.mulf %cf1279, %m1269 : vector<1x2xf64>
    %x1283 = arith.addf %x1280, %x1281 : vector<1x2xf64>
    %fl1284 = arith.addf %x1283, %x1282 : vector<1x2xf64>
    %lo1285 = vector.extract %fl1284[0, 0] : f64 from vector<1x2xf64>
    %hi1286 = vector.extract %fl1284[0, 1] : f64 from vector<1x2xf64>
    %t1288 = arith.andi %li, %j28 : i32
    %t1289 = arith.andi %li, %j3 : i32
    %t1290 = arith.shrui %t1289, %j1 : i32
    %t1291 = arith.ori %t1288, %t1290 : i32
    %s1292, %unused1292 = gpu.shuffle idx %lo1285, %t1291, %j32 : f64
    %s1294, %unused1294 = gpu.shuffle idx %hi1286, %t1291, %j32 : f64
    %t1296 = arith.andi %li, %j1 : i32
    %t1297 = arith.cmpi eq, %t1296, %j1 : i32
    %sv1298 = arith.select %t1297, %s1294, %s1292 : f64
    %a1299 = vector.insert %sv1298, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f1300 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1301 = nvgpu.mma.sync(%a1299, %f1300, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %lo1302 = vector.extract %fl1284[0, 0] : f64 from vector<1x2xf64>
    %hi1303 = vector.extract %fl1284[0, 1] : f64 from vector<1x2xf64>
    %t1305 = arith.andi %li, %j28 : i32
    %t1306 = arith.andi %li, %j3 : i32
    %t1307 = arith.shrui %t1306, %j1 : i32
    %t1308 = arith.ori %t1305, %t1307 : i32
    %t1309 = arith.ori %t1308, %j2 : i32
    %s1310, %unused1310 = gpu.shuffle idx %lo1302, %t1309, %j32 : f64
    %s1312, %unused1312 = gpu.shuffle idx %hi1303, %t1309, %j32 : f64
    %t1314 = arith.andi %li, %j1 : i32
    %t1315 = arith.cmpi eq, %t1314, %j1 : i32
    %sv1316 = arith.select %t1315, %s1312, %s1310 : f64
    %a1317 = vector.insert %sv1316, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f1318 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1319 = nvgpu.mma.sync(%a1317, %f1318, %m1301) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1320 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo1321 = vector.extract %m1319[0, 0] : f64 from vector<1x2xf64>
    %hi1322 = vector.extract %m1319[0, 1] : f64 from vector<1x2xf64>
    %t1323 = arith.andi %li, %j3 : i32
    %t1324 = arith.muli %t1323, %j4 : i32
    %t1325 = arith.shrui %li, %j3 : i32
    %t1326 = arith.addi %t1324, %t1325 : i32
    %s1327, %unused1327 = gpu.shuffle idx %lo1321, %t1326, %j32 : f64
    %s1329, %unused1329 = gpu.shuffle idx %hi1322, %t1326, %j32 : f64
    %t1331 = arith.shrui %li, %j2 : i32
    %t1332 = arith.andi %t1331, %j1 : i32
    %t1333 = arith.cmpi eq, %t1332, %j1 : i32
    %sv1334 = arith.select %t1333, %s1329, %s1327 : f64
    %b1335 = vector.insert %sv1334, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m1336 = nvgpu.mma.sync(%f1320, %b1335, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1337 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo1338 = vector.extract %m1319[0, 0] : f64 from vector<1x2xf64>
    %hi1339 = vector.extract %m1319[0, 1] : f64 from vector<1x2xf64>
    %t1340 = arith.andi %li, %j3 : i32
    %t1341 = arith.muli %t1340, %j4 : i32
    %t1342 = arith.shrui %li, %j3 : i32
    %t1343 = arith.addi %t1341, %t1342 : i32
    %t1344 = arith.addi %t1343, %j16 : i32
    %s1345, %unused1345 = gpu.shuffle idx %lo1338, %t1344, %j32 : f64
    %s1347, %unused1347 = gpu.shuffle idx %hi1339, %t1344, %j32 : f64
    %t1349 = arith.shrui %li, %j2 : i32
    %t1350 = arith.andi %t1349, %j1 : i32
    %t1351 = arith.cmpi eq, %t1350, %j1 : i32
    %sv1352 = arith.select %t1351, %s1347, %s1345 : f64
    %b1353 = vector.insert %sv1352, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m1354 = nvgpu.mma.sync(%f1337, %b1353, %m1336) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %sv1355 = memref.subview %el2[0, 3, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv1356 = memref.subview %el2[0, 4, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %y1357 = vector.transfer_read %sv1355[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<1x2xf64>
    %y1358 = arith.subf %y1357, %m1354 : vector<1x2xf64>
    vector.transfer_write %y1358, %sv1355[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    %y1359 = vector.transfer_read %sv1356[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<1x2xf64>
    %y1360 = arith.addf %y1359, %m1354 : vector<1x2xf64>
    vector.transfer_write %y1360, %sv1356[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv1361 = memref.subview %v3963[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv1362 = memref.subview %v3964[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %f1363 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f1364 = vector.transfer_read %sv1361[%k, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m1365 = nvgpu.mma.sync(%f1363, %f1364, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1366 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f1367 = vector.transfer_read %sv1361[%k4, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m1368 = nvgpu.mma.sync(%f1366, %f1367, %m1365) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1369 = vector.transfer_read %sv1361[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f1370 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1371 = nvgpu.mma.sync(%f1369, %f1370, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1372 = vector.transfer_read %sv1361[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f1373 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1374 = nvgpu.mma.sync(%f1372, %f1373, %m1371) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %cf1375 = vector.transfer_read %sv1362[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<1x2xf64>
    %cf1376 = vector.transfer_read %el8[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf1377 = vector.transfer_read %el6[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf1378 = vector.transfer_read %el7[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %x1379 = arith.mulf %cf1376, %cf1375 : vector<1x2xf64>
    %x1380 = arith.mulf %cf1377, %m1374 : vector<1x2xf64>
    %x1381 = arith.mulf %cf1378, %m1368 : vector<1x2xf64>
    %x1382 = arith.addf %x1379, %x1380 : vector<1x2xf64>
    %fl1383 = arith.addf %x1382, %x1381 : vector<1x2xf64>
    %lo1384 = vector.extract %fl1383[0, 0] : f64 from vector<1x2xf64>
    %hi1385 = vector.extract %fl1383[0, 1] : f64 from vector<1x2xf64>
    %t1387 = arith.andi %li, %j28 : i32
    %t1388 = arith.andi %li, %j3 : i32
    %t1389 = arith.shrui %t1388, %j1 : i32
    %t1390 = arith.ori %t1387, %t1389 : i32
    %s1391, %unused1391 = gpu.shuffle idx %lo1384, %t1390, %j32 : f64
    %s1393, %unused1393 = gpu.shuffle idx %hi1385, %t1390, %j32 : f64
    %t1395 = arith.andi %li, %j1 : i32
    %t1396 = arith.cmpi eq, %t1395, %j1 : i32
    %sv1397 = arith.select %t1396, %s1393, %s1391 : f64
    %a1398 = vector.insert %sv1397, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f1399 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1400 = nvgpu.mma.sync(%a1398, %f1399, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %lo1401 = vector.extract %fl1383[0, 0] : f64 from vector<1x2xf64>
    %hi1402 = vector.extract %fl1383[0, 1] : f64 from vector<1x2xf64>
    %t1404 = arith.andi %li, %j28 : i32
    %t1405 = arith.andi %li, %j3 : i32
    %t1406 = arith.shrui %t1405, %j1 : i32
    %t1407 = arith.ori %t1404, %t1406 : i32
    %t1408 = arith.ori %t1407, %j2 : i32
    %s1409, %unused1409 = gpu.shuffle idx %lo1401, %t1408, %j32 : f64
    %s1411, %unused1411 = gpu.shuffle idx %hi1402, %t1408, %j32 : f64
    %t1413 = arith.andi %li, %j1 : i32
    %t1414 = arith.cmpi eq, %t1413, %j1 : i32
    %sv1415 = arith.select %t1414, %s1411, %s1409 : f64
    %a1416 = vector.insert %sv1415, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f1417 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1418 = nvgpu.mma.sync(%a1416, %f1417, %m1400) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1419 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo1420 = vector.extract %m1418[0, 0] : f64 from vector<1x2xf64>
    %hi1421 = vector.extract %m1418[0, 1] : f64 from vector<1x2xf64>
    %t1422 = arith.andi %li, %j3 : i32
    %t1423 = arith.muli %t1422, %j4 : i32
    %t1424 = arith.shrui %li, %j3 : i32
    %t1425 = arith.addi %t1423, %t1424 : i32
    %s1426, %unused1426 = gpu.shuffle idx %lo1420, %t1425, %j32 : f64
    %s1428, %unused1428 = gpu.shuffle idx %hi1421, %t1425, %j32 : f64
    %t1430 = arith.shrui %li, %j2 : i32
    %t1431 = arith.andi %t1430, %j1 : i32
    %t1432 = arith.cmpi eq, %t1431, %j1 : i32
    %sv1433 = arith.select %t1432, %s1428, %s1426 : f64
    %b1434 = vector.insert %sv1433, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m1435 = nvgpu.mma.sync(%f1419, %b1434, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1436 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo1437 = vector.extract %m1418[0, 0] : f64 from vector<1x2xf64>
    %hi1438 = vector.extract %m1418[0, 1] : f64 from vector<1x2xf64>
    %t1439 = arith.andi %li, %j3 : i32
    %t1440 = arith.muli %t1439, %j4 : i32
    %t1441 = arith.shrui %li, %j3 : i32
    %t1442 = arith.addi %t1440, %t1441 : i32
    %t1443 = arith.addi %t1442, %j16 : i32
    %s1444, %unused1444 = gpu.shuffle idx %lo1437, %t1443, %j32 : f64
    %s1446, %unused1446 = gpu.shuffle idx %hi1438, %t1443, %j32 : f64
    %t1448 = arith.shrui %li, %j2 : i32
    %t1449 = arith.andi %t1448, %j1 : i32
    %t1450 = arith.cmpi eq, %t1449, %j1 : i32
    %sv1451 = arith.select %t1450, %s1446, %s1444 : f64
    %b1452 = vector.insert %sv1451, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m1453 = nvgpu.mma.sync(%f1436, %b1452, %m1435) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %sv1454 = memref.subview %el2[0, 4, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv1455 = memref.subview %el2[0, 5, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %y1456 = vector.transfer_read %sv1454[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<1x2xf64>
    %y1457 = arith.subf %y1456, %m1453 : vector<1x2xf64>
    vector.transfer_write %y1457, %sv1454[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    %y1458 = vector.transfer_read %sv1455[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<1x2xf64>
    %y1459 = arith.addf %y1458, %m1453 : vector<1x2xf64>
    vector.transfer_write %y1459, %sv1455[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv1460 = memref.subview %v3963[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv1461 = memref.subview %v3964[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %f1462 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f1463 = vector.transfer_read %sv1460[%k, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m1464 = nvgpu.mma.sync(%f1462, %f1463, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1465 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f1466 = vector.transfer_read %sv1460[%k4, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m1467 = nvgpu.mma.sync(%f1465, %f1466, %m1464) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1468 = vector.transfer_read %sv1460[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f1469 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1470 = nvgpu.mma.sync(%f1468, %f1469, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1471 = vector.transfer_read %sv1460[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f1472 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1473 = nvgpu.mma.sync(%f1471, %f1472, %m1470) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %cf1474 = vector.transfer_read %sv1461[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<1x2xf64>
    %cf1475 = vector.transfer_read %el8[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf1476 = vector.transfer_read %el6[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf1477 = vector.transfer_read %el7[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %x1478 = arith.mulf %cf1475, %cf1474 : vector<1x2xf64>
    %x1479 = arith.mulf %cf1476, %m1473 : vector<1x2xf64>
    %x1480 = arith.mulf %cf1477, %m1467 : vector<1x2xf64>
    %x1481 = arith.addf %x1478, %x1479 : vector<1x2xf64>
    %fl1482 = arith.addf %x1481, %x1480 : vector<1x2xf64>
    %lo1483 = vector.extract %fl1482[0, 0] : f64 from vector<1x2xf64>
    %hi1484 = vector.extract %fl1482[0, 1] : f64 from vector<1x2xf64>
    %t1486 = arith.andi %li, %j28 : i32
    %t1487 = arith.andi %li, %j3 : i32
    %t1488 = arith.shrui %t1487, %j1 : i32
    %t1489 = arith.ori %t1486, %t1488 : i32
    %s1490, %unused1490 = gpu.shuffle idx %lo1483, %t1489, %j32 : f64
    %s1492, %unused1492 = gpu.shuffle idx %hi1484, %t1489, %j32 : f64
    %t1494 = arith.andi %li, %j1 : i32
    %t1495 = arith.cmpi eq, %t1494, %j1 : i32
    %sv1496 = arith.select %t1495, %s1492, %s1490 : f64
    %a1497 = vector.insert %sv1496, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f1498 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1499 = nvgpu.mma.sync(%a1497, %f1498, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %lo1500 = vector.extract %fl1482[0, 0] : f64 from vector<1x2xf64>
    %hi1501 = vector.extract %fl1482[0, 1] : f64 from vector<1x2xf64>
    %t1503 = arith.andi %li, %j28 : i32
    %t1504 = arith.andi %li, %j3 : i32
    %t1505 = arith.shrui %t1504, %j1 : i32
    %t1506 = arith.ori %t1503, %t1505 : i32
    %t1507 = arith.ori %t1506, %j2 : i32
    %s1508, %unused1508 = gpu.shuffle idx %lo1500, %t1507, %j32 : f64
    %s1510, %unused1510 = gpu.shuffle idx %hi1501, %t1507, %j32 : f64
    %t1512 = arith.andi %li, %j1 : i32
    %t1513 = arith.cmpi eq, %t1512, %j1 : i32
    %sv1514 = arith.select %t1513, %s1510, %s1508 : f64
    %a1515 = vector.insert %sv1514, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f1516 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1517 = nvgpu.mma.sync(%a1515, %f1516, %m1499) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1518 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo1519 = vector.extract %m1517[0, 0] : f64 from vector<1x2xf64>
    %hi1520 = vector.extract %m1517[0, 1] : f64 from vector<1x2xf64>
    %t1521 = arith.andi %li, %j3 : i32
    %t1522 = arith.muli %t1521, %j4 : i32
    %t1523 = arith.shrui %li, %j3 : i32
    %t1524 = arith.addi %t1522, %t1523 : i32
    %s1525, %unused1525 = gpu.shuffle idx %lo1519, %t1524, %j32 : f64
    %s1527, %unused1527 = gpu.shuffle idx %hi1520, %t1524, %j32 : f64
    %t1529 = arith.shrui %li, %j2 : i32
    %t1530 = arith.andi %t1529, %j1 : i32
    %t1531 = arith.cmpi eq, %t1530, %j1 : i32
    %sv1532 = arith.select %t1531, %s1527, %s1525 : f64
    %b1533 = vector.insert %sv1532, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m1534 = nvgpu.mma.sync(%f1518, %b1533, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1535 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo1536 = vector.extract %m1517[0, 0] : f64 from vector<1x2xf64>
    %hi1537 = vector.extract %m1517[0, 1] : f64 from vector<1x2xf64>
    %t1538 = arith.andi %li, %j3 : i32
    %t1539 = arith.muli %t1538, %j4 : i32
    %t1540 = arith.shrui %li, %j3 : i32
    %t1541 = arith.addi %t1539, %t1540 : i32
    %t1542 = arith.addi %t1541, %j16 : i32
    %s1543, %unused1543 = gpu.shuffle idx %lo1536, %t1542, %j32 : f64
    %s1545, %unused1545 = gpu.shuffle idx %hi1537, %t1542, %j32 : f64
    %t1547 = arith.shrui %li, %j2 : i32
    %t1548 = arith.andi %t1547, %j1 : i32
    %t1549 = arith.cmpi eq, %t1548, %j1 : i32
    %sv1550 = arith.select %t1549, %s1545, %s1543 : f64
    %b1551 = vector.insert %sv1550, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m1552 = nvgpu.mma.sync(%f1535, %b1551, %m1534) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %sv1553 = memref.subview %el2[0, 5, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv1554 = memref.subview %el2[0, 6, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %y1555 = vector.transfer_read %sv1553[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<1x2xf64>
    %y1556 = arith.subf %y1555, %m1552 : vector<1x2xf64>
    vector.transfer_write %y1556, %sv1553[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    %y1557 = vector.transfer_read %sv1554[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<1x2xf64>
    %y1558 = arith.addf %y1557, %m1552 : vector<1x2xf64>
    vector.transfer_write %y1558, %sv1554[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv1559 = memref.subview %v3963[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv1560 = memref.subview %v3964[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %f1561 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f1562 = vector.transfer_read %sv1559[%k, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m1563 = nvgpu.mma.sync(%f1561, %f1562, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1564 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f1565 = vector.transfer_read %sv1559[%k4, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m1566 = nvgpu.mma.sync(%f1564, %f1565, %m1563) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1567 = vector.transfer_read %sv1559[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f1568 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1569 = nvgpu.mma.sync(%f1567, %f1568, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1570 = vector.transfer_read %sv1559[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f1571 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1572 = nvgpu.mma.sync(%f1570, %f1571, %m1569) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %cf1573 = vector.transfer_read %sv1560[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<1x2xf64>
    %cf1574 = vector.transfer_read %el8[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf1575 = vector.transfer_read %el6[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf1576 = vector.transfer_read %el7[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %x1577 = arith.mulf %cf1574, %cf1573 : vector<1x2xf64>
    %x1578 = arith.mulf %cf1575, %m1572 : vector<1x2xf64>
    %x1579 = arith.mulf %cf1576, %m1566 : vector<1x2xf64>
    %x1580 = arith.addf %x1577, %x1578 : vector<1x2xf64>
    %fl1581 = arith.addf %x1580, %x1579 : vector<1x2xf64>
    %lo1582 = vector.extract %fl1581[0, 0] : f64 from vector<1x2xf64>
    %hi1583 = vector.extract %fl1581[0, 1] : f64 from vector<1x2xf64>
    %t1585 = arith.andi %li, %j28 : i32
    %t1586 = arith.andi %li, %j3 : i32
    %t1587 = arith.shrui %t1586, %j1 : i32
    %t1588 = arith.ori %t1585, %t1587 : i32
    %s1589, %unused1589 = gpu.shuffle idx %lo1582, %t1588, %j32 : f64
    %s1591, %unused1591 = gpu.shuffle idx %hi1583, %t1588, %j32 : f64
    %t1593 = arith.andi %li, %j1 : i32
    %t1594 = arith.cmpi eq, %t1593, %j1 : i32
    %sv1595 = arith.select %t1594, %s1591, %s1589 : f64
    %a1596 = vector.insert %sv1595, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f1597 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1598 = nvgpu.mma.sync(%a1596, %f1597, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %lo1599 = vector.extract %fl1581[0, 0] : f64 from vector<1x2xf64>
    %hi1600 = vector.extract %fl1581[0, 1] : f64 from vector<1x2xf64>
    %t1602 = arith.andi %li, %j28 : i32
    %t1603 = arith.andi %li, %j3 : i32
    %t1604 = arith.shrui %t1603, %j1 : i32
    %t1605 = arith.ori %t1602, %t1604 : i32
    %t1606 = arith.ori %t1605, %j2 : i32
    %s1607, %unused1607 = gpu.shuffle idx %lo1599, %t1606, %j32 : f64
    %s1609, %unused1609 = gpu.shuffle idx %hi1600, %t1606, %j32 : f64
    %t1611 = arith.andi %li, %j1 : i32
    %t1612 = arith.cmpi eq, %t1611, %j1 : i32
    %sv1613 = arith.select %t1612, %s1609, %s1607 : f64
    %a1614 = vector.insert %sv1613, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f1615 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1616 = nvgpu.mma.sync(%a1614, %f1615, %m1598) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1617 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo1618 = vector.extract %m1616[0, 0] : f64 from vector<1x2xf64>
    %hi1619 = vector.extract %m1616[0, 1] : f64 from vector<1x2xf64>
    %t1620 = arith.andi %li, %j3 : i32
    %t1621 = arith.muli %t1620, %j4 : i32
    %t1622 = arith.shrui %li, %j3 : i32
    %t1623 = arith.addi %t1621, %t1622 : i32
    %s1624, %unused1624 = gpu.shuffle idx %lo1618, %t1623, %j32 : f64
    %s1626, %unused1626 = gpu.shuffle idx %hi1619, %t1623, %j32 : f64
    %t1628 = arith.shrui %li, %j2 : i32
    %t1629 = arith.andi %t1628, %j1 : i32
    %t1630 = arith.cmpi eq, %t1629, %j1 : i32
    %sv1631 = arith.select %t1630, %s1626, %s1624 : f64
    %b1632 = vector.insert %sv1631, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m1633 = nvgpu.mma.sync(%f1617, %b1632, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1634 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo1635 = vector.extract %m1616[0, 0] : f64 from vector<1x2xf64>
    %hi1636 = vector.extract %m1616[0, 1] : f64 from vector<1x2xf64>
    %t1637 = arith.andi %li, %j3 : i32
    %t1638 = arith.muli %t1637, %j4 : i32
    %t1639 = arith.shrui %li, %j3 : i32
    %t1640 = arith.addi %t1638, %t1639 : i32
    %t1641 = arith.addi %t1640, %j16 : i32
    %s1642, %unused1642 = gpu.shuffle idx %lo1635, %t1641, %j32 : f64
    %s1644, %unused1644 = gpu.shuffle idx %hi1636, %t1641, %j32 : f64
    %t1646 = arith.shrui %li, %j2 : i32
    %t1647 = arith.andi %t1646, %j1 : i32
    %t1648 = arith.cmpi eq, %t1647, %j1 : i32
    %sv1649 = arith.select %t1648, %s1644, %s1642 : f64
    %b1650 = vector.insert %sv1649, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m1651 = nvgpu.mma.sync(%f1634, %b1650, %m1633) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %sv1652 = memref.subview %el2[0, 6, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv1653 = memref.subview %el2[0, 7, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %y1654 = vector.transfer_read %sv1652[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<1x2xf64>
    %y1655 = arith.subf %y1654, %m1651 : vector<1x2xf64>
    vector.transfer_write %y1655, %sv1652[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    %y1656 = vector.transfer_read %sv1653[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<1x2xf64>
    %y1657 = arith.addf %y1656, %m1651 : vector<1x2xf64>
    vector.transfer_write %y1657, %sv1653[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    gpu.barrier
    %sv1658 = memref.subview %el1[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1659 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1660 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1661 = vector.transfer_read %sv1658[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1662 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1660, %b1661, %acc1659 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1663 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1664 = vector.transfer_read %sv1658[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1665 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1663, %b1664, %r1662 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1665, %interp_all[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1666 = memref.subview %el1[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1667 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1668 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1669 = vector.transfer_read %sv1666[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1670 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1668, %b1669, %acc1667 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1671 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1672 = vector.transfer_read %sv1666[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1673 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1671, %b1672, %r1670 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1673, %interp_all[%c0, %c8] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1674 = memref.subview %el1[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1675 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1676 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1677 = vector.transfer_read %sv1674[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1678 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1676, %b1677, %acc1675 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1679 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1680 = vector.transfer_read %sv1674[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1681 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1679, %b1680, %r1678 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1681, %interp_all[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1682 = memref.subview %el1[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1683 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1684 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1685 = vector.transfer_read %sv1682[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1686 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1684, %b1685, %acc1683 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1687 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1688 = vector.transfer_read %sv1682[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1689 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1687, %b1688, %r1686 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1689, %interp_all[%c0, %c24] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1690 = memref.subview %el1[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1691 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1692 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1693 = vector.transfer_read %sv1690[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1694 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1692, %b1693, %acc1691 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1695 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1696 = vector.transfer_read %sv1690[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1697 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1695, %b1696, %r1694 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1697, %interp_all[%c0, %c32] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1698 = memref.subview %el1[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1699 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1700 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1701 = vector.transfer_read %sv1698[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1702 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1700, %b1701, %acc1699 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1703 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1704 = vector.transfer_read %sv1698[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1705 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1703, %b1704, %r1702 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1705, %interp_all[%c0, %c40] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1706 = memref.subview %el1[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1707 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1708 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1709 = vector.transfer_read %sv1706[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1710 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1708, %b1709, %acc1707 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1711 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1712 = vector.transfer_read %sv1706[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1713 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1711, %b1712, %r1710 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1713, %interp_all[%c0, %c48] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1714 = memref.subview %el1[7, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1715 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1716 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1717 = vector.transfer_read %sv1714[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1718 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1716, %b1717, %acc1715 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1719 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1720 = vector.transfer_read %sv1714[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1721 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1719, %b1720, %r1718 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1721, %interp_all[%c0, %c56] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1722 = memref.subview %el1[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1723 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1724 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1725 = vector.transfer_read %sv1722[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1726 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1724, %b1725, %acc1723 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1727 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1728 = vector.transfer_read %sv1722[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1729 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1727, %b1728, %r1726 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1729, %deriv_all[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1730 = memref.subview %el1[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1731 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1732 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1733 = vector.transfer_read %sv1730[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1734 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1732, %b1733, %acc1731 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1735 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1736 = vector.transfer_read %sv1730[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1737 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1735, %b1736, %r1734 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1737, %deriv_all[%c0, %c8] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1738 = memref.subview %el1[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1739 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1740 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1741 = vector.transfer_read %sv1738[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1742 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1740, %b1741, %acc1739 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1743 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1744 = vector.transfer_read %sv1738[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1745 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1743, %b1744, %r1742 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1745, %deriv_all[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1746 = memref.subview %el1[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1747 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1748 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1749 = vector.transfer_read %sv1746[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1750 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1748, %b1749, %acc1747 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1751 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1752 = vector.transfer_read %sv1746[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1753 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1751, %b1752, %r1750 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1753, %deriv_all[%c0, %c24] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1754 = memref.subview %el1[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1755 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1756 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1757 = vector.transfer_read %sv1754[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1758 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1756, %b1757, %acc1755 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1759 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1760 = vector.transfer_read %sv1754[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1761 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1759, %b1760, %r1758 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1761, %deriv_all[%c0, %c32] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1762 = memref.subview %el1[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1763 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1764 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1765 = vector.transfer_read %sv1762[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1766 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1764, %b1765, %acc1763 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1767 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1768 = vector.transfer_read %sv1762[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1769 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1767, %b1768, %r1766 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1769, %deriv_all[%c0, %c40] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1770 = memref.subview %el1[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1771 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1772 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1773 = vector.transfer_read %sv1770[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1774 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1772, %b1773, %acc1771 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1775 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1776 = vector.transfer_read %sv1770[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1777 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1775, %b1776, %r1774 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1777, %deriv_all[%c0, %c48] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1778 = memref.subview %el1[7, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1779 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1780 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1781 = vector.transfer_read %sv1778[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1782 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1780, %b1781, %acc1779 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1783 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1784 = vector.transfer_read %sv1778[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1785 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1783, %b1784, %r1782 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1785, %deriv_all[%c0, %c56] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %v31786 = memref.reinterpret_cast %interp_all to offset: [0], sizes: [8, 8, 8], strides: [64, 8, 1] : memref<8x64xf64, #gpu.address_space<workgroup>> to memref<8x8x8xf64, #gpu.address_space<workgroup>>
    %v31787 = memref.reinterpret_cast %deriv_all to offset: [0], sizes: [8, 8, 8], strides: [64, 8, 1] : memref<8x64xf64, #gpu.address_space<workgroup>> to memref<8x8x8xf64, #gpu.address_space<workgroup>>
    %sv1788 = memref.subview %v31786[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>
    %sv1789 = memref.subview %v31787[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>
    %f1790 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f1791 = vector.transfer_read %sv1788[%k, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m1792 = nvgpu.mma.sync(%f1790, %f1791, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1793 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f1794 = vector.transfer_read %sv1788[%k4, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m1795 = nvgpu.mma.sync(%f1793, %f1794, %m1792) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1796 = vector.transfer_read %sv1788[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f1797 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1798 = nvgpu.mma.sync(%f1796, %f1797, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1799 = vector.transfer_read %sv1788[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f1800 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1801 = nvgpu.mma.sync(%f1799, %f1800, %m1798) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %cf1802 = vector.transfer_read %sv1789[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<1x2xf64>
    %cf1803 = vector.transfer_read %el11[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf1804 = vector.transfer_read %el9[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf1805 = vector.transfer_read %el10[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %x1806 = arith.mulf %cf1803, %cf1802 : vector<1x2xf64>
    %x1807 = arith.mulf %cf1804, %m1801 : vector<1x2xf64>
    %x1808 = arith.mulf %cf1805, %m1795 : vector<1x2xf64>
    %x1809 = arith.addf %x1806, %x1807 : vector<1x2xf64>
    %fl1810 = arith.addf %x1809, %x1808 : vector<1x2xf64>
    %lo1811 = vector.extract %fl1810[0, 0] : f64 from vector<1x2xf64>
    %hi1812 = vector.extract %fl1810[0, 1] : f64 from vector<1x2xf64>
    %t1814 = arith.andi %li, %j28 : i32
    %t1815 = arith.andi %li, %j3 : i32
    %t1816 = arith.shrui %t1815, %j1 : i32
    %t1817 = arith.ori %t1814, %t1816 : i32
    %s1818, %unused1818 = gpu.shuffle idx %lo1811, %t1817, %j32 : f64
    %s1820, %unused1820 = gpu.shuffle idx %hi1812, %t1817, %j32 : f64
    %t1822 = arith.andi %li, %j1 : i32
    %t1823 = arith.cmpi eq, %t1822, %j1 : i32
    %sv1824 = arith.select %t1823, %s1820, %s1818 : f64
    %a1825 = vector.insert %sv1824, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f1826 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1827 = nvgpu.mma.sync(%a1825, %f1826, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %lo1828 = vector.extract %fl1810[0, 0] : f64 from vector<1x2xf64>
    %hi1829 = vector.extract %fl1810[0, 1] : f64 from vector<1x2xf64>
    %t1831 = arith.andi %li, %j28 : i32
    %t1832 = arith.andi %li, %j3 : i32
    %t1833 = arith.shrui %t1832, %j1 : i32
    %t1834 = arith.ori %t1831, %t1833 : i32
    %t1835 = arith.ori %t1834, %j2 : i32
    %s1836, %unused1836 = gpu.shuffle idx %lo1828, %t1835, %j32 : f64
    %s1838, %unused1838 = gpu.shuffle idx %hi1829, %t1835, %j32 : f64
    %t1840 = arith.andi %li, %j1 : i32
    %t1841 = arith.cmpi eq, %t1840, %j1 : i32
    %sv1842 = arith.select %t1841, %s1838, %s1836 : f64
    %a1843 = vector.insert %sv1842, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f1844 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1845 = nvgpu.mma.sync(%a1843, %f1844, %m1827) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1846 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo1847 = vector.extract %m1845[0, 0] : f64 from vector<1x2xf64>
    %hi1848 = vector.extract %m1845[0, 1] : f64 from vector<1x2xf64>
    %t1849 = arith.andi %li, %j3 : i32
    %t1850 = arith.muli %t1849, %j4 : i32
    %t1851 = arith.shrui %li, %j3 : i32
    %t1852 = arith.addi %t1850, %t1851 : i32
    %s1853, %unused1853 = gpu.shuffle idx %lo1847, %t1852, %j32 : f64
    %s1855, %unused1855 = gpu.shuffle idx %hi1848, %t1852, %j32 : f64
    %t1857 = arith.shrui %li, %j2 : i32
    %t1858 = arith.andi %t1857, %j1 : i32
    %t1859 = arith.cmpi eq, %t1858, %j1 : i32
    %sv1860 = arith.select %t1859, %s1855, %s1853 : f64
    %b1861 = vector.insert %sv1860, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m1862 = nvgpu.mma.sync(%f1846, %b1861, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1863 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo1864 = vector.extract %m1845[0, 0] : f64 from vector<1x2xf64>
    %hi1865 = vector.extract %m1845[0, 1] : f64 from vector<1x2xf64>
    %t1866 = arith.andi %li, %j3 : i32
    %t1867 = arith.muli %t1866, %j4 : i32
    %t1868 = arith.shrui %li, %j3 : i32
    %t1869 = arith.addi %t1867, %t1868 : i32
    %t1870 = arith.addi %t1869, %j16 : i32
    %s1871, %unused1871 = gpu.shuffle idx %lo1864, %t1870, %j32 : f64
    %s1873, %unused1873 = gpu.shuffle idx %hi1865, %t1870, %j32 : f64
    %t1875 = arith.shrui %li, %j2 : i32
    %t1876 = arith.andi %t1875, %j1 : i32
    %t1877 = arith.cmpi eq, %t1876, %j1 : i32
    %sv1878 = arith.select %t1877, %s1873, %s1871 : f64
    %b1879 = vector.insert %sv1878, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m1880 = nvgpu.mma.sync(%f1863, %b1879, %m1862) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %sv1881 = memref.subview %el2[0, 0, 0] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv1882 = memref.subview %el2[0, 0, 1] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %y1883 = vector.transfer_read %sv1881[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<1x2xf64>
    %y1884 = arith.subf %y1883, %m1880 : vector<1x2xf64>
    vector.transfer_write %y1884, %sv1881[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    %y1885 = vector.transfer_read %sv1882[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<1x2xf64>
    %y1886 = arith.addf %y1885, %m1880 : vector<1x2xf64>
    vector.transfer_write %y1886, %sv1882[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv1887 = memref.subview %v31786[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv1888 = memref.subview %v31787[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %f1889 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f1890 = vector.transfer_read %sv1887[%k, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m1891 = nvgpu.mma.sync(%f1889, %f1890, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1892 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f1893 = vector.transfer_read %sv1887[%k4, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m1894 = nvgpu.mma.sync(%f1892, %f1893, %m1891) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1895 = vector.transfer_read %sv1887[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f1896 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1897 = nvgpu.mma.sync(%f1895, %f1896, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1898 = vector.transfer_read %sv1887[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f1899 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1900 = nvgpu.mma.sync(%f1898, %f1899, %m1897) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %cf1901 = vector.transfer_read %sv1888[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<1x2xf64>
    %cf1902 = vector.transfer_read %el11[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf1903 = vector.transfer_read %el9[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf1904 = vector.transfer_read %el10[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %x1905 = arith.mulf %cf1902, %cf1901 : vector<1x2xf64>
    %x1906 = arith.mulf %cf1903, %m1900 : vector<1x2xf64>
    %x1907 = arith.mulf %cf1904, %m1894 : vector<1x2xf64>
    %x1908 = arith.addf %x1905, %x1906 : vector<1x2xf64>
    %fl1909 = arith.addf %x1908, %x1907 : vector<1x2xf64>
    %lo1910 = vector.extract %fl1909[0, 0] : f64 from vector<1x2xf64>
    %hi1911 = vector.extract %fl1909[0, 1] : f64 from vector<1x2xf64>
    %t1913 = arith.andi %li, %j28 : i32
    %t1914 = arith.andi %li, %j3 : i32
    %t1915 = arith.shrui %t1914, %j1 : i32
    %t1916 = arith.ori %t1913, %t1915 : i32
    %s1917, %unused1917 = gpu.shuffle idx %lo1910, %t1916, %j32 : f64
    %s1919, %unused1919 = gpu.shuffle idx %hi1911, %t1916, %j32 : f64
    %t1921 = arith.andi %li, %j1 : i32
    %t1922 = arith.cmpi eq, %t1921, %j1 : i32
    %sv1923 = arith.select %t1922, %s1919, %s1917 : f64
    %a1924 = vector.insert %sv1923, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f1925 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1926 = nvgpu.mma.sync(%a1924, %f1925, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %lo1927 = vector.extract %fl1909[0, 0] : f64 from vector<1x2xf64>
    %hi1928 = vector.extract %fl1909[0, 1] : f64 from vector<1x2xf64>
    %t1930 = arith.andi %li, %j28 : i32
    %t1931 = arith.andi %li, %j3 : i32
    %t1932 = arith.shrui %t1931, %j1 : i32
    %t1933 = arith.ori %t1930, %t1932 : i32
    %t1934 = arith.ori %t1933, %j2 : i32
    %s1935, %unused1935 = gpu.shuffle idx %lo1927, %t1934, %j32 : f64
    %s1937, %unused1937 = gpu.shuffle idx %hi1928, %t1934, %j32 : f64
    %t1939 = arith.andi %li, %j1 : i32
    %t1940 = arith.cmpi eq, %t1939, %j1 : i32
    %sv1941 = arith.select %t1940, %s1937, %s1935 : f64
    %a1942 = vector.insert %sv1941, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f1943 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1944 = nvgpu.mma.sync(%a1942, %f1943, %m1926) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1945 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo1946 = vector.extract %m1944[0, 0] : f64 from vector<1x2xf64>
    %hi1947 = vector.extract %m1944[0, 1] : f64 from vector<1x2xf64>
    %t1948 = arith.andi %li, %j3 : i32
    %t1949 = arith.muli %t1948, %j4 : i32
    %t1950 = arith.shrui %li, %j3 : i32
    %t1951 = arith.addi %t1949, %t1950 : i32
    %s1952, %unused1952 = gpu.shuffle idx %lo1946, %t1951, %j32 : f64
    %s1954, %unused1954 = gpu.shuffle idx %hi1947, %t1951, %j32 : f64
    %t1956 = arith.shrui %li, %j2 : i32
    %t1957 = arith.andi %t1956, %j1 : i32
    %t1958 = arith.cmpi eq, %t1957, %j1 : i32
    %sv1959 = arith.select %t1958, %s1954, %s1952 : f64
    %b1960 = vector.insert %sv1959, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m1961 = nvgpu.mma.sync(%f1945, %b1960, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1962 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo1963 = vector.extract %m1944[0, 0] : f64 from vector<1x2xf64>
    %hi1964 = vector.extract %m1944[0, 1] : f64 from vector<1x2xf64>
    %t1965 = arith.andi %li, %j3 : i32
    %t1966 = arith.muli %t1965, %j4 : i32
    %t1967 = arith.shrui %li, %j3 : i32
    %t1968 = arith.addi %t1966, %t1967 : i32
    %t1969 = arith.addi %t1968, %j16 : i32
    %s1970, %unused1970 = gpu.shuffle idx %lo1963, %t1969, %j32 : f64
    %s1972, %unused1972 = gpu.shuffle idx %hi1964, %t1969, %j32 : f64
    %t1974 = arith.shrui %li, %j2 : i32
    %t1975 = arith.andi %t1974, %j1 : i32
    %t1976 = arith.cmpi eq, %t1975, %j1 : i32
    %sv1977 = arith.select %t1976, %s1972, %s1970 : f64
    %b1978 = vector.insert %sv1977, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m1979 = nvgpu.mma.sync(%f1962, %b1978, %m1961) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %sv1980 = memref.subview %el2[0, 0, 1] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv1981 = memref.subview %el2[0, 0, 2] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %y1982 = vector.transfer_read %sv1980[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<1x2xf64>
    %y1983 = arith.subf %y1982, %m1979 : vector<1x2xf64>
    vector.transfer_write %y1983, %sv1980[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    %y1984 = vector.transfer_read %sv1981[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<1x2xf64>
    %y1985 = arith.addf %y1984, %m1979 : vector<1x2xf64>
    vector.transfer_write %y1985, %sv1981[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv1986 = memref.subview %v31786[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv1987 = memref.subview %v31787[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %f1988 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f1989 = vector.transfer_read %sv1986[%k, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m1990 = nvgpu.mma.sync(%f1988, %f1989, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1991 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f1992 = vector.transfer_read %sv1986[%k4, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m1993 = nvgpu.mma.sync(%f1991, %f1992, %m1990) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1994 = vector.transfer_read %sv1986[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f1995 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1996 = nvgpu.mma.sync(%f1994, %f1995, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f1997 = vector.transfer_read %sv1986[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f1998 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m1999 = nvgpu.mma.sync(%f1997, %f1998, %m1996) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %cf2000 = vector.transfer_read %sv1987[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<1x2xf64>
    %cf2001 = vector.transfer_read %el11[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf2002 = vector.transfer_read %el9[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf2003 = vector.transfer_read %el10[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %x2004 = arith.mulf %cf2001, %cf2000 : vector<1x2xf64>
    %x2005 = arith.mulf %cf2002, %m1999 : vector<1x2xf64>
    %x2006 = arith.mulf %cf2003, %m1993 : vector<1x2xf64>
    %x2007 = arith.addf %x2004, %x2005 : vector<1x2xf64>
    %fl2008 = arith.addf %x2007, %x2006 : vector<1x2xf64>
    %lo2009 = vector.extract %fl2008[0, 0] : f64 from vector<1x2xf64>
    %hi2010 = vector.extract %fl2008[0, 1] : f64 from vector<1x2xf64>
    %t2012 = arith.andi %li, %j28 : i32
    %t2013 = arith.andi %li, %j3 : i32
    %t2014 = arith.shrui %t2013, %j1 : i32
    %t2015 = arith.ori %t2012, %t2014 : i32
    %s2016, %unused2016 = gpu.shuffle idx %lo2009, %t2015, %j32 : f64
    %s2018, %unused2018 = gpu.shuffle idx %hi2010, %t2015, %j32 : f64
    %t2020 = arith.andi %li, %j1 : i32
    %t2021 = arith.cmpi eq, %t2020, %j1 : i32
    %sv2022 = arith.select %t2021, %s2018, %s2016 : f64
    %a2023 = vector.insert %sv2022, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f2024 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m2025 = nvgpu.mma.sync(%a2023, %f2024, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %lo2026 = vector.extract %fl2008[0, 0] : f64 from vector<1x2xf64>
    %hi2027 = vector.extract %fl2008[0, 1] : f64 from vector<1x2xf64>
    %t2029 = arith.andi %li, %j28 : i32
    %t2030 = arith.andi %li, %j3 : i32
    %t2031 = arith.shrui %t2030, %j1 : i32
    %t2032 = arith.ori %t2029, %t2031 : i32
    %t2033 = arith.ori %t2032, %j2 : i32
    %s2034, %unused2034 = gpu.shuffle idx %lo2026, %t2033, %j32 : f64
    %s2036, %unused2036 = gpu.shuffle idx %hi2027, %t2033, %j32 : f64
    %t2038 = arith.andi %li, %j1 : i32
    %t2039 = arith.cmpi eq, %t2038, %j1 : i32
    %sv2040 = arith.select %t2039, %s2036, %s2034 : f64
    %a2041 = vector.insert %sv2040, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f2042 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m2043 = nvgpu.mma.sync(%a2041, %f2042, %m2025) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f2044 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo2045 = vector.extract %m2043[0, 0] : f64 from vector<1x2xf64>
    %hi2046 = vector.extract %m2043[0, 1] : f64 from vector<1x2xf64>
    %t2047 = arith.andi %li, %j3 : i32
    %t2048 = arith.muli %t2047, %j4 : i32
    %t2049 = arith.shrui %li, %j3 : i32
    %t2050 = arith.addi %t2048, %t2049 : i32
    %s2051, %unused2051 = gpu.shuffle idx %lo2045, %t2050, %j32 : f64
    %s2053, %unused2053 = gpu.shuffle idx %hi2046, %t2050, %j32 : f64
    %t2055 = arith.shrui %li, %j2 : i32
    %t2056 = arith.andi %t2055, %j1 : i32
    %t2057 = arith.cmpi eq, %t2056, %j1 : i32
    %sv2058 = arith.select %t2057, %s2053, %s2051 : f64
    %b2059 = vector.insert %sv2058, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m2060 = nvgpu.mma.sync(%f2044, %b2059, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f2061 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo2062 = vector.extract %m2043[0, 0] : f64 from vector<1x2xf64>
    %hi2063 = vector.extract %m2043[0, 1] : f64 from vector<1x2xf64>
    %t2064 = arith.andi %li, %j3 : i32
    %t2065 = arith.muli %t2064, %j4 : i32
    %t2066 = arith.shrui %li, %j3 : i32
    %t2067 = arith.addi %t2065, %t2066 : i32
    %t2068 = arith.addi %t2067, %j16 : i32
    %s2069, %unused2069 = gpu.shuffle idx %lo2062, %t2068, %j32 : f64
    %s2071, %unused2071 = gpu.shuffle idx %hi2063, %t2068, %j32 : f64
    %t2073 = arith.shrui %li, %j2 : i32
    %t2074 = arith.andi %t2073, %j1 : i32
    %t2075 = arith.cmpi eq, %t2074, %j1 : i32
    %sv2076 = arith.select %t2075, %s2071, %s2069 : f64
    %b2077 = vector.insert %sv2076, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m2078 = nvgpu.mma.sync(%f2061, %b2077, %m2060) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %sv2079 = memref.subview %el2[0, 0, 2] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv2080 = memref.subview %el2[0, 0, 3] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %y2081 = vector.transfer_read %sv2079[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<1x2xf64>
    %y2082 = arith.subf %y2081, %m2078 : vector<1x2xf64>
    vector.transfer_write %y2082, %sv2079[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    %y2083 = vector.transfer_read %sv2080[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<1x2xf64>
    %y2084 = arith.addf %y2083, %m2078 : vector<1x2xf64>
    vector.transfer_write %y2084, %sv2080[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv2085 = memref.subview %v31786[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv2086 = memref.subview %v31787[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %f2087 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f2088 = vector.transfer_read %sv2085[%k, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m2089 = nvgpu.mma.sync(%f2087, %f2088, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f2090 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f2091 = vector.transfer_read %sv2085[%k4, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m2092 = nvgpu.mma.sync(%f2090, %f2091, %m2089) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f2093 = vector.transfer_read %sv2085[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f2094 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m2095 = nvgpu.mma.sync(%f2093, %f2094, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f2096 = vector.transfer_read %sv2085[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f2097 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m2098 = nvgpu.mma.sync(%f2096, %f2097, %m2095) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %cf2099 = vector.transfer_read %sv2086[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<1x2xf64>
    %cf2100 = vector.transfer_read %el11[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf2101 = vector.transfer_read %el9[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf2102 = vector.transfer_read %el10[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %x2103 = arith.mulf %cf2100, %cf2099 : vector<1x2xf64>
    %x2104 = arith.mulf %cf2101, %m2098 : vector<1x2xf64>
    %x2105 = arith.mulf %cf2102, %m2092 : vector<1x2xf64>
    %x2106 = arith.addf %x2103, %x2104 : vector<1x2xf64>
    %fl2107 = arith.addf %x2106, %x2105 : vector<1x2xf64>
    %lo2108 = vector.extract %fl2107[0, 0] : f64 from vector<1x2xf64>
    %hi2109 = vector.extract %fl2107[0, 1] : f64 from vector<1x2xf64>
    %t2111 = arith.andi %li, %j28 : i32
    %t2112 = arith.andi %li, %j3 : i32
    %t2113 = arith.shrui %t2112, %j1 : i32
    %t2114 = arith.ori %t2111, %t2113 : i32
    %s2115, %unused2115 = gpu.shuffle idx %lo2108, %t2114, %j32 : f64
    %s2117, %unused2117 = gpu.shuffle idx %hi2109, %t2114, %j32 : f64
    %t2119 = arith.andi %li, %j1 : i32
    %t2120 = arith.cmpi eq, %t2119, %j1 : i32
    %sv2121 = arith.select %t2120, %s2117, %s2115 : f64
    %a2122 = vector.insert %sv2121, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f2123 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m2124 = nvgpu.mma.sync(%a2122, %f2123, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %lo2125 = vector.extract %fl2107[0, 0] : f64 from vector<1x2xf64>
    %hi2126 = vector.extract %fl2107[0, 1] : f64 from vector<1x2xf64>
    %t2128 = arith.andi %li, %j28 : i32
    %t2129 = arith.andi %li, %j3 : i32
    %t2130 = arith.shrui %t2129, %j1 : i32
    %t2131 = arith.ori %t2128, %t2130 : i32
    %t2132 = arith.ori %t2131, %j2 : i32
    %s2133, %unused2133 = gpu.shuffle idx %lo2125, %t2132, %j32 : f64
    %s2135, %unused2135 = gpu.shuffle idx %hi2126, %t2132, %j32 : f64
    %t2137 = arith.andi %li, %j1 : i32
    %t2138 = arith.cmpi eq, %t2137, %j1 : i32
    %sv2139 = arith.select %t2138, %s2135, %s2133 : f64
    %a2140 = vector.insert %sv2139, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f2141 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m2142 = nvgpu.mma.sync(%a2140, %f2141, %m2124) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f2143 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo2144 = vector.extract %m2142[0, 0] : f64 from vector<1x2xf64>
    %hi2145 = vector.extract %m2142[0, 1] : f64 from vector<1x2xf64>
    %t2146 = arith.andi %li, %j3 : i32
    %t2147 = arith.muli %t2146, %j4 : i32
    %t2148 = arith.shrui %li, %j3 : i32
    %t2149 = arith.addi %t2147, %t2148 : i32
    %s2150, %unused2150 = gpu.shuffle idx %lo2144, %t2149, %j32 : f64
    %s2152, %unused2152 = gpu.shuffle idx %hi2145, %t2149, %j32 : f64
    %t2154 = arith.shrui %li, %j2 : i32
    %t2155 = arith.andi %t2154, %j1 : i32
    %t2156 = arith.cmpi eq, %t2155, %j1 : i32
    %sv2157 = arith.select %t2156, %s2152, %s2150 : f64
    %b2158 = vector.insert %sv2157, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m2159 = nvgpu.mma.sync(%f2143, %b2158, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f2160 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo2161 = vector.extract %m2142[0, 0] : f64 from vector<1x2xf64>
    %hi2162 = vector.extract %m2142[0, 1] : f64 from vector<1x2xf64>
    %t2163 = arith.andi %li, %j3 : i32
    %t2164 = arith.muli %t2163, %j4 : i32
    %t2165 = arith.shrui %li, %j3 : i32
    %t2166 = arith.addi %t2164, %t2165 : i32
    %t2167 = arith.addi %t2166, %j16 : i32
    %s2168, %unused2168 = gpu.shuffle idx %lo2161, %t2167, %j32 : f64
    %s2170, %unused2170 = gpu.shuffle idx %hi2162, %t2167, %j32 : f64
    %t2172 = arith.shrui %li, %j2 : i32
    %t2173 = arith.andi %t2172, %j1 : i32
    %t2174 = arith.cmpi eq, %t2173, %j1 : i32
    %sv2175 = arith.select %t2174, %s2170, %s2168 : f64
    %b2176 = vector.insert %sv2175, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m2177 = nvgpu.mma.sync(%f2160, %b2176, %m2159) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %sv2178 = memref.subview %el2[0, 0, 3] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv2179 = memref.subview %el2[0, 0, 4] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %y2180 = vector.transfer_read %sv2178[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<1x2xf64>
    %y2181 = arith.subf %y2180, %m2177 : vector<1x2xf64>
    vector.transfer_write %y2181, %sv2178[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    %y2182 = vector.transfer_read %sv2179[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<1x2xf64>
    %y2183 = arith.addf %y2182, %m2177 : vector<1x2xf64>
    vector.transfer_write %y2183, %sv2179[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv2184 = memref.subview %v31786[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv2185 = memref.subview %v31787[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %f2186 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f2187 = vector.transfer_read %sv2184[%k, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m2188 = nvgpu.mma.sync(%f2186, %f2187, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f2189 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f2190 = vector.transfer_read %sv2184[%k4, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m2191 = nvgpu.mma.sync(%f2189, %f2190, %m2188) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f2192 = vector.transfer_read %sv2184[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f2193 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m2194 = nvgpu.mma.sync(%f2192, %f2193, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f2195 = vector.transfer_read %sv2184[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f2196 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m2197 = nvgpu.mma.sync(%f2195, %f2196, %m2194) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %cf2198 = vector.transfer_read %sv2185[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<1x2xf64>
    %cf2199 = vector.transfer_read %el11[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf2200 = vector.transfer_read %el9[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf2201 = vector.transfer_read %el10[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %x2202 = arith.mulf %cf2199, %cf2198 : vector<1x2xf64>
    %x2203 = arith.mulf %cf2200, %m2197 : vector<1x2xf64>
    %x2204 = arith.mulf %cf2201, %m2191 : vector<1x2xf64>
    %x2205 = arith.addf %x2202, %x2203 : vector<1x2xf64>
    %fl2206 = arith.addf %x2205, %x2204 : vector<1x2xf64>
    %lo2207 = vector.extract %fl2206[0, 0] : f64 from vector<1x2xf64>
    %hi2208 = vector.extract %fl2206[0, 1] : f64 from vector<1x2xf64>
    %t2210 = arith.andi %li, %j28 : i32
    %t2211 = arith.andi %li, %j3 : i32
    %t2212 = arith.shrui %t2211, %j1 : i32
    %t2213 = arith.ori %t2210, %t2212 : i32
    %s2214, %unused2214 = gpu.shuffle idx %lo2207, %t2213, %j32 : f64
    %s2216, %unused2216 = gpu.shuffle idx %hi2208, %t2213, %j32 : f64
    %t2218 = arith.andi %li, %j1 : i32
    %t2219 = arith.cmpi eq, %t2218, %j1 : i32
    %sv2220 = arith.select %t2219, %s2216, %s2214 : f64
    %a2221 = vector.insert %sv2220, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f2222 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m2223 = nvgpu.mma.sync(%a2221, %f2222, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %lo2224 = vector.extract %fl2206[0, 0] : f64 from vector<1x2xf64>
    %hi2225 = vector.extract %fl2206[0, 1] : f64 from vector<1x2xf64>
    %t2227 = arith.andi %li, %j28 : i32
    %t2228 = arith.andi %li, %j3 : i32
    %t2229 = arith.shrui %t2228, %j1 : i32
    %t2230 = arith.ori %t2227, %t2229 : i32
    %t2231 = arith.ori %t2230, %j2 : i32
    %s2232, %unused2232 = gpu.shuffle idx %lo2224, %t2231, %j32 : f64
    %s2234, %unused2234 = gpu.shuffle idx %hi2225, %t2231, %j32 : f64
    %t2236 = arith.andi %li, %j1 : i32
    %t2237 = arith.cmpi eq, %t2236, %j1 : i32
    %sv2238 = arith.select %t2237, %s2234, %s2232 : f64
    %a2239 = vector.insert %sv2238, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f2240 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m2241 = nvgpu.mma.sync(%a2239, %f2240, %m2223) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f2242 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo2243 = vector.extract %m2241[0, 0] : f64 from vector<1x2xf64>
    %hi2244 = vector.extract %m2241[0, 1] : f64 from vector<1x2xf64>
    %t2245 = arith.andi %li, %j3 : i32
    %t2246 = arith.muli %t2245, %j4 : i32
    %t2247 = arith.shrui %li, %j3 : i32
    %t2248 = arith.addi %t2246, %t2247 : i32
    %s2249, %unused2249 = gpu.shuffle idx %lo2243, %t2248, %j32 : f64
    %s2251, %unused2251 = gpu.shuffle idx %hi2244, %t2248, %j32 : f64
    %t2253 = arith.shrui %li, %j2 : i32
    %t2254 = arith.andi %t2253, %j1 : i32
    %t2255 = arith.cmpi eq, %t2254, %j1 : i32
    %sv2256 = arith.select %t2255, %s2251, %s2249 : f64
    %b2257 = vector.insert %sv2256, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m2258 = nvgpu.mma.sync(%f2242, %b2257, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f2259 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo2260 = vector.extract %m2241[0, 0] : f64 from vector<1x2xf64>
    %hi2261 = vector.extract %m2241[0, 1] : f64 from vector<1x2xf64>
    %t2262 = arith.andi %li, %j3 : i32
    %t2263 = arith.muli %t2262, %j4 : i32
    %t2264 = arith.shrui %li, %j3 : i32
    %t2265 = arith.addi %t2263, %t2264 : i32
    %t2266 = arith.addi %t2265, %j16 : i32
    %s2267, %unused2267 = gpu.shuffle idx %lo2260, %t2266, %j32 : f64
    %s2269, %unused2269 = gpu.shuffle idx %hi2261, %t2266, %j32 : f64
    %t2271 = arith.shrui %li, %j2 : i32
    %t2272 = arith.andi %t2271, %j1 : i32
    %t2273 = arith.cmpi eq, %t2272, %j1 : i32
    %sv2274 = arith.select %t2273, %s2269, %s2267 : f64
    %b2275 = vector.insert %sv2274, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m2276 = nvgpu.mma.sync(%f2259, %b2275, %m2258) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %sv2277 = memref.subview %el2[0, 0, 4] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv2278 = memref.subview %el2[0, 0, 5] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %y2279 = vector.transfer_read %sv2277[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<1x2xf64>
    %y2280 = arith.subf %y2279, %m2276 : vector<1x2xf64>
    vector.transfer_write %y2280, %sv2277[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    %y2281 = vector.transfer_read %sv2278[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<1x2xf64>
    %y2282 = arith.addf %y2281, %m2276 : vector<1x2xf64>
    vector.transfer_write %y2282, %sv2278[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv2283 = memref.subview %v31786[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv2284 = memref.subview %v31787[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %f2285 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f2286 = vector.transfer_read %sv2283[%k, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m2287 = nvgpu.mma.sync(%f2285, %f2286, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f2288 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f2289 = vector.transfer_read %sv2283[%k4, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m2290 = nvgpu.mma.sync(%f2288, %f2289, %m2287) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f2291 = vector.transfer_read %sv2283[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f2292 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m2293 = nvgpu.mma.sync(%f2291, %f2292, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f2294 = vector.transfer_read %sv2283[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f2295 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m2296 = nvgpu.mma.sync(%f2294, %f2295, %m2293) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %cf2297 = vector.transfer_read %sv2284[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<1x2xf64>
    %cf2298 = vector.transfer_read %el11[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf2299 = vector.transfer_read %el9[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf2300 = vector.transfer_read %el10[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %x2301 = arith.mulf %cf2298, %cf2297 : vector<1x2xf64>
    %x2302 = arith.mulf %cf2299, %m2296 : vector<1x2xf64>
    %x2303 = arith.mulf %cf2300, %m2290 : vector<1x2xf64>
    %x2304 = arith.addf %x2301, %x2302 : vector<1x2xf64>
    %fl2305 = arith.addf %x2304, %x2303 : vector<1x2xf64>
    %lo2306 = vector.extract %fl2305[0, 0] : f64 from vector<1x2xf64>
    %hi2307 = vector.extract %fl2305[0, 1] : f64 from vector<1x2xf64>
    %t2309 = arith.andi %li, %j28 : i32
    %t2310 = arith.andi %li, %j3 : i32
    %t2311 = arith.shrui %t2310, %j1 : i32
    %t2312 = arith.ori %t2309, %t2311 : i32
    %s2313, %unused2313 = gpu.shuffle idx %lo2306, %t2312, %j32 : f64
    %s2315, %unused2315 = gpu.shuffle idx %hi2307, %t2312, %j32 : f64
    %t2317 = arith.andi %li, %j1 : i32
    %t2318 = arith.cmpi eq, %t2317, %j1 : i32
    %sv2319 = arith.select %t2318, %s2315, %s2313 : f64
    %a2320 = vector.insert %sv2319, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f2321 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m2322 = nvgpu.mma.sync(%a2320, %f2321, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %lo2323 = vector.extract %fl2305[0, 0] : f64 from vector<1x2xf64>
    %hi2324 = vector.extract %fl2305[0, 1] : f64 from vector<1x2xf64>
    %t2326 = arith.andi %li, %j28 : i32
    %t2327 = arith.andi %li, %j3 : i32
    %t2328 = arith.shrui %t2327, %j1 : i32
    %t2329 = arith.ori %t2326, %t2328 : i32
    %t2330 = arith.ori %t2329, %j2 : i32
    %s2331, %unused2331 = gpu.shuffle idx %lo2323, %t2330, %j32 : f64
    %s2333, %unused2333 = gpu.shuffle idx %hi2324, %t2330, %j32 : f64
    %t2335 = arith.andi %li, %j1 : i32
    %t2336 = arith.cmpi eq, %t2335, %j1 : i32
    %sv2337 = arith.select %t2336, %s2333, %s2331 : f64
    %a2338 = vector.insert %sv2337, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f2339 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m2340 = nvgpu.mma.sync(%a2338, %f2339, %m2322) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f2341 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo2342 = vector.extract %m2340[0, 0] : f64 from vector<1x2xf64>
    %hi2343 = vector.extract %m2340[0, 1] : f64 from vector<1x2xf64>
    %t2344 = arith.andi %li, %j3 : i32
    %t2345 = arith.muli %t2344, %j4 : i32
    %t2346 = arith.shrui %li, %j3 : i32
    %t2347 = arith.addi %t2345, %t2346 : i32
    %s2348, %unused2348 = gpu.shuffle idx %lo2342, %t2347, %j32 : f64
    %s2350, %unused2350 = gpu.shuffle idx %hi2343, %t2347, %j32 : f64
    %t2352 = arith.shrui %li, %j2 : i32
    %t2353 = arith.andi %t2352, %j1 : i32
    %t2354 = arith.cmpi eq, %t2353, %j1 : i32
    %sv2355 = arith.select %t2354, %s2350, %s2348 : f64
    %b2356 = vector.insert %sv2355, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m2357 = nvgpu.mma.sync(%f2341, %b2356, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f2358 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo2359 = vector.extract %m2340[0, 0] : f64 from vector<1x2xf64>
    %hi2360 = vector.extract %m2340[0, 1] : f64 from vector<1x2xf64>
    %t2361 = arith.andi %li, %j3 : i32
    %t2362 = arith.muli %t2361, %j4 : i32
    %t2363 = arith.shrui %li, %j3 : i32
    %t2364 = arith.addi %t2362, %t2363 : i32
    %t2365 = arith.addi %t2364, %j16 : i32
    %s2366, %unused2366 = gpu.shuffle idx %lo2359, %t2365, %j32 : f64
    %s2368, %unused2368 = gpu.shuffle idx %hi2360, %t2365, %j32 : f64
    %t2370 = arith.shrui %li, %j2 : i32
    %t2371 = arith.andi %t2370, %j1 : i32
    %t2372 = arith.cmpi eq, %t2371, %j1 : i32
    %sv2373 = arith.select %t2372, %s2368, %s2366 : f64
    %b2374 = vector.insert %sv2373, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m2375 = nvgpu.mma.sync(%f2358, %b2374, %m2357) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %sv2376 = memref.subview %el2[0, 0, 5] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv2377 = memref.subview %el2[0, 0, 6] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %y2378 = vector.transfer_read %sv2376[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<1x2xf64>
    %y2379 = arith.subf %y2378, %m2375 : vector<1x2xf64>
    vector.transfer_write %y2379, %sv2376[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    %y2380 = vector.transfer_read %sv2377[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<1x2xf64>
    %y2381 = arith.addf %y2380, %m2375 : vector<1x2xf64>
    vector.transfer_write %y2381, %sv2377[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv2382 = memref.subview %v31786[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv2383 = memref.subview %v31787[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %f2384 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f2385 = vector.transfer_read %sv2382[%k, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m2386 = nvgpu.mma.sync(%f2384, %f2385, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f2387 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f2388 = vector.transfer_read %sv2382[%k4, %i], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %m2389 = nvgpu.mma.sync(%f2387, %f2388, %m2386) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f2390 = vector.transfer_read %sv2382[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f2391 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m2392 = nvgpu.mma.sync(%f2390, %f2391, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f2393 = vector.transfer_read %sv2382[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<1x1xf64>
    %f2394 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m2395 = nvgpu.mma.sync(%f2393, %f2394, %m2392) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %cf2396 = vector.transfer_read %sv2383[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<1x2xf64>
    %cf2397 = vector.transfer_read %el11[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf2398 = vector.transfer_read %el9[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %cf2399 = vector.transfer_read %el10[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
    %x2400 = arith.mulf %cf2397, %cf2396 : vector<1x2xf64>
    %x2401 = arith.mulf %cf2398, %m2395 : vector<1x2xf64>
    %x2402 = arith.mulf %cf2399, %m2389 : vector<1x2xf64>
    %x2403 = arith.addf %x2400, %x2401 : vector<1x2xf64>
    %fl2404 = arith.addf %x2403, %x2402 : vector<1x2xf64>
    %lo2405 = vector.extract %fl2404[0, 0] : f64 from vector<1x2xf64>
    %hi2406 = vector.extract %fl2404[0, 1] : f64 from vector<1x2xf64>
    %t2408 = arith.andi %li, %j28 : i32
    %t2409 = arith.andi %li, %j3 : i32
    %t2410 = arith.shrui %t2409, %j1 : i32
    %t2411 = arith.ori %t2408, %t2410 : i32
    %s2412, %unused2412 = gpu.shuffle idx %lo2405, %t2411, %j32 : f64
    %s2414, %unused2414 = gpu.shuffle idx %hi2406, %t2411, %j32 : f64
    %t2416 = arith.andi %li, %j1 : i32
    %t2417 = arith.cmpi eq, %t2416, %j1 : i32
    %sv2418 = arith.select %t2417, %s2414, %s2412 : f64
    %a2419 = vector.insert %sv2418, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f2420 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m2421 = nvgpu.mma.sync(%a2419, %f2420, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %lo2422 = vector.extract %fl2404[0, 0] : f64 from vector<1x2xf64>
    %hi2423 = vector.extract %fl2404[0, 1] : f64 from vector<1x2xf64>
    %t2425 = arith.andi %li, %j28 : i32
    %t2426 = arith.andi %li, %j3 : i32
    %t2427 = arith.shrui %t2426, %j1 : i32
    %t2428 = arith.ori %t2425, %t2427 : i32
    %t2429 = arith.ori %t2428, %j2 : i32
    %s2430, %unused2430 = gpu.shuffle idx %lo2422, %t2429, %j32 : f64
    %s2432, %unused2432 = gpu.shuffle idx %hi2423, %t2429, %j32 : f64
    %t2434 = arith.andi %li, %j1 : i32
    %t2435 = arith.cmpi eq, %t2434, %j1 : i32
    %sv2436 = arith.select %t2435, %s2432, %s2430 : f64
    %a2437 = vector.insert %sv2436, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f2438 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m2439 = nvgpu.mma.sync(%a2437, %f2438, %m2421) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f2440 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo2441 = vector.extract %m2439[0, 0] : f64 from vector<1x2xf64>
    %hi2442 = vector.extract %m2439[0, 1] : f64 from vector<1x2xf64>
    %t2443 = arith.andi %li, %j3 : i32
    %t2444 = arith.muli %t2443, %j4 : i32
    %t2445 = arith.shrui %li, %j3 : i32
    %t2446 = arith.addi %t2444, %t2445 : i32
    %s2447, %unused2447 = gpu.shuffle idx %lo2441, %t2446, %j32 : f64
    %s2449, %unused2449 = gpu.shuffle idx %hi2442, %t2446, %j32 : f64
    %t2451 = arith.shrui %li, %j2 : i32
    %t2452 = arith.andi %t2451, %j1 : i32
    %t2453 = arith.cmpi eq, %t2452, %j1 : i32
    %sv2454 = arith.select %t2453, %s2449, %s2447 : f64
    %b2455 = vector.insert %sv2454, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m2456 = nvgpu.mma.sync(%f2440, %b2455, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f2457 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo2458 = vector.extract %m2439[0, 0] : f64 from vector<1x2xf64>
    %hi2459 = vector.extract %m2439[0, 1] : f64 from vector<1x2xf64>
    %t2460 = arith.andi %li, %j3 : i32
    %t2461 = arith.muli %t2460, %j4 : i32
    %t2462 = arith.shrui %li, %j3 : i32
    %t2463 = arith.addi %t2461, %t2462 : i32
    %t2464 = arith.addi %t2463, %j16 : i32
    %s2465, %unused2465 = gpu.shuffle idx %lo2458, %t2464, %j32 : f64
    %s2467, %unused2467 = gpu.shuffle idx %hi2459, %t2464, %j32 : f64
    %t2469 = arith.shrui %li, %j2 : i32
    %t2470 = arith.andi %t2469, %j1 : i32
    %t2471 = arith.cmpi eq, %t2470, %j1 : i32
    %sv2472 = arith.select %t2471, %s2467, %s2465 : f64
    %b2473 = vector.insert %sv2472, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m2474 = nvgpu.mma.sync(%f2457, %b2473, %m2456) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %sv2475 = memref.subview %el2[0, 0, 6] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv2476 = memref.subview %el2[0, 0, 7] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %y2477 = vector.transfer_read %sv2475[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<1x2xf64>
    %y2478 = arith.subf %y2477, %m2474 : vector<1x2xf64>
    vector.transfer_write %y2478, %sv2475[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    %y2479 = vector.transfer_read %sv2476[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<1x2xf64>
    %y2480 = arith.addf %y2479, %m2474 : vector<1x2xf64>
    vector.transfer_write %y2480, %sv2476[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    gpu.barrier
    gpu.return
  }
}
