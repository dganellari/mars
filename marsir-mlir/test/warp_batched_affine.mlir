#a = affine_map<(m, n, k) -> (m, k)>
#b = affine_map<(m, n, k) -> (k, n)>
#bt = affine_map<(m, n, k) -> (n, k)>
#c = affine_map<(m, n, k) -> (m, n)>
gpu.module @mir_kernels [#nvvm.target<chip = "sm_90", O = 3>] {
  gpu.func @full_batched_affine(%U: memref<?x8x8x8xf64>, %Btil: memref<8x8xf64>, %Dtil: memref<8x8xf64>, %Dm: memref<8x8xf64>, %W: memref<8x8xf64>, %g00all: memref<?x8x8xf64>, %g01all: memref<?x8x8xf64>, %g02all: memref<?x8x8xf64>, %g10all: memref<?x8x8xf64>, %g11all: memref<?x8x8xf64>, %g12all: memref<?x8x8xf64>, %g20all: memref<?x8x8xf64>, %g21all: memref<?x8x8xf64>, %g22all: memref<?x8x8xf64>, %Y: memref<?x8x8x8xf64>) workgroup(%interp_all: memref<8x64xf64, #gpu.address_space<workgroup>>, %deriv_all: memref<8x64xf64, #gpu.address_space<workgroup>>, %flux: memref<8x8xf64, #gpu.address_space<workgroup>>, %tmp: memref<8x8xf64, #gpu.address_space<workgroup>>, %intf: memref<8x8xf64, #gpu.address_space<workgroup>>) kernel {
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
      %dt152 = arith.mulf %r150, %g151 : vector<8x8xf64>
      %acc153 = arith.constant dense<0.0> : vector<8x8xf64>
      %a154 = vector.transfer_read %sv142[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b155 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r156 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a154, %b155, %acc153 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a157 = vector.transfer_read %sv142[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b158 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r159 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a157, %b158, %r156 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g160 = vector.transfer_read %el3[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt161 = arith.mulf %r159, %g160 : vector<8x8xf64>
      %dv162 = vector.transfer_read %sv143[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g163 = vector.transfer_read %el5[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t164 = arith.mulf %dv162, %g163 : vector<8x8xf64>
      %t2165 = arith.addf %t164, %dt161 : vector<8x8xf64>
      %fl166 = arith.addf %t2165, %dt152 : vector<8x8xf64>
      vector.transfer_write %fl166, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc167 = arith.constant dense<0.0> : vector<8x8xf64>
      %a168 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b169 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r170 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a168, %b169, %acc167 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a171 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b172 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r173 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a171, %b172, %r170 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r173, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc174 = arith.constant dense<0.0> : vector<8x8xf64>
      %a175 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b176 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r177 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a175, %b176, %acc174 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a178 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b179 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r180 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a178, %b179, %r177 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r180, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv181 = memref.subview %el2[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv182 = memref.subview %el2[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa183 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r184 = vector.transfer_read %sv181[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %m185 = arith.subf %r184, %fa183 : vector<8x8xf64>
      vector.transfer_write %m185, %sv181[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      %r186 = vector.transfer_read %sv182[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %p187 = arith.addf %r186, %fa183 : vector<8x8xf64>
      vector.transfer_write %p187, %sv182[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    }
    gpu.barrier
    %sv188 = memref.subview %v3140[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv189 = memref.subview %v3141[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc190 = arith.constant dense<0.0> : vector<8x8xf64>
      %a191 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b192 = vector.transfer_read %sv188[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r193 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a191, %b192, %acc190 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a194 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b195 = vector.transfer_read %sv188[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r196 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a194, %b195, %r193 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g197 = vector.transfer_read %el4[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt198 = arith.mulf %r196, %g197 : vector<8x8xf64>
      %acc199 = arith.constant dense<0.0> : vector<8x8xf64>
      %a200 = vector.transfer_read %sv188[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b201 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r202 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a200, %b201, %acc199 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a203 = vector.transfer_read %sv188[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b204 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r205 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a203, %b204, %r202 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g206 = vector.transfer_read %el3[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt207 = arith.mulf %r205, %g206 : vector<8x8xf64>
      %dv208 = vector.transfer_read %sv189[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g209 = vector.transfer_read %el5[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t210 = arith.mulf %dv208, %g209 : vector<8x8xf64>
      %t2211 = arith.addf %t210, %dt207 : vector<8x8xf64>
      %fl212 = arith.addf %t2211, %dt198 : vector<8x8xf64>
      vector.transfer_write %fl212, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc213 = arith.constant dense<0.0> : vector<8x8xf64>
      %a214 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b215 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r216 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a214, %b215, %acc213 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a217 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b218 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r219 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a217, %b218, %r216 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r219, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc220 = arith.constant dense<0.0> : vector<8x8xf64>
      %a221 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b222 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r223 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a221, %b222, %acc220 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a224 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b225 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r226 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a224, %b225, %r223 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r226, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv227 = memref.subview %el2[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv228 = memref.subview %el2[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa229 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r230 = vector.transfer_read %sv227[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %m231 = arith.subf %r230, %fa229 : vector<8x8xf64>
      vector.transfer_write %m231, %sv227[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      %r232 = vector.transfer_read %sv228[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %p233 = arith.addf %r232, %fa229 : vector<8x8xf64>
      vector.transfer_write %p233, %sv228[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    }
    gpu.barrier
    %sv234 = memref.subview %v3140[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv235 = memref.subview %v3141[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc236 = arith.constant dense<0.0> : vector<8x8xf64>
      %a237 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b238 = vector.transfer_read %sv234[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r239 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a237, %b238, %acc236 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a240 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b241 = vector.transfer_read %sv234[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r242 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a240, %b241, %r239 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g243 = vector.transfer_read %el4[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt244 = arith.mulf %r242, %g243 : vector<8x8xf64>
      %acc245 = arith.constant dense<0.0> : vector<8x8xf64>
      %a246 = vector.transfer_read %sv234[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b247 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r248 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a246, %b247, %acc245 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a249 = vector.transfer_read %sv234[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b250 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r251 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a249, %b250, %r248 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g252 = vector.transfer_read %el3[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt253 = arith.mulf %r251, %g252 : vector<8x8xf64>
      %dv254 = vector.transfer_read %sv235[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g255 = vector.transfer_read %el5[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t256 = arith.mulf %dv254, %g255 : vector<8x8xf64>
      %t2257 = arith.addf %t256, %dt253 : vector<8x8xf64>
      %fl258 = arith.addf %t2257, %dt244 : vector<8x8xf64>
      vector.transfer_write %fl258, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc259 = arith.constant dense<0.0> : vector<8x8xf64>
      %a260 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b261 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r262 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a260, %b261, %acc259 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a263 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b264 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r265 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a263, %b264, %r262 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r265, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc266 = arith.constant dense<0.0> : vector<8x8xf64>
      %a267 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b268 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r269 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a267, %b268, %acc266 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a270 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b271 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r272 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a270, %b271, %r269 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r272, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv273 = memref.subview %el2[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv274 = memref.subview %el2[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa275 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r276 = vector.transfer_read %sv273[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %m277 = arith.subf %r276, %fa275 : vector<8x8xf64>
      vector.transfer_write %m277, %sv273[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      %r278 = vector.transfer_read %sv274[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %p279 = arith.addf %r278, %fa275 : vector<8x8xf64>
      vector.transfer_write %p279, %sv274[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    }
    gpu.barrier
    %sv280 = memref.subview %v3140[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv281 = memref.subview %v3141[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc282 = arith.constant dense<0.0> : vector<8x8xf64>
      %a283 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b284 = vector.transfer_read %sv280[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r285 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a283, %b284, %acc282 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a286 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b287 = vector.transfer_read %sv280[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r288 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a286, %b287, %r285 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g289 = vector.transfer_read %el4[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt290 = arith.mulf %r288, %g289 : vector<8x8xf64>
      %acc291 = arith.constant dense<0.0> : vector<8x8xf64>
      %a292 = vector.transfer_read %sv280[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b293 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r294 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a292, %b293, %acc291 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a295 = vector.transfer_read %sv280[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b296 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r297 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a295, %b296, %r294 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g298 = vector.transfer_read %el3[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt299 = arith.mulf %r297, %g298 : vector<8x8xf64>
      %dv300 = vector.transfer_read %sv281[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g301 = vector.transfer_read %el5[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t302 = arith.mulf %dv300, %g301 : vector<8x8xf64>
      %t2303 = arith.addf %t302, %dt299 : vector<8x8xf64>
      %fl304 = arith.addf %t2303, %dt290 : vector<8x8xf64>
      vector.transfer_write %fl304, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc305 = arith.constant dense<0.0> : vector<8x8xf64>
      %a306 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b307 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r308 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a306, %b307, %acc305 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a309 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b310 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r311 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a309, %b310, %r308 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r311, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc312 = arith.constant dense<0.0> : vector<8x8xf64>
      %a313 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b314 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r315 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a313, %b314, %acc312 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a316 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b317 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r318 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a316, %b317, %r315 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r318, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv319 = memref.subview %el2[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv320 = memref.subview %el2[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa321 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r322 = vector.transfer_read %sv319[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %m323 = arith.subf %r322, %fa321 : vector<8x8xf64>
      vector.transfer_write %m323, %sv319[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      %r324 = vector.transfer_read %sv320[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %p325 = arith.addf %r324, %fa321 : vector<8x8xf64>
      vector.transfer_write %p325, %sv320[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    }
    gpu.barrier
    %sv326 = memref.subview %v3140[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv327 = memref.subview %v3141[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc328 = arith.constant dense<0.0> : vector<8x8xf64>
      %a329 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b330 = vector.transfer_read %sv326[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r331 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a329, %b330, %acc328 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a332 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b333 = vector.transfer_read %sv326[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r334 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a332, %b333, %r331 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g335 = vector.transfer_read %el4[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt336 = arith.mulf %r334, %g335 : vector<8x8xf64>
      %acc337 = arith.constant dense<0.0> : vector<8x8xf64>
      %a338 = vector.transfer_read %sv326[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b339 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r340 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a338, %b339, %acc337 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a341 = vector.transfer_read %sv326[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b342 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r343 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a341, %b342, %r340 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g344 = vector.transfer_read %el3[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt345 = arith.mulf %r343, %g344 : vector<8x8xf64>
      %dv346 = vector.transfer_read %sv327[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g347 = vector.transfer_read %el5[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t348 = arith.mulf %dv346, %g347 : vector<8x8xf64>
      %t2349 = arith.addf %t348, %dt345 : vector<8x8xf64>
      %fl350 = arith.addf %t2349, %dt336 : vector<8x8xf64>
      vector.transfer_write %fl350, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc351 = arith.constant dense<0.0> : vector<8x8xf64>
      %a352 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b353 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r354 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a352, %b353, %acc351 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a355 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b356 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r357 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a355, %b356, %r354 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r357, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc358 = arith.constant dense<0.0> : vector<8x8xf64>
      %a359 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b360 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r361 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a359, %b360, %acc358 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a362 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b363 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r364 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a362, %b363, %r361 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r364, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv365 = memref.subview %el2[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv366 = memref.subview %el2[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa367 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r368 = vector.transfer_read %sv365[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %m369 = arith.subf %r368, %fa367 : vector<8x8xf64>
      vector.transfer_write %m369, %sv365[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      %r370 = vector.transfer_read %sv366[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %p371 = arith.addf %r370, %fa367 : vector<8x8xf64>
      vector.transfer_write %p371, %sv366[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    }
    gpu.barrier
    %sv372 = memref.subview %v3140[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv373 = memref.subview %v3141[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc374 = arith.constant dense<0.0> : vector<8x8xf64>
      %a375 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b376 = vector.transfer_read %sv372[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r377 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a375, %b376, %acc374 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a378 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b379 = vector.transfer_read %sv372[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r380 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a378, %b379, %r377 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g381 = vector.transfer_read %el4[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt382 = arith.mulf %r380, %g381 : vector<8x8xf64>
      %acc383 = arith.constant dense<0.0> : vector<8x8xf64>
      %a384 = vector.transfer_read %sv372[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b385 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r386 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a384, %b385, %acc383 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a387 = vector.transfer_read %sv372[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b388 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r389 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a387, %b388, %r386 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g390 = vector.transfer_read %el3[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt391 = arith.mulf %r389, %g390 : vector<8x8xf64>
      %dv392 = vector.transfer_read %sv373[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g393 = vector.transfer_read %el5[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t394 = arith.mulf %dv392, %g393 : vector<8x8xf64>
      %t2395 = arith.addf %t394, %dt391 : vector<8x8xf64>
      %fl396 = arith.addf %t2395, %dt382 : vector<8x8xf64>
      vector.transfer_write %fl396, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc397 = arith.constant dense<0.0> : vector<8x8xf64>
      %a398 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b399 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r400 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a398, %b399, %acc397 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a401 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b402 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r403 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a401, %b402, %r400 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r403, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc404 = arith.constant dense<0.0> : vector<8x8xf64>
      %a405 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b406 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r407 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a405, %b406, %acc404 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a408 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b409 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r410 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a408, %b409, %r407 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r410, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv411 = memref.subview %el2[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv412 = memref.subview %el2[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa413 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r414 = vector.transfer_read %sv411[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %m415 = arith.subf %r414, %fa413 : vector<8x8xf64>
      vector.transfer_write %m415, %sv411[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      %r416 = vector.transfer_read %sv412[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %p417 = arith.addf %r416, %fa413 : vector<8x8xf64>
      vector.transfer_write %p417, %sv412[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    }
    gpu.barrier
    %sv418 = memref.subview %v3140[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv419 = memref.subview %v3141[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc420 = arith.constant dense<0.0> : vector<8x8xf64>
      %a421 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b422 = vector.transfer_read %sv418[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r423 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a421, %b422, %acc420 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a424 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b425 = vector.transfer_read %sv418[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r426 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a424, %b425, %r423 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g427 = vector.transfer_read %el4[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt428 = arith.mulf %r426, %g427 : vector<8x8xf64>
      %acc429 = arith.constant dense<0.0> : vector<8x8xf64>
      %a430 = vector.transfer_read %sv418[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b431 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r432 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a430, %b431, %acc429 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a433 = vector.transfer_read %sv418[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b434 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r435 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a433, %b434, %r432 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g436 = vector.transfer_read %el3[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt437 = arith.mulf %r435, %g436 : vector<8x8xf64>
      %dv438 = vector.transfer_read %sv419[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g439 = vector.transfer_read %el5[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t440 = arith.mulf %dv438, %g439 : vector<8x8xf64>
      %t2441 = arith.addf %t440, %dt437 : vector<8x8xf64>
      %fl442 = arith.addf %t2441, %dt428 : vector<8x8xf64>
      vector.transfer_write %fl442, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc443 = arith.constant dense<0.0> : vector<8x8xf64>
      %a444 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b445 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r446 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a444, %b445, %acc443 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a447 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b448 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r449 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a447, %b448, %r446 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r449, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc450 = arith.constant dense<0.0> : vector<8x8xf64>
      %a451 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b452 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r453 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a451, %b452, %acc450 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a454 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b455 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r456 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a454, %b455, %r453 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r456, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv457 = memref.subview %el2[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv458 = memref.subview %el2[7, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa459 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r460 = vector.transfer_read %sv457[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %m461 = arith.subf %r460, %fa459 : vector<8x8xf64>
      vector.transfer_write %m461, %sv457[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      %r462 = vector.transfer_read %sv458[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %p463 = arith.addf %r462, %fa459 : vector<8x8xf64>
      vector.transfer_write %p463, %sv458[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    }
    gpu.barrier
    %sv464 = memref.subview %el1[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc465 = arith.constant dense<0.0> : vector<8x8xf64>
      %a466 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b467 = vector.transfer_read %sv464[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r468 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a466, %b467, %acc465 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a469 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b470 = vector.transfer_read %sv464[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r471 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a469, %b470, %r468 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r471, %interp_all[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv472 = memref.subview %el1[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc473 = arith.constant dense<0.0> : vector<8x8xf64>
      %a474 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b475 = vector.transfer_read %sv472[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r476 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a474, %b475, %acc473 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a477 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b478 = vector.transfer_read %sv472[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r479 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a477, %b478, %r476 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r479, %interp_all[%c0, %c8] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv480 = memref.subview %el1[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc481 = arith.constant dense<0.0> : vector<8x8xf64>
      %a482 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b483 = vector.transfer_read %sv480[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r484 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a482, %b483, %acc481 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a485 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b486 = vector.transfer_read %sv480[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r487 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a485, %b486, %r484 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r487, %interp_all[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv488 = memref.subview %el1[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc489 = arith.constant dense<0.0> : vector<8x8xf64>
      %a490 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b491 = vector.transfer_read %sv488[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r492 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a490, %b491, %acc489 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a493 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b494 = vector.transfer_read %sv488[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r495 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a493, %b494, %r492 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r495, %interp_all[%c0, %c24] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv496 = memref.subview %el1[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc497 = arith.constant dense<0.0> : vector<8x8xf64>
      %a498 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b499 = vector.transfer_read %sv496[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r500 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a498, %b499, %acc497 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a501 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b502 = vector.transfer_read %sv496[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r503 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a501, %b502, %r500 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r503, %interp_all[%c0, %c32] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv504 = memref.subview %el1[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc505 = arith.constant dense<0.0> : vector<8x8xf64>
      %a506 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b507 = vector.transfer_read %sv504[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r508 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a506, %b507, %acc505 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a509 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b510 = vector.transfer_read %sv504[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r511 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a509, %b510, %r508 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r511, %interp_all[%c0, %c40] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv512 = memref.subview %el1[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc513 = arith.constant dense<0.0> : vector<8x8xf64>
      %a514 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b515 = vector.transfer_read %sv512[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r516 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a514, %b515, %acc513 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a517 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b518 = vector.transfer_read %sv512[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r519 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a517, %b518, %r516 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r519, %interp_all[%c0, %c48] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv520 = memref.subview %el1[7, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc521 = arith.constant dense<0.0> : vector<8x8xf64>
      %a522 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b523 = vector.transfer_read %sv520[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r524 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a522, %b523, %acc521 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a525 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b526 = vector.transfer_read %sv520[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r527 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a525, %b526, %r524 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r527, %interp_all[%c0, %c56] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv528 = memref.subview %el1[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc529 = arith.constant dense<0.0> : vector<8x8xf64>
      %a530 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b531 = vector.transfer_read %sv528[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r532 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a530, %b531, %acc529 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a533 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b534 = vector.transfer_read %sv528[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r535 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a533, %b534, %r532 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r535, %deriv_all[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv536 = memref.subview %el1[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc537 = arith.constant dense<0.0> : vector<8x8xf64>
      %a538 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b539 = vector.transfer_read %sv536[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r540 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a538, %b539, %acc537 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a541 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b542 = vector.transfer_read %sv536[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r543 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a541, %b542, %r540 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r543, %deriv_all[%c0, %c8] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv544 = memref.subview %el1[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc545 = arith.constant dense<0.0> : vector<8x8xf64>
      %a546 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b547 = vector.transfer_read %sv544[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r548 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a546, %b547, %acc545 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a549 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b550 = vector.transfer_read %sv544[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r551 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a549, %b550, %r548 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r551, %deriv_all[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv552 = memref.subview %el1[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc553 = arith.constant dense<0.0> : vector<8x8xf64>
      %a554 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b555 = vector.transfer_read %sv552[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r556 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a554, %b555, %acc553 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a557 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b558 = vector.transfer_read %sv552[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r559 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a557, %b558, %r556 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r559, %deriv_all[%c0, %c24] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv560 = memref.subview %el1[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc561 = arith.constant dense<0.0> : vector<8x8xf64>
      %a562 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b563 = vector.transfer_read %sv560[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r564 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a562, %b563, %acc561 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a565 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b566 = vector.transfer_read %sv560[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r567 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a565, %b566, %r564 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r567, %deriv_all[%c0, %c32] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv568 = memref.subview %el1[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc569 = arith.constant dense<0.0> : vector<8x8xf64>
      %a570 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b571 = vector.transfer_read %sv568[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r572 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a570, %b571, %acc569 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a573 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b574 = vector.transfer_read %sv568[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r575 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a573, %b574, %r572 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r575, %deriv_all[%c0, %c40] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv576 = memref.subview %el1[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc577 = arith.constant dense<0.0> : vector<8x8xf64>
      %a578 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b579 = vector.transfer_read %sv576[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r580 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a578, %b579, %acc577 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a581 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b582 = vector.transfer_read %sv576[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r583 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a581, %b582, %r580 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r583, %deriv_all[%c0, %c48] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv584 = memref.subview %el1[7, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc585 = arith.constant dense<0.0> : vector<8x8xf64>
      %a586 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b587 = vector.transfer_read %sv584[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r588 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a586, %b587, %acc585 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a589 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b590 = vector.transfer_read %sv584[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r591 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a589, %b590, %r588 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r591, %deriv_all[%c0, %c56] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %v3592 = memref.reinterpret_cast %interp_all to offset: [0], sizes: [8, 8, 8], strides: [64, 8, 1] : memref<8x64xf64, #gpu.address_space<workgroup>> to memref<8x8x8xf64, #gpu.address_space<workgroup>>
    %v3593 = memref.reinterpret_cast %deriv_all to offset: [0], sizes: [8, 8, 8], strides: [64, 8, 1] : memref<8x64xf64, #gpu.address_space<workgroup>> to memref<8x8x8xf64, #gpu.address_space<workgroup>>
    %sv594 = memref.subview %v3592[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>
    %sv595 = memref.subview %v3593[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc596 = arith.constant dense<0.0> : vector<8x8xf64>
      %a597 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b598 = vector.transfer_read %sv594[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r599 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a597, %b598, %acc596 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a600 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b601 = vector.transfer_read %sv594[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r602 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a600, %b601, %r599 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g603 = vector.transfer_read %el7[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt604 = arith.mulf %r602, %g603 : vector<8x8xf64>
      %acc605 = arith.constant dense<0.0> : vector<8x8xf64>
      %a606 = vector.transfer_read %sv594[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b607 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r608 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a606, %b607, %acc605 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a609 = vector.transfer_read %sv594[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b610 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r611 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a609, %b610, %r608 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g612 = vector.transfer_read %el6[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt613 = arith.mulf %r611, %g612 : vector<8x8xf64>
      %dv614 = vector.transfer_read %sv595[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g615 = vector.transfer_read %el8[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t616 = arith.mulf %dv614, %g615 : vector<8x8xf64>
      %t2617 = arith.addf %t616, %dt613 : vector<8x8xf64>
      %fl618 = arith.addf %t2617, %dt604 : vector<8x8xf64>
      vector.transfer_write %fl618, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc619 = arith.constant dense<0.0> : vector<8x8xf64>
      %a620 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b621 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r622 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a620, %b621, %acc619 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a623 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b624 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r625 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a623, %b624, %r622 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r625, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc626 = arith.constant dense<0.0> : vector<8x8xf64>
      %a627 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b628 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r629 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a627, %b628, %acc626 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a630 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b631 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r632 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a630, %b631, %r629 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r632, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv633 = memref.subview %el2[0, 0, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv634 = memref.subview %el2[0, 1, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa635 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r636 = vector.transfer_read %sv633[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %m637 = arith.subf %r636, %fa635 : vector<8x8xf64>
      vector.transfer_write %m637, %sv633[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
      %r638 = vector.transfer_read %sv634[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %p639 = arith.addf %r638, %fa635 : vector<8x8xf64>
      vector.transfer_write %p639, %sv634[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    }
    gpu.barrier
    %sv640 = memref.subview %v3592[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv641 = memref.subview %v3593[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc642 = arith.constant dense<0.0> : vector<8x8xf64>
      %a643 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b644 = vector.transfer_read %sv640[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r645 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a643, %b644, %acc642 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a646 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b647 = vector.transfer_read %sv640[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r648 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a646, %b647, %r645 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g649 = vector.transfer_read %el7[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt650 = arith.mulf %r648, %g649 : vector<8x8xf64>
      %acc651 = arith.constant dense<0.0> : vector<8x8xf64>
      %a652 = vector.transfer_read %sv640[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b653 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r654 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a652, %b653, %acc651 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a655 = vector.transfer_read %sv640[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b656 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r657 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a655, %b656, %r654 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g658 = vector.transfer_read %el6[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt659 = arith.mulf %r657, %g658 : vector<8x8xf64>
      %dv660 = vector.transfer_read %sv641[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g661 = vector.transfer_read %el8[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t662 = arith.mulf %dv660, %g661 : vector<8x8xf64>
      %t2663 = arith.addf %t662, %dt659 : vector<8x8xf64>
      %fl664 = arith.addf %t2663, %dt650 : vector<8x8xf64>
      vector.transfer_write %fl664, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc665 = arith.constant dense<0.0> : vector<8x8xf64>
      %a666 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b667 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r668 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a666, %b667, %acc665 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a669 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b670 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r671 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a669, %b670, %r668 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r671, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc672 = arith.constant dense<0.0> : vector<8x8xf64>
      %a673 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b674 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r675 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a673, %b674, %acc672 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a676 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b677 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r678 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a676, %b677, %r675 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r678, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv679 = memref.subview %el2[0, 1, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv680 = memref.subview %el2[0, 2, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa681 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r682 = vector.transfer_read %sv679[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %m683 = arith.subf %r682, %fa681 : vector<8x8xf64>
      vector.transfer_write %m683, %sv679[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
      %r684 = vector.transfer_read %sv680[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %p685 = arith.addf %r684, %fa681 : vector<8x8xf64>
      vector.transfer_write %p685, %sv680[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    }
    gpu.barrier
    %sv686 = memref.subview %v3592[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv687 = memref.subview %v3593[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc688 = arith.constant dense<0.0> : vector<8x8xf64>
      %a689 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b690 = vector.transfer_read %sv686[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r691 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a689, %b690, %acc688 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a692 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b693 = vector.transfer_read %sv686[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r694 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a692, %b693, %r691 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g695 = vector.transfer_read %el7[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt696 = arith.mulf %r694, %g695 : vector<8x8xf64>
      %acc697 = arith.constant dense<0.0> : vector<8x8xf64>
      %a698 = vector.transfer_read %sv686[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b699 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r700 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a698, %b699, %acc697 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a701 = vector.transfer_read %sv686[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b702 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r703 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a701, %b702, %r700 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g704 = vector.transfer_read %el6[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt705 = arith.mulf %r703, %g704 : vector<8x8xf64>
      %dv706 = vector.transfer_read %sv687[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g707 = vector.transfer_read %el8[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t708 = arith.mulf %dv706, %g707 : vector<8x8xf64>
      %t2709 = arith.addf %t708, %dt705 : vector<8x8xf64>
      %fl710 = arith.addf %t2709, %dt696 : vector<8x8xf64>
      vector.transfer_write %fl710, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc711 = arith.constant dense<0.0> : vector<8x8xf64>
      %a712 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b713 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r714 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a712, %b713, %acc711 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a715 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b716 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r717 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a715, %b716, %r714 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r717, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc718 = arith.constant dense<0.0> : vector<8x8xf64>
      %a719 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b720 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r721 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a719, %b720, %acc718 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a722 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b723 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r724 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a722, %b723, %r721 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r724, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv725 = memref.subview %el2[0, 2, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv726 = memref.subview %el2[0, 3, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa727 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r728 = vector.transfer_read %sv725[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %m729 = arith.subf %r728, %fa727 : vector<8x8xf64>
      vector.transfer_write %m729, %sv725[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
      %r730 = vector.transfer_read %sv726[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %p731 = arith.addf %r730, %fa727 : vector<8x8xf64>
      vector.transfer_write %p731, %sv726[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    }
    gpu.barrier
    %sv732 = memref.subview %v3592[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv733 = memref.subview %v3593[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc734 = arith.constant dense<0.0> : vector<8x8xf64>
      %a735 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b736 = vector.transfer_read %sv732[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r737 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a735, %b736, %acc734 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a738 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b739 = vector.transfer_read %sv732[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r740 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a738, %b739, %r737 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g741 = vector.transfer_read %el7[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt742 = arith.mulf %r740, %g741 : vector<8x8xf64>
      %acc743 = arith.constant dense<0.0> : vector<8x8xf64>
      %a744 = vector.transfer_read %sv732[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b745 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r746 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a744, %b745, %acc743 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a747 = vector.transfer_read %sv732[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b748 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r749 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a747, %b748, %r746 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g750 = vector.transfer_read %el6[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt751 = arith.mulf %r749, %g750 : vector<8x8xf64>
      %dv752 = vector.transfer_read %sv733[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g753 = vector.transfer_read %el8[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t754 = arith.mulf %dv752, %g753 : vector<8x8xf64>
      %t2755 = arith.addf %t754, %dt751 : vector<8x8xf64>
      %fl756 = arith.addf %t2755, %dt742 : vector<8x8xf64>
      vector.transfer_write %fl756, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc757 = arith.constant dense<0.0> : vector<8x8xf64>
      %a758 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b759 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r760 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a758, %b759, %acc757 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a761 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b762 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r763 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a761, %b762, %r760 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r763, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc764 = arith.constant dense<0.0> : vector<8x8xf64>
      %a765 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b766 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r767 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a765, %b766, %acc764 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a768 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b769 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r770 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a768, %b769, %r767 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r770, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv771 = memref.subview %el2[0, 3, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv772 = memref.subview %el2[0, 4, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa773 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r774 = vector.transfer_read %sv771[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %m775 = arith.subf %r774, %fa773 : vector<8x8xf64>
      vector.transfer_write %m775, %sv771[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
      %r776 = vector.transfer_read %sv772[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %p777 = arith.addf %r776, %fa773 : vector<8x8xf64>
      vector.transfer_write %p777, %sv772[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    }
    gpu.barrier
    %sv778 = memref.subview %v3592[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv779 = memref.subview %v3593[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc780 = arith.constant dense<0.0> : vector<8x8xf64>
      %a781 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b782 = vector.transfer_read %sv778[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r783 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a781, %b782, %acc780 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a784 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b785 = vector.transfer_read %sv778[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r786 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a784, %b785, %r783 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g787 = vector.transfer_read %el7[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt788 = arith.mulf %r786, %g787 : vector<8x8xf64>
      %acc789 = arith.constant dense<0.0> : vector<8x8xf64>
      %a790 = vector.transfer_read %sv778[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b791 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r792 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a790, %b791, %acc789 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a793 = vector.transfer_read %sv778[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b794 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r795 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a793, %b794, %r792 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g796 = vector.transfer_read %el6[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt797 = arith.mulf %r795, %g796 : vector<8x8xf64>
      %dv798 = vector.transfer_read %sv779[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g799 = vector.transfer_read %el8[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t800 = arith.mulf %dv798, %g799 : vector<8x8xf64>
      %t2801 = arith.addf %t800, %dt797 : vector<8x8xf64>
      %fl802 = arith.addf %t2801, %dt788 : vector<8x8xf64>
      vector.transfer_write %fl802, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc803 = arith.constant dense<0.0> : vector<8x8xf64>
      %a804 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b805 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r806 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a804, %b805, %acc803 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a807 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b808 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r809 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a807, %b808, %r806 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r809, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc810 = arith.constant dense<0.0> : vector<8x8xf64>
      %a811 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b812 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r813 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a811, %b812, %acc810 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a814 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b815 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r816 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a814, %b815, %r813 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r816, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv817 = memref.subview %el2[0, 4, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv818 = memref.subview %el2[0, 5, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa819 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r820 = vector.transfer_read %sv817[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %m821 = arith.subf %r820, %fa819 : vector<8x8xf64>
      vector.transfer_write %m821, %sv817[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
      %r822 = vector.transfer_read %sv818[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %p823 = arith.addf %r822, %fa819 : vector<8x8xf64>
      vector.transfer_write %p823, %sv818[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    }
    gpu.barrier
    %sv824 = memref.subview %v3592[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv825 = memref.subview %v3593[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc826 = arith.constant dense<0.0> : vector<8x8xf64>
      %a827 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b828 = vector.transfer_read %sv824[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r829 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a827, %b828, %acc826 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a830 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b831 = vector.transfer_read %sv824[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r832 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a830, %b831, %r829 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g833 = vector.transfer_read %el7[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt834 = arith.mulf %r832, %g833 : vector<8x8xf64>
      %acc835 = arith.constant dense<0.0> : vector<8x8xf64>
      %a836 = vector.transfer_read %sv824[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b837 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r838 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a836, %b837, %acc835 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a839 = vector.transfer_read %sv824[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b840 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r841 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a839, %b840, %r838 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g842 = vector.transfer_read %el6[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt843 = arith.mulf %r841, %g842 : vector<8x8xf64>
      %dv844 = vector.transfer_read %sv825[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g845 = vector.transfer_read %el8[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t846 = arith.mulf %dv844, %g845 : vector<8x8xf64>
      %t2847 = arith.addf %t846, %dt843 : vector<8x8xf64>
      %fl848 = arith.addf %t2847, %dt834 : vector<8x8xf64>
      vector.transfer_write %fl848, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc849 = arith.constant dense<0.0> : vector<8x8xf64>
      %a850 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b851 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r852 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a850, %b851, %acc849 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a853 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b854 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r855 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a853, %b854, %r852 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r855, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc856 = arith.constant dense<0.0> : vector<8x8xf64>
      %a857 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b858 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r859 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a857, %b858, %acc856 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a860 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b861 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r862 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a860, %b861, %r859 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r862, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv863 = memref.subview %el2[0, 5, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv864 = memref.subview %el2[0, 6, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa865 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r866 = vector.transfer_read %sv863[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %m867 = arith.subf %r866, %fa865 : vector<8x8xf64>
      vector.transfer_write %m867, %sv863[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
      %r868 = vector.transfer_read %sv864[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %p869 = arith.addf %r868, %fa865 : vector<8x8xf64>
      vector.transfer_write %p869, %sv864[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    }
    gpu.barrier
    %sv870 = memref.subview %v3592[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv871 = memref.subview %v3593[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc872 = arith.constant dense<0.0> : vector<8x8xf64>
      %a873 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b874 = vector.transfer_read %sv870[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r875 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a873, %b874, %acc872 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a876 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b877 = vector.transfer_read %sv870[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r878 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a876, %b877, %r875 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g879 = vector.transfer_read %el7[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt880 = arith.mulf %r878, %g879 : vector<8x8xf64>
      %acc881 = arith.constant dense<0.0> : vector<8x8xf64>
      %a882 = vector.transfer_read %sv870[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b883 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r884 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a882, %b883, %acc881 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a885 = vector.transfer_read %sv870[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b886 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r887 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a885, %b886, %r884 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g888 = vector.transfer_read %el6[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt889 = arith.mulf %r887, %g888 : vector<8x8xf64>
      %dv890 = vector.transfer_read %sv871[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g891 = vector.transfer_read %el8[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t892 = arith.mulf %dv890, %g891 : vector<8x8xf64>
      %t2893 = arith.addf %t892, %dt889 : vector<8x8xf64>
      %fl894 = arith.addf %t2893, %dt880 : vector<8x8xf64>
      vector.transfer_write %fl894, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc895 = arith.constant dense<0.0> : vector<8x8xf64>
      %a896 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b897 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r898 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a896, %b897, %acc895 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a899 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b900 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r901 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a899, %b900, %r898 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r901, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc902 = arith.constant dense<0.0> : vector<8x8xf64>
      %a903 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b904 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r905 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a903, %b904, %acc902 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a906 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b907 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r908 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a906, %b907, %r905 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r908, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv909 = memref.subview %el2[0, 6, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv910 = memref.subview %el2[0, 7, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa911 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r912 = vector.transfer_read %sv909[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %m913 = arith.subf %r912, %fa911 : vector<8x8xf64>
      vector.transfer_write %m913, %sv909[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
      %r914 = vector.transfer_read %sv910[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %p915 = arith.addf %r914, %fa911 : vector<8x8xf64>
      vector.transfer_write %p915, %sv910[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    }
    gpu.barrier
    %sv916 = memref.subview %el1[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc917 = arith.constant dense<0.0> : vector<8x8xf64>
      %a918 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b919 = vector.transfer_read %sv916[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r920 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a918, %b919, %acc917 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a921 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b922 = vector.transfer_read %sv916[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r923 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a921, %b922, %r920 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r923, %interp_all[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv924 = memref.subview %el1[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc925 = arith.constant dense<0.0> : vector<8x8xf64>
      %a926 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b927 = vector.transfer_read %sv924[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r928 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a926, %b927, %acc925 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a929 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b930 = vector.transfer_read %sv924[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r931 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a929, %b930, %r928 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r931, %interp_all[%c0, %c8] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv932 = memref.subview %el1[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc933 = arith.constant dense<0.0> : vector<8x8xf64>
      %a934 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b935 = vector.transfer_read %sv932[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r936 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a934, %b935, %acc933 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a937 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b938 = vector.transfer_read %sv932[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r939 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a937, %b938, %r936 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r939, %interp_all[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv940 = memref.subview %el1[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc941 = arith.constant dense<0.0> : vector<8x8xf64>
      %a942 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b943 = vector.transfer_read %sv940[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r944 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a942, %b943, %acc941 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a945 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b946 = vector.transfer_read %sv940[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r947 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a945, %b946, %r944 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r947, %interp_all[%c0, %c24] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv948 = memref.subview %el1[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc949 = arith.constant dense<0.0> : vector<8x8xf64>
      %a950 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b951 = vector.transfer_read %sv948[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r952 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a950, %b951, %acc949 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a953 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b954 = vector.transfer_read %sv948[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r955 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a953, %b954, %r952 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r955, %interp_all[%c0, %c32] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv956 = memref.subview %el1[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc957 = arith.constant dense<0.0> : vector<8x8xf64>
      %a958 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b959 = vector.transfer_read %sv956[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r960 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a958, %b959, %acc957 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a961 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b962 = vector.transfer_read %sv956[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r963 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a961, %b962, %r960 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r963, %interp_all[%c0, %c40] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv964 = memref.subview %el1[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc965 = arith.constant dense<0.0> : vector<8x8xf64>
      %a966 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b967 = vector.transfer_read %sv964[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r968 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a966, %b967, %acc965 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a969 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b970 = vector.transfer_read %sv964[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r971 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a969, %b970, %r968 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r971, %interp_all[%c0, %c48] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv972 = memref.subview %el1[7, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc973 = arith.constant dense<0.0> : vector<8x8xf64>
      %a974 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b975 = vector.transfer_read %sv972[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r976 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a974, %b975, %acc973 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a977 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b978 = vector.transfer_read %sv972[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r979 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a977, %b978, %r976 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r979, %interp_all[%c0, %c56] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv980 = memref.subview %el1[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc981 = arith.constant dense<0.0> : vector<8x8xf64>
      %a982 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b983 = vector.transfer_read %sv980[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r984 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a982, %b983, %acc981 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a985 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b986 = vector.transfer_read %sv980[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r987 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a985, %b986, %r984 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r987, %deriv_all[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv988 = memref.subview %el1[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc989 = arith.constant dense<0.0> : vector<8x8xf64>
      %a990 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b991 = vector.transfer_read %sv988[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r992 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a990, %b991, %acc989 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a993 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b994 = vector.transfer_read %sv988[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r995 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a993, %b994, %r992 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r995, %deriv_all[%c0, %c8] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv996 = memref.subview %el1[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc997 = arith.constant dense<0.0> : vector<8x8xf64>
      %a998 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b999 = vector.transfer_read %sv996[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1000 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a998, %b999, %acc997 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1001 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1002 = vector.transfer_read %sv996[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1003 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1001, %b1002, %r1000 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1003, %deriv_all[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1004 = memref.subview %el1[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1005 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1006 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1007 = vector.transfer_read %sv1004[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1008 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1006, %b1007, %acc1005 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1009 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1010 = vector.transfer_read %sv1004[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1011 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1009, %b1010, %r1008 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1011, %deriv_all[%c0, %c24] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1012 = memref.subview %el1[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1013 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1014 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1015 = vector.transfer_read %sv1012[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1016 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1014, %b1015, %acc1013 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1017 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1018 = vector.transfer_read %sv1012[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1019 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1017, %b1018, %r1016 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1019, %deriv_all[%c0, %c32] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1020 = memref.subview %el1[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1021 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1022 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1023 = vector.transfer_read %sv1020[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1024 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1022, %b1023, %acc1021 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1025 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1026 = vector.transfer_read %sv1020[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1027 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1025, %b1026, %r1024 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1027, %deriv_all[%c0, %c40] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1028 = memref.subview %el1[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1029 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1030 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1031 = vector.transfer_read %sv1028[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1032 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1030, %b1031, %acc1029 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1033 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1034 = vector.transfer_read %sv1028[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1035 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1033, %b1034, %r1032 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1035, %deriv_all[%c0, %c48] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1036 = memref.subview %el1[7, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1037 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1038 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1039 = vector.transfer_read %sv1036[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1040 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1038, %b1039, %acc1037 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1041 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1042 = vector.transfer_read %sv1036[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1043 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1041, %b1042, %r1040 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1043, %deriv_all[%c0, %c56] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %v31044 = memref.reinterpret_cast %interp_all to offset: [0], sizes: [8, 8, 8], strides: [64, 8, 1] : memref<8x64xf64, #gpu.address_space<workgroup>> to memref<8x8x8xf64, #gpu.address_space<workgroup>>
    %v31045 = memref.reinterpret_cast %deriv_all to offset: [0], sizes: [8, 8, 8], strides: [64, 8, 1] : memref<8x64xf64, #gpu.address_space<workgroup>> to memref<8x8x8xf64, #gpu.address_space<workgroup>>
    %sv1046 = memref.subview %v31044[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>
    %sv1047 = memref.subview %v31045[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1048 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1049 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1050 = vector.transfer_read %sv1046[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1051 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1049, %b1050, %acc1048 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1052 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1053 = vector.transfer_read %sv1046[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1054 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1052, %b1053, %r1051 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1055 = vector.transfer_read %el10[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt1056 = arith.mulf %r1054, %g1055 : vector<8x8xf64>
      %acc1057 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1058 = vector.transfer_read %sv1046[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1059 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1060 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1058, %b1059, %acc1057 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1061 = vector.transfer_read %sv1046[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1062 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1063 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1061, %b1062, %r1060 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1064 = vector.transfer_read %el9[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt1065 = arith.mulf %r1063, %g1064 : vector<8x8xf64>
      %dv1066 = vector.transfer_read %sv1047[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g1067 = vector.transfer_read %el11[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t1068 = arith.mulf %dv1066, %g1067 : vector<8x8xf64>
      %t21069 = arith.addf %t1068, %dt1065 : vector<8x8xf64>
      %fl1070 = arith.addf %t21069, %dt1056 : vector<8x8xf64>
      vector.transfer_write %fl1070, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1071 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1072 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1073 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1074 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1072, %b1073, %acc1071 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1075 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1076 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1077 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1075, %b1076, %r1074 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1077, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1078 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1079 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1080 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1081 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1079, %b1080, %acc1078 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1082 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1083 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1084 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1082, %b1083, %r1081 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1084, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1085 = memref.subview %el2[0, 0, 0] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv1086 = memref.subview %el2[0, 0, 1] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1087 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1088 = vector.transfer_read %sv1085[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %m1089 = arith.subf %r1088, %fa1087 : vector<8x8xf64>
      vector.transfer_write %m1089, %sv1085[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
      %r1090 = vector.transfer_read %sv1086[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %p1091 = arith.addf %r1090, %fa1087 : vector<8x8xf64>
      vector.transfer_write %p1091, %sv1086[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    }
    gpu.barrier
    %sv1092 = memref.subview %v31044[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv1093 = memref.subview %v31045[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1094 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1095 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1096 = vector.transfer_read %sv1092[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1097 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1095, %b1096, %acc1094 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1098 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1099 = vector.transfer_read %sv1092[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1100 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1098, %b1099, %r1097 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1101 = vector.transfer_read %el10[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt1102 = arith.mulf %r1100, %g1101 : vector<8x8xf64>
      %acc1103 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1104 = vector.transfer_read %sv1092[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1105 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1106 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1104, %b1105, %acc1103 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1107 = vector.transfer_read %sv1092[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1108 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1109 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1107, %b1108, %r1106 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1110 = vector.transfer_read %el9[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt1111 = arith.mulf %r1109, %g1110 : vector<8x8xf64>
      %dv1112 = vector.transfer_read %sv1093[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g1113 = vector.transfer_read %el11[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t1114 = arith.mulf %dv1112, %g1113 : vector<8x8xf64>
      %t21115 = arith.addf %t1114, %dt1111 : vector<8x8xf64>
      %fl1116 = arith.addf %t21115, %dt1102 : vector<8x8xf64>
      vector.transfer_write %fl1116, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1117 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1118 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1119 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1120 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1118, %b1119, %acc1117 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1121 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1122 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1123 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1121, %b1122, %r1120 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1123, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1124 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1125 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1126 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1127 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1125, %b1126, %acc1124 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1128 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1129 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1130 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1128, %b1129, %r1127 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1130, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1131 = memref.subview %el2[0, 0, 1] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv1132 = memref.subview %el2[0, 0, 2] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1133 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1134 = vector.transfer_read %sv1131[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %m1135 = arith.subf %r1134, %fa1133 : vector<8x8xf64>
      vector.transfer_write %m1135, %sv1131[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
      %r1136 = vector.transfer_read %sv1132[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %p1137 = arith.addf %r1136, %fa1133 : vector<8x8xf64>
      vector.transfer_write %p1137, %sv1132[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    }
    gpu.barrier
    %sv1138 = memref.subview %v31044[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv1139 = memref.subview %v31045[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1140 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1141 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1142 = vector.transfer_read %sv1138[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1143 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1141, %b1142, %acc1140 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1144 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1145 = vector.transfer_read %sv1138[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1146 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1144, %b1145, %r1143 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1147 = vector.transfer_read %el10[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt1148 = arith.mulf %r1146, %g1147 : vector<8x8xf64>
      %acc1149 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1150 = vector.transfer_read %sv1138[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1151 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1152 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1150, %b1151, %acc1149 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1153 = vector.transfer_read %sv1138[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1154 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1155 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1153, %b1154, %r1152 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1156 = vector.transfer_read %el9[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt1157 = arith.mulf %r1155, %g1156 : vector<8x8xf64>
      %dv1158 = vector.transfer_read %sv1139[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g1159 = vector.transfer_read %el11[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t1160 = arith.mulf %dv1158, %g1159 : vector<8x8xf64>
      %t21161 = arith.addf %t1160, %dt1157 : vector<8x8xf64>
      %fl1162 = arith.addf %t21161, %dt1148 : vector<8x8xf64>
      vector.transfer_write %fl1162, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1163 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1164 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1165 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1166 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1164, %b1165, %acc1163 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1167 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1168 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1169 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1167, %b1168, %r1166 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1169, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1170 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1171 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1172 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1173 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1171, %b1172, %acc1170 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1174 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1175 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1176 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1174, %b1175, %r1173 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1176, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1177 = memref.subview %el2[0, 0, 2] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv1178 = memref.subview %el2[0, 0, 3] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1179 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1180 = vector.transfer_read %sv1177[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %m1181 = arith.subf %r1180, %fa1179 : vector<8x8xf64>
      vector.transfer_write %m1181, %sv1177[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
      %r1182 = vector.transfer_read %sv1178[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %p1183 = arith.addf %r1182, %fa1179 : vector<8x8xf64>
      vector.transfer_write %p1183, %sv1178[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    }
    gpu.barrier
    %sv1184 = memref.subview %v31044[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv1185 = memref.subview %v31045[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1186 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1187 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1188 = vector.transfer_read %sv1184[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1189 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1187, %b1188, %acc1186 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1190 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1191 = vector.transfer_read %sv1184[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1192 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1190, %b1191, %r1189 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1193 = vector.transfer_read %el10[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt1194 = arith.mulf %r1192, %g1193 : vector<8x8xf64>
      %acc1195 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1196 = vector.transfer_read %sv1184[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1197 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1198 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1196, %b1197, %acc1195 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1199 = vector.transfer_read %sv1184[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1200 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1201 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1199, %b1200, %r1198 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1202 = vector.transfer_read %el9[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt1203 = arith.mulf %r1201, %g1202 : vector<8x8xf64>
      %dv1204 = vector.transfer_read %sv1185[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g1205 = vector.transfer_read %el11[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t1206 = arith.mulf %dv1204, %g1205 : vector<8x8xf64>
      %t21207 = arith.addf %t1206, %dt1203 : vector<8x8xf64>
      %fl1208 = arith.addf %t21207, %dt1194 : vector<8x8xf64>
      vector.transfer_write %fl1208, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1209 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1210 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1211 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1212 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1210, %b1211, %acc1209 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1213 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1214 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1215 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1213, %b1214, %r1212 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1215, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1216 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1217 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1218 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1219 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1217, %b1218, %acc1216 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1220 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1221 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1222 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1220, %b1221, %r1219 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1222, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1223 = memref.subview %el2[0, 0, 3] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv1224 = memref.subview %el2[0, 0, 4] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1225 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1226 = vector.transfer_read %sv1223[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %m1227 = arith.subf %r1226, %fa1225 : vector<8x8xf64>
      vector.transfer_write %m1227, %sv1223[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
      %r1228 = vector.transfer_read %sv1224[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %p1229 = arith.addf %r1228, %fa1225 : vector<8x8xf64>
      vector.transfer_write %p1229, %sv1224[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    }
    gpu.barrier
    %sv1230 = memref.subview %v31044[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv1231 = memref.subview %v31045[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1232 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1233 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1234 = vector.transfer_read %sv1230[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1235 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1233, %b1234, %acc1232 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1236 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1237 = vector.transfer_read %sv1230[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1238 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1236, %b1237, %r1235 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1239 = vector.transfer_read %el10[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt1240 = arith.mulf %r1238, %g1239 : vector<8x8xf64>
      %acc1241 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1242 = vector.transfer_read %sv1230[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1243 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1244 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1242, %b1243, %acc1241 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1245 = vector.transfer_read %sv1230[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1246 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1247 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1245, %b1246, %r1244 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1248 = vector.transfer_read %el9[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt1249 = arith.mulf %r1247, %g1248 : vector<8x8xf64>
      %dv1250 = vector.transfer_read %sv1231[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g1251 = vector.transfer_read %el11[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t1252 = arith.mulf %dv1250, %g1251 : vector<8x8xf64>
      %t21253 = arith.addf %t1252, %dt1249 : vector<8x8xf64>
      %fl1254 = arith.addf %t21253, %dt1240 : vector<8x8xf64>
      vector.transfer_write %fl1254, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1255 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1256 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1257 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1258 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1256, %b1257, %acc1255 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1259 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1260 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1261 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1259, %b1260, %r1258 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1261, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1262 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1263 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1264 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1265 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1263, %b1264, %acc1262 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1266 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1267 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1268 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1266, %b1267, %r1265 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1268, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1269 = memref.subview %el2[0, 0, 4] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv1270 = memref.subview %el2[0, 0, 5] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1271 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1272 = vector.transfer_read %sv1269[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %m1273 = arith.subf %r1272, %fa1271 : vector<8x8xf64>
      vector.transfer_write %m1273, %sv1269[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
      %r1274 = vector.transfer_read %sv1270[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %p1275 = arith.addf %r1274, %fa1271 : vector<8x8xf64>
      vector.transfer_write %p1275, %sv1270[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    }
    gpu.barrier
    %sv1276 = memref.subview %v31044[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv1277 = memref.subview %v31045[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1278 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1279 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1280 = vector.transfer_read %sv1276[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1281 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1279, %b1280, %acc1278 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1282 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1283 = vector.transfer_read %sv1276[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1284 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1282, %b1283, %r1281 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1285 = vector.transfer_read %el10[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt1286 = arith.mulf %r1284, %g1285 : vector<8x8xf64>
      %acc1287 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1288 = vector.transfer_read %sv1276[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1289 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1290 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1288, %b1289, %acc1287 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1291 = vector.transfer_read %sv1276[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1292 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1293 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1291, %b1292, %r1290 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1294 = vector.transfer_read %el9[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt1295 = arith.mulf %r1293, %g1294 : vector<8x8xf64>
      %dv1296 = vector.transfer_read %sv1277[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g1297 = vector.transfer_read %el11[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t1298 = arith.mulf %dv1296, %g1297 : vector<8x8xf64>
      %t21299 = arith.addf %t1298, %dt1295 : vector<8x8xf64>
      %fl1300 = arith.addf %t21299, %dt1286 : vector<8x8xf64>
      vector.transfer_write %fl1300, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1301 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1302 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1303 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1304 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1302, %b1303, %acc1301 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1305 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1306 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1307 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1305, %b1306, %r1304 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1307, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1308 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1309 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1310 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1311 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1309, %b1310, %acc1308 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1312 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1313 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1314 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1312, %b1313, %r1311 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1314, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1315 = memref.subview %el2[0, 0, 5] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv1316 = memref.subview %el2[0, 0, 6] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1317 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1318 = vector.transfer_read %sv1315[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %m1319 = arith.subf %r1318, %fa1317 : vector<8x8xf64>
      vector.transfer_write %m1319, %sv1315[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
      %r1320 = vector.transfer_read %sv1316[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %p1321 = arith.addf %r1320, %fa1317 : vector<8x8xf64>
      vector.transfer_write %p1321, %sv1316[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    }
    gpu.barrier
    %sv1322 = memref.subview %v31044[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv1323 = memref.subview %v31045[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1324 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1325 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1326 = vector.transfer_read %sv1322[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1327 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1325, %b1326, %acc1324 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1328 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1329 = vector.transfer_read %sv1322[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1330 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1328, %b1329, %r1327 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1331 = vector.transfer_read %el10[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt1332 = arith.mulf %r1330, %g1331 : vector<8x8xf64>
      %acc1333 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1334 = vector.transfer_read %sv1322[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1335 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1336 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1334, %b1335, %acc1333 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1337 = vector.transfer_read %sv1322[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1338 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1339 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1337, %b1338, %r1336 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1340 = vector.transfer_read %el9[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %dt1341 = arith.mulf %r1339, %g1340 : vector<8x8xf64>
      %dv1342 = vector.transfer_read %sv1323[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g1343 = vector.transfer_read %el11[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t1344 = arith.mulf %dv1342, %g1343 : vector<8x8xf64>
      %t21345 = arith.addf %t1344, %dt1341 : vector<8x8xf64>
      %fl1346 = arith.addf %t21345, %dt1332 : vector<8x8xf64>
      vector.transfer_write %fl1346, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1347 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1348 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1349 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1350 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1348, %b1349, %acc1347 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1351 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1352 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1353 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1351, %b1352, %r1350 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1353, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1354 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1355 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1356 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1357 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1355, %b1356, %acc1354 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1358 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1359 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1360 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1358, %b1359, %r1357 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1360, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1361 = memref.subview %el2[0, 0, 6] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv1362 = memref.subview %el2[0, 0, 7] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1363 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1364 = vector.transfer_read %sv1361[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %m1365 = arith.subf %r1364, %fa1363 : vector<8x8xf64>
      vector.transfer_write %m1365, %sv1361[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
      %r1366 = vector.transfer_read %sv1362[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %p1367 = arith.addf %r1366, %fa1363 : vector<8x8xf64>
      vector.transfer_write %p1367, %sv1362[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    }
    gpu.barrier
    gpu.return
  }
}
