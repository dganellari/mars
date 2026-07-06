#a = affine_map<(m, n, k) -> (m, k)>
#b = affine_map<(m, n, k) -> (k, n)>
#bt = affine_map<(m, n, k) -> (n, k)>
#c = affine_map<(m, n, k) -> (m, n)>
gpu.module @mir_kernels [#nvvm.target<chip = "sm_90", O = 3>] {
  gpu.func @full_batched(%U: memref<?x8x8x8xf64>, %Btil: memref<8x8xf64>, %Dtil: memref<8x8xf64>, %Dm: memref<8x8xf64>, %W: memref<8x8xf64>, %g00all: memref<?x8x8x8xf64>, %g01all: memref<?x8x8x8xf64>, %g02all: memref<?x8x8x8xf64>, %g10all: memref<?x8x8x8xf64>, %g11all: memref<?x8x8x8xf64>, %g12all: memref<?x8x8x8xf64>, %g20all: memref<?x8x8x8xf64>, %g21all: memref<?x8x8x8xf64>, %g22all: memref<?x8x8x8xf64>, %Y: memref<?x8x8x8xf64>) workgroup(%interp_all: memref<8x64xf64, #gpu.address_space<workgroup>>, %deriv_all: memref<8x64xf64, #gpu.address_space<workgroup>>, %dt1g: memref<8x8xf64, #gpu.address_space<workgroup>>, %dt2g: memref<8x8xf64, #gpu.address_space<workgroup>>, %flux: memref<8x8xf64, #gpu.address_space<workgroup>>, %tmp: memref<8x8xf64, #gpu.address_space<workgroup>>, %intf: memref<8x8xf64, #gpu.address_space<workgroup>>) kernel {
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
      %sc155 = arith.mulf %r153, %g154 : vector<8x8xf64>
      vector.transfer_write %sc155, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc156 = arith.constant dense<0.0> : vector<8x8xf64>
      %a157 = vector.transfer_read %sv142[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b158 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r159 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a157, %b158, %acc156 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a160 = vector.transfer_read %sv142[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b161 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r162 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a160, %b161, %r159 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g163 = vector.transfer_read %sv144[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc164 = arith.mulf %r162, %g163 : vector<8x8xf64>
      vector.transfer_write %sc164, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d165 = vector.transfer_read %sv143[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2166 = vector.transfer_read %sv146[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t167 = arith.mulf %d165, %g2166 : vector<8x8xf64>
      %a168 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2169 = arith.addf %t167, %a168 : vector<8x8xf64>
      %b170 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl171 = arith.addf %t2169, %b170 : vector<8x8xf64>
      vector.transfer_write %fl171, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc172 = arith.constant dense<0.0> : vector<8x8xf64>
      %a173 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b174 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r175 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a173, %b174, %acc172 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a176 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b177 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r178 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a176, %b177, %r175 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r178, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc179 = arith.constant dense<0.0> : vector<8x8xf64>
      %a180 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b181 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r182 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a180, %b181, %acc179 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a183 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b184 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r185 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a183, %b184, %r182 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r185, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv186 = memref.subview %el2[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv187 = memref.subview %el2[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa188 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r189 = vector.transfer_read %sv186[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %m190 = arith.subf %r189, %fa188 : vector<8x8xf64>
      vector.transfer_write %m190, %sv186[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      %fa191 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r192 = vector.transfer_read %sv187[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %p193 = arith.addf %r192, %fa191 : vector<8x8xf64>
      vector.transfer_write %p193, %sv187[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    }
    gpu.barrier
    %sv194 = memref.subview %v3140[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv195 = memref.subview %v3141[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv196 = memref.subview %el3[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv197 = memref.subview %el4[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv198 = memref.subview %el5[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc199 = arith.constant dense<0.0> : vector<8x8xf64>
      %a200 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b201 = vector.transfer_read %sv194[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r202 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a200, %b201, %acc199 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a203 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b204 = vector.transfer_read %sv194[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r205 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a203, %b204, %r202 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g206 = vector.transfer_read %sv197[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc207 = arith.mulf %r205, %g206 : vector<8x8xf64>
      vector.transfer_write %sc207, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc208 = arith.constant dense<0.0> : vector<8x8xf64>
      %a209 = vector.transfer_read %sv194[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b210 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r211 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a209, %b210, %acc208 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a212 = vector.transfer_read %sv194[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b213 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r214 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a212, %b213, %r211 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g215 = vector.transfer_read %sv196[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc216 = arith.mulf %r214, %g215 : vector<8x8xf64>
      vector.transfer_write %sc216, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d217 = vector.transfer_read %sv195[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2218 = vector.transfer_read %sv198[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t219 = arith.mulf %d217, %g2218 : vector<8x8xf64>
      %a220 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2221 = arith.addf %t219, %a220 : vector<8x8xf64>
      %b222 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl223 = arith.addf %t2221, %b222 : vector<8x8xf64>
      vector.transfer_write %fl223, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc224 = arith.constant dense<0.0> : vector<8x8xf64>
      %a225 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b226 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r227 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a225, %b226, %acc224 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a228 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b229 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r230 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a228, %b229, %r227 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r230, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc231 = arith.constant dense<0.0> : vector<8x8xf64>
      %a232 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b233 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r234 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a232, %b233, %acc231 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a235 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b236 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r237 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a235, %b236, %r234 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r237, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv238 = memref.subview %el2[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv239 = memref.subview %el2[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa240 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r241 = vector.transfer_read %sv238[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %m242 = arith.subf %r241, %fa240 : vector<8x8xf64>
      vector.transfer_write %m242, %sv238[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      %fa243 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r244 = vector.transfer_read %sv239[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %p245 = arith.addf %r244, %fa243 : vector<8x8xf64>
      vector.transfer_write %p245, %sv239[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    }
    gpu.barrier
    %sv246 = memref.subview %v3140[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv247 = memref.subview %v3141[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv248 = memref.subview %el3[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv249 = memref.subview %el4[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv250 = memref.subview %el5[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc251 = arith.constant dense<0.0> : vector<8x8xf64>
      %a252 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b253 = vector.transfer_read %sv246[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r254 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a252, %b253, %acc251 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a255 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b256 = vector.transfer_read %sv246[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r257 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a255, %b256, %r254 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g258 = vector.transfer_read %sv249[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc259 = arith.mulf %r257, %g258 : vector<8x8xf64>
      vector.transfer_write %sc259, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc260 = arith.constant dense<0.0> : vector<8x8xf64>
      %a261 = vector.transfer_read %sv246[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b262 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r263 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a261, %b262, %acc260 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a264 = vector.transfer_read %sv246[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b265 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r266 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a264, %b265, %r263 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g267 = vector.transfer_read %sv248[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc268 = arith.mulf %r266, %g267 : vector<8x8xf64>
      vector.transfer_write %sc268, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d269 = vector.transfer_read %sv247[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2270 = vector.transfer_read %sv250[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t271 = arith.mulf %d269, %g2270 : vector<8x8xf64>
      %a272 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2273 = arith.addf %t271, %a272 : vector<8x8xf64>
      %b274 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl275 = arith.addf %t2273, %b274 : vector<8x8xf64>
      vector.transfer_write %fl275, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc276 = arith.constant dense<0.0> : vector<8x8xf64>
      %a277 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b278 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r279 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a277, %b278, %acc276 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a280 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b281 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r282 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a280, %b281, %r279 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r282, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc283 = arith.constant dense<0.0> : vector<8x8xf64>
      %a284 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b285 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r286 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a284, %b285, %acc283 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a287 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b288 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r289 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a287, %b288, %r286 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r289, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv290 = memref.subview %el2[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv291 = memref.subview %el2[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa292 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r293 = vector.transfer_read %sv290[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %m294 = arith.subf %r293, %fa292 : vector<8x8xf64>
      vector.transfer_write %m294, %sv290[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      %fa295 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r296 = vector.transfer_read %sv291[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %p297 = arith.addf %r296, %fa295 : vector<8x8xf64>
      vector.transfer_write %p297, %sv291[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    }
    gpu.barrier
    %sv298 = memref.subview %v3140[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv299 = memref.subview %v3141[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv300 = memref.subview %el3[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv301 = memref.subview %el4[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv302 = memref.subview %el5[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc303 = arith.constant dense<0.0> : vector<8x8xf64>
      %a304 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b305 = vector.transfer_read %sv298[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r306 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a304, %b305, %acc303 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a307 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b308 = vector.transfer_read %sv298[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r309 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a307, %b308, %r306 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g310 = vector.transfer_read %sv301[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc311 = arith.mulf %r309, %g310 : vector<8x8xf64>
      vector.transfer_write %sc311, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc312 = arith.constant dense<0.0> : vector<8x8xf64>
      %a313 = vector.transfer_read %sv298[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b314 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r315 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a313, %b314, %acc312 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a316 = vector.transfer_read %sv298[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b317 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r318 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a316, %b317, %r315 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g319 = vector.transfer_read %sv300[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc320 = arith.mulf %r318, %g319 : vector<8x8xf64>
      vector.transfer_write %sc320, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d321 = vector.transfer_read %sv299[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2322 = vector.transfer_read %sv302[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t323 = arith.mulf %d321, %g2322 : vector<8x8xf64>
      %a324 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2325 = arith.addf %t323, %a324 : vector<8x8xf64>
      %b326 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl327 = arith.addf %t2325, %b326 : vector<8x8xf64>
      vector.transfer_write %fl327, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc328 = arith.constant dense<0.0> : vector<8x8xf64>
      %a329 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b330 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r331 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a329, %b330, %acc328 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a332 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b333 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r334 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a332, %b333, %r331 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r334, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc335 = arith.constant dense<0.0> : vector<8x8xf64>
      %a336 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b337 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r338 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a336, %b337, %acc335 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a339 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b340 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r341 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a339, %b340, %r338 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r341, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv342 = memref.subview %el2[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv343 = memref.subview %el2[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa344 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r345 = vector.transfer_read %sv342[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %m346 = arith.subf %r345, %fa344 : vector<8x8xf64>
      vector.transfer_write %m346, %sv342[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      %fa347 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r348 = vector.transfer_read %sv343[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %p349 = arith.addf %r348, %fa347 : vector<8x8xf64>
      vector.transfer_write %p349, %sv343[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    }
    gpu.barrier
    %sv350 = memref.subview %v3140[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv351 = memref.subview %v3141[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv352 = memref.subview %el3[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv353 = memref.subview %el4[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv354 = memref.subview %el5[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc355 = arith.constant dense<0.0> : vector<8x8xf64>
      %a356 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b357 = vector.transfer_read %sv350[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r358 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a356, %b357, %acc355 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a359 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b360 = vector.transfer_read %sv350[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r361 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a359, %b360, %r358 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g362 = vector.transfer_read %sv353[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc363 = arith.mulf %r361, %g362 : vector<8x8xf64>
      vector.transfer_write %sc363, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc364 = arith.constant dense<0.0> : vector<8x8xf64>
      %a365 = vector.transfer_read %sv350[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b366 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r367 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a365, %b366, %acc364 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a368 = vector.transfer_read %sv350[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b369 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r370 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a368, %b369, %r367 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g371 = vector.transfer_read %sv352[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc372 = arith.mulf %r370, %g371 : vector<8x8xf64>
      vector.transfer_write %sc372, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d373 = vector.transfer_read %sv351[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2374 = vector.transfer_read %sv354[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t375 = arith.mulf %d373, %g2374 : vector<8x8xf64>
      %a376 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2377 = arith.addf %t375, %a376 : vector<8x8xf64>
      %b378 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl379 = arith.addf %t2377, %b378 : vector<8x8xf64>
      vector.transfer_write %fl379, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc380 = arith.constant dense<0.0> : vector<8x8xf64>
      %a381 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b382 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r383 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a381, %b382, %acc380 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a384 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b385 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r386 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a384, %b385, %r383 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r386, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc387 = arith.constant dense<0.0> : vector<8x8xf64>
      %a388 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b389 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r390 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a388, %b389, %acc387 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a391 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b392 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r393 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a391, %b392, %r390 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r393, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv394 = memref.subview %el2[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv395 = memref.subview %el2[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa396 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r397 = vector.transfer_read %sv394[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %m398 = arith.subf %r397, %fa396 : vector<8x8xf64>
      vector.transfer_write %m398, %sv394[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      %fa399 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r400 = vector.transfer_read %sv395[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %p401 = arith.addf %r400, %fa399 : vector<8x8xf64>
      vector.transfer_write %p401, %sv395[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    }
    gpu.barrier
    %sv402 = memref.subview %v3140[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv403 = memref.subview %v3141[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv404 = memref.subview %el3[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv405 = memref.subview %el4[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv406 = memref.subview %el5[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc407 = arith.constant dense<0.0> : vector<8x8xf64>
      %a408 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b409 = vector.transfer_read %sv402[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r410 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a408, %b409, %acc407 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a411 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b412 = vector.transfer_read %sv402[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r413 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a411, %b412, %r410 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g414 = vector.transfer_read %sv405[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc415 = arith.mulf %r413, %g414 : vector<8x8xf64>
      vector.transfer_write %sc415, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc416 = arith.constant dense<0.0> : vector<8x8xf64>
      %a417 = vector.transfer_read %sv402[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b418 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r419 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a417, %b418, %acc416 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a420 = vector.transfer_read %sv402[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b421 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r422 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a420, %b421, %r419 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g423 = vector.transfer_read %sv404[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc424 = arith.mulf %r422, %g423 : vector<8x8xf64>
      vector.transfer_write %sc424, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d425 = vector.transfer_read %sv403[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2426 = vector.transfer_read %sv406[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t427 = arith.mulf %d425, %g2426 : vector<8x8xf64>
      %a428 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2429 = arith.addf %t427, %a428 : vector<8x8xf64>
      %b430 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl431 = arith.addf %t2429, %b430 : vector<8x8xf64>
      vector.transfer_write %fl431, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc432 = arith.constant dense<0.0> : vector<8x8xf64>
      %a433 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b434 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r435 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a433, %b434, %acc432 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a436 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b437 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r438 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a436, %b437, %r435 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r438, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc439 = arith.constant dense<0.0> : vector<8x8xf64>
      %a440 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b441 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r442 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a440, %b441, %acc439 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a443 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b444 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r445 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a443, %b444, %r442 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r445, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv446 = memref.subview %el2[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv447 = memref.subview %el2[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa448 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r449 = vector.transfer_read %sv446[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %m450 = arith.subf %r449, %fa448 : vector<8x8xf64>
      vector.transfer_write %m450, %sv446[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      %fa451 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r452 = vector.transfer_read %sv447[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %p453 = arith.addf %r452, %fa451 : vector<8x8xf64>
      vector.transfer_write %p453, %sv447[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    }
    gpu.barrier
    %sv454 = memref.subview %v3140[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv455 = memref.subview %v3141[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv456 = memref.subview %el3[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv457 = memref.subview %el4[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv458 = memref.subview %el5[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc459 = arith.constant dense<0.0> : vector<8x8xf64>
      %a460 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b461 = vector.transfer_read %sv454[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r462 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a460, %b461, %acc459 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a463 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b464 = vector.transfer_read %sv454[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r465 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a463, %b464, %r462 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g466 = vector.transfer_read %sv457[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc467 = arith.mulf %r465, %g466 : vector<8x8xf64>
      vector.transfer_write %sc467, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc468 = arith.constant dense<0.0> : vector<8x8xf64>
      %a469 = vector.transfer_read %sv454[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b470 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r471 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a469, %b470, %acc468 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a472 = vector.transfer_read %sv454[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b473 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r474 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a472, %b473, %r471 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g475 = vector.transfer_read %sv456[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc476 = arith.mulf %r474, %g475 : vector<8x8xf64>
      vector.transfer_write %sc476, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d477 = vector.transfer_read %sv455[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2478 = vector.transfer_read %sv458[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t479 = arith.mulf %d477, %g2478 : vector<8x8xf64>
      %a480 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2481 = arith.addf %t479, %a480 : vector<8x8xf64>
      %b482 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl483 = arith.addf %t2481, %b482 : vector<8x8xf64>
      vector.transfer_write %fl483, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc484 = arith.constant dense<0.0> : vector<8x8xf64>
      %a485 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b486 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r487 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a485, %b486, %acc484 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a488 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b489 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r490 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a488, %b489, %r487 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r490, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc491 = arith.constant dense<0.0> : vector<8x8xf64>
      %a492 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b493 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r494 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a492, %b493, %acc491 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a495 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b496 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r497 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a495, %b496, %r494 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r497, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv498 = memref.subview %el2[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv499 = memref.subview %el2[7, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa500 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r501 = vector.transfer_read %sv498[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %m502 = arith.subf %r501, %fa500 : vector<8x8xf64>
      vector.transfer_write %m502, %sv498[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      %fa503 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r504 = vector.transfer_read %sv499[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %p505 = arith.addf %r504, %fa503 : vector<8x8xf64>
      vector.transfer_write %p505, %sv499[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    }
    gpu.barrier
    %sv506 = memref.subview %el1[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc507 = arith.constant dense<0.0> : vector<8x8xf64>
      %a508 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b509 = vector.transfer_read %sv506[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r510 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a508, %b509, %acc507 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a511 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b512 = vector.transfer_read %sv506[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r513 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a511, %b512, %r510 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r513, %interp_all[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv514 = memref.subview %el1[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc515 = arith.constant dense<0.0> : vector<8x8xf64>
      %a516 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b517 = vector.transfer_read %sv514[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r518 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a516, %b517, %acc515 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a519 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b520 = vector.transfer_read %sv514[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r521 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a519, %b520, %r518 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r521, %interp_all[%c0, %c8] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv522 = memref.subview %el1[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc523 = arith.constant dense<0.0> : vector<8x8xf64>
      %a524 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b525 = vector.transfer_read %sv522[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r526 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a524, %b525, %acc523 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a527 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b528 = vector.transfer_read %sv522[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r529 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a527, %b528, %r526 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r529, %interp_all[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv530 = memref.subview %el1[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc531 = arith.constant dense<0.0> : vector<8x8xf64>
      %a532 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b533 = vector.transfer_read %sv530[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r534 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a532, %b533, %acc531 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a535 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b536 = vector.transfer_read %sv530[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r537 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a535, %b536, %r534 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r537, %interp_all[%c0, %c24] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv538 = memref.subview %el1[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc539 = arith.constant dense<0.0> : vector<8x8xf64>
      %a540 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b541 = vector.transfer_read %sv538[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r542 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a540, %b541, %acc539 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a543 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b544 = vector.transfer_read %sv538[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r545 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a543, %b544, %r542 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r545, %interp_all[%c0, %c32] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv546 = memref.subview %el1[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc547 = arith.constant dense<0.0> : vector<8x8xf64>
      %a548 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b549 = vector.transfer_read %sv546[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r550 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a548, %b549, %acc547 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a551 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b552 = vector.transfer_read %sv546[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r553 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a551, %b552, %r550 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r553, %interp_all[%c0, %c40] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv554 = memref.subview %el1[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc555 = arith.constant dense<0.0> : vector<8x8xf64>
      %a556 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b557 = vector.transfer_read %sv554[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r558 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a556, %b557, %acc555 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a559 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b560 = vector.transfer_read %sv554[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r561 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a559, %b560, %r558 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r561, %interp_all[%c0, %c48] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv562 = memref.subview %el1[7, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc563 = arith.constant dense<0.0> : vector<8x8xf64>
      %a564 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b565 = vector.transfer_read %sv562[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r566 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a564, %b565, %acc563 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a567 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b568 = vector.transfer_read %sv562[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r569 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a567, %b568, %r566 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r569, %interp_all[%c0, %c56] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv570 = memref.subview %el1[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc571 = arith.constant dense<0.0> : vector<8x8xf64>
      %a572 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b573 = vector.transfer_read %sv570[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r574 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a572, %b573, %acc571 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a575 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b576 = vector.transfer_read %sv570[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r577 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a575, %b576, %r574 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r577, %deriv_all[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv578 = memref.subview %el1[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc579 = arith.constant dense<0.0> : vector<8x8xf64>
      %a580 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b581 = vector.transfer_read %sv578[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r582 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a580, %b581, %acc579 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a583 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b584 = vector.transfer_read %sv578[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r585 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a583, %b584, %r582 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r585, %deriv_all[%c0, %c8] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv586 = memref.subview %el1[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc587 = arith.constant dense<0.0> : vector<8x8xf64>
      %a588 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b589 = vector.transfer_read %sv586[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r590 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a588, %b589, %acc587 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a591 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b592 = vector.transfer_read %sv586[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r593 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a591, %b592, %r590 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r593, %deriv_all[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv594 = memref.subview %el1[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc595 = arith.constant dense<0.0> : vector<8x8xf64>
      %a596 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b597 = vector.transfer_read %sv594[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r598 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a596, %b597, %acc595 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a599 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b600 = vector.transfer_read %sv594[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r601 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a599, %b600, %r598 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r601, %deriv_all[%c0, %c24] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv602 = memref.subview %el1[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc603 = arith.constant dense<0.0> : vector<8x8xf64>
      %a604 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b605 = vector.transfer_read %sv602[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r606 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a604, %b605, %acc603 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a607 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b608 = vector.transfer_read %sv602[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r609 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a607, %b608, %r606 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r609, %deriv_all[%c0, %c32] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv610 = memref.subview %el1[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc611 = arith.constant dense<0.0> : vector<8x8xf64>
      %a612 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b613 = vector.transfer_read %sv610[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r614 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a612, %b613, %acc611 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a615 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b616 = vector.transfer_read %sv610[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r617 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a615, %b616, %r614 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r617, %deriv_all[%c0, %c40] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv618 = memref.subview %el1[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc619 = arith.constant dense<0.0> : vector<8x8xf64>
      %a620 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b621 = vector.transfer_read %sv618[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r622 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a620, %b621, %acc619 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a623 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b624 = vector.transfer_read %sv618[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r625 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a623, %b624, %r622 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r625, %deriv_all[%c0, %c48] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv626 = memref.subview %el1[7, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc627 = arith.constant dense<0.0> : vector<8x8xf64>
      %a628 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b629 = vector.transfer_read %sv626[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r630 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a628, %b629, %acc627 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a631 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b632 = vector.transfer_read %sv626[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<4x8xf64>
      %r633 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a631, %b632, %r630 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r633, %deriv_all[%c0, %c56] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %v3634 = memref.reinterpret_cast %interp_all to offset: [0], sizes: [8, 8, 8], strides: [64, 8, 1] : memref<8x64xf64, #gpu.address_space<workgroup>> to memref<8x8x8xf64, #gpu.address_space<workgroup>>
    %v3635 = memref.reinterpret_cast %deriv_all to offset: [0], sizes: [8, 8, 8], strides: [64, 8, 1] : memref<8x64xf64, #gpu.address_space<workgroup>> to memref<8x8x8xf64, #gpu.address_space<workgroup>>
    %sv636 = memref.subview %v3634[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>
    %sv637 = memref.subview %v3635[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>
    %sv638 = memref.subview %el6[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv639 = memref.subview %el7[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv640 = memref.subview %el8[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc641 = arith.constant dense<0.0> : vector<8x8xf64>
      %a642 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b643 = vector.transfer_read %sv636[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r644 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a642, %b643, %acc641 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a645 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b646 = vector.transfer_read %sv636[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r647 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a645, %b646, %r644 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g648 = vector.transfer_read %sv639[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc649 = arith.mulf %r647, %g648 : vector<8x8xf64>
      vector.transfer_write %sc649, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc650 = arith.constant dense<0.0> : vector<8x8xf64>
      %a651 = vector.transfer_read %sv636[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b652 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r653 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a651, %b652, %acc650 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a654 = vector.transfer_read %sv636[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b655 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r656 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a654, %b655, %r653 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g657 = vector.transfer_read %sv638[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc658 = arith.mulf %r656, %g657 : vector<8x8xf64>
      vector.transfer_write %sc658, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d659 = vector.transfer_read %sv637[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2660 = vector.transfer_read %sv640[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t661 = arith.mulf %d659, %g2660 : vector<8x8xf64>
      %a662 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2663 = arith.addf %t661, %a662 : vector<8x8xf64>
      %b664 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl665 = arith.addf %t2663, %b664 : vector<8x8xf64>
      vector.transfer_write %fl665, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc666 = arith.constant dense<0.0> : vector<8x8xf64>
      %a667 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b668 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r669 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a667, %b668, %acc666 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a670 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b671 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r672 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a670, %b671, %r669 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r672, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc673 = arith.constant dense<0.0> : vector<8x8xf64>
      %a674 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b675 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r676 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a674, %b675, %acc673 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a677 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b678 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r679 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a677, %b678, %r676 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r679, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv680 = memref.subview %el2[0, 0, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv681 = memref.subview %el2[0, 1, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa682 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r683 = vector.transfer_read %sv680[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %m684 = arith.subf %r683, %fa682 : vector<8x8xf64>
      vector.transfer_write %m684, %sv680[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
      %fa685 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r686 = vector.transfer_read %sv681[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %p687 = arith.addf %r686, %fa685 : vector<8x8xf64>
      vector.transfer_write %p687, %sv681[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    }
    gpu.barrier
    %sv688 = memref.subview %v3634[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv689 = memref.subview %v3635[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv690 = memref.subview %el6[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv691 = memref.subview %el7[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv692 = memref.subview %el8[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc693 = arith.constant dense<0.0> : vector<8x8xf64>
      %a694 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b695 = vector.transfer_read %sv688[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r696 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a694, %b695, %acc693 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a697 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b698 = vector.transfer_read %sv688[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r699 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a697, %b698, %r696 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g700 = vector.transfer_read %sv691[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc701 = arith.mulf %r699, %g700 : vector<8x8xf64>
      vector.transfer_write %sc701, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc702 = arith.constant dense<0.0> : vector<8x8xf64>
      %a703 = vector.transfer_read %sv688[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b704 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r705 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a703, %b704, %acc702 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a706 = vector.transfer_read %sv688[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b707 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r708 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a706, %b707, %r705 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g709 = vector.transfer_read %sv690[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc710 = arith.mulf %r708, %g709 : vector<8x8xf64>
      vector.transfer_write %sc710, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d711 = vector.transfer_read %sv689[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2712 = vector.transfer_read %sv692[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t713 = arith.mulf %d711, %g2712 : vector<8x8xf64>
      %a714 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2715 = arith.addf %t713, %a714 : vector<8x8xf64>
      %b716 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl717 = arith.addf %t2715, %b716 : vector<8x8xf64>
      vector.transfer_write %fl717, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc718 = arith.constant dense<0.0> : vector<8x8xf64>
      %a719 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b720 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r721 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a719, %b720, %acc718 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a722 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b723 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r724 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a722, %b723, %r721 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r724, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc725 = arith.constant dense<0.0> : vector<8x8xf64>
      %a726 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b727 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r728 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a726, %b727, %acc725 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a729 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b730 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r731 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a729, %b730, %r728 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r731, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv732 = memref.subview %el2[0, 1, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv733 = memref.subview %el2[0, 2, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa734 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r735 = vector.transfer_read %sv732[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %m736 = arith.subf %r735, %fa734 : vector<8x8xf64>
      vector.transfer_write %m736, %sv732[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
      %fa737 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r738 = vector.transfer_read %sv733[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %p739 = arith.addf %r738, %fa737 : vector<8x8xf64>
      vector.transfer_write %p739, %sv733[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    }
    gpu.barrier
    %sv740 = memref.subview %v3634[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv741 = memref.subview %v3635[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv742 = memref.subview %el6[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv743 = memref.subview %el7[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv744 = memref.subview %el8[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc745 = arith.constant dense<0.0> : vector<8x8xf64>
      %a746 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b747 = vector.transfer_read %sv740[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r748 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a746, %b747, %acc745 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a749 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b750 = vector.transfer_read %sv740[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r751 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a749, %b750, %r748 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g752 = vector.transfer_read %sv743[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc753 = arith.mulf %r751, %g752 : vector<8x8xf64>
      vector.transfer_write %sc753, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc754 = arith.constant dense<0.0> : vector<8x8xf64>
      %a755 = vector.transfer_read %sv740[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b756 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r757 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a755, %b756, %acc754 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a758 = vector.transfer_read %sv740[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b759 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r760 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a758, %b759, %r757 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g761 = vector.transfer_read %sv742[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc762 = arith.mulf %r760, %g761 : vector<8x8xf64>
      vector.transfer_write %sc762, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d763 = vector.transfer_read %sv741[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2764 = vector.transfer_read %sv744[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t765 = arith.mulf %d763, %g2764 : vector<8x8xf64>
      %a766 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2767 = arith.addf %t765, %a766 : vector<8x8xf64>
      %b768 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl769 = arith.addf %t2767, %b768 : vector<8x8xf64>
      vector.transfer_write %fl769, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc770 = arith.constant dense<0.0> : vector<8x8xf64>
      %a771 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b772 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r773 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a771, %b772, %acc770 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a774 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b775 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r776 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a774, %b775, %r773 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r776, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc777 = arith.constant dense<0.0> : vector<8x8xf64>
      %a778 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b779 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r780 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a778, %b779, %acc777 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a781 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b782 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r783 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a781, %b782, %r780 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r783, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv784 = memref.subview %el2[0, 2, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv785 = memref.subview %el2[0, 3, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa786 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r787 = vector.transfer_read %sv784[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %m788 = arith.subf %r787, %fa786 : vector<8x8xf64>
      vector.transfer_write %m788, %sv784[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
      %fa789 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r790 = vector.transfer_read %sv785[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %p791 = arith.addf %r790, %fa789 : vector<8x8xf64>
      vector.transfer_write %p791, %sv785[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    }
    gpu.barrier
    %sv792 = memref.subview %v3634[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv793 = memref.subview %v3635[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv794 = memref.subview %el6[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv795 = memref.subview %el7[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv796 = memref.subview %el8[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc797 = arith.constant dense<0.0> : vector<8x8xf64>
      %a798 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b799 = vector.transfer_read %sv792[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r800 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a798, %b799, %acc797 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a801 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b802 = vector.transfer_read %sv792[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r803 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a801, %b802, %r800 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g804 = vector.transfer_read %sv795[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc805 = arith.mulf %r803, %g804 : vector<8x8xf64>
      vector.transfer_write %sc805, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc806 = arith.constant dense<0.0> : vector<8x8xf64>
      %a807 = vector.transfer_read %sv792[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b808 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r809 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a807, %b808, %acc806 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a810 = vector.transfer_read %sv792[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b811 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r812 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a810, %b811, %r809 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g813 = vector.transfer_read %sv794[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc814 = arith.mulf %r812, %g813 : vector<8x8xf64>
      vector.transfer_write %sc814, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d815 = vector.transfer_read %sv793[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2816 = vector.transfer_read %sv796[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t817 = arith.mulf %d815, %g2816 : vector<8x8xf64>
      %a818 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2819 = arith.addf %t817, %a818 : vector<8x8xf64>
      %b820 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl821 = arith.addf %t2819, %b820 : vector<8x8xf64>
      vector.transfer_write %fl821, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc822 = arith.constant dense<0.0> : vector<8x8xf64>
      %a823 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b824 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r825 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a823, %b824, %acc822 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a826 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b827 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r828 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a826, %b827, %r825 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r828, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc829 = arith.constant dense<0.0> : vector<8x8xf64>
      %a830 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b831 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r832 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a830, %b831, %acc829 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a833 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b834 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r835 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a833, %b834, %r832 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r835, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv836 = memref.subview %el2[0, 3, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv837 = memref.subview %el2[0, 4, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa838 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r839 = vector.transfer_read %sv836[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %m840 = arith.subf %r839, %fa838 : vector<8x8xf64>
      vector.transfer_write %m840, %sv836[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
      %fa841 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r842 = vector.transfer_read %sv837[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %p843 = arith.addf %r842, %fa841 : vector<8x8xf64>
      vector.transfer_write %p843, %sv837[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    }
    gpu.barrier
    %sv844 = memref.subview %v3634[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv845 = memref.subview %v3635[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv846 = memref.subview %el6[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv847 = memref.subview %el7[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv848 = memref.subview %el8[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc849 = arith.constant dense<0.0> : vector<8x8xf64>
      %a850 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b851 = vector.transfer_read %sv844[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r852 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a850, %b851, %acc849 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a853 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b854 = vector.transfer_read %sv844[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r855 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a853, %b854, %r852 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g856 = vector.transfer_read %sv847[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc857 = arith.mulf %r855, %g856 : vector<8x8xf64>
      vector.transfer_write %sc857, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc858 = arith.constant dense<0.0> : vector<8x8xf64>
      %a859 = vector.transfer_read %sv844[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b860 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r861 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a859, %b860, %acc858 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a862 = vector.transfer_read %sv844[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b863 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r864 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a862, %b863, %r861 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g865 = vector.transfer_read %sv846[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc866 = arith.mulf %r864, %g865 : vector<8x8xf64>
      vector.transfer_write %sc866, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d867 = vector.transfer_read %sv845[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2868 = vector.transfer_read %sv848[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t869 = arith.mulf %d867, %g2868 : vector<8x8xf64>
      %a870 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2871 = arith.addf %t869, %a870 : vector<8x8xf64>
      %b872 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl873 = arith.addf %t2871, %b872 : vector<8x8xf64>
      vector.transfer_write %fl873, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc874 = arith.constant dense<0.0> : vector<8x8xf64>
      %a875 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b876 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r877 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a875, %b876, %acc874 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a878 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b879 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r880 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a878, %b879, %r877 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r880, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc881 = arith.constant dense<0.0> : vector<8x8xf64>
      %a882 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b883 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r884 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a882, %b883, %acc881 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a885 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b886 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r887 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a885, %b886, %r884 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r887, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv888 = memref.subview %el2[0, 4, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv889 = memref.subview %el2[0, 5, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa890 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r891 = vector.transfer_read %sv888[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %m892 = arith.subf %r891, %fa890 : vector<8x8xf64>
      vector.transfer_write %m892, %sv888[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
      %fa893 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r894 = vector.transfer_read %sv889[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %p895 = arith.addf %r894, %fa893 : vector<8x8xf64>
      vector.transfer_write %p895, %sv889[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    }
    gpu.barrier
    %sv896 = memref.subview %v3634[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv897 = memref.subview %v3635[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv898 = memref.subview %el6[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv899 = memref.subview %el7[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv900 = memref.subview %el8[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc901 = arith.constant dense<0.0> : vector<8x8xf64>
      %a902 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b903 = vector.transfer_read %sv896[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r904 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a902, %b903, %acc901 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a905 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b906 = vector.transfer_read %sv896[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r907 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a905, %b906, %r904 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g908 = vector.transfer_read %sv899[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc909 = arith.mulf %r907, %g908 : vector<8x8xf64>
      vector.transfer_write %sc909, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc910 = arith.constant dense<0.0> : vector<8x8xf64>
      %a911 = vector.transfer_read %sv896[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b912 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r913 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a911, %b912, %acc910 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a914 = vector.transfer_read %sv896[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b915 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r916 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a914, %b915, %r913 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g917 = vector.transfer_read %sv898[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc918 = arith.mulf %r916, %g917 : vector<8x8xf64>
      vector.transfer_write %sc918, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d919 = vector.transfer_read %sv897[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2920 = vector.transfer_read %sv900[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t921 = arith.mulf %d919, %g2920 : vector<8x8xf64>
      %a922 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2923 = arith.addf %t921, %a922 : vector<8x8xf64>
      %b924 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl925 = arith.addf %t2923, %b924 : vector<8x8xf64>
      vector.transfer_write %fl925, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc926 = arith.constant dense<0.0> : vector<8x8xf64>
      %a927 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b928 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r929 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a927, %b928, %acc926 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a930 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b931 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r932 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a930, %b931, %r929 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r932, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc933 = arith.constant dense<0.0> : vector<8x8xf64>
      %a934 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b935 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r936 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a934, %b935, %acc933 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a937 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b938 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r939 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a937, %b938, %r936 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r939, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv940 = memref.subview %el2[0, 5, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv941 = memref.subview %el2[0, 6, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa942 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r943 = vector.transfer_read %sv940[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %m944 = arith.subf %r943, %fa942 : vector<8x8xf64>
      vector.transfer_write %m944, %sv940[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
      %fa945 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r946 = vector.transfer_read %sv941[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %p947 = arith.addf %r946, %fa945 : vector<8x8xf64>
      vector.transfer_write %p947, %sv941[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    }
    gpu.barrier
    %sv948 = memref.subview %v3634[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv949 = memref.subview %v3635[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv950 = memref.subview %el6[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv951 = memref.subview %el7[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv952 = memref.subview %el8[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc953 = arith.constant dense<0.0> : vector<8x8xf64>
      %a954 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b955 = vector.transfer_read %sv948[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r956 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a954, %b955, %acc953 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a957 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b958 = vector.transfer_read %sv948[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r959 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a957, %b958, %r956 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g960 = vector.transfer_read %sv951[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc961 = arith.mulf %r959, %g960 : vector<8x8xf64>
      vector.transfer_write %sc961, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc962 = arith.constant dense<0.0> : vector<8x8xf64>
      %a963 = vector.transfer_read %sv948[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b964 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r965 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a963, %b964, %acc962 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a966 = vector.transfer_read %sv948[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b967 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r968 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a966, %b967, %r965 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g969 = vector.transfer_read %sv950[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc970 = arith.mulf %r968, %g969 : vector<8x8xf64>
      vector.transfer_write %sc970, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d971 = vector.transfer_read %sv949[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g2972 = vector.transfer_read %sv952[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t973 = arith.mulf %d971, %g2972 : vector<8x8xf64>
      %a974 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t2975 = arith.addf %t973, %a974 : vector<8x8xf64>
      %b976 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl977 = arith.addf %t2975, %b976 : vector<8x8xf64>
      vector.transfer_write %fl977, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc978 = arith.constant dense<0.0> : vector<8x8xf64>
      %a979 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b980 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r981 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a979, %b980, %acc978 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a982 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b983 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r984 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a982, %b983, %r981 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r984, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc985 = arith.constant dense<0.0> : vector<8x8xf64>
      %a986 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b987 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r988 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a986, %b987, %acc985 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a989 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b990 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r991 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a989, %b990, %r988 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r991, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv992 = memref.subview %el2[0, 6, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    %sv993 = memref.subview %el2[0, 7, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa994 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r995 = vector.transfer_read %sv992[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %m996 = arith.subf %r995, %fa994 : vector<8x8xf64>
      vector.transfer_write %m996, %sv992[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
      %fa997 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r998 = vector.transfer_read %sv993[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %p999 = arith.addf %r998, %fa997 : vector<8x8xf64>
      vector.transfer_write %p999, %sv993[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    }
    gpu.barrier
    %sv1000 = memref.subview %el1[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1001 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1002 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1003 = vector.transfer_read %sv1000[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1004 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1002, %b1003, %acc1001 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1005 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1006 = vector.transfer_read %sv1000[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1007 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1005, %b1006, %r1004 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1007, %interp_all[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1008 = memref.subview %el1[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1009 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1010 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1011 = vector.transfer_read %sv1008[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1012 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1010, %b1011, %acc1009 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1013 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1014 = vector.transfer_read %sv1008[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1015 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1013, %b1014, %r1012 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1015, %interp_all[%c0, %c8] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1016 = memref.subview %el1[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1017 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1018 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1019 = vector.transfer_read %sv1016[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1020 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1018, %b1019, %acc1017 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1021 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1022 = vector.transfer_read %sv1016[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1023 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1021, %b1022, %r1020 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1023, %interp_all[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1024 = memref.subview %el1[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1025 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1026 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1027 = vector.transfer_read %sv1024[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1028 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1026, %b1027, %acc1025 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1029 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1030 = vector.transfer_read %sv1024[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1031 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1029, %b1030, %r1028 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1031, %interp_all[%c0, %c24] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1032 = memref.subview %el1[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1033 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1034 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1035 = vector.transfer_read %sv1032[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1036 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1034, %b1035, %acc1033 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1037 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1038 = vector.transfer_read %sv1032[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1039 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1037, %b1038, %r1036 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1039, %interp_all[%c0, %c32] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1040 = memref.subview %el1[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1041 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1042 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1043 = vector.transfer_read %sv1040[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1044 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1042, %b1043, %acc1041 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1045 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1046 = vector.transfer_read %sv1040[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1047 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1045, %b1046, %r1044 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1047, %interp_all[%c0, %c40] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1048 = memref.subview %el1[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1049 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1050 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1051 = vector.transfer_read %sv1048[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1052 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1050, %b1051, %acc1049 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1053 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1054 = vector.transfer_read %sv1048[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1055 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1053, %b1054, %r1052 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1055, %interp_all[%c0, %c48] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1056 = memref.subview %el1[7, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1057 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1058 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1059 = vector.transfer_read %sv1056[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1060 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1058, %b1059, %acc1057 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1061 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1062 = vector.transfer_read %sv1056[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1063 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1061, %b1062, %r1060 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1063, %interp_all[%c0, %c56] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1064 = memref.subview %el1[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1065 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1066 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1067 = vector.transfer_read %sv1064[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1068 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1066, %b1067, %acc1065 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1069 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1070 = vector.transfer_read %sv1064[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1071 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1069, %b1070, %r1068 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1071, %deriv_all[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1072 = memref.subview %el1[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1073 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1074 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1075 = vector.transfer_read %sv1072[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1076 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1074, %b1075, %acc1073 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1077 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1078 = vector.transfer_read %sv1072[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1079 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1077, %b1078, %r1076 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1079, %deriv_all[%c0, %c8] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1080 = memref.subview %el1[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1081 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1082 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1083 = vector.transfer_read %sv1080[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1084 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1082, %b1083, %acc1081 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1085 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1086 = vector.transfer_read %sv1080[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1087 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1085, %b1086, %r1084 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1087, %deriv_all[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1088 = memref.subview %el1[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1089 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1090 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1091 = vector.transfer_read %sv1088[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1092 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1090, %b1091, %acc1089 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1093 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1094 = vector.transfer_read %sv1088[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1095 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1093, %b1094, %r1092 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1095, %deriv_all[%c0, %c24] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1096 = memref.subview %el1[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1097 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1098 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1099 = vector.transfer_read %sv1096[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1100 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1098, %b1099, %acc1097 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1101 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1102 = vector.transfer_read %sv1096[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1103 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1101, %b1102, %r1100 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1103, %deriv_all[%c0, %c32] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1104 = memref.subview %el1[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1105 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1106 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1107 = vector.transfer_read %sv1104[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1108 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1106, %b1107, %acc1105 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1109 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1110 = vector.transfer_read %sv1104[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1111 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1109, %b1110, %r1108 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1111, %deriv_all[%c0, %c40] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1112 = memref.subview %el1[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1113 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1114 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1115 = vector.transfer_read %sv1112[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1116 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1114, %b1115, %acc1113 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1117 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1118 = vector.transfer_read %sv1112[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1119 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1117, %b1118, %r1116 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1119, %deriv_all[%c0, %c48] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1120 = memref.subview %el1[7, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1121 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1122 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1123 = vector.transfer_read %sv1120[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1124 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1122, %b1123, %acc1121 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1125 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1126 = vector.transfer_read %sv1120[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x4xf64>
      %r1127 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1125, %b1126, %r1124 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1127, %deriv_all[%c0, %c56] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %v31128 = memref.reinterpret_cast %interp_all to offset: [0], sizes: [8, 8, 8], strides: [64, 8, 1] : memref<8x64xf64, #gpu.address_space<workgroup>> to memref<8x8x8xf64, #gpu.address_space<workgroup>>
    %v31129 = memref.reinterpret_cast %deriv_all to offset: [0], sizes: [8, 8, 8], strides: [64, 8, 1] : memref<8x64xf64, #gpu.address_space<workgroup>> to memref<8x8x8xf64, #gpu.address_space<workgroup>>
    %sv1130 = memref.subview %v31128[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>
    %sv1131 = memref.subview %v31129[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>
    %sv1132 = memref.subview %el9[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv1133 = memref.subview %el10[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv1134 = memref.subview %el11[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1135 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1136 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1137 = vector.transfer_read %sv1130[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1138 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1136, %b1137, %acc1135 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1139 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1140 = vector.transfer_read %sv1130[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1141 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1139, %b1140, %r1138 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1142 = vector.transfer_read %sv1133[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc1143 = arith.mulf %r1141, %g1142 : vector<8x8xf64>
      vector.transfer_write %sc1143, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1144 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1145 = vector.transfer_read %sv1130[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1146 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1147 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1145, %b1146, %acc1144 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1148 = vector.transfer_read %sv1130[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1149 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1150 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1148, %b1149, %r1147 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1151 = vector.transfer_read %sv1132[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc1152 = arith.mulf %r1150, %g1151 : vector<8x8xf64>
      vector.transfer_write %sc1152, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d1153 = vector.transfer_read %sv1131[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g21154 = vector.transfer_read %sv1134[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t1155 = arith.mulf %d1153, %g21154 : vector<8x8xf64>
      %a1156 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t21157 = arith.addf %t1155, %a1156 : vector<8x8xf64>
      %b1158 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl1159 = arith.addf %t21157, %b1158 : vector<8x8xf64>
      vector.transfer_write %fl1159, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1160 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1161 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1162 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1163 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1161, %b1162, %acc1160 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1164 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1165 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1166 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1164, %b1165, %r1163 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1166, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1167 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1168 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1169 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1170 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1168, %b1169, %acc1167 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1171 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1172 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1173 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1171, %b1172, %r1170 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1173, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1174 = memref.subview %el2[0, 0, 0] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv1175 = memref.subview %el2[0, 0, 1] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1176 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1177 = vector.transfer_read %sv1174[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %m1178 = arith.subf %r1177, %fa1176 : vector<8x8xf64>
      vector.transfer_write %m1178, %sv1174[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
      %fa1179 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1180 = vector.transfer_read %sv1175[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %p1181 = arith.addf %r1180, %fa1179 : vector<8x8xf64>
      vector.transfer_write %p1181, %sv1175[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    }
    gpu.barrier
    %sv1182 = memref.subview %v31128[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv1183 = memref.subview %v31129[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv1184 = memref.subview %el9[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv1185 = memref.subview %el10[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv1186 = memref.subview %el11[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1187 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1188 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1189 = vector.transfer_read %sv1182[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1190 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1188, %b1189, %acc1187 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1191 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1192 = vector.transfer_read %sv1182[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1193 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1191, %b1192, %r1190 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1194 = vector.transfer_read %sv1185[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc1195 = arith.mulf %r1193, %g1194 : vector<8x8xf64>
      vector.transfer_write %sc1195, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1196 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1197 = vector.transfer_read %sv1182[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1198 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1199 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1197, %b1198, %acc1196 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1200 = vector.transfer_read %sv1182[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1201 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1202 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1200, %b1201, %r1199 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1203 = vector.transfer_read %sv1184[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc1204 = arith.mulf %r1202, %g1203 : vector<8x8xf64>
      vector.transfer_write %sc1204, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d1205 = vector.transfer_read %sv1183[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g21206 = vector.transfer_read %sv1186[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t1207 = arith.mulf %d1205, %g21206 : vector<8x8xf64>
      %a1208 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t21209 = arith.addf %t1207, %a1208 : vector<8x8xf64>
      %b1210 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl1211 = arith.addf %t21209, %b1210 : vector<8x8xf64>
      vector.transfer_write %fl1211, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1212 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1213 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1214 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1215 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1213, %b1214, %acc1212 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1216 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1217 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1218 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1216, %b1217, %r1215 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1218, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1219 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1220 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1221 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1222 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1220, %b1221, %acc1219 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1223 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1224 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1225 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1223, %b1224, %r1222 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1225, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1226 = memref.subview %el2[0, 0, 1] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv1227 = memref.subview %el2[0, 0, 2] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1228 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1229 = vector.transfer_read %sv1226[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %m1230 = arith.subf %r1229, %fa1228 : vector<8x8xf64>
      vector.transfer_write %m1230, %sv1226[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
      %fa1231 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1232 = vector.transfer_read %sv1227[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %p1233 = arith.addf %r1232, %fa1231 : vector<8x8xf64>
      vector.transfer_write %p1233, %sv1227[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    }
    gpu.barrier
    %sv1234 = memref.subview %v31128[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv1235 = memref.subview %v31129[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv1236 = memref.subview %el9[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv1237 = memref.subview %el10[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv1238 = memref.subview %el11[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1239 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1240 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1241 = vector.transfer_read %sv1234[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1242 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1240, %b1241, %acc1239 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1243 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1244 = vector.transfer_read %sv1234[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1245 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1243, %b1244, %r1242 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1246 = vector.transfer_read %sv1237[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc1247 = arith.mulf %r1245, %g1246 : vector<8x8xf64>
      vector.transfer_write %sc1247, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1248 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1249 = vector.transfer_read %sv1234[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1250 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1251 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1249, %b1250, %acc1248 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1252 = vector.transfer_read %sv1234[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1253 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1254 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1252, %b1253, %r1251 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1255 = vector.transfer_read %sv1236[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc1256 = arith.mulf %r1254, %g1255 : vector<8x8xf64>
      vector.transfer_write %sc1256, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d1257 = vector.transfer_read %sv1235[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g21258 = vector.transfer_read %sv1238[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t1259 = arith.mulf %d1257, %g21258 : vector<8x8xf64>
      %a1260 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t21261 = arith.addf %t1259, %a1260 : vector<8x8xf64>
      %b1262 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl1263 = arith.addf %t21261, %b1262 : vector<8x8xf64>
      vector.transfer_write %fl1263, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1264 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1265 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1266 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1267 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1265, %b1266, %acc1264 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1268 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1269 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1270 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1268, %b1269, %r1267 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1270, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1271 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1272 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1273 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1274 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1272, %b1273, %acc1271 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1275 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1276 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1277 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1275, %b1276, %r1274 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1277, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1278 = memref.subview %el2[0, 0, 2] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv1279 = memref.subview %el2[0, 0, 3] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1280 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1281 = vector.transfer_read %sv1278[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %m1282 = arith.subf %r1281, %fa1280 : vector<8x8xf64>
      vector.transfer_write %m1282, %sv1278[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
      %fa1283 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1284 = vector.transfer_read %sv1279[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %p1285 = arith.addf %r1284, %fa1283 : vector<8x8xf64>
      vector.transfer_write %p1285, %sv1279[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    }
    gpu.barrier
    %sv1286 = memref.subview %v31128[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv1287 = memref.subview %v31129[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv1288 = memref.subview %el9[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv1289 = memref.subview %el10[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv1290 = memref.subview %el11[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1291 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1292 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1293 = vector.transfer_read %sv1286[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1294 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1292, %b1293, %acc1291 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1295 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1296 = vector.transfer_read %sv1286[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1297 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1295, %b1296, %r1294 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1298 = vector.transfer_read %sv1289[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc1299 = arith.mulf %r1297, %g1298 : vector<8x8xf64>
      vector.transfer_write %sc1299, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1300 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1301 = vector.transfer_read %sv1286[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1302 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1303 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1301, %b1302, %acc1300 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1304 = vector.transfer_read %sv1286[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1305 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1306 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1304, %b1305, %r1303 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1307 = vector.transfer_read %sv1288[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc1308 = arith.mulf %r1306, %g1307 : vector<8x8xf64>
      vector.transfer_write %sc1308, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d1309 = vector.transfer_read %sv1287[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g21310 = vector.transfer_read %sv1290[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t1311 = arith.mulf %d1309, %g21310 : vector<8x8xf64>
      %a1312 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t21313 = arith.addf %t1311, %a1312 : vector<8x8xf64>
      %b1314 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl1315 = arith.addf %t21313, %b1314 : vector<8x8xf64>
      vector.transfer_write %fl1315, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1316 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1317 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1318 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1319 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1317, %b1318, %acc1316 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1320 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1321 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1322 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1320, %b1321, %r1319 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1322, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1323 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1324 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1325 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1326 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1324, %b1325, %acc1323 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1327 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1328 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1329 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1327, %b1328, %r1326 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1329, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1330 = memref.subview %el2[0, 0, 3] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv1331 = memref.subview %el2[0, 0, 4] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1332 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1333 = vector.transfer_read %sv1330[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %m1334 = arith.subf %r1333, %fa1332 : vector<8x8xf64>
      vector.transfer_write %m1334, %sv1330[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
      %fa1335 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1336 = vector.transfer_read %sv1331[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %p1337 = arith.addf %r1336, %fa1335 : vector<8x8xf64>
      vector.transfer_write %p1337, %sv1331[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    }
    gpu.barrier
    %sv1338 = memref.subview %v31128[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv1339 = memref.subview %v31129[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv1340 = memref.subview %el9[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv1341 = memref.subview %el10[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv1342 = memref.subview %el11[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1343 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1344 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1345 = vector.transfer_read %sv1338[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1346 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1344, %b1345, %acc1343 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1347 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1348 = vector.transfer_read %sv1338[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1349 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1347, %b1348, %r1346 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1350 = vector.transfer_read %sv1341[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc1351 = arith.mulf %r1349, %g1350 : vector<8x8xf64>
      vector.transfer_write %sc1351, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1352 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1353 = vector.transfer_read %sv1338[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1354 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1355 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1353, %b1354, %acc1352 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1356 = vector.transfer_read %sv1338[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1357 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1358 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1356, %b1357, %r1355 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1359 = vector.transfer_read %sv1340[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc1360 = arith.mulf %r1358, %g1359 : vector<8x8xf64>
      vector.transfer_write %sc1360, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d1361 = vector.transfer_read %sv1339[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g21362 = vector.transfer_read %sv1342[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t1363 = arith.mulf %d1361, %g21362 : vector<8x8xf64>
      %a1364 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t21365 = arith.addf %t1363, %a1364 : vector<8x8xf64>
      %b1366 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl1367 = arith.addf %t21365, %b1366 : vector<8x8xf64>
      vector.transfer_write %fl1367, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1368 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1369 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1370 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1371 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1369, %b1370, %acc1368 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1372 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1373 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1374 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1372, %b1373, %r1371 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1374, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1375 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1376 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1377 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1378 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1376, %b1377, %acc1375 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1379 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1380 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1381 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1379, %b1380, %r1378 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1381, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1382 = memref.subview %el2[0, 0, 4] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv1383 = memref.subview %el2[0, 0, 5] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1384 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1385 = vector.transfer_read %sv1382[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %m1386 = arith.subf %r1385, %fa1384 : vector<8x8xf64>
      vector.transfer_write %m1386, %sv1382[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
      %fa1387 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1388 = vector.transfer_read %sv1383[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %p1389 = arith.addf %r1388, %fa1387 : vector<8x8xf64>
      vector.transfer_write %p1389, %sv1383[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    }
    gpu.barrier
    %sv1390 = memref.subview %v31128[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv1391 = memref.subview %v31129[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv1392 = memref.subview %el9[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv1393 = memref.subview %el10[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv1394 = memref.subview %el11[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1395 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1396 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1397 = vector.transfer_read %sv1390[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1398 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1396, %b1397, %acc1395 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1399 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1400 = vector.transfer_read %sv1390[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1401 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1399, %b1400, %r1398 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1402 = vector.transfer_read %sv1393[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc1403 = arith.mulf %r1401, %g1402 : vector<8x8xf64>
      vector.transfer_write %sc1403, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1404 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1405 = vector.transfer_read %sv1390[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1406 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1407 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1405, %b1406, %acc1404 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1408 = vector.transfer_read %sv1390[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1409 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1410 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1408, %b1409, %r1407 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1411 = vector.transfer_read %sv1392[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc1412 = arith.mulf %r1410, %g1411 : vector<8x8xf64>
      vector.transfer_write %sc1412, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d1413 = vector.transfer_read %sv1391[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g21414 = vector.transfer_read %sv1394[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t1415 = arith.mulf %d1413, %g21414 : vector<8x8xf64>
      %a1416 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t21417 = arith.addf %t1415, %a1416 : vector<8x8xf64>
      %b1418 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl1419 = arith.addf %t21417, %b1418 : vector<8x8xf64>
      vector.transfer_write %fl1419, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1420 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1421 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1422 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1423 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1421, %b1422, %acc1420 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1424 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1425 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1426 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1424, %b1425, %r1423 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1426, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1427 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1428 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1429 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1430 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1428, %b1429, %acc1427 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1431 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1432 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1433 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1431, %b1432, %r1430 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1433, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1434 = memref.subview %el2[0, 0, 5] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv1435 = memref.subview %el2[0, 0, 6] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1436 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1437 = vector.transfer_read %sv1434[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %m1438 = arith.subf %r1437, %fa1436 : vector<8x8xf64>
      vector.transfer_write %m1438, %sv1434[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
      %fa1439 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1440 = vector.transfer_read %sv1435[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %p1441 = arith.addf %r1440, %fa1439 : vector<8x8xf64>
      vector.transfer_write %p1441, %sv1435[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    }
    gpu.barrier
    %sv1442 = memref.subview %v31128[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv1443 = memref.subview %v31129[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv1444 = memref.subview %el9[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv1445 = memref.subview %el10[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    %sv1446 = memref.subview %el11[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[8, 1], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1447 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1448 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1449 = vector.transfer_read %sv1442[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1450 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1448, %b1449, %acc1447 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1451 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1452 = vector.transfer_read %sv1442[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1453 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1451, %b1452, %r1450 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1454 = vector.transfer_read %sv1445[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc1455 = arith.mulf %r1453, %g1454 : vector<8x8xf64>
      vector.transfer_write %sc1455, %dt1g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1456 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1457 = vector.transfer_read %sv1442[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1458 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1459 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1457, %b1458, %acc1456 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1460 = vector.transfer_read %sv1442[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1461 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1462 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1460, %b1461, %r1459 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1463 = vector.transfer_read %sv1444[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %sc1464 = arith.mulf %r1462, %g1463 : vector<8x8xf64>
      vector.transfer_write %sc1464, %dt2g[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d1465 = vector.transfer_read %sv1443[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g21466 = vector.transfer_read %sv1446[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %t1467 = arith.mulf %d1465, %g21466 : vector<8x8xf64>
      %a1468 = vector.transfer_read %dt2g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %t21469 = arith.addf %t1467, %a1468 : vector<8x8xf64>
      %b1470 = vector.transfer_read %dt1g[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %fl1471 = arith.addf %t21469, %b1470 : vector<8x8xf64>
      vector.transfer_write %fl1471, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1472 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1473 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1474 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1475 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1473, %b1474, %acc1472 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1476 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1477 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1478 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1476, %b1477, %r1475 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1478, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1479 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1480 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1481 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1482 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1480, %b1481, %acc1479 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1483 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1484 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1485 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1483, %b1484, %r1482 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1485, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1486 = memref.subview %el2[0, 0, 6] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    %sv1487 = memref.subview %el2[0, 0, 7] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64, strided<[64, 8, 1], offset: ?>> to memref<8x8xf64, strided<[64, 8], offset: ?>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1488 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1489 = vector.transfer_read %sv1486[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %m1490 = arith.subf %r1489, %fa1488 : vector<8x8xf64>
      vector.transfer_write %m1490, %sv1486[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
      %fa1491 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1492 = vector.transfer_read %sv1487[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %p1493 = arith.addf %r1492, %fa1491 : vector<8x8xf64>
      vector.transfer_write %p1493, %sv1487[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    }
    gpu.barrier
    gpu.return
  }
}
