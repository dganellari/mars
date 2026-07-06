#a = affine_map<(m, n, k) -> (m, k)>
#b = affine_map<(m, n, k) -> (k, n)>
#bt = affine_map<(m, n, k) -> (n, k)>
#c = affine_map<(m, n, k) -> (m, n)>
gpu.module @mir_kernels [#nvvm.target<chip = "sm_90", O = 3>] {
  gpu.func @full(%u: memref<8x8x8xf64>, %Btil: memref<8x8xf64>, %Dtil: memref<8x8xf64>, %Dm: memref<8x8xf64>, %W: memref<8x8xf64>, %g00all: memref<8x8x8xf64>, %g01all: memref<8x8x8xf64>, %g02all: memref<8x8x8xf64>, %g10all: memref<8x8x8xf64>, %g11all: memref<8x8x8xf64>, %g12all: memref<8x8x8xf64>, %g20all: memref<8x8x8xf64>, %g21all: memref<8x8x8xf64>, %g22all: memref<8x8x8xf64>, %y3: memref<8x8x8xf64>) workgroup(%interp_all: memref<8x64xf64, #gpu.address_space<workgroup>>, %deriv_all: memref<8x64xf64, #gpu.address_space<workgroup>>, %flux: memref<8x8xf64, #gpu.address_space<workgroup>>, %tmp: memref<8x8xf64, #gpu.address_space<workgroup>>, %intf: memref<8x8xf64, #gpu.address_space<workgroup>>) kernel {
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
    %sv1 = memref.subview %u[0, 0, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 0>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc2 = arith.constant dense<0.0> : vector<8x8xf64>
      %a3 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b4 = vector.transfer_read %sv1[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 0>>, vector<4x8xf64>
      %r5 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a3, %b4, %acc2 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a6 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b7 = vector.transfer_read %sv1[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 0>>, vector<4x8xf64>
      %r8 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a6, %b7, %r5 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r8, %interp_all[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv9 = memref.subview %u[0, 1, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 8>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc10 = arith.constant dense<0.0> : vector<8x8xf64>
      %a11 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b12 = vector.transfer_read %sv9[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 8>>, vector<4x8xf64>
      %r13 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a11, %b12, %acc10 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a14 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b15 = vector.transfer_read %sv9[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 8>>, vector<4x8xf64>
      %r16 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a14, %b15, %r13 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r16, %interp_all[%c0, %c8] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv17 = memref.subview %u[0, 2, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 16>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc18 = arith.constant dense<0.0> : vector<8x8xf64>
      %a19 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b20 = vector.transfer_read %sv17[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 16>>, vector<4x8xf64>
      %r21 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a19, %b20, %acc18 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a22 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b23 = vector.transfer_read %sv17[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 16>>, vector<4x8xf64>
      %r24 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a22, %b23, %r21 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r24, %interp_all[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv25 = memref.subview %u[0, 3, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 24>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc26 = arith.constant dense<0.0> : vector<8x8xf64>
      %a27 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b28 = vector.transfer_read %sv25[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 24>>, vector<4x8xf64>
      %r29 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a27, %b28, %acc26 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a30 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b31 = vector.transfer_read %sv25[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 24>>, vector<4x8xf64>
      %r32 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a30, %b31, %r29 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r32, %interp_all[%c0, %c24] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv33 = memref.subview %u[0, 4, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 32>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc34 = arith.constant dense<0.0> : vector<8x8xf64>
      %a35 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b36 = vector.transfer_read %sv33[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 32>>, vector<4x8xf64>
      %r37 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a35, %b36, %acc34 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a38 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b39 = vector.transfer_read %sv33[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 32>>, vector<4x8xf64>
      %r40 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a38, %b39, %r37 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r40, %interp_all[%c0, %c32] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv41 = memref.subview %u[0, 5, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 40>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc42 = arith.constant dense<0.0> : vector<8x8xf64>
      %a43 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b44 = vector.transfer_read %sv41[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 40>>, vector<4x8xf64>
      %r45 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a43, %b44, %acc42 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a46 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b47 = vector.transfer_read %sv41[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 40>>, vector<4x8xf64>
      %r48 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a46, %b47, %r45 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r48, %interp_all[%c0, %c40] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv49 = memref.subview %u[0, 6, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 48>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc50 = arith.constant dense<0.0> : vector<8x8xf64>
      %a51 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b52 = vector.transfer_read %sv49[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 48>>, vector<4x8xf64>
      %r53 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a51, %b52, %acc50 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a54 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b55 = vector.transfer_read %sv49[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 48>>, vector<4x8xf64>
      %r56 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a54, %b55, %r53 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r56, %interp_all[%c0, %c48] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv57 = memref.subview %u[0, 7, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 56>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc58 = arith.constant dense<0.0> : vector<8x8xf64>
      %a59 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b60 = vector.transfer_read %sv57[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 56>>, vector<4x8xf64>
      %r61 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a59, %b60, %acc58 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a62 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b63 = vector.transfer_read %sv57[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 56>>, vector<4x8xf64>
      %r64 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a62, %b63, %r61 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r64, %interp_all[%c0, %c56] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv65 = memref.subview %u[0, 0, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 0>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc66 = arith.constant dense<0.0> : vector<8x8xf64>
      %a67 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b68 = vector.transfer_read %sv65[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 0>>, vector<4x8xf64>
      %r69 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a67, %b68, %acc66 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a70 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b71 = vector.transfer_read %sv65[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 0>>, vector<4x8xf64>
      %r72 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a70, %b71, %r69 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r72, %deriv_all[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv73 = memref.subview %u[0, 1, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 8>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc74 = arith.constant dense<0.0> : vector<8x8xf64>
      %a75 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b76 = vector.transfer_read %sv73[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 8>>, vector<4x8xf64>
      %r77 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a75, %b76, %acc74 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a78 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b79 = vector.transfer_read %sv73[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 8>>, vector<4x8xf64>
      %r80 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a78, %b79, %r77 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r80, %deriv_all[%c0, %c8] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv81 = memref.subview %u[0, 2, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 16>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc82 = arith.constant dense<0.0> : vector<8x8xf64>
      %a83 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b84 = vector.transfer_read %sv81[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 16>>, vector<4x8xf64>
      %r85 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a83, %b84, %acc82 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a86 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b87 = vector.transfer_read %sv81[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 16>>, vector<4x8xf64>
      %r88 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a86, %b87, %r85 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r88, %deriv_all[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv89 = memref.subview %u[0, 3, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 24>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc90 = arith.constant dense<0.0> : vector<8x8xf64>
      %a91 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b92 = vector.transfer_read %sv89[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 24>>, vector<4x8xf64>
      %r93 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a91, %b92, %acc90 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a94 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b95 = vector.transfer_read %sv89[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 24>>, vector<4x8xf64>
      %r96 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a94, %b95, %r93 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r96, %deriv_all[%c0, %c24] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv97 = memref.subview %u[0, 4, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 32>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc98 = arith.constant dense<0.0> : vector<8x8xf64>
      %a99 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b100 = vector.transfer_read %sv97[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 32>>, vector<4x8xf64>
      %r101 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a99, %b100, %acc98 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a102 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b103 = vector.transfer_read %sv97[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 32>>, vector<4x8xf64>
      %r104 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a102, %b103, %r101 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r104, %deriv_all[%c0, %c32] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv105 = memref.subview %u[0, 5, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 40>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc106 = arith.constant dense<0.0> : vector<8x8xf64>
      %a107 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b108 = vector.transfer_read %sv105[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 40>>, vector<4x8xf64>
      %r109 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a107, %b108, %acc106 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a110 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b111 = vector.transfer_read %sv105[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 40>>, vector<4x8xf64>
      %r112 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a110, %b111, %r109 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r112, %deriv_all[%c0, %c40] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv113 = memref.subview %u[0, 6, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 48>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc114 = arith.constant dense<0.0> : vector<8x8xf64>
      %a115 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b116 = vector.transfer_read %sv113[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 48>>, vector<4x8xf64>
      %r117 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a115, %b116, %acc114 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a118 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b119 = vector.transfer_read %sv113[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 48>>, vector<4x8xf64>
      %r120 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a118, %b119, %r117 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r120, %deriv_all[%c0, %c48] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv121 = memref.subview %u[0, 7, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 56>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc122 = arith.constant dense<0.0> : vector<8x8xf64>
      %a123 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b124 = vector.transfer_read %sv121[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 56>>, vector<4x8xf64>
      %r125 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a123, %b124, %acc122 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a126 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b127 = vector.transfer_read %sv121[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 56>>, vector<4x8xf64>
      %r128 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a126, %b127, %r125 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r128, %deriv_all[%c0, %c56] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %v3129 = memref.reinterpret_cast %interp_all to offset: [0], sizes: [8, 8, 8], strides: [64, 8, 1] : memref<8x64xf64, #gpu.address_space<workgroup>> to memref<8x8x8xf64, #gpu.address_space<workgroup>>
    %v3130 = memref.reinterpret_cast %deriv_all to offset: [0], sizes: [8, 8, 8], strides: [64, 8, 1] : memref<8x64xf64, #gpu.address_space<workgroup>> to memref<8x8x8xf64, #gpu.address_space<workgroup>>
    %sv131 = memref.subview %v3129[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>
    %sv132 = memref.subview %v3130[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>
    %sv133 = memref.subview %g00all[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 0>>
    %sv134 = memref.subview %g01all[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 0>>
    %sv135 = memref.subview %g02all[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 0>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc136 = arith.constant dense<0.0> : vector<8x8xf64>
      %a137 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b138 = vector.transfer_read %sv131[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r139 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a137, %b138, %acc136 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a140 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b141 = vector.transfer_read %sv131[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r142 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a140, %b141, %r139 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g143 = vector.transfer_read %sv134[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>>, vector<8x8xf64>
      %dt144 = arith.mulf %r142, %g143 : vector<8x8xf64>
      %acc145 = arith.constant dense<0.0> : vector<8x8xf64>
      %a146 = vector.transfer_read %sv131[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b147 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r148 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a146, %b147, %acc145 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a149 = vector.transfer_read %sv131[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b150 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r151 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a149, %b150, %r148 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g152 = vector.transfer_read %sv133[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>>, vector<8x8xf64>
      %dt153 = arith.mulf %r151, %g152 : vector<8x8xf64>
      %dv154 = vector.transfer_read %sv132[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g155 = vector.transfer_read %sv135[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>>, vector<8x8xf64>
      %t156 = arith.mulf %dv154, %g155 : vector<8x8xf64>
      %t2157 = arith.addf %t156, %dt153 : vector<8x8xf64>
      %fl158 = arith.addf %t2157, %dt144 : vector<8x8xf64>
      vector.transfer_write %fl158, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc159 = arith.constant dense<0.0> : vector<8x8xf64>
      %a160 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b161 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r162 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a160, %b161, %acc159 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a163 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b164 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r165 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a163, %b164, %r162 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r165, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc166 = arith.constant dense<0.0> : vector<8x8xf64>
      %a167 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b168 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r169 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a167, %b168, %acc166 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a170 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b171 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r172 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a170, %b171, %r169 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r172, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv173 = memref.subview %y3[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 0>>
    %sv174 = memref.subview %y3[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 64>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa175 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r176 = vector.transfer_read %sv173[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>>, vector<8x8xf64>
      %m177 = arith.subf %r176, %fa175 : vector<8x8xf64>
      vector.transfer_write %m177, %sv173[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 0>>
      %r178 = vector.transfer_read %sv174[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<8x8xf64>
      %p179 = arith.addf %r178, %fa175 : vector<8x8xf64>
      vector.transfer_write %p179, %sv174[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 64>>
    }
    gpu.barrier
    %sv180 = memref.subview %v3129[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv181 = memref.subview %v3130[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv182 = memref.subview %g00all[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 64>>
    %sv183 = memref.subview %g01all[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 64>>
    %sv184 = memref.subview %g02all[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 64>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc185 = arith.constant dense<0.0> : vector<8x8xf64>
      %a186 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b187 = vector.transfer_read %sv180[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r188 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a186, %b187, %acc185 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a189 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b190 = vector.transfer_read %sv180[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r191 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a189, %b190, %r188 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g192 = vector.transfer_read %sv183[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<8x8xf64>
      %dt193 = arith.mulf %r191, %g192 : vector<8x8xf64>
      %acc194 = arith.constant dense<0.0> : vector<8x8xf64>
      %a195 = vector.transfer_read %sv180[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b196 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r197 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a195, %b196, %acc194 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a198 = vector.transfer_read %sv180[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b199 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r200 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a198, %b199, %r197 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g201 = vector.transfer_read %sv182[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<8x8xf64>
      %dt202 = arith.mulf %r200, %g201 : vector<8x8xf64>
      %dv203 = vector.transfer_read %sv181[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g204 = vector.transfer_read %sv184[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<8x8xf64>
      %t205 = arith.mulf %dv203, %g204 : vector<8x8xf64>
      %t2206 = arith.addf %t205, %dt202 : vector<8x8xf64>
      %fl207 = arith.addf %t2206, %dt193 : vector<8x8xf64>
      vector.transfer_write %fl207, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc208 = arith.constant dense<0.0> : vector<8x8xf64>
      %a209 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b210 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r211 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a209, %b210, %acc208 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a212 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b213 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r214 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a212, %b213, %r211 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r214, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc215 = arith.constant dense<0.0> : vector<8x8xf64>
      %a216 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b217 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r218 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a216, %b217, %acc215 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a219 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b220 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r221 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a219, %b220, %r218 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r221, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv222 = memref.subview %y3[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 64>>
    %sv223 = memref.subview %y3[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 128>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa224 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r225 = vector.transfer_read %sv222[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<8x8xf64>
      %m226 = arith.subf %r225, %fa224 : vector<8x8xf64>
      vector.transfer_write %m226, %sv222[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 64>>
      %r227 = vector.transfer_read %sv223[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<8x8xf64>
      %p228 = arith.addf %r227, %fa224 : vector<8x8xf64>
      vector.transfer_write %p228, %sv223[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 128>>
    }
    gpu.barrier
    %sv229 = memref.subview %v3129[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv230 = memref.subview %v3130[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv231 = memref.subview %g00all[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 128>>
    %sv232 = memref.subview %g01all[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 128>>
    %sv233 = memref.subview %g02all[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 128>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc234 = arith.constant dense<0.0> : vector<8x8xf64>
      %a235 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b236 = vector.transfer_read %sv229[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r237 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a235, %b236, %acc234 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a238 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b239 = vector.transfer_read %sv229[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r240 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a238, %b239, %r237 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g241 = vector.transfer_read %sv232[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<8x8xf64>
      %dt242 = arith.mulf %r240, %g241 : vector<8x8xf64>
      %acc243 = arith.constant dense<0.0> : vector<8x8xf64>
      %a244 = vector.transfer_read %sv229[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b245 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r246 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a244, %b245, %acc243 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a247 = vector.transfer_read %sv229[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b248 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r249 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a247, %b248, %r246 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g250 = vector.transfer_read %sv231[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<8x8xf64>
      %dt251 = arith.mulf %r249, %g250 : vector<8x8xf64>
      %dv252 = vector.transfer_read %sv230[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g253 = vector.transfer_read %sv233[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<8x8xf64>
      %t254 = arith.mulf %dv252, %g253 : vector<8x8xf64>
      %t2255 = arith.addf %t254, %dt251 : vector<8x8xf64>
      %fl256 = arith.addf %t2255, %dt242 : vector<8x8xf64>
      vector.transfer_write %fl256, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc257 = arith.constant dense<0.0> : vector<8x8xf64>
      %a258 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b259 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r260 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a258, %b259, %acc257 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a261 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b262 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r263 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a261, %b262, %r260 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r263, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc264 = arith.constant dense<0.0> : vector<8x8xf64>
      %a265 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b266 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r267 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a265, %b266, %acc264 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a268 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b269 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r270 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a268, %b269, %r267 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r270, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv271 = memref.subview %y3[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 128>>
    %sv272 = memref.subview %y3[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 192>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa273 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r274 = vector.transfer_read %sv271[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<8x8xf64>
      %m275 = arith.subf %r274, %fa273 : vector<8x8xf64>
      vector.transfer_write %m275, %sv271[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 128>>
      %r276 = vector.transfer_read %sv272[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<8x8xf64>
      %p277 = arith.addf %r276, %fa273 : vector<8x8xf64>
      vector.transfer_write %p277, %sv272[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 192>>
    }
    gpu.barrier
    %sv278 = memref.subview %v3129[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv279 = memref.subview %v3130[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv280 = memref.subview %g00all[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 192>>
    %sv281 = memref.subview %g01all[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 192>>
    %sv282 = memref.subview %g02all[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 192>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc283 = arith.constant dense<0.0> : vector<8x8xf64>
      %a284 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b285 = vector.transfer_read %sv278[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r286 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a284, %b285, %acc283 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a287 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b288 = vector.transfer_read %sv278[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r289 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a287, %b288, %r286 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g290 = vector.transfer_read %sv281[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<8x8xf64>
      %dt291 = arith.mulf %r289, %g290 : vector<8x8xf64>
      %acc292 = arith.constant dense<0.0> : vector<8x8xf64>
      %a293 = vector.transfer_read %sv278[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b294 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r295 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a293, %b294, %acc292 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a296 = vector.transfer_read %sv278[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b297 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r298 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a296, %b297, %r295 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g299 = vector.transfer_read %sv280[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<8x8xf64>
      %dt300 = arith.mulf %r298, %g299 : vector<8x8xf64>
      %dv301 = vector.transfer_read %sv279[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g302 = vector.transfer_read %sv282[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<8x8xf64>
      %t303 = arith.mulf %dv301, %g302 : vector<8x8xf64>
      %t2304 = arith.addf %t303, %dt300 : vector<8x8xf64>
      %fl305 = arith.addf %t2304, %dt291 : vector<8x8xf64>
      vector.transfer_write %fl305, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc306 = arith.constant dense<0.0> : vector<8x8xf64>
      %a307 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b308 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r309 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a307, %b308, %acc306 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a310 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b311 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r312 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a310, %b311, %r309 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r312, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc313 = arith.constant dense<0.0> : vector<8x8xf64>
      %a314 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b315 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r316 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a314, %b315, %acc313 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a317 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b318 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r319 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a317, %b318, %r316 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r319, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv320 = memref.subview %y3[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 192>>
    %sv321 = memref.subview %y3[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 256>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa322 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r323 = vector.transfer_read %sv320[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<8x8xf64>
      %m324 = arith.subf %r323, %fa322 : vector<8x8xf64>
      vector.transfer_write %m324, %sv320[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 192>>
      %r325 = vector.transfer_read %sv321[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<8x8xf64>
      %p326 = arith.addf %r325, %fa322 : vector<8x8xf64>
      vector.transfer_write %p326, %sv321[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 256>>
    }
    gpu.barrier
    %sv327 = memref.subview %v3129[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv328 = memref.subview %v3130[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv329 = memref.subview %g00all[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 256>>
    %sv330 = memref.subview %g01all[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 256>>
    %sv331 = memref.subview %g02all[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 256>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc332 = arith.constant dense<0.0> : vector<8x8xf64>
      %a333 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b334 = vector.transfer_read %sv327[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r335 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a333, %b334, %acc332 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a336 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b337 = vector.transfer_read %sv327[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r338 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a336, %b337, %r335 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g339 = vector.transfer_read %sv330[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<8x8xf64>
      %dt340 = arith.mulf %r338, %g339 : vector<8x8xf64>
      %acc341 = arith.constant dense<0.0> : vector<8x8xf64>
      %a342 = vector.transfer_read %sv327[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b343 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r344 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a342, %b343, %acc341 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a345 = vector.transfer_read %sv327[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b346 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r347 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a345, %b346, %r344 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g348 = vector.transfer_read %sv329[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<8x8xf64>
      %dt349 = arith.mulf %r347, %g348 : vector<8x8xf64>
      %dv350 = vector.transfer_read %sv328[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g351 = vector.transfer_read %sv331[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<8x8xf64>
      %t352 = arith.mulf %dv350, %g351 : vector<8x8xf64>
      %t2353 = arith.addf %t352, %dt349 : vector<8x8xf64>
      %fl354 = arith.addf %t2353, %dt340 : vector<8x8xf64>
      vector.transfer_write %fl354, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc355 = arith.constant dense<0.0> : vector<8x8xf64>
      %a356 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b357 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r358 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a356, %b357, %acc355 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a359 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b360 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r361 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a359, %b360, %r358 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r361, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc362 = arith.constant dense<0.0> : vector<8x8xf64>
      %a363 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b364 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r365 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a363, %b364, %acc362 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a366 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b367 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r368 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a366, %b367, %r365 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r368, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv369 = memref.subview %y3[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 256>>
    %sv370 = memref.subview %y3[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 320>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa371 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r372 = vector.transfer_read %sv369[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<8x8xf64>
      %m373 = arith.subf %r372, %fa371 : vector<8x8xf64>
      vector.transfer_write %m373, %sv369[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 256>>
      %r374 = vector.transfer_read %sv370[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<8x8xf64>
      %p375 = arith.addf %r374, %fa371 : vector<8x8xf64>
      vector.transfer_write %p375, %sv370[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 320>>
    }
    gpu.barrier
    %sv376 = memref.subview %v3129[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv377 = memref.subview %v3130[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv378 = memref.subview %g00all[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 320>>
    %sv379 = memref.subview %g01all[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 320>>
    %sv380 = memref.subview %g02all[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 320>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc381 = arith.constant dense<0.0> : vector<8x8xf64>
      %a382 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b383 = vector.transfer_read %sv376[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r384 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a382, %b383, %acc381 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a385 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b386 = vector.transfer_read %sv376[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r387 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a385, %b386, %r384 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g388 = vector.transfer_read %sv379[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<8x8xf64>
      %dt389 = arith.mulf %r387, %g388 : vector<8x8xf64>
      %acc390 = arith.constant dense<0.0> : vector<8x8xf64>
      %a391 = vector.transfer_read %sv376[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b392 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r393 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a391, %b392, %acc390 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a394 = vector.transfer_read %sv376[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b395 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r396 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a394, %b395, %r393 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g397 = vector.transfer_read %sv378[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<8x8xf64>
      %dt398 = arith.mulf %r396, %g397 : vector<8x8xf64>
      %dv399 = vector.transfer_read %sv377[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g400 = vector.transfer_read %sv380[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<8x8xf64>
      %t401 = arith.mulf %dv399, %g400 : vector<8x8xf64>
      %t2402 = arith.addf %t401, %dt398 : vector<8x8xf64>
      %fl403 = arith.addf %t2402, %dt389 : vector<8x8xf64>
      vector.transfer_write %fl403, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc404 = arith.constant dense<0.0> : vector<8x8xf64>
      %a405 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b406 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r407 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a405, %b406, %acc404 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a408 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b409 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r410 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a408, %b409, %r407 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r410, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc411 = arith.constant dense<0.0> : vector<8x8xf64>
      %a412 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b413 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r414 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a412, %b413, %acc411 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a415 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b416 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r417 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a415, %b416, %r414 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r417, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv418 = memref.subview %y3[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 320>>
    %sv419 = memref.subview %y3[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 384>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa420 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r421 = vector.transfer_read %sv418[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<8x8xf64>
      %m422 = arith.subf %r421, %fa420 : vector<8x8xf64>
      vector.transfer_write %m422, %sv418[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 320>>
      %r423 = vector.transfer_read %sv419[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<8x8xf64>
      %p424 = arith.addf %r423, %fa420 : vector<8x8xf64>
      vector.transfer_write %p424, %sv419[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 384>>
    }
    gpu.barrier
    %sv425 = memref.subview %v3129[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv426 = memref.subview %v3130[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv427 = memref.subview %g00all[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 384>>
    %sv428 = memref.subview %g01all[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 384>>
    %sv429 = memref.subview %g02all[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 384>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc430 = arith.constant dense<0.0> : vector<8x8xf64>
      %a431 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b432 = vector.transfer_read %sv425[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r433 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a431, %b432, %acc430 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a434 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b435 = vector.transfer_read %sv425[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r436 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a434, %b435, %r433 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g437 = vector.transfer_read %sv428[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<8x8xf64>
      %dt438 = arith.mulf %r436, %g437 : vector<8x8xf64>
      %acc439 = arith.constant dense<0.0> : vector<8x8xf64>
      %a440 = vector.transfer_read %sv425[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b441 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r442 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a440, %b441, %acc439 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a443 = vector.transfer_read %sv425[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b444 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r445 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a443, %b444, %r442 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g446 = vector.transfer_read %sv427[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<8x8xf64>
      %dt447 = arith.mulf %r445, %g446 : vector<8x8xf64>
      %dv448 = vector.transfer_read %sv426[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g449 = vector.transfer_read %sv429[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<8x8xf64>
      %t450 = arith.mulf %dv448, %g449 : vector<8x8xf64>
      %t2451 = arith.addf %t450, %dt447 : vector<8x8xf64>
      %fl452 = arith.addf %t2451, %dt438 : vector<8x8xf64>
      vector.transfer_write %fl452, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc453 = arith.constant dense<0.0> : vector<8x8xf64>
      %a454 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b455 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r456 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a454, %b455, %acc453 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a457 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b458 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r459 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a457, %b458, %r456 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r459, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc460 = arith.constant dense<0.0> : vector<8x8xf64>
      %a461 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b462 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r463 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a461, %b462, %acc460 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a464 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b465 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r466 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a464, %b465, %r463 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r466, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv467 = memref.subview %y3[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 384>>
    %sv468 = memref.subview %y3[7, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 448>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa469 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r470 = vector.transfer_read %sv467[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<8x8xf64>
      %m471 = arith.subf %r470, %fa469 : vector<8x8xf64>
      vector.transfer_write %m471, %sv467[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 384>>
      %r472 = vector.transfer_read %sv468[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 448>>, vector<8x8xf64>
      %p473 = arith.addf %r472, %fa469 : vector<8x8xf64>
      vector.transfer_write %p473, %sv468[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: 448>>
    }
    gpu.barrier
    %sv474 = memref.subview %u[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 0>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc475 = arith.constant dense<0.0> : vector<8x8xf64>
      %a476 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b477 = vector.transfer_read %sv474[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>>, vector<4x8xf64>
      %r478 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a476, %b477, %acc475 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a479 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b480 = vector.transfer_read %sv474[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>>, vector<4x8xf64>
      %r481 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a479, %b480, %r478 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r481, %interp_all[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv482 = memref.subview %u[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 64>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc483 = arith.constant dense<0.0> : vector<8x8xf64>
      %a484 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b485 = vector.transfer_read %sv482[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<4x8xf64>
      %r486 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a484, %b485, %acc483 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a487 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b488 = vector.transfer_read %sv482[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<4x8xf64>
      %r489 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a487, %b488, %r486 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r489, %interp_all[%c0, %c8] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv490 = memref.subview %u[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 128>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc491 = arith.constant dense<0.0> : vector<8x8xf64>
      %a492 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b493 = vector.transfer_read %sv490[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<4x8xf64>
      %r494 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a492, %b493, %acc491 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a495 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b496 = vector.transfer_read %sv490[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<4x8xf64>
      %r497 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a495, %b496, %r494 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r497, %interp_all[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv498 = memref.subview %u[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 192>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc499 = arith.constant dense<0.0> : vector<8x8xf64>
      %a500 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b501 = vector.transfer_read %sv498[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<4x8xf64>
      %r502 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a500, %b501, %acc499 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a503 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b504 = vector.transfer_read %sv498[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<4x8xf64>
      %r505 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a503, %b504, %r502 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r505, %interp_all[%c0, %c24] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv506 = memref.subview %u[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 256>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc507 = arith.constant dense<0.0> : vector<8x8xf64>
      %a508 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b509 = vector.transfer_read %sv506[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<4x8xf64>
      %r510 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a508, %b509, %acc507 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a511 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b512 = vector.transfer_read %sv506[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<4x8xf64>
      %r513 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a511, %b512, %r510 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r513, %interp_all[%c0, %c32] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv514 = memref.subview %u[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 320>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc515 = arith.constant dense<0.0> : vector<8x8xf64>
      %a516 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b517 = vector.transfer_read %sv514[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<4x8xf64>
      %r518 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a516, %b517, %acc515 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a519 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b520 = vector.transfer_read %sv514[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<4x8xf64>
      %r521 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a519, %b520, %r518 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r521, %interp_all[%c0, %c40] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv522 = memref.subview %u[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 384>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc523 = arith.constant dense<0.0> : vector<8x8xf64>
      %a524 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b525 = vector.transfer_read %sv522[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<4x8xf64>
      %r526 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a524, %b525, %acc523 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a527 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b528 = vector.transfer_read %sv522[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<4x8xf64>
      %r529 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a527, %b528, %r526 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r529, %interp_all[%c0, %c48] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv530 = memref.subview %u[7, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 448>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc531 = arith.constant dense<0.0> : vector<8x8xf64>
      %a532 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b533 = vector.transfer_read %sv530[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 448>>, vector<4x8xf64>
      %r534 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a532, %b533, %acc531 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a535 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b536 = vector.transfer_read %sv530[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 448>>, vector<4x8xf64>
      %r537 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a535, %b536, %r534 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r537, %interp_all[%c0, %c56] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv538 = memref.subview %u[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 0>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc539 = arith.constant dense<0.0> : vector<8x8xf64>
      %a540 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b541 = vector.transfer_read %sv538[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>>, vector<4x8xf64>
      %r542 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a540, %b541, %acc539 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a543 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b544 = vector.transfer_read %sv538[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>>, vector<4x8xf64>
      %r545 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a543, %b544, %r542 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r545, %deriv_all[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv546 = memref.subview %u[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 64>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc547 = arith.constant dense<0.0> : vector<8x8xf64>
      %a548 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b549 = vector.transfer_read %sv546[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<4x8xf64>
      %r550 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a548, %b549, %acc547 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a551 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b552 = vector.transfer_read %sv546[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<4x8xf64>
      %r553 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a551, %b552, %r550 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r553, %deriv_all[%c0, %c8] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv554 = memref.subview %u[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 128>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc555 = arith.constant dense<0.0> : vector<8x8xf64>
      %a556 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b557 = vector.transfer_read %sv554[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<4x8xf64>
      %r558 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a556, %b557, %acc555 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a559 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b560 = vector.transfer_read %sv554[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<4x8xf64>
      %r561 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a559, %b560, %r558 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r561, %deriv_all[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv562 = memref.subview %u[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 192>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc563 = arith.constant dense<0.0> : vector<8x8xf64>
      %a564 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b565 = vector.transfer_read %sv562[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<4x8xf64>
      %r566 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a564, %b565, %acc563 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a567 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b568 = vector.transfer_read %sv562[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<4x8xf64>
      %r569 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a567, %b568, %r566 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r569, %deriv_all[%c0, %c24] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv570 = memref.subview %u[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 256>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc571 = arith.constant dense<0.0> : vector<8x8xf64>
      %a572 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b573 = vector.transfer_read %sv570[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<4x8xf64>
      %r574 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a572, %b573, %acc571 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a575 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b576 = vector.transfer_read %sv570[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<4x8xf64>
      %r577 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a575, %b576, %r574 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r577, %deriv_all[%c0, %c32] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv578 = memref.subview %u[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 320>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc579 = arith.constant dense<0.0> : vector<8x8xf64>
      %a580 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b581 = vector.transfer_read %sv578[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<4x8xf64>
      %r582 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a580, %b581, %acc579 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a583 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b584 = vector.transfer_read %sv578[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<4x8xf64>
      %r585 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a583, %b584, %r582 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r585, %deriv_all[%c0, %c40] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv586 = memref.subview %u[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 384>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc587 = arith.constant dense<0.0> : vector<8x8xf64>
      %a588 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b589 = vector.transfer_read %sv586[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<4x8xf64>
      %r590 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a588, %b589, %acc587 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a591 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b592 = vector.transfer_read %sv586[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<4x8xf64>
      %r593 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a591, %b592, %r590 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r593, %deriv_all[%c0, %c48] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv594 = memref.subview %u[7, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 448>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc595 = arith.constant dense<0.0> : vector<8x8xf64>
      %a596 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b597 = vector.transfer_read %sv594[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 448>>, vector<4x8xf64>
      %r598 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a596, %b597, %acc595 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a599 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b600 = vector.transfer_read %sv594[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 448>>, vector<4x8xf64>
      %r601 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a599, %b600, %r598 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r601, %deriv_all[%c0, %c56] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %v3602 = memref.reinterpret_cast %interp_all to offset: [0], sizes: [8, 8, 8], strides: [64, 8, 1] : memref<8x64xf64, #gpu.address_space<workgroup>> to memref<8x8x8xf64, #gpu.address_space<workgroup>>
    %v3603 = memref.reinterpret_cast %deriv_all to offset: [0], sizes: [8, 8, 8], strides: [64, 8, 1] : memref<8x64xf64, #gpu.address_space<workgroup>> to memref<8x8x8xf64, #gpu.address_space<workgroup>>
    %sv604 = memref.subview %v3602[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>
    %sv605 = memref.subview %v3603[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>
    %sv606 = memref.subview %g10all[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 0>>
    %sv607 = memref.subview %g11all[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 0>>
    %sv608 = memref.subview %g12all[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 0>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc609 = arith.constant dense<0.0> : vector<8x8xf64>
      %a610 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b611 = vector.transfer_read %sv604[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r612 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a610, %b611, %acc609 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a613 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b614 = vector.transfer_read %sv604[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r615 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a613, %b614, %r612 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g616 = vector.transfer_read %sv607[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>>, vector<8x8xf64>
      %dt617 = arith.mulf %r615, %g616 : vector<8x8xf64>
      %acc618 = arith.constant dense<0.0> : vector<8x8xf64>
      %a619 = vector.transfer_read %sv604[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b620 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r621 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a619, %b620, %acc618 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a622 = vector.transfer_read %sv604[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b623 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r624 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a622, %b623, %r621 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g625 = vector.transfer_read %sv606[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>>, vector<8x8xf64>
      %dt626 = arith.mulf %r624, %g625 : vector<8x8xf64>
      %dv627 = vector.transfer_read %sv605[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g628 = vector.transfer_read %sv608[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>>, vector<8x8xf64>
      %t629 = arith.mulf %dv627, %g628 : vector<8x8xf64>
      %t2630 = arith.addf %t629, %dt626 : vector<8x8xf64>
      %fl631 = arith.addf %t2630, %dt617 : vector<8x8xf64>
      vector.transfer_write %fl631, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc632 = arith.constant dense<0.0> : vector<8x8xf64>
      %a633 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b634 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r635 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a633, %b634, %acc632 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a636 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b637 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r638 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a636, %b637, %r635 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r638, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc639 = arith.constant dense<0.0> : vector<8x8xf64>
      %a640 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b641 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r642 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a640, %b641, %acc639 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a643 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b644 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r645 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a643, %b644, %r642 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r645, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv646 = memref.subview %y3[0, 0, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 0>>
    %sv647 = memref.subview %y3[0, 1, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 8>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa648 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r649 = vector.transfer_read %sv646[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 0>>, vector<8x8xf64>
      %m650 = arith.subf %r649, %fa648 : vector<8x8xf64>
      vector.transfer_write %m650, %sv646[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: 0>>
      %r651 = vector.transfer_read %sv647[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 8>>, vector<8x8xf64>
      %p652 = arith.addf %r651, %fa648 : vector<8x8xf64>
      vector.transfer_write %p652, %sv647[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: 8>>
    }
    gpu.barrier
    %sv653 = memref.subview %v3602[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv654 = memref.subview %v3603[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv655 = memref.subview %g10all[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 64>>
    %sv656 = memref.subview %g11all[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 64>>
    %sv657 = memref.subview %g12all[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 64>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc658 = arith.constant dense<0.0> : vector<8x8xf64>
      %a659 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b660 = vector.transfer_read %sv653[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r661 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a659, %b660, %acc658 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a662 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b663 = vector.transfer_read %sv653[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r664 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a662, %b663, %r661 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g665 = vector.transfer_read %sv656[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<8x8xf64>
      %dt666 = arith.mulf %r664, %g665 : vector<8x8xf64>
      %acc667 = arith.constant dense<0.0> : vector<8x8xf64>
      %a668 = vector.transfer_read %sv653[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b669 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r670 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a668, %b669, %acc667 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a671 = vector.transfer_read %sv653[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b672 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r673 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a671, %b672, %r670 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g674 = vector.transfer_read %sv655[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<8x8xf64>
      %dt675 = arith.mulf %r673, %g674 : vector<8x8xf64>
      %dv676 = vector.transfer_read %sv654[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g677 = vector.transfer_read %sv657[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<8x8xf64>
      %t678 = arith.mulf %dv676, %g677 : vector<8x8xf64>
      %t2679 = arith.addf %t678, %dt675 : vector<8x8xf64>
      %fl680 = arith.addf %t2679, %dt666 : vector<8x8xf64>
      vector.transfer_write %fl680, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc681 = arith.constant dense<0.0> : vector<8x8xf64>
      %a682 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b683 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r684 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a682, %b683, %acc681 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a685 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b686 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r687 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a685, %b686, %r684 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r687, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc688 = arith.constant dense<0.0> : vector<8x8xf64>
      %a689 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b690 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r691 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a689, %b690, %acc688 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a692 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b693 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r694 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a692, %b693, %r691 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r694, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv695 = memref.subview %y3[0, 1, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 8>>
    %sv696 = memref.subview %y3[0, 2, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 16>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa697 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r698 = vector.transfer_read %sv695[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 8>>, vector<8x8xf64>
      %m699 = arith.subf %r698, %fa697 : vector<8x8xf64>
      vector.transfer_write %m699, %sv695[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: 8>>
      %r700 = vector.transfer_read %sv696[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 16>>, vector<8x8xf64>
      %p701 = arith.addf %r700, %fa697 : vector<8x8xf64>
      vector.transfer_write %p701, %sv696[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: 16>>
    }
    gpu.barrier
    %sv702 = memref.subview %v3602[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv703 = memref.subview %v3603[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv704 = memref.subview %g10all[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 128>>
    %sv705 = memref.subview %g11all[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 128>>
    %sv706 = memref.subview %g12all[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 128>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc707 = arith.constant dense<0.0> : vector<8x8xf64>
      %a708 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b709 = vector.transfer_read %sv702[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r710 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a708, %b709, %acc707 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a711 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b712 = vector.transfer_read %sv702[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r713 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a711, %b712, %r710 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g714 = vector.transfer_read %sv705[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<8x8xf64>
      %dt715 = arith.mulf %r713, %g714 : vector<8x8xf64>
      %acc716 = arith.constant dense<0.0> : vector<8x8xf64>
      %a717 = vector.transfer_read %sv702[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b718 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r719 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a717, %b718, %acc716 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a720 = vector.transfer_read %sv702[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b721 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r722 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a720, %b721, %r719 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g723 = vector.transfer_read %sv704[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<8x8xf64>
      %dt724 = arith.mulf %r722, %g723 : vector<8x8xf64>
      %dv725 = vector.transfer_read %sv703[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g726 = vector.transfer_read %sv706[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<8x8xf64>
      %t727 = arith.mulf %dv725, %g726 : vector<8x8xf64>
      %t2728 = arith.addf %t727, %dt724 : vector<8x8xf64>
      %fl729 = arith.addf %t2728, %dt715 : vector<8x8xf64>
      vector.transfer_write %fl729, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc730 = arith.constant dense<0.0> : vector<8x8xf64>
      %a731 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b732 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r733 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a731, %b732, %acc730 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a734 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b735 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r736 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a734, %b735, %r733 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r736, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc737 = arith.constant dense<0.0> : vector<8x8xf64>
      %a738 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b739 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r740 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a738, %b739, %acc737 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a741 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b742 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r743 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a741, %b742, %r740 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r743, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv744 = memref.subview %y3[0, 2, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 16>>
    %sv745 = memref.subview %y3[0, 3, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 24>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa746 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r747 = vector.transfer_read %sv744[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 16>>, vector<8x8xf64>
      %m748 = arith.subf %r747, %fa746 : vector<8x8xf64>
      vector.transfer_write %m748, %sv744[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: 16>>
      %r749 = vector.transfer_read %sv745[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 24>>, vector<8x8xf64>
      %p750 = arith.addf %r749, %fa746 : vector<8x8xf64>
      vector.transfer_write %p750, %sv745[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: 24>>
    }
    gpu.barrier
    %sv751 = memref.subview %v3602[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv752 = memref.subview %v3603[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv753 = memref.subview %g10all[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 192>>
    %sv754 = memref.subview %g11all[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 192>>
    %sv755 = memref.subview %g12all[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 192>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc756 = arith.constant dense<0.0> : vector<8x8xf64>
      %a757 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b758 = vector.transfer_read %sv751[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r759 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a757, %b758, %acc756 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a760 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b761 = vector.transfer_read %sv751[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r762 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a760, %b761, %r759 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g763 = vector.transfer_read %sv754[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<8x8xf64>
      %dt764 = arith.mulf %r762, %g763 : vector<8x8xf64>
      %acc765 = arith.constant dense<0.0> : vector<8x8xf64>
      %a766 = vector.transfer_read %sv751[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b767 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r768 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a766, %b767, %acc765 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a769 = vector.transfer_read %sv751[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b770 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r771 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a769, %b770, %r768 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g772 = vector.transfer_read %sv753[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<8x8xf64>
      %dt773 = arith.mulf %r771, %g772 : vector<8x8xf64>
      %dv774 = vector.transfer_read %sv752[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g775 = vector.transfer_read %sv755[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<8x8xf64>
      %t776 = arith.mulf %dv774, %g775 : vector<8x8xf64>
      %t2777 = arith.addf %t776, %dt773 : vector<8x8xf64>
      %fl778 = arith.addf %t2777, %dt764 : vector<8x8xf64>
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
      %r798 = vector.transfer_read %sv794[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 32>>, vector<8x8xf64>
      %p799 = arith.addf %r798, %fa795 : vector<8x8xf64>
      vector.transfer_write %p799, %sv794[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: 32>>
    }
    gpu.barrier
    %sv800 = memref.subview %v3602[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv801 = memref.subview %v3603[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv802 = memref.subview %g10all[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 256>>
    %sv803 = memref.subview %g11all[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 256>>
    %sv804 = memref.subview %g12all[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 256>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc805 = arith.constant dense<0.0> : vector<8x8xf64>
      %a806 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b807 = vector.transfer_read %sv800[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r808 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a806, %b807, %acc805 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a809 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b810 = vector.transfer_read %sv800[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r811 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a809, %b810, %r808 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g812 = vector.transfer_read %sv803[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<8x8xf64>
      %dt813 = arith.mulf %r811, %g812 : vector<8x8xf64>
      %acc814 = arith.constant dense<0.0> : vector<8x8xf64>
      %a815 = vector.transfer_read %sv800[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b816 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r817 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a815, %b816, %acc814 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a818 = vector.transfer_read %sv800[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b819 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r820 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a818, %b819, %r817 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g821 = vector.transfer_read %sv802[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<8x8xf64>
      %dt822 = arith.mulf %r820, %g821 : vector<8x8xf64>
      %dv823 = vector.transfer_read %sv801[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g824 = vector.transfer_read %sv804[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<8x8xf64>
      %t825 = arith.mulf %dv823, %g824 : vector<8x8xf64>
      %t2826 = arith.addf %t825, %dt822 : vector<8x8xf64>
      %fl827 = arith.addf %t2826, %dt813 : vector<8x8xf64>
      vector.transfer_write %fl827, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc828 = arith.constant dense<0.0> : vector<8x8xf64>
      %a829 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b830 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r831 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a829, %b830, %acc828 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a832 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b833 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r834 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a832, %b833, %r831 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r834, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc835 = arith.constant dense<0.0> : vector<8x8xf64>
      %a836 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b837 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r838 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a836, %b837, %acc835 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a839 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b840 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r841 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a839, %b840, %r838 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r841, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv842 = memref.subview %y3[0, 4, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 32>>
    %sv843 = memref.subview %y3[0, 5, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 40>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa844 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r845 = vector.transfer_read %sv842[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 32>>, vector<8x8xf64>
      %m846 = arith.subf %r845, %fa844 : vector<8x8xf64>
      vector.transfer_write %m846, %sv842[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: 32>>
      %r847 = vector.transfer_read %sv843[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 40>>, vector<8x8xf64>
      %p848 = arith.addf %r847, %fa844 : vector<8x8xf64>
      vector.transfer_write %p848, %sv843[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: 40>>
    }
    gpu.barrier
    %sv849 = memref.subview %v3602[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv850 = memref.subview %v3603[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv851 = memref.subview %g10all[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 320>>
    %sv852 = memref.subview %g11all[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 320>>
    %sv853 = memref.subview %g12all[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 320>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc854 = arith.constant dense<0.0> : vector<8x8xf64>
      %a855 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b856 = vector.transfer_read %sv849[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r857 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a855, %b856, %acc854 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a858 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b859 = vector.transfer_read %sv849[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r860 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a858, %b859, %r857 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g861 = vector.transfer_read %sv852[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<8x8xf64>
      %dt862 = arith.mulf %r860, %g861 : vector<8x8xf64>
      %acc863 = arith.constant dense<0.0> : vector<8x8xf64>
      %a864 = vector.transfer_read %sv849[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b865 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r866 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a864, %b865, %acc863 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a867 = vector.transfer_read %sv849[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b868 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r869 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a867, %b868, %r866 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g870 = vector.transfer_read %sv851[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<8x8xf64>
      %dt871 = arith.mulf %r869, %g870 : vector<8x8xf64>
      %dv872 = vector.transfer_read %sv850[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g873 = vector.transfer_read %sv853[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<8x8xf64>
      %t874 = arith.mulf %dv872, %g873 : vector<8x8xf64>
      %t2875 = arith.addf %t874, %dt871 : vector<8x8xf64>
      %fl876 = arith.addf %t2875, %dt862 : vector<8x8xf64>
      vector.transfer_write %fl876, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc877 = arith.constant dense<0.0> : vector<8x8xf64>
      %a878 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b879 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r880 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a878, %b879, %acc877 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a881 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b882 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r883 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a881, %b882, %r880 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r883, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc884 = arith.constant dense<0.0> : vector<8x8xf64>
      %a885 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b886 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r887 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a885, %b886, %acc884 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a888 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b889 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r890 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a888, %b889, %r887 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r890, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv891 = memref.subview %y3[0, 5, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 40>>
    %sv892 = memref.subview %y3[0, 6, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 48>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa893 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r894 = vector.transfer_read %sv891[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 40>>, vector<8x8xf64>
      %m895 = arith.subf %r894, %fa893 : vector<8x8xf64>
      vector.transfer_write %m895, %sv891[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: 40>>
      %r896 = vector.transfer_read %sv892[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 48>>, vector<8x8xf64>
      %p897 = arith.addf %r896, %fa893 : vector<8x8xf64>
      vector.transfer_write %p897, %sv892[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: 48>>
    }
    gpu.barrier
    %sv898 = memref.subview %v3602[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv899 = memref.subview %v3603[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv900 = memref.subview %g10all[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 384>>
    %sv901 = memref.subview %g11all[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 384>>
    %sv902 = memref.subview %g12all[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 384>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc903 = arith.constant dense<0.0> : vector<8x8xf64>
      %a904 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b905 = vector.transfer_read %sv898[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r906 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a904, %b905, %acc903 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a907 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b908 = vector.transfer_read %sv898[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r909 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a907, %b908, %r906 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g910 = vector.transfer_read %sv901[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<8x8xf64>
      %dt911 = arith.mulf %r909, %g910 : vector<8x8xf64>
      %acc912 = arith.constant dense<0.0> : vector<8x8xf64>
      %a913 = vector.transfer_read %sv898[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b914 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r915 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a913, %b914, %acc912 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a916 = vector.transfer_read %sv898[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b917 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r918 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a916, %b917, %r915 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g919 = vector.transfer_read %sv900[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<8x8xf64>
      %dt920 = arith.mulf %r918, %g919 : vector<8x8xf64>
      %dv921 = vector.transfer_read %sv899[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g922 = vector.transfer_read %sv902[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<8x8xf64>
      %t923 = arith.mulf %dv921, %g922 : vector<8x8xf64>
      %t2924 = arith.addf %t923, %dt920 : vector<8x8xf64>
      %fl925 = arith.addf %t2924, %dt911 : vector<8x8xf64>
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
    %sv940 = memref.subview %y3[0, 6, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 48>>
    %sv941 = memref.subview %y3[0, 7, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: 56>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa942 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r943 = vector.transfer_read %sv940[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 48>>, vector<8x8xf64>
      %m944 = arith.subf %r943, %fa942 : vector<8x8xf64>
      vector.transfer_write %m944, %sv940[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: 48>>
      %r945 = vector.transfer_read %sv941[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 1], offset: 56>>, vector<8x8xf64>
      %p946 = arith.addf %r945, %fa942 : vector<8x8xf64>
      vector.transfer_write %p946, %sv941[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: 56>>
    }
    gpu.barrier
    %sv947 = memref.subview %u[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 0>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc948 = arith.constant dense<0.0> : vector<8x8xf64>
      %a949 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b950 = vector.transfer_read %sv947[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>>, vector<8x4xf64>
      %r951 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a949, %b950, %acc948 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a952 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b953 = vector.transfer_read %sv947[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>>, vector<8x4xf64>
      %r954 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a952, %b953, %r951 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r954, %interp_all[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv955 = memref.subview %u[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 64>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc956 = arith.constant dense<0.0> : vector<8x8xf64>
      %a957 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b958 = vector.transfer_read %sv955[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<8x4xf64>
      %r959 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a957, %b958, %acc956 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a960 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b961 = vector.transfer_read %sv955[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<8x4xf64>
      %r962 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a960, %b961, %r959 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r962, %interp_all[%c0, %c8] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv963 = memref.subview %u[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 128>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc964 = arith.constant dense<0.0> : vector<8x8xf64>
      %a965 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b966 = vector.transfer_read %sv963[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<8x4xf64>
      %r967 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a965, %b966, %acc964 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a968 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b969 = vector.transfer_read %sv963[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<8x4xf64>
      %r970 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a968, %b969, %r967 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r970, %interp_all[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv971 = memref.subview %u[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 192>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc972 = arith.constant dense<0.0> : vector<8x8xf64>
      %a973 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b974 = vector.transfer_read %sv971[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<8x4xf64>
      %r975 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a973, %b974, %acc972 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a976 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b977 = vector.transfer_read %sv971[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<8x4xf64>
      %r978 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a976, %b977, %r975 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r978, %interp_all[%c0, %c24] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv979 = memref.subview %u[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 256>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc980 = arith.constant dense<0.0> : vector<8x8xf64>
      %a981 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b982 = vector.transfer_read %sv979[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<8x4xf64>
      %r983 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a981, %b982, %acc980 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a984 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b985 = vector.transfer_read %sv979[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<8x4xf64>
      %r986 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a984, %b985, %r983 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r986, %interp_all[%c0, %c32] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv987 = memref.subview %u[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 320>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc988 = arith.constant dense<0.0> : vector<8x8xf64>
      %a989 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b990 = vector.transfer_read %sv987[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<8x4xf64>
      %r991 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a989, %b990, %acc988 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a992 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b993 = vector.transfer_read %sv987[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<8x4xf64>
      %r994 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a992, %b993, %r991 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r994, %interp_all[%c0, %c40] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv995 = memref.subview %u[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 384>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc996 = arith.constant dense<0.0> : vector<8x8xf64>
      %a997 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b998 = vector.transfer_read %sv995[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<8x4xf64>
      %r999 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a997, %b998, %acc996 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1000 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1001 = vector.transfer_read %sv995[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<8x4xf64>
      %r1002 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1000, %b1001, %r999 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1002, %interp_all[%c0, %c48] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1003 = memref.subview %u[7, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 448>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1004 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1005 = vector.transfer_read %Btil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1006 = vector.transfer_read %sv1003[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 448>>, vector<8x4xf64>
      %r1007 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1005, %b1006, %acc1004 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1008 = vector.transfer_read %Btil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1009 = vector.transfer_read %sv1003[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 448>>, vector<8x4xf64>
      %r1010 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1008, %b1009, %r1007 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1010, %interp_all[%c0, %c56] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1011 = memref.subview %u[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 0>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1012 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1013 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1014 = vector.transfer_read %sv1011[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>>, vector<8x4xf64>
      %r1015 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1013, %b1014, %acc1012 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1016 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1017 = vector.transfer_read %sv1011[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>>, vector<8x4xf64>
      %r1018 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1016, %b1017, %r1015 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1018, %deriv_all[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1019 = memref.subview %u[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 64>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1020 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1021 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1022 = vector.transfer_read %sv1019[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<8x4xf64>
      %r1023 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1021, %b1022, %acc1020 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1024 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1025 = vector.transfer_read %sv1019[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<8x4xf64>
      %r1026 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1024, %b1025, %r1023 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1026, %deriv_all[%c0, %c8] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1027 = memref.subview %u[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 128>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1028 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1029 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1030 = vector.transfer_read %sv1027[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<8x4xf64>
      %r1031 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1029, %b1030, %acc1028 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1032 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1033 = vector.transfer_read %sv1027[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<8x4xf64>
      %r1034 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1032, %b1033, %r1031 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1034, %deriv_all[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1035 = memref.subview %u[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 192>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1036 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1037 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1038 = vector.transfer_read %sv1035[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<8x4xf64>
      %r1039 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1037, %b1038, %acc1036 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1040 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1041 = vector.transfer_read %sv1035[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<8x4xf64>
      %r1042 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1040, %b1041, %r1039 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1042, %deriv_all[%c0, %c24] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1043 = memref.subview %u[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 256>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1044 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1045 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1046 = vector.transfer_read %sv1043[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<8x4xf64>
      %r1047 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1045, %b1046, %acc1044 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1048 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1049 = vector.transfer_read %sv1043[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<8x4xf64>
      %r1050 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1048, %b1049, %r1047 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1050, %deriv_all[%c0, %c32] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1051 = memref.subview %u[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 320>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1052 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1053 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1054 = vector.transfer_read %sv1051[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<8x4xf64>
      %r1055 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1053, %b1054, %acc1052 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1056 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1057 = vector.transfer_read %sv1051[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<8x4xf64>
      %r1058 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1056, %b1057, %r1055 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1058, %deriv_all[%c0, %c40] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1059 = memref.subview %u[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 384>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1060 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1061 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1062 = vector.transfer_read %sv1059[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<8x4xf64>
      %r1063 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1061, %b1062, %acc1060 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1064 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1065 = vector.transfer_read %sv1059[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<8x4xf64>
      %r1066 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1064, %b1065, %r1063 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1066, %deriv_all[%c0, %c48] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    %sv1067 = memref.subview %u[7, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 448>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1068 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1069 = vector.transfer_read %Dtil[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1070 = vector.transfer_read %sv1067[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 448>>, vector<8x4xf64>
      %r1071 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1069, %b1070, %acc1068 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1072 = vector.transfer_read %Dtil[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1073 = vector.transfer_read %sv1067[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 448>>, vector<8x4xf64>
      %r1074 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1072, %b1073, %r1071 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1074, %deriv_all[%c0, %c56] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %v31075 = memref.reinterpret_cast %interp_all to offset: [0], sizes: [8, 8, 8], strides: [64, 8, 1] : memref<8x64xf64, #gpu.address_space<workgroup>> to memref<8x8x8xf64, #gpu.address_space<workgroup>>
    %v31076 = memref.reinterpret_cast %deriv_all to offset: [0], sizes: [8, 8, 8], strides: [64, 8, 1] : memref<8x64xf64, #gpu.address_space<workgroup>> to memref<8x8x8xf64, #gpu.address_space<workgroup>>
    %sv1077 = memref.subview %v31075[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>
    %sv1078 = memref.subview %v31076[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>
    %sv1079 = memref.subview %g20all[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 0>>
    %sv1080 = memref.subview %g21all[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 0>>
    %sv1081 = memref.subview %g22all[0, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 0>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1082 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1083 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1084 = vector.transfer_read %sv1077[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1085 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1083, %b1084, %acc1082 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1086 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1087 = vector.transfer_read %sv1077[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1088 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1086, %b1087, %r1085 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1089 = vector.transfer_read %sv1080[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>>, vector<8x8xf64>
      %dt1090 = arith.mulf %r1088, %g1089 : vector<8x8xf64>
      %acc1091 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1092 = vector.transfer_read %sv1077[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1093 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1094 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1092, %b1093, %acc1091 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1095 = vector.transfer_read %sv1077[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1096 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1097 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1095, %b1096, %r1094 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1098 = vector.transfer_read %sv1079[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>>, vector<8x8xf64>
      %dt1099 = arith.mulf %r1097, %g1098 : vector<8x8xf64>
      %dv1100 = vector.transfer_read %sv1078[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g1101 = vector.transfer_read %sv1081[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 0>>, vector<8x8xf64>
      %t1102 = arith.mulf %dv1100, %g1101 : vector<8x8xf64>
      %t21103 = arith.addf %t1102, %dt1099 : vector<8x8xf64>
      %fl1104 = arith.addf %t21103, %dt1090 : vector<8x8xf64>
      vector.transfer_write %fl1104, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1105 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1106 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1107 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1108 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1106, %b1107, %acc1105 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1109 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1110 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1111 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1109, %b1110, %r1108 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1111, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1112 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1113 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1114 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1115 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1113, %b1114, %acc1112 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1116 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1117 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1118 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1116, %b1117, %r1115 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1118, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1119 = memref.subview %y3[0, 0, 0] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: 0>>
    %sv1120 = memref.subview %y3[0, 0, 1] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: 1>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1121 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1122 = vector.transfer_read %sv1119[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: 0>>, vector<8x8xf64>
      %m1123 = arith.subf %r1122, %fa1121 : vector<8x8xf64>
      vector.transfer_write %m1123, %sv1119[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: 0>>
      %r1124 = vector.transfer_read %sv1120[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: 1>>, vector<8x8xf64>
      %p1125 = arith.addf %r1124, %fa1121 : vector<8x8xf64>
      vector.transfer_write %p1125, %sv1120[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: 1>>
    }
    gpu.barrier
    %sv1126 = memref.subview %v31075[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv1127 = memref.subview %v31076[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>
    %sv1128 = memref.subview %g20all[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 64>>
    %sv1129 = memref.subview %g21all[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 64>>
    %sv1130 = memref.subview %g22all[1, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 64>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1131 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1132 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1133 = vector.transfer_read %sv1126[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1134 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1132, %b1133, %acc1131 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1135 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1136 = vector.transfer_read %sv1126[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1137 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1135, %b1136, %r1134 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1138 = vector.transfer_read %sv1129[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<8x8xf64>
      %dt1139 = arith.mulf %r1137, %g1138 : vector<8x8xf64>
      %acc1140 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1141 = vector.transfer_read %sv1126[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1142 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1143 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1141, %b1142, %acc1140 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1144 = vector.transfer_read %sv1126[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1145 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1146 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1144, %b1145, %r1143 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1147 = vector.transfer_read %sv1128[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<8x8xf64>
      %dt1148 = arith.mulf %r1146, %g1147 : vector<8x8xf64>
      %dv1149 = vector.transfer_read %sv1127[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g1150 = vector.transfer_read %sv1130[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 64>>, vector<8x8xf64>
      %t1151 = arith.mulf %dv1149, %g1150 : vector<8x8xf64>
      %t21152 = arith.addf %t1151, %dt1148 : vector<8x8xf64>
      %fl1153 = arith.addf %t21152, %dt1139 : vector<8x8xf64>
      vector.transfer_write %fl1153, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1154 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1155 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1156 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1157 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1155, %b1156, %acc1154 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1158 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1159 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1160 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1158, %b1159, %r1157 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1160, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1161 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1162 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1163 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1164 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1162, %b1163, %acc1161 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1165 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1166 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1167 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1165, %b1166, %r1164 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1167, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1168 = memref.subview %y3[0, 0, 1] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: 1>>
    %sv1169 = memref.subview %y3[0, 0, 2] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: 2>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1170 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1171 = vector.transfer_read %sv1168[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: 1>>, vector<8x8xf64>
      %m1172 = arith.subf %r1171, %fa1170 : vector<8x8xf64>
      vector.transfer_write %m1172, %sv1168[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: 1>>
      %r1173 = vector.transfer_read %sv1169[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: 2>>, vector<8x8xf64>
      %p1174 = arith.addf %r1173, %fa1170 : vector<8x8xf64>
      vector.transfer_write %p1174, %sv1169[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: 2>>
    }
    gpu.barrier
    %sv1175 = memref.subview %v31075[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
    %sv1176 = memref.subview %v31076[2, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>
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
      %dt1188 = arith.mulf %r1186, %g1187 : vector<8x8xf64>
      %acc1189 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1190 = vector.transfer_read %sv1175[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1191 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1192 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1190, %b1191, %acc1189 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1193 = vector.transfer_read %sv1175[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1194 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1195 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1193, %b1194, %r1192 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1196 = vector.transfer_read %sv1177[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<8x8xf64>
      %dt1197 = arith.mulf %r1195, %g1196 : vector<8x8xf64>
      %dv1198 = vector.transfer_read %sv1176[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g1199 = vector.transfer_read %sv1179[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 128>>, vector<8x8xf64>
      %t1200 = arith.mulf %dv1198, %g1199 : vector<8x8xf64>
      %t21201 = arith.addf %t1200, %dt1197 : vector<8x8xf64>
      %fl1202 = arith.addf %t21201, %dt1188 : vector<8x8xf64>
      vector.transfer_write %fl1202, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1203 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1204 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1205 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1206 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1204, %b1205, %acc1203 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1207 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1208 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1209 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1207, %b1208, %r1206 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1209, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1210 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1211 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1212 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1213 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1211, %b1212, %acc1210 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1214 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1215 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1216 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1214, %b1215, %r1213 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1216, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1217 = memref.subview %y3[0, 0, 2] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: 2>>
    %sv1218 = memref.subview %y3[0, 0, 3] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: 3>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1219 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1220 = vector.transfer_read %sv1217[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: 2>>, vector<8x8xf64>
      %m1221 = arith.subf %r1220, %fa1219 : vector<8x8xf64>
      vector.transfer_write %m1221, %sv1217[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: 2>>
      %r1222 = vector.transfer_read %sv1218[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: 3>>, vector<8x8xf64>
      %p1223 = arith.addf %r1222, %fa1219 : vector<8x8xf64>
      vector.transfer_write %p1223, %sv1218[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: 3>>
    }
    gpu.barrier
    %sv1224 = memref.subview %v31075[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv1225 = memref.subview %v31076[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>
    %sv1226 = memref.subview %g20all[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 192>>
    %sv1227 = memref.subview %g21all[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 192>>
    %sv1228 = memref.subview %g22all[3, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 192>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1229 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1230 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1231 = vector.transfer_read %sv1224[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1232 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1230, %b1231, %acc1229 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1233 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1234 = vector.transfer_read %sv1224[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1235 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1233, %b1234, %r1232 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1236 = vector.transfer_read %sv1227[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<8x8xf64>
      %dt1237 = arith.mulf %r1235, %g1236 : vector<8x8xf64>
      %acc1238 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1239 = vector.transfer_read %sv1224[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1240 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1241 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1239, %b1240, %acc1238 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1242 = vector.transfer_read %sv1224[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1243 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1244 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1242, %b1243, %r1241 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1245 = vector.transfer_read %sv1226[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<8x8xf64>
      %dt1246 = arith.mulf %r1244, %g1245 : vector<8x8xf64>
      %dv1247 = vector.transfer_read %sv1225[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g1248 = vector.transfer_read %sv1228[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 192>>, vector<8x8xf64>
      %t1249 = arith.mulf %dv1247, %g1248 : vector<8x8xf64>
      %t21250 = arith.addf %t1249, %dt1246 : vector<8x8xf64>
      %fl1251 = arith.addf %t21250, %dt1237 : vector<8x8xf64>
      vector.transfer_write %fl1251, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1252 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1253 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1254 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1255 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1253, %b1254, %acc1252 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1256 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1257 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1258 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1256, %b1257, %r1255 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1258, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1259 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1260 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1261 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1262 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1260, %b1261, %acc1259 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1263 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1264 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1265 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1263, %b1264, %r1262 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1265, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1266 = memref.subview %y3[0, 0, 3] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: 3>>
    %sv1267 = memref.subview %y3[0, 0, 4] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: 4>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1268 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1269 = vector.transfer_read %sv1266[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: 3>>, vector<8x8xf64>
      %m1270 = arith.subf %r1269, %fa1268 : vector<8x8xf64>
      vector.transfer_write %m1270, %sv1266[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: 3>>
      %r1271 = vector.transfer_read %sv1267[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: 4>>, vector<8x8xf64>
      %p1272 = arith.addf %r1271, %fa1268 : vector<8x8xf64>
      vector.transfer_write %p1272, %sv1267[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: 4>>
    }
    gpu.barrier
    %sv1273 = memref.subview %v31075[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv1274 = memref.subview %v31076[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>
    %sv1275 = memref.subview %g20all[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 256>>
    %sv1276 = memref.subview %g21all[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 256>>
    %sv1277 = memref.subview %g22all[4, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 256>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1278 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1279 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1280 = vector.transfer_read %sv1273[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1281 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1279, %b1280, %acc1278 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1282 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1283 = vector.transfer_read %sv1273[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1284 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1282, %b1283, %r1281 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1285 = vector.transfer_read %sv1276[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<8x8xf64>
      %dt1286 = arith.mulf %r1284, %g1285 : vector<8x8xf64>
      %acc1287 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1288 = vector.transfer_read %sv1273[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1289 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1290 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1288, %b1289, %acc1287 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1291 = vector.transfer_read %sv1273[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1292 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1293 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1291, %b1292, %r1290 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1294 = vector.transfer_read %sv1275[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<8x8xf64>
      %dt1295 = arith.mulf %r1293, %g1294 : vector<8x8xf64>
      %dv1296 = vector.transfer_read %sv1274[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g1297 = vector.transfer_read %sv1277[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 256>>, vector<8x8xf64>
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
    %sv1315 = memref.subview %y3[0, 0, 4] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: 4>>
    %sv1316 = memref.subview %y3[0, 0, 5] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: 5>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1317 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1318 = vector.transfer_read %sv1315[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: 4>>, vector<8x8xf64>
      %m1319 = arith.subf %r1318, %fa1317 : vector<8x8xf64>
      vector.transfer_write %m1319, %sv1315[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: 4>>
      %r1320 = vector.transfer_read %sv1316[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: 5>>, vector<8x8xf64>
      %p1321 = arith.addf %r1320, %fa1317 : vector<8x8xf64>
      vector.transfer_write %p1321, %sv1316[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: 5>>
    }
    gpu.barrier
    %sv1322 = memref.subview %v31075[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv1323 = memref.subview %v31076[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>
    %sv1324 = memref.subview %g20all[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 320>>
    %sv1325 = memref.subview %g21all[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 320>>
    %sv1326 = memref.subview %g22all[5, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 320>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1327 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1328 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1329 = vector.transfer_read %sv1322[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1330 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1328, %b1329, %acc1327 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1331 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1332 = vector.transfer_read %sv1322[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1333 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1331, %b1332, %r1330 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1334 = vector.transfer_read %sv1325[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<8x8xf64>
      %dt1335 = arith.mulf %r1333, %g1334 : vector<8x8xf64>
      %acc1336 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1337 = vector.transfer_read %sv1322[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1338 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1339 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1337, %b1338, %acc1336 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1340 = vector.transfer_read %sv1322[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1341 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1342 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1340, %b1341, %r1339 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1343 = vector.transfer_read %sv1324[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<8x8xf64>
      %dt1344 = arith.mulf %r1342, %g1343 : vector<8x8xf64>
      %dv1345 = vector.transfer_read %sv1323[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g1346 = vector.transfer_read %sv1326[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 320>>, vector<8x8xf64>
      %t1347 = arith.mulf %dv1345, %g1346 : vector<8x8xf64>
      %t21348 = arith.addf %t1347, %dt1344 : vector<8x8xf64>
      %fl1349 = arith.addf %t21348, %dt1335 : vector<8x8xf64>
      vector.transfer_write %fl1349, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1350 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1351 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1352 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1353 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1351, %b1352, %acc1350 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1354 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1355 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1356 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1354, %b1355, %r1353 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1356, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1357 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1358 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1359 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1360 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1358, %b1359, %acc1357 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1361 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1362 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1363 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1361, %b1362, %r1360 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1363, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1364 = memref.subview %y3[0, 0, 5] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: 5>>
    %sv1365 = memref.subview %y3[0, 0, 6] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: 6>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1366 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1367 = vector.transfer_read %sv1364[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: 5>>, vector<8x8xf64>
      %m1368 = arith.subf %r1367, %fa1366 : vector<8x8xf64>
      vector.transfer_write %m1368, %sv1364[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: 5>>
      %r1369 = vector.transfer_read %sv1365[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: 6>>, vector<8x8xf64>
      %p1370 = arith.addf %r1369, %fa1366 : vector<8x8xf64>
      vector.transfer_write %p1370, %sv1365[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: 6>>
    }
    gpu.barrier
    %sv1371 = memref.subview %v31075[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv1372 = memref.subview %v31076[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64, #gpu.address_space<workgroup>> to memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>
    %sv1373 = memref.subview %g20all[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 384>>
    %sv1374 = memref.subview %g21all[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 384>>
    %sv1375 = memref.subview %g22all[6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: 384>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1376 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1377 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1378 = vector.transfer_read %sv1371[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1379 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1377, %b1378, %acc1376 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1380 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1381 = vector.transfer_read %sv1371[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1382 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1380, %b1381, %r1379 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g1383 = vector.transfer_read %sv1374[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<8x8xf64>
      %dt1384 = arith.mulf %r1382, %g1383 : vector<8x8xf64>
      %acc1385 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1386 = vector.transfer_read %sv1371[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1387 = vector.transfer_read %Dm[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1388 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1386, %b1387, %acc1385 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1389 = vector.transfer_read %sv1371[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1390 = vector.transfer_read %Dm[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1391 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1389, %b1390, %r1388 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %g1392 = vector.transfer_read %sv1373[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<8x8xf64>
      %dt1393 = arith.mulf %r1391, %g1392 : vector<8x8xf64>
      %dv1394 = vector.transfer_read %sv1372[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %g1395 = vector.transfer_read %sv1375[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8, 1], offset: 384>>, vector<8x8xf64>
      %t1396 = arith.mulf %dv1394, %g1395 : vector<8x8xf64>
      %t21397 = arith.addf %t1396, %dt1393 : vector<8x8xf64>
      %fl1398 = arith.addf %t21397, %dt1384 : vector<8x8xf64>
      vector.transfer_write %fl1398, %flux[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1399 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1400 = vector.transfer_read %flux[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1401 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1402 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1400, %b1401, %acc1399 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %a1403 = vector.transfer_read %flux[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x4xf64>
      %b1404 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %r1405 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1403, %b1404, %r1402 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %r1405, %tmp[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %acc1406 = arith.constant dense<0.0> : vector<8x8xf64>
      %a1407 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1408 = vector.transfer_read %tmp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1409 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1407, %b1408, %acc1406 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %a1410 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1411 = vector.transfer_read %tmp[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<4x8xf64>
      %r1412 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1410, %b1411, %r1409 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %r1412, %intf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, #gpu.address_space<workgroup>>
    }
    gpu.barrier
    %sv1413 = memref.subview %y3[0, 0, 6] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: 6>>
    %sv1414 = memref.subview %y3[0, 0, 7] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: 7>>
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa1415 = vector.transfer_read %intf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, #gpu.address_space<workgroup>>, vector<8x8xf64>
      %r1416 = vector.transfer_read %sv1413[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: 6>>, vector<8x8xf64>
      %m1417 = arith.subf %r1416, %fa1415 : vector<8x8xf64>
      vector.transfer_write %m1417, %sv1413[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: 6>>
      %r1418 = vector.transfer_read %sv1414[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[64, 8], offset: 7>>, vector<8x8xf64>
      %p1419 = arith.addf %r1418, %fa1415 : vector<8x8xf64>
      vector.transfer_write %p1419, %sv1414[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: 7>>
    }
    gpu.barrier
    gpu.return
  }
}
