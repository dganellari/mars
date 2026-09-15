#map = affine_map<(d0, d1, d2) -> (d0, d2)>
#map1 = affine_map<(d0, d1, d2) -> (d2, d1)>
#map2 = affine_map<(d0, d1, d2) -> (d0, d1)>
#map3 = affine_map<(d0, d1, d2, d3) -> (d1, d3)>
#map4 = affine_map<(d0, d1, d2, d3) -> (d0, d3, d2)>
#map5 = affine_map<(d0, d1, d2, d3) -> (d0, d1, d2)>
#map6 = affine_map<(d0, d1, d2, d3) -> (d2, d3)>
#map7 = affine_map<(d0, d1, d2, d3) -> (d0, d1, d3)>
module {
  func.func @laplacian_apply(%arg0: memref<8x8x8xf64>, %arg1: memref<7x8xf64>, %arg2: memref<7x8xf64>, %arg3: memref<8x8xf64>, %arg4: memref<8x8xf64>, %arg5: memref<3x7x8x8x3xf64>) -> memref<8x8x8xf64> {
    %cst = arith.constant 0.000000e+00 : f64
    %cst_0 = arith.constant dense<0.000000e+00> : vector<1x1xf64>
    %c2 = arith.constant 2 : index
    %c4 = arith.constant 4 : index
    %thread_id_x = gpu.thread_id  x
    %0 = arith.divui %thread_id_x, %c4 : index
    %1 = arith.remui %thread_id_x, %c4 : index
    %2 = arith.index_cast %thread_id_x : index to i32
    %cst_1 = arith.constant dense<0.000000e+00> : vector<8x8x7xf64>
    %cst_2 = arith.constant dense<0.000000e+00> : vector<8x7x8xf64>
    %cst_3 = arith.constant dense<0.000000e+00> : vector<8x8xf64>
    %cst_4 = arith.constant 0.000000e+00 : f64
    %cst_5 = arith.constant dense<0.000000e+00> : vector<7x64xf64>
    %cst_6 = arith.constant dense<0.000000e+00> : vector<8x8x8xf64>
    %c4_7 = arith.constant 4 : index
    %c8 = arith.constant 8 : index
    %c64 = arith.constant 64 : index
    %c7 = arith.constant 7 : index
    %c1 = arith.constant 1 : index
    %c0 = arith.constant 0 : index
    %alloc = memref.alloc() {alignment = 64 : i64} : memref<8x8x8xf64>
    vector.transfer_write %cst_6, %alloc[%c0, %c0, %c0] {in_bounds = [true, true, true]} : vector<8x8x8xf64>, memref<8x8x8xf64>
    %collapse_shape = memref.collapse_shape %arg0 [[0], [1, 2]] : memref<8x8x8xf64> into memref<8x64xf64>
    %alloc_8 = memref.alloc() {alignment = 64 : i64} : memref<7x64xf64>
    vector.transfer_write %cst_5, %alloc_8[%c0, %c0] {in_bounds = [true, true]} : vector<7x64xf64>, memref<7x64xf64>
    scf.for %arg6 = %c0 to %c64 step %c8 {
      scf.for %arg7 = %c0 to %c8 step %c4_7 {
        %subview = memref.subview %arg1[0, %arg7] [7, 4] [1, 1] : memref<7x8xf64> to memref<7x4xf64, strided<[8, 1], offset: ?>>
        %subview_16 = memref.subview %collapse_shape[%arg7, %arg6] [4, 8] [1, 1] : memref<8x64xf64> to memref<4x8xf64, strided<[64, 1], offset: ?>>
        %subview_17 = memref.subview %alloc_8[0, %arg6] [7, 8] [1, 1] : memref<7x64xf64> to memref<7x8xf64, strided<[64, 1], offset: ?>>
        %22 = vector.transfer_read %subview[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<7x4xf64, strided<[8, 1], offset: ?>>, vector<7x4xf64>
        %23 = vector.transfer_read %subview_16[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<4x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
        %24 = vector.transfer_read %subview_17[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<7x8xf64, strided<[64, 1], offset: ?>>, vector<7x8xf64>
        %25 = vector.contract {indexing_maps = [#map, #map1, #map2], iterator_types = ["parallel", "parallel", "reduction"], kind = #vector.kind<add>} %22, %23, %24 : vector<7x4xf64>, vector<4x8xf64> into vector<7x8xf64>
        vector.transfer_write %25, %subview_17[%c0, %c0] {in_bounds = [true, true]} : vector<7x8xf64>, memref<7x8xf64, strided<[64, 1], offset: ?>>
      }
    }
    %expand_shape = memref.expand_shape %alloc_8 [[0], [1, 2]] output_shape [7, 8, 8] : memref<7x64xf64> into memref<7x8x8xf64>
    %collapse_shape_9 = memref.collapse_shape %arg0 [[0], [1, 2]] : memref<8x8x8xf64> into memref<8x64xf64>
    %alloc_10 = memref.alloc() {alignment = 64 : i64} : memref<7x64xf64>
    vector.transfer_write %cst_5, %alloc_10[%c0, %c0] {in_bounds = [true, true]} : vector<7x64xf64>, memref<7x64xf64>
    scf.for %arg6 = %c0 to %c64 step %c8 {
      scf.for %arg7 = %c0 to %c8 step %c4_7 {
        %subview = memref.subview %arg2[0, %arg7] [7, 4] [1, 1] : memref<7x8xf64> to memref<7x4xf64, strided<[8, 1], offset: ?>>
        %subview_16 = memref.subview %collapse_shape_9[%arg7, %arg6] [4, 8] [1, 1] : memref<8x64xf64> to memref<4x8xf64, strided<[64, 1], offset: ?>>
        %subview_17 = memref.subview %alloc_10[0, %arg6] [7, 8] [1, 1] : memref<7x64xf64> to memref<7x8xf64, strided<[64, 1], offset: ?>>
        %22 = vector.transfer_read %subview[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<7x4xf64, strided<[8, 1], offset: ?>>, vector<7x4xf64>
        %23 = vector.transfer_read %subview_16[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<4x8xf64, strided<[64, 1], offset: ?>>, vector<4x8xf64>
        %24 = vector.transfer_read %subview_17[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<7x8xf64, strided<[64, 1], offset: ?>>, vector<7x8xf64>
        %25 = vector.contract {indexing_maps = [#map, #map1, #map2], iterator_types = ["parallel", "parallel", "reduction"], kind = #vector.kind<add>} %22, %23, %24 : vector<7x4xf64>, vector<4x8xf64> into vector<7x8xf64>
        vector.transfer_write %25, %subview_17[%c0, %c0] {in_bounds = [true, true]} : vector<7x8xf64>, memref<7x8xf64, strided<[64, 1], offset: ?>>
      }
    }
    %expand_shape_11 = memref.expand_shape %alloc_10 [[0], [1, 2]] output_shape [7, 8, 8] : memref<7x64xf64> into memref<7x8x8xf64>
    %3 = scf.for %arg6 = %c0 to %c7 step %c1 iter_args(%arg7 = %alloc) -> (memref<8x8x8xf64>) {
      %subview = memref.subview %expand_shape_11[%arg6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<7x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: ?>>
      %subview_16 = memref.subview %expand_shape[%arg6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<7x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: ?>>
      %alloc_17 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_3, %alloc_17[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %22 = vector.transfer_read %subview_16[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %23 = vector.transfer_read %arg4[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %24 = vector.transfer_read %alloc_17[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %25 = arith.muli %1, %c2 : index
      %26 = arith.addi %c0, %0 : index
      %27 = arith.addi %c0, %25 : index
      %28 = vector.transfer_read %alloc_17[%26, %27], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %29 = vector.transfer_read %arg4[%0, %1], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %30 = vector.transfer_read %subview_16[%1, %0], %cst {in_bounds = [true, true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x1xf64>
      %31 = nvgpu.mma.sync(%29, %30, %28) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %c4_18 = arith.constant 4 : index
      %32 = arith.addi %1, %c4_18 : index
      %33 = vector.transfer_read %arg4[%0, %32], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %c4_19 = arith.constant 4 : index
      %34 = arith.addi %1, %c4_19 : index
      %35 = vector.transfer_read %subview_16[%34, %0], %cst {in_bounds = [true, true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x1xf64>
      %36 = nvgpu.mma.sync(%33, %35, %31) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %37 = arith.muli %1, %c2 : index
      %38 = arith.addi %c0, %0 : index
      %39 = arith.addi %c0, %37 : index
      vector.transfer_write %36, %alloc_17[%38, %39] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %alloc_20 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_3, %alloc_20[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %40 = vector.transfer_read %subview_16[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %41 = vector.transfer_read %arg4[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %42 = vector.transfer_read %alloc_20[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %43 = arith.muli %1, %c2 : index
      %44 = arith.addi %c0, %0 : index
      %45 = arith.addi %c0, %43 : index
      %46 = vector.transfer_read %alloc_20[%44, %45], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %47 = vector.transfer_read %subview_16[%0, %1], %cst {in_bounds = [true, true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x1xf64>
      %48 = vector.transfer_read %arg4[%0, %1], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %49 = nvgpu.mma.sync(%47, %48, %46) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %c4_21 = arith.constant 4 : index
      %50 = arith.addi %1, %c4_21 : index
      %51 = vector.transfer_read %subview_16[%0, %50], %cst {in_bounds = [true, true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x1xf64>
      %c4_22 = arith.constant 4 : index
      %52 = arith.addi %1, %c4_22 : index
      %53 = vector.transfer_read %arg4[%0, %52], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %54 = nvgpu.mma.sync(%51, %53, %49) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %55 = arith.muli %1, %c2 : index
      %56 = arith.addi %c0, %0 : index
      %57 = arith.addi %c0, %55 : index
      vector.transfer_write %54, %alloc_20[%56, %57] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %subview_23 = memref.subview %arg5[0, %arg6, 0, 0, 0] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %subview_24 = memref.subview %arg5[0, %arg6, 0, 0, 1] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %subview_25 = memref.subview %arg5[0, %arg6, 0, 0, 2] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %alloc_26 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %58 = vector.transfer_read %subview[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %59 = vector.transfer_read %alloc_17[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %60 = vector.transfer_read %alloc_20[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %61 = vector.transfer_read %subview_23[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<8x8xf64>
      %62 = vector.transfer_read %subview_24[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<8x8xf64>
      %63 = vector.transfer_read %subview_25[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<8x8xf64>
      %64 = arith.mulf %63, %58 : vector<8x8xf64>
      %65 = arith.mulf %61, %60 : vector<8x8xf64>
      %66 = arith.addf %64, %65 : vector<8x8xf64>
      %67 = arith.mulf %62, %59 : vector<8x8xf64>
      %68 = arith.addf %66, %67 : vector<8x8xf64>
      vector.transfer_write %68, %alloc_26[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %alloc_27 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_3, %alloc_27[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %69 = vector.transfer_read %alloc_26[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %70 = vector.transfer_read %arg3[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %71 = vector.transfer_read %alloc_27[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %72 = arith.muli %1, %c2 : index
      %73 = arith.addi %c0, %0 : index
      %74 = arith.addi %c0, %72 : index
      %75 = vector.transfer_read %alloc_27[%73, %74], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %76 = vector.transfer_read %alloc_26[%0, %1], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %77 = vector.transfer_read %arg3[%0, %1], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %78 = nvgpu.mma.sync(%76, %77, %75) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %c4_28 = arith.constant 4 : index
      %79 = arith.addi %1, %c4_28 : index
      %80 = vector.transfer_read %alloc_26[%0, %79], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %c4_29 = arith.constant 4 : index
      %81 = arith.addi %1, %c4_29 : index
      %82 = vector.transfer_read %arg3[%0, %81], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %83 = nvgpu.mma.sync(%80, %82, %78) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %84 = arith.muli %1, %c2 : index
      %85 = arith.addi %c0, %0 : index
      %86 = arith.addi %c0, %84 : index
      vector.transfer_write %83, %alloc_27[%85, %86] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %alloc_30 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_3, %alloc_30[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %87 = vector.transfer_read %alloc_27[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %88 = vector.transfer_read %arg3[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %89 = vector.transfer_read %alloc_30[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %90 = arith.muli %1, %c2 : index
      %91 = arith.addi %c0, %0 : index
      %92 = arith.addi %c0, %90 : index
      %93 = vector.transfer_read %alloc_30[%91, %92], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %94 = vector.transfer_read %arg3[%0, %1], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %95 = vector.transfer_read %alloc_27[%1, %0], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %96 = nvgpu.mma.sync(%94, %95, %93) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %c4_31 = arith.constant 4 : index
      %97 = arith.addi %1, %c4_31 : index
      %98 = vector.transfer_read %arg3[%0, %97], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %c4_32 = arith.constant 4 : index
      %99 = arith.addi %1, %c4_32 : index
      %100 = vector.transfer_read %alloc_27[%99, %0], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %101 = nvgpu.mma.sync(%98, %100, %96) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %102 = arith.muli %1, %c2 : index
      %103 = arith.addi %c0, %0 : index
      %104 = arith.addi %c0, %102 : index
      vector.transfer_write %101, %alloc_30[%103, %104] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %105 = arith.addi %arg6, %c1 : index
      %subview_33 = memref.subview %arg7[%arg6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: ?>>
      %alloc_34 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %106 = vector.transfer_read %subview_33[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %107 = vector.transfer_read %alloc_30[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %108 = arith.subf %106, %107 : vector<8x8xf64>
      vector.transfer_write %108, %alloc_34[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %subview_35 = memref.subview %arg7[%arg6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: ?>>
      %109 = vector.transfer_read %alloc_34[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      vector.transfer_write %109, %subview_35[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      %subview_36 = memref.subview %arg7[%105, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: ?>>
      %alloc_37 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %110 = vector.transfer_read %subview_36[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %111 = vector.transfer_read %alloc_30[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %112 = arith.addf %110, %111 : vector<8x8xf64>
      vector.transfer_write %112, %alloc_37[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %subview_38 = memref.subview %arg7[%105, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: ?>>
      %113 = vector.transfer_read %alloc_37[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      vector.transfer_write %113, %subview_38[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      scf.yield %arg7 : memref<8x8x8xf64>
    }
    %alloc_12 = memref.alloc() {alignment = 64 : i64} : memref<8x7x8xf64>
    vector.transfer_write %cst_2, %alloc_12[%c0, %c0, %c0] {in_bounds = [true, true, true]} : vector<8x7x8xf64>, memref<8x7x8xf64>
    %4 = vector.transfer_read %arg0[%c0, %c0, %c0], %cst_4 {in_bounds = [true, true, true]} : memref<8x8x8xf64>, vector<8x8x8xf64>
    %5 = vector.transfer_read %arg1[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<7x8xf64>, vector<7x8xf64>
    %6 = vector.transfer_read %alloc_12[%c0, %c0, %c0], %cst_4 {in_bounds = [true, true, true]} : memref<8x7x8xf64>, vector<8x7x8xf64>
    %7 = vector.contract {indexing_maps = [#map3, #map4, #map5], iterator_types = ["parallel", "parallel", "parallel", "reduction"], kind = #vector.kind<add>} %5, %4, %6 : vector<7x8xf64>, vector<8x8x8xf64> into vector<8x7x8xf64>
    vector.transfer_write %7, %alloc_12[%c0, %c0, %c0] {in_bounds = [true, true, true]} : vector<8x7x8xf64>, memref<8x7x8xf64>
    %alloc_13 = memref.alloc() {alignment = 64 : i64} : memref<8x7x8xf64>
    vector.transfer_write %cst_2, %alloc_13[%c0, %c0, %c0] {in_bounds = [true, true, true]} : vector<8x7x8xf64>, memref<8x7x8xf64>
    %8 = vector.transfer_read %arg0[%c0, %c0, %c0], %cst_4 {in_bounds = [true, true, true]} : memref<8x8x8xf64>, vector<8x8x8xf64>
    %9 = vector.transfer_read %arg2[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<7x8xf64>, vector<7x8xf64>
    %10 = vector.transfer_read %alloc_13[%c0, %c0, %c0], %cst_4 {in_bounds = [true, true, true]} : memref<8x7x8xf64>, vector<8x7x8xf64>
    %11 = vector.contract {indexing_maps = [#map3, #map4, #map5], iterator_types = ["parallel", "parallel", "parallel", "reduction"], kind = #vector.kind<add>} %9, %8, %10 : vector<7x8xf64>, vector<8x8x8xf64> into vector<8x7x8xf64>
    vector.transfer_write %11, %alloc_13[%c0, %c0, %c0] {in_bounds = [true, true, true]} : vector<8x7x8xf64>, memref<8x7x8xf64>
    %12 = scf.for %arg6 = %c0 to %c7 step %c1 iter_args(%arg7 = %3) -> (memref<8x8x8xf64>) {
      %subview = memref.subview %alloc_13[0, %arg6, 0] [8, 1, 8] [1, 1, 1] : memref<8x7x8xf64> to memref<8x8xf64, strided<[56, 1], offset: ?>>
      %subview_16 = memref.subview %alloc_12[0, %arg6, 0] [8, 1, 8] [1, 1, 1] : memref<8x7x8xf64> to memref<8x8xf64, strided<[56, 1], offset: ?>>
      %alloc_17 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_3, %alloc_17[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %22 = vector.transfer_read %subview_16[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 1], offset: ?>>, vector<8x8xf64>
      %23 = vector.transfer_read %arg4[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %24 = vector.transfer_read %alloc_17[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %25 = arith.muli %1, %c2 : index
      %26 = arith.addi %c0, %0 : index
      %27 = arith.addi %c0, %25 : index
      %28 = vector.transfer_read %alloc_17[%26, %27], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %29 = vector.transfer_read %arg4[%0, %1], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %30 = vector.transfer_read %subview_16[%1, %0], %cst {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 1], offset: ?>>, vector<1x1xf64>
      %31 = nvgpu.mma.sync(%29, %30, %28) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %c4_18 = arith.constant 4 : index
      %32 = arith.addi %1, %c4_18 : index
      %33 = vector.transfer_read %arg4[%0, %32], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %c4_19 = arith.constant 4 : index
      %34 = arith.addi %1, %c4_19 : index
      %35 = vector.transfer_read %subview_16[%34, %0], %cst {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 1], offset: ?>>, vector<1x1xf64>
      %36 = nvgpu.mma.sync(%33, %35, %31) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %37 = arith.muli %1, %c2 : index
      %38 = arith.addi %c0, %0 : index
      %39 = arith.addi %c0, %37 : index
      vector.transfer_write %36, %alloc_17[%38, %39] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %alloc_20 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_3, %alloc_20[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %40 = vector.transfer_read %subview_16[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 1], offset: ?>>, vector<8x8xf64>
      %41 = vector.transfer_read %arg4[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %42 = vector.transfer_read %alloc_20[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %43 = arith.muli %1, %c2 : index
      %44 = arith.addi %c0, %0 : index
      %45 = arith.addi %c0, %43 : index
      %46 = vector.transfer_read %alloc_20[%44, %45], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %47 = vector.transfer_read %subview_16[%0, %1], %cst {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 1], offset: ?>>, vector<1x1xf64>
      %48 = vector.transfer_read %arg4[%0, %1], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %49 = nvgpu.mma.sync(%47, %48, %46) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %c4_21 = arith.constant 4 : index
      %50 = arith.addi %1, %c4_21 : index
      %51 = vector.transfer_read %subview_16[%0, %50], %cst {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 1], offset: ?>>, vector<1x1xf64>
      %c4_22 = arith.constant 4 : index
      %52 = arith.addi %1, %c4_22 : index
      %53 = vector.transfer_read %arg4[%0, %52], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %54 = nvgpu.mma.sync(%51, %53, %49) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %55 = arith.muli %1, %c2 : index
      %56 = arith.addi %c0, %0 : index
      %57 = arith.addi %c0, %55 : index
      vector.transfer_write %54, %alloc_20[%56, %57] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %subview_23 = memref.subview %arg5[1, %arg6, 0, 0, 0] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %subview_24 = memref.subview %arg5[1, %arg6, 0, 0, 1] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %subview_25 = memref.subview %arg5[1, %arg6, 0, 0, 2] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %alloc_26 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %58 = vector.transfer_read %subview[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 1], offset: ?>>, vector<8x8xf64>
      %59 = vector.transfer_read %alloc_17[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %60 = vector.transfer_read %alloc_20[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %61 = vector.transfer_read %subview_23[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<8x8xf64>
      %62 = vector.transfer_read %subview_24[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<8x8xf64>
      %63 = vector.transfer_read %subview_25[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<8x8xf64>
      %64 = arith.mulf %63, %58 : vector<8x8xf64>
      %65 = arith.mulf %61, %60 : vector<8x8xf64>
      %66 = arith.addf %64, %65 : vector<8x8xf64>
      %67 = arith.mulf %62, %59 : vector<8x8xf64>
      %68 = arith.addf %66, %67 : vector<8x8xf64>
      vector.transfer_write %68, %alloc_26[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %alloc_27 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_3, %alloc_27[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %69 = vector.transfer_read %alloc_26[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %70 = vector.transfer_read %arg3[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %71 = vector.transfer_read %alloc_27[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %72 = arith.muli %1, %c2 : index
      %73 = arith.addi %c0, %0 : index
      %74 = arith.addi %c0, %72 : index
      %75 = vector.transfer_read %alloc_27[%73, %74], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %76 = vector.transfer_read %alloc_26[%0, %1], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %77 = vector.transfer_read %arg3[%0, %1], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %78 = nvgpu.mma.sync(%76, %77, %75) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %c4_28 = arith.constant 4 : index
      %79 = arith.addi %1, %c4_28 : index
      %80 = vector.transfer_read %alloc_26[%0, %79], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %c4_29 = arith.constant 4 : index
      %81 = arith.addi %1, %c4_29 : index
      %82 = vector.transfer_read %arg3[%0, %81], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %83 = nvgpu.mma.sync(%80, %82, %78) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %84 = arith.muli %1, %c2 : index
      %85 = arith.addi %c0, %0 : index
      %86 = arith.addi %c0, %84 : index
      vector.transfer_write %83, %alloc_27[%85, %86] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %alloc_30 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_3, %alloc_30[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %87 = vector.transfer_read %alloc_27[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %88 = vector.transfer_read %arg3[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %89 = vector.transfer_read %alloc_30[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %90 = arith.muli %1, %c2 : index
      %91 = arith.addi %c0, %0 : index
      %92 = arith.addi %c0, %90 : index
      %93 = vector.transfer_read %alloc_30[%91, %92], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %94 = vector.transfer_read %arg3[%0, %1], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %95 = vector.transfer_read %alloc_27[%1, %0], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %96 = nvgpu.mma.sync(%94, %95, %93) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %c4_31 = arith.constant 4 : index
      %97 = arith.addi %1, %c4_31 : index
      %98 = vector.transfer_read %arg3[%0, %97], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %c4_32 = arith.constant 4 : index
      %99 = arith.addi %1, %c4_32 : index
      %100 = vector.transfer_read %alloc_27[%99, %0], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %101 = nvgpu.mma.sync(%98, %100, %96) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %102 = arith.muli %1, %c2 : index
      %103 = arith.addi %c0, %0 : index
      %104 = arith.addi %c0, %102 : index
      vector.transfer_write %101, %alloc_30[%103, %104] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %105 = arith.addi %arg6, %c1 : index
      %subview_33 = memref.subview %arg7[0, %arg6, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: ?>>
      %alloc_34 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %106 = vector.transfer_read %subview_33[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %107 = vector.transfer_read %alloc_30[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %108 = arith.subf %106, %107 : vector<8x8xf64>
      vector.transfer_write %108, %alloc_34[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %subview_35 = memref.subview %arg7[0, %arg6, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: ?>>
      %109 = vector.transfer_read %alloc_34[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      vector.transfer_write %109, %subview_35[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
      %subview_36 = memref.subview %arg7[0, %105, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: ?>>
      %alloc_37 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %110 = vector.transfer_read %subview_36[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %111 = vector.transfer_read %alloc_30[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %112 = arith.addf %110, %111 : vector<8x8xf64>
      vector.transfer_write %112, %alloc_37[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %subview_38 = memref.subview %arg7[0, %105, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: ?>>
      %113 = vector.transfer_read %alloc_37[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      vector.transfer_write %113, %subview_38[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
      scf.yield %arg7 : memref<8x8x8xf64>
    }
    %alloc_14 = memref.alloc() {alignment = 64 : i64} : memref<8x8x7xf64>
    vector.transfer_write %cst_1, %alloc_14[%c0, %c0, %c0] {in_bounds = [true, true, true]} : vector<8x8x7xf64>, memref<8x8x7xf64>
    %13 = vector.transfer_read %arg0[%c0, %c0, %c0], %cst_4 {in_bounds = [true, true, true]} : memref<8x8x8xf64>, vector<8x8x8xf64>
    %14 = vector.transfer_read %arg1[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<7x8xf64>, vector<7x8xf64>
    %15 = vector.transfer_read %alloc_14[%c0, %c0, %c0], %cst_4 {in_bounds = [true, true, true]} : memref<8x8x7xf64>, vector<8x8x7xf64>
    %16 = vector.contract {indexing_maps = [#map6, #map7, #map5], iterator_types = ["parallel", "parallel", "parallel", "reduction"], kind = #vector.kind<add>} %14, %13, %15 : vector<7x8xf64>, vector<8x8x8xf64> into vector<8x8x7xf64>
    vector.transfer_write %16, %alloc_14[%c0, %c0, %c0] {in_bounds = [true, true, true]} : vector<8x8x7xf64>, memref<8x8x7xf64>
    %alloc_15 = memref.alloc() {alignment = 64 : i64} : memref<8x8x7xf64>
    vector.transfer_write %cst_1, %alloc_15[%c0, %c0, %c0] {in_bounds = [true, true, true]} : vector<8x8x7xf64>, memref<8x8x7xf64>
    %17 = vector.transfer_read %arg0[%c0, %c0, %c0], %cst_4 {in_bounds = [true, true, true]} : memref<8x8x8xf64>, vector<8x8x8xf64>
    %18 = vector.transfer_read %arg2[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<7x8xf64>, vector<7x8xf64>
    %19 = vector.transfer_read %alloc_15[%c0, %c0, %c0], %cst_4 {in_bounds = [true, true, true]} : memref<8x8x7xf64>, vector<8x8x7xf64>
    %20 = vector.contract {indexing_maps = [#map6, #map7, #map5], iterator_types = ["parallel", "parallel", "parallel", "reduction"], kind = #vector.kind<add>} %18, %17, %19 : vector<7x8xf64>, vector<8x8x8xf64> into vector<8x8x7xf64>
    vector.transfer_write %20, %alloc_15[%c0, %c0, %c0] {in_bounds = [true, true, true]} : vector<8x8x7xf64>, memref<8x8x7xf64>
    %21 = scf.for %arg6 = %c0 to %c7 step %c1 iter_args(%arg7 = %12) -> (memref<8x8x8xf64>) {
      %subview = memref.subview %alloc_15[0, 0, %arg6] [8, 8, 1] [1, 1, 1] : memref<8x8x7xf64> to memref<8x8xf64, strided<[56, 7], offset: ?>>
      %subview_16 = memref.subview %alloc_14[0, 0, %arg6] [8, 8, 1] [1, 1, 1] : memref<8x8x7xf64> to memref<8x8xf64, strided<[56, 7], offset: ?>>
      %alloc_17 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_3, %alloc_17[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %22 = vector.transfer_read %subview_16[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 7], offset: ?>>, vector<8x8xf64>
      %23 = vector.transfer_read %arg4[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %24 = vector.transfer_read %alloc_17[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %25 = arith.muli %1, %c2 : index
      %26 = arith.addi %c0, %0 : index
      %27 = arith.addi %c0, %25 : index
      %28 = vector.transfer_read %alloc_17[%26, %27], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %29 = vector.transfer_read %arg4[%0, %1], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %30 = vector.transfer_read %subview_16[%1, %0], %cst {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 7], offset: ?>>, vector<1x1xf64>
      %31 = nvgpu.mma.sync(%29, %30, %28) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %c4_18 = arith.constant 4 : index
      %32 = arith.addi %1, %c4_18 : index
      %33 = vector.transfer_read %arg4[%0, %32], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %c4_19 = arith.constant 4 : index
      %34 = arith.addi %1, %c4_19 : index
      %35 = vector.transfer_read %subview_16[%34, %0], %cst {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 7], offset: ?>>, vector<1x1xf64>
      %36 = nvgpu.mma.sync(%33, %35, %31) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %37 = arith.muli %1, %c2 : index
      %38 = arith.addi %c0, %0 : index
      %39 = arith.addi %c0, %37 : index
      vector.transfer_write %36, %alloc_17[%38, %39] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %alloc_20 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_3, %alloc_20[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %40 = vector.transfer_read %subview_16[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 7], offset: ?>>, vector<8x8xf64>
      %41 = vector.transfer_read %arg4[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %42 = vector.transfer_read %alloc_20[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %43 = arith.muli %1, %c2 : index
      %44 = arith.addi %c0, %0 : index
      %45 = arith.addi %c0, %43 : index
      %46 = vector.transfer_read %alloc_20[%44, %45], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %47 = vector.transfer_read %subview_16[%0, %1], %cst {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 7], offset: ?>>, vector<1x1xf64>
      %48 = vector.transfer_read %arg4[%0, %1], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %49 = nvgpu.mma.sync(%47, %48, %46) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %c4_21 = arith.constant 4 : index
      %50 = arith.addi %1, %c4_21 : index
      %51 = vector.transfer_read %subview_16[%0, %50], %cst {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 7], offset: ?>>, vector<1x1xf64>
      %c4_22 = arith.constant 4 : index
      %52 = arith.addi %1, %c4_22 : index
      %53 = vector.transfer_read %arg4[%0, %52], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %54 = nvgpu.mma.sync(%51, %53, %49) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %55 = arith.muli %1, %c2 : index
      %56 = arith.addi %c0, %0 : index
      %57 = arith.addi %c0, %55 : index
      vector.transfer_write %54, %alloc_20[%56, %57] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %subview_23 = memref.subview %arg5[2, %arg6, 0, 0, 0] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %subview_24 = memref.subview %arg5[2, %arg6, 0, 0, 1] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %subview_25 = memref.subview %arg5[2, %arg6, 0, 0, 2] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %alloc_26 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %58 = vector.transfer_read %subview[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 7], offset: ?>>, vector<8x8xf64>
      %59 = vector.transfer_read %alloc_17[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %60 = vector.transfer_read %alloc_20[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %61 = vector.transfer_read %subview_23[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<8x8xf64>
      %62 = vector.transfer_read %subview_24[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<8x8xf64>
      %63 = vector.transfer_read %subview_25[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<8x8xf64>
      %64 = arith.mulf %63, %58 : vector<8x8xf64>
      %65 = arith.mulf %61, %60 : vector<8x8xf64>
      %66 = arith.addf %64, %65 : vector<8x8xf64>
      %67 = arith.mulf %62, %59 : vector<8x8xf64>
      %68 = arith.addf %66, %67 : vector<8x8xf64>
      vector.transfer_write %68, %alloc_26[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %alloc_27 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_3, %alloc_27[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %69 = vector.transfer_read %alloc_26[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %70 = vector.transfer_read %arg3[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %71 = vector.transfer_read %alloc_27[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %72 = arith.muli %1, %c2 : index
      %73 = arith.addi %c0, %0 : index
      %74 = arith.addi %c0, %72 : index
      %75 = vector.transfer_read %alloc_27[%73, %74], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %76 = vector.transfer_read %alloc_26[%0, %1], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %77 = vector.transfer_read %arg3[%0, %1], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %78 = nvgpu.mma.sync(%76, %77, %75) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %c4_28 = arith.constant 4 : index
      %79 = arith.addi %1, %c4_28 : index
      %80 = vector.transfer_read %alloc_26[%0, %79], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %c4_29 = arith.constant 4 : index
      %81 = arith.addi %1, %c4_29 : index
      %82 = vector.transfer_read %arg3[%0, %81], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %83 = nvgpu.mma.sync(%80, %82, %78) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %84 = arith.muli %1, %c2 : index
      %85 = arith.addi %c0, %0 : index
      %86 = arith.addi %c0, %84 : index
      vector.transfer_write %83, %alloc_27[%85, %86] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %alloc_30 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_3, %alloc_30[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %87 = vector.transfer_read %alloc_27[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %88 = vector.transfer_read %arg3[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %89 = vector.transfer_read %alloc_30[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %90 = arith.muli %1, %c2 : index
      %91 = arith.addi %c0, %0 : index
      %92 = arith.addi %c0, %90 : index
      %93 = vector.transfer_read %alloc_30[%91, %92], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %94 = vector.transfer_read %arg3[%0, %1], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %95 = vector.transfer_read %alloc_27[%1, %0], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %96 = nvgpu.mma.sync(%94, %95, %93) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %c4_31 = arith.constant 4 : index
      %97 = arith.addi %1, %c4_31 : index
      %98 = vector.transfer_read %arg3[%0, %97], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %c4_32 = arith.constant 4 : index
      %99 = arith.addi %1, %c4_32 : index
      %100 = vector.transfer_read %alloc_27[%99, %0], %cst {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %101 = nvgpu.mma.sync(%98, %100, %96) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %102 = arith.muli %1, %c2 : index
      %103 = arith.addi %c0, %0 : index
      %104 = arith.addi %c0, %102 : index
      vector.transfer_write %101, %alloc_30[%103, %104] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %105 = arith.addi %arg6, %c1 : index
      %subview_33 = memref.subview %arg7[0, 0, %arg6] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: ?>>
      %alloc_34 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %106 = vector.transfer_read %subview_33[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %107 = vector.transfer_read %alloc_30[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %108 = arith.subf %106, %107 : vector<8x8xf64>
      vector.transfer_write %108, %alloc_34[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %subview_35 = memref.subview %arg7[0, 0, %arg6] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: ?>>
      %109 = vector.transfer_read %alloc_34[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      vector.transfer_write %109, %subview_35[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
      %subview_36 = memref.subview %arg7[0, 0, %105] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: ?>>
      %alloc_37 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %110 = vector.transfer_read %subview_36[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %111 = vector.transfer_read %alloc_30[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %112 = arith.addf %110, %111 : vector<8x8xf64>
      vector.transfer_write %112, %alloc_37[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %subview_38 = memref.subview %arg7[0, 0, %105] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: ?>>
      %113 = vector.transfer_read %alloc_37[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      vector.transfer_write %113, %subview_38[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
      scf.yield %arg7 : memref<8x8x8xf64>
    }
    return %21 : memref<8x8x8xf64>
  }
}

