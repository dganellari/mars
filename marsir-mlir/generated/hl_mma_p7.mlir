#map = affine_map<(d0, d1, d2, d3) -> (d1, d3)>
#map1 = affine_map<(d0, d1, d2, d3) -> (d0, d3, d2)>
#map2 = affine_map<(d0, d1, d2, d3) -> (d0, d1, d2)>
module {
  func.func @laplacian_apply(%arg0: memref<8x8x8xf64>, %arg1: memref<7x8xf64>, %arg2: memref<7x8xf64>, %arg3: memref<8x8xf64>, %arg4: memref<8x8xf64>, %arg5: memref<3x7x8x8x3xf64>) -> memref<8x8x8xf64> {
    %c16_i32 = arith.constant 16 : i32
    %c4_i32 = arith.constant 4 : i32
    %c2_i32 = arith.constant 2 : i32
    %c32_i32 = arith.constant 32 : i32
    %c1_i32 = arith.constant 1 : i32
    %c3_i32 = arith.constant 3 : i32
    %c28_i32 = arith.constant 28 : i32
    %cst = arith.constant dense<0.000000e+00> : vector<1x2xf64>
    %c0 = arith.constant 0 : index
    %c1 = arith.constant 1 : index
    %c7 = arith.constant 7 : index
    %c64 = arith.constant 64 : index
    %c8 = arith.constant 8 : index
    %cst_0 = arith.constant dense<0.000000e+00> : vector<8x8x8xf64>
    %cst_1 = arith.constant dense<0.000000e+00> : vector<7x64xf64>
    %cst_2 = arith.constant dense<0.000000e+00> : vector<8x8xf64>
    %cst_3 = arith.constant dense<0.000000e+00> : vector<8x7x8xf64>
    %cst_4 = arith.constant dense<0.000000e+00> : vector<64x7xf64>
    %cst_5 = arith.constant 0.000000e+00 : f64
    %cst_6 = arith.constant dense<0.000000e+00> : vector<1x1xf64>
    %c2 = arith.constant 2 : index
    %c4 = arith.constant 4 : index
    %thread_id_x = gpu.thread_id  x
    %0 = arith.divui %thread_id_x, %c4 : index
    %1 = arith.remui %thread_id_x, %c4 : index
    %2 = arith.index_cast %thread_id_x : index to i32
    %alloc = memref.alloc() {alignment = 64 : i64} : memref<8x8x8xf64>
    vector.transfer_write %cst_0, %alloc[%c0, %c0, %c0] {in_bounds = [true, true, true]} : vector<8x8x8xf64>, memref<8x8x8xf64>
    %collapse_shape = memref.collapse_shape %arg0 [[0], [1, 2]] : memref<8x8x8xf64> into memref<8x64xf64>
    %alloc_7 = memref.alloc() {alignment = 64 : i64} : memref<7x64xf64>
    vector.transfer_write %cst_1, %alloc_7[%c0, %c0] {in_bounds = [true, true]} : vector<7x64xf64>, memref<7x64xf64>
    %3 = arith.muli %1, %c2 : index
    scf.for %arg6 = %c0 to %c64 step %c8 {
      %subview = memref.subview %alloc_7[0, %arg6] [7, 8] [1, 1] : memref<7x64xf64> to memref<7x8xf64, strided<[64, 1], offset: ?>>
      %65 = vector.transfer_read %subview[%0, %3], %cst_5 {in_bounds = [false, true]} : memref<7x8xf64, strided<[64, 1], offset: ?>>, vector<1x2xf64>
      %66 = scf.for %arg7 = %c0 to %c8 step %c4 iter_args(%arg8 = %65) -> (vector<1x2xf64>) {
        %subview_19 = memref.subview %arg1[0, %arg7] [7, 4] [1, 1] : memref<7x8xf64> to memref<7x4xf64, strided<[8, 1], offset: ?>>
        %subview_20 = memref.subview %collapse_shape[%arg7, %arg6] [4, 8] [1, 1] : memref<8x64xf64> to memref<4x8xf64, strided<[64, 1], offset: ?>>
        %67 = vector.transfer_read %subview_19[%0, %1], %cst_5 {in_bounds = [false, true]} : memref<7x4xf64, strided<[8, 1], offset: ?>>, vector<1x1xf64>
        %68 = vector.transfer_read %subview_20[%1, %0], %cst_5 {in_bounds = [true, true]} : memref<4x8xf64, strided<[64, 1], offset: ?>>, vector<1x1xf64>
        %69 = nvgpu.mma.sync(%67, %68, %arg8) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
        scf.yield %69 : vector<1x2xf64>
      }
      vector.transfer_write %66, %subview[%0, %3] : vector<1x2xf64>, memref<7x8xf64, strided<[64, 1], offset: ?>>
    }
    %expand_shape = memref.expand_shape %alloc_7 [[0], [1, 2]] output_shape [7, 8, 8] : memref<7x64xf64> into memref<7x8x8xf64>
    %alloc_8 = memref.alloc() {alignment = 64 : i64} : memref<7x64xf64>
    vector.transfer_write %cst_1, %alloc_8[%c0, %c0] {in_bounds = [true, true]} : vector<7x64xf64>, memref<7x64xf64>
    %4 = arith.muli %1, %c2 : index
    scf.for %arg6 = %c0 to %c64 step %c8 {
      %subview = memref.subview %alloc_8[0, %arg6] [7, 8] [1, 1] : memref<7x64xf64> to memref<7x8xf64, strided<[64, 1], offset: ?>>
      %65 = vector.transfer_read %subview[%0, %4], %cst_5 {in_bounds = [false, true]} : memref<7x8xf64, strided<[64, 1], offset: ?>>, vector<1x2xf64>
      %66 = scf.for %arg7 = %c0 to %c8 step %c4 iter_args(%arg8 = %65) -> (vector<1x2xf64>) {
        %subview_19 = memref.subview %arg2[0, %arg7] [7, 4] [1, 1] : memref<7x8xf64> to memref<7x4xf64, strided<[8, 1], offset: ?>>
        %subview_20 = memref.subview %collapse_shape[%arg7, %arg6] [4, 8] [1, 1] : memref<8x64xf64> to memref<4x8xf64, strided<[64, 1], offset: ?>>
        %67 = vector.transfer_read %subview_19[%0, %1], %cst_5 {in_bounds = [false, true]} : memref<7x4xf64, strided<[8, 1], offset: ?>>, vector<1x1xf64>
        %68 = vector.transfer_read %subview_20[%1, %0], %cst_5 {in_bounds = [true, true]} : memref<4x8xf64, strided<[64, 1], offset: ?>>, vector<1x1xf64>
        %69 = nvgpu.mma.sync(%67, %68, %arg8) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
        scf.yield %69 : vector<1x2xf64>
      }
      vector.transfer_write %66, %subview[%0, %4] : vector<1x2xf64>, memref<7x8xf64, strided<[64, 1], offset: ?>>
    }
    %expand_shape_9 = memref.expand_shape %alloc_8 [[0], [1, 2]] output_shape [7, 8, 8] : memref<7x64xf64> into memref<7x8x8xf64>
    %5 = arith.addi %1, %c4 : index
    %6 = arith.muli %1, %c2 : index
    %7 = arith.andi %2, %c28_i32 : i32
    %8 = arith.andi %2, %c3_i32 : i32
    %9 = arith.shrui %8, %c1_i32 : i32
    %10 = arith.ori %7, %9 : i32
    %11 = arith.andi %2, %c1_i32 : i32
    %12 = arith.cmpi eq, %11, %c1_i32 : i32
    %13 = arith.ori %10, %c2_i32 : i32
    %14 = arith.muli %8, %c4_i32 : i32
    %15 = arith.shrui %2, %c3_i32 : i32
    %16 = arith.addi %14, %15 : i32
    %17 = arith.shrui %2, %c2_i32 : i32
    %18 = arith.andi %17, %c1_i32 : i32
    %19 = arith.cmpi eq, %18, %c1_i32 : i32
    %20 = arith.addi %16, %c16_i32 : i32
    scf.for %arg6 = %c0 to %c7 step %c1 {
      %subview = memref.subview %expand_shape_9[%arg6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<7x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: ?>>
      %subview_19 = memref.subview %expand_shape[%arg6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<7x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: ?>>
      %alloc_20 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_2, %alloc_20[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %65 = vector.transfer_read %arg4[%0, %1], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %66 = vector.transfer_read %subview_19[%1, %0], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x1xf64>
      %67 = nvgpu.mma.sync(%65, %66, %cst) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %68 = vector.transfer_read %arg4[%0, %5], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %69 = vector.transfer_read %subview_19[%5, %0], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x1xf64>
      %70 = nvgpu.mma.sync(%68, %69, %67) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %70, %alloc_20[%0, %6] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %alloc_21 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_2, %alloc_21[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %71 = vector.transfer_read %subview_19[%0, %1], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x1xf64>
      %72 = vector.transfer_read %arg4[%0, %1], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %73 = nvgpu.mma.sync(%71, %72, %cst) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %74 = vector.transfer_read %subview_19[%0, %5], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x1xf64>
      %75 = vector.transfer_read %arg4[%0, %5], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %76 = nvgpu.mma.sync(%74, %75, %73) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %76, %alloc_21[%0, %6] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %subview_22 = memref.subview %arg5[0, %arg6, 0, 0, 0] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %subview_23 = memref.subview %arg5[0, %arg6, 0, 0, 1] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %subview_24 = memref.subview %arg5[0, %arg6, 0, 0, 2] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %alloc_25 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %77 = vector.transfer_read %subview_24[%0, %6], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<1x2xf64>
      %78 = vector.transfer_read %subview[%0, %6], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
      %79 = arith.mulf %77, %78 : vector<1x2xf64>
      %80 = vector.transfer_read %subview_22[%0, %6], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<1x2xf64>
      %81 = arith.mulf %80, %76 : vector<1x2xf64>
      %82 = arith.addf %79, %81 : vector<1x2xf64>
      %83 = vector.transfer_read %subview_23[%0, %6], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<1x2xf64>
      %84 = arith.mulf %83, %70 : vector<1x2xf64>
      %85 = arith.addf %82, %84 : vector<1x2xf64>
      vector.transfer_write %85, %alloc_25[%0, %6] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %alloc_26 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_2, %alloc_26[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %86 = vector.extract %85[0, 0] : f64 from vector<1x2xf64>
      %87 = vector.extract %85[0, 1] : f64 from vector<1x2xf64>
      %shuffleResult, %valid = gpu.shuffle  idx %86, %10, %c32_i32 : f64
      %shuffleResult_27, %valid_28 = gpu.shuffle  idx %87, %10, %c32_i32 : f64
      %88 = arith.select %12, %shuffleResult_27, %shuffleResult : f64
      %89 = vector.insert %88, %cst_6 [0, 0] : f64 into vector<1x1xf64>
      %90 = vector.transfer_read %arg3[%0, %1], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %91 = nvgpu.mma.sync(%89, %90, %cst) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %shuffleResult_29, %valid_30 = gpu.shuffle  idx %86, %13, %c32_i32 : f64
      %shuffleResult_31, %valid_32 = gpu.shuffle  idx %87, %13, %c32_i32 : f64
      %92 = arith.select %12, %shuffleResult_31, %shuffleResult_29 : f64
      %93 = vector.insert %92, %cst_6 [0, 0] : f64 into vector<1x1xf64>
      %94 = vector.transfer_read %arg3[%0, %5], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %95 = nvgpu.mma.sync(%93, %94, %91) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %95, %alloc_26[%0, %6] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %alloc_33 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_2, %alloc_33[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %96 = vector.transfer_read %arg3[%0, %1], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %97 = vector.extract %95[0, 0] : f64 from vector<1x2xf64>
      %98 = vector.extract %95[0, 1] : f64 from vector<1x2xf64>
      %shuffleResult_34, %valid_35 = gpu.shuffle  idx %97, %16, %c32_i32 : f64
      %shuffleResult_36, %valid_37 = gpu.shuffle  idx %98, %16, %c32_i32 : f64
      %99 = arith.select %19, %shuffleResult_36, %shuffleResult_34 : f64
      %100 = vector.insert %99, %cst_6 [0, 0] : f64 into vector<1x1xf64>
      %101 = nvgpu.mma.sync(%96, %100, %cst) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %102 = vector.transfer_read %arg3[%0, %5], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %shuffleResult_38, %valid_39 = gpu.shuffle  idx %97, %20, %c32_i32 : f64
      %shuffleResult_40, %valid_41 = gpu.shuffle  idx %98, %20, %c32_i32 : f64
      %103 = arith.select %19, %shuffleResult_40, %shuffleResult_38 : f64
      %104 = vector.insert %103, %cst_6 [0, 0] : f64 into vector<1x1xf64>
      %105 = nvgpu.mma.sync(%102, %104, %101) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %105, %alloc_33[%0, %6] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %106 = arith.addi %arg6, %c1 : index
      %subview_42 = memref.subview %alloc[%arg6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: ?>>
      %alloc_43 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %107 = vector.transfer_read %subview_42[%0, %6], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
      %108 = arith.subf %107, %105 : vector<1x2xf64>
      vector.transfer_write %108, %alloc_43[%0, %6] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      vector.transfer_write %108, %subview_42[%0, %6] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      %subview_44 = memref.subview %alloc[%106, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: ?>>
      %alloc_45 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %109 = vector.transfer_read %subview_44[%0, %6], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x2xf64>
      %110 = arith.addf %109, %105 : vector<1x2xf64>
      vector.transfer_write %110, %alloc_45[%0, %6] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      vector.transfer_write %110, %subview_44[%0, %6] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    }
    %alloc_10 = memref.alloc() {alignment = 64 : i64} : memref<8x7x8xf64>
    vector.transfer_write %cst_3, %alloc_10[%c0, %c0, %c0] {in_bounds = [true, true, true]} : vector<8x7x8xf64>, memref<8x7x8xf64>
    %21 = vector.transfer_read %arg0[%c0, %c0, %c0], %cst_5 {in_bounds = [true, true, true]} : memref<8x8x8xf64>, vector<8x8x8xf64>
    %22 = vector.transfer_read %arg1[%c0, %c0], %cst_5 {in_bounds = [true, true]} : memref<7x8xf64>, vector<7x8xf64>
    %23 = vector.transfer_read %alloc_10[%c0, %c0, %c0], %cst_5 {in_bounds = [true, true, true]} : memref<8x7x8xf64>, vector<8x7x8xf64>
    %24 = vector.contract {indexing_maps = [#map, #map1, #map2], iterator_types = ["parallel", "parallel", "parallel", "reduction"], kind = #vector.kind<add>} %22, %21, %23 : vector<7x8xf64>, vector<8x8x8xf64> into vector<8x7x8xf64>
    vector.transfer_write %24, %alloc_10[%c0, %c0, %c0] {in_bounds = [true, true, true]} : vector<8x7x8xf64>, memref<8x7x8xf64>
    %alloc_11 = memref.alloc() {alignment = 64 : i64} : memref<8x7x8xf64>
    vector.transfer_write %cst_3, %alloc_11[%c0, %c0, %c0] {in_bounds = [true, true, true]} : vector<8x7x8xf64>, memref<8x7x8xf64>
    %25 = vector.transfer_read %arg0[%c0, %c0, %c0], %cst_5 {in_bounds = [true, true, true]} : memref<8x8x8xf64>, vector<8x8x8xf64>
    %26 = vector.transfer_read %arg2[%c0, %c0], %cst_5 {in_bounds = [true, true]} : memref<7x8xf64>, vector<7x8xf64>
    %27 = vector.transfer_read %alloc_11[%c0, %c0, %c0], %cst_5 {in_bounds = [true, true, true]} : memref<8x7x8xf64>, vector<8x7x8xf64>
    %28 = vector.contract {indexing_maps = [#map, #map1, #map2], iterator_types = ["parallel", "parallel", "parallel", "reduction"], kind = #vector.kind<add>} %26, %25, %27 : vector<7x8xf64>, vector<8x8x8xf64> into vector<8x7x8xf64>
    vector.transfer_write %28, %alloc_11[%c0, %c0, %c0] {in_bounds = [true, true, true]} : vector<8x7x8xf64>, memref<8x7x8xf64>
    %29 = arith.addi %1, %c4 : index
    %30 = arith.muli %1, %c2 : index
    %31 = arith.andi %2, %c28_i32 : i32
    %32 = arith.andi %2, %c3_i32 : i32
    %33 = arith.shrui %32, %c1_i32 : i32
    %34 = arith.ori %31, %33 : i32
    %35 = arith.andi %2, %c1_i32 : i32
    %36 = arith.cmpi eq, %35, %c1_i32 : i32
    %37 = arith.ori %34, %c2_i32 : i32
    %38 = arith.muli %32, %c4_i32 : i32
    %39 = arith.shrui %2, %c3_i32 : i32
    %40 = arith.addi %38, %39 : i32
    %41 = arith.shrui %2, %c2_i32 : i32
    %42 = arith.andi %41, %c1_i32 : i32
    %43 = arith.cmpi eq, %42, %c1_i32 : i32
    %44 = arith.addi %40, %c16_i32 : i32
    scf.for %arg6 = %c0 to %c7 step %c1 {
      %subview = memref.subview %alloc_11[0, %arg6, 0] [8, 1, 8] [1, 1, 1] : memref<8x7x8xf64> to memref<8x8xf64, strided<[56, 1], offset: ?>>
      %subview_19 = memref.subview %alloc_10[0, %arg6, 0] [8, 1, 8] [1, 1, 1] : memref<8x7x8xf64> to memref<8x8xf64, strided<[56, 1], offset: ?>>
      %alloc_20 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_2, %alloc_20[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %65 = vector.transfer_read %arg4[%0, %1], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %66 = vector.transfer_read %subview_19[%1, %0], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 1], offset: ?>>, vector<1x1xf64>
      %67 = nvgpu.mma.sync(%65, %66, %cst) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %68 = vector.transfer_read %arg4[%0, %29], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %69 = vector.transfer_read %subview_19[%29, %0], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 1], offset: ?>>, vector<1x1xf64>
      %70 = nvgpu.mma.sync(%68, %69, %67) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %70, %alloc_20[%0, %30] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %alloc_21 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_2, %alloc_21[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %71 = vector.transfer_read %subview_19[%0, %1], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 1], offset: ?>>, vector<1x1xf64>
      %72 = vector.transfer_read %arg4[%0, %1], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %73 = nvgpu.mma.sync(%71, %72, %cst) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %74 = vector.transfer_read %subview_19[%0, %29], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 1], offset: ?>>, vector<1x1xf64>
      %75 = vector.transfer_read %arg4[%0, %29], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %76 = nvgpu.mma.sync(%74, %75, %73) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %76, %alloc_21[%0, %30] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %subview_22 = memref.subview %arg5[1, %arg6, 0, 0, 0] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %subview_23 = memref.subview %arg5[1, %arg6, 0, 0, 1] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %subview_24 = memref.subview %arg5[1, %arg6, 0, 0, 2] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %alloc_25 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %77 = vector.transfer_read %subview_24[%0, %30], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<1x2xf64>
      %78 = vector.transfer_read %subview[%0, %30], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 1], offset: ?>>, vector<1x2xf64>
      %79 = arith.mulf %77, %78 : vector<1x2xf64>
      %80 = vector.transfer_read %subview_22[%0, %30], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<1x2xf64>
      %81 = arith.mulf %80, %76 : vector<1x2xf64>
      %82 = arith.addf %79, %81 : vector<1x2xf64>
      %83 = vector.transfer_read %subview_23[%0, %30], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<1x2xf64>
      %84 = arith.mulf %83, %70 : vector<1x2xf64>
      %85 = arith.addf %82, %84 : vector<1x2xf64>
      vector.transfer_write %85, %alloc_25[%0, %30] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %alloc_26 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_2, %alloc_26[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %86 = vector.extract %85[0, 0] : f64 from vector<1x2xf64>
      %87 = vector.extract %85[0, 1] : f64 from vector<1x2xf64>
      %shuffleResult, %valid = gpu.shuffle  idx %86, %34, %c32_i32 : f64
      %shuffleResult_27, %valid_28 = gpu.shuffle  idx %87, %34, %c32_i32 : f64
      %88 = arith.select %36, %shuffleResult_27, %shuffleResult : f64
      %89 = vector.insert %88, %cst_6 [0, 0] : f64 into vector<1x1xf64>
      %90 = vector.transfer_read %arg3[%0, %1], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %91 = nvgpu.mma.sync(%89, %90, %cst) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %shuffleResult_29, %valid_30 = gpu.shuffle  idx %86, %37, %c32_i32 : f64
      %shuffleResult_31, %valid_32 = gpu.shuffle  idx %87, %37, %c32_i32 : f64
      %92 = arith.select %36, %shuffleResult_31, %shuffleResult_29 : f64
      %93 = vector.insert %92, %cst_6 [0, 0] : f64 into vector<1x1xf64>
      %94 = vector.transfer_read %arg3[%0, %29], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %95 = nvgpu.mma.sync(%93, %94, %91) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %95, %alloc_26[%0, %30] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %alloc_33 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_2, %alloc_33[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %96 = vector.transfer_read %arg3[%0, %1], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %97 = vector.extract %95[0, 0] : f64 from vector<1x2xf64>
      %98 = vector.extract %95[0, 1] : f64 from vector<1x2xf64>
      %shuffleResult_34, %valid_35 = gpu.shuffle  idx %97, %40, %c32_i32 : f64
      %shuffleResult_36, %valid_37 = gpu.shuffle  idx %98, %40, %c32_i32 : f64
      %99 = arith.select %43, %shuffleResult_36, %shuffleResult_34 : f64
      %100 = vector.insert %99, %cst_6 [0, 0] : f64 into vector<1x1xf64>
      %101 = nvgpu.mma.sync(%96, %100, %cst) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %102 = vector.transfer_read %arg3[%0, %29], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %shuffleResult_38, %valid_39 = gpu.shuffle  idx %97, %44, %c32_i32 : f64
      %shuffleResult_40, %valid_41 = gpu.shuffle  idx %98, %44, %c32_i32 : f64
      %103 = arith.select %43, %shuffleResult_40, %shuffleResult_38 : f64
      %104 = vector.insert %103, %cst_6 [0, 0] : f64 into vector<1x1xf64>
      %105 = nvgpu.mma.sync(%102, %104, %101) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %105, %alloc_33[%0, %30] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %106 = arith.addi %arg6, %c1 : index
      %subview_42 = memref.subview %alloc[0, %arg6, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: ?>>
      %alloc_43 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %107 = vector.transfer_read %subview_42[%0, %30], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<1x2xf64>
      %108 = arith.subf %107, %105 : vector<1x2xf64>
      vector.transfer_write %108, %alloc_43[%0, %30] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      vector.transfer_write %108, %subview_42[%0, %30] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
      %subview_44 = memref.subview %alloc[0, %106, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: ?>>
      %alloc_45 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %109 = vector.transfer_read %subview_44[%0, %30], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<1x2xf64>
      %110 = arith.addf %109, %105 : vector<1x2xf64>
      vector.transfer_write %110, %alloc_45[%0, %30] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      vector.transfer_write %110, %subview_44[%0, %30] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    }
    %collapse_shape_12 = memref.collapse_shape %arg0 [[0, 1], [2]] : memref<8x8x8xf64> into memref<64x8xf64>
    %alloc_13 = memref.alloc() {alignment = 64 : i64} : memref<8x7xf64>
    %45 = vector.transfer_read %arg1[%c0, %c0], %cst_5 {in_bounds = [true, true]} : memref<7x8xf64>, vector<7x8xf64>
    %46 = vector.transpose %45, [1, 0] : vector<7x8xf64> to vector<8x7xf64>
    vector.transfer_write %46, %alloc_13[%c0, %c0] {in_bounds = [true, true]} : vector<8x7xf64>, memref<8x7xf64>
    %alloc_14 = memref.alloc() {alignment = 64 : i64} : memref<64x7xf64>
    vector.transfer_write %cst_4, %alloc_14[%c0, %c0] {in_bounds = [true, true]} : vector<64x7xf64>, memref<64x7xf64>
    scf.for %arg6 = %c0 to %c8 step %c4 {
      %subview = memref.subview %collapse_shape_12[0, %arg6] [64, 4] [1, 1] : memref<64x8xf64> to memref<64x4xf64, strided<[8, 1], offset: ?>>
      %subview_19 = memref.subview %alloc_13[%arg6, 0] [4, 7] [1, 1] : memref<8x7xf64> to memref<4x7xf64, strided<[7, 1], offset: ?>>
      scf.for %arg7 = %c0 to %c64 step %c1 {
        scf.for %arg8 = %c0 to %c7 step %c1 {
          scf.for %arg9 = %c0 to %c4 step %c1 {
            %65 = memref.load %subview[%arg7, %arg9] : memref<64x4xf64, strided<[8, 1], offset: ?>>
            %66 = memref.load %subview_19[%arg9, %arg8] : memref<4x7xf64, strided<[7, 1], offset: ?>>
            %67 = memref.load %alloc_14[%arg7, %arg8] : memref<64x7xf64>
            %68 = arith.mulf %65, %66 : f64
            %69 = arith.addf %67, %68 : f64
            memref.store %69, %alloc_14[%arg7, %arg8] : memref<64x7xf64>
          }
        }
      }
    }
    %expand_shape_15 = memref.expand_shape %alloc_14 [[0, 1], [2]] output_shape [8, 8, 7] : memref<64x7xf64> into memref<8x8x7xf64>
    %alloc_16 = memref.alloc() {alignment = 64 : i64} : memref<8x7xf64>
    %47 = vector.transfer_read %arg2[%c0, %c0], %cst_5 {in_bounds = [true, true]} : memref<7x8xf64>, vector<7x8xf64>
    %48 = vector.transpose %47, [1, 0] : vector<7x8xf64> to vector<8x7xf64>
    vector.transfer_write %48, %alloc_16[%c0, %c0] {in_bounds = [true, true]} : vector<8x7xf64>, memref<8x7xf64>
    %alloc_17 = memref.alloc() {alignment = 64 : i64} : memref<64x7xf64>
    vector.transfer_write %cst_4, %alloc_17[%c0, %c0] {in_bounds = [true, true]} : vector<64x7xf64>, memref<64x7xf64>
    scf.for %arg6 = %c0 to %c8 step %c4 {
      %subview = memref.subview %collapse_shape_12[0, %arg6] [64, 4] [1, 1] : memref<64x8xf64> to memref<64x4xf64, strided<[8, 1], offset: ?>>
      %subview_19 = memref.subview %alloc_16[%arg6, 0] [4, 7] [1, 1] : memref<8x7xf64> to memref<4x7xf64, strided<[7, 1], offset: ?>>
      scf.for %arg7 = %c0 to %c64 step %c1 {
        scf.for %arg8 = %c0 to %c7 step %c1 {
          scf.for %arg9 = %c0 to %c4 step %c1 {
            %65 = memref.load %subview[%arg7, %arg9] : memref<64x4xf64, strided<[8, 1], offset: ?>>
            %66 = memref.load %subview_19[%arg9, %arg8] : memref<4x7xf64, strided<[7, 1], offset: ?>>
            %67 = memref.load %alloc_17[%arg7, %arg8] : memref<64x7xf64>
            %68 = arith.mulf %65, %66 : f64
            %69 = arith.addf %67, %68 : f64
            memref.store %69, %alloc_17[%arg7, %arg8] : memref<64x7xf64>
          }
        }
      }
    }
    %expand_shape_18 = memref.expand_shape %alloc_17 [[0, 1], [2]] output_shape [8, 8, 7] : memref<64x7xf64> into memref<8x8x7xf64>
    %49 = arith.addi %1, %c4 : index
    %50 = arith.muli %1, %c2 : index
    %51 = arith.andi %2, %c28_i32 : i32
    %52 = arith.andi %2, %c3_i32 : i32
    %53 = arith.shrui %52, %c1_i32 : i32
    %54 = arith.ori %51, %53 : i32
    %55 = arith.andi %2, %c1_i32 : i32
    %56 = arith.cmpi eq, %55, %c1_i32 : i32
    %57 = arith.ori %54, %c2_i32 : i32
    %58 = arith.muli %52, %c4_i32 : i32
    %59 = arith.shrui %2, %c3_i32 : i32
    %60 = arith.addi %58, %59 : i32
    %61 = arith.shrui %2, %c2_i32 : i32
    %62 = arith.andi %61, %c1_i32 : i32
    %63 = arith.cmpi eq, %62, %c1_i32 : i32
    %64 = arith.addi %60, %c16_i32 : i32
    scf.for %arg6 = %c0 to %c7 step %c1 {
      %subview = memref.subview %expand_shape_18[0, 0, %arg6] [8, 8, 1] [1, 1, 1] : memref<8x8x7xf64> to memref<8x8xf64, strided<[56, 7], offset: ?>>
      %subview_19 = memref.subview %expand_shape_15[0, 0, %arg6] [8, 8, 1] [1, 1, 1] : memref<8x8x7xf64> to memref<8x8xf64, strided<[56, 7], offset: ?>>
      %alloc_20 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_2, %alloc_20[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %65 = vector.transfer_read %arg4[%0, %1], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %66 = vector.transfer_read %subview_19[%1, %0], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 7], offset: ?>>, vector<1x1xf64>
      %67 = nvgpu.mma.sync(%65, %66, %cst) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %68 = vector.transfer_read %arg4[%0, %49], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %69 = vector.transfer_read %subview_19[%49, %0], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 7], offset: ?>>, vector<1x1xf64>
      %70 = nvgpu.mma.sync(%68, %69, %67) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %70, %alloc_20[%0, %50] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %alloc_21 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_2, %alloc_21[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %71 = vector.transfer_read %subview_19[%0, %1], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 7], offset: ?>>, vector<1x1xf64>
      %72 = vector.transfer_read %arg4[%0, %1], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %73 = nvgpu.mma.sync(%71, %72, %cst) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %74 = vector.transfer_read %subview_19[%0, %49], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 7], offset: ?>>, vector<1x1xf64>
      %75 = vector.transfer_read %arg4[%0, %49], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %76 = nvgpu.mma.sync(%74, %75, %73) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %76, %alloc_21[%0, %50] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %subview_22 = memref.subview %arg5[2, %arg6, 0, 0, 0] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %subview_23 = memref.subview %arg5[2, %arg6, 0, 0, 1] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %subview_24 = memref.subview %arg5[2, %arg6, 0, 0, 2] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %alloc_25 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %77 = vector.transfer_read %subview_24[%0, %50], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<1x2xf64>
      %78 = vector.transfer_read %subview[%0, %50], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 7], offset: ?>>, vector<1x2xf64>
      %79 = arith.mulf %77, %78 : vector<1x2xf64>
      %80 = vector.transfer_read %subview_22[%0, %50], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<1x2xf64>
      %81 = arith.mulf %80, %76 : vector<1x2xf64>
      %82 = arith.addf %79, %81 : vector<1x2xf64>
      %83 = vector.transfer_read %subview_23[%0, %50], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<1x2xf64>
      %84 = arith.mulf %83, %70 : vector<1x2xf64>
      %85 = arith.addf %82, %84 : vector<1x2xf64>
      vector.transfer_write %85, %alloc_25[%0, %50] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %alloc_26 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_2, %alloc_26[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %86 = vector.extract %85[0, 0] : f64 from vector<1x2xf64>
      %87 = vector.extract %85[0, 1] : f64 from vector<1x2xf64>
      %shuffleResult, %valid = gpu.shuffle  idx %86, %54, %c32_i32 : f64
      %shuffleResult_27, %valid_28 = gpu.shuffle  idx %87, %54, %c32_i32 : f64
      %88 = arith.select %56, %shuffleResult_27, %shuffleResult : f64
      %89 = vector.insert %88, %cst_6 [0, 0] : f64 into vector<1x1xf64>
      %90 = vector.transfer_read %arg3[%0, %1], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %91 = nvgpu.mma.sync(%89, %90, %cst) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %shuffleResult_29, %valid_30 = gpu.shuffle  idx %86, %57, %c32_i32 : f64
      %shuffleResult_31, %valid_32 = gpu.shuffle  idx %87, %57, %c32_i32 : f64
      %92 = arith.select %56, %shuffleResult_31, %shuffleResult_29 : f64
      %93 = vector.insert %92, %cst_6 [0, 0] : f64 into vector<1x1xf64>
      %94 = vector.transfer_read %arg3[%0, %49], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %95 = nvgpu.mma.sync(%93, %94, %91) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %95, %alloc_26[%0, %50] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %alloc_33 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_2, %alloc_33[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %96 = vector.transfer_read %arg3[%0, %1], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %97 = vector.extract %95[0, 0] : f64 from vector<1x2xf64>
      %98 = vector.extract %95[0, 1] : f64 from vector<1x2xf64>
      %shuffleResult_34, %valid_35 = gpu.shuffle  idx %97, %60, %c32_i32 : f64
      %shuffleResult_36, %valid_37 = gpu.shuffle  idx %98, %60, %c32_i32 : f64
      %99 = arith.select %63, %shuffleResult_36, %shuffleResult_34 : f64
      %100 = vector.insert %99, %cst_6 [0, 0] : f64 into vector<1x1xf64>
      %101 = nvgpu.mma.sync(%96, %100, %cst) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %102 = vector.transfer_read %arg3[%0, %49], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %shuffleResult_38, %valid_39 = gpu.shuffle  idx %97, %64, %c32_i32 : f64
      %shuffleResult_40, %valid_41 = gpu.shuffle  idx %98, %64, %c32_i32 : f64
      %103 = arith.select %63, %shuffleResult_40, %shuffleResult_38 : f64
      %104 = vector.insert %103, %cst_6 [0, 0] : f64 into vector<1x1xf64>
      %105 = nvgpu.mma.sync(%102, %104, %101) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %105, %alloc_33[%0, %50] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %106 = arith.addi %arg6, %c1 : index
      %subview_42 = memref.subview %alloc[0, 0, %arg6] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: ?>>
      %alloc_43 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %107 = vector.transfer_read %subview_42[%0, %50], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<1x2xf64>
      %108 = arith.subf %107, %105 : vector<1x2xf64>
      vector.transfer_write %108, %alloc_43[%0, %50] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      vector.transfer_write %108, %subview_42[%0, %50] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
      %subview_44 = memref.subview %alloc[0, 0, %106] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: ?>>
      %alloc_45 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %109 = vector.transfer_read %subview_44[%0, %50], %cst_5 {in_bounds = [true, true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<1x2xf64>
      %110 = arith.addf %109, %105 : vector<1x2xf64>
      vector.transfer_write %110, %alloc_45[%0, %50] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      vector.transfer_write %110, %subview_44[%0, %50] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    }
    return %alloc : memref<8x8x8xf64>
  }
}

