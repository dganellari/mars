#map = affine_map<(d0, d1, d2, d3) -> (d1, d3)>
#map1 = affine_map<(d0, d1, d2, d3) -> (d0, d3, d2)>
#map2 = affine_map<(d0, d1, d2, d3) -> (d0, d1, d2)>
#map3 = affine_map<(d0, d1, d2, d3) -> (d2, d3)>
#map4 = affine_map<(d0, d1, d2, d3) -> (d0, d1, d3)>
module {
  func.func @laplacian_apply(%arg0: memref<8x8x8xf64>, %arg1: memref<7x8xf64>, %arg2: memref<7x8xf64>, %arg3: memref<8x8xf64>, %arg4: memref<8x8xf64>, %arg5: memref<3x7x8x8x3xf64>) -> memref<8x8x8xf64> {
    %c0 = arith.constant 0 : index
    %c1 = arith.constant 1 : index
    %c7 = arith.constant 7 : index
    %c64 = arith.constant 64 : index
    %c8 = arith.constant 8 : index
    %cst = arith.constant dense<0.000000e+00> : vector<8x8x8xf64>
    %cst_0 = arith.constant dense<0.000000e+00> : vector<7x64xf64>
    %cst_1 = arith.constant dense<0.000000e+00> : vector<8x8xf64>
    %cst_2 = arith.constant dense<0.000000e+00> : vector<8x7x8xf64>
    %cst_3 = arith.constant dense<0.000000e+00> : vector<8x8x7xf64>
    %cst_4 = arith.constant 0.000000e+00 : f64
    %c2 = arith.constant 2 : index
    %c4 = arith.constant 4 : index
    %thread_id_x = gpu.thread_id  x
    %0 = arith.divui %thread_id_x, %c4 : index
    %1 = arith.remui %thread_id_x, %c4 : index
    %alloc = memref.alloc() {alignment = 64 : i64} : memref<8x8x8xf64>
    vector.transfer_write %cst, %alloc[%c0, %c0, %c0] {in_bounds = [true, true, true]} : vector<8x8x8xf64>, memref<8x8x8xf64>
    %collapse_shape = memref.collapse_shape %arg0 [[0], [1, 2]] : memref<8x8x8xf64> into memref<8x64xf64>
    %alloc_5 = memref.alloc() {alignment = 64 : i64} : memref<7x64xf64>
    vector.transfer_write %cst_0, %alloc_5[%c0, %c0] {in_bounds = [true, true]} : vector<7x64xf64>, memref<7x64xf64>
    %2 = arith.muli %1, %c2 : index
    scf.for %arg6 = %c0 to %c64 step %c8 {
      %subview = memref.subview %alloc_5[0, %arg6] [7, 8] [1, 1] : memref<7x64xf64> to memref<7x8xf64, strided<[64, 1], offset: ?>>
      %26 = vector.transfer_read %subview[%0, %2], %cst_4 {in_bounds = [false, true]} : memref<7x8xf64, strided<[64, 1], offset: ?>>, vector<1x2xf64>
      %27 = scf.for %arg7 = %c0 to %c8 step %c4 iter_args(%arg8 = %26) -> (vector<1x2xf64>) {
        %subview_12 = memref.subview %arg1[0, %arg7] [7, 4] [1, 1] : memref<7x8xf64> to memref<7x4xf64, strided<[8, 1], offset: ?>>
        %subview_13 = memref.subview %collapse_shape[%arg7, %arg6] [4, 8] [1, 1] : memref<8x64xf64> to memref<4x8xf64, strided<[64, 1], offset: ?>>
        %28 = vector.transfer_read %subview_12[%0, %1], %cst_4 {in_bounds = [false, true]} : memref<7x4xf64, strided<[8, 1], offset: ?>>, vector<1x1xf64>
        %29 = vector.transfer_read %subview_13[%1, %0], %cst_4 {in_bounds = [true, true]} : memref<4x8xf64, strided<[64, 1], offset: ?>>, vector<1x1xf64>
        %30 = nvgpu.mma.sync(%28, %29, %arg8) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
        scf.yield %30 : vector<1x2xf64>
      }
      vector.transfer_write %27, %subview[%0, %2] : vector<1x2xf64>, memref<7x8xf64, strided<[64, 1], offset: ?>>
    }
    %expand_shape = memref.expand_shape %alloc_5 [[0], [1, 2]] output_shape [7, 8, 8] : memref<7x64xf64> into memref<7x8x8xf64>
    %alloc_6 = memref.alloc() {alignment = 64 : i64} : memref<7x64xf64>
    vector.transfer_write %cst_0, %alloc_6[%c0, %c0] {in_bounds = [true, true]} : vector<7x64xf64>, memref<7x64xf64>
    %3 = arith.muli %1, %c2 : index
    scf.for %arg6 = %c0 to %c64 step %c8 {
      %subview = memref.subview %alloc_6[0, %arg6] [7, 8] [1, 1] : memref<7x64xf64> to memref<7x8xf64, strided<[64, 1], offset: ?>>
      %26 = vector.transfer_read %subview[%0, %3], %cst_4 {in_bounds = [false, true]} : memref<7x8xf64, strided<[64, 1], offset: ?>>, vector<1x2xf64>
      %27 = scf.for %arg7 = %c0 to %c8 step %c4 iter_args(%arg8 = %26) -> (vector<1x2xf64>) {
        %subview_12 = memref.subview %arg2[0, %arg7] [7, 4] [1, 1] : memref<7x8xf64> to memref<7x4xf64, strided<[8, 1], offset: ?>>
        %subview_13 = memref.subview %collapse_shape[%arg7, %arg6] [4, 8] [1, 1] : memref<8x64xf64> to memref<4x8xf64, strided<[64, 1], offset: ?>>
        %28 = vector.transfer_read %subview_12[%0, %1], %cst_4 {in_bounds = [false, true]} : memref<7x4xf64, strided<[8, 1], offset: ?>>, vector<1x1xf64>
        %29 = vector.transfer_read %subview_13[%1, %0], %cst_4 {in_bounds = [true, true]} : memref<4x8xf64, strided<[64, 1], offset: ?>>, vector<1x1xf64>
        %30 = nvgpu.mma.sync(%28, %29, %arg8) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
        scf.yield %30 : vector<1x2xf64>
      }
      vector.transfer_write %27, %subview[%0, %3] : vector<1x2xf64>, memref<7x8xf64, strided<[64, 1], offset: ?>>
    }
    %expand_shape_7 = memref.expand_shape %alloc_6 [[0], [1, 2]] output_shape [7, 8, 8] : memref<7x64xf64> into memref<7x8x8xf64>
    %4 = arith.muli %1, %c2 : index
    %5 = arith.addi %1, %c4 : index
    scf.for %arg6 = %c0 to %c7 step %c1 {
      %subview = memref.subview %expand_shape_7[%arg6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<7x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: ?>>
      %subview_12 = memref.subview %expand_shape[%arg6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<7x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: ?>>
      %alloc_13 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_1, %alloc_13[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %26 = vector.transfer_read %alloc_13[%0, %4], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %27 = vector.transfer_read %arg4[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %28 = vector.transfer_read %subview_12[%1, %0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x1xf64>
      %29 = nvgpu.mma.sync(%27, %28, %26) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %30 = vector.transfer_read %arg4[%0, %5], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %31 = vector.transfer_read %subview_12[%5, %0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x1xf64>
      %32 = nvgpu.mma.sync(%30, %31, %29) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %32, %alloc_13[%0, %4] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %alloc_14 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_1, %alloc_14[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %33 = vector.transfer_read %alloc_14[%0, %4], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %34 = vector.transfer_read %subview_12[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x1xf64>
      %35 = vector.transfer_read %arg4[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %36 = nvgpu.mma.sync(%34, %35, %33) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %37 = vector.transfer_read %subview_12[%0, %5], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x1xf64>
      %38 = vector.transfer_read %arg4[%0, %5], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %39 = nvgpu.mma.sync(%37, %38, %36) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %39, %alloc_14[%0, %4] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %subview_15 = memref.subview %arg5[0, %arg6, 0, 0, 0] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %subview_16 = memref.subview %arg5[0, %arg6, 0, 0, 1] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %subview_17 = memref.subview %arg5[0, %arg6, 0, 0, 2] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %alloc_18 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %40 = vector.transfer_read %subview[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %41 = vector.transfer_read %alloc_13[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %42 = vector.transfer_read %alloc_14[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %43 = vector.transfer_read %subview_15[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<8x8xf64>
      %44 = vector.transfer_read %subview_16[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<8x8xf64>
      %45 = vector.transfer_read %subview_17[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<8x8xf64>
      %46 = arith.mulf %45, %40 : vector<8x8xf64>
      %47 = arith.mulf %43, %42 : vector<8x8xf64>
      %48 = arith.addf %46, %47 : vector<8x8xf64>
      %49 = arith.mulf %44, %41 : vector<8x8xf64>
      %50 = arith.addf %48, %49 : vector<8x8xf64>
      vector.transfer_write %50, %alloc_18[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %alloc_19 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_1, %alloc_19[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %51 = vector.transfer_read %alloc_19[%0, %4], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %52 = vector.transfer_read %alloc_18[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %53 = vector.transfer_read %arg3[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %54 = nvgpu.mma.sync(%52, %53, %51) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %55 = vector.transfer_read %alloc_18[%0, %5], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %56 = vector.transfer_read %arg3[%0, %5], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %57 = nvgpu.mma.sync(%55, %56, %54) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %57, %alloc_19[%0, %4] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %alloc_20 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_1, %alloc_20[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %58 = vector.transfer_read %alloc_20[%0, %4], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %59 = vector.transfer_read %arg3[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %60 = vector.transfer_read %alloc_19[%1, %0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %61 = nvgpu.mma.sync(%59, %60, %58) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %62 = vector.transfer_read %arg3[%0, %5], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %63 = vector.transfer_read %alloc_19[%5, %0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %64 = nvgpu.mma.sync(%62, %63, %61) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %64, %alloc_20[%0, %4] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %65 = arith.addi %arg6, %c1 : index
      %subview_21 = memref.subview %alloc[%arg6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: ?>>
      %alloc_22 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %66 = vector.transfer_read %subview_21[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %67 = vector.transfer_read %alloc_20[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %68 = arith.subf %66, %67 : vector<8x8xf64>
      vector.transfer_write %68, %alloc_22[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %69 = vector.transfer_read %alloc_22[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      vector.transfer_write %69, %subview_21[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      %subview_23 = memref.subview %alloc[%65, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: ?>>
      %alloc_24 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %70 = vector.transfer_read %subview_23[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %71 = vector.transfer_read %alloc_20[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %72 = arith.addf %70, %71 : vector<8x8xf64>
      vector.transfer_write %72, %alloc_24[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %73 = vector.transfer_read %alloc_24[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      vector.transfer_write %73, %subview_23[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
    }
    %alloc_8 = memref.alloc() {alignment = 64 : i64} : memref<8x7x8xf64>
    vector.transfer_write %cst_2, %alloc_8[%c0, %c0, %c0] {in_bounds = [true, true, true]} : vector<8x7x8xf64>, memref<8x7x8xf64>
    %6 = vector.transfer_read %arg0[%c0, %c0, %c0], %cst_4 {in_bounds = [true, true, true]} : memref<8x8x8xf64>, vector<8x8x8xf64>
    %7 = vector.transfer_read %arg1[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<7x8xf64>, vector<7x8xf64>
    %8 = vector.transfer_read %alloc_8[%c0, %c0, %c0], %cst_4 {in_bounds = [true, true, true]} : memref<8x7x8xf64>, vector<8x7x8xf64>
    %9 = vector.contract {indexing_maps = [#map, #map1, #map2], iterator_types = ["parallel", "parallel", "parallel", "reduction"], kind = #vector.kind<add>} %7, %6, %8 : vector<7x8xf64>, vector<8x8x8xf64> into vector<8x7x8xf64>
    vector.transfer_write %9, %alloc_8[%c0, %c0, %c0] {in_bounds = [true, true, true]} : vector<8x7x8xf64>, memref<8x7x8xf64>
    %alloc_9 = memref.alloc() {alignment = 64 : i64} : memref<8x7x8xf64>
    vector.transfer_write %cst_2, %alloc_9[%c0, %c0, %c0] {in_bounds = [true, true, true]} : vector<8x7x8xf64>, memref<8x7x8xf64>
    %10 = vector.transfer_read %arg0[%c0, %c0, %c0], %cst_4 {in_bounds = [true, true, true]} : memref<8x8x8xf64>, vector<8x8x8xf64>
    %11 = vector.transfer_read %arg2[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<7x8xf64>, vector<7x8xf64>
    %12 = vector.transfer_read %alloc_9[%c0, %c0, %c0], %cst_4 {in_bounds = [true, true, true]} : memref<8x7x8xf64>, vector<8x7x8xf64>
    %13 = vector.contract {indexing_maps = [#map, #map1, #map2], iterator_types = ["parallel", "parallel", "parallel", "reduction"], kind = #vector.kind<add>} %11, %10, %12 : vector<7x8xf64>, vector<8x8x8xf64> into vector<8x7x8xf64>
    vector.transfer_write %13, %alloc_9[%c0, %c0, %c0] {in_bounds = [true, true, true]} : vector<8x7x8xf64>, memref<8x7x8xf64>
    %14 = arith.muli %1, %c2 : index
    %15 = arith.addi %1, %c4 : index
    scf.for %arg6 = %c0 to %c7 step %c1 {
      %subview = memref.subview %alloc_9[0, %arg6, 0] [8, 1, 8] [1, 1, 1] : memref<8x7x8xf64> to memref<8x8xf64, strided<[56, 1], offset: ?>>
      %subview_12 = memref.subview %alloc_8[0, %arg6, 0] [8, 1, 8] [1, 1, 1] : memref<8x7x8xf64> to memref<8x8xf64, strided<[56, 1], offset: ?>>
      %alloc_13 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_1, %alloc_13[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %26 = vector.transfer_read %alloc_13[%0, %14], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %27 = vector.transfer_read %arg4[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %28 = vector.transfer_read %subview_12[%1, %0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 1], offset: ?>>, vector<1x1xf64>
      %29 = nvgpu.mma.sync(%27, %28, %26) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %30 = vector.transfer_read %arg4[%0, %15], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %31 = vector.transfer_read %subview_12[%15, %0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 1], offset: ?>>, vector<1x1xf64>
      %32 = nvgpu.mma.sync(%30, %31, %29) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %32, %alloc_13[%0, %14] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %alloc_14 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_1, %alloc_14[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %33 = vector.transfer_read %alloc_14[%0, %14], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %34 = vector.transfer_read %subview_12[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 1], offset: ?>>, vector<1x1xf64>
      %35 = vector.transfer_read %arg4[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %36 = nvgpu.mma.sync(%34, %35, %33) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %37 = vector.transfer_read %subview_12[%0, %15], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 1], offset: ?>>, vector<1x1xf64>
      %38 = vector.transfer_read %arg4[%0, %15], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %39 = nvgpu.mma.sync(%37, %38, %36) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %39, %alloc_14[%0, %14] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %subview_15 = memref.subview %arg5[1, %arg6, 0, 0, 0] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %subview_16 = memref.subview %arg5[1, %arg6, 0, 0, 1] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %subview_17 = memref.subview %arg5[1, %arg6, 0, 0, 2] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %alloc_18 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %40 = vector.transfer_read %subview[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 1], offset: ?>>, vector<8x8xf64>
      %41 = vector.transfer_read %alloc_13[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %42 = vector.transfer_read %alloc_14[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %43 = vector.transfer_read %subview_15[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<8x8xf64>
      %44 = vector.transfer_read %subview_16[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<8x8xf64>
      %45 = vector.transfer_read %subview_17[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<8x8xf64>
      %46 = arith.mulf %45, %40 : vector<8x8xf64>
      %47 = arith.mulf %43, %42 : vector<8x8xf64>
      %48 = arith.addf %46, %47 : vector<8x8xf64>
      %49 = arith.mulf %44, %41 : vector<8x8xf64>
      %50 = arith.addf %48, %49 : vector<8x8xf64>
      vector.transfer_write %50, %alloc_18[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %alloc_19 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_1, %alloc_19[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %51 = vector.transfer_read %alloc_19[%0, %14], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %52 = vector.transfer_read %alloc_18[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %53 = vector.transfer_read %arg3[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %54 = nvgpu.mma.sync(%52, %53, %51) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %55 = vector.transfer_read %alloc_18[%0, %15], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %56 = vector.transfer_read %arg3[%0, %15], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %57 = nvgpu.mma.sync(%55, %56, %54) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %57, %alloc_19[%0, %14] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %alloc_20 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_1, %alloc_20[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %58 = vector.transfer_read %alloc_20[%0, %14], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %59 = vector.transfer_read %arg3[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %60 = vector.transfer_read %alloc_19[%1, %0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %61 = nvgpu.mma.sync(%59, %60, %58) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %62 = vector.transfer_read %arg3[%0, %15], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %63 = vector.transfer_read %alloc_19[%15, %0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %64 = nvgpu.mma.sync(%62, %63, %61) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %64, %alloc_20[%0, %14] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %65 = arith.addi %arg6, %c1 : index
      %subview_21 = memref.subview %alloc[0, %arg6, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: ?>>
      %alloc_22 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %66 = vector.transfer_read %subview_21[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %67 = vector.transfer_read %alloc_20[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %68 = arith.subf %66, %67 : vector<8x8xf64>
      vector.transfer_write %68, %alloc_22[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %69 = vector.transfer_read %alloc_22[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      vector.transfer_write %69, %subview_21[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
      %subview_23 = memref.subview %alloc[0, %65, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: ?>>
      %alloc_24 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %70 = vector.transfer_read %subview_23[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %71 = vector.transfer_read %alloc_20[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %72 = arith.addf %70, %71 : vector<8x8xf64>
      vector.transfer_write %72, %alloc_24[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %73 = vector.transfer_read %alloc_24[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      vector.transfer_write %73, %subview_23[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    }
    %alloc_10 = memref.alloc() {alignment = 64 : i64} : memref<8x8x7xf64>
    vector.transfer_write %cst_3, %alloc_10[%c0, %c0, %c0] {in_bounds = [true, true, true]} : vector<8x8x7xf64>, memref<8x8x7xf64>
    %16 = vector.transfer_read %arg0[%c0, %c0, %c0], %cst_4 {in_bounds = [true, true, true]} : memref<8x8x8xf64>, vector<8x8x8xf64>
    %17 = vector.transfer_read %arg1[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<7x8xf64>, vector<7x8xf64>
    %18 = vector.transfer_read %alloc_10[%c0, %c0, %c0], %cst_4 {in_bounds = [true, true, true]} : memref<8x8x7xf64>, vector<8x8x7xf64>
    %19 = vector.contract {indexing_maps = [#map3, #map4, #map2], iterator_types = ["parallel", "parallel", "parallel", "reduction"], kind = #vector.kind<add>} %17, %16, %18 : vector<7x8xf64>, vector<8x8x8xf64> into vector<8x8x7xf64>
    vector.transfer_write %19, %alloc_10[%c0, %c0, %c0] {in_bounds = [true, true, true]} : vector<8x8x7xf64>, memref<8x8x7xf64>
    %alloc_11 = memref.alloc() {alignment = 64 : i64} : memref<8x8x7xf64>
    vector.transfer_write %cst_3, %alloc_11[%c0, %c0, %c0] {in_bounds = [true, true, true]} : vector<8x8x7xf64>, memref<8x8x7xf64>
    %20 = vector.transfer_read %arg0[%c0, %c0, %c0], %cst_4 {in_bounds = [true, true, true]} : memref<8x8x8xf64>, vector<8x8x8xf64>
    %21 = vector.transfer_read %arg2[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<7x8xf64>, vector<7x8xf64>
    %22 = vector.transfer_read %alloc_11[%c0, %c0, %c0], %cst_4 {in_bounds = [true, true, true]} : memref<8x8x7xf64>, vector<8x8x7xf64>
    %23 = vector.contract {indexing_maps = [#map3, #map4, #map2], iterator_types = ["parallel", "parallel", "parallel", "reduction"], kind = #vector.kind<add>} %21, %20, %22 : vector<7x8xf64>, vector<8x8x8xf64> into vector<8x8x7xf64>
    vector.transfer_write %23, %alloc_11[%c0, %c0, %c0] {in_bounds = [true, true, true]} : vector<8x8x7xf64>, memref<8x8x7xf64>
    %24 = arith.muli %1, %c2 : index
    %25 = arith.addi %1, %c4 : index
    scf.for %arg6 = %c0 to %c7 step %c1 {
      %subview = memref.subview %alloc_11[0, 0, %arg6] [8, 8, 1] [1, 1, 1] : memref<8x8x7xf64> to memref<8x8xf64, strided<[56, 7], offset: ?>>
      %subview_12 = memref.subview %alloc_10[0, 0, %arg6] [8, 8, 1] [1, 1, 1] : memref<8x8x7xf64> to memref<8x8xf64, strided<[56, 7], offset: ?>>
      %alloc_13 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_1, %alloc_13[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %26 = vector.transfer_read %alloc_13[%0, %24], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %27 = vector.transfer_read %arg4[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %28 = vector.transfer_read %subview_12[%1, %0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 7], offset: ?>>, vector<1x1xf64>
      %29 = nvgpu.mma.sync(%27, %28, %26) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %30 = vector.transfer_read %arg4[%0, %25], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %31 = vector.transfer_read %subview_12[%25, %0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 7], offset: ?>>, vector<1x1xf64>
      %32 = nvgpu.mma.sync(%30, %31, %29) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %32, %alloc_13[%0, %24] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %alloc_14 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_1, %alloc_14[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %33 = vector.transfer_read %alloc_14[%0, %24], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %34 = vector.transfer_read %subview_12[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 7], offset: ?>>, vector<1x1xf64>
      %35 = vector.transfer_read %arg4[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %36 = nvgpu.mma.sync(%34, %35, %33) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %37 = vector.transfer_read %subview_12[%0, %25], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 7], offset: ?>>, vector<1x1xf64>
      %38 = vector.transfer_read %arg4[%0, %25], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %39 = nvgpu.mma.sync(%37, %38, %36) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %39, %alloc_14[%0, %24] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %subview_15 = memref.subview %arg5[2, %arg6, 0, 0, 0] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %subview_16 = memref.subview %arg5[2, %arg6, 0, 0, 1] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %subview_17 = memref.subview %arg5[2, %arg6, 0, 0, 2] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %alloc_18 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %40 = vector.transfer_read %subview[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 7], offset: ?>>, vector<8x8xf64>
      %41 = vector.transfer_read %alloc_13[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %42 = vector.transfer_read %alloc_14[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %43 = vector.transfer_read %subview_15[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<8x8xf64>
      %44 = vector.transfer_read %subview_16[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<8x8xf64>
      %45 = vector.transfer_read %subview_17[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<8x8xf64>
      %46 = arith.mulf %45, %40 : vector<8x8xf64>
      %47 = arith.mulf %43, %42 : vector<8x8xf64>
      %48 = arith.addf %46, %47 : vector<8x8xf64>
      %49 = arith.mulf %44, %41 : vector<8x8xf64>
      %50 = arith.addf %48, %49 : vector<8x8xf64>
      vector.transfer_write %50, %alloc_18[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %alloc_19 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_1, %alloc_19[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %51 = vector.transfer_read %alloc_19[%0, %24], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %52 = vector.transfer_read %alloc_18[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %53 = vector.transfer_read %arg3[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %54 = nvgpu.mma.sync(%52, %53, %51) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %55 = vector.transfer_read %alloc_18[%0, %25], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %56 = vector.transfer_read %arg3[%0, %25], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %57 = nvgpu.mma.sync(%55, %56, %54) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %57, %alloc_19[%0, %24] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %alloc_20 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_1, %alloc_20[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %58 = vector.transfer_read %alloc_20[%0, %24], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %59 = vector.transfer_read %arg3[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %60 = vector.transfer_read %alloc_19[%1, %0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %61 = nvgpu.mma.sync(%59, %60, %58) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %62 = vector.transfer_read %arg3[%0, %25], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %63 = vector.transfer_read %alloc_19[%25, %0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %64 = nvgpu.mma.sync(%62, %63, %61) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %64, %alloc_20[%0, %24] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %65 = arith.addi %arg6, %c1 : index
      %subview_21 = memref.subview %alloc[0, 0, %arg6] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: ?>>
      %alloc_22 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %66 = vector.transfer_read %subview_21[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %67 = vector.transfer_read %alloc_20[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %68 = arith.subf %66, %67 : vector<8x8xf64>
      vector.transfer_write %68, %alloc_22[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %69 = vector.transfer_read %alloc_22[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      vector.transfer_write %69, %subview_21[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
      %subview_23 = memref.subview %alloc[0, 0, %65] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: ?>>
      %alloc_24 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %70 = vector.transfer_read %subview_23[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %71 = vector.transfer_read %alloc_20[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %72 = arith.addf %70, %71 : vector<8x8xf64>
      vector.transfer_write %72, %alloc_24[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %73 = vector.transfer_read %alloc_24[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      vector.transfer_write %73, %subview_23[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    }
    return %alloc : memref<8x8x8xf64>
  }
}

