#map = affine_map<(d0, d1, d2, d3) -> (d1, d3)>
#map1 = affine_map<(d0, d1, d2, d3) -> (d0, d3, d2)>
#map2 = affine_map<(d0, d1, d2, d3) -> (d0, d1, d2)>
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
    %cst_3 = arith.constant dense<0.000000e+00> : vector<64x7xf64>
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
      %22 = vector.transfer_read %subview[%0, %2], %cst_4 {in_bounds = [false, true]} : memref<7x8xf64, strided<[64, 1], offset: ?>>, vector<1x2xf64>
      %23 = scf.for %arg7 = %c0 to %c8 step %c4 iter_args(%arg8 = %22) -> (vector<1x2xf64>) {
        %subview_17 = memref.subview %arg1[0, %arg7] [7, 4] [1, 1] : memref<7x8xf64> to memref<7x4xf64, strided<[8, 1], offset: ?>>
        %subview_18 = memref.subview %collapse_shape[%arg7, %arg6] [4, 8] [1, 1] : memref<8x64xf64> to memref<4x8xf64, strided<[64, 1], offset: ?>>
        %24 = vector.transfer_read %subview_17[%0, %1], %cst_4 {in_bounds = [false, true]} : memref<7x4xf64, strided<[8, 1], offset: ?>>, vector<1x1xf64>
        %25 = vector.transfer_read %subview_18[%1, %0], %cst_4 {in_bounds = [true, true]} : memref<4x8xf64, strided<[64, 1], offset: ?>>, vector<1x1xf64>
        %26 = nvgpu.mma.sync(%24, %25, %arg8) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
        scf.yield %26 : vector<1x2xf64>
      }
      vector.transfer_write %23, %subview[%0, %2] : vector<1x2xf64>, memref<7x8xf64, strided<[64, 1], offset: ?>>
    }
    %expand_shape = memref.expand_shape %alloc_5 [[0], [1, 2]] output_shape [7, 8, 8] : memref<7x64xf64> into memref<7x8x8xf64>
    %alloc_6 = memref.alloc() {alignment = 64 : i64} : memref<7x64xf64>
    vector.transfer_write %cst_0, %alloc_6[%c0, %c0] {in_bounds = [true, true]} : vector<7x64xf64>, memref<7x64xf64>
    %3 = arith.muli %1, %c2 : index
    scf.for %arg6 = %c0 to %c64 step %c8 {
      %subview = memref.subview %alloc_6[0, %arg6] [7, 8] [1, 1] : memref<7x64xf64> to memref<7x8xf64, strided<[64, 1], offset: ?>>
      %22 = vector.transfer_read %subview[%0, %3], %cst_4 {in_bounds = [false, true]} : memref<7x8xf64, strided<[64, 1], offset: ?>>, vector<1x2xf64>
      %23 = scf.for %arg7 = %c0 to %c8 step %c4 iter_args(%arg8 = %22) -> (vector<1x2xf64>) {
        %subview_17 = memref.subview %arg2[0, %arg7] [7, 4] [1, 1] : memref<7x8xf64> to memref<7x4xf64, strided<[8, 1], offset: ?>>
        %subview_18 = memref.subview %collapse_shape[%arg7, %arg6] [4, 8] [1, 1] : memref<8x64xf64> to memref<4x8xf64, strided<[64, 1], offset: ?>>
        %24 = vector.transfer_read %subview_17[%0, %1], %cst_4 {in_bounds = [false, true]} : memref<7x4xf64, strided<[8, 1], offset: ?>>, vector<1x1xf64>
        %25 = vector.transfer_read %subview_18[%1, %0], %cst_4 {in_bounds = [true, true]} : memref<4x8xf64, strided<[64, 1], offset: ?>>, vector<1x1xf64>
        %26 = nvgpu.mma.sync(%24, %25, %arg8) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
        scf.yield %26 : vector<1x2xf64>
      }
      vector.transfer_write %23, %subview[%0, %3] : vector<1x2xf64>, memref<7x8xf64, strided<[64, 1], offset: ?>>
    }
    %expand_shape_7 = memref.expand_shape %alloc_6 [[0], [1, 2]] output_shape [7, 8, 8] : memref<7x64xf64> into memref<7x8x8xf64>
    %4 = arith.muli %1, %c2 : index
    %5 = arith.addi %1, %c4 : index
    scf.for %arg6 = %c0 to %c7 step %c1 {
      %subview = memref.subview %expand_shape_7[%arg6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<7x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: ?>>
      %subview_17 = memref.subview %expand_shape[%arg6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<7x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: ?>>
      %alloc_18 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_1, %alloc_18[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %22 = vector.transfer_read %alloc_18[%0, %4], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %23 = vector.transfer_read %arg4[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %24 = vector.transfer_read %subview_17[%1, %0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x1xf64>
      %25 = nvgpu.mma.sync(%23, %24, %22) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %26 = vector.transfer_read %arg4[%0, %5], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %27 = vector.transfer_read %subview_17[%5, %0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x1xf64>
      %28 = nvgpu.mma.sync(%26, %27, %25) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %28, %alloc_18[%0, %4] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %alloc_19 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_1, %alloc_19[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %29 = vector.transfer_read %alloc_19[%0, %4], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %30 = vector.transfer_read %subview_17[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x1xf64>
      %31 = vector.transfer_read %arg4[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %32 = nvgpu.mma.sync(%30, %31, %29) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %33 = vector.transfer_read %subview_17[%0, %5], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<1x1xf64>
      %34 = vector.transfer_read %arg4[%0, %5], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %35 = nvgpu.mma.sync(%33, %34, %32) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %35, %alloc_19[%0, %4] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %subview_20 = memref.subview %arg5[0, %arg6, 0, 0, 0] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %subview_21 = memref.subview %arg5[0, %arg6, 0, 0, 1] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %subview_22 = memref.subview %arg5[0, %arg6, 0, 0, 2] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %alloc_23 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %36 = vector.transfer_read %subview[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %37 = vector.transfer_read %alloc_18[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %38 = vector.transfer_read %alloc_19[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %39 = vector.transfer_read %subview_20[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<8x8xf64>
      %40 = vector.transfer_read %subview_21[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<8x8xf64>
      %41 = vector.transfer_read %subview_22[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<8x8xf64>
      %42 = arith.mulf %41, %36 : vector<8x8xf64>
      %43 = arith.mulf %39, %38 : vector<8x8xf64>
      %44 = arith.addf %42, %43 : vector<8x8xf64>
      %45 = arith.mulf %40, %37 : vector<8x8xf64>
      %46 = arith.addf %44, %45 : vector<8x8xf64>
      vector.transfer_write %46, %alloc_23[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %alloc_24 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_1, %alloc_24[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %47 = vector.transfer_read %alloc_24[%0, %4], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %48 = vector.transfer_read %alloc_23[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %49 = vector.transfer_read %arg3[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %50 = nvgpu.mma.sync(%48, %49, %47) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %51 = vector.transfer_read %alloc_23[%0, %5], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %52 = vector.transfer_read %arg3[%0, %5], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %53 = nvgpu.mma.sync(%51, %52, %50) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %53, %alloc_24[%0, %4] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %alloc_25 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_1, %alloc_25[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %54 = vector.transfer_read %alloc_25[%0, %4], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %55 = vector.transfer_read %arg3[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %56 = vector.transfer_read %alloc_24[%1, %0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %57 = nvgpu.mma.sync(%55, %56, %54) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %58 = vector.transfer_read %arg3[%0, %5], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %59 = vector.transfer_read %alloc_24[%5, %0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %60 = nvgpu.mma.sync(%58, %59, %57) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %60, %alloc_25[%0, %4] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %61 = arith.addi %arg6, %c1 : index
      %subview_26 = memref.subview %alloc[%arg6, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: ?>>
      %alloc_27 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %62 = vector.transfer_read %subview_26[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %63 = vector.transfer_read %alloc_25[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %64 = arith.subf %62, %63 : vector<8x8xf64>
      vector.transfer_write %64, %alloc_27[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %65 = vector.transfer_read %alloc_27[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      vector.transfer_write %65, %subview_26[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      %subview_28 = memref.subview %alloc[%61, 0, 0] [1, 8, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: ?>>
      %alloc_29 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %66 = vector.transfer_read %subview_28[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %67 = vector.transfer_read %alloc_25[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %68 = arith.addf %66, %67 : vector<8x8xf64>
      vector.transfer_write %68, %alloc_29[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %69 = vector.transfer_read %alloc_29[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      vector.transfer_write %69, %subview_28[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
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
      %subview_17 = memref.subview %alloc_8[0, %arg6, 0] [8, 1, 8] [1, 1, 1] : memref<8x7x8xf64> to memref<8x8xf64, strided<[56, 1], offset: ?>>
      %alloc_18 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_1, %alloc_18[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %22 = vector.transfer_read %alloc_18[%0, %14], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %23 = vector.transfer_read %arg4[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %24 = vector.transfer_read %subview_17[%1, %0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 1], offset: ?>>, vector<1x1xf64>
      %25 = nvgpu.mma.sync(%23, %24, %22) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %26 = vector.transfer_read %arg4[%0, %15], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %27 = vector.transfer_read %subview_17[%15, %0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 1], offset: ?>>, vector<1x1xf64>
      %28 = nvgpu.mma.sync(%26, %27, %25) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %28, %alloc_18[%0, %14] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %alloc_19 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_1, %alloc_19[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %29 = vector.transfer_read %alloc_19[%0, %14], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %30 = vector.transfer_read %subview_17[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 1], offset: ?>>, vector<1x1xf64>
      %31 = vector.transfer_read %arg4[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %32 = nvgpu.mma.sync(%30, %31, %29) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %33 = vector.transfer_read %subview_17[%0, %15], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 1], offset: ?>>, vector<1x1xf64>
      %34 = vector.transfer_read %arg4[%0, %15], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %35 = nvgpu.mma.sync(%33, %34, %32) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %35, %alloc_19[%0, %14] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %subview_20 = memref.subview %arg5[1, %arg6, 0, 0, 0] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %subview_21 = memref.subview %arg5[1, %arg6, 0, 0, 1] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %subview_22 = memref.subview %arg5[1, %arg6, 0, 0, 2] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %alloc_23 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %36 = vector.transfer_read %subview[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 1], offset: ?>>, vector<8x8xf64>
      %37 = vector.transfer_read %alloc_18[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %38 = vector.transfer_read %alloc_19[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %39 = vector.transfer_read %subview_20[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<8x8xf64>
      %40 = vector.transfer_read %subview_21[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<8x8xf64>
      %41 = vector.transfer_read %subview_22[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<8x8xf64>
      %42 = arith.mulf %41, %36 : vector<8x8xf64>
      %43 = arith.mulf %39, %38 : vector<8x8xf64>
      %44 = arith.addf %42, %43 : vector<8x8xf64>
      %45 = arith.mulf %40, %37 : vector<8x8xf64>
      %46 = arith.addf %44, %45 : vector<8x8xf64>
      vector.transfer_write %46, %alloc_23[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %alloc_24 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_1, %alloc_24[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %47 = vector.transfer_read %alloc_24[%0, %14], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %48 = vector.transfer_read %alloc_23[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %49 = vector.transfer_read %arg3[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %50 = nvgpu.mma.sync(%48, %49, %47) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %51 = vector.transfer_read %alloc_23[%0, %15], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %52 = vector.transfer_read %arg3[%0, %15], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %53 = nvgpu.mma.sync(%51, %52, %50) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %53, %alloc_24[%0, %14] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %alloc_25 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_1, %alloc_25[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %54 = vector.transfer_read %alloc_25[%0, %14], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %55 = vector.transfer_read %arg3[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %56 = vector.transfer_read %alloc_24[%1, %0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %57 = nvgpu.mma.sync(%55, %56, %54) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %58 = vector.transfer_read %arg3[%0, %15], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %59 = vector.transfer_read %alloc_24[%15, %0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %60 = nvgpu.mma.sync(%58, %59, %57) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %60, %alloc_25[%0, %14] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %61 = arith.addi %arg6, %c1 : index
      %subview_26 = memref.subview %alloc[0, %arg6, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: ?>>
      %alloc_27 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %62 = vector.transfer_read %subview_26[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %63 = vector.transfer_read %alloc_25[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %64 = arith.subf %62, %63 : vector<8x8xf64>
      vector.transfer_write %64, %alloc_27[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %65 = vector.transfer_read %alloc_27[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      vector.transfer_write %65, %subview_26[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
      %subview_28 = memref.subview %alloc[0, %61, 0] [8, 1, 8] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 1], offset: ?>>
      %alloc_29 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %66 = vector.transfer_read %subview_28[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[64, 1], offset: ?>>, vector<8x8xf64>
      %67 = vector.transfer_read %alloc_25[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %68 = arith.addf %66, %67 : vector<8x8xf64>
      vector.transfer_write %68, %alloc_29[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %69 = vector.transfer_read %alloc_29[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      vector.transfer_write %69, %subview_28[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 1], offset: ?>>
    }
    %collapse_shape_10 = memref.collapse_shape %arg0 [[0, 1], [2]] : memref<8x8x8xf64> into memref<64x8xf64>
    %alloc_11 = memref.alloc() {alignment = 64 : i64} : memref<8x7xf64>
    %16 = vector.transfer_read %arg1[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<7x8xf64>, vector<7x8xf64>
    %17 = vector.transpose %16, [1, 0] : vector<7x8xf64> to vector<8x7xf64>
    vector.transfer_write %17, %alloc_11[%c0, %c0] {in_bounds = [true, true]} : vector<8x7xf64>, memref<8x7xf64>
    %alloc_12 = memref.alloc() {alignment = 64 : i64} : memref<64x7xf64>
    vector.transfer_write %cst_3, %alloc_12[%c0, %c0] {in_bounds = [true, true]} : vector<64x7xf64>, memref<64x7xf64>
    scf.for %arg6 = %c0 to %c8 step %c4 {
      %subview = memref.subview %collapse_shape_10[0, %arg6] [64, 4] [1, 1] : memref<64x8xf64> to memref<64x4xf64, strided<[8, 1], offset: ?>>
      %subview_17 = memref.subview %alloc_11[%arg6, 0] [4, 7] [1, 1] : memref<8x7xf64> to memref<4x7xf64, strided<[7, 1], offset: ?>>
      scf.for %arg7 = %c0 to %c64 step %c1 {
        scf.for %arg8 = %c0 to %c7 step %c1 {
          scf.for %arg9 = %c0 to %c4 step %c1 {
            %22 = memref.load %subview[%arg7, %arg9] : memref<64x4xf64, strided<[8, 1], offset: ?>>
            %23 = memref.load %subview_17[%arg9, %arg8] : memref<4x7xf64, strided<[7, 1], offset: ?>>
            %24 = memref.load %alloc_12[%arg7, %arg8] : memref<64x7xf64>
            %25 = arith.mulf %22, %23 : f64
            %26 = arith.addf %24, %25 : f64
            memref.store %26, %alloc_12[%arg7, %arg8] : memref<64x7xf64>
          }
        }
      }
    }
    %expand_shape_13 = memref.expand_shape %alloc_12 [[0, 1], [2]] output_shape [8, 8, 7] : memref<64x7xf64> into memref<8x8x7xf64>
    %alloc_14 = memref.alloc() {alignment = 64 : i64} : memref<8x7xf64>
    %18 = vector.transfer_read %arg2[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<7x8xf64>, vector<7x8xf64>
    %19 = vector.transpose %18, [1, 0] : vector<7x8xf64> to vector<8x7xf64>
    vector.transfer_write %19, %alloc_14[%c0, %c0] {in_bounds = [true, true]} : vector<8x7xf64>, memref<8x7xf64>
    %alloc_15 = memref.alloc() {alignment = 64 : i64} : memref<64x7xf64>
    vector.transfer_write %cst_3, %alloc_15[%c0, %c0] {in_bounds = [true, true]} : vector<64x7xf64>, memref<64x7xf64>
    scf.for %arg6 = %c0 to %c8 step %c4 {
      %subview = memref.subview %collapse_shape_10[0, %arg6] [64, 4] [1, 1] : memref<64x8xf64> to memref<64x4xf64, strided<[8, 1], offset: ?>>
      %subview_17 = memref.subview %alloc_14[%arg6, 0] [4, 7] [1, 1] : memref<8x7xf64> to memref<4x7xf64, strided<[7, 1], offset: ?>>
      scf.for %arg7 = %c0 to %c64 step %c1 {
        scf.for %arg8 = %c0 to %c7 step %c1 {
          scf.for %arg9 = %c0 to %c4 step %c1 {
            %22 = memref.load %subview[%arg7, %arg9] : memref<64x4xf64, strided<[8, 1], offset: ?>>
            %23 = memref.load %subview_17[%arg9, %arg8] : memref<4x7xf64, strided<[7, 1], offset: ?>>
            %24 = memref.load %alloc_15[%arg7, %arg8] : memref<64x7xf64>
            %25 = arith.mulf %22, %23 : f64
            %26 = arith.addf %24, %25 : f64
            memref.store %26, %alloc_15[%arg7, %arg8] : memref<64x7xf64>
          }
        }
      }
    }
    %expand_shape_16 = memref.expand_shape %alloc_15 [[0, 1], [2]] output_shape [8, 8, 7] : memref<64x7xf64> into memref<8x8x7xf64>
    %20 = arith.muli %1, %c2 : index
    %21 = arith.addi %1, %c4 : index
    scf.for %arg6 = %c0 to %c7 step %c1 {
      %subview = memref.subview %expand_shape_16[0, 0, %arg6] [8, 8, 1] [1, 1, 1] : memref<8x8x7xf64> to memref<8x8xf64, strided<[56, 7], offset: ?>>
      %subview_17 = memref.subview %expand_shape_13[0, 0, %arg6] [8, 8, 1] [1, 1, 1] : memref<8x8x7xf64> to memref<8x8xf64, strided<[56, 7], offset: ?>>
      %alloc_18 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_1, %alloc_18[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %22 = vector.transfer_read %alloc_18[%0, %20], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %23 = vector.transfer_read %arg4[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %24 = vector.transfer_read %subview_17[%1, %0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 7], offset: ?>>, vector<1x1xf64>
      %25 = nvgpu.mma.sync(%23, %24, %22) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %26 = vector.transfer_read %arg4[%0, %21], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %27 = vector.transfer_read %subview_17[%21, %0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 7], offset: ?>>, vector<1x1xf64>
      %28 = nvgpu.mma.sync(%26, %27, %25) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %28, %alloc_18[%0, %20] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %alloc_19 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_1, %alloc_19[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %29 = vector.transfer_read %alloc_19[%0, %20], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %30 = vector.transfer_read %subview_17[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 7], offset: ?>>, vector<1x1xf64>
      %31 = vector.transfer_read %arg4[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %32 = nvgpu.mma.sync(%30, %31, %29) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %33 = vector.transfer_read %subview_17[%0, %21], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 7], offset: ?>>, vector<1x1xf64>
      %34 = vector.transfer_read %arg4[%0, %21], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %35 = nvgpu.mma.sync(%33, %34, %32) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %35, %alloc_19[%0, %20] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %subview_20 = memref.subview %arg5[2, %arg6, 0, 0, 0] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %subview_21 = memref.subview %arg5[2, %arg6, 0, 0, 1] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %subview_22 = memref.subview %arg5[2, %arg6, 0, 0, 2] [1, 1, 8, 8, 1] [1, 1, 1, 1, 1] : memref<3x7x8x8x3xf64> to memref<8x8xf64, strided<[24, 3], offset: ?>>
      %alloc_23 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %36 = vector.transfer_read %subview[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[56, 7], offset: ?>>, vector<8x8xf64>
      %37 = vector.transfer_read %alloc_18[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %38 = vector.transfer_read %alloc_19[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %39 = vector.transfer_read %subview_20[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<8x8xf64>
      %40 = vector.transfer_read %subview_21[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<8x8xf64>
      %41 = vector.transfer_read %subview_22[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[24, 3], offset: ?>>, vector<8x8xf64>
      %42 = arith.mulf %41, %36 : vector<8x8xf64>
      %43 = arith.mulf %39, %38 : vector<8x8xf64>
      %44 = arith.addf %42, %43 : vector<8x8xf64>
      %45 = arith.mulf %40, %37 : vector<8x8xf64>
      %46 = arith.addf %44, %45 : vector<8x8xf64>
      vector.transfer_write %46, %alloc_23[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %alloc_24 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_1, %alloc_24[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %47 = vector.transfer_read %alloc_24[%0, %20], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %48 = vector.transfer_read %alloc_23[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %49 = vector.transfer_read %arg3[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %50 = nvgpu.mma.sync(%48, %49, %47) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %51 = vector.transfer_read %alloc_23[%0, %21], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %52 = vector.transfer_read %arg3[%0, %21], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %53 = nvgpu.mma.sync(%51, %52, %50) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %53, %alloc_24[%0, %20] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %alloc_25 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      vector.transfer_write %cst_1, %alloc_25[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %54 = vector.transfer_read %alloc_25[%0, %20], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x2xf64>
      %55 = vector.transfer_read %arg3[%0, %1], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %56 = vector.transfer_read %alloc_24[%1, %0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %57 = nvgpu.mma.sync(%55, %56, %54) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      %58 = vector.transfer_read %arg3[%0, %21], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %59 = vector.transfer_read %alloc_24[%21, %0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<1x1xf64>
      %60 = nvgpu.mma.sync(%58, %59, %57) {mmaShape = [8, 8, 4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
      vector.transfer_write %60, %alloc_25[%0, %20] {in_bounds = [true, true]} : vector<1x2xf64>, memref<8x8xf64>
      %61 = arith.addi %arg6, %c1 : index
      %subview_26 = memref.subview %alloc[0, 0, %arg6] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: ?>>
      %alloc_27 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %62 = vector.transfer_read %subview_26[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %63 = vector.transfer_read %alloc_25[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %64 = arith.subf %62, %63 : vector<8x8xf64>
      vector.transfer_write %64, %alloc_27[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %65 = vector.transfer_read %alloc_27[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      vector.transfer_write %65, %subview_26[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
      %subview_28 = memref.subview %alloc[0, 0, %61] [8, 8, 1] [1, 1, 1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[64, 8], offset: ?>>
      %alloc_29 = memref.alloc() {alignment = 64 : i64} : memref<8x8xf64>
      %66 = vector.transfer_read %subview_28[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64, strided<[64, 8], offset: ?>>, vector<8x8xf64>
      %67 = vector.transfer_read %alloc_25[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      %68 = arith.addf %66, %67 : vector<8x8xf64>
      vector.transfer_write %68, %alloc_29[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
      %69 = vector.transfer_read %alloc_29[%c0, %c0], %cst_4 {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
      vector.transfer_write %69, %subview_28[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64, strided<[64, 8], offset: ?>>
    }
    return %alloc : memref<8x8x8xf64>
  }
}

