gpu.module @mir_kernels [#nvvm.target<chip = "sm_90", O = 3>] {
  gpu.func @face_reg(%interp: memref<8x8xf64>, %deriv: memref<8x8xf64>, %Dm: memref<8x8xf64>, %W: memref<8x8xf64>, %g0: memref<8x8xf64>, %g1: memref<8x8xf64>, %g2: memref<8x8xf64>, %intf: memref<8x8xf64>) kernel {
    %z = arith.constant 0.0 : f64
    %zc1 = arith.constant dense<0.0> : vector<1x1xf64>
    %zc2 = arith.constant dense<0.0> : vector<1x2xf64>
    %c2i = arith.constant 2 : index
    %c4i = arith.constant 4 : index
    %lane = gpu.thread_id x
    %i = arith.divui %lane, %c4i : index
    %k = arith.remui %lane, %c4i : index
    %k4 = arith.addi %k, %c4i : index
    %tc = arith.muli %k, %c2i : index
    %li = arith.index_cast %lane : index to i32
    %c1 = arith.constant 1 : i32
    %c2 = arith.constant 2 : i32
    %c3 = arith.constant 3 : i32
    %c4 = arith.constant 4 : i32
    %c16 = arith.constant 16 : i32
    %c28 = arith.constant 28 : i32
    %c32 = arith.constant 32 : i32
    %f1 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f2 = vector.transfer_read %interp[%k, %i], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m3 = nvgpu.mma.sync(%f1, %f2, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f4 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f5 = vector.transfer_read %interp[%k4, %i], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m6 = nvgpu.mma.sync(%f4, %f5, %m3) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f7 = vector.transfer_read %interp[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f8 = vector.transfer_read %Dm[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m9 = nvgpu.mma.sync(%f7, %f8, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f10 = vector.transfer_read %interp[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %f11 = vector.transfer_read %Dm[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m12 = nvgpu.mma.sync(%f10, %f11, %m9) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %cf13 = vector.transfer_read %deriv[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x2xf64>
    %cf14 = vector.transfer_read %g2[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x2xf64>
    %cf15 = vector.transfer_read %g0[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x2xf64>
    %cf16 = vector.transfer_read %g1[%i, %tc], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x2xf64>
    %x17 = arith.mulf %cf14, %cf13 : vector<1x2xf64>
    %x18 = arith.mulf %cf15, %m12 : vector<1x2xf64>
    %x19 = arith.mulf %cf16, %m6 : vector<1x2xf64>
    %x20 = arith.addf %x17, %x18 : vector<1x2xf64>
    %fl21 = arith.addf %x20, %x19 : vector<1x2xf64>
    %lo22 = vector.extract %fl21[0, 0] : f64 from vector<1x2xf64>
    %hi23 = vector.extract %fl21[0, 1] : f64 from vector<1x2xf64>
    %t25 = arith.andi %li, %c28 : i32
    %t26 = arith.andi %li, %c3 : i32
    %t27 = arith.shrui %t26, %c1 : i32
    %t28 = arith.ori %t25, %t27 : i32
    %s29, %unused29 = gpu.shuffle idx %lo22, %t28, %c32 : f64
    %s31, %unused31 = gpu.shuffle idx %hi23, %t28, %c32 : f64
    %t33 = arith.andi %li, %c1 : i32
    %t34 = arith.cmpi eq, %t33, %c1 : i32
    %sv35 = arith.select %t34, %s31, %s29 : f64
    %a36 = vector.insert %sv35, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f37 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m38 = nvgpu.mma.sync(%a36, %f37, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %lo39 = vector.extract %fl21[0, 0] : f64 from vector<1x2xf64>
    %hi40 = vector.extract %fl21[0, 1] : f64 from vector<1x2xf64>
    %t42 = arith.andi %li, %c28 : i32
    %t43 = arith.andi %li, %c3 : i32
    %t44 = arith.shrui %t43, %c1 : i32
    %t45 = arith.ori %t42, %t44 : i32
    %t46 = arith.ori %t45, %c2 : i32
    %s47, %unused47 = gpu.shuffle idx %lo39, %t46, %c32 : f64
    %s49, %unused49 = gpu.shuffle idx %hi40, %t46, %c32 : f64
    %t51 = arith.andi %li, %c1 : i32
    %t52 = arith.cmpi eq, %t51, %c1 : i32
    %sv53 = arith.select %t52, %s49, %s47 : f64
    %a54 = vector.insert %sv53, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %f55 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %m56 = nvgpu.mma.sync(%a54, %f55, %m38) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f57 = vector.transfer_read %W[%i, %k], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo58 = vector.extract %m56[0, 0] : f64 from vector<1x2xf64>
    %hi59 = vector.extract %m56[0, 1] : f64 from vector<1x2xf64>
    %t60 = arith.andi %li, %c3 : i32
    %t61 = arith.muli %t60, %c4 : i32
    %t62 = arith.shrui %li, %c3 : i32
    %t63 = arith.addi %t61, %t62 : i32
    %s64, %unused64 = gpu.shuffle idx %lo58, %t63, %c32 : f64
    %s66, %unused66 = gpu.shuffle idx %hi59, %t63, %c32 : f64
    %t68 = arith.shrui %li, %c2 : i32
    %t69 = arith.andi %t68, %c1 : i32
    %t70 = arith.cmpi eq, %t69, %c1 : i32
    %sv71 = arith.select %t70, %s66, %s64 : f64
    %b72 = vector.insert %sv71, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m73 = nvgpu.mma.sync(%f57, %b72, %zc2) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    %f74 = vector.transfer_read %W[%i, %k4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<1x1xf64>
    %lo75 = vector.extract %m56[0, 0] : f64 from vector<1x2xf64>
    %hi76 = vector.extract %m56[0, 1] : f64 from vector<1x2xf64>
    %t77 = arith.andi %li, %c3 : i32
    %t78 = arith.muli %t77, %c4 : i32
    %t79 = arith.shrui %li, %c3 : i32
    %t80 = arith.addi %t78, %t79 : i32
    %t81 = arith.addi %t80, %c16 : i32
    %s82, %unused82 = gpu.shuffle idx %lo75, %t81, %c32 : f64
    %s84, %unused84 = gpu.shuffle idx %hi76, %t81, %c32 : f64
    %t86 = arith.shrui %li, %c2 : i32
    %t87 = arith.andi %t86, %c1 : i32
    %t88 = arith.cmpi eq, %t87, %c1 : i32
    %sv89 = arith.select %t88, %s84, %s82 : f64
    %b90 = vector.insert %sv89, %zc1 [0, 0] : f64 into vector<1x1xf64>
    %m91 = nvgpu.mma.sync(%f74, %b90, %m73) {mmaShape=[8,8,4]} : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
    vector.transfer_write %m91, %intf[%i, %tc] {in_bounds=[true,true]} : vector<1x2xf64>, memref<8x8xf64>
    gpu.return
  }
}
