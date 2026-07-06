// PROBE 1 experiment: two m8n8k4 f64 nvgpu.mma.sync chained REGISTER-RESIDENT.
// mma2's A operand is a gpu.shuffle relayout of mma1's C-fragment result.
// No shared memory. Launch: grid 1, block (32,1,1).
//
// C-frag lane L holds C[L/4, 2*(L%4)], C[L/4, 2*(L%4)+1]  (vector<1x2>)
// A-frag lane L holds A[L/4, L%4]                          (vector<1x1>)
// want A2[i,k] = D1[i,k], i=L/4, k=L%4 (k in 0..3 -> D1 cols 0..3):
//   src lane = 4*(L/4) + ((L%4)>>1) = (L & 28) | ((L & 3) >> 1)
//   element  = L & 1   (0 -> D1 lo, 1 -> D1 hi)
gpu.module @mir_kernels [#nvvm.target<chip = "sm_90", O = 3>] {
  gpu.func @chain_mma(%A1: memref<8x4xf64>, %B1: memref<4x8xf64>,
                      %B2: memref<4x8xf64>, %C: memref<8x8xf64>) kernel {
    %z  = arith.constant 0.0 : f64
    %c2 = arith.constant 2 : index
    %c4 = arith.constant 4 : index
    %lane = gpu.thread_id x
    %i = arith.divui %lane, %c4 : index      // L/4
    %k = arith.remui %lane, %c4 : index      // L%4

    // fragment reads for mma1
    %a1 = vector.transfer_read %A1[%i, %k], %z {in_bounds=[true,true]}
          : memref<8x4xf64>, vector<1x1xf64>
    %b1 = vector.transfer_read %B1[%k, %i], %z {in_bounds=[true,true]}
          : memref<4x8xf64>, vector<1x1xf64>
    %j0 = arith.muli %k, %c2 : index
    %cacc = vector.transfer_read %C[%i, %j0], %z {in_bounds=[true,true]}
          : memref<8x8xf64>, vector<1x2xf64>

    %d1 = nvgpu.mma.sync(%a1, %b1, %cacc) {mmaShape=[8,8,4]}
          : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>

    // ---- register-resident C-frag -> A-frag relayout via gpu.shuffle ----
    %li  = arith.index_cast %lane : index to i32
    %c28 = arith.constant 28 : i32
    %c3  = arith.constant 3  : i32
    %c1  = arith.constant 1  : i32
    %c32 = arith.constant 32 : i32
    %hi   = arith.andi %li, %c28 : i32       // 4*(L/4)
    %klo  = arith.andi %li, %c3  : i32       // L%4
    %kh   = arith.shrui %klo, %c1 : i32      // (L%4)>>1
    %src  = arith.ori %hi, %kh : i32         // source lane
    %d1lo = vector.extract %d1[0, 0] : f64 from vector<1x2xf64>
    %d1hi = vector.extract %d1[0, 1] : f64 from vector<1x2xf64>
    %s0, %v0 = gpu.shuffle idx %d1lo, %src, %c32 : f64
    %s1, %v1 = gpu.shuffle idx %d1hi, %src, %c32 : f64
    %par  = arith.andi %li, %c1 : i32        // L&1
    %odd  = arith.cmpi eq, %par, %c1 : i32
    %sel  = arith.select %odd, %s1, %s0 : f64
    %zero1 = arith.constant dense<0.0> : vector<1x1xf64>
    %a2   = vector.insert %sel, %zero1 [0, 0] : f64 into vector<1x1xf64>

    // mma2 consumes the shuffled A2; accumulate onto D1
    %b2 = vector.transfer_read %B2[%k, %i], %z {in_bounds=[true,true]}
          : memref<4x8xf64>, vector<1x1xf64>
    %d2 = nvgpu.mma.sync(%a2, %b2, %d1) {mmaShape=[8,8,4]}
          : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>

    vector.transfer_write %d2, %C[%i, %j0] {in_bounds=[true,true]}
          : vector<1x2xf64>, memref<8x8xf64>
    gpu.return
  }
}
