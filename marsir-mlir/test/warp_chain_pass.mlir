#a  = affine_map<(m, n, k) -> (m, k)>
#b  = affine_map<(m, n, k) -> (k, n)>
#c  = affine_map<(m, n, k) -> (m, n)>
gpu.module @mir_kernels [#nvvm.target<chip = "sm_90", O = 3>] {
  gpu.func @chain2(%A: memref<8x8xf64>, %B1: memref<8x8xf64>,
                   %B2: memref<8x8xf64>, %Out: memref<8x8xf64>) kernel {
    %c0 = arith.constant 0 : index
    %z  = arith.constant 0.0 : f64
    %zc = arith.constant dense<0.0> : vector<8x8xf64>
    %va  = vector.transfer_read %A[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x8xf64>
    %vb1 = vector.transfer_read %B1[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x8xf64>
    %vb2 = vector.transfer_read %B2[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x8xf64>
    %t = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>}
         %va, %vb1, %zc : vector<8x8xf64>, vector<8x8xf64> into vector<8x8xf64>
    %r = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>}
         %t, %vb2, %zc : vector<8x8xf64>, vector<8x8xf64> into vector<8x8xf64>
    vector.transfer_write %r, %Out[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64>
    gpu.return
  }
}
