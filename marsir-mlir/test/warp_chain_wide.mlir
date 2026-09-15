// A WIDE sum-factorization sweep: Btil(8x8) @ U(8x64) -> 8x64, the B-sweep shape
// that mlir_warp.py's emit_matmul hand-tiles into 8 column tiles of 8 in Python.
// m8n8k4 fixes m = 8 in hardware, but K splits into slabs of 4 and N into column
// tiles of 8, so this is 8 tiles x 2 slabs = 16 mma, no shuffles (both operands
// are leaf reads, nothing is relayed out of a C-fragment).
// RUN: mir-opt %s --mir-chain-contracts
#a = affine_map<(m, n, k) -> (m, k)>
#b = affine_map<(m, n, k) -> (k, n)>
#c = affine_map<(m, n, k) -> (m, n)>
gpu.module @mir_kernels [#nvvm.target<chip = "sm_90", O = 3>] {
  gpu.func @bsweep(%Bt: memref<8x8xf64>, %U: memref<8x64xf64>,
                   %Out: memref<8x64xf64>) kernel {
    %c0 = arith.constant 0 : index
    %z  = arith.constant 0.0 : f64
    %zc = arith.constant dense<0.0> : vector<8x64xf64>
    %va = vector.transfer_read %Bt[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x8xf64>
    %vb = vector.transfer_read %U[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<8x64xf64>
    %r = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>}
         %va, %vb, %zc : vector<8x8xf64>, vector<8x64xf64> into vector<8x64xf64>
    vector.transfer_write %r, %Out[%c0, %c0] {in_bounds=[true,true]} : vector<8x64xf64>, memref<8x64xf64>
    gpu.return
  }
}
