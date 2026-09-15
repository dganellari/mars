// The +/- PLANE SCATTER, the last piece of the face chain, expressed as plain
// 8x8 vector ops rather than hand-emitted per-lane code:
//   y[plane_a] -= flux ; y[plane_b] += flux
// The flux comes out of a contract, so it is already a C-fragment; the two
// accumulations are elementwise and fuse lane-locally; the two writes land at
// the planes' own base offsets. Nothing here is scatter-specific in the pass --
// it falls out of contract + elementwise + transfer_write.
// RUN: mir-opt %s --mir-chain-contracts
#a = affine_map<(m, n, k) -> (m, k)>
#b = affine_map<(m, n, k) -> (k, n)>
#c = affine_map<(m, n, k) -> (m, n)>
gpu.module @mir_kernels [#nvvm.target<chip = "sm_90", O = 3>] {
  gpu.func @scatter(%A: memref<8x8xf64>, %B: memref<8x8xf64>,
                    %Y: memref<8x64xf64>) kernel {
    %c0  = arith.constant 0 : index
    %c8  = arith.constant 8 : index
    %c16 = arith.constant 16 : index
    %z   = arith.constant 0.0 : f64
    %zc  = arith.constant dense<0.0> : vector<8x8xf64>
    %va = vector.transfer_read %A[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x8xf64>
    %vb = vector.transfer_read %B[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x8xf64>
    %fx = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>}
          %va, %vb, %zc : vector<8x8xf64>, vector<8x8xf64> into vector<8x8xf64>
    // minus plane at column 8, plus plane at column 16
    %ya = vector.transfer_read %Y[%c0, %c8],  %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<8x8xf64>
    %yb = vector.transfer_read %Y[%c0, %c16], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<8x8xf64>
    %na = arith.subf %ya, %fx : vector<8x8xf64>
    %nb = arith.addf %yb, %fx : vector<8x8xf64>
    vector.transfer_write %na, %Y[%c0, %c8]  {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64>
    vector.transfer_write %nb, %Y[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64>
    gpu.return
  }
}
