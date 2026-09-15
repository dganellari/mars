// THE FULL REGISTER-RESIDENT FACE CHAIN, as plain 8x8 vector ops -- the thing
// mlir_warp_reg.py hand-writes per lane today:
//   deriv = D @ interp ; dt1 = D @ t1 ; dt2 = D @ t2      (3 contracts)
//   flux  = g2*deriv + g0*dt2 + g1*dt1                    (the authored line)
//   intf  = W @ flux                                      (integrate back)
//   y[-plane] -= intf ; y[+plane] += intf                 (scatter)
// 4 contracts x 2 k-slabs = 8 mma. flux is elementwise so it stays lane-local on
// the vector<1x2> C-fragment; only intf's operand needs a C->B shuffle relayout.
// RUN: mir-opt %s --mir-chain-contracts
#a = affine_map<(m, n, k) -> (m, k)>
#b = affine_map<(m, n, k) -> (k, n)>
#c = affine_map<(m, n, k) -> (m, n)>
gpu.module @mir_kernels [#nvvm.target<chip = "sm_90", O = 3>] {
  gpu.func @face(%D: memref<8x8xf64>, %W: memref<8x8xf64>,
                 %interp: memref<8x8xf64>, %t1: memref<8x8xf64>, %t2: memref<8x8xf64>,
                 %g0: memref<8x8xf64>, %g1: memref<8x8xf64>, %g2: memref<8x8xf64>,
                 %Y: memref<8x64xf64>) kernel {
    %c0 = arith.constant 0 : index
    %c8 = arith.constant 8 : index
    %c16 = arith.constant 16 : index
    %z  = arith.constant 0.0 : f64
    %zc = arith.constant dense<0.0> : vector<8x8xf64>
    %vD = vector.transfer_read %D[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x8xf64>
    %vW = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x8xf64>
    %vi = vector.transfer_read %interp[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x8xf64>
    %v1 = vector.transfer_read %t1[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x8xf64>
    %v2 = vector.transfer_read %t2[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x8xf64>
    %m0 = vector.transfer_read %g0[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x8xf64>
    %m1 = vector.transfer_read %g1[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x8xf64>
    %m2 = vector.transfer_read %g2[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x8xf64>

    %deriv = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>}
             %vD, %vi, %zc : vector<8x8xf64>, vector<8x8xf64> into vector<8x8xf64>
    %dt1 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>}
           %vD, %v1, %zc : vector<8x8xf64>, vector<8x8xf64> into vector<8x8xf64>
    %dt2 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>}
           %vD, %v2, %zc : vector<8x8xf64>, vector<8x8xf64> into vector<8x8xf64>

    %p2 = arith.mulf %m2, %deriv : vector<8x8xf64>
    %p0 = arith.mulf %m0, %dt2 : vector<8x8xf64>
    %p1 = arith.mulf %m1, %dt1 : vector<8x8xf64>
    %s0 = arith.addf %p2, %p0 : vector<8x8xf64>
    %flux = arith.addf %s0, %p1 : vector<8x8xf64>

    %intf = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>}
            %vW, %flux, %zc : vector<8x8xf64>, vector<8x8xf64> into vector<8x8xf64>

    %ya = vector.transfer_read %Y[%c0, %c8],  %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<8x8xf64>
    %yb = vector.transfer_read %Y[%c0, %c16], %z {in_bounds=[true,true]} : memref<8x64xf64>, vector<8x8xf64>
    %na = arith.subf %ya, %intf : vector<8x8xf64>
    %nb = arith.addf %yb, %intf : vector<8x8xf64>
    vector.transfer_write %na, %Y[%c0, %c8]  {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64>
    vector.transfer_write %nb, %Y[%c0, %c16] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x64xf64>
    gpu.return
  }
}
