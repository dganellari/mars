// A flux-like pointwise stage (g2*deriv + g0*dt2) on an 8x8 face inside a warp
// region. After --mir-warp-distribute each lane should own a 1x2 slice, the
// warp op dissolving into per-lane vector code.
// RUN: mir-opt %s --mir-warp-distribute
func.func @flux_stage(%lane: index, %A: memref<8x8xf64>, %B: memref<8x8xf64>,
                      %C: memref<8x8xf64>, %O: memref<8x8xf64>) {
  %c0 = arith.constant 0 : index
  %z = arith.constant 0.0 : f64
  vector.warp_execute_on_lane_0(%lane)[32] {
    %a = vector.transfer_read %A[%c0, %c0], %z {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
    %b = vector.transfer_read %B[%c0, %c0], %z {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
    %c = vector.transfer_read %C[%c0, %c0], %z {in_bounds = [true, true]} : memref<8x8xf64>, vector<8x8xf64>
    %m1 = arith.mulf %a, %b : vector<8x8xf64>
    %m2 = arith.mulf %c, %m1 : vector<8x8xf64>
    %s = arith.addf %m1, %m2 : vector<8x8xf64>
    vector.transfer_write %s, %O[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64>
  }
  return
}
