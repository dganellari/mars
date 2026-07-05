// A high-order matrix-free PARTIAL-ASSEMBLY operator in the mir dialect:
//   y = B^T . D . B . u        (sum-factorized)
//   B   : mir.contract   -- evaluate to the quadrature grid (one direction shown)
//   D   : mir.flux       -- the pointwise per-quad action (Galerkin metric / CVFEM flux)
//   B^T : mir.contract   -- integrate back to modes
// This is the shape MARSIR targets for Galerkin AND CVFEM high-order matrix-free.
// RUN: mir-opt %s | mir-opt
func.func @pa_laplacian(%u: tensor<8x8x8xf64>, %D: tensor<8x8xf64>, %g: tensor<8x8x8xf64>)
    -> tensor<8x8x8xf64> {
  // B: sum-factorization sweep to the quadrature grid.
  %grad = mir.contract %u, %D {axis = 0 : i64}
          : (tensor<8x8x8xf64>, tensor<8x8xf64>) -> tensor<8x8x8xf64>

  // D: the authored flux -- here the pointwise metric action g * grad.
  %flux = mir.flux ins(%grad, %g) : (tensor<8x8x8xf64>, tensor<8x8x8xf64>) -> tensor<8x8x8xf64> {
  ^bb0(%gr: f64, %gm: f64):
    %f = arith.mulf %gm, %gr : f64
    mir.yield %f : f64
  }

  // B^T: integrate back to modal coefficients.
  %y = mir.contract %flux, %D {axis = 0 : i64}
       : (tensor<8x8x8xf64>, tensor<8x8xf64>) -> tensor<8x8x8xf64>
  return %y : tensor<8x8x8xf64>
}
