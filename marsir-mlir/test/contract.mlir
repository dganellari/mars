// Round-trip: mir-opt must parse and re-print the mir dialect.
// RUN: mir-opt %s | mir-opt

// One direction of a p=7 (n=8) tensor-product sum-factorization: apply the
// derivative operator D (8x8) along axis 0 of the n^3 element field.
func.func @sumfac_one_axis(%u: tensor<8x8x8xf64>, %D: tensor<8x8xf64>)
    -> tensor<8x8x8xf64> {
  %g = mir.contract %u, %D {axis = 0 : i64}
       : (tensor<8x8x8xf64>, tensor<8x8xf64>) -> tensor<8x8x8xf64>
  return %g : tensor<8x8x8xf64>
}
