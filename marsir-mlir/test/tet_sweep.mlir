// The collapsed (Duffy/PKD) tet ragged sweep in high-level mir:
//   out[p,q,k] = sum_{r=0}^{D-p-q} u[p,q,r] * C[p,q,r,k]
// D=2, so the dense modal cube is W=D+1=3, with n=3 quadrature points. The
// reduction bound depends on the OUTPUT indices, which is what keeps tets off
// the hyperrectangular linalg path and out of the m8n8k4 tile.
//
// The whole tet path, as verified end to end:
//   mir-opt %s --convert-mir-to-linalg                      ragged scf.for nest
//   | mlir-opt --one-shot-bufferize="bufferize-function-boundaries \
//              function-boundary-type-conversion=identity-layout-map"
//   | mir-opt --mir-batch-elements                          warp-per-element
// (mark %u {mir.element} before the last step; %C is the reference table and is
// shared by every element, so it must NOT be batched.)
// Numerics for this op are gated by test/exec_gate.py gate 4.
// RUN: mir-opt %s --convert-mir-to-linalg
func.func @tet_sweep(%u: tensor<3x3x3xf64>, %C: tensor<3x3x3x3xf64>) -> tensor<3x3x3xf64> {
  %o = mir.simplex_contract %u, %C {degree = 2 : i64}
       : (tensor<3x3x3xf64>, tensor<3x3x3x3xf64>) -> tensor<3x3x3xf64>
  return %o : tensor<3x3x3xf64>
}
