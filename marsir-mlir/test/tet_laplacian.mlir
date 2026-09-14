// The collapsed (Duffy/PKD) tet Laplacian  y = B^T G_c B u, sum-factorized, in
// high-level mir. D=3 so the modal cube is W=4 and there are n=4 Gauss-Jacobi
// points per direction. This is the tet analogue of test/pa_laplacian.mlir.
//
// Per gradient component the forward half is three sweeps, of which two are
// RAGGED (the reduction bound depends on the output indices, which is why they
// are mir.simplex_contract and not linalg):
//   r -> k   ragged, bound D-p-q     mir.simplex_contract axis=2
//   q -> j   ragged, bound D-p       mir.simplex_contract axis=1
//   p -> i   full range              mir.contract axis=0, with A transposed
// then the 3x3 collapse metric maps gradient to flux pointwise (mir.flux), and
// the transposed sweeps integrate back to modal coefficients.
//
// Component c differentiates factor c: (Ad,B,C), (A,Bd,C), (A,B,Cd).
// Numerics are gated by test/exec_gate.py gate 8, against tet_galerkin.py.
// RUN: mir-opt %s --convert-mir-to-linalg
func.func @tet_laplacian(
    %u: tensor<4x4x4xf64>,
    %A: tensor<4x4xf64>, %Ad: tensor<4x4xf64>,      // A[p,i]  (integrate back)
    %AT: tensor<4x4xf64>, %AdT: tensor<4x4xf64>,    // A^T[i,p] (evaluate)
    %B: tensor<4x4x4xf64>, %Bd: tensor<4x4x4xf64>,
    %C: tensor<4x4x4x4xf64>, %Cd: tensor<4x4x4x4xf64>,
    %g00: tensor<4x4x4xf64>, %g01: tensor<4x4x4xf64>, %g02: tensor<4x4x4xf64>,
    %g10: tensor<4x4x4xf64>, %g11: tensor<4x4x4xf64>, %g12: tensor<4x4x4xf64>,
    %g20: tensor<4x4x4xf64>, %g21: tensor<4x4x4xf64>, %g22: tensor<4x4x4xf64>)
    -> tensor<4x4x4xf64> {
  // component 0: d/dr1 -> Ad
  %f1_0 = mir.simplex_contract %u, %C {degree = 3 : i64, axis = 2 : i64}
        : (tensor<4x4x4xf64>, tensor<4x4x4x4xf64>) -> tensor<4x4x4xf64>
  %f2_0 = mir.simplex_contract %f1_0, %B {degree = 3 : i64, axis = 1 : i64}
        : (tensor<4x4x4xf64>, tensor<4x4x4xf64>) -> tensor<4x4x4xf64>
  %gr0 = mir.contract %f2_0, %AdT {axis = 0 : i64}
        : (tensor<4x4x4xf64>, tensor<4x4xf64>) -> tensor<4x4x4xf64>
  // component 1: d/dr2 -> Bd
  %f1_1 = mir.simplex_contract %u, %C {degree = 3 : i64, axis = 2 : i64}
        : (tensor<4x4x4xf64>, tensor<4x4x4x4xf64>) -> tensor<4x4x4xf64>
  %f2_1 = mir.simplex_contract %f1_1, %Bd {degree = 3 : i64, axis = 1 : i64}
        : (tensor<4x4x4xf64>, tensor<4x4x4xf64>) -> tensor<4x4x4xf64>
  %gr1 = mir.contract %f2_1, %AT {axis = 0 : i64}
        : (tensor<4x4x4xf64>, tensor<4x4xf64>) -> tensor<4x4x4xf64>
  // component 2: d/dr3 -> Cd
  %f1_2 = mir.simplex_contract %u, %Cd {degree = 3 : i64, axis = 2 : i64}
        : (tensor<4x4x4xf64>, tensor<4x4x4x4xf64>) -> tensor<4x4x4xf64>
  %f2_2 = mir.simplex_contract %f1_2, %B {degree = 3 : i64, axis = 1 : i64}
        : (tensor<4x4x4xf64>, tensor<4x4x4xf64>) -> tensor<4x4x4xf64>
  %gr2 = mir.contract %f2_2, %AT {axis = 0 : i64}
        : (tensor<4x4x4xf64>, tensor<4x4xf64>) -> tensor<4x4x4xf64>

  // flux_a = sum_b Gc[a,b] * grad_b   -- the one line a user authors
  %fx0 = mir.flux ins(%gr0, %gr1, %gr2, %g00, %g01, %g02)
       : (tensor<4x4x4xf64>, tensor<4x4x4xf64>, tensor<4x4x4xf64>,
          tensor<4x4x4xf64>, tensor<4x4x4xf64>, tensor<4x4x4xf64>) -> tensor<4x4x4xf64> {
  ^bb0(%d0: f64, %d1: f64, %d2: f64, %m0: f64, %m1: f64, %m2: f64):
    %a0 = arith.mulf %m0, %d0 : f64
    %a1 = arith.mulf %m1, %d1 : f64
    %a2 = arith.mulf %m2, %d2 : f64
    %s0 = arith.addf %a0, %a1 : f64
    %s1 = arith.addf %s0, %a2 : f64
    mir.yield %s1 : f64
  }
  %fx1 = mir.flux ins(%gr0, %gr1, %gr2, %g10, %g11, %g12)
       : (tensor<4x4x4xf64>, tensor<4x4x4xf64>, tensor<4x4x4xf64>,
          tensor<4x4x4xf64>, tensor<4x4x4xf64>, tensor<4x4x4xf64>) -> tensor<4x4x4xf64> {
  ^bb0(%d0: f64, %d1: f64, %d2: f64, %m0: f64, %m1: f64, %m2: f64):
    %a0 = arith.mulf %m0, %d0 : f64
    %a1 = arith.mulf %m1, %d1 : f64
    %a2 = arith.mulf %m2, %d2 : f64
    %s0 = arith.addf %a0, %a1 : f64
    %s1 = arith.addf %s0, %a2 : f64
    mir.yield %s1 : f64
  }
  %fx2 = mir.flux ins(%gr0, %gr1, %gr2, %g20, %g21, %g22)
       : (tensor<4x4x4xf64>, tensor<4x4x4xf64>, tensor<4x4x4xf64>,
          tensor<4x4x4xf64>, tensor<4x4x4xf64>, tensor<4x4x4xf64>) -> tensor<4x4x4xf64> {
  ^bb0(%d0: f64, %d1: f64, %d2: f64, %m0: f64, %m1: f64, %m2: f64):
    %a0 = arith.mulf %m0, %d0 : f64
    %a1 = arith.mulf %m1, %d1 : f64
    %a2 = arith.mulf %m2, %d2 : f64
    %s0 = arith.addf %a0, %a1 : f64
    %s1 = arith.addf %s0, %a2 : f64
    mir.yield %s1 : f64
  }

  // integrate back: full p sweep, then the two ragged transposes
  %h1_0 = mir.contract %fx0, %Ad {axis = 0 : i64}
        : (tensor<4x4x4xf64>, tensor<4x4xf64>) -> tensor<4x4x4xf64>
  %h2_0 = mir.simplex_contract %h1_0, %B {degree = 3 : i64, axis = 1 : i64, transposed = true}
        : (tensor<4x4x4xf64>, tensor<4x4x4xf64>) -> tensor<4x4x4xf64>
  %y0 = mir.simplex_contract %h2_0, %C {degree = 3 : i64, axis = 2 : i64, transposed = true}
        : (tensor<4x4x4xf64>, tensor<4x4x4x4xf64>) -> tensor<4x4x4xf64>
  %h1_1 = mir.contract %fx1, %A {axis = 0 : i64}
        : (tensor<4x4x4xf64>, tensor<4x4xf64>) -> tensor<4x4x4xf64>
  %h2_1 = mir.simplex_contract %h1_1, %Bd {degree = 3 : i64, axis = 1 : i64, transposed = true}
        : (tensor<4x4x4xf64>, tensor<4x4x4xf64>) -> tensor<4x4x4xf64>
  %y1 = mir.simplex_contract %h2_1, %C {degree = 3 : i64, axis = 2 : i64, transposed = true}
        : (tensor<4x4x4xf64>, tensor<4x4x4x4xf64>) -> tensor<4x4x4xf64>
  %h1_2 = mir.contract %fx2, %A {axis = 0 : i64}
        : (tensor<4x4x4xf64>, tensor<4x4xf64>) -> tensor<4x4x4xf64>
  %h2_2 = mir.simplex_contract %h1_2, %B {degree = 3 : i64, axis = 1 : i64, transposed = true}
        : (tensor<4x4x4xf64>, tensor<4x4x4xf64>) -> tensor<4x4x4xf64>
  %y2 = mir.simplex_contract %h2_2, %Cd {degree = 3 : i64, axis = 2 : i64, transposed = true}
        : (tensor<4x4x4xf64>, tensor<4x4x4x4xf64>) -> tensor<4x4x4xf64>

  %y = mir.flux ins(%y0, %y1, %y2)
     : (tensor<4x4x4xf64>, tensor<4x4x4xf64>, tensor<4x4x4xf64>) -> tensor<4x4x4xf64> {
  ^bb0(%b0: f64, %b1: f64, %b2: f64):
    %t0 = arith.addf %b0, %b1 : f64
    %t1 = arith.addf %t0, %b2 : f64
    mir.yield %t1 : f64
  }
  return %y : tensor<4x4x4xf64>
}
