// The !mir.block type: carries FEM semantics (shape/order/deformed/batch) as a
// compile-time type. mir.unbind(mir.bind(x)) folds to x under canonicalize, so
// the semantic annotation erases before lowering.
// RUN: mir-opt %s | mir-opt
// RUN: mir-opt %s --canonicalize   (the bind/unbind pair must fold away)
func.func @annotate(%u: tensor<1024x8x8x8xf64>) -> tensor<1024x8x8x8xf64> {
  %b = mir.bind %u : tensor<1024x8x8x8xf64>
       -> !mir.block<shape = "hex", p = 7, deformed = 0, batch = 1024>
  %r = mir.unbind %b : !mir.block<shape = "hex", p = 7, deformed = 0, batch = 1024>
       -> tensor<1024x8x8x8xf64>
  return %r : tensor<1024x8x8x8xf64>
}
