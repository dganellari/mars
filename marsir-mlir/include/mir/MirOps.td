#ifndef MIR_OPS
#define MIR_OPS

include "mlir/IR/AttrTypeBase.td"
include "mlir/IR/OpBase.td"
include "mlir/Interfaces/SideEffectInterfaces.td"

//===----------------------------------------------------------------------===//
// The mir dialect: MARS FEM intermediate representation.
//===----------------------------------------------------------------------===//

def Mir_Dialect : Dialect {
  let name = "mir";
  let summary = "MARS FEM intermediate representation (sum-factorization)";
  let description = [{
    High-level dialect for matrix-free finite-element operators. Encodes the
    tensor-product sum-factorization structure so it can be lowered
    progressively: mir -> linalg -> gpu/nvgpu (tensor core) -> nvvm -> PTX.
  }];
  let cppNamespace = "::mir";
  let useDefaultTypePrinterParser = 1;
}

class Mir_Op<string mnemonic, list<Trait> traits = []>
    : Op<Mir_Dialect, mnemonic, traits>;

//===----------------------------------------------------------------------===//
// !mir.block : a typed batch of alike elements (the FEM semantics carrier).
//===----------------------------------------------------------------------===//

def Mir_BlockType : TypeDef<Mir_Dialect, "Block"> {
  let mnemonic = "block";
  let summary = "a batch of alike finite elements";
  let description = [{
    Carries the FEM metadata of a batch of elements as a compile-time type:
    element shape ("hex"/"tet"), polynomial order p, deformed (1 = curved
    geometry, per-point metric; 0 = affine), and the batch size E (number of
    elements; -1 = dynamic). Ops bind element data to a block (mir.bind) so
    schedules can key on FEM semantics; lowering erases the pair
    (mir.unbind(mir.bind(x)) folds to x).

    Example: !mir.block<shape = "hex", p = 7, deformed = 0, batch = 1024>
  }];
  let parameters = (ins StringRefParameter<"element shape">:$shape,
                        "int64_t":$order,
                        "int64_t":$deformed,
                        "int64_t":$batch);
  let assemblyFormat =
      "`<` `shape` `=` $shape `,` `p` `=` $order `,` `deformed` `=` $deformed "
      "`,` `batch` `=` $batch `>`";
}

def Mir_BindOp : Mir_Op<"bind", [Pure]> {
  let summary = "bind raw element data to a typed block (attach FEM semantics)";
  let arguments = (ins AnyRankedTensor:$data);
  let results = (outs Mir_BlockType:$block);
  let assemblyFormat =
      "$data attr-dict `:` type($data) `->` qualified(type($block))";
}

def Mir_UnbindOp : Mir_Op<"unbind", [Pure]> {
  let summary = "recover the raw tensor from a block (folds away over mir.bind)";
  let arguments = (ins Mir_BlockType:$block);
  let results = (outs AnyRankedTensor:$data);
  let assemblyFormat =
      "$block attr-dict `:` qualified(type($block)) `->` type($data)";
  let hasFolder = 1;
}

//===----------------------------------------------------------------------===//
// mir.contract : the atom of tensor-product sum-factorization.
//===----------------------------------------------------------------------===//

def Mir_ContractOp : Mir_Op<"contract", [Pure]> {
  let summary = "1D reference-operator contraction along one axis";
  let description = [{
    Contracts the 1D reference operator matrix `op_matrix` against `input`
    along dimension `axis`. This is the single 1D pass of tensor-product
    sum-factorization (e.g. Btil/Dtil/D/W applied along x, y, or z). It lowers
    to a `linalg.matmul` (the unfold-along-axis matrix multiply), which is
    MLIR's on-ramp to GPU tensor-core codegen.

    Example:
      %g = mir.contract %u, %D {axis = 0 : i64}
             : (tensor<8x8x8xf64>, tensor<8x8xf64>) -> tensor<8x8x8xf64>
  }];

  let arguments = (ins AnyRankedTensor:$input,
                       AnyRankedTensor:$op_matrix,
                       I64Attr:$axis);
  let results = (outs AnyRankedTensor:$result);

  let assemblyFormat = [{
    $input `,` $op_matrix attr-dict `:`
      functional-type(operands, results)
  }];
}

//===----------------------------------------------------------------------===//
// mir.simplex_contract : the ragged (collapsed-tet) sum-factorization sweep.
//===----------------------------------------------------------------------===//

def Mir_SimplexContractOp : Mir_Op<"simplex_contract", [Pure]> {
  let summary = "ragged modal sweep of the collapsed (Duffy/PKD) tetrahedron";
  let description = [{
    The innermost sweep of collapsed-coordinate tet sum-factorization: for the
    simplex modal set {p+q+r <= D},

      out[p, q, k] = sum_{r = 0}^{D - p - q}  u[p, q, r] * table[p, q, r, k]

    Unlike mir.contract this is NOT hyperrectangular: the reduction bound
    depends on the outer (output) indices, and the 1D factor differs per (p, q)
    -- hence the 4-D operator `table` (the PKD C-factor evaluated at the n
    Gauss-Jacobi points). u and out use the dense W = D+1 cube layout with the
    ragged region valid (entries beyond p+q+r <= D are ignored / zeroed).
    Lowers to explicit scf.for loops (triangular nests do not fit linalg).
  }];

  let arguments = (ins AnyRankedTensor:$input,     // u[W, W, W]
                       AnyRankedTensor:$table,     // C[W, W, W, n]
                       I64Attr:$degree);           // D (order of the simplex set)
  let results = (outs AnyRankedTensor:$result);    // out[W, W, n]

  let assemblyFormat = [{
    $input `,` $table attr-dict `:` functional-type(operands, results)
  }];
}

//===----------------------------------------------------------------------===//
// mir.flux : the pointwise per-quadrature-point action (the "D" of B^T D B).
// Named for what you author in a MARSIR .op spec -- the flux.
//===----------------------------------------------------------------------===//

def Mir_FluxOp : Mir_Op<"flux",
    [SingleBlockImplicitTerminator<"YieldOp">, RecursiveMemoryEffects]> {
  let summary = "pointwise per-quadrature-point action (the D of B^T D B)";
  let description = [{
    The operator-specific pointwise action applied at every quadrature point --
    the "D" of matrix-free partial assembly, and exactly the `flux` a user writes
    in a MARSIR .op spec. The B / B^T sweeps around it are `mir.contract`; only
    this op differs between operators:

      * Galerkin Laplacian: flux = G_c . grad  (the collapse-metric contraction)
      * CVFEM (Knaus SCS):  flux = g2*deriv + g0*dt2 + g1*dt1

    Its region receives the per-point input quantities (gradient components,
    metric, ...) and yields the per-point output(s). Lowers to a `linalg.generic`
    elementwise map over the quadrature grid, and fuses with the surrounding
    contractions into one loop.
  }];

  let arguments = (ins Variadic<AnyRankedTensor>:$inputs);
  let results = (outs Variadic<AnyRankedTensor>:$results);
  let regions = (region SizedRegion<1>:$body);

  let assemblyFormat = [{
    attr-dict `ins` `(` $inputs `)` `:` functional-type($inputs, $results) $body
  }];
}

def Mir_YieldOp : Mir_Op<"yield", [Pure, Terminator, HasParent<"FluxOp">]> {
  let summary = "terminator yielding the per-point result of a mir.qfunction";
  let arguments = (ins Variadic<AnyType>:$operands);
  let assemblyFormat = "attr-dict ($operands^ `:` type($operands))?";
  // Zero-operand builder (defined in MirOps.cpp) so SingleBlockImplicitTerminator
  // can insert the implicit terminator.
  let builders = [OpBuilder<(ins)>];
}

#endif // MIR_OPS
