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

    `axis` selects the modal axis being swept and `transposed` the direction
    (modal->quadrature, or quadrature->modal for integrating back). The four
    combinations are the four ragged stages of the collapsed tet operator; the
    p-axis sweep is full-range and is a plain `mir.contract`:

      axis=2               out[p,q,k] = sum_{r<=D-p-q} in[p,q,r] * T[p,q,r,k]
        in[W,W,W] T[W,W,W,n] -> out[W,W,n]
      axis=2 transposed    out[p,q,r] = sum_k in[p,q,k] * T[p,q,r,k],  r<=D-p-q
        in[W,W,n] T[W,W,W,n] -> out[W,W,W]
      axis=1               out[p,j,k] = sum_{q<=D-p} T[p,q,j] * in[p,q,k]
        in[W,W,n] T[W,W,n]   -> out[W,n,n]
      axis=1 transposed    out[p,q,k] = sum_j T[p,q,j] * in[p,j,k],    q<=D-p
        in[W,n,n] T[W,W,n]   -> out[W,W,n]
  }];

  let arguments = (ins AnyRankedTensor:$input,
                       AnyRankedTensor:$table,
                       I64Attr:$degree,                            // D
                       DefaultValuedAttr<I64Attr, "2">:$axis,      // modal axis: 1=q, 2=r
                       DefaultValuedAttr<BoolAttr, "false">:$transposed);
  let results = (outs AnyRankedTensor:$result);

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

def Mir_YieldOp : Mir_Op<"yield", [Pure, Terminator, ParentOneOf<["FluxOp", "DssOp", "CellwiseOp"]>]> {
  let summary = "terminator of mir.flux (per-point result), mir.dss and mir.cellwise (per-copy results)";
  let arguments = (ins Variadic<AnyType>:$operands);
  let assemblyFormat = "attr-dict ($operands^ `:` type($operands))?";
  // Zero-operand builder (defined in MirOps.cpp) so SingleBlockImplicitTerminator
  // can insert the implicit terminator.
  let builders = [OpBuilder<(ins)>];
}

//===----------------------------------------------------------------------===//
// mir.dss : direct stiffness summation on element-local (cell-wise) vectors.
//===----------------------------------------------------------------------===//

def Mir_DssOp : Mir_Op<"dss", [AttrSizedOperandSegments]> {
  let summary = "direct stiffness summation on element-local hex vectors";
  let description = [{
    Every element-local copy of a node gets the sum over all copies of that
    node: the assembly G G^T of a cell-wise solver, without forming the global
    vector. Edge and vertex copies are summed in the order of their star (the
    canonical order of the tables), and a face adds its two copies, so every
    copy of a node gets the same bits. Only the copies of the elements listed
    in `elements` are written; nodes are numbered (a, b, c) -> (a n + b) n + c
    with n = degree + 1.

    Copy i reads values[i] when i < size(values), else ghost[i - size(values)]:
    the copies of ghost elements, received from other ranks. With `minus`, a
    local copy reads values[i] - minus[i] (the DSS of a residual b - y; the
    ghost copies then hold the exchanged difference).

    The topology is that of MARS's UnstructuredTopology
    (backend/distributed/unstructured/solvers/mars_cellwise_topology.hpp):
      face_nbr[e*6+f]     the other copy of face f (element*6 + face), or -1
      face_code[e*6+f]    own frame | other frame << 3 (bit 6: Dirichlet)
      edge_of[e*12+k]     global edge of local edge k; bit 31: reversed
      edge_off, edge_ent  edge stars (CSR); entries element*12 + edge, bit 31:
                          reversed
      vert_of[e*8+c], vert_off, vert_ent
                          corner stars; entries element*8 + corner

    Without an epilogue, outs is one array that receives the sums.

    With an epilogue (what the hand kernels pass to for_each_copy_sum as
    `visit`), the region runs for every written copy t with
      (sum, dirichlet: i1, counted: i1, copies[0][t], copies[1][t], ...)
    and yields one value per `outs` array (stored at t) followed by
    `reductions` terms. Each term is summed over the copies of an element, and
    the sums go to partials[w * reductions + j] for the w-th listed element:
    reductions without atomics, in a fixed order. `dirichlet` and `counted`
    need flags(counting, edge_dirichlet, vert_dirichlet): counting[e] has bit b
    set when element e holds the counting copy of its entity with bit b (faces
    0-5, edges 8-19, corners 20-27, interior 31), so a sum over counted copies
    takes every node once.

    Example (the Jacobi preconditioner with <a, z> and <z, z>):
      mir.dss ins(%q, %ghost) elements(%ids) faces(%fn, %fc)
              edges(%eo, %eoff, %eent) vertices(%vo, %voff, %vent)
              flags(%cnt, %edir, %vdir) copies(%diag, %a) outs(%z)
              partials(%part) {degree = 7 : i64, reductions = 2 : i64}
              : memref<?xf64>, ... epilogue {
      ^bb0(%s: f64, %d: i1, %c: i1, %dg: f64, %a: f64):
        ... mir.yield %z, %az, %zz : f64, f64, f64
      }
  }];

  let arguments = (ins
      Arg<MemRefRankOf<[F64], [1]>, "local copies", [MemRead]>:$values,
      Arg<MemRefRankOf<[F64], [1]>, "ghost copies", [MemRead]>:$ghost,
      Arg<Optional<MemRefRankOf<[F64], [1]>>, "subtracted from the local copies", [MemRead]>:$minus,
      Arg<MemRefRankOf<[I32], [1]>, "elements to write", [MemRead]>:$elements,
      Arg<MemRefRankOf<[I32], [1]>, "", [MemRead]>:$face_nbr,
      Arg<MemRefRankOf<[I8], [1]>, "", [MemRead]>:$face_code,
      Arg<MemRefRankOf<[I32], [1]>, "", [MemRead]>:$edge_of,
      Arg<MemRefRankOf<[I32], [1]>, "", [MemRead]>:$edge_off,
      Arg<MemRefRankOf<[I32], [1]>, "", [MemRead]>:$edge_ent,
      Arg<MemRefRankOf<[I32], [1]>, "", [MemRead]>:$vert_of,
      Arg<MemRefRankOf<[I32], [1]>, "", [MemRead]>:$vert_off,
      Arg<MemRefRankOf<[I32], [1]>, "", [MemRead]>:$vert_ent,
      Arg<Variadic<AnyMemRef>, "counting, edge_dirichlet, vert_dirichlet", [MemRead]>:$flags,
      Arg<Variadic<MemRefRankOf<[F64], [1]>>, "per-copy epilogue inputs", [MemRead]>:$copies,
      Arg<Variadic<MemRefRankOf<[F64], [1]>>, "per-copy results", [MemWrite]>:$outs,
      Arg<Optional<MemRefRankOf<[F64], [1]>>, "per-element reduction partials", [MemWrite]>:$partials,
      I64Attr:$degree,
      DefaultValuedAttr<I64Attr, "0">:$reductions);
  let regions = (region AnyRegion:$epilogue);

  let assemblyFormat = [{
    `ins` `(` $values `,` $ghost `)` (`minus` `(` $minus^ `)`)? `elements` `(` $elements `)`
    `faces` `(` $face_nbr `,` $face_code `)`
    `edges` `(` $edge_of `,` $edge_off `,` $edge_ent `)`
    `vertices` `(` $vert_of `,` $vert_off `,` $vert_ent `)`
    (`flags` `(` $flags^ `)`)? (`copies` `(` $copies^ `)`)?
    `outs` `(` $outs `)` (`partials` `(` $partials^ `)`)?
    attr-dict `:` type(operands) (`epilogue` $epilogue^)?
  }];
  let hasVerifier = 1;
}

//===----------------------------------------------------------------------===//
// mir.cellwise : a pointwise map over element-local copies, with reductions.
//===----------------------------------------------------------------------===//

def Mir_CellwiseOp : Mir_Op<"cellwise", [AttrSizedOperandSegments]> {
  let summary = "per-copy map over the element-local vectors of listed elements";
  let description = [{
    The vector operations of a cell-wise Krylov solver: for every copy t of the
    elements in `elements` (nodes numbered as for mir.dss), the region runs on
      (counted: i1, inputs[0][t], inputs[1][t], ...)
    and yields one value per `outs` array (stored at t) followed by `reductions`
    terms, summed over the element's copies into partials[w * reductions + j]
    for the w-th listed element. `counted` is bit b of counting[e] for the
    copy's entity bit b (see mir.dss), so a sum over counted copies takes every
    node once; without `counting` it is true. With
    dirichlet(face_code, edge_of, vert_of, edge_dirichlet, vert_dirichlet) (the
    tables of mir.dss) the region takes a second flag, (counted, dirichlet, ...):
    whether the copy's node is a Dirichlet node.

    An array may be both an input and an output (x = x + ...): each copy is
    read and written by the same iteration. Values from above the op (e.g.
    Krylov scalars loaded once) may be used in the region.

    Example (s = r - alpha v):
      mir.cellwise elements(%ids) inputs(%r, %v) outs(%s) {degree = 7 : i64}
                   : memref<?xi32>, ... body {
      ^bb0(%c: i1, %ri: f64, %vi: f64):
        %av = arith.mulf %alpha, %vi : f64
        %si = arith.subf %ri, %av : f64
        mir.yield %si : f64
      }
  }];

  let arguments = (ins
      Arg<MemRefRankOf<[I32], [1]>, "elements to visit", [MemRead]>:$elements,
      Arg<Optional<MemRefRankOf<[I32], [1]>>, "counting bits per element", [MemRead]>:$counting,
      Arg<Variadic<AnyMemRef>, "face_code, edge_of, vert_of, edge_dirichlet, vert_dirichlet", [MemRead]>:$dirichlet,
      Arg<Variadic<MemRefRankOf<[F64], [1]>>, "per-copy inputs", [MemRead]>:$inputs,
      Arg<Variadic<MemRefRankOf<[F64], [1]>>, "per-copy results", [MemWrite]>:$outs,
      Arg<Optional<MemRefRankOf<[F64], [1]>>, "per-element reduction partials", [MemWrite]>:$partials,
      I64Attr:$degree,
      DefaultValuedAttr<I64Attr, "0">:$reductions);
  let regions = (region SizedRegion<1>:$body);

  let assemblyFormat = [{
    `elements` `(` $elements `)` (`counting` `(` $counting^ `)`)? (`dirichlet` `(` $dirichlet^ `)`)?
    (`inputs` `(` $inputs^ `)`)? (`outs` `(` $outs^ `)`)? (`partials` `(` $partials^ `)`)?
    attr-dict `:` type(operands) `body` $body
  }];
  let hasVerifier = 1;
}

//===----------------------------------------------------------------------===//
// mir.interpolate : M (x) M (x) M on element-local vectors (p-multigrid transfers).
//===----------------------------------------------------------------------===//

def Mir_InterpolateOp : Mir_Op<"interpolate"> {
  let summary = "per-element tensor-product map between two nodal degrees";
  let description = [{
    For every listed element e:
      out[e][a][b][c] (+)= sum_{j,k,m} M[a][j] M[b][k] M[c][m] in[e][j][k][m]
    with M the row-major `rows` x `cols` matrix (rows = nodes per axis of the
    output, cols = of the input). Prolongation passes the GLL interpolation
    matrix I (fine x coarse), restriction its transpose. `accumulate` adds to
    out instead of overwriting it.

    Example (prolong p = 3 -> 7, x += I u):
      mir.interpolate elements(%ids) matrix(%I) ins(%u) outs(%x)
                      {rows = 8 : i64, cols = 4 : i64, accumulate}
                      : memref<?xi32>, memref<?xf64>, memref<?xf64>, memref<?xf64>
  }];
  let arguments = (ins
      Arg<MemRefRankOf<[I32], [1]>, "elements to visit", [MemRead]>:$elements,
      Arg<MemRefRankOf<[F64], [1]>, "rows x cols, row-major", [MemRead]>:$matrix,
      Arg<MemRefRankOf<[F64], [1]>, "", [MemRead]>:$input,
      Arg<MemRefRankOf<[F64], [1]>, "", [MemRead, MemWrite]>:$output,
      I64Attr:$rows,
      I64Attr:$cols,
      UnitAttr:$accumulate);
  let assemblyFormat = [{
    `elements` `(` $elements `)` `matrix` `(` $matrix `)` `ins` `(` $input `)` `outs` `(` $output `)`
    attr-dict `:` type(operands)
  }];
}

#endif // MIR_OPS
