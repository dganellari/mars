// The bridge: lower the high-level mir dialect to linalg -- MLIR's on-ramp to
// GPU + tensor-core codegen (linalg -> gpu/nvgpu(mma) -> nvvm -> PTX).
//
// mir.contract %input, %op_matrix {axis} : (tensor, tensor) -> tensor
//   becomes a linalg.generic contraction:
//     out[d0..d_{R-1}] = sum_p op_matrix[d_axis, p] * input[.. p at axis ..]
// which is exactly the 1D reference-operator sweep of tensor-product
// sum-factorization, expressed as a linalg reduction that MLIR can vectorize to
// tensor cores.

#include "mir/MirDialect.h"
#include "mir/MirOps.h"
#include "mir/MirPasses.h"

#include "mlir/Dialect/Arith/IR/Arith.h"
#include "mlir/Dialect/Func/IR/FuncOps.h"
#include "mlir/Dialect/Linalg/IR/Linalg.h"
#include "mlir/Dialect/SCF/IR/SCF.h"
#include "mlir/Dialect/Tensor/IR/Tensor.h"
#include "mlir/IR/AffineMap.h"
#include "mlir/IR/Builders.h"
#include "mlir/IR/BuiltinTypes.h"
#include "mlir/IR/IRMapping.h"
#include "mlir/IR/PatternMatch.h"
#include "mlir/Pass/Pass.h"
#include "mlir/Transforms/GreedyPatternRewriteDriver.h"

using namespace mlir;

namespace {

struct ContractLowering : public OpRewritePattern<mir::ContractOp> {
  using OpRewritePattern<mir::ContractOp>::OpRewritePattern;

  LogicalResult matchAndRewrite(mir::ContractOp op,
                                PatternRewriter &rewriter) const override {
    Location loc = op.getLoc();
    MLIRContext *ctx = op.getContext();
    Value input = op.getInput();
    Value opMatrix = op.getOpMatrix();

    auto resType = dyn_cast<RankedTensorType>(op.getResult().getType());
    auto inType = dyn_cast<RankedTensorType>(input.getType());
    if (!resType || !inType)
      return failure();

    const int64_t rank = inType.getRank();
    const int64_t axis = op.getAxis();
    if (axis < 0 || axis >= rank)
      return op.emitOpError("axis out of range");

    // Axis-0 static case: emit the UNFOLDED MATMUL form instead of a rank-N
    // generic. out[i,(jk)] = sum_p D[i,p] * U[p,(jk)] -- collapse the
    // non-contracted axes into one column dim, matmul, expand back. Same math,
    // but linalg.matmul is what the tensor-core schedule (tile m8n8k4 ->
    // vectorize -> nvgpu.mma.sync) pattern-matches.
    if (axis == 0 && rank == 3 && inType.hasStaticShape() &&
        resType.hasStaticShape()) {
      SmallVector<ReassociationIndices> reassoc = {{0}, {1, 2}};
      auto elemTy = resType.getElementType();
      int64_t rows = resType.getDimSize(0);
      int64_t cols = resType.getDimSize(1) * resType.getDimSize(2);
      Value u2 = rewriter.create<tensor::CollapseShapeOp>(loc, input, reassoc);
      Value zero = rewriter.create<arith::ConstantOp>(
          loc, rewriter.getZeroAttr(elemTy));
      Value empty2 = rewriter.create<tensor::EmptyOp>(
          loc, ArrayRef<int64_t>{rows, cols}, elemTy);
      Value init2 = rewriter
                        .create<linalg::FillOp>(loc, ValueRange{zero},
                                                ValueRange{empty2})
                        .getResult(0);
      Value mm = rewriter
                     .create<linalg::MatmulOp>(loc, TypeRange{init2.getType()},
                                               ValueRange{opMatrix, u2},
                                               ValueRange{init2})
                     .getResult(0);
      rewriter.replaceOpWithNewOp<tensor::ExpandShapeOp>(op, resType, mm,
                                                         reassoc);
      return success();
    }

    const int64_t nLoops = rank + 1;  // R parallel output dims + 1 reduction (p)

    // Affine maps. Loops are (d0..d_{R-1}, p) with p = dim rank.
    AffineExpr p = getAffineDimExpr(rank, ctx);
    SmallVector<AffineExpr> outExprs, inExprs, opmExprs;
    for (int64_t i = 0; i < rank; ++i)
      outExprs.push_back(getAffineDimExpr(i, ctx));
    for (int64_t i = 0; i < rank; ++i)
      inExprs.push_back(i == axis ? p : getAffineDimExpr(i, ctx));
    opmExprs.push_back(getAffineDimExpr(axis, ctx));  // op_matrix[out_axis, p]
    opmExprs.push_back(p);

    SmallVector<AffineMap> maps{
        AffineMap::get(nLoops, 0, inExprs, ctx),
        AffineMap::get(nLoops, 0, opmExprs, ctx),
        AffineMap::get(nLoops, 0, outExprs, ctx)};

    SmallVector<utils::IteratorType> iters(rank, utils::IteratorType::parallel);
    iters.push_back(utils::IteratorType::reduction);

    // Zero-initialized output (the reduction accumulates into it).
    Value zero = rewriter.create<arith::ConstantOp>(
        loc, rewriter.getZeroAttr(resType.getElementType()));
    Value empty = rewriter.create<tensor::EmptyOp>(
        loc, resType.getShape(), resType.getElementType());
    Value init =
        rewriter.create<linalg::FillOp>(loc, ValueRange{zero}, ValueRange{empty})
            .getResult(0);

    auto generic = rewriter.create<linalg::GenericOp>(
        loc, TypeRange{resType}, ValueRange{input, opMatrix}, ValueRange{init},
        maps, iters,
        [&](OpBuilder &b, Location l, ValueRange args) {
          // args = [input_elem, op_matrix_elem, out_elem]
          Value m = b.create<arith::MulFOp>(l, args[1], args[0]);
          Value a = b.create<arith::AddFOp>(l, args[2], m);
          b.create<linalg::YieldOp>(l, a);
        });

    rewriter.replaceOp(op, generic.getResults());
    return success();
  }
};

// mir.simplex_contract -> explicit scf.for nest. The triangular bounds
// (q <= D-p, r <= D-p-q) do NOT fit linalg's hyperrectangular model, so this
// lowers to loops directly, threading the result tensor through iter_args:
//   out[p,q,k] = sum_{r=0}^{D-p-q} u[p,q,r] * table[p,q,r,k]
struct SimplexContractLowering
    : public OpRewritePattern<mir::SimplexContractOp> {
  using OpRewritePattern<mir::SimplexContractOp>::OpRewritePattern;

  LogicalResult matchAndRewrite(mir::SimplexContractOp op,
                                PatternRewriter &rewriter) const override {
    Location loc = op.getLoc();
    auto resType = dyn_cast<RankedTensorType>(op.getResult().getType());
    auto inType = dyn_cast<RankedTensorType>(op.getInput().getType());
    if (!resType || !inType || !resType.hasStaticShape())
      return failure();
    const int64_t W = inType.getDimSize(0);   // = D+1 (dense cube layout)
    const int64_t n = resType.getDimSize(2);  // quadrature points
    if (op.getDegree() + 1 != W)
      return op.emitOpError("degree+1 must equal the modal cube extent");
    Type elemTy = resType.getElementType();

    Value fzero = rewriter.create<arith::ConstantOp>(
        loc, rewriter.getZeroAttr(elemTy));
    Value empty = rewriter.create<tensor::EmptyOp>(
        loc, resType.getShape(), elemTy);
    Value init =
        rewriter.create<linalg::FillOp>(loc, ValueRange{fzero},
                                        ValueRange{empty})
            .getResult(0);

    auto idx = [&](int64_t v) {
      return rewriter.create<arith::ConstantIndexOp>(loc, v).getResult();
    };
    Value c0 = idx(0), c1 = idx(1), cW = idx(W), cn = idx(n);

    // for p in [0, W): triangular bounds derive from p and q below.
    auto pLoop = rewriter.create<scf::ForOp>(
        loc, c0, cW, c1, ValueRange{init},
        [&](OpBuilder &bp, Location lp, Value p, ValueRange op0) {
          Value ubq = bp.create<arith::SubIOp>(lp, cW, p);   // q < W - p
          auto qLoop = bp.create<scf::ForOp>(
              lp, c0, ubq, c1, op0,
              [&](OpBuilder &bq, Location lq, Value q, ValueRange oq0) {
                Value ubr = bq.create<arith::SubIOp>(lq, ubq, q);  // r < W-p-q
                auto kLoop = bq.create<scf::ForOp>(
                    lq, c0, cn, c1, oq0,
                    [&](OpBuilder &bk, Location lk, Value k, ValueRange ok0) {
                      auto rLoop = bk.create<scf::ForOp>(
                          lk, c0, ubr, c1, ValueRange{fzero},
                          [&](OpBuilder &br, Location lr, Value r,
                              ValueRange acc) {
                            Value uv = br.create<tensor::ExtractOp>(
                                lr, op.getInput(), ValueRange{p, q, r});
                            Value tv = br.create<tensor::ExtractOp>(
                                lr, op.getTable(), ValueRange{p, q, r, k});
                            Value m = br.create<arith::MulFOp>(lr, uv, tv);
                            Value s =
                                br.create<arith::AddFOp>(lr, acc[0], m);
                            br.create<scf::YieldOp>(lr, s);
                          });
                      Value updated = bk.create<tensor::InsertOp>(
                          lk, rLoop.getResult(0), ok0[0], ValueRange{p, q, k});
                      bk.create<scf::YieldOp>(lk, updated);
                    });
                bq.create<scf::YieldOp>(lq, kLoop.getResults());
              });
          bp.create<scf::YieldOp>(lp, qLoop.getResults());
        });

    rewriter.replaceOp(op, pLoop.getResults());
    return success();
  }
};

// mir.flux -> linalg.generic, but as an ELEMENTWISE MAP (no reduction): every
// quadrature point is independent, so all indexing maps are the identity and all
// iterators are parallel. The flux's own region (the authored pointwise math) is
// cloned into the linalg body.
struct FluxLowering : public OpRewritePattern<mir::FluxOp> {
  using OpRewritePattern<mir::FluxOp>::OpRewritePattern;

  LogicalResult matchAndRewrite(mir::FluxOp op,
                                PatternRewriter &rewriter) const override {
    if (op.getResults().size() != 1)
      return failure();  // single-output flux for now
    auto resType = dyn_cast<RankedTensorType>(op.getResult(0).getType());
    if (!resType)
      return failure();
    Location loc = op.getLoc();
    int64_t rank = resType.getRank();

    // Identity map for every input and the output; all loops parallel.
    AffineMap id = rewriter.getMultiDimIdentityMap(rank);
    SmallVector<AffineMap> maps(op.getInputs().size() + 1, id);
    SmallVector<utils::IteratorType> iters(rank, utils::IteratorType::parallel);

    Value init = rewriter.create<tensor::EmptyOp>(loc, resType.getShape(),
                                                  resType.getElementType());

    Block &fluxBody = op.getBody().front();
    auto generic = rewriter.create<linalg::GenericOp>(
        loc, TypeRange{resType}, op.getInputs(), ValueRange{init}, maps, iters,
        [&](OpBuilder &b, Location l, ValueRange args) {
          // Wire flux block args -> linalg per-element args, clone the pointwise
          // ops, then translate mir.yield into linalg.yield.
          IRMapping map;
          for (auto [blockArg, elem] :
               llvm::zip(fluxBody.getArguments(),
                         args.take_front(fluxBody.getNumArguments())))
            map.map(blockArg, elem);
          for (Operation &inner : fluxBody.without_terminator())
            b.clone(inner, map);
          auto yieldOp = cast<mir::YieldOp>(fluxBody.getTerminator());
          SmallVector<Value> yielded;
          for (Value v : yieldOp.getOperands())
            yielded.push_back(map.lookupOrDefault(v));
          b.create<linalg::YieldOp>(l, yielded);
        });

    rewriter.replaceOp(op, generic.getResults());
    return success();
  }
};

struct ConvertMirToLinalgPass
    : public PassWrapper<ConvertMirToLinalgPass, OperationPass<func::FuncOp>> {
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(ConvertMirToLinalgPass)

  StringRef getArgument() const final { return "convert-mir-to-linalg"; }
  StringRef getDescription() const final {
    return "Lower mir.contract to a linalg.generic contraction";
  }
  void getDependentDialects(DialectRegistry &registry) const override {
    registry.insert<linalg::LinalgDialect, tensor::TensorDialect, scf::SCFDialect,
                    arith::ArithDialect>();
  }
  void runOnOperation() override {
    RewritePatternSet patterns(&getContext());
    patterns.add<ContractLowering, FluxLowering, SimplexContractLowering>(&getContext());
    if (failed(applyPatternsAndFoldGreedily(getOperation(),
                                            std::move(patterns))))
      signalPassFailure();
  }
};

}  // namespace

namespace mir {
void registerMirPasses() { PassRegistration<ConvertMirToLinalgPass>(); }
}  // namespace mir
