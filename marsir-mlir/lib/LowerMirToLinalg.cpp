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

    // A rank-3 contraction along a NON-leading axis cannot be folded into one
    // 2-D matmul without moving the field, but it IS a batch of 2-D ones: one
    // per index of the leading axis. Tag it so the tensor-core schedule tiles
    // that axis by 1 and each tile becomes a plain 2-D contraction.
    if (rank == 3 && axis != 0)
      generic->setAttr("mir.batch_contract", rewriter.getUnitAttr());

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
    auto tabType = dyn_cast<RankedTensorType>(op.getTable().getType());
    if (!resType || !inType || !tabType || !resType.hasStaticShape())
      return failure();
    const int64_t D = op.getDegree();
    const int64_t W = D + 1;
    const int64_t axis = op.getAxis();
    const bool tr = op.getTransposed();
    if (axis != 1 && axis != 2)
      return op.emitOpError("axis must be 1 (q) or 2 (r); the p sweep is "
                            "full-range and is a plain mir.contract");
    if (inType.getDimSize(0) != W)
      return op.emitOpError("degree+1 must equal the modal cube extent");
    // n = quadrature points. It is the trailing table extent in every form.
    const int64_t n = tabType.getDimSize(tabType.getRank() - 1);
    Type elemTy = resType.getElementType();

    // The four ragged stages share one nest: three output loops o0/o1/o2 over
    // the result, and one reduction. Only the bounds and the index order of the
    // input/table reads differ. o0 is always the p axis.
    //   A = axis 2 forward, B = axis 2 transposed,
    //   C = axis 1 forward, D_ = axis 1 transposed
    const bool A = (axis == 2 && !tr), B_ = (axis == 2 && tr),
               C_ = (axis == 1 && !tr), D_ = (axis == 1 && tr);
    if (tabType.getRank() != (axis == 2 ? 4 : 3))
      return op.emitOpError("table rank must be 4 for axis=2, 3 for axis=1");

    Value fzero = rewriter.create<arith::ConstantOp>(
        loc, rewriter.getZeroAttr(elemTy));
    Value empty = rewriter.create<tensor::EmptyOp>(
        loc, resType.getShape(), elemTy);
    Value init = rewriter
                     .create<linalg::FillOp>(loc, ValueRange{fzero},
                                             ValueRange{empty})
                     .getResult(0);

    auto idx = [&](int64_t v) {
      return rewriter.create<arith::ConstantIndexOp>(loc, v).getResult();
    };
    Value c0 = idx(0), c1 = idx(1), cW = idx(W), cn = idx(n);

    auto pLoop = rewriter.create<scf::ForOp>(
        loc, c0, cW, c1, ValueRange{init},
        [&](OpBuilder &bp, Location lp, Value o0, ValueRange it0) {
          Value wmp = bp.create<arith::SubIOp>(lp, cW, o0);   // W - p
          Value ub1 = (C_ ? cn : wmp);
          auto l1 = bp.create<scf::ForOp>(
              lp, c0, ub1, c1, it0,
              [&](OpBuilder &b1, Location l1loc, Value o1, ValueRange it1) {
                // W - p - o1, only meaningful where o1 is the q axis.
                Value wmpq = b1.create<arith::SubIOp>(l1loc, wmp, o1);
                Value ub2 = (B_ ? wmpq : cn);
                auto l2 = b1.create<scf::ForOp>(
                    l1loc, c0, ub2, c1, it1,
                    [&](OpBuilder &b2, Location l2loc, Value o2,
                        ValueRange it2) {
                      Value ubR = A ? wmpq : (C_ ? wmp : cn);
                      auto rLoop = b2.create<scf::ForOp>(
                          l2loc, c0, ubR, c1, ValueRange{fzero},
                          [&](OpBuilder &br, Location lr, Value red,
                              ValueRange acc) {
                            SmallVector<Value> inIdx, tabIdx;
                            if (A || B_)
                              inIdx = {o0, o1, red};
                            else
                              inIdx = {o0, red, o2};
                            if (A)
                              tabIdx = {o0, o1, red, o2};
                            else if (B_)
                              tabIdx = {o0, o1, o2, red};
                            else if (C_)
                              tabIdx = {o0, red, o1};
                            else
                              tabIdx = {o0, o1, red};
                            Value uv = br.create<tensor::ExtractOp>(
                                lr, op.getInput(), inIdx);
                            Value tv = br.create<tensor::ExtractOp>(
                                lr, op.getTable(), tabIdx);
                            Value m = br.create<arith::MulFOp>(lr, uv, tv);
                            Value sum = br.create<arith::AddFOp>(lr, acc[0], m);
                            br.create<scf::YieldOp>(lr, sum);
                          });
                      Value updated = b2.create<tensor::InsertOp>(
                          l2loc, rLoop.getResult(0), it2[0],
                          ValueRange{o0, o1, o2});
                      b2.create<scf::YieldOp>(l2loc, updated);
                    });
                b1.create<scf::YieldOp>(l1loc, l2.getResults());
              });
          bp.create<scf::YieldOp>(lp, l1.getResults());
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
