// --mir-lower-copies: rewrite memref.copy into plain scf.for loops. GPU
// kernels cannot call the memrefCopy host-runtime fallback that
// finalize-memref-to-llvm emits for strided layouts; explicit loops are legal
// anywhere and these copies are tiny thread-local tiles.

#include "mir/MirPasses.h"

#include "mlir/Dialect/Arith/IR/Arith.h"
#include "mlir/Dialect/MemRef/IR/MemRef.h"
#include "mlir/Dialect/SCF/IR/SCF.h"
#include "mlir/IR/Builders.h"
#include "mlir/IR/PatternMatch.h"
#include "mlir/Pass/Pass.h"
#include "mlir/Transforms/GreedyPatternRewriteDriver.h"

using namespace mlir;

namespace {

struct CopyToLoops : public OpRewritePattern<memref::CopyOp> {
  using OpRewritePattern<memref::CopyOp>::OpRewritePattern;

  LogicalResult matchAndRewrite(memref::CopyOp op,
                                PatternRewriter &rewriter) const override {
    auto ty = dyn_cast<MemRefType>(op.getSource().getType());
    if (!ty || !ty.hasStaticShape())
      return failure();
    Location loc = op.getLoc();
    int64_t rank = ty.getRank();
    Value c0 = rewriter.create<arith::ConstantIndexOp>(loc, 0);
    Value c1 = rewriter.create<arith::ConstantIndexOp>(loc, 1);
    SmallVector<Value> ivs;
    OpBuilder b = rewriter;
    for (int64_t d = 0; d < rank; ++d) {
      Value ub = b.create<arith::ConstantIndexOp>(loc, ty.getDimSize(d));
      auto loop = b.create<scf::ForOp>(loc, c0, ub, c1);
      ivs.push_back(loop.getInductionVar());
      b.setInsertionPointToStart(loop.getBody());
    }
    Value v = b.create<memref::LoadOp>(loc, op.getSource(), ivs);
    b.create<memref::StoreOp>(loc, v, op.getTarget(), ivs);
    rewriter.eraseOp(op);
    return success();
  }
};

struct LowerCopiesPass
    : public PassWrapper<LowerCopiesPass, OperationPass<>> {
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(LowerCopiesPass)
  StringRef getArgument() const final { return "mir-lower-copies"; }
  StringRef getDescription() const final {
    return "Rewrite static-shape memref.copy into scf.for loops (GPU-legal)";
  }
  void runOnOperation() override {
    RewritePatternSet patterns(&getContext());
    patterns.add<CopyToLoops>(&getContext());
    if (failed(applyPatternsAndFoldGreedily(getOperation(),
                                            std::move(patterns))))
      signalPassFailure();
  }
};

}  // namespace

namespace mir {
void registerLowerCopiesPass() { PassRegistration<LowerCopiesPass>(); }
}  // namespace mir
