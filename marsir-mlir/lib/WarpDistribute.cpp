// --mir-warp-distribute: distribute vector code inside
// vector.warp_execute_on_lane_0 regions across the 32 lanes of a warp, using
// upstream's VectorDistribution pattern library (which ships in the MLIR
// install but is not exposed by any release pass -- same situation as the
// hoisting utilities, same remedy: link it from our own pass).
//
// What distributes today: transfer_write at the region boundary (each lane
// writes its slice), then propagation through elementwise ops, constants,
// broadcasts and transfer_read (each lane computes only its slice; cross-lane
// values become gpu.shuffle). vector.contract has NO upstream propagation
// pattern -- contracts must be handled separately (fragment-level mma
// distribution), which is the next increment.
//
// Distribution map: greedily split the TRAILING dims of a vector so their
// product covers the 32 lanes (8x8 -> lanes (8,4), per-lane vector<1x2>;
// 8x8x8 -> (8,4), per-lane vector<1x2x8> ... trailing-major delinearization).

#include "mir/MirPasses.h"

#include "mlir/Dialect/Arith/IR/Arith.h"
#include "mlir/Dialect/GPU/IR/GPUDialect.h"
#include "mlir/Dialect/MemRef/IR/MemRef.h"
#include "mlir/Dialect/SCF/IR/SCF.h"
#include "mlir/Dialect/Vector/IR/VectorOps.h"
#include "mlir/Dialect/Vector/Transforms/VectorDistribution.h"
#include "mlir/IR/AffineMap.h"
#include "mlir/IR/Builders.h"
#include "mlir/IR/PatternMatch.h"
#include "mlir/Pass/Pass.h"
#include "mlir/Transforms/GreedyPatternRewriteDriver.h"

#include <numeric>

using namespace mlir;

namespace {

// Split warpSize over the leading dims of `shape` (largest dims first is not
// needed: we walk from dim 0) such that each chosen dim divides evenly.
// Returns the affine map whose results are the distributed dims.
static AffineMap distributionMap(Value v) {
  auto ty = dyn_cast<VectorType>(v.getType());
  MLIRContext *ctx = v.getContext();
  if (!ty)
    return AffineMap::get(ctx);
  int64_t lanes = 32;
  SmallVector<AffineExpr> dims;
  for (int64_t d = 0; d < ty.getRank() && lanes > 1; ++d) {
    int64_t sz = ty.getDimSize(d);
    int64_t take = std::gcd(sz, lanes);
    if (take > 1) {
      dims.push_back(getAffineDimExpr(d, ctx));
      lanes /= take;
    }
  }
  return AffineMap::get(ty.getRank(), 0, dims, ctx);
}

static Value shuffleFromIdx(Location loc, OpBuilder &b, Value val, Value srcIdx,
                            int64_t warpSz) {
  Value idx32 = b.create<arith::IndexCastOp>(loc, b.getI32Type(), srcIdx);
  Value width = b.create<arith::ConstantOp>(
      loc, b.getI32IntegerAttr((int32_t)warpSz));
  auto shuf = b.create<gpu::ShuffleOp>(loc, val, idx32, width,
                                       gpu::ShuffleMode::IDX);
  return shuf.getShuffleResult();
}

struct WarpDistributePass
    : public PassWrapper<WarpDistributePass, OperationPass<>> {
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(WarpDistributePass)

  StringRef getArgument() const final { return "mir-warp-distribute"; }
  StringRef getDescription() const final {
    return "Distribute vector.warp_execute_on_lane_0 regions across lanes "
           "(upstream VectorDistribution patterns)";
  }
  void getDependentDialects(DialectRegistry &registry) const override {
    registry.insert<gpu::GPUDialect, vector::VectorDialect,
                    memref::MemRefDialect, arith::ArithDialect,
                    scf::SCFDialect>();
  }

  void runOnOperation() override {
    MLIRContext *ctx = &getContext();
    {
      RewritePatternSet patterns(ctx);
      vector::populateDistributeTransferWriteOpPatterns(
          patterns, distributionMap, /*maxNumElementsToExtract=*/1);
      vector::populatePropagateWarpVectorDistributionPatterns(
          patterns, distributionMap, shuffleFromIdx);
      if (failed(applyPatternsAndFoldGreedily(getOperation(),
                                              std::move(patterns))))
        signalPassFailure();
    }
  }
};

}  // namespace

namespace mir {
void registerWarpDistributePass() { PassRegistration<WarpDistributePass>(); }
}  // namespace mir
