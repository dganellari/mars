// --mir-warp-wrap: wrap the body of each innermost scf.parallel (the
// per-element loop of the batched operator) in vector.warp_execute_on_lane_0,
// so --mir-warp-distribute can spread its vector ops across the 32 lanes.
// The lane id is materialized with gpu.lane_id (legal pre-outlining; becomes
// %laneid in PTX). After wrapping, the launch contract changes: blockDim.x
// becomes the 32 lanes and the element index moves to thread_y -- the payload
// emitter's mapping attributes must match (thread_y for the element loop).

#include "mir/MirPasses.h"

#include "mlir/Dialect/GPU/IR/GPUDialect.h"
#include "mlir/Dialect/SCF/IR/SCF.h"
#include "mlir/Dialect/Vector/IR/VectorOps.h"
#include "mlir/Dialect/Vector/Transforms/VectorDistribution.h"
#include "mlir/IR/Builders.h"
#include "mlir/Pass/Pass.h"

using namespace mlir;

namespace {

struct WarpWrapPass : public PassWrapper<WarpWrapPass, OperationPass<>> {
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(WarpWrapPass)

  StringRef getArgument() const final { return "mir-warp-wrap"; }
  StringRef getDescription() const final {
    return "Wrap innermost scf.parallel bodies in "
           "vector.warp_execute_on_lane_0 (32 lanes)";
  }
  void getDependentDialects(DialectRegistry &registry) const override {
    registry.insert<gpu::GPUDialect, vector::VectorDialect, scf::SCFDialect>();
  }

  void runOnOperation() override {
    SmallVector<scf::ParallelOp> innermost;
    getOperation()->walk([&](scf::ParallelOp op) {
      bool hasNested = false;
      op.getBody()->walk([&](scf::ParallelOp) { hasNested = true; });
      if (!hasNested)
        innermost.push_back(op);
    });

    for (scf::ParallelOp loop : innermost) {
      Block *body = loop.getBody();
      // Ops to move: everything except the terminator (scf.reduce).
      SmallVector<Operation *> toMove;
      for (Operation &op : body->without_terminator())
        toMove.push_back(&op);
      if (toMove.empty())
        continue;

      OpBuilder b(toMove.front());
      Location loc = loop.getLoc();
      Value lane =
          b.create<gpu::LaneIdOp>(loc, b.getIndexType(), /*upperBound=*/nullptr);
      auto warp = b.create<vector::WarpExecuteOnLane0Op>(
          loc, TypeRange{}, lane, /*warpSize=*/32);
      // The generic builder creates the region with an EMPTY block (or none):
      // ensure a block exists AND that it ends with vector.yield.
      if (warp.getWarpRegion().empty())
        warp.getWarpRegion().push_back(new Block());
      Block &wbody = warp.getWarpRegion().front();
      if (wbody.empty() ||
          !wbody.back().hasTrait<OpTrait::IsTerminator>()) {
        OpBuilder tb = OpBuilder::atBlockEnd(&wbody);
        tb.create<vector::YieldOp>(loc);
      }
      // Splice everything after the warp op (up to the loop terminator) into
      // the warp body in one shot -- per-op moves invalidate iteration order.
      Block::OpListType &srcOps = body->getOperations();
      Block::OpListType &dstOps = wbody.getOperations();
      auto first = std::next(Block::iterator(warp.getOperation()));
      dstOps.splice(Block::iterator(wbody.getTerminator()), srcOps, first,
                    Block::iterator(body->getTerminator()));
      // Hoist uniform scalar/index/view ops back out so boundary vector ops
      // (transfer_write) can distribute -- their memref operands must be
      // defined outside the region.
      vector::moveScalarUniformCode(warp);
    }
  }
};

}  // namespace

namespace mir {
void registerWarpWrapPass() { PassRegistration<WarpWrapPass>(); }
}  // namespace mir
