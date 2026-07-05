// --mir-hoist-transfer-pairs: promote a loop-invariant vector.transfer_read /
// vector.transfer_write pair on the same memref view to a loop-carried
// iter_arg, so the accumulator of a tiled reduction (e.g. the C tile of the
// k-loop in the DMMA schedule) lives in registers instead of round-tripping
// global memory every iteration.
//
// WHY THIS EXISTS: upstream linalg::hoistRedundantVectorTransfers declines the
// pattern whenever other transfers (the A/B operand reads) share the loop --
// it cannot prove distinct function-arg memrefs don't alias. This is a DOMAIN
// pass: mir kernels take D/U/Y as distinct buffers by construction, so we key
// the safety check on the SSA view (the pair must be the view's only in-loop
// users) rather than global alias analysis. Do not reuse on IR where the
// operand buffers may alias the accumulator.

#include "mir/MirPasses.h"

#include "mlir/Dialect/Arith/IR/Arith.h"
#include "mlir/Dialect/SCF/IR/SCF.h"
#include "mlir/Dialect/Vector/IR/VectorOps.h"
#include "mlir/IR/Builders.h"
#include "mlir/IR/PatternMatch.h"
#include "mlir/Pass/Pass.h"

using namespace mlir;

namespace {

struct HoistTransferPairsPass
    : public PassWrapper<HoistTransferPairsPass, OperationPass<>> {
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(HoistTransferPairsPass)

  HoistTransferPairsPass() = default;
  HoistTransferPairsPass(const HoistTransferPairsPass &other)
      : PassWrapper(other) {}

  StringRef getArgument() const final { return "mir-hoist-transfer-pairs"; }
  StringRef getDescription() const final {
    return "Promote loop-invariant transfer_read/write pairs to iter_args "
           "(accumulator registerization for tiled reductions)";
  }

  // The mir apply contract pre-zeroes the output, so the accumulator's initial
  // read would only load zeros: start from a zero splat instead and save the
  // read traffic. Changes semantics from Y += A*B to Y = A*B -- only enable
  // when the surrounding operator guarantees zero-initialized output.
  Option<bool> assumeZeroInit{
      *this, "assume-zero-init",
      llvm::cl::desc("replace the hoisted initial read with a zero vector "
                     "(output is guaranteed pre-zeroed)"),
      llvm::cl::init(false)};

  void runOnOperation() override {
    // Repeat until no more pairs hoist (a pair may become hoistable after an
    // inner loop was rewritten).
    bool changed = true;
    while (changed) {
      changed = false;
      getOperation()->walk([&](scf::ForOp loop) {
        if (succeeded(hoistOnePair(loop))) {
          changed = true;
          return WalkResult::interrupt();
        }
        return WalkResult::advance();
      });
    }
  }

  LogicalResult hoistOnePair(scf::ForOp loop) {
    vector::TransferReadOp read;
    vector::TransferWriteOp write;

    // Find a read/write pair on the same source with loop-invariant operands.
    for (Operation &op : loop.getBody()->getOperations()) {
      auto r = dyn_cast<vector::TransferReadOp>(op);
      if (!r)
        continue;
      Value src = r.getSource();
      if (!loop.isDefinedOutsideOfLoop(src))
        continue;
      if (llvm::any_of(r.getIndices(), [&](Value v) {
            return !loop.isDefinedOutsideOfLoop(v);
          }))
        continue;
      // The pair must be the ONLY in-loop users of this view (SSA identity --
      // see header note on why this is safe for mir kernels).
      vector::TransferWriteOp w;
      bool clean = true;
      for (Operation *user : src.getUsers()) {
        if (!loop->isAncestor(user))
          continue;
        if (user == r.getOperation())
          continue;
        if (auto cw = dyn_cast<vector::TransferWriteOp>(user)) {
          if (w || cw.getSource() != src ||
              cw.getIndices() != r.getIndices() ||
              cw.getVectorType() != r.getVectorType()) {
            clean = false;
            break;
          }
          w = cw;
          continue;
        }
        clean = false;
        break;
      }
      if (!clean || !w)
        continue;
      // The read must happen before the write in the block.
      if (!r->isBeforeInBlock(w))
        continue;
      read = r;
      write = w;
      break;
    }
    if (!read)
      return failure();

    // Capture everything needed AFTER the rewrite before any op is erased.
    Value src = read.getSource();
    SmallVector<Value> indices(read.getIndices());
    Location wloc = write.getLoc();
    Value storedVal = write.getVector();

    // 1. Initial accumulator: re-read the tile before the loop, or a zero
    //    splat when the output is guaranteed pre-zeroed (assume-zero-init).
    IRRewriter rewriter(loop.getContext());
    rewriter.setInsertionPoint(loop);
    Value init;
    if (assumeZeroInit) {
      init = rewriter.create<arith::ConstantOp>(
          read.getLoc(), rewriter.getZeroAttr(read.getVectorType()));
    } else {
      init = cast<vector::TransferReadOp>(rewriter.clone(*read.getOperation()))
                 .getResult();
    }

    // 2. Add an iter_arg initialized with it; yield what the write stored.
    // replaceWithAdditionalYields MOVES the body into a new ForOp, so the
    // original read/write handles stay valid inside the new loop.
    auto newLoopOr = loop.replaceWithAdditionalYields(
        rewriter, ValueRange{init},
        /*replaceInitOperandUsesInLoop=*/false,
        [&](OpBuilder &, Location, ArrayRef<BlockArgument>) {
          return SmallVector<Value>{storedVal};
        });
    if (failed(newLoopOr))
      return failure();
    auto newLoop = cast<scf::ForOp>(*newLoopOr);

    // 3. Inside the loop: the read becomes the new iter_arg; the write is gone.
    BlockArgument acc = newLoop.getRegionIterArgs().back();
    rewriter.replaceAllUsesWith(read.getResult(), acc);
    rewriter.eraseOp(read);
    rewriter.eraseOp(write);

    // 4. Write the final accumulator back AFTER the loop.
    rewriter.setInsertionPointAfter(newLoop);
    rewriter.create<vector::TransferWriteOp>(
        wloc, newLoop.getResults().back(), src, indices);

    return success();
  }
};

}  // namespace

namespace mir {
void registerHoistTransferPairsPass() {
  PassRegistration<HoistTransferPairsPass>();
}
}  // namespace mir
