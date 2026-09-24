// --mir-unroll-loops: fully unroll every scf.for with a small constant trip
// count, innermost first.
//
// In a warp kernel the loops left after tiling are short (7 faces, 8 column
// tiles, 2 k-slabs) and their iterations touch different planes and tiles of the
// same buffers. Unrolled, each access has constant offsets, so passes can tell
// which planes two accesses touch (--mir-forward-owned keeps a plane's fragment
// in registers between the faces that update it), and ptxas can schedule loads
// of the next face under the math of the current one.
//
// Loops whose bounds are not constant (a lane-strided fill) are left alone.

#include "mir/MirPasses.h"

#include "mlir/Dialect/SCF/IR/SCF.h"
#include "mlir/Dialect/SCF/Utils/Utils.h"
#include "mlir/Dialect/Utils/StaticValueUtils.h"
#include "mlir/IR/PatternMatch.h"
#include "mlir/Pass/Pass.h"

using namespace mlir;

namespace {

struct UnrollLoopsPass : public PassWrapper<UnrollLoopsPass, OperationPass<>> {
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(UnrollLoopsPass)

  UnrollLoopsPass() = default;
  UnrollLoopsPass(const UnrollLoopsPass &other) : PassWrapper(other) {}

  Option<unsigned> maxTrip{*this, "max-trip",
      llvm::cl::desc("unroll loops with at most this many iterations"),
      llvm::cl::init(8)};

  StringRef getArgument() const final { return "mir-unroll-loops"; }
  StringRef getDescription() const final {
    return "Fully unroll scf.for loops with a small constant trip count";
  }
  void getDependentDialects(DialectRegistry &registry) const override {
    registry.insert<scf::SCFDialect>();
  }

  void runOnOperation() override {
    SmallVector<scf::ForOp> loops;
    getOperation()->walk([&](scf::ForOp f) { loops.push_back(f); });  // inner first
    IRRewriter rewriter(&getContext());
    for (scf::ForOp f : loops) {
      std::optional<int64_t> lb = getConstantIntValue(f.getLowerBound());
      std::optional<int64_t> ub = getConstantIntValue(f.getUpperBound());
      std::optional<int64_t> step = getConstantIntValue(f.getStep());
      if (!lb || !ub || !step || *step <= 0 || *ub <= *lb)
        continue;
      const int64_t trip = (*ub - *lb + *step - 1) / *step;
      if (trip > (int64_t)maxTrip)
        continue;
      // Unrolling by the full trip count leaves one iteration, which
      // loopUnrollByFactor itself promotes (erasing the loop).
      if (trip == 1) {
        (void)f.promoteIfSingleIteration(rewriter);
      } else if (failed(loopUnrollByFactor(f, (uint64_t)trip))) {
        f.emitOpError("mir-unroll-loops: could not unroll a constant-trip loop");
        return signalPassFailure();
      }
    }
  }
};

}  // namespace

namespace mir {
void registerUnrollLoopsPass() { PassRegistration<UnrollLoopsPass>(); }
}  // namespace mir
