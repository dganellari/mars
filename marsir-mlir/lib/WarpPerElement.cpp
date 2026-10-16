// --mir-warp-per-element: map per-element parallel loops onto one warp per element.
//
//   scf.parallel (w) = (0) to (N) {mir.elements}      ->  w = gpu.block_id x
//     %s = memref.alloca {mir.scratch}                     scf.if w < N {
//     scf.parallel (i) = (a) to (b) {mir.lanes} {A}          %s = the warp's shared buffer
//     scf.parallel (j) = (c) to (d) {mir.lanes} {B}          scf.for i = a + lane to b step 32 {A}
//                                                             gpu.barrier
//                                                             scf.for j = c + lane to d step 32 {B}
//                                                           }
//
// Launch contract: grid.x >= N, blockDim.x = 32, as for --mir-batch-elements.
//
// Why this is correct: iterations of an scf.parallel are independent by
// definition, so any lane may run any of them, and a lane-cyclic split needs no
// synchronisation inside a phase. Consecutive phases of one element do depend
// on each other (one writes the scratch the next reads), and the barrier
// between them is exactly that order. Every lane reaches every barrier: the
// guard depends only on the block, and the phases sit directly in the element
// body. The scratch becomes a {mir.workgroup} argument, which --mir-gpu-wrap
// turns into shared memory: one copy per warp, which the lanes share.

#include "mir/MirPasses.h"

#include "mlir/Dialect/Arith/IR/Arith.h"
#include "mlir/Dialect/Func/IR/FuncOps.h"
#include "mlir/Dialect/GPU/IR/GPUDialect.h"
#include "mlir/Dialect/MemRef/IR/MemRef.h"
#include "mlir/Dialect/SCF/IR/SCF.h"
#include "mlir/Dialect/Utils/StaticValueUtils.h"
#include "mlir/IR/Builders.h"
#include "mlir/IR/BuiltinOps.h"
#include "mlir/Pass/Pass.h"

using namespace mlir;

namespace {

constexpr int kLanes = 32;

// Moves the body of a one-dimensional scf.parallel in front of `dest`'s
// terminator, with its induction variable replaced by `iv`.
void moveBody(scf::ParallelOp par, Block *dest, Value iv) {
  Block *src = par.getBody();
  src->getArgument(0).replaceAllUsesWith(iv);
  dest->getOperations().splice(dest->getTerminator()->getIterator(), src->getOperations(), src->begin(),
                               std::prev(src->end()));   // all but the scf.reduce terminator
}

LogicalResult mapElements(func::FuncOp func, scf::ParallelOp elements) {
  if (elements.getNumLoops() != 1 || !elements.getInitVals().empty())
    return elements.emitOpError("{mir.elements} must be a one-dimensional loop without reductions");
  if (getConstantIntValue(elements.getLowerBound()[0]) != 0 || getConstantIntValue(elements.getStep()[0]) != 1)
    return elements.emitOpError("{mir.elements} must run from 0 in steps of 1");
  OpBuilder b(elements);
  Location loc = elements.getLoc();

  // Lanes first, while their positions in the element body are still easy to see.
  Value lane = b.create<gpu::ThreadIdOp>(loc, gpu::Dimension::x);
  bool first = true;
  for (Operation &op : llvm::make_early_inc_range(*elements.getBody())) {
    auto par = dyn_cast<scf::ParallelOp>(op);
    if (!par || !par->hasAttr("mir.lanes"))
      continue;
    if (par.getNumLoops() != 1 || !par.getInitVals().empty())
      return par.emitOpError("{mir.lanes} must be a one-dimensional loop without reductions");
    OpBuilder pb(par);
    Location pl = par.getLoc();
    if (!first)
      pb.create<gpu::BarrierOp>(pl);
    first = false;
    Value step = par.getStep()[0];
    Value start = pb.create<arith::AddIOp>(pl, par.getLowerBound()[0], pb.create<arith::MulIOp>(pl, lane, step));
    Value stride = pb.create<arith::MulIOp>(pl, step, pb.create<arith::ConstantIndexOp>(pl, kLanes));
    auto loop = pb.create<scf::ForOp>(pl, start, par.getUpperBound()[0], stride);
    moveBody(par, loop.getBody(), loop.getInductionVar());
    par.erase();
  }
  WalkResult nested = elements.walk([](scf::ParallelOp p) {
    return p->hasAttr("mir.lanes") ? WalkResult::interrupt() : WalkResult::advance();
  });
  if (nested.wasInterrupted())
    return elements.emitOpError("{mir.lanes} loops must sit directly in the element body");

  // The scratch: one buffer per warp, in shared memory.
  for (Operation &op : llvm::make_early_inc_range(*elements.getBody())) {
    auto alloca = dyn_cast<memref::AllocaOp>(op);
    if (!alloca || !alloca->hasAttr("mir.scratch"))
      continue;
    if (!alloca.getType().hasStaticShape())
      return alloca.emitOpError("{mir.scratch} must have a static shape");
    const unsigned at = func.getNumArguments();
    func.insertArgument(at, alloca.getType(),
                        b.getDictionaryAttr(b.getNamedAttr("mir.workgroup", b.getUnitAttr())), loc);
    alloca.replaceAllUsesWith(func.getArgument(at));
    alloca.erase();
  }

  Value w = b.create<gpu::BlockIdOp>(loc, gpu::Dimension::x);
  Value inRange = b.create<arith::CmpIOp>(loc, arith::CmpIPredicate::ult, w, elements.getUpperBound()[0]);
  auto guard = b.create<scf::IfOp>(loc, inRange, /*withElseRegion=*/false);
  moveBody(elements, guard.thenBlock(), w);
  elements.erase();
  return success();
}

struct WarpPerElementPass : public PassWrapper<WarpPerElementPass, OperationPass<ModuleOp>> {
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(WarpPerElementPass)

  StringRef getArgument() const final { return "mir-warp-per-element"; }
  StringRef getDescription() const final {
    return "Map a {mir.elements} scf.parallel onto gpu.block_id x and its "
           "{mir.lanes} phases onto the 32 lanes, with a barrier between phases";
  }
  void getDependentDialects(DialectRegistry &registry) const override {
    registry.insert<arith::ArithDialect, gpu::GPUDialect, memref::MemRefDialect, scf::SCFDialect>();
  }

  void runOnOperation() override {
    for (auto func : getOperation().getOps<func::FuncOp>()) {
      SmallVector<scf::ParallelOp> loops;
      func.walk([&](scf::ParallelOp p) {
        if (p->hasAttr("mir.elements"))
          loops.push_back(p);
      });
      if (loops.size() > 1) {
        func.emitOpError("expects at most one {mir.elements} loop");
        return signalPassFailure();
      }
      for (scf::ParallelOp p : loops)
        if (failed(mapElements(func, p)))
          return signalPassFailure();
    }
  }
};

}  // namespace

namespace mir {
void registerWarpPerElementPass() { PassRegistration<WarpPerElementPass>(); }
}  // namespace mir
