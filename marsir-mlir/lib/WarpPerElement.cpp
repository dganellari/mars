// --mir-warp-per-element: map per-element parallel loops onto one warp per element.
//
//   scf.parallel (w) = (0) to (N) {mir.elements}      ->  scf.for w = block_id to N step grid_dim {
//     %s = memref.alloca {mir.scratch}                       %s = the warp's shared buffer
//     scf.parallel (i) = (a) to (b) {mir.lanes} {A}          scf.for i = a + lane to b step 32 {A}
//     scf.parallel (j) = (c) to (d) {mir.lanes} {B}          gpu.barrier
//                                                             scf.for j = c + lane to d step 32 {B}
//                                                             gpu.barrier
//                                                           }
//
// Launch contract: blockDim.x = 32, any grid.x >= 1. With grid.x = N each block
// does one element; with fewer, each block walks the elements in grid-stride
// order, so per-block setup (see --mir-workgroup-constants) is paid once per
// block rather than once per element.
//
// Reductions (scf.parallel with init values and scf.reduce): each lane folds its
// iterations in order, then an xor-shuffle tree over the 32 lanes, in a fixed
// order, so results are reproducible bit for bit; every lane ends with the same
// value. A 32-iteration OR reduction this way is a ballot. Stores in the element
// body itself (outside the phases, e.g. reduction results) are done by lane 0.
//
// Barriers: phases are assumed to pass data only through the {mir.scratch}
// buffers, so a barrier goes between two phases only when both touch scratch and
// one writes it, and at the end of the element body only if some phase touches
// it.
//
// Why this is correct: iterations of an scf.parallel are independent by
// definition, so any lane may run any of them, and a lane-cyclic split needs no
// synchronisation inside a phase. Consecutive phases of one element do depend
// on each other (one writes the scratch the next reads), and the barrier
// between them is exactly that order; the barrier at the end of the body keeps
// the next element from overwriting scratch that a slower lane still reads.
// Every lane reaches every barrier: the element loop depends only on the block,
// and the phases sit directly in the element body. The scratch becomes a {mir.workgroup} argument, which --mir-gpu-wrap
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
#include "mlir/IR/IRMapping.h"
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

  // Phases communicate only through the {mir.scratch} buffers (a phase does not
  // read global memory that an earlier phase of the same element wrote), so a
  // barrier is needed between two phases only when both touch scratch and one of
  // them writes it.
  DenseSet<Value> scratch;
  for (Operation &op : *elements.getBody())
    if (auto a = dyn_cast<memref::AllocaOp>(op); a && a->hasAttr("mir.scratch"))
      scratch.insert(a.getResult());
  enum { kNone = 0, kRead = 1, kWrite = 2 };
  auto scratchUse = [&](Operation *phase) {
    int use = kNone;
    phase->walk([&](Operation *op) {
      if (auto ld = dyn_cast<memref::LoadOp>(op); ld && scratch.contains(ld.getMemRef()))
        use |= kRead;
      if (auto st = dyn_cast<memref::StoreOp>(op); st && scratch.contains(st.getMemRef()))
        use |= kWrite;
    });
    return use;
  };

  Value lane = b.create<gpu::ThreadIdOp>(loc, gpu::Dimension::x);
  int previous = kNone;
  bool anyScratch = false;
  for (Operation &op : llvm::make_early_inc_range(*elements.getBody())) {
    auto par = dyn_cast<scf::ParallelOp>(op);
    if (!par || !par->hasAttr("mir.lanes"))
      continue;
    if (par.getNumLoops() != 1)
      return par.emitOpError("{mir.lanes} must be a one-dimensional loop");
    OpBuilder pb(par);
    Location pl = par.getLoc();
    const int use = scratchUse(par);
    if (use != kNone && previous != kNone && ((use | previous) & kWrite))
      pb.create<gpu::BarrierOp>(pl);
    if (use != kNone) {
      previous = use;
      anyScratch = true;
    }

    // Each lane runs the iterations lane, lane + 32, ...; reductions accumulate in
    // that order per lane, then an xor-shuffle tree combines the 32 lane values in
    // a fixed order (every lane ends with the same bits: combine is commutative).
    Value step = par.getStep()[0];
    Value start = pb.create<arith::AddIOp>(pl, par.getLowerBound()[0], pb.create<arith::MulIOp>(pl, lane, step));
    Value stride = pb.create<arith::MulIOp>(pl, step, pb.create<arith::ConstantIndexOp>(pl, kLanes));
    // The yield starts as the identity; with reductions it gets the accumulators below.
    auto loop = pb.create<scf::ForOp>(pl, start, par.getUpperBound()[0], stride, par.getInitVals(),
                                      [](OpBuilder &lb, Location ll, Value, ValueRange acc) {
                                        lb.create<scf::YieldOp>(ll, acc);
                                      });
    auto reduce = cast<scf::ReduceOp>(par.getBody()->getTerminator());
    auto combine = [&](OpBuilder &cb, unsigned i, Value acc, Value term) {
      Block &blk = reduce.getReductions()[i].front();
      IRMapping map;
      map.map(blk.getArgument(0), acc);
      map.map(blk.getArgument(1), term);
      for (Operation &rop : blk.without_terminator()) cb.clone(rop, map);
      return map.lookup(cast<scf::ReduceReturnOp>(blk.getTerminator()).getResult());
    };
    if (!par.getInitVals().empty()) {
      // The terms now become per-lane accumulators.
      OpBuilder yb(reduce);
      SmallVector<Value> next;
      for (unsigned i = 0; i < reduce.getNumOperands(); ++i)
        next.push_back(combine(yb, i, loop.getRegionIterArgs()[i], reduce.getOperands()[i]));
      loop.getBody()->getTerminator()->setOperands(next);
    }
    Block *src = par.getBody();
    src->getArgument(0).replaceAllUsesWith(loop.getInductionVar());
    loop.getBody()->getOperations().splice(loop.getBody()->getTerminator()->getIterator(), src->getOperations(),
                                           src->begin(), std::prev(src->end()));
    // (the scf.reduce stays behind in the parallel op, which is erased)
    SmallVector<Value> results(loop.getResults());
    for (int offset = kLanes / 2; offset >= 1; offset /= 2)
      for (unsigned i = 0; i < results.size(); ++i) {
        Value other = pb.create<gpu::ShuffleOp>(pl, results[i], offset, kLanes, gpu::ShuffleMode::XOR)
                          .getShuffleResult();
        results[i] = combine(pb, i, results[i], other);
      }
    par.replaceAllUsesWith(results);
    par.erase();
  }
  WalkResult nested = elements.walk([](scf::ParallelOp p) {
    return p->hasAttr("mir.lanes") ? WalkResult::interrupt() : WalkResult::advance();
  });
  if (nested.wasInterrupted())
    return elements.emitOpError("{mir.lanes} loops must sit directly in the element body");

  // Element-level work runs on every lane alike; its stores (e.g. the reduction
  // results) happen once, from lane 0.
  Value lane0 = b.create<arith::CmpIOp>(loc, arith::CmpIPredicate::eq, lane, b.create<arith::ConstantIndexOp>(loc, 0));
  for (Operation &op : llvm::make_early_inc_range(*elements.getBody())) {
    auto st = dyn_cast<memref::StoreOp>(op);
    if (!st || scratch.contains(st.getMemRef()))
      continue;
    OpBuilder sb(st);
    auto guard = sb.create<scf::IfOp>(st.getLoc(), lane0, /*withElseRegion=*/false);
    st->moveBefore(guard.thenBlock()->getTerminator());
  }

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

  Value start = b.create<gpu::BlockIdOp>(loc, gpu::Dimension::x);
  Value stride = b.create<gpu::GridDimOp>(loc, gpu::Dimension::x);
  auto loop = b.create<scf::ForOp>(loc, start, elements.getUpperBound()[0], stride);
  moveBody(elements, loop.getBody(), loop.getInductionVar());
  if (anyScratch)   // the next element must not overwrite scratch a slower lane still reads
    OpBuilder(loop.getBody()->getTerminator()).create<gpu::BarrierOp>(loc);
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
