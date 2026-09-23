// --mir-hoist-invariant-reads: move a read out of a loop when neither its address
// nor what it reads can change between iterations.
//
// Upstream LICM hoists only side-effect-free ops, so a fragment read of a small
// operator matrix (D, W) inside the face loop runs 7 times per direction for the
// same value. On the GPU that is load/store-unit traffic, and the fast kernels are
// bound by exactly that unit, not by DRAM.
//
// A vector.transfer_read / memref.load is hoisted when its source and indices
// (and mask) are defined outside the loop and its source is a view of a function
// argument that NOTHING in the function writes: then every iteration reads the
// same memory holding the same values. Loops are visited innermost first, so a
// read can climb as far as it stays invariant.

#include "mir/MirPasses.h"

#include "mlir/Dialect/MemRef/IR/MemRef.h"
#include "mlir/Dialect/SCF/IR/SCF.h"
#include "mlir/Dialect/Vector/IR/VectorOps.h"
#include "mlir/Interfaces/FunctionInterfaces.h"
#include "mlir/Interfaces/SideEffectInterfaces.h"
#include "mlir/Interfaces/ViewLikeInterface.h"
#include "mlir/Pass/Pass.h"

using namespace mlir;

namespace {

// The value a memref is a view of: a block argument, or whatever op made it.
static Value rootOf(Value mem) {
  while (Operation *d = mem.getDefiningOp()) {
    auto v = dyn_cast<ViewLikeOpInterface>(d);
    if (!v)
      break;
    mem = v.getViewSource();
  }
  return mem;
}

struct HoistInvariantReadsPass
    : public PassWrapper<HoistInvariantReadsPass, OperationPass<>> {
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(HoistInvariantReadsPass)

  StringRef getArgument() const final { return "mir-hoist-invariant-reads"; }
  StringRef getDescription() const final {
    return "Hoist loop-invariant reads of never-written function arguments out "
           "of scf.for loops";
  }

  void runOnOperation() override {
    getOperation()->walk([&](FunctionOpInterface fn) {
      if (fn.isExternal())
        return;
      // Every root something in the function may write. An op that touches a
      // memref without declaring its effects counts as writing it.
      DenseSet<Value> written;
      fn->walk([&](Operation *op) {
        auto me = dyn_cast<MemoryEffectOpInterface>(op);
        if (!me) {
          if (op->getNumRegions() == 0)
            for (Value o : op->getOperands())
              if (isa<BaseMemRefType>(o.getType()))
                written.insert(rootOf(o));
          return;
        }
        SmallVector<MemoryEffects::EffectInstance> effs;
        me.getEffects(effs);
        for (auto &e : effs)
          if (isa<MemoryEffects::Write>(e.getEffect())) {
            if (Value v = e.getValue())
              written.insert(rootOf(v));
            else
              for (Value o : op->getOperands())
                if (isa<BaseMemRefType>(o.getType()))
                  written.insert(rootOf(o));
          }
      });
      auto readOnlyArg = [&](Value mem) {
        Value r = rootOf(mem);
        auto ba = dyn_cast<BlockArgument>(r);
        return ba && ba.getOwner()->isEntryBlock() &&
               isa<FunctionOpInterface>(ba.getOwner()->getParentOp()) &&
               !written.contains(r);
      };

      SmallVector<scf::ForOp> loops;
      fn->walk([&](scf::ForOp f) { loops.push_back(f); });   // post-order: inner first
      for (scf::ForOp loop : loops) {
        SmallVector<Operation *> hoist;
        for (Operation &op : *loop.getBody()) {
          Value src;
          if (auto r = dyn_cast<vector::TransferReadOp>(op))
            src = r.getSource();
          else if (auto l = dyn_cast<memref::LoadOp>(op))
            src = l.getMemRef();
          else
            continue;
          if (!isa<BaseMemRefType>(src.getType()) || !readOnlyArg(src) ||
              !llvm::all_of(op.getOperands(), [&](Value v) {
                return loop.isDefinedOutsideOfLoop(v);
              }))
            continue;
          hoist.push_back(&op);
        }
        for (Operation *op : hoist)
          op->moveBefore(loop);
      }
    });
  }
};

}  // namespace

namespace mir {
void registerHoistInvariantReadsPass() {
  PassRegistration<HoistInvariantReadsPass>();
}
}  // namespace mir
