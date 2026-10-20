// --mir-fuse-phases: merge consecutive phases of an element body whose
// dependences stay with one lane, then keep the values they hand over in
// registers.
//
// Phases are the scf.parallel {mir.lanes} loops of a {mir.elements} body; between
// two of them --mir-warp-per-element puts a barrier when they share a buffer. But
// when both phases have the same iteration space, point q of each runs on the
// same lane (the lane-cyclic split depends only on the space), and if every
// conflicting pair of accesses (one of them a write, same buffer) touches the same
// point q in both, no other lane is involved: the phases can run as one, point by
// point, with no barrier. A pointwise flux after the contractions it reads, or a
// plane update after the contraction it adds, is such a pair.
//
// "Same point" is decided syntactically: the same memref value, and indices that
// are the same induction variables of the two loops, position by position. Views
// of one buffer through different values are treated as conflicting.
//
// A phase may also join an earlier phase it can merge with when no phase in
// between conflicts with it (shares a buffer that one of them writes): its work
// then simply runs earlier, past phases it has nothing to do with. This groups the
// same step of different faces into one phase, as a hand kernel's thread loops
// over the faces.
//
// In a merged phase, reduction loops over the same range (two contractions of the
// same input) become one loop when neither stores to memory, so their loads of the
// shared input are made once.
//
// After merging, a load of what the same phase just stored at the same indices
// takes the stored value, and scratch buffers nothing loads any more are dropped
// (with their stores), which frees shared memory and its traffic.
//
// Values defined between phases (views, constants, scratch allocations) are first
// hoisted to the top of the element body, in order, so merged code still sees
// them; none of them depends on a phase.

#include "mir/MirPasses.h"

#include "mlir/Dialect/Func/IR/FuncOps.h"
#include "mlir/Dialect/MemRef/IR/MemRef.h"
#include "mlir/Dialect/SCF/IR/SCF.h"
#include "mlir/Dialect/Utils/StaticValueUtils.h"
#include "mlir/IR/Builders.h"
#include "mlir/IR/BuiltinOps.h"
#include "mlir/Interfaces/SideEffectInterfaces.h"
#include "mlir/Interfaces/ViewLikeInterface.h"
#include "mlir/Pass/Pass.h"

using namespace mlir;

namespace {

Value rootOf(Value m) {
  while (Operation *def = m.getDefiningOp()) {
    auto view = dyn_cast<ViewLikeOpInterface>(def);
    if (!view)
      break;
    m = view.getViewSource();
  }
  return m;
}

bool isPhase(Operation *op) {
  auto p = dyn_cast<scf::ParallelOp>(op);
  return p && p->hasAttr("mir.lanes") && p.getInitVals().empty();
}

struct Access {
  Operation *op;
  Value memref;
  SmallVector<Value> indices;
  bool write;
};

void accesses(scf::ParallelOp phase, SmallVector<Access> &out) {
  phase.walk([&](Operation *op) {
    if (auto ld = dyn_cast<memref::LoadOp>(op))
      out.push_back({op, ld.getMemRef(), SmallVector<Value>(ld.getIndices()), false});
    else if (auto st = dyn_cast<memref::StoreOp>(op))
      out.push_back({op, st.getMemRef(), SmallVector<Value>(st.getIndices()), true});
  });
}

bool sameSpace(scf::ParallelOp a, scf::ParallelOp b) {
  if (a.getNumLoops() != b.getNumLoops())
    return false;
  for (unsigned d = 0; d < a.getNumLoops(); ++d) {
    auto eq = [](Value x, Value y) {
      auto cx = getConstantIntValue(x), cy = getConstantIntValue(y);
      return cx && cy && *cx == *cy;
    };
    if (!eq(a.getLowerBound()[d], b.getLowerBound()[d]) || !eq(a.getUpperBound()[d], b.getUpperBound()[d]) ||
        !eq(a.getStep()[d], b.getStep()[d]))
      return false;
  }
  return true;
}

// Position of v among the loop's induction variables, or -1.
int ivPosition(scf::ParallelOp p, Value v) {
  for (unsigned k = 0; k < p.getNumLoops(); ++k)
    if (p.getInductionVars()[k] == v)
      return (int)k;
  return -1;
}

bool samePoint(const Access &x, scf::ParallelOp a, const Access &y, scf::ParallelOp b) {
  if (x.memref != y.memref || x.indices.size() != y.indices.size())
    return false;
  for (size_t j = 0; j < x.indices.size(); ++j) {
    const int k = ivPosition(a, x.indices[j]);
    if (k < 0 || k != ivPosition(b, y.indices[j]))
      return false;
  }
  return true;
}

// Whether the phases share a buffer that one of them writes.
bool conflict(scf::ParallelOp a, scf::ParallelOp b) {
  SmallVector<Access> xa, xb;
  accesses(a, xa);
  accesses(b, xb);
  for (const Access &x : xa)
    for (const Access &y : xb)
      if ((x.write || y.write) && rootOf(x.memref) == rootOf(y.memref))
        return true;
  return false;
}

bool canMerge(scf::ParallelOp a, scf::ParallelOp b) {
  if (!sameSpace(a, b))
    return false;
  SmallVector<Access> xa, xb;
  accesses(a, xa);
  accesses(b, xb);
  for (const Access &x : xa)
    for (const Access &y : xb)
      if ((x.write || y.write) && rootOf(x.memref) == rootOf(y.memref) && !samePoint(x, a, y, b))
        return false;
  return true;
}

// b's body runs after a's, point by point.
void merge(scf::ParallelOp a, scf::ParallelOp b) {
  Block *src = b.getBody(), *dst = a.getBody();
  for (unsigned k = 0; k < b.getNumLoops(); ++k) src->getArgument(k).replaceAllUsesWith(dst->getArgument(k));
  dst->getOperations().splice(dst->getTerminator()->getIterator(), src->getOperations(), src->begin(),
                              std::prev(src->end()));
  b.erase();
}

// Within a phase body (top level): a load of what was just stored at the same
// indices takes the stored value.
void forward(scf::ParallelOp phase) {
  struct Stored {
    Value memref;
    SmallVector<Value> indices;
    Value value;
  };
  SmallVector<Stored> live;
  auto kill = [&](Value root) {
    llvm::erase_if(live, [&](const Stored &s) { return rootOf(s.memref) == root; });
  };
  for (Operation &op : llvm::make_early_inc_range(*phase.getBody())) {
    if (auto st = dyn_cast<memref::StoreOp>(op)) {
      kill(rootOf(st.getMemRef()));
      live.push_back({st.getMemRef(), SmallVector<Value>(st.getIndices()), st.getValueToStore()});
    } else if (auto ld = dyn_cast<memref::LoadOp>(op)) {
      for (const Stored &s : live)
        if (s.memref == ld.getMemRef() && llvm::equal(s.indices, ld.getIndices())) {
          ld.getResult().replaceAllUsesWith(s.value);
          ld->erase();
          break;
        }
    } else if (op.getNumRegions() > 0) {   // a nested loop that may write
      op.walk([&](memref::StoreOp st) { kill(rootOf(st.getMemRef())); });
    }
  }
}

bool storesInside(scf::ForOp loop) {
  return loop.walk([](memref::StoreOp) { return WalkResult::interrupt(); }).wasInterrupted();
}

bool sameRange(scf::ForOp a, scf::ForOp b) {
  auto eq = [](Value x, Value y) {
    auto cx = getConstantIntValue(x), cy = getConstantIntValue(y);
    return x == y || (cx && cy && *cx == *cy);
  };
  return eq(a.getLowerBound(), b.getLowerBound()) && eq(a.getUpperBound(), b.getUpperBound()) &&
         eq(a.getStep(), b.getStep());
}

// Reduction loops of a phase body (top level) over the same range: a later one
// joins an earlier one when neither stores to memory, nothing in between defines a
// value it uses, and it reads no buffer stored to in between.
void fuseReductions(scf::ParallelOp phase) {
  Block *body = phase.getBody();
  bool changed = true;
  while (changed) {
    changed = false;
    SmallVector<scf::ForOp> loops;
    for (Operation &op : *body)
      if (auto f = dyn_cast<scf::ForOp>(op); f && !storesInside(f))
        loops.push_back(f);
    for (size_t j = 1; j < loops.size() && !changed; ++j)
      for (size_t i = 0; i < j && !changed; ++i) {
        scf::ForOp a = loops[i], b = loops[j];
        if (!sameRange(a, b))
          continue;
        // What b uses from outside itself must be defined before a; b's loads must not
        // see a store between a and b.
        DenseSet<Operation *> between;
        for (Operation *op = a->getNextNode(); op != b.getOperation(); op = op->getNextNode()) between.insert(op);
        bool ok = true;
        b.walk([&](Operation *op) {
          for (Value v : op->getOperands())
            if (Operation *def = v.getDefiningOp(); def && (between.count(def) || def == a.getOperation()))
              ok = false;
        });
        DenseSet<Value> storedBetween;
        for (Operation *op : between)
          op->walk([&](memref::StoreOp st) { storedBetween.insert(rootOf(st.getMemRef())); });
        b.walk([&](memref::LoadOp ld) {
          if (storedBetween.count(rootOf(ld.getMemRef())))
            ok = false;
        });
        if (!ok)
          continue;
        // One loop: a's iter_args then b's; a's body then b's.
        OpBuilder ob(a);
        SmallVector<Value> inits(a.getInitArgs());
        llvm::append_range(inits, b.getInitArgs());
        auto fused = ob.create<scf::ForOp>(a.getLoc(), a.getLowerBound(), a.getUpperBound(), a.getStep(), inits);
        Block *dst = fused.getBody();
        if (!dst->empty())
          dst->getTerminator()->erase();
        const unsigned na = a.getNumRegionIterArgs();
        a.getInductionVar().replaceAllUsesWith(fused.getInductionVar());
        b.getInductionVar().replaceAllUsesWith(fused.getInductionVar());
        for (unsigned k = 0; k < na; ++k) a.getRegionIterArgs()[k].replaceAllUsesWith(fused.getRegionIterArgs()[k]);
        for (unsigned k = 0; k < b.getNumRegionIterArgs(); ++k)
          b.getRegionIterArgs()[k].replaceAllUsesWith(fused.getRegionIterArgs()[na + k]);
        auto ya = cast<scf::YieldOp>(a.getBody()->getTerminator());
        auto yb = cast<scf::YieldOp>(b.getBody()->getTerminator());
        SmallVector<Value> yields(ya.getOperands());
        llvm::append_range(yields, yb.getOperands());
        dst->getOperations().splice(dst->end(), a.getBody()->getOperations(), a.getBody()->begin(),
                                    std::prev(a.getBody()->end()));
        dst->getOperations().splice(dst->end(), b.getBody()->getOperations(), b.getBody()->begin(),
                                    std::prev(b.getBody()->end()));
        OpBuilder::atBlockEnd(dst).create<scf::YieldOp>(a.getLoc(), yields);
        for (unsigned k = 0; k < na; ++k) a.getResult(k).replaceAllUsesWith(fused.getResult(k));
        for (unsigned k = 0; k < b.getNumResults(); ++k) b.getResult(k).replaceAllUsesWith(fused.getResult(na + k));
        a.erase();
        b.erase();
        changed = true;
      }
  }
}

void fuse(scf::ParallelOp elements) {
  Block *body = elements.getBody();
  // Hoist what sits between phases (it does not depend on them).
  SmallVector<Operation *> hoist;
  bool seenPhase = false;
  for (Operation &op : *body) {
    if (isPhase(&op)) {
      seenPhase = true;
      continue;
    }
    if (seenPhase && &op != body->getTerminator() && op.getNumRegions() == 0 && isMemoryEffectFree(&op))
      hoist.push_back(&op);
    else if (seenPhase && isa<memref::AllocaOp>(op))
      hoist.push_back(&op);
  }
  Operation *firstPhase = nullptr;
  for (Operation &op : *body)
    if (isPhase(&op)) {
      firstPhase = &op;
      break;
    }
  for (Operation *op : hoist) op->moveBefore(firstPhase);

  // Each phase, in order, joins the earliest phase it can merge with that it can
  // reach: it moves up past phases it does not conflict with, and stops at the
  // first one it conflicts with (which it may still join, when that conflict is
  // lane-local; moving past it would reorder same-point accesses). Joining the
  // earliest groups the same step of the different faces.
  SmallVector<scf::ParallelOp> done;   // the phases so far, in order
  for (Operation &op : llvm::make_early_inc_range(*body)) {
    if (!isPhase(&op)) {
      if (!isMemoryEffectFree(&op) && op.getNumRegions() == 0)
        done.clear();   // nothing moves past this
      continue;
    }
    auto cur = cast<scf::ParallelOp>(op);
    scf::ParallelOp into;
    for (int i = (int)done.size() - 1; i >= 0; --i) {
      const bool mergeable = canMerge(done[i], cur);
      if (mergeable)
        into = done[i];
      if (conflict(done[i], cur))
        break;
    }
    if (into)
      merge(into, cur);
    else
      done.push_back(cur);
  }
  for (Operation &op : *body)
    if (isPhase(&op)) {
      fuseReductions(cast<scf::ParallelOp>(&op));
      forward(cast<scf::ParallelOp>(&op));
    }

  // Scratch that nothing loads any more: drop its stores and the buffer.
  SmallVector<memref::AllocaOp> scratch;
  body->walk([&](memref::AllocaOp a) {
    if (a->hasAttr("mir.scratch"))
      scratch.push_back(a);
  });
  for (memref::AllocaOp a : scratch) {
    SmallVector<Operation *> stores, views;
    bool loaded = false;
    SmallVector<Value> work{a.getResult()};
    while (!work.empty()) {
      Value v = work.pop_back_val();
      for (Operation *u : v.getUsers()) {
        if (auto view = dyn_cast<ViewLikeOpInterface>(u); view && view.getViewSource() == v) {
          for (Value r : u->getResults()) work.push_back(r);
          views.push_back(u);
        } else if (isa<memref::StoreOp>(u)) {
          stores.push_back(u);
        } else {
          loaded = true;
        }
      }
    }
    if (loaded)
      continue;
    for (Operation *st : stores) st->erase();
    for (Operation *v : llvm::reverse(views))   // innermost views first
      if (v->use_empty())
        v->erase();
    if (a.use_empty())
      a->erase();
  }
}

struct FusePhasesPass : public PassWrapper<FusePhasesPass, OperationPass<ModuleOp>> {
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(FusePhasesPass)

  StringRef getArgument() const final { return "mir-fuse-phases"; }
  StringRef getDescription() const final {
    return "Merge consecutive element phases whose dependences stay on one lane, forward "
           "the values they hand over, drop scratch nothing reads";
  }

  void runOnOperation() override {
    SmallVector<scf::ParallelOp> loops;   // first collect: fusing erases phases
    getOperation()->walk([&](scf::ParallelOp p) {
      if (p->hasAttr("mir.elements"))
        loops.push_back(p);
    });
    for (scf::ParallelOp p : loops) fuse(p);
  }
};

}  // namespace

namespace mir {
void registerFusePhasesPass() { PassRegistration<FusePhasesPass>(); }
}  // namespace mir
