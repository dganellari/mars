// --mir-forward-owned: keep a lane's own fragment of a shared buffer in
// registers between the accesses that update it.
//
// After --mir-unroll-loops the face chain updates Y plane by plane: face f reads
// plane f, subtracts its flux and writes it back, then does the same with a plus
// sign on plane f+1 -- which face f+1 reads again. All of these accesses are
// lane-owned (tagged by --mir-chain-contracts: lane L touches only its own
// entries of the window), so between two of them only the lane itself can have
// changed its entries. Per straight-line block this pass
//   * forwards an owned read of the window and indices an earlier owned write
//     stored to: the read IS the stored fragment;
//   * forwards an owned read inside a view an earlier splat write filled, when
//     the read's window is disjoint from every owned write since: the read IS
//     the splat (the first direction's Y planes and every sweep accumulator
//     start at zero);
//   * drops an owned write that a later owned write to the same window and
//     indices overwrites before any memory read could see it.
// A Y plane then costs one read and one write per direction instead of two of
// each per face, and the first direction reads nothing.
//
// Only in-bounds reads are forwarded: a masked read returns padding in its
// out-of-bounds entries, which the stored fragment need not hold.
//
// Assumes distinct function arguments do not alias -- the destination-passing
// kernel contract: the out-parameter is written, the inputs are only read.

#include "mir/MirPasses.h"
#include "mir/Windows.h"

#include "mlir/Dialect/Arith/IR/Arith.h"
#include "mlir/Dialect/Vector/IR/VectorOps.h"
#include "mlir/IR/Builders.h"
#include "mlir/IR/Matchers.h"
#include "mlir/Interfaces/FunctionInterfaces.h"
#include "mlir/Interfaces/SideEffectInterfaces.h"
#include "mlir/Interfaces/ViewLikeInterface.h"
#include "mlir/Pass/Pass.h"

using namespace mlir;

namespace {

// The function argument a memref views, or null for per-thread memory.
static Value sharedRoot(Value mem) {
  while (Operation *d = mem.getDefiningOp()) {
    auto v = dyn_cast<ViewLikeOpInterface>(d);
    if (!v)
      return Value();
    mem = v.getViewSource();
  }
  auto ba = cast<BlockArgument>(mem);
  if (!ba.getOwner()->isEntryBlock() ||
      !isa<FunctionOpInterface>(ba.getOwner()->getParentOp()))
    return Value();
  return mem;
}

// Same memory in the same lane: same window, same index values, same shape.
static bool sameAccess(VectorTransferOpInterface a, VectorTransferOpInterface b) {
  if (mir::stripCasts(a.getSource()) != mir::stripCasts(b.getSource()) ||
      a.getVectorType() != b.getVectorType() ||
      a.getPermutationMap() != b.getPermutationMap() ||
      a.getIndices().size() != b.getIndices().size())
    return false;
  for (auto [x, y] : llvm::zip(a.getIndices(), b.getIndices()))
    if (x != y)
      return false;
  return true;
}

// A splat written over the whole of its memref (a buffer, or a view such as a
// workgroup slot's typed view): its value, or null.
static TypedAttr wholeViewSplat(vector::TransferWriteOp w) {
  auto cst = w.getVector().getDefiningOp<arith::ConstantOp>();
  auto dv = cst ? dyn_cast<DenseElementsAttr>(cst.getValue()) : nullptr;
  auto mt = dyn_cast<MemRefType>(w.getSource().getType());
  if (!dv || !dv.isSplat() || !mt || w.getMask() ||
      w.getVectorType().getShape() != mt.getShape() ||
      !w.getPermutationMap().isMinorIdentity() ||
      llvm::any_of(w.getInBoundsValues(), [](bool b) { return !b; }) ||
      llvm::any_of(w.getIndices(), [](Value i) {
        std::optional<int64_t> c = getConstantIntValue(i);
        return !c || *c != 0;
      }))
    return nullptr;
  return cast<TypedAttr>(dv.getSplatValue<Attribute>());
}

// Is `mem` the view `win`, or a subview carved out of it?
static bool within(Value mem, Value win) {
  win = mir::stripCasts(win);
  for (mem = mir::stripCasts(mem);; ) {
    if (mem == win)
      return true;
    auto sv = mem.getDefiningOp<memref::SubViewOp>();
    if (!sv)
      return false;
    mem = mir::stripCasts(sv.getSource());
  }
}

struct ForwardOwnedPass : public PassWrapper<ForwardOwnedPass, OperationPass<>> {
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(ForwardOwnedPass)

  StringRef getArgument() const final { return "mir-forward-owned"; }
  StringRef getDescription() const final {
    return "Forward lane-owned fragments of shared buffers between accesses, "
           "and drop overwritten owned writes";
  }
  void getDependentDialects(DialectRegistry &registry) const override {
    registry.insert<arith::ArithDialect>();
  }

  // An owned write whose value is still what memory holds, and whether any
  // memory read may have seen it since.
  struct Pending {
    vector::TransferWriteOp write;
    bool observed;
  };
  struct Fill {
    TypedAttr splat;
    Value window;   // the view the splat covers entirely
    // Windows owned writes stored to since the fill. Windows, not the writes:
    // a write here can be erased later as overwritten.
    SmallVector<Value> since;
  };

  // Shared roots an op reads / writes, by its declared effects (an op that
  // declares none is assumed to do both with every memref operand).
  static void effectsOf(Operation *op, SmallVectorImpl<std::pair<Value, bool>> &out) {
    auto add = [&](Value mem, bool w) {
      if (Value r = sharedRoot(mem))
        out.push_back({r, w});
    };
    auto all = [&](bool w) {
      for (Value o : op->getOperands())
        if (isa<BaseMemRefType>(o.getType()))
          add(o, w);
    };
    auto me = dyn_cast<MemoryEffectOpInterface>(op);
    if (!me) {
      all(false);
      all(true);
      return;
    }
    SmallVector<MemoryEffects::EffectInstance> effs;
    me.getEffects(effs);
    for (auto &e : effs) {
      const bool w = isa<MemoryEffects::Write>(e.getEffect());
      if (!w && !isa<MemoryEffects::Read>(e.getEffect()))
        continue;
      if (Value v = e.getValue())
        add(v, w);
      else
        all(w);
    }
  }

  void runOnBlock(Block &blk) {
    DenseMap<Value, SmallVector<Pending>> owned;
    DenseMap<Value, Fill> fills;

    for (Operation &opRef : llvm::make_early_inc_range(blk)) {
      Operation *op = &opRef;
      const bool isOwned = op->hasAttr(mir::kLaneOwnedAttr);

      if (auto w = dyn_cast<vector::TransferWriteOp>(op); w && isOwned) {
        Value root = sharedRoot(w.getSource());
        if (!root)
          continue;
        auto &list = owned[root];
        for (auto it = list.begin(); it != list.end();) {
          const bool same = sameAccess(it->write, w);
          if (same && !it->observed)
            it->write.erase();   // overwritten before anything read it
          if (same || !mir::disjointWindows(it->write.getSource(), w.getSource()))
            it = list.erase(it);   // superseded, or possibly overlapped: forget it
          else
            ++it;
        }
        list.push_back({w, false});
        if (auto f = fills.find(root); f != fills.end())
          f->second.since.push_back(w.getSource());
        continue;
      }

      if (auto r = dyn_cast<vector::TransferReadOp>(op); r && isOwned) {
        Value root = sharedRoot(r.getSource());
        if (!root)
          continue;
        const bool inBounds =
            !r.getMask() &&
            llvm::all_of(r.getInBoundsValues(), [](bool b) { return b; });
        Value fwd;
        if (inBounds)
          for (Pending &p : owned[root])
            if (sameAccess(p.write, r)) {
              fwd = p.write.getVector();
              break;
            }
        // From a fill, a masked read is exact too when the padding it returns
        // out of bounds equals the splat (a short tile's zero accumulator).
        auto padIsSplat = [&](TypedAttr splat) {
          Attribute pad;
          return !r.getMask() && matchPattern(r.getPadding(), m_Constant(&pad)) &&
                 pad == splat;
        };
        if (!fwd)
          if (auto f = fills.find(root);
              f != fills.end() && (inBounds || padIsSplat(f->second.splat)) &&
              within(r.getSource(), f->second.window) &&
              llvm::all_of(f->second.since, [&](Value win) {
                return mir::disjointWindows(win, r.getSource());
              }))
            fwd = OpBuilder(r).create<arith::ConstantOp>(
                r.getLoc(), r.getVectorType(),
                SplatElementsAttr::get(r.getVectorType(), f->second.splat));
        if (fwd) {
          r.getResult().replaceAllUsesWith(fwd);
          r.erase();
          continue;
        }
        for (Pending &p : owned[root])
          if (!mir::disjointWindows(p.write.getSource(), r.getSource()))
            p.observed = true;
        continue;
      }

      // Everything else, including region ops as a whole: a read may observe
      // any pending write of its root, a write invalidates what is known.
      SmallVector<std::pair<Value, bool>> eff;
      if (op->getNumRegions() == 0)
        effectsOf(op, eff);
      else
        op->walk([&](Operation *inner) {
          if (inner != op && inner->getNumRegions() == 0)
            effectsOf(inner, eff);
        });
      for (auto [root, write] : eff) {
        for (Pending &p : owned[root])
          p.observed = true;
        if (!write)
          continue;
        owned.erase(root);
        fills.erase(root);
        if (auto w = dyn_cast<vector::TransferWriteOp>(op))
          if (TypedAttr splat = wholeViewSplat(w))
            fills[root] = {splat, w.getSource(), {}};
      }
    }
  }

  void runOnOperation() override {
    SmallVector<Block *> blocks;
    getOperation()->walk([&](Block *b) { blocks.push_back(b); });
    for (Block *b : blocks)
      runOnBlock(*b);
  }
};

}  // namespace

namespace mir {
void registerForwardOwnedPass() { PassRegistration<ForwardOwnedPass>(); }
}  // namespace mir
