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
// A fill that owned writes overwrite entirely before any memory read is
// erased: the sweep tiles cover each scratch slot and the first direction's
// planes cover Y, so their zero fills -- and the barriers ordering them -- go.
// An owned write's footprint is read off its indices: the chain pass builds a
// C-fragment index as constant + lane/4 (row) and constant + 2*(lane%4)
// (column), so the 32 lanes together write an 8x8 box at the constants.
//
// Only in-bounds reads are forwarded: a masked read returns padding in its
// out-of-bounds entries, which the stored fragment need not hold.
//
// Assumes distinct function arguments do not alias -- the destination-passing
// kernel contract: the out-parameter is written, the inputs are only read.

#include "mir/MirPasses.h"
#include "mir/Windows.h"

#include "mlir/Dialect/Arith/IR/Arith.h"
#include "mlir/Dialect/GPU/IR/GPUDialect.h"
#include "mlir/Dialect/Utils/StaticValueUtils.h"
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

// The lane id, possibly through an index cast.
static bool isLane(Value v) {
  if (auto c = v.getDefiningOp<arith::IndexCastOp>())
    v = c.getIn();
  auto t = v.getDefiningOp<gpu::ThreadIdOp>();
  return t && t.getDimension() == gpu::Dimension::x;
}
static bool isConst(Value v, int64_t c) {
  std::optional<int64_t> x = getConstantIntValue(v);
  return x && *x == c;
}
// lane / 4: the fragment row.
static bool isLaneRow(Value v) {
  if (auto d = v.getDefiningOp<arith::DivUIOp>())
    return isLane(d.getLhs()) && isConst(d.getRhs(), 4);
  if (auto d = v.getDefiningOp<arith::ShRUIOp>())
    return isLane(d.getLhs()) && isConst(d.getRhs(), 2);
  return false;
}
// 2 * (lane % 4): the fragment's first column.
static bool isLaneCol(Value v) {
  auto laneMod4 = [](Value k) {
    if (auto r = k.getDefiningOp<arith::RemUIOp>())
      return isLane(r.getLhs()) && isConst(r.getRhs(), 4);
    if (auto a = k.getDefiningOp<arith::AndIOp>())
      return isLane(a.getLhs()) && isConst(a.getRhs(), 3);
    return false;
  };
  if (auto m = v.getDefiningOp<arith::MulIOp>())
    return (laneMod4(m.getLhs()) && isConst(m.getRhs(), 2)) ||
           (laneMod4(m.getRhs()) && isConst(m.getLhs(), 2));
  if (auto sh = v.getDefiningOp<arith::ShLIOp>())
    return laneMod4(sh.getLhs()) && isConst(sh.getRhs(), 1);
  return false;
}
// The constant in `lanePart + constant`, or nullopt.
static std::optional<int64_t> laneBase(Value idx, bool row) {
  auto isPart = [&](Value v) { return row ? isLaneRow(v) : isLaneCol(v); };
  if (isPart(idx))
    return 0;
  if (auto add = idx.getDefiningOp<arith::AddIOp>())
    for (int s = 0; s < 2; ++s)
      if (isPart(add->getOperand(s)))
        if (std::optional<int64_t> c = getConstantIntValue(add->getOperand(1 - s)))
          return *c;
  return std::nullopt;
}

// The box of `w`'s memref the 32 lanes of a C-fragment write cover together:
// the leading indices fixed, 8 rows and 8 columns at the index constants,
// clipped to the memref (a masked short tile writes only what exists).
static bool ownedFootprint(vector::TransferWriteOp w, SmallVectorImpl<int64_t> &lo,
                           SmallVectorImpl<int64_t> &hi) {
  auto mt = dyn_cast<MemRefType>(w.getSource().getType());
  auto vt = w.getVectorType();
  if (!mt || !mt.hasStaticShape() || w.getMask() || vt.getRank() != 2 ||
      vt.getDimSize(0) != 1 || vt.getDimSize(1) != 2 ||
      !w.getPermutationMap().isMinorIdentity())
    return false;
  const int64_t R = mt.getRank();
  ValueRange idx = w.getIndices();
  lo.assign(R, 0);
  hi.assign(R, 0);
  for (int64_t d = 0; d + 2 < R; ++d) {
    std::optional<int64_t> c = getConstantIntValue(idx[d]);
    if (!c)
      return false;
    lo[d] = *c;
    hi[d] = *c + 1;
  }
  std::optional<int64_t> r0 = laneBase(idx[R - 2], true);
  std::optional<int64_t> c0 = laneBase(idx[R - 1], false);
  if (!r0 || !c0)
    return false;
  lo[R - 2] = *r0;
  hi[R - 2] = std::min(*r0 + 8, mt.getDimSize(R - 2));
  lo[R - 1] = *c0;
  hi[R - 1] = std::min(*c0 + 8, mt.getDimSize(R - 1));
  return lo[R - 2] < hi[R - 2] && lo[R - 1] < hi[R - 1];
}

// Map a box in `mem` coordinates into `win` coordinates through static,
// unit-stride subviews (rank-reducing ones included). False if `mem` is not
// carved out of `win` that way.
static bool boxInWindow(Value mem, Value win, SmallVector<int64_t> &lo,
                        SmallVector<int64_t> &hi) {
  win = mir::stripCasts(win);
  for (mem = mir::stripCasts(mem); mem != win;) {
    auto sv = mem.getDefiningOp<memref::SubViewOp>();
    if (!sv)
      return false;
    ArrayRef<int64_t> offs = sv.getStaticOffsets(), strides = sv.getStaticStrides();
    llvm::SmallBitVector dropped = sv.getDroppedDims();
    SmallVector<int64_t> nlo, nhi;
    unsigned k = 0;
    for (unsigned d = 0; d < offs.size(); ++d) {
      if (ShapedType::isDynamic(offs[d]) || strides[d] != 1)
        return false;
      if (dropped.test(d)) {
        nlo.push_back(offs[d]);
        nhi.push_back(offs[d] + 1);
        continue;
      }
      nlo.push_back(offs[d] + lo[k]);
      nhi.push_back(offs[d] + hi[k]);
      ++k;
    }
    lo = std::move(nlo);
    hi = std::move(nhi);
    mem = mir::stripCasts(sv.getSource());
  }
  return true;
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
    vector::TransferWriteOp op;
    SmallVector<int64_t> shape;   // of `window`
    std::vector<bool> covered;    // elements owned writes have overwritten
    int64_t left = 0;             // elements still holding the splat
    bool observed = false;        // a memory read may have seen the splat
  };

  // Mark what an owned write overwrites; true once the whole fill is covered.
  static bool cover(Fill &f, vector::TransferWriteOp w) {
    SmallVector<int64_t> lo, hi;
    if (f.shape.empty() || !ownedFootprint(w, lo, hi) ||
        !boxInWindow(w.getSource(), f.window, lo, hi) || lo.size() != f.shape.size())
      return false;
    for (size_t d = 0; d < lo.size(); ++d) {
      lo[d] = std::max<int64_t>(lo[d], 0);
      hi[d] = std::min(hi[d], f.shape[d]);
      if (lo[d] >= hi[d])
        return false;
    }
    SmallVector<int64_t> at(lo.begin(), lo.end());
    while (true) {
      int64_t lin = 0;
      for (size_t d = 0; d < at.size(); ++d)
        lin = lin * f.shape[d] + at[d];
      if (!f.covered[lin]) {
        f.covered[lin] = true;
        --f.left;
      }
      size_t d = at.size();
      while (d > 0 && ++at[d - 1] == hi[d - 1]) {
        at[d - 1] = lo[d - 1];
        --d;
      }
      if (d == 0)
        break;
    }
    return f.left == 0;
  }

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
        if (auto f = fills.find(root); f != fills.end()) {
          f->second.since.push_back(w.getSource());
          if (cover(f->second, w) && !f->second.observed) {
            f->second.op.erase();   // every element overwritten, none read
            fills.erase(f);
          }
        }
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
        if (auto f = fills.find(root); f != fills.end())
          f->second.observed = true;
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
        if (!write) {
          if (auto f = fills.find(root); f != fills.end())
            f->second.observed = true;
          continue;
        }
        owned.erase(root);
        fills.erase(root);
        if (auto w = dyn_cast<vector::TransferWriteOp>(op))
          if (TypedAttr splat = wholeViewSplat(w)) {
            Fill f;
            f.splat = splat;
            f.window = w.getSource();
            f.op = w;
            auto mt = cast<MemRefType>(w.getSource().getType());
            if (mt.hasStaticShape()) {
              f.shape.assign(mt.getShape().begin(), mt.getShape().end());
              f.left = mt.getNumElements();
              f.covered.assign(f.left, false);
            }
            fills[root] = std::move(f);
          }
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
