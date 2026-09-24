// --mir-warp-barriers: order a warp kernel's accesses to memory the warp shares.
//
// After --mir-chain-contracts / --mir-workgroup-buffers / --mir-distribute-fills
// a kernel touches shared memory (kernel arguments and workgroup slots) in two
// ways: lane-owned C-fragment accesses (tagged mir.lane_owned: lane L touches
// only its own entries of the window), and everything else -- operand fragments
// read across lanes, distributed fills, redundant whole-value accesses. Where one
// lane may read or overwrite what another lane wrote, or overwrite what another
// lane still reads, a gpu.barrier goes between the two.
//
// Two accesses to the same root argument CONFLICT when at least one writes,
// except in three provable cases:
//   * LANE-CONSISTENT: both lane-owned, same vector type and map, index values
//     equal, and windows that map (row, col) to the same elements -- the same
//     view, or subviews of one source that differ only in a DROPPED (unit,
//     rank-reduced) dimension. Every element is then touched by the same lane in
//     both: the plane-by-plane read-modify-write of Y within one direction.
//   * DISJOINT TILES: a carried access from an earlier iteration of an scf.for
//     and a current one, both through the same subview whose offset in a kept
//     dimension IS the induction variable with size <= step. Iterations touch
//     disjoint tiles: the column tiles of a sweep.
//   * DISJOINT WINDOWS: two transfers through subviews of one buffer whose
//     static ranges do not meet -- the same tiles once the loop is unrolled.
// "Equal" for index values means equal at run time: the same SSA value (unless
// it is defined inside the loop that carries the earlier access -- a new value
// every iteration), or the same pure expression over such values.
//
// Loops are scanned with their own accesses already pending, so a body that
// conflicts with its previous iteration gets its barrier inside. A barrier
// before a loop is preferred (it runs once). A loop or branch whose control
// depends on the lane (a distributed fill) is one access as a whole: a barrier
// inside it would be divergent.

#include "mir/MirPasses.h"
#include "mir/Windows.h"

#include "mlir/Dialect/Arith/IR/Arith.h"
#include "mlir/Dialect/GPU/IR/GPUDialect.h"
#include "mlir/Dialect/MemRef/IR/MemRef.h"
#include "mlir/Dialect/SCF/IR/SCF.h"
#include "mlir/Dialect/Utils/StaticValueUtils.h"
#include "mlir/Dialect/Vector/IR/VectorOps.h"
#include "mlir/IR/Builders.h"
#include "mlir/Interfaces/FunctionInterfaces.h"
#include "mlir/Interfaces/SideEffectInterfaces.h"
#include "mlir/Interfaces/ViewLikeInterface.h"
#include "mlir/Pass/Pass.h"

using namespace mlir;

namespace {

// The function argument a memref is a view of, or null if it is not one
// (an alloc, an alloca: per-thread memory, never shared).
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

static bool dependsOnLane(Value v, int depth = 0) {
  if (depth > 16)
    return true;
  Operation *d = v.getDefiningOp();
  if (!d)
    return false;
  if (isa<gpu::ThreadIdOp, gpu::LaneIdOp>(d))
    return true;
  return llvm::any_of(d->getOperands(),
                      [&](Value o) { return dependsOnLane(o, depth + 1); });
}

// Does this region op's control flow differ between lanes?
static bool laneDivergent(Operation *op) {
  if (auto f = dyn_cast<scf::ForOp>(op))
    return dependsOnLane(f.getLowerBound()) || dependsOnLane(f.getUpperBound()) ||
           dependsOnLane(f.getStep());
  if (auto i = dyn_cast<scf::IfOp>(op))
    return dependsOnLane(i.getCondition());
  return !isa<LoopLikeOpInterface>(op);   // unknown region op: do not enter
}

// Equal at run time, for the access pair being compared. `carried` is the loop
// whose earlier iteration produced one of them: a value defined inside it is a
// different value in the current iteration.
static bool sameValue(Value a, Value b, Operation *carried, int depth = 0) {
  if (a == b) {
    if (!carried)
      return true;
    Operation *scope = a.getParentRegion()->getParentOp();
    return !carried->isAncestor(scope);
  }
  if (depth > 8)
    return false;
  std::optional<int64_t> ca = getConstantIntValue(a), cb = getConstantIntValue(b);
  if (ca && cb)
    return *ca == *cb;
  Operation *da = a.getDefiningOp(), *db = b.getDefiningOp();
  if (!da || !db || da->getName() != db->getName() ||
      da->getAttrDictionary() != db->getAttrDictionary() ||
      da->getNumOperands() != db->getNumOperands() || da->getNumResults() != 1 ||
      !isPure(da))
    return false;
  for (auto [x, y] : llvm::zip(da->getOperands(), db->getOperands()))
    if (!sameValue(x, y, carried, depth + 1))
      return false;
  return true;
}

static bool sameOffset(OpFoldResult a, OpFoldResult b, Operation *carried) {
  std::optional<int64_t> ca = getConstantIntValue(a), cb = getConstantIntValue(b);
  if (ca || cb)
    return ca && cb && *ca == *cb;
  return sameValue(cast<Value>(a), cast<Value>(b), carried);
}

// Do two windows map every (row, col) of a lane-owned access to the same element
// of the root, up to a dropped dimension?
static bool sameWindow(Value a, Value b, Operation *carried) {
  if (sameValue(a, b, carried))
    return true;
  auto sa = a.getDefiningOp<memref::SubViewOp>();
  auto sb = b.getDefiningOp<memref::SubViewOp>();
  if (!sa || !sb || sa.getType() != sb.getType() ||
      !sameWindow(sa.getSource(), sb.getSource(), carried) ||
      sa.getStaticSizes() != sb.getStaticSizes() ||
      sa.getStaticStrides() != sb.getStaticStrides() ||
      ShapedType::isDynamicShape(sa.getStaticSizes()) ||
      ShapedType::isDynamicShape(sa.getStaticStrides()))
    return false;
  llvm::SmallBitVector dropped = sa.getDroppedDims();
  if (dropped != sb.getDroppedDims())
    return false;
  auto oa = sa.getMixedOffsets(), ob = sb.getMixedOffsets();
  for (unsigned d = 0; d < oa.size(); ++d)
    if (!dropped.test(d) && !sameOffset(oa[d], ob[d], carried))
      return false;
  return true;
}

static bool laneConsistent(Operation *x, Operation *y, Operation *carried) {
  if (!x->hasAttr(mir::kLaneOwnedAttr) || !y->hasAttr(mir::kLaneOwnedAttr))
    return false;
  auto tx = dyn_cast<VectorTransferOpInterface>(x);
  auto ty = dyn_cast<VectorTransferOpInterface>(y);
  if (!tx || !ty || tx.getVectorType() != ty.getVectorType() ||
      tx.getPermutationMap() != ty.getPermutationMap() ||
      tx.getIndices().size() != ty.getIndices().size())
    return false;
  for (auto [a, b] : llvm::zip(tx.getIndices(), ty.getIndices()))
    if (!sameValue(a, b, carried))
      return false;
  return sameWindow(tx.getSource(), ty.getSource(), carried);
}

static bool disjointTiles(Operation *x, Operation *y, Operation *carried) {
  auto loop = dyn_cast_or_null<scf::ForOp>(carried);
  auto tx = dyn_cast<VectorTransferOpInterface>(x);
  auto ty = dyn_cast<VectorTransferOpInterface>(y);
  if (!loop || !tx || !ty || tx.getSource() != ty.getSource())
    return false;
  auto sv = tx.getSource().getDefiningOp<memref::SubViewOp>();
  std::optional<int64_t> step = getConstantIntValue(loop.getStep());
  if (!sv || !step)
    return false;
  llvm::SmallBitVector dropped = sv.getDroppedDims();
  auto offs = sv.getMixedOffsets();
  for (unsigned d = 0; d < offs.size(); ++d) {
    if (dropped.test(d))
      continue;
    auto off = dyn_cast<Value>(offs[d]);
    const int64_t size = sv.getStaticSizes()[d];
    if (off && off == loop.getInductionVar() && !ShapedType::isDynamic(size) &&
        size <= *step)
      return true;
  }
  return false;
}

// Two transfers through statically disjoint windows of the same buffer (the
// column tiles of an unrolled sweep) never touch the same element.
static bool disjointAccesses(Operation *x, Operation *y) {
  auto tx = dyn_cast<VectorTransferOpInterface>(x);
  auto ty = dyn_cast<VectorTransferOpInterface>(y);
  return tx && ty && mir::disjointWindows(tx.getSource(), ty.getSource());
}

struct Rec {
  Value root;
  bool write;
  Operation *op;
  Operation *carried;   // region op this came out of, or null: same scope
  bool operator==(const Rec &o) const {
    return root == o.root && write == o.write && op == o.op && carried == o.carried;
  }
};

static void accessesOf(Operation *op, SmallVectorImpl<Rec> &out) {
  auto add = [&](Value mem, bool w) {
    if (Value r = sharedRoot(mem))
      out.push_back({r, w, op, nullptr});
  };
  auto allMemrefs = [&](bool w) {
    for (Value o : op->getOperands())
      if (isa<BaseMemRefType>(o.getType()))
        add(o, w);
  };
  auto me = dyn_cast<MemoryEffectOpInterface>(op);
  if (!me) {
    allMemrefs(false);
    allMemrefs(true);
    return;
  }
  SmallVector<MemoryEffects::EffectInstance> effs;
  me.getEffects(effs);
  for (auto &e : effs) {
    const bool w = isa<MemoryEffects::Write>(e.getEffect());
    if (!w && !isa<MemoryEffects::Read>(e.getEffect()))
      continue;   // allocate / free
    if (Value v = e.getValue()) {
      if (isa<BaseMemRefType>(v.getType()))
        add(v, w);
    } else {
      allMemrefs(w);
    }
  }
}

struct WarpBarriersPass : public PassWrapper<WarpBarriersPass, OperationPass<>> {
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(WarpBarriersPass)

  StringRef getArgument() const final { return "mir-warp-barriers"; }
  StringRef getDescription() const final {
    return "Insert gpu.barrier between conflicting warp-shared memory accesses";
  }
  void getDependentDialects(DialectRegistry &registry) const override {
    registry.insert<gpu::GPUDialect>();
  }

  static bool conflicts(ArrayRef<Rec> st, const Rec &y) {
    for (const Rec &x : st) {
      if (x.root != y.root || (!x.write && !y.write))
        continue;
      if (laneConsistent(x.op, y.op, x.carried) ||
          disjointTiles(x.op, y.op, x.carried) || disjointAccesses(x.op, y.op))
        continue;
      return true;
    }
    return false;
  }

  static void merge(SmallVector<Rec> &st, ArrayRef<Rec> add) {
    for (const Rec &r : add)
      if (!llvm::is_contained(st, r))
        st.push_back(r);
  }

  // Every access inside a region op, attributed to the op itself.
  static void summarize(Operation *op, SmallVectorImpl<Rec> &out) {
    op->walk([&](Operation *inner) {
      if (inner != op && inner->getNumRegions() == 0)
        accessesOf(inner, out);
    });
    for (Rec &r : out)
      r.carried = op;
  }

  void barrierBefore(Operation *op, SmallVector<Rec> &st) {
    OpBuilder(op).create<gpu::BarrierOp>(op->getLoc());
    st.clear();
  }

  void scan(Block &blk, SmallVector<Rec> &st) {
    for (Operation &opRef : llvm::make_early_inc_range(blk)) {
      Operation *op = &opRef;
      if (isa<gpu::BarrierOp>(op)) {
        st.clear();
        continue;
      }
      SmallVector<Rec> acc;
      const bool leaf = op->getNumRegions() == 0 || laneDivergent(op);
      if (leaf)
        accessesOf(op, acc);
      if (op->getNumRegions() > 0)
        summarize(op, acc);
      if (acc.empty())
        continue;
      if (llvm::any_of(acc, [&](const Rec &y) { return conflicts(st, y); }))
        barrierBefore(op, st);
      if (leaf) {
        merge(st, acc);
        continue;
      }
      // Earlier iterations of a loop are pending when its body starts.
      const bool loop = isa<LoopLikeOpInterface>(op);
      SmallVector<Rec> out = st;
      for (Region &r : op->getRegions())
        for (Block &b : r) {
          SmallVector<Rec> in = st;
          if (loop)
            merge(in, acc);
          scan(b, in);
          for (Rec &x : in)
            if (!x.carried)
              x.carried = op;
          merge(out, in);
        }
      st = std::move(out);
    }
  }

  void runOnOperation() override {
    getOperation()->walk([&](FunctionOpInterface fn) {
      if (fn.isExternal())
        return;
      SmallVector<Rec> st;
      scan(fn.getFunctionBody().front(), st);
    });
  }
};

}  // namespace

namespace mir {
void registerWarpBarriersPass() { PassRegistration<WarpBarriersPass>(); }
}  // namespace mir
