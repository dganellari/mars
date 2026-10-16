// --mir-warp-per-element: map per-element parallel loops onto the lanes of a block.
//
//   scf.parallel (w) = (0) to (N) {mir.elements}      ->  scf.for c = block_id to N / G step grid_dim {
//     %s = memref.alloca {mir.scratch}                       w = c G + group, %s = the group's shared buffer
//     scf.parallel (i) = (a) to (b) {mir.lanes} {A}          scf.for i = a + glane to b step L {A}
//     scf.parallel (j) = (c) to (d) {mir.lanes} {B}          gpu.barrier
//                                                             scf.for j = c + glane to d step L {B}
//                                                             gpu.barrier
//                                                           }
//
// Shape: L lanes per element (a power of two, at most 32) and W warps per block,
// so G = 32 W / L elements per block; group = thread / L, glane = thread % L. L
// comes from the element loop's {mir.lanes_per_element} (default 32), W from
// {mir.warps_per_block} (default 1 for L = 32, else 4). With L = 32 and W = 1 a
// block is one warp on one element at a time. The element loop of the last
// block may run past N: those groups compute on a valid element and store
// nothing. The function gets {mir.block_threads} and {mir.elements_per_block}
// for the launcher (--mir-gpu-wrap publishes them). Launch: blockDim.x = 32 W,
// any grid.x >= 1 (blocks stride over the element groups).
//
// Phases may be multi-dimensional scf.parallel loops; their iterations are
// numbered row-major and dealt to the lanes cyclically.
//
// Reductions (scf.parallel with init values and scf.reduce): each lane folds its
// iterations in order, then an xor-shuffle tree over the group's L lanes, in a
// fixed order, so results are reproducible bit for bit; every lane of the group
// ends with the same value. A 32-iteration OR reduction this way is a ballot.
// Stores in the element body itself (outside the phases, e.g. reduction results)
// are done by the group's lane 0.
//
// Barriers: a barrier goes between two phases when some buffer (traced to its
// root through views) is written by one and accessed by the other with no
// barrier in between, and at the end of the element body when some phase uses
// scratch (the next element reuses it).
//
// Why this is correct: iterations of an scf.parallel are independent by
// definition, so any lane may run any of them, and a lane-cyclic split needs no
// synchronisation inside a phase. Phases of one element run in order, and the
// barriers enforce exactly the orders the memory accesses need. Every thread
// reaches every barrier: the element loop depends only on the block, the phases
// sit directly in the element body, and out-of-range groups are masked only at
// their stores. The scratch becomes a {mir.workgroup} argument, which
// --mir-gpu-wrap turns into shared memory: one copy per group.

#include "mir/MirPasses.h"
#include "mir/ViewTypes.h"

#include "mlir/Dialect/Arith/IR/Arith.h"
#include "mlir/Dialect/Func/IR/FuncOps.h"
#include "mlir/Dialect/GPU/IR/GPUDialect.h"
#include "mlir/Dialect/MemRef/IR/MemRef.h"
#include "mlir/Dialect/SCF/IR/SCF.h"
#include "mlir/Dialect/Utils/StaticValueUtils.h"
#include "mlir/IR/Builders.h"
#include "mlir/IR/BuiltinOps.h"
#include "mlir/IR/IRMapping.h"
#include "mlir/Interfaces/SideEffectInterfaces.h"
#include "mlir/Interfaces/ViewLikeInterface.h"
#include "mlir/Pass/Pass.h"

using namespace mlir;

namespace {

constexpr int kWarp = 32;

// Moves the body of a one-dimensional scf.parallel in front of `dest`'s
// terminator, with its induction variable replaced by `iv`.
void moveBody(scf::ParallelOp par, Block *dest, Value iv) {
  Block *src = par.getBody();
  src->getArgument(0).replaceAllUsesWith(iv);
  dest->getOperations().splice(dest->getTerminator()->getIterator(), src->getOperations(), src->begin(),
                               std::prev(src->end()));   // all but the scf.reduce terminator
}

// The buffer a memref value views: a function argument, an alloca or an alloc.
Value rootOf(Value m) {
  while (Operation *def = m.getDefiningOp()) {
    auto view = dyn_cast<ViewLikeOpInterface>(def);
    if (!view)
      break;
    m = view.getViewSource();
  }
  return m;
}

// ceil(x / y) for x >= 0, y > 0, written out (arith.ceildivui needs a separate
// expansion before the LLVM lowering).
Value ceilDiv(OpBuilder &b, Location loc, Value x, Value y) {
  Value one = b.create<arith::ConstantIndexOp>(loc, 1);
  return b.create<arith::DivUIOp>(loc, b.create<arith::AddIOp>(loc, x, b.create<arith::SubIOp>(loc, y, one)), y);
}

enum { kRead = 1, kWrite = 2 };
using Accesses = llvm::MapVector<Value, int>;   // root -> kRead | kWrite

void collect(Operation *phase, Accesses &acc) {
  phase->walk([&](MemoryEffectOpInterface op) {
    SmallVector<MemoryEffects::EffectInstance> effects;
    op.getEffects(effects);
    for (auto &e : effects) {
      Value v = e.getValue();
      if (!v || !isa<MemRefType>(v.getType()))
        continue;
      int &flags = acc[rootOf(v)];
      if (isa<MemoryEffects::Read>(e.getEffect()))
        flags |= kRead;
      if (isa<MemoryEffects::Write>(e.getEffect()))
        flags |= kWrite;
    }
  });
}

bool conflict(const Accesses &a, const Accesses &b) {
  for (auto &[root, fa] : a) {
    auto it = b.find(root);
    if (it != b.end() && ((fa | it->second) & kWrite))
      return true;
  }
  return false;
}

// Packs an element body's scratch buffers into one shared buffer per element
// type, G group parts of S slots each: buffers whose lifetimes (first to last
// access, in the body's top-level order) do not overlap share slots, first fit
// in order of first access. Where a buffer reuses the slots of an earlier one, a
// barrier must separate the last access of the one from the first of the next
// (a lane may still read what another lane would overwrite); one is added if the
// body has none in between. Each buffer becomes a view at group * S + its offset
// (`group` null: one element per block, offset static).
LogicalResult packScratch(func::FuncOp func, Block *body, ArrayRef<memref::AllocaOp> allocas, Value group,
                          int64_t G) {
  DenseMap<Operation *, int> pos;
  int k = 0;
  for (Operation &op : *body) pos[&op] = k++;
  auto top = [&](Operation *op) {
    while (op->getBlock() != body) op = op->getParentOp();
    return op;
  };
  struct Item {
    memref::AllocaOp alloca;
    int64_t size = 0, offset = 0;
    int first = INT_MAX, last = -1;
    Operation *firstOp = nullptr;
  };
  SmallVector<Item> items;
  for (memref::AllocaOp a : allocas) {
    Item it{a, a.getType().getNumElements()};
    SmallVector<Value> work{a.getResult()};
    while (!work.empty()) {
      Value v = work.pop_back_val();
      for (Operation *user : v.getUsers()) {
        if (auto view = dyn_cast<ViewLikeOpInterface>(user); view && view.getViewSource() == v) {
          for (Value r : user->getResults()) work.push_back(r);
          continue;
        }
        Operation *t = top(user);
        if (pos[t] < it.first) {
          it.first = pos[t];
          it.firstOp = t;
        }
        it.last = std::max(it.last, pos[t]);
      }
    }
    if (it.last < 0) {   // never accessed
      it.first = it.last = pos[a];
      it.firstOp = a;
    }
    items.push_back(it);
  }
  SmallVector<int> barriers;
  for (Operation &op : *body)
    if (isa<gpu::BarrierOp>(op))
      barriers.push_back(pos[&op]);
  auto live = [](const Item &x, const Item &y) { return !(x.last < y.first || y.last < x.first); };
  OpBuilder b(func.getContext());
  Location loc = func.getLoc();
  SmallVector<Operation *> needBarrier;
  DenseMap<Type, SmallVector<int>> byType;
  for (int i = 0; i < (int)items.size(); ++i) byType[items[i].alloca.getType().getElementType()].push_back(i);
  for (auto &[type, idx] : byType) {
    llvm::stable_sort(idx, [&](int x, int y) { return items[x].first < items[y].first; });
    SmallVector<int> placed;
    int64_t S = 0;
    for (int i : idx) {
      Item &it = items[i];
      SmallVector<int64_t> candidates{0};
      for (int j : placed)
        if (live(items[j], it)) candidates.push_back(items[j].offset + items[j].size);
      llvm::sort(candidates);
      for (int64_t c : candidates) {
        bool fits = true;
        for (int j : placed)
          if (live(items[j], it) && c < items[j].offset + items[j].size && items[j].offset < c + it.size)
            fits = false;
        if (fits) {
          it.offset = c;
          break;
        }
      }
      for (int j : placed) {   // reused slots: is there a barrier between the two lifetimes?
        const Item &prev = items[j];
        if (live(prev, it) || !(it.offset < prev.offset + prev.size && prev.offset < it.offset + it.size))
          continue;
        const int from = std::min(prev.last, it.last), to = std::max(prev.first, it.first);
        if (!llvm::any_of(barriers, [&](int p) { return from < p && p < to; }))
          needBarrier.push_back(prev.last < it.first ? it.firstOp : prev.firstOp);
      }
      placed.push_back(i);
      S = std::max(S, it.offset + it.size);
    }
    const unsigned at = func.getNumArguments();
    func.insertArgument(at, MemRefType::get({G * S}, type),
                        b.getDictionaryAttr(b.getNamedAttr("mir.workgroup", b.getUnitAttr())), loc);
    Value buffer = func.getArgument(at);
    for (int i : idx) {
      Item &it = items[i];
      MemRefType t = it.alloca.getType();
      OpBuilder vb(it.alloca);
      SmallVector<int64_t> strides(t.getRank(), 1);
      for (int d = (int)t.getRank() - 2; d >= 0; --d) strides[d] = strides[d + 1] * t.getDimSize(d + 1);
      SmallVector<OpFoldResult> sizes, strideFolds;
      for (int d = 0; d < t.getRank(); ++d) {
        sizes.push_back(vb.getIndexAttr(t.getDimSize(d)));
        strideFolds.push_back(vb.getIndexAttr(strides[d]));
      }
      OpFoldResult offset = vb.getIndexAttr(it.offset);
      int64_t staticOffset = it.offset;
      if (group) {
        offset = vb.create<arith::AddIOp>(loc, vb.create<arith::MulIOp>(loc, group, vb.create<arith::ConstantIndexOp>(loc, S)),
                                          vb.create<arith::ConstantIndexOp>(loc, it.offset))
                     .getResult();
        staticOffset = ShapedType::kDynamic;
      }
      auto viewType = MemRefType::get(t.getShape(), t.getElementType(),
                                      StridedLayoutAttr::get(vb.getContext(), staticOffset, strides));
      Value view = vb.create<memref::ReinterpretCastOp>(loc, viewType, buffer, offset, sizes, strideFolds);
      it.alloca.replaceAllUsesWith(view);
      it.alloca.erase();
    }
  }
  for (Operation *op : needBarrier) OpBuilder(op).create<gpu::BarrierOp>(loc);
  return success();
}

LogicalResult mapElements(func::FuncOp func, scf::ParallelOp elements) {
  if (elements.getNumLoops() != 1 || !elements.getInitVals().empty())
    return elements.emitOpError("{mir.elements} must be a one-dimensional loop without reductions");
  if (getConstantIntValue(elements.getLowerBound()[0]) != 0 || getConstantIntValue(elements.getStep()[0]) != 1)
    return elements.emitOpError("{mir.elements} must run from 0 in steps of 1");
  auto attr = [&](StringRef name, int64_t fallback) {
    auto a = elements->getAttrOfType<IntegerAttr>(name);
    return a ? a.getInt() : fallback;
  };
  const int64_t L = attr("mir.lanes_per_element", kWarp);
  const int64_t W = attr("mir.warps_per_block", L == kWarp ? 1 : 4);
  if (L < 1 || L > kWarp || (L & (L - 1)) != 0 || W < 1)
    return elements.emitOpError("lanes per element must be a power of two <= 32, warps per block >= 1");
  const int64_t G = kWarp * W / L;
  const bool single = G == 1;   // one element per block: the original form
  OpBuilder b(elements);
  Location loc = elements.getLoc();

  DenseSet<Value> scratch;
  for (Operation &op : *elements.getBody())
    if (auto a = dyn_cast<memref::AllocaOp>(op); a && a->hasAttr("mir.scratch"))
      scratch.insert(a.getResult());

  Value tid = b.create<gpu::ThreadIdOp>(loc, gpu::Dimension::x);
  Value glane = tid, group;
  if (!single || L != kWarp) {
    glane = b.create<arith::RemUIOp>(loc, tid, b.create<arith::ConstantIndexOp>(loc, L));
    group = b.create<arith::DivUIOp>(loc, tid, b.create<arith::ConstantIndexOp>(loc, L));
  }

  Accesses pending;   // accesses since the last barrier
  bool anyScratch = false;
  for (Operation &op : llvm::make_early_inc_range(*elements.getBody())) {
    auto par = dyn_cast<scf::ParallelOp>(op);
    if (!par || !par->hasAttr("mir.lanes"))
      continue;
    OpBuilder pb(par);
    Location pl = par.getLoc();
    Accesses mine;
    collect(par, mine);
    for (auto &[root, f] : mine) anyScratch |= scratch.contains(root);
    if (conflict(pending, mine)) {
      pb.create<gpu::BarrierOp>(pl);
      pending.clear();
    }
    for (auto &[root, f] : mine) pending[root] |= f;

    // Each lane runs the iterations glane, glane + L, ... (row-major over the
    // phase's dimensions); reductions accumulate in that order per lane, then an
    // xor-shuffle tree combines the group's L lane values in a fixed order.
    const unsigned dims = par.getNumLoops();
    Value start, ub, stride;
    if (dims == 1) {
      Value step = par.getStep()[0];
      start = pb.create<arith::AddIOp>(pl, par.getLowerBound()[0], pb.create<arith::MulIOp>(pl, glane, step));
      stride = pb.create<arith::MulIOp>(pl, step, pb.create<arith::ConstantIndexOp>(pl, L));
      ub = par.getUpperBound()[0];
    } else {
      Value total = pb.create<arith::ConstantIndexOp>(pl, 1);
      for (unsigned d = 0; d < dims; ++d) {
        Value extent = ceilDiv(pb, pl, pb.create<arith::SubIOp>(pl, par.getUpperBound()[d], par.getLowerBound()[d]),
                               par.getStep()[d]);
        total = pb.create<arith::MulIOp>(pl, total, extent);
      }
      start = glane;
      stride = pb.create<arith::ConstantIndexOp>(pl, L);
      ub = total;
    }
    auto loop = pb.create<scf::ForOp>(pl, start, ub, stride, par.getInitVals(),
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
      OpBuilder yb(reduce);
      SmallVector<Value> next;
      for (unsigned i = 0; i < reduce.getNumOperands(); ++i)
        next.push_back(combine(yb, i, loop.getRegionIterArgs()[i], reduce.getOperands()[i]));
      loop.getBody()->getTerminator()->setOperands(next);
    }
    Block *src = par.getBody();
    if (dims == 1) {
      src->getArgument(0).replaceAllUsesWith(loop.getInductionVar());
    } else {   // the row-major position of the linear index q
      OpBuilder ib = OpBuilder::atBlockBegin(loop.getBody());
      Value q = loop.getInductionVar();
      for (int d = (int)dims - 1; d >= 0; --d) {
        Value extent = ceilDiv(ib, pl, ib.create<arith::SubIOp>(pl, par.getUpperBound()[d], par.getLowerBound()[d]),
                               par.getStep()[d]);
        Value k = d == 0 ? q : ib.create<arith::RemUIOp>(pl, q, extent).getResult();
        if (d > 0)
          q = ib.create<arith::DivUIOp>(pl, q, extent);
        Value iv = ib.create<arith::AddIOp>(pl, par.getLowerBound()[d], ib.create<arith::MulIOp>(pl, k, par.getStep()[d]));
        src->getArgument(d).replaceAllUsesWith(iv);
      }
    }
    loop.getBody()->getOperations().splice(loop.getBody()->getTerminator()->getIterator(), src->getOperations(),
                                           src->begin(), std::prev(src->end()));
    SmallVector<Value> results(loop.getResults());
    for (int64_t offset = L / 2; offset >= 1; offset /= 2)
      for (unsigned i = 0; i < results.size(); ++i) {
        Value other = pb.create<gpu::ShuffleOp>(pl, results[i], (int32_t)offset, kWarp, gpu::ShuffleMode::XOR)
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

  // The element loop: one element per block in grid-stride order, or groups of G.
  Value N = elements.getUpperBound()[0];
  Value valid;
  scf::ForOp loop;
  if (single) {
    loop = b.create<scf::ForOp>(loc, b.create<gpu::BlockIdOp>(loc, gpu::Dimension::x), N,
                                b.create<gpu::GridDimOp>(loc, gpu::Dimension::x));
    moveBody(elements, loop.getBody(), loop.getInductionVar());
  } else {
    Value g = b.create<arith::ConstantIndexOp>(loc, G);
    loop = b.create<scf::ForOp>(loc, b.create<gpu::BlockIdOp>(loc, gpu::Dimension::x),
                                ceilDiv(b, loc, N, g),
                                b.create<gpu::GridDimOp>(loc, gpu::Dimension::x));
    OpBuilder ib = OpBuilder::atBlockBegin(loop.getBody());
    Value w = ib.create<arith::AddIOp>(loc, ib.create<arith::MulIOp>(loc, loop.getInductionVar(), g), group);
    valid = ib.create<arith::CmpIOp>(loc, arith::CmpIPredicate::ult, w, N);
    Value last = ib.create<arith::SubIOp>(loc, N, ib.create<arith::ConstantIndexOp>(loc, 1));
    moveBody(elements, loop.getBody(), ib.create<arith::SelectOp>(loc, valid, w, last));
  }
  Block *body = loop.getBody();

  // Stores to anything but scratch: element-level ones once per group (lane 0),
  // and none at all from a group past the last element.
  auto guard = [](memref::StoreOp st, Value cond) {
    OpBuilder sb(st);
    auto ifOp = sb.create<scf::IfOp>(st.getLoc(), cond, /*withElseRegion=*/false);
    st->moveBefore(ifOp.thenBlock()->getTerminator());
  };
  OpBuilder gb = OpBuilder::atBlockBegin(body);
  if (valid)
    gb.setInsertionPointAfter(valid.getDefiningOp());
  Value lane0 = gb.create<arith::CmpIOp>(loc, arith::CmpIPredicate::eq, glane, gb.create<arith::ConstantIndexOp>(loc, 0));
  if (valid)
    lane0 = gb.create<arith::AndIOp>(loc, lane0, valid);
  SmallVector<memref::StoreOp> stores;
  body->walk([&](memref::StoreOp st) {
    if (!scratch.contains(rootOf(st.getMemRef())))
      stores.push_back(st);
  });
  for (memref::StoreOp st : stores) {
    if (st->getBlock() == body)
      guard(st, lane0);
    else if (valid)
      guard(st, valid);
  }

  // The scratch, in shared memory: one buffer per scratch alloca for one element
  // per block with flat buffers (the DSS kernels), else packed (below).
  SmallVector<memref::AllocaOp> allocas;
  for (Operation &op : *body)
    if (auto a = dyn_cast<memref::AllocaOp>(op); a && a->hasAttr("mir.scratch"))
      allocas.push_back(a);
  for (memref::AllocaOp a : allocas)
    if (!a.getType().hasStaticShape() || !a.getType().getLayout().isIdentity())
      return a.emitOpError("{mir.scratch} must be a static buffer with the identity layout");
  const bool flat = llvm::all_of(allocas, [](memref::AllocaOp a) { return a.getType().getRank() == 1; });
  if (single && flat) {
    for (memref::AllocaOp alloca : allocas) {
      const unsigned at = func.getNumArguments();
      func.insertArgument(at, alloca.getType(),
                          b.getDictionaryAttr(b.getNamedAttr("mir.workgroup", b.getUnitAttr())), loc);
      alloca.replaceAllUsesWith(func.getArgument(at));
      alloca.erase();
    }
  } else if (failed(packScratch(func, body, allocas, single ? Value() : group, G))) {
    return failure();
  }
  if (anyScratch)   // the next element must not overwrite scratch a slower lane still reads
    OpBuilder(body->getTerminator()).create<gpu::BarrierOp>(loc);
  elements.erase();
  mir::refreshViewTypes(func.getOperation());
  if (!single || L != kWarp) {
    func->setAttr("mir.block_threads", b.getI64IntegerAttr(kWarp * W));
    func->setAttr("mir.elements_per_block", b.getI64IntegerAttr(G));
  }
  return success();
}

struct WarpPerElementPass : public PassWrapper<WarpPerElementPass, OperationPass<ModuleOp>> {
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(WarpPerElementPass)

  StringRef getArgument() const final { return "mir-warp-per-element"; }
  StringRef getDescription() const final {
    return "Map a {mir.elements} scf.parallel onto groups of lanes (one element "
           "each) and its {mir.lanes} phases onto the group's lanes, with barriers "
           "where phases share memory";
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
