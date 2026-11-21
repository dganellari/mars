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
// Shape: L lanes per element (a power of two up to 32, or a multiple of 32: several
// warps on one element) and W warps per block, so G = 32 W / L elements per block;
// group = thread / L, glane = thread % L. L
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
// are done by the group's lane 0, except to a lane's own buffer (a memref.alloca
// without {mir.scratch}: a scalar every lane computes for itself).
//
// Barriers: a barrier goes between two phases when some buffer (traced to its
// root through views) is written by one and accessed by the other with no
// barrier in between, and at the end of the element body when some phase uses
// scratch (the next element reuses it).
//
// Loops of phases: an scf.for in the element body may hold phases (the planes of a
// z-march schedule) when its trip count is a constant, at least one, and it holds
// nothing else that touches memory. Barriers inside it follow the same rule; at the
// end of each trip one goes in when the accesses since the last barrier conflict
// with any access of the loop, since the next trip's phases may touch those buffers.
//
// Why this is correct: iterations of an scf.parallel are independent by
// definition, so any lane may run any of them, and a lane-cyclic split needs no
// synchronisation inside a phase. Phases of one element run in order, and the
// barriers enforce exactly the orders the memory accesses need. Every thread
// reaches every barrier: the element loop depends only on the block, the phases
// sit in the element body or in a loop of phases whose trip count is the same for
// every lane, and out-of-range groups are masked only at their stores. The scratch becomes a {mir.workgroup} argument, which
// --mir-gpu-wrap turns into shared memory: one copy per group.

#include "mir/MirPasses.h"
#include "mir/Target.h"
#include "mir/ViewTypes.h"

#include "mlir/Dialect/Arith/IR/Arith.h"
#include "mlir/Dialect/Func/IR/FuncOps.h"
#include "mlir/Dialect/GPU/IR/GPUDialect.h"
#include "mlir/Dialect/MemRef/IR/MemRef.h"
#include "mlir/Dialect/SCF/IR/SCF.h"
#include "mlir/Dialect/SCF/Utils/Utils.h"
#include "mlir/Dialect/Utils/StaticValueUtils.h"
#include "mlir/IR/Builders.h"
#include "mlir/IR/BuiltinOps.h"
#include "mlir/IR/IRMapping.h"
#include "mlir/Interfaces/SideEffectInterfaces.h"
#include "mlir/Interfaces/ViewLikeInterface.h"
#include "mlir/Pass/Pass.h"

#include "llvm/ADT/SmallBitVector.h"

#include <map>
#include <optional>
#include <set>

using namespace mlir;

namespace {

// The bank model below is NVIDIA's (32 banks of 4 bytes); the lanes per warp come from the
// target (mir::getGpuTarget).
constexpr int kBanks = 32;

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

// Scratch layout: the strides each buffer gets in shared memory. Row-major is the
// default; padding grows the stride of one dimension (the dimensions after it stay
// one contiguous block). Register blocking keeps one axis inside the thread, so a
// warp's lanes can spread over axes whose power-of-two strides land on a few of the
// 32 banks (measured: 4- and 8-way conflicts at p = 7 and 3), while other buffers
// are walked row by row, where the same padding would cause conflicts. So each
// buffer's layout is chosen from its own accesses: every load and store of the
// buffer is evaluated for the lanes of one element at the first iteration of its
// loops (indices are integer expressions of the thread id and loop counters),
// counted in wavefronts as the hardware serves them, and weighted by how often it
// runs; the cheapest candidate wins, row-major on ties. A collapse_shape of the
// buffer fixes the boundaries it merges.
struct Layout {
  SmallVector<int64_t> strides;
  int64_t size = 0;
};

Layout rowMajor(ArrayRef<int64_t> shape) {
  Layout l;
  l.strides.assign(shape.size(), 1);
  for (int d = (int)shape.size() - 2; d >= 0; --d) l.strides[d] = l.strides[d + 1] * shape[d + 1];
  l.size = shape.empty() ? 1 : l.strides[0] * shape[0];
  return l;
}

// Dimension d's stride grows by pad; the dimensions after it stay row-major.
Layout padded(ArrayRef<int64_t> shape, int64_t d, int64_t pad) {
  Layout l = rowMajor(shape);
  l.strides[d] += pad;
  for (int k = (int)d - 1; k >= 0; --k) l.strides[k] = l.strides[k + 1] * shape[k + 1];
  l.size = l.strides[0] * shape[0];
  return l;
}

// The value of an index expression for thread `tid`, loop counters at their first
// iteration; none when it depends on anything else.
struct IndexEval {
  explicit IndexEval(int64_t t) : tid(t) {}
  int64_t tid;
  DenseMap<Value, std::optional<int64_t>> memo;

  std::optional<int64_t> operator()(OpFoldResult f) {
    if (auto a = dyn_cast<Attribute>(f))
      return cast<IntegerAttr>(a).getInt();
    return (*this)(cast<Value>(f));
  }

  std::optional<int64_t> operator()(Value v) {
    if (auto it = memo.find(v); it != memo.end())
      return it->second;
    std::optional<int64_t> r = compute(v);
    memo[v] = r;
    return r;
  }

  std::optional<int64_t> compute(Value v) {
    if (auto arg = dyn_cast<BlockArgument>(v)) {
      auto loop = dyn_cast<scf::ForOp>(arg.getOwner()->getParentOp());
      if (loop && arg == loop.getInductionVar())
        return (*this)(loop.getLowerBound());
      return std::nullopt;
    }
    Operation *op = v.getDefiningOp();
    if (auto c = dyn_cast<arith::ConstantOp>(op)) {
      if (auto i = dyn_cast<IntegerAttr>(c.getValue()))
        return i.getValue().getSExtValue();
      return std::nullopt;
    }
    if (isa<gpu::ThreadIdOp>(op))
      return tid;
    if (isa<gpu::BlockIdOp>(op))
      return 0;
    if (isa<arith::IndexCastOp, arith::IndexCastUIOp, arith::ExtSIOp, arith::ExtUIOp, arith::TruncIOp>(op))
      return (*this)(op->getOperand(0));
    if (auto sel = dyn_cast<arith::SelectOp>(op)) {
      std::optional<int64_t> c = (*this)(sel.getCondition());
      if (!c)
        return std::nullopt;
      return (*this)(*c ? sel.getTrueValue() : sel.getFalseValue());
    }
    if (op->getNumOperands() != 2)
      return std::nullopt;
    std::optional<int64_t> a = (*this)(op->getOperand(0)), b = (*this)(op->getOperand(1));
    if (!a || !b)
      return std::nullopt;
    if (auto cmp = dyn_cast<arith::CmpIOp>(op)) {
      switch (cmp.getPredicate()) {
      case arith::CmpIPredicate::eq: return *a == *b;
      case arith::CmpIPredicate::ne: return *a != *b;
      case arith::CmpIPredicate::slt: case arith::CmpIPredicate::ult: return *a < *b;
      case arith::CmpIPredicate::sle: case arith::CmpIPredicate::ule: return *a <= *b;
      case arith::CmpIPredicate::sgt: case arith::CmpIPredicate::ugt: return *a > *b;
      case arith::CmpIPredicate::sge: case arith::CmpIPredicate::uge: return *a >= *b;
      }
    }
    if (isa<arith::AddIOp>(op)) return *a + *b;
    if (isa<arith::SubIOp>(op)) return *a - *b;
    if (isa<arith::MulIOp>(op)) return *a * *b;
    if (isa<arith::MinSIOp, arith::MinUIOp>(op)) return std::min(*a, *b);
    if (isa<arith::MaxSIOp, arith::MaxUIOp>(op)) return std::max(*a, *b);
    if (*b == 0)
      return std::nullopt;
    if (isa<arith::DivUIOp, arith::DivSIOp>(op)) return *a / *b;
    if (isa<arith::RemUIOp, arith::RemSIOp>(op)) return *a % *b;
    return std::nullopt;
  }
};

// The indices into `root` that an access to its view `m` at `idx` touches.
std::optional<SmallVector<int64_t>> rootIndices(Value m, SmallVector<int64_t> idx, Value root, IndexEval &eval) {
  while (m != root) {
    Operation *def = m.getDefiningOp();
    if (auto sv = dyn_cast_or_null<memref::SubViewOp>(def)) {
      const llvm::SmallBitVector dropped = sv.getDroppedDims();
      SmallVector<int64_t> src;
      unsigned k = 0;
      for (unsigned d = 0; d < (unsigned)sv.getSourceType().getRank(); ++d) {
        std::optional<int64_t> off = eval(sv.getMixedOffsets()[d]), st = eval(sv.getMixedStrides()[d]);
        if (!off || !st)
          return std::nullopt;
        src.push_back(*off + (dropped.test(d) ? 0 : idx[k++]) * *st);
      }
      idx = src;
      m = sv.getSource();
    } else if (auto ex = dyn_cast_or_null<memref::ExpandShapeOp>(def)) {
      ArrayRef<int64_t> shape = ex.getResultType().getShape();
      SmallVector<int64_t> src;
      for (const ReassociationIndices &group : ex.getReassociationIndices()) {
        int64_t lin = 0;
        for (int64_t d : group) lin = lin * shape[d] + idx[d];
        src.push_back(lin);
      }
      idx = src;
      m = ex.getSrc();
    } else if (auto co = dyn_cast_or_null<memref::CollapseShapeOp>(def)) {
      ArrayRef<int64_t> shape = co.getSrcType().getShape();
      SmallVector<int64_t> src(shape.size());
      SmallVector<ReassociationIndices> groups = co.getReassociationIndices();
      for (size_t g = 0; g < groups.size(); ++g) {
        int64_t lin = idx[g];
        for (int k = (int)groups[g].size() - 1; k >= 0; --k) {
          src[groups[g][k]] = lin % shape[groups[g][k]];
          lin /= shape[groups[g][k]];
        }
      }
      idx = src;
      m = co.getSrc();
    } else if (auto cast = dyn_cast_or_null<memref::CastOp>(def)) {
      m = cast.getSource();
    } else {
      return std::nullopt;
    }
  }
  return idx;
}

// Shared-memory wavefronts of one warp instruction: 32 banks of 4 bytes, an access
// wider than 4 bytes served per half (8 bytes) or quarter warp (16), lanes reading
// the same word sharing it.
int64_t wavefronts(ArrayRef<std::pair<int64_t, int64_t>> laneBytes, int64_t bytes) {
  const int64_t span = bytes <= 4 ? kBanks : 128 / bytes, words = std::max<int64_t>(1, bytes / 4);
  std::map<int64_t, std::map<int64_t, std::set<int64_t>>> passes;   // pass -> bank -> words
  for (auto [lane, byte] : laneBytes)
    for (int64_t w = 0; w < words; ++w) passes[lane % kBanks / span][(byte / 4 + w) % kBanks].insert(byte / 4 + w);
  int64_t total = 0;
  for (auto &[p, banks] : passes) {
    size_t worst = 0;
    for (auto &[bank, set] : banks) worst = std::max(worst, set.size());
    total += (int64_t)worst;
  }
  return total;
}

// The layout of one scratch buffer (see Layout above). `lanes`: lanes per element.
Layout chooseLayout(memref::AllocaOp alloca, Block *body, int64_t lanes) {
  MemRefType type = alloca.getType();
  ArrayRef<int64_t> shape = type.getShape();
  const int64_t rank = shape.size();
  const int64_t bytes = std::max<int64_t>(1, type.getElementTypeBitWidth() / 8);
  Layout best = rowMajor(shape);
  if (rank < 2)
    return best;
  // Accesses (loads and stores through views) and the boundaries a collapse fixes.
  SmallVector<bool> fixed(rank - 1, false);
  struct Use {
    Value memref;
    Operation *op;
    ValueRange indices;
  };
  SmallVector<Use> uses;
  SmallVector<Value> work{alloca.getResult()};
  while (!work.empty()) {
    Value v = work.pop_back_val();
    for (Operation *user : v.getUsers()) {
      if (auto c = dyn_cast<memref::CollapseShapeOp>(user)) {
        if (v == alloca.getResult())
          for (const ReassociationIndices &group : c.getReassociationIndices())
            for (size_t k = 0; k + 1 < group.size(); ++k) fixed[group[k]] = true;
        else
          fixed.assign(fixed.size(), true);   // a reshape behind another view: not followed
      }
      if (auto view = dyn_cast<ViewLikeOpInterface>(user); view && view.getViewSource() == v) {
        for (Value r : user->getResults()) work.push_back(r);
      } else if (auto l = dyn_cast<memref::LoadOp>(user)) {
        uses.push_back({v, user, l.getIndices()});
      } else if (auto st = dyn_cast<memref::StoreOp>(user)) {
        uses.push_back({v, user, st.getIndices()});
      }
    }
  }
  // Root indices of every access for the lanes of one element, and its weight.
  struct Sample {
    SmallVector<std::pair<int64_t, SmallVector<int64_t>>> lanes;   // (tid, root indices)
    int64_t weight = 1;
  };
  SmallVector<Sample> samples;
  for (const Use &u : uses) {
    Sample s;
    for (int64_t tid = 0; tid < std::min<int64_t>(lanes, kBanks); ++tid) {
      IndexEval eval(tid);
      SmallVector<int64_t> idx;
      for (Value i : u.indices)
        if (std::optional<int64_t> x = eval(i))
          idx.push_back(*x);
      if (idx.size() != u.indices.size())
        break;
      if (std::optional<SmallVector<int64_t>> r = rootIndices(u.memref, idx, alloca.getResult(), eval))
        s.lanes.push_back({tid, *r});
    }
    IndexEval eval0(0);
    for (Operation *p = u.op->getParentOp(); p && p != body->getParentOp(); p = p->getParentOp())
      if (auto loop = dyn_cast<scf::ForOp>(p)) {
        std::optional<int64_t> lb = eval0(loop.getLowerBound()), ub = eval0(loop.getUpperBound()),
                               st = eval0(loop.getStep());
        if (lb && ub && st && *st > 0 && *ub > *lb)
          s.weight *= (*ub - *lb + *st - 1) / *st;
      }
    if (!s.lanes.empty())
      samples.push_back(std::move(s));
  }
  auto cost = [&](const Layout &l) {
    int64_t total = 0;
    for (const Sample &s : samples) {
      SmallVector<std::pair<int64_t, int64_t>> laneBytes;
      for (auto &[tid, idx] : s.lanes) {
        int64_t lin = 0;
        for (int64_t d = 0; d < rank; ++d) lin += idx[d] * l.strides[d];
        laneBytes.push_back({tid, lin * bytes});
      }
      total += s.weight * wavefronts(laneBytes, bytes);
    }
    return total;
  };
  int64_t bestCost = cost(best);
  for (int64_t d = rank - 2; d >= 0; --d) {
    int64_t inner = 1;
    for (int64_t k = d + 1; k < rank; ++k) inner *= shape[k];
    if (fixed[d] || inner < 4)
      continue;
    for (int64_t pad : {1, 2}) {
      Layout l = padded(shape, d, pad);
      const int64_t c = cost(l);
      if (c < bestCost || (c == bestCost && l.size < best.size)) {
        best = l;
        bestCost = c;
      }
    }
  }
  return best;
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
                          int64_t G, int64_t lanes, bool pad, bool largestFirstAllowed) {
  DenseMap<Operation *, int> pos;
  int k = 0;
  for (Operation &op : *body) pos[&op] = k++;
  auto top = [&](Operation *op) {
    while (op->getBlock() != body) op = op->getParentOp();
    return op;
  };
  struct Item {
    memref::AllocaOp alloca;
    Layout layout;
    int64_t size = 0, offset = 0;
    int first = INT_MAX, last = -1;
    Operation *firstOp = nullptr;
  };
  SmallVector<Item> items;
  for (memref::AllocaOp a : allocas) {
    Layout layout = pad ? chooseLayout(a, body, lanes) : rowMajor(a.getType().getShape());
    Item it{a, layout, layout.size};
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
  // First fit of the buffers in the given order; returns the slots it needs.
  auto pack = [&](ArrayRef<int> order) {
    SmallVector<int> placed;
    int64_t S = 0;
    for (int i : order) {
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
      placed.push_back(i);
      S = std::max(S, it.offset + it.size);
    }
    return S;
  };
  for (auto &[type, idx] : byType) {
    llvm::stable_sort(idx, [&](int x, int y) { return items[x].first < items[y].first; });
    // Order of first access, or largest first when that needs fewer slots: a
    // buffer placed by first access can sit on slots a larger, later one then
    // has to go around (a Galerkin operator's interpolated u after a dead sweep).
    SmallVector<int> bySize(idx);
    llvm::stable_sort(bySize, [&](int x, int y) { return items[x].size > items[y].size; });
    const bool largestFirst = largestFirstAllowed && pack(bySize) < pack(idx);
    ArrayRef<int> order = largestFirst ? ArrayRef<int>(bySize) : ArrayRef<int>(idx);
    int64_t S = pack(order);   // leaves the chosen order's offsets
    for (int n = 0; n < (int)order.size(); ++n) {   // reused slots: is there a barrier between the two lifetimes?
      const Item &it = items[order[n]];
      for (int m = 0; m < n; ++m) {
        const Item &prev = items[order[m]];
        if (live(prev, it) || !(it.offset < prev.offset + prev.size && prev.offset < it.offset + it.size))
          continue;
        const int from = std::min(prev.last, it.last), to = std::max(prev.first, it.first);
        if (!llvm::any_of(barriers, [&](int p) { return from < p && p < to; }))
          needBarrier.push_back(prev.last < it.first ? it.firstOp : prev.firstOp);
      }
    }
    // Several elements' lanes are served in one pass over the 32 banks (a warp for
    // 4-byte values, a half warp for 8-byte ones) when lanes x words < 32. Their
    // slots then must not start on the same bank, or the same access of each element
    // conflicts (measured: 2-way on every access at p = 1). S = lanes (mod 32 /
    // words) puts them lanes x words words apart, so they tile the banks.
    const int64_t words = std::max<int64_t>(1, (int64_t)type.getIntOrFloatBitWidth() / 32);
    if (group && lanes * words < 32) {
      const int64_t period = 32 / words;
      S += ((lanes - S) % period + period) % period;
    }
    const unsigned at = func.getNumArguments();
    func.insertArgument(at, MemRefType::get({G * S}, type),
                        b.getDictionaryAttr(b.getNamedAttr("mir.workgroup", b.getUnitAttr())), loc);
    Value buffer = func.getArgument(at);
    for (int i : idx) {
      Item &it = items[i];
      MemRefType t = it.alloca.getType();
      OpBuilder vb(it.alloca);
      ArrayRef<int64_t> strides = it.layout.strides;
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


// Lane-owned scratch (option lane-owned): a buffer every phase touches only at points the
// accessing lane owns lives in the lane's registers instead of shared memory. It is lane-owned
// when every access sits in a phase (no reductions) over one common iteration space and indexes
// the buffer, through its views, with each of the phase's induction variables exactly once, in
// the same positions in every phase, plus constants (dimensions kept in the thread): point q of
// the space then runs on lane q mod L in all of them and owns its own elements. A phase that
// only stores one value at its induction variables over the buffer's whole shape is a fill. Each
// lane keeps its points in J = ceil(points / L) slots: a private array [J][constant dims...],
// indexed by constants once the slots are unrolled, which LLVM keeps in registers.
std::optional<int64_t> staticIndex(Value v, int depth = 0) {
  if (auto c = getConstantIntValue(v))
    return *c;
  Operation *d = v.getDefiningOp();
  if (!d || depth > 8 || d->getNumOperands() != 2)
    return std::nullopt;
  auto a = staticIndex(d->getOperand(0), depth + 1), b = staticIndex(d->getOperand(1), depth + 1);
  if (!a || !b)
    return std::nullopt;
  if (isa<arith::AddIOp>(d)) return *a + *b;
  if (isa<arith::SubIOp>(d)) return *a - *b;
  if (isa<arith::MulIOp>(d)) return *a * *b;
  return std::nullopt;
}

struct OwnedIndex {
  bool iv;         // the phase's induction variable at position `value`, else the constant `value`
  int64_t value;
};

// Indices of an access through views of `root`, as indices of `root`.
std::optional<SmallVector<OwnedIndex>> toRoot(Value m, SmallVector<OwnedIndex> idx, Value root) {
  if (m == root)
    return idx;
  Operation *def = m.getDefiningOp();
  if (auto sv = dyn_cast_or_null<memref::SubViewOp>(def)) {
    auto src = cast<MemRefType>(sv.getSource().getType());
    llvm::SmallBitVector dropped = sv.getDroppedDims();
    SmallVector<OpFoldResult> offs = sv.getMixedOffsets(), strides = sv.getMixedStrides();
    SmallVector<OwnedIndex> out;
    unsigned k = 0;
    for (int64_t d = 0; d < src.getRank(); ++d) {
      std::optional<int64_t> off = isa<Attribute>(offs[d]) ? getConstantIntValue(offs[d])
                                                           : staticIndex(cast<Value>(offs[d]));
      if (!off || getConstantIntValue(strides[d]) != std::optional<int64_t>(1))
        return std::nullopt;
      if (dropped.test(d)) {
        out.push_back({false, *off});
        continue;
      }
      OwnedIndex r = idx[k++];
      if (r.iv && *off != 0)
        return std::nullopt;
      out.push_back(r.iv ? r : OwnedIndex{false, r.value + *off});
    }
    return toRoot(sv.getSource(), out, root);
  }
  if (auto cs = dyn_cast_or_null<memref::CollapseShapeOp>(def)) {
    ArrayRef<int64_t> shape = cs.getSrcType().getShape();
    SmallVector<OwnedIndex> out(shape.size(), OwnedIndex{false, 0});
    for (auto [g, group] : llvm::enumerate(cs.getReassociationIndices())) {
      SmallVector<int64_t> wide;
      for (int64_t d : group) if (shape[d] != 1) wide.push_back(d);
      if (wide.size() > 1)
        return std::nullopt;
      if (wide.size() == 1)
        out[wide[0]] = idx[g];
    }
    return toRoot(cs.getSrc(), out, root);
  }
  if (auto es = dyn_cast_or_null<memref::ExpandShapeOp>(def)) {
    ArrayRef<int64_t> shape = es.getResultType().getShape();
    SmallVector<OwnedIndex> out;
    for (const ReassociationIndices &group : es.getReassociationIndices()) {
      std::optional<OwnedIndex> one;
      for (int64_t d : group) {
        if (shape[d] == 1)
          continue;
        if (one)
          return std::nullopt;
        one = idx[d];
      }
      out.push_back(one.value_or(OwnedIndex{false, 0}));
    }
    return toRoot(es.getSrc(), out, root);
  }
  return std::nullopt;
}

struct Owned {
  SmallVector<int64_t> space;             // the common iteration space
  SmallVector<int64_t> constDims;         // root dimensions indexed by constants
  DenseMap<Operation *, SmallVector<int64_t>> at;   // access -> its constants, in constDims order
  DenseMap<Operation *, Value> fills;     // fill phase -> the value it stores
  int64_t slots = 0;                      // J
  Value priv;
};

std::optional<Owned> laneOwned(memref::AllocaOp alloca, Block *body, int64_t L) {
  MemRefType type = alloca.getType();
  if (!type.hasStaticShape() || type.getRank() == 0)
    return std::nullopt;
  Owned o;
  std::optional<SmallVector<OwnedIndex>> signature;
  SmallVector<Value> work{alloca.getResult()};
  while (!work.empty()) {
    Value v = work.pop_back_val();
    for (Operation *user : v.getUsers()) {
      if (isa<memref::SubViewOp, memref::CollapseShapeOp, memref::ExpandShapeOp>(user)) {
        work.push_back(user->getResult(0));
        continue;
      }
      Value mem;
      SmallVector<Value> indices;
      if (auto ld = dyn_cast<memref::LoadOp>(user)) mem = ld.getMemRef(), indices = ld.getIndices();
      else if (auto st = dyn_cast<memref::StoreOp>(user)) mem = st.getMemRef(), indices = st.getIndices();
      else return std::nullopt;
      auto phase = dyn_cast<scf::ParallelOp>(user->getParentOp());
      if (!phase || !phase->hasAttr("mir.lanes") || phase->getBlock() != body || !phase.getInitVals().empty())
        return std::nullopt;
      SmallVector<int64_t> space;
      for (unsigned d = 0; d < phase.getNumLoops(); ++d) {
        auto lb = getConstantIntValue(phase.getLowerBound()[d]), ub = getConstantIntValue(phase.getUpperBound()[d]),
             st = getConstantIntValue(phase.getStep()[d]);
        if (lb != std::optional<int64_t>(0) || st != std::optional<int64_t>(1) || !ub)
          return std::nullopt;
        space.push_back(*ub);
      }
      SmallVector<OwnedIndex> idx;
      for (Value i : indices) {
        int pos = -1;
        for (unsigned d = 0; d < phase.getNumLoops(); ++d)
          if (phase.getInductionVars()[d] == i) pos = (int)d;
        if (pos >= 0) {
          idx.push_back({true, pos});
        } else if (auto c = staticIndex(i)) {
          idx.push_back({false, *c});
        } else {
          return std::nullopt;
        }
      }
      std::optional<SmallVector<OwnedIndex>> root = toRoot(mem, idx, alloca.getResult());
      if (!root)
        return std::nullopt;
      // A fill: the phase stores one value at its induction variables, over the whole buffer.
      if (auto st = dyn_cast<memref::StoreOp>(user);
          st && std::distance(phase.getBody()->begin(), phase.getBody()->end()) == 2 &&
          st.getValueToStore().getDefiningOp() && !phase->isAncestor(st.getValueToStore().getDefiningOp()) &&
          llvm::equal(space, type.getShape()) &&
          llvm::all_of(llvm::enumerate(*root), [](auto e) { return e.value().iv && e.value().value == (int64_t)e.index(); })) {
        o.fills[phase] = st.getValueToStore();
        continue;
      }
      // Each induction variable exactly once; the same positions and space in every phase.
      SmallVector<int> seen(phase.getNumLoops(), 0);
      for (const OwnedIndex &r : *root) if (r.iv) ++seen[r.value];
      if (llvm::any_of(seen, [](int n) { return n != 1; }))
        return std::nullopt;
      if (!signature) {
        signature = root;
        o.space = space;
        for (auto [d, r] : llvm::enumerate(*root)) if (!r.iv) o.constDims.push_back(d);
      } else {
        if (space != o.space)
          return std::nullopt;
        for (auto [a, b] : llvm::zip(*root, *signature))
          if (a.iv != b.iv || (a.iv && a.value != b.value))
            return std::nullopt;
      }
      SmallVector<int64_t> consts;
      for (int64_t d : o.constDims) consts.push_back((*root)[d].value);
      o.at[user] = consts;
    }
  }
  if (!signature)
    return std::nullopt;
  int64_t points = 1;
  for (int64_t s : o.space) points *= s;
  o.slots = (points + L - 1) / L;
  return o;
}

LogicalResult mapElements(func::FuncOp func, scf::ParallelOp elements, bool padScratch, bool largestFirst,
                          bool promoteOwned = false) {
  if (elements.getNumLoops() != 1 || !elements.getInitVals().empty())
    return elements.emitOpError("{mir.elements} must be a one-dimensional loop without reductions");
  if (getConstantIntValue(elements.getLowerBound()[0]) != 0 || getConstantIntValue(elements.getStep()[0]) != 1)
    return elements.emitOpError("{mir.elements} must run from 0 in steps of 1");
  auto attr = [&](StringRef name, int64_t fallback) {
    auto a = elements->getAttrOfType<IntegerAttr>(name);
    return a ? a.getInt() : fallback;
  };
  const int64_t wave = mir::getGpuTarget(func).waveSize;
  const int64_t L = attr("mir.lanes_per_element", wave);
  const int64_t W = attr("mir.warps_per_block", L == wave ? 1 : 4);
  // Several warps per element keep more warps resident when scratch per element limits
  // the blocks per SM (q = 9: four q^3 buffers per BP3 element).
  const bool multiWarp = L > wave;
  if (L < 1 || (!multiWarp && (L & (L - 1)) != 0) || (multiWarp && L % wave != 0) || W < 1 ||
      (wave * W) % L != 0)
    return elements.emitOpError("lanes per element must be a power of two <= the wave size or a multiple "
                                "of it dividing the block, warps per block >= 1");
  const int64_t G = wave * W / L;
  const bool single = G == 1;   // one element per block: the original form
  OpBuilder b(elements);
  Location loc = elements.getLoc();

  // Lane-owned scratch into registers: loops of phases are unrolled first, so a plane index
  // becomes a constant, then each scratch buffer every phase touches lane by lane gets a
  // private array (laneOwned above).
  DenseMap<Value, Owned> owned;   // scratch alloca -> its private form
  if (promoteOwned) {
    SmallVector<scf::ForOp> loops;
    for (Operation &op : *elements.getBody())
      if (auto f = dyn_cast<scf::ForOp>(op);
          f && llvm::any_of(f.getBody()->getOperations(), [](Operation &o) { return o.hasAttr("mir.lanes"); }))
        loops.push_back(f);
    for (scf::ForOp f : loops) {
      auto lo = getConstantIntValue(f.getLowerBound()), hi = getConstantIntValue(f.getUpperBound()),
           st = getConstantIntValue(f.getStep());
      if (!lo || !hi || !st || failed(loopUnrollByFactor(f, (uint64_t)((*hi - *lo + *st - 1) / *st))))
        return f.emitOpError("cannot unroll a loop of phases for lane-owned scratch");
    }
    // The private arrays go before the element loop: an alloca inside it is a dynamic stack
    // allocation on every element (on AMDGPU it overruns the private segment, and LLVM keeps
    // only entry-block allocas in registers). Each element's fill re-initializes them.
    OpBuilder ab(elements);
    for (Operation &op : *elements.getBody())
      if (auto a = dyn_cast<memref::AllocaOp>(op); a && a->hasAttr("mir.scratch"))
        if (std::optional<Owned> o = laneOwned(a, elements.getBody(), L)) {
          SmallVector<int64_t> shape{o->slots};
          for (int64_t d : o->constDims) shape.push_back(a.getType().getDimSize(d));
          o->priv = ab.create<memref::AllocaOp>(a.getLoc(), MemRefType::get(shape, a.getType().getElementType()));
          owned[a.getResult()] = std::move(*o);
        }
  }

  // A lane's own buffers (allocas without {mir.scratch}, e.g. the per-element scalars of a
  // trilinear map) go before the element loop on AMD: inside it they are dynamic stack
  // allocations on every element, which overrun AMDGPU's private segment. Each element writes
  // them before reading. (NVIDIA's lowering copes, and its PTX stays as it is.)
  if (mir::getGpuTarget(func).isAmd())
    for (Operation &op : llvm::make_early_inc_range(*elements.getBody()))
      if (auto a = dyn_cast<memref::AllocaOp>(op); a && !a->hasAttr("mir.scratch") && a->getNumOperands() == 0)
        a->moveBefore(elements);

  DenseSet<Value> scratch;
  for (Operation &op : *elements.getBody())
    if (auto a = dyn_cast<memref::AllocaOp>(op); a && a->hasAttr("mir.scratch"))
      scratch.insert(a.getResult());

  Value tid = b.create<gpu::ThreadIdOp>(loc, gpu::Dimension::x);
  Value glane = tid, group;
  if (!single || L != wave) {
    glane = b.create<arith::RemUIOp>(loc, tid, b.create<arith::ConstantIndexOp>(loc, L));
    group = b.create<arith::DivUIOp>(loc, tid, b.create<arith::ConstantIndexOp>(loc, L));
  }

  Accesses pending;   // accesses since the last barrier
  bool anyScratch = false;
  // A barrier before a phase that conflicts with what ran since the last one.
  auto order = [&](scf::ParallelOp par) {
    Accesses mine;
    collect(par, mine);
    for (auto &[root, o] : owned) mine.erase(root);   // lane-local by construction
    for (auto &[root, f] : mine) anyScratch |= scratch.contains(root);
    if (conflict(pending, mine)) {
      OpBuilder(par).create<gpu::BarrierOp>(par.getLoc());
      pending.clear();
    }
    for (auto &[root, f] : mine) pending[root] |= f;
  };
  // One phase onto the group's lanes.
  auto mapPhase = [&](scf::ParallelOp par) -> LogicalResult {
    OpBuilder pb(par);
    Location pl = par.getLoc();
    // Each lane runs the iterations glane, glane + L, ... (row-major over the
    // phase's dimensions); reductions accumulate in that order per lane, then an
    // xor-shuffle tree combines the group's L lane values in a fixed order.
    const unsigned dims = par.getNumLoops();
    if (multiWarp && !par.getInitVals().empty())
      return par.emitOpError("reductions need the element's lanes in one warp");
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
        Value other = pb.create<gpu::ShuffleOp>(pl, results[i], (int32_t)offset, (int32_t)wave, gpu::ShuffleMode::XOR)
                          .getShuffleResult();
        results[i] = combine(pb, i, results[i], other);
      }
    par.replaceAllUsesWith(results);
    par.erase();
    return success();
  };
  auto isPhase = [](Operation &op) {
    auto p = dyn_cast<scf::ParallelOp>(op);
    return p && p->hasAttr("mir.lanes");
  };
  // A phase touching lane-owned buffers: its J slots unrolled, slot j running point
  // q = glane + j L when q is in range, the owned accesses on the private arrays at [j, ...].
  auto mapOwned = [&](scf::ParallelOp par) -> LogicalResult {
    OpBuilder pb(par);
    Location pl = par.getLoc();
    SmallVector<int64_t> space;
    for (Value ub : par.getUpperBound()) space.push_back(*getConstantIntValue(ub));
    int64_t points = 1;
    for (int64_t s : space) points *= s;
    const int64_t slots = (points + L - 1) / L;
    for (int64_t j = 0; j < slots; ++j) {
      Value q = pb.create<arith::AddIOp>(pl, glane, pb.create<arith::ConstantIndexOp>(pl, j * L));
      Value in = pb.create<arith::CmpIOp>(pl, arith::CmpIPredicate::ult, q, pb.create<arith::ConstantIndexOp>(pl, points));
      auto guard = pb.create<scf::IfOp>(pl, in, /*withElseRegion=*/false);
      OpBuilder gb = guard.getThenBodyBuilder();
      IRMapping map;
      Value rest = q;
      for (int d = (int)space.size() - 1; d >= 0; --d) {   // row-major position of q
        Value extent = gb.create<arith::ConstantIndexOp>(pl, space[d]);
        map.map(par.getInductionVars()[d], d == 0 ? rest : gb.create<arith::RemUIOp>(pl, rest, extent).getResult());
        if (d > 0) rest = gb.create<arith::DivUIOp>(pl, rest, extent);
      }
      for (Operation &o : par.getBody()->without_terminator()) {
        const Owned *own = nullptr;
        const SmallVector<int64_t> *consts = nullptr;
        for (auto &[root, ow] : owned)
          if (auto it = ow.at.find(&o); it != ow.at.end()) own = &ow, consts = &it->second;
        if (!own) {
          gb.clone(o, map);
          continue;
        }
        SmallVector<Value> idx{gb.create<arith::ConstantIndexOp>(pl, j)};
        for (int64_t c : *consts) idx.push_back(gb.create<arith::ConstantIndexOp>(pl, c));
        if (auto ld = dyn_cast<memref::LoadOp>(o))
          map.map(ld.getResult(), gb.create<memref::LoadOp>(pl, own->priv, idx).getResult());
        else
          gb.create<memref::StoreOp>(pl, map.lookupOrDefault(cast<memref::StoreOp>(o).getValueToStore()), own->priv, idx);
      }
    }
    // The accesses die with the phase: no stale pointer may match a later op.
    for (Operation &o : par.getBody()->without_terminator())
      for (auto &[root, ow] : owned) ow.at.erase(&o);
    par.erase();
    return success();
  };
  // Phases with lane-owned accesses, fixed before any phase is mapped (and erased).
  DenseSet<Operation *> ownedPhases;
  for (auto &[root, ow] : owned)
    for (auto &[access, consts] : ow.at) ownedPhases.insert(access->getParentOp());
  for (Operation &op : llvm::make_early_inc_range(*elements.getBody())) {
    if (isPhase(op)) {
      auto par = cast<scf::ParallelOp>(op);
      // A fill of a lane-owned buffer: every private entry of every lane.
      const Owned *filled = nullptr;
      for (auto &[root, ow] : owned) if (ow.fills.count(par)) filled = &ow;
      if (filled) {
        const Value fillValue = filled->fills.lookup(par);
        for (auto &[root, ow] : owned) ow.fills.erase(par);
        OpBuilder fb(par);
        auto pt = cast<MemRefType>(filled->priv.getType());
        int64_t entries = pt.getNumElements();
        for (int64_t e = 0; e < entries; ++e) {
          SmallVector<Value> idx;
          int64_t rest = e;
          SmallVector<int64_t> pos(pt.getRank());
          for (int d = (int)pt.getRank() - 1; d >= 0; --d) pos[d] = rest % pt.getDimSize(d), rest /= pt.getDimSize(d);
          for (int64_t v : pos) idx.push_back(fb.create<arith::ConstantIndexOp>(par.getLoc(), v));
          fb.create<memref::StoreOp>(par.getLoc(), fillValue, filled->priv, idx);
        }
        par.erase();
        continue;
      }
      const bool touchesOwned = ownedPhases.contains(par.getOperation());
      order(par);
      if (failed(touchesOwned ? mapOwned(par) : mapPhase(par)))
        return failure();
      continue;
    }
    // A loop of phases (the planes of a z-march). Every lane must reach the same barriers, so
    // the trip count is a constant, at least one, and the loop holds only phases and pure ops.
    auto loop = dyn_cast<scf::ForOp>(op);
    if (!loop || !llvm::any_of(loop.getBody()->getOperations(), isPhase))
      continue;
    auto lo = getConstantIntValue(loop.getLowerBound()), hi = getConstantIntValue(loop.getUpperBound()),
         st = getConstantIntValue(loop.getStep());
    if (!lo || !hi || !st || *hi <= *lo)
      return loop.emitOpError("a loop holding {mir.lanes} phases needs constant bounds and at least one trip");
    SmallVector<scf::ParallelOp> phases;
    for (Operation &inner : loop.getBody()->without_terminator()) {
      if (isPhase(inner))
        phases.push_back(cast<scf::ParallelOp>(inner));
      else if (!isMemoryEffectFree(&inner))
        return inner.emitOpError("a loop holding {mir.lanes} phases may hold only phases and pure ops");
    }
    Accesses all;
    collect(loop, all);
    for (scf::ParallelOp par : phases) {
      order(par);
      if (failed(mapPhase(par)))
        return failure();
    }
    // The next trip's phases must not touch a buffer before this trip's last ones are done with it.
    if (conflict(pending, all)) {
      OpBuilder(loop.getBody()->getTerminator()).create<gpu::BarrierOp>(loop.getLoc());
      pending.clear();
    }
  }
  WalkResult nested = elements.walk([](scf::ParallelOp p) {
    return p->hasAttr("mir.lanes") ? WalkResult::interrupt() : WalkResult::advance();
  });
  if (nested.wasInterrupted())
    return elements.emitOpError("{mir.lanes} loops must sit in the element body or in a loop of phases there");

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
  auto isPrivate = [](Value m) {
    auto a = m.getDefiningOp<memref::AllocaOp>();
    return a && !a->hasAttr("mir.scratch");
  };
  body->walk([&](memref::StoreOp st) {
    Value root = rootOf(st.getMemRef());
    if (!scratch.contains(root) && !isPrivate(root))
      stores.push_back(st);
  });
  for (memref::StoreOp st : stores) {
    if (st->getBlock() == body)
      guard(st, lane0);
    else if (valid)
      guard(st, valid);
  }

  // Lane-owned scratch now lives in the private arrays: drop the buffer and its views.
  for (auto &[root, ow] : owned) {
    SmallVector<Operation *> views;
    SmallVector<Value> work{root};
    while (!work.empty()) {
      Value v = work.pop_back_val();
      for (Operation *user : v.getUsers()) { views.push_back(user); work.push_back(user->getResult(0)); }
    }
    for (Operation *view : llvm::reverse(views)) view->erase();
    root.getDefiningOp()->erase();
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
  } else if (failed(packScratch(func, body, allocas, single ? Value() : group, G, L, padScratch, largestFirst))) {
    return failure();
  }
  if (anyScratch)   // the next element must not overwrite scratch a slower lane still reads
    OpBuilder(body->getTerminator()).create<gpu::BarrierOp>(loc);
  // A launch with no ghost buffer can only read local copies. Branch once
  // around the scheduled loop so every lane still reaches the same barriers.
  if (auto ghostArg = elements->getAttrOfType<IntegerAttr>("mir.interior_ghost_arg")) {
    OpBuilder at(loop);
    Value address = at.create<memref::ExtractAlignedPointerAsIndexOp>(loc, func.getArgument(ghostArg.getInt()));
    Value noGhost = at.create<arith::CmpIOp>(loc, arith::CmpIPredicate::eq, address,
                                             at.create<arith::ConstantIndexOp>(loc, 0));
    auto dispatch = at.create<scf::IfOp>(loc, noGhost, /*withElseRegion=*/true);
    OpBuilder inside = OpBuilder::atBlockBegin(dispatch.thenBlock());
    Value local = inside.create<arith::ConstantIntOp>(loc, 1, 1);
    Operation *specialized = inside.clone(*loop);
    specialized->walk([&](scf::IfOp branch) {
      if (branch->removeAttr("mir.local_copy")) branch.getConditionMutable().assign(local);
    });
    loop.walk([](scf::IfOp branch) { branch->removeAttr("mir.local_copy"); });
    loop->moveBefore(dispatch.elseBlock()->getTerminator());
  }
  elements.erase();
  mir::refreshViewTypes(func.getOperation());
  // The launch shape, unless it is NVIDIA's original form (one warp, one element), which
  // the loaders assume when none is published; on AMD that default would be half a wave.
  if (!single || L != wave || mir::getGpuTarget(func).isAmd()) {
    func->setAttr("mir.block_threads", b.getI64IntegerAttr(wave * W));
    func->setAttr("mir.elements_per_block", b.getI64IntegerAttr(G));
  }
  return success();
}

struct WarpPerElementPass : public PassWrapper<WarpPerElementPass, OperationPass<ModuleOp>> {
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(WarpPerElementPass)

  WarpPerElementPass() = default;
  WarpPerElementPass(const WarpPerElementPass &other) : PassWrapper(other) {}
  Option<bool> padScratch{*this, "pad-scratch",
                          llvm::cl::desc("choose each scratch buffer's padding by its bank conflicts"),
                          llvm::cl::init(true)};
  Option<bool> largestFirst{*this, "largest-first",
                            llvm::cl::desc("pack scratch largest first when that needs fewer slots"),
                            llvm::cl::init(true)};
  Option<bool> laneOwned{*this, "lane-owned",
                         llvm::cl::desc("keep scratch every phase touches lane by lane in registers "
                                        "(loops of phases are unrolled first)"),
                         llvm::cl::init(false)};
  Option<bool> fitShared{*this, "fit-shared",
                         llvm::cl::desc("halve the warps per block until the shared memory fits the "
                                        "target's limit (for a warp count that is a default, not a choice)"),
                         llvm::cl::init(false)};

  StringRef getArgument() const final { return "mir-warp-per-element"; }
  StringRef getDescription() const final {
    return "Map a {mir.elements} scf.parallel onto groups of lanes (one element "
           "each) and its {mir.lanes} phases onto the group's lanes, with barriers "
           "where phases share memory";
  }
  void getDependentDialects(DialectRegistry &registry) const override {
    registry.insert<arith::ArithDialect, gpu::GPUDialect, memref::MemRefDialect, scf::SCFDialect>();
  }

  // Static shared memory a mapping with the loop's current warps per block declares: map a
  // copy of the function (in the module, so it sees the target) and sum its buffers.
  FailureOr<int64_t> sharedBytes(func::FuncOp func, scf::ParallelOp loop) {
    OpBuilder b(func);
    b.setInsertionPointAfter(func);
    IRMapping map;
    auto copy = cast<func::FuncOp>(b.clone(*func, map));
    copy.setSymName((func.getSymName() + "__fit").str());
    auto copyLoop = cast<scf::ParallelOp>(map.lookup(loop.getOperation()));
    int64_t bytes = 0;
    const bool ok = succeeded(mapElements(copy, copyLoop, padScratch, largestFirst, laneOwned));
    for (unsigned i = 0; ok && i < copy.getNumArguments(); ++i)
      if (copy.getArgAttr(i, "mir.workgroup"))
        if (auto t = dyn_cast<MemRefType>(copy.getArgument(i).getType()))
          bytes += t.getNumElements() * t.getElementTypeBitWidth() / 8;
    copy.erase();
    if (!ok)
      return failure();
    return bytes;
  }

  // The warps per block, halved until the block's shared memory fits the target's static
  // limit, or one element per block remains (then the kernel cannot be built: an error).
  LogicalResult fitWarps(func::FuncOp func, scf::ParallelOp loop) {
    const mir::GpuTarget target = mir::getGpuTarget(func);
    auto attr = [&](StringRef name) { return loop->getAttrOfType<IntegerAttr>(name); };
    const int64_t L = attr("mir.lanes_per_element") ? attr("mir.lanes_per_element").getInt() : target.waveSize;
    int64_t W = attr("mir.warps_per_block") ? attr("mir.warps_per_block").getInt() : (L == target.waveSize ? 1 : 4);
    for (;;) {
      loop->setAttr("mir.warps_per_block", IntegerAttr::get(IntegerType::get(func.getContext(), 64), W));
      FailureOr<int64_t> bytes = sharedBytes(func, loop);
      if (failed(bytes))
        return failure();
      if (*bytes <= target.maxSharedBytes)
        return success();
      if (W == 1 || target.waveSize * W <= L)
        return func.emitOpError("shared memory ") << *bytes << " B per block exceeds the target's "
               << target.maxSharedBytes << " B even with " << W << " warp(s) per block";
      W /= 2;
    }
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
      for (scf::ParallelOp p : loops) {
        if (fitShared && failed(fitWarps(func, p)))
          return signalPassFailure();
        if (failed(mapElements(func, p, padScratch, largestFirst, laneOwned)))
          return signalPassFailure();
      }
    }
  }
};

}  // namespace

namespace mir {
void registerWarpPerElementPass() { PassRegistration<WarpPerElementPass>(); }
}  // namespace mir
