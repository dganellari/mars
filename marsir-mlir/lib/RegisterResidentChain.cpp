// --mir-chain-contracts: lower a chain of m8n8k4 f64 vector.contract ops to
// DIRECT per-lane nvgpu.mma.sync + gpu.shuffle, REGISTER-RESIDENT -- no shared
// memory, no barriers. This is the pass form of what mlir_warp_reg.py emits by
// hand: when one contract's result feeds the next, the C-fragment is repacked
// into the next contract's A or B fragment with a warp shuffle instead of a
// shared-memory round-trip.
//
// Input (register-resident intent): a gpu.func / func.func whose body is a
// straight-line chain of 8x8x8 vector.contract ops (standard A@B or transposed
// A@B^T, identified by indexing maps), where each operand is EITHER a
// transfer_read of an 8x8 memref (a leaf, read from memory as a fragment) OR
// the result of an earlier contract in the chain (an internal edge, relayed via
// shuffle). Pointwise ops (mulf/addf on the 8x8 C-fragment) between contracts
// are preserved -- they act per lane on the vector<1x2> and need no relayout.
//
// COHERENCE RULE: a C-fragment is only complete across the whole warp. It may be
// stored piecewise (each lane its two entries) only into memory the warp shares.
// Before it goes into per-thread memory -- a function-local alloc/alloca -- it is
// materialized: every lane gathers the full value. Missing this left each lane's
// private buffers holding two real entries and zeros, which the rest of the
// kernel then read as complete.
//
// m8n8k4 f64 fragment conventions (row.col), lane L, i=L/4, k=L%4, slab s:
//   A-frag        = A[i, 4s+k]           (read mem [i, 4s+k])
//   B-frag std    = B[4s+k, i]           (read mem [4s+k, i])
//   B-frag transp = read mem [i, 4s+k]   (B stored [n,k])
//   C-frag        = C[i, 2k], C[i, 2k+1] (vector<1x2>, read [i, 2k])
// Relayouts (register-only, gpu.shuffle idx width 32):
//   C->A slab s: src=(L&28)|(2s)|((L&3)>>1),          elem=L&1
//   C->B std   s: src=16s+4*(L&3)+(L>>3),             elem=(L>>2)&1

#include "mir/MirPasses.h"

#include "mlir/Dialect/Arith/IR/Arith.h"
#include "mlir/Dialect/Func/IR/FuncOps.h"
#include "mlir/Dialect/GPU/IR/GPUDialect.h"
#include "mlir/Dialect/MemRef/IR/MemRef.h"
#include "mlir/Interfaces/FunctionInterfaces.h"
#include "mlir/Dialect/NVGPU/IR/NVGPUDialect.h"
#include "mlir/Dialect/Vector/IR/VectorOps.h"
#include "mlir/IR/AffineMap.h"
#include "mlir/IR/Builders.h"
#include "mlir/Pass/Pass.h"

using namespace mlir;

namespace {

// Per-lane index + constant bundle, materialized once at the top of the kernel.
struct Lane {
  Value lane, i, k, li;                 // lane id, L/4, L%4 (index), lane as i32
  Value c2idx, c4idx;                   // 2, 4 : index
  Value f0;                             // 0.0 : f64
  VectorType frag1, frag2;              // vector<1x1>, vector<1x2>
  Value zc1;                            // dense<0.0> : vector<1x1>
  Type i32;
  Value ci(OpBuilder &b, Location loc, int v) {
    return b.create<arith::ConstantIntOp>(loc, v, 32);
  }
  Value idx(OpBuilder &b, Location loc, int v) {
    return b.create<arith::ConstantIndexOp>(loc, v);
  }
};

static Lane makeLane(OpBuilder &b, Location loc) {
  Lane L;
  auto f64 = b.getF64Type();
  L.i32 = b.getIntegerType(32);
  L.frag1 = VectorType::get({1, 1}, f64);
  L.frag2 = VectorType::get({1, 2}, f64);
  L.f0 = b.create<arith::ConstantOp>(loc, f64, b.getF64FloatAttr(0.0));
  L.zc1 = b.create<arith::ConstantOp>(
      loc, L.frag1, DenseElementsAttr::get(L.frag1, b.getF64FloatAttr(0.0)));
  L.c2idx = b.create<arith::ConstantIndexOp>(loc, 2);
  L.c4idx = b.create<arith::ConstantIndexOp>(loc, 4);
  L.lane = b.create<gpu::ThreadIdOp>(loc, gpu::Dimension::x);
  L.i = b.create<arith::DivUIOp>(loc, L.lane, L.c4idx);
  L.k = b.create<arith::RemUIOp>(loc, L.lane, L.c4idx);
  L.li = b.create<arith::IndexCastOp>(loc, L.i32, L.lane);
  return L;
}

// col = 4*slab + k (index).
static Value colOf(OpBuilder &b, Location loc, Lane &L, int slab) {
  if (slab == 0)
    return L.k;
  Value c = b.create<arith::ConstantIndexOp>(loc, 4 * slab);
  return b.create<arith::AddIOp>(loc, L.k, c);
}

static Value readFrag1(OpBuilder &b, Location loc, Lane &L, Value mem, Value row,
                       Value col, bool rowInBounds = true) {
  return b.create<vector::TransferReadOp>(
      loc, L.frag1, mem, ValueRange{row, col},
      AffineMap::getMinorIdentityMap(2, 2, b.getContext()), L.f0, Value(),
      b.getBoolArrayAttr({rowInBounds, true}));
}

static Value readCTile(OpBuilder &b, Location loc, Lane &L, Value mem, int tile,
                       ValueRange base = {}, bool rowInBounds = true) {
  Value col2k = b.create<arith::MulIOp>(loc, L.k, L.c2idx);
  if (tile) {
    Value off = b.create<arith::ConstantIndexOp>(loc, 8 * tile);
    col2k = b.create<arith::AddIOp>(loc, col2k, off);
  }
  Value row = L.i;
  if (base.size() == 2) {   // the read's own window origin
    row = b.create<arith::AddIOp>(loc, base[0], row);
    col2k = b.create<arith::AddIOp>(loc, base[1], col2k);
  }
  return b.create<vector::TransferReadOp>(
      loc, L.frag2, mem, ValueRange{row, col2k},
      AffineMap::getMinorIdentityMap(2, 2, b.getContext()), L.f0, Value(),
      b.getBoolArrayAttr({rowInBounds, true}));
}

static Value mma(OpBuilder &b, Location loc, Lane &L, Value a, Value bfrag,
                 Value acc) {
  return b.create<nvgpu::MmaSyncOp>(loc, a, bfrag, acc,
                                    b.getDenseI64ArrayAttr({8, 8, 4}));
}

// Repack a C-fragment (cfrag, vector<1x2>) into an A- or B-fragment for k-slab
// `slab` via gpu.shuffle. `toB` selects the C->B(std) formula; else C->A.
static Value relayout(OpBuilder &b, Location loc, Lane &L, Value cfrag, int slab,
                      bool toB) {
  Value lo = b.create<vector::ExtractOp>(loc, cfrag, ArrayRef<int64_t>{0, 0});
  Value hi = b.create<vector::ExtractOp>(loc, cfrag, ArrayRef<int64_t>{0, 1});
  Value li = L.li, src, elemOdd;
  auto ci = [&](int v) { return L.ci(b, loc, v); };
  if (!toB) {
    // C->A: src=(L&28)|(2s)|((L&3)>>1), elem=L&1
    Value m28 = b.create<arith::AndIOp>(loc, li, ci(28));
    Value m3 = b.create<arith::AndIOp>(loc, li, ci(3));
    Value sh = b.create<arith::ShRUIOp>(loc, m3, ci(1));
    src = b.create<arith::OrIOp>(loc, m28, sh);
    if (slab)
      src = b.create<arith::OrIOp>(loc, src, ci(2 * slab));
    Value par = b.create<arith::AndIOp>(loc, li, ci(1));
    elemOdd = b.create<arith::CmpIOp>(loc, arith::CmpIPredicate::eq, par, ci(1));
  } else {
    // C->B std: src=16s+4*(L&3)+(L>>3), elem=(L>>2)&1
    Value m3 = b.create<arith::AndIOp>(loc, li, ci(3));
    Value q4 = b.create<arith::MulIOp>(loc, m3, ci(4));
    Value hi3 = b.create<arith::ShRUIOp>(loc, li, ci(3));
    src = b.create<arith::AddIOp>(loc, q4, hi3);
    if (slab)
      src = b.create<arith::AddIOp>(loc, src, ci(16 * slab));
    Value sh2 = b.create<arith::ShRUIOp>(loc, li, ci(2));
    Value par = b.create<arith::AndIOp>(loc, sh2, ci(1));
    elemOdd = b.create<arith::CmpIOp>(loc, arith::CmpIPredicate::eq, par, ci(1));
  }
  Value width = ci(32);
  auto s0 = b.create<gpu::ShuffleOp>(loc, lo, src, width, gpu::ShuffleMode::IDX);
  auto s1 = b.create<gpu::ShuffleOp>(loc, hi, src, width, gpu::ShuffleMode::IDX);
  Value sel = b.create<arith::SelectOp>(loc, elemOdd, s1.getShuffleResult(),
                                        s0.getShuffleResult());
  return b.create<vector::InsertOp>(loc, sel, L.zc1, ArrayRef<int64_t>{0, 0});
}

// Is `mem` private to each thread? A function-local buffer is: memref.alloc
// becomes a per-thread malloc on the GPU and memref.alloca per-thread stack.
// Kernel arguments are global memory the whole warp shares. Anything this walk
// cannot trace to a kernel argument is treated as private -- that only costs a
// gather, while the opposite mistake silently loses data.
static bool isThreadPrivate(Value mem) {
  for (int guard = 0; guard < 64; ++guard) {
    if (auto ba = dyn_cast<BlockArgument>(mem)) {
      Operation *owner = ba.getOwner()->getParentOp();
      return !(owner && isa<FunctionOpInterface>(owner) &&
               ba.getOwner()->isEntryBlock());
    }
    Operation *d = mem.getDefiningOp();
    if (!d || isa<memref::AllocOp, memref::AllocaOp>(d))
      return true;
    if (auto sv = dyn_cast<memref::SubViewOp>(d)) { mem = sv.getSource(); continue; }
    if (auto cs = dyn_cast<memref::CollapseShapeOp>(d)) { mem = cs.getSrc(); continue; }
    if (auto es = dyn_cast<memref::ExpandShapeOp>(d)) { mem = es.getSrc(); continue; }
    if (auto rc = dyn_cast<memref::ReinterpretCastOp>(d)) { mem = rc.getSource(); continue; }
    if (auto ca = dyn_cast<memref::CastOp>(d)) { mem = ca.getSource(); continue; }
    return true;
  }
  return true;
}

// Rebuild the full value of `ty` (R x 8*nTiles) in EVERY lane from its C-fragments.
// Element (r, 8t + c) lives in lane 4r + c/2 as component c%2 of tile t, so each
// element is one broadcast shuffle. Needed before a fragment is stored into
// per-thread memory: there each lane has its own copy of the buffer, and a lane
// that stored only its two entries would leave the rest of its copy stale.
static Value materialize(OpBuilder &b, Location loc, Lane &L,
                         ArrayRef<Value> tiles, VectorType ty) {
  const int64_t R = ty.getDimSize(0), C = ty.getDimSize(1);
  Value width = L.ci(b, loc, 32);
  SmallVector<Value> lo, hi;
  for (Value t : tiles) {
    lo.push_back(b.create<vector::ExtractOp>(loc, t, ArrayRef<int64_t>{0, 0}));
    hi.push_back(b.create<vector::ExtractOp>(loc, t, ArrayRef<int64_t>{0, 1}));
  }
  Value acc = b.create<arith::ConstantOp>(
      loc, ty, DenseElementsAttr::get(ty, b.getF64FloatAttr(0.0)));
  for (int64_t r = 0; r < R; ++r)
    for (int64_t col = 0; col < C; ++col) {
      const int t = (int)(col / 8), c = (int)(col % 8);
      Value src = L.ci(b, loc, (int)(4 * r + c / 2));
      Value v = (c % 2) ? hi[t] : lo[t];
      Value e = b.create<gpu::ShuffleOp>(loc, v, src, width, gpu::ShuffleMode::IDX)
                    .getShuffleResult();
      acc = b.create<vector::InsertOp>(loc, e, acc, ArrayRef<int64_t>{r, col});
    }
  return acc;
}

// A transfer the fragment path can address directly: rank-2 memref, identity
// permutation map, no mask. Anything else is declined (reads) or materialized
// and handed to the original op (writes), which keeps its own map and indices.
static bool isPlain2D(Operation *op) {
  if (auto r = dyn_cast<vector::TransferReadOp>(op))
    return r.getShapedType().getRank() == 2 && r.getPermutationMap().isIdentity() &&
           !r.getMask();
  if (auto w = dyn_cast<vector::TransferWriteOp>(op))
    return w.getShapedType().getRank() == 2 && w.getPermutationMap().isIdentity() &&
           !w.getMask();
  return false;
}

// m8n8k4 iteration-space maps (m=d0, n=d1, k=d2).
struct Maps {
  AffineMap mk, kn, nk, mn;
  Maps(MLIRContext *ctx) {
    auto d0 = getAffineDimExpr(0, ctx), d1 = getAffineDimExpr(1, ctx),
         d2 = getAffineDimExpr(2, ctx);
    mk = AffineMap::get(3, 0, {d0, d2}, ctx);
    kn = AffineMap::get(3, 0, {d2, d1}, ctx);
    nk = AffineMap::get(3, 0, {d1, d2}, ctx);
    mn = AffineMap::get(3, 0, {d0, d1}, ctx);
  }
};

// Classify a vector.contract: is it m8n8k4 f64, and which operand is A vs B,
// and is B transposed? Returns false if it is not a supported contract.
static bool classify(vector::ContractionOp c, Maps &M, Value &A, Value &B,
                     bool &bTransp, int64_t &K, int64_t &N, int64_t &M_) {
  auto maps = c.getIndexingMapsArray();
  if (maps.size() != 3 || maps[2] != M.mn)
    return false;
  auto iters = c.getIteratorTypesArray();
  if (iters.size() != 3 || iters[0] != vector::IteratorType::parallel ||
      iters[1] != vector::IteratorType::parallel ||
      iters[2] != vector::IteratorType::reduction)
    return false;
  Value lhs = c.getLhs(), rhs = c.getRhs();
  if (maps[0] == M.mk && (maps[1] == M.kn || maps[1] == M.nk)) {
    A = lhs; B = rhs; bTransp = (maps[1] == M.nk);
  } else if (maps[1] == M.mk && (maps[0] == M.kn || maps[0] == M.nk)) {
    A = rhs; B = lhs; bTransp = (maps[0] == M.nk);
  } else {
    return false;
  }
  auto shp = [](Type t) { return cast<VectorType>(t).getShape(); };
  if (!cast<VectorType>(A.getType()).getElementType().isF64())
    return false;
  auto as = shp(A.getType()), bs = shp(B.getType()),
       cs = shp(c.getResultType());
  if (as.size() != 2 || bs.size() != 2 || cs.size() != 2)
    return false;
  // m8n8k4 fixes the TILE at m = 8, but the operator's m may be smaller -- the
  // Knaus B-sweep is Pxn with P = 7 faces. A shorter m rides in the same tile:
  // the tail rows read out of bounds (transfer_read pads them with 0) and their
  // writes are dropped, so the result is exact with one wasted row.
  if (as[0] != cs[0] || as[0] < 1 || as[0] > 8)
    return false;
  M_ = as[0];
  K = as[1];
  N = cs[1];
  if (K % 4 != 0 || N % 8 != 0)
    return false;
  if (bTransp)
    return bs[0] == N && bs[1] == K;
  return bs[0] == K && bs[1] == N;
}

struct ChainContractsPass
    : public PassWrapper<ChainContractsPass, OperationPass<>> {
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(ChainContractsPass)

  ChainContractsPass() = default;
  ChainContractsPass(const ChainContractsPass &other) : PassWrapper(other) {}

  // Feature switches, all on by default. Turning one off makes the pass DECLINE
  // what that feature would have lowered, so those ops fall back to ordinary
  // lowering. Used to bisect a numerical failure to one feature.
  Option<bool> fusePointwise{*this, "fuse-pointwise",
      llvm::cl::desc("lower elementwise ops onto C-fragments"), llvm::cl::init(true)};
  Option<bool> shortM{*this, "short-m",
      llvm::cl::desc("lower contractions with m < 8 via out-of-bounds padding"),
      llvm::cl::init(true)};
  Option<bool> wideN{*this, "wide-n",
      llvm::cl::desc("lower contractions with N > 8 as several column tiles"),
      llvm::cl::init(true)};
  Option<bool> memAcc{*this, "mem-acc",
      llvm::cl::desc("accept an accumulator read from memory"), llvm::cl::init(true)};

  StringRef getArgument() const final { return "mir-chain-contracts"; }
  StringRef getDescription() const final {
    return "Lower a chain of m8n8k4 vector.contract to register-resident "
           "nvgpu.mma.sync + gpu.shuffle (no shared memory)";
  }
  void getDependentDialects(DialectRegistry &registry) const override {
    registry.insert<gpu::GPUDialect, nvgpu::NVGPUDialect, vector::VectorDialect,
                    arith::ArithDialect>();
  }

  void runOnOperation() override {
    Operation *root = getOperation();
    MLIRContext *ctx = &getContext();
    Maps M(ctx);

    // EVERY block holding contracts, not just the first: a real operator spreads
    // them over several scf.for bodies, and taking only the first silently left
    // the rest unlowered.
    SmallVector<Block *> blocks;
    DenseSet<Block *> seenBlocks;
    root->walk([&](vector::ContractionOp c) {
      if (seenBlocks.insert(c->getBlock()).second)
        blocks.push_back(c->getBlock());
    });
    if (blocks.empty())
      return;

    // Lane values go at the enclosing function's entry block so they dominate
    // every block below; putting them in a loop body would not dominate a
    // sibling loop.
    Operation *fnOp = blocks.front()->getParentOp();
    while (fnOp && !isa<FunctionOpInterface>(fnOp))
      fnOp = fnOp->getParentOp();
    Block *entry = (fnOp && fnOp->getNumRegions() && !fnOp->getRegion(0).empty())
                       ? &fnOp->getRegion(0).front()
                       : blocks.front();

    OpBuilder b(ctx);
    b.setInsertionPointToStart(entry);
    Lane L = makeLane(b, root->getLoc());

    // A value of type vector<8xN> lives as N/8 per-lane C-fragments, one per
    // column tile, each vector<1x2>. Contract results and pointwise results are
    // BOTH in C layout, so a pointwise op reads its operands from the same map a
    // following contract relayouts out of.
    DenseMap<Value, SmallVector<Value>> frag;
    // How many of the tile's 8 rows a lowered value actually occupies, so its
    // write knows whether to clip the tail.
    DenseMap<Value, int64_t> fragRows;
    SmallVector<Operation *> dead;

    // Column offset of tile t, as an index value.
    auto tileCol = [&](Location loc, int t) -> Value {
      return b.create<arith::ConstantIndexOp>(loc, 8 * t);
    };
    auto addCol = [&](Location loc, Value base, int t) -> Value {
      return t ? b.create<arith::AddIOp>(loc, base, tileCol(loc, t)).getResult()
               : base;
    };

    // A is 8xK: the same fragment feeds every column tile.
    auto operandFragA = [&](Value v, int slab, bool rowIB) -> Value {
      Location loc = v.getLoc();
      if (auto it = frag.find(v); it != frag.end()) {
        if (it->second.size() != 1)
          return Value();   // a multi-tile value cannot be an A operand
        return relayout(b, loc, L, it->second[0], slab, /*toB=*/false);
      }
      auto rd = v.getDefiningOp<vector::TransferReadOp>();
      if (!rd || !isPlain2D(rd))
        return Value();  // neither a plain leaf read nor a value we lowered
      Value col = colOf(b, loc, L, slab);
      Value r0 = b.create<arith::AddIOp>(loc, rd.getIndices()[0], L.i);
      Value c0 = b.create<arith::AddIOp>(loc, rd.getIndices()[1], col);
      return readFrag1(b, loc, L, rd.getSource(), r0, c0, rowIB);  // A[i,4s+k]
    };
    // B supplies columns 8t..8t+7 for tile t.
    auto operandFragB = [&](Value v, int slab, bool transp, int tile) -> Value {
      Location loc = v.getLoc();
      if (auto it = frag.find(v); it != frag.end()) {
        if (it->second.size() != 1)
          return Value();
        // A TRANSPOSED B operand X (indexed [n,k]) supplies B[k][n] = X[n][k] to
        // lane L as X[L/4][4s + L%4] -- which is X's A-fragment, not its
        // B-fragment. Only the standard (k,n) form takes the C->B relayout.
        return relayout(b, loc, L, it->second[0], slab, /*toB=*/!transp);
      }
      auto rd = v.getDefiningOp<vector::TransferReadOp>();
      if (!rd || !isPlain2D(rd))
        return Value();
      Value col = colOf(b, loc, L, slab);
      Value nIdx = addCol(loc, L.i, tile);
      Value b0 = rd.getIndices()[0], b1 = rd.getIndices()[1];
      if (transp)   // B stored [n,k]: read [8t+i, 4s+k]
        return readFrag1(b, loc, L, rd.getSource(),
                         b.create<arith::AddIOp>(loc, b0, nIdx),
                         b.create<arith::AddIOp>(loc, b1, col));
      // B stored [k,n]: read [4s+k, 8t+i]
      return readFrag1(b, loc, L, rd.getSource(),
                       b.create<arith::AddIOp>(loc, b0, col),
                       b.create<arith::AddIOp>(loc, b1, nIdx));
    };

    // Shape of a value in tiles, or 0 if it is not an f64 8xN vector.
    auto tilesOf = [](Value v) -> int {
      auto t = dyn_cast<VectorType>(v.getType());
      if (!t || !t.getElementType().isF64() || t.getRank() != 2) return 0;
      if (t.getDimSize(0) != 8 || t.getDimSize(1) % 8 != 0) return 0;
      return (int)(t.getDimSize(1) / 8);
    };

    // Fragment tile `tile` of a POINTWISE operand, which needs C layout.
    // Already lowered -> reuse. Leaf transfer_read -> read [i, 8t + 2k]. Splat
    // constant -> a vector<1x2> splat (a uniform value has no layout).
    auto operandFragC = [&](Value v, int tile) -> Value {
      if (auto it = frag.find(v); it != frag.end())
        return tile < (int)it->second.size() ? it->second[tile] : Value();
      Location loc = v.getLoc();
      if (auto rd = v.getDefiningOp<vector::TransferReadOp>())
        return readCTile(b, loc, L, rd.getSource(), tile, rd.getIndices());
      if (auto cst = v.getDefiningOp<arith::ConstantOp>())
        if (auto d = dyn_cast<DenseElementsAttr>(cst.getValue()))
          if (d.isSplat())
            return b.create<arith::ConstantOp>(
                loc, L.frag2,
                DenseElementsAttr::get(L.frag2, d.getSplatValue<APFloat>()));
      return Value();
    };

    // ---- Phase 1: the LOWERABLE CLOSURE ----------------------------------
    // An op can be lowered only if every operand it needs is available as a
    // fragment AND every consumer of its result can itself consume a fragment.
    // Both directions matter: a fragment cannot be materialized back into a full
    // vector, so lowering a value whose consumer we cannot handle would strand
    // it. Neither condition is local, hence the fixpoint.
    auto isLeafRead = [](Value v) {
      auto rd = v.getDefiningOp<vector::TransferReadOp>();
      return rd && isPlain2D(rd);
    };
    auto isSplatCst = [](Value v) {
      if (auto cst = v.getDefiningOp<arith::ConstantOp>())
        if (auto dv = dyn_cast<DenseElementsAttr>(cst.getValue()))
          return dv.isSplat();
      return false;
    };

    SmallVector<Operation *> cand;
    DenseSet<Operation *> inCand;
    for (Block *blk : blocks)
      for (Operation &o : *blk) {
        Operation *op = &o;
        bool ok = false;
        if (auto c = dyn_cast<vector::ContractionOp>(op)) {
          Value A2, B2; bool bt2; int64_t K2, N2, m2;
          ok = classify(c, M, A2, B2, bt2, K2, N2, m2) &&
               (shortM || m2 == 8) && (wideN || N2 == 8);
        } else if (fusePointwise && op->hasTrait<OpTrait::Elementwise>() &&
                   op->getNumResults() == 1 && tilesOf(op->getResult(0)) > 0) {
          ok = llvm::all_of(op->getOperands(), [&](Value v) {
            return tilesOf(v) == tilesOf(op->getResult(0));
          });
        }
        if (ok) { cand.push_back(op); inCand.insert(op); }
      }

    for (bool changed = true; changed;) {
      changed = false;
      for (Operation *op : cand) {
        if (!inCand.count(op))
          continue;
        // A contract's A/B must be a leaf read or a SINGLE-tile lowered value:
        // the relayout formulas are per-tile, so a wide value cannot be an
        // operand. This has to mirror the emit-time check exactly, or the
        // closure promises something emit then declines and the chain breaks.
        auto operandOk = [&](Value v, bool singleTileOnly) {
          if (isLeafRead(v) || isSplatCst(v))
            return true;
          Operation *d = v.getDefiningOp();
          if (!d || !inCand.count(d))
            return false;
          return !singleTileOnly || tilesOf(v) == 1;
        };
        bool good;
        if (auto c = dyn_cast<vector::ContractionOp>(op)) {
          Value A2, B2; bool bt2; int64_t K2, N2, m2;
          good = classify(c, M, A2, B2, bt2, K2, N2, m2) &&
                 operandOk(A2, /*singleTileOnly=*/true) &&
                 operandOk(B2, /*singleTileOnly=*/true) &&
                 operandOk(c.getAcc(), /*singleTileOnly=*/false) &&
                 (memAcc || !isLeafRead(c.getAcc()));
        } else {
          good = llvm::all_of(op->getOperands(), [&](Value v) {
            return operandOk(v, /*singleTileOnly=*/false);
          });
        }
        if (good)
          for (Operation *u : op->getResult(0).getUsers())
            if (!inCand.count(u) && !isa<vector::TransferWriteOp>(u)) {
              good = false;
              break;
            }
        if (!good) { inCand.erase(op); changed = true; }
      }
    }

    // One in-order pass per block: within a block the chain is straight-line, so
    // lowering each op as it is reached keeps `frag` populated before any
    // consumer needs it.
    for (Block *body : blocks)
    for (Operation &opRef : *body) {
      Operation *op = &opRef;

      if (auto c = dyn_cast<vector::ContractionOp>(op)) {
        if (!inCand.count(op))
          continue;   // outside the lowerable closure
        Value A, B;
        bool bt;
        int64_t K = 0, N = 0, mDim = 8;
        // DECLINE rather than fail: a real operator mixes shapes, and m is fixed
        // at 8 by the hardware tile. A contraction that does not fit (the Knaus
        // B-sweep is PxN with P = 7 faces) is left alone for another lowering,
        // not erased and not silently mangled.
        if (!classify(c, M, A, B, bt, K, N, mDim))
          continue;
        const bool rowIB = (mDim == 8);   // else the tail rows ride OOB
        const int nTiles = (int)(N / 8), nSlabs = (int)(K / 4);

        // Check every operand BEFORE emitting anything, so declining leaves no
        // half-lowered contract behind. A multi-tile value cannot be an A or B
        // operand (its layout is per-tile, the relayout formulas are not).
        auto usable = [&](Value v) {
          auto it = frag.find(v);
          if (it != frag.end()) return it->second.size() == 1;
          return (bool)v.getDefiningOp<vector::TransferReadOp>();
        };
        if (!usable(A) || !usable(B))
          continue;
        {
          Value av = c.getAcc();
          bool accOk = frag.count(av) || av.getDefiningOp<vector::TransferReadOp>();
          if (!accOk)
            if (auto cst = av.getDefiningOp<arith::ConstantOp>())
              if (auto dv = dyn_cast<DenseElementsAttr>(cst.getValue()))
                accOk = dv.isSplat() && dv.getSplatValue<APFloat>().isZero();
          if (!accOk) {
            c.emitOpError("mir-chain-contracts: accumulator is neither a zero "
                          "splat, a memory read, nor a value this pass lowered, "
                          "so it would be dropped");
            signalPassFailure();
            return;
          }
        }
        b.setInsertionPoint(c);

        // The accumulator decides where each tile's chain starts. A zero splat
        // starts a fresh one; a value this pass already lowered continues one,
        // in C layout. Anything else would be SILENTLY DROPPED -- refuse it.
        SmallVector<Value> accFrags;
        if (auto it = frag.find(c.getAcc()); it != frag.end()) {
          if ((int)it->second.size() != nTiles) {
            c.emitOpError("mir-chain-contracts: accumulator tile count does not "
                          "match the result");
            signalPassFailure();
            return;
          }
          accFrags = it->second;
        } else if (auto accRd =
                       c.getAcc().getDefiningOp<vector::TransferReadOp>()) {
          // An accumulator staged in memory (what tiling a matmul produces):
          // take it in C layout and accumulate straight into it.
          for (int t = 0; t < nTiles; ++t)
            accFrags.push_back(readCTile(b, c.getLoc(), L, accRd.getSource(), t,
                                         accRd.getIndices(), rowIB));
        } else {
          bool zeroAcc = false;
          if (auto cst = c.getAcc().getDefiningOp<arith::ConstantOp>())
            if (auto dv = dyn_cast<DenseElementsAttr>(cst.getValue()))
              zeroAcc = dv.isSplat() && dv.getSplatValue<APFloat>().isZero();
          if (!zeroAcc) {
            c.emitOpError("mir-chain-contracts: accumulator is neither a zero "
                          "splat nor a value this pass lowered, so it would be "
                          "dropped");
            signalPassFailure();
            return;
          }
          Value z = b.create<arith::ConstantOp>(
              c.getLoc(), L.frag2,
              DenseElementsAttr::get(L.frag2, b.getF64FloatAttr(0.0)));
          accFrags.assign(nTiles, z);
        }

        SmallVector<Value> out;
        for (int t = 0; t < nTiles; ++t) {
          Value cfrag = accFrags[t];
          for (int s2 = 0; s2 < nSlabs; ++s2) {
            Value af = operandFragA(A, s2, rowIB);
            Value bf = operandFragB(B, s2, bt, t);
            if (!af || !bf) {   // pre-checked above; defensive
              c.emitOpError("mir-chain-contracts: operand became unfragmentable");
              signalPassFailure();
              return;
            }
            cfrag = mma(b, c.getLoc(), L, af, bf, cfrag);
          }
          out.push_back(cfrag);
        }
        frag[c.getResult()] = out;
        fragRows[c.getResult()] = mDim;
        dead.push_back(op);
        continue;
      }

      // Pointwise (the flux): elementwise on 8x8 needs NO relayout -- it acts
      // lane-locally on the vector<1x2>. Only fuse ops that actually touch the
      // chain; an elementwise op on two leaf reads is left alone.
      if (op->hasTrait<OpTrait::Elementwise>() && inCand.count(op) &&
          op->getNumResults() == 1 && tilesOf(op->getResult(0)) > 0 &&
          llvm::all_of(op->getOperands(),
                       [&](Value v) {
                         return tilesOf(v) == tilesOf(op->getResult(0));
                       })) {
        // No "must already touch the chain" test here: the closure decided, and
        // an op it kept may legitimately have only leaf-read operands (a flux
        // term built from two metric reads) while its CONSUMER is on the chain.
        // Skipping it would strand that consumer with one fragment operand and
        // one full vector.
        b.setInsertionPoint(op);
        const int nTiles = tilesOf(op->getResult(0));
        SmallVector<Value> out;
        for (int t = 0; t < nTiles; ++t) {
          SmallVector<Value> fops;
          for (Value v : op->getOperands()) {
            Value f = operandFragC(v, t);
            if (!f) {
              op->emitOpError("mir-chain-contracts: pointwise operand is not a "
                              "fragment, a leaf transfer_read or a splat");
              signalPassFailure();
              return;
            }
            fops.push_back(f);
          }
          OperationState st(op->getLoc(), op->getName());
          st.addOperands(fops);
          st.addTypes({L.frag2});
          st.addAttributes(op->getAttrs());
          out.push_back(b.create(st)->getResult(0));
        }
        frag[op->getResult(0)] = out;
        {   // a pointwise result is as tall as its tallest fragmented operand
          int64_t rows = 8;
          for (Value v : op->getOperands())
            if (auto it = fragRows.find(v); it != fragRows.end()) rows = it->second;
          fragRows[op->getResult(0)] = rows;
        }
        dead.push_back(op);
        continue;
      }

      // A write of a lowered value. Into memory the warp SHARES (a kernel
      // argument), each lane stores its own two entries: tile t's C-fragment at
      // [i, 8t + 2k]. Into PER-THREAD memory that would leave every lane's copy
      // missing the other lanes' entries, so the full value is gathered first and
      // handed to the original write, which keeps its own map and indices.
      if (auto w = dyn_cast<vector::TransferWriteOp>(op)) {
        auto it = frag.find(w.getVector());
        if (it == frag.end())
          continue;
        b.setInsertionPoint(w);
        if (isThreadPrivate(w.getSource()) || !isPlain2D(w)) {
          Value full = materialize(b, w.getLoc(), L, it->second,
                                   cast<VectorType>(w.getVector().getType()));
          w->setOperand(0, full);
          continue;   // the write stays; it now stores the gathered value
        }
        int64_t wRows = 8;
        if (auto it = fragRows.find(w.getVector()); it != fragRows.end())
          wRows = it->second;
        Value col2k = b.create<arith::MulIOp>(w.getLoc(), L.k, L.c2idx);
        // The write's own indices are the base of the destination window.
        SmallVector<Value> base(w.getIndices().begin(), w.getIndices().end());
        Value row = base.size() == 2 ? b.create<arith::AddIOp>(w.getLoc(), base[0], L.i)
                                     : L.i;
        Value colBase = base.size() == 2 ? base[1] : Value();
        for (int t = 0; t < (int)it->second.size(); ++t) {
          Value col = col2k;
          if (t)
            col = b.create<arith::AddIOp>(
                w.getLoc(), col,
                b.create<arith::ConstantIndexOp>(w.getLoc(), 8 * t));
          if (colBase)
            col = b.create<arith::AddIOp>(w.getLoc(), colBase, col);
          b.create<vector::TransferWriteOp>(
              w.getLoc(), it->second[t], w.getSource(), ValueRange{row, col},
              AffineMapAttr::get(AffineMap::getMinorIdentityMap(2, 2, ctx)),
              /*mask=*/Value(),
              b.getBoolArrayAttr({wRows == 8, true}));
        }
        dead.push_back(op);
        continue;
      }
    }

    // Erasing a lowered op whose result is still consumed by something we did
    // not rewrite leaves a dangling operand -- that used to crash the verifier
    // rather than report anything. Refuse instead.
    DenseSet<Operation *> deadSet(dead.begin(), dead.end());
    for (Operation *op : dead)
      for (Value r : op->getResults())
        for (Operation *u : r.getUsers())
          if (!deadSet.count(u)) {
            u->emitOpError("mir-chain-contracts: unhandled consumer of a "
                           "register-resident value (only contracts, "
                           "elementwise 8x8 ops and transfer_write are lowered)");
            signalPassFailure();
            return;
          }
    for (Operation *op : llvm::reverse(dead))
      op->erase();
  }
};

}  // namespace

namespace mir {
void registerChainContractsPass() { PassRegistration<ChainContractsPass>(); }
}  // namespace mir
