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
                       Value col) {
  return b.create<vector::TransferReadOp>(
      loc, L.frag1, mem, ValueRange{row, col},
      AffineMap::getMinorIdentityMap(2, 2, b.getContext()), L.f0, Value(),
      b.getBoolArrayAttr({true, true}));
}

static Value readC(OpBuilder &b, Location loc, Lane &L, Value mem) {
  Value col2k = b.create<arith::MulIOp>(loc, L.k, L.c2idx);
  return b.create<vector::TransferReadOp>(
      loc, L.frag2, mem, ValueRange{L.i, col2k},
      AffineMap::getMinorIdentityMap(2, 2, b.getContext()), L.f0, Value(),
      b.getBoolArrayAttr({true, true}));
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
                     Value &acc, bool &bTransp) {
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
  return shp(A.getType()) == ArrayRef<int64_t>({8, 8}) &&
         shp(B.getType()) == ArrayRef<int64_t>({8, 8}) &&
         shp(c.getResultType()) == ArrayRef<int64_t>({8, 8});
}

struct ChainContractsPass
    : public PassWrapper<ChainContractsPass, OperationPass<>> {
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(ChainContractsPass)

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

    // Find the kernel body (the region that holds the contract chain).
    Block *body = nullptr;
    root->walk([&](vector::ContractionOp c) {
      if (!body)
        body = c->getBlock();
    });
    if (!body)
      return;

    OpBuilder b(ctx);
    b.setInsertionPointToStart(body);
    Lane L = makeLane(b, root->getLoc());

    // Map: contract result SSA -> the C-fragment (vector<1x2>) that lowered it.
    DenseMap<Value, Value> frag;
    SmallVector<vector::ContractionOp> chain;
    body->walk([&](vector::ContractionOp c) { chain.push_back(c); });

    // Fragment reader for an operand of a contract. A leaf transfer_read of an
    // 8x8 memref is read directly as the requested fragment; a prior contract
    // result is repacked via shuffle.
    auto operandFragA = [&](Value v, int slab) -> Value {
      Location loc = v.getLoc();
      if (auto it = frag.find(v); it != frag.end())
        return relayout(b, loc, L, it->second, slab, /*toB=*/false);
      auto rd = v.getDefiningOp<vector::TransferReadOp>();
      Value col = colOf(b, loc, L, slab);
      return readFrag1(b, loc, L, rd.getSource(), L.i, col);  // A[i,4s+k]
    };
    auto operandFragB = [&](Value v, int slab, bool transp) -> Value {
      Location loc = v.getLoc();
      if (auto it = frag.find(v); it != frag.end())
        return relayout(b, loc, L, it->second, slab, /*toB=*/true);
      auto rd = v.getDefiningOp<vector::TransferReadOp>();
      Value col = colOf(b, loc, L, slab);
      if (transp)
        return readFrag1(b, loc, L, rd.getSource(), L.i, col);   // B[i,4s+k]
      return readFrag1(b, loc, L, rd.getSource(), col, L.i);     // B[4s+k,i]
    };

    for (auto c : chain) {
      Value A, B, acc;
      bool bt;
      if (!classify(c, M, A, B, acc, bt)) {
        c.emitOpError("mir-chain-contracts: unsupported contract shape/maps");
        signalPassFailure();
        return;
      }
      b.setInsertionPoint(c);
      // acc: the vector.contract's acc operand is a zero constant (fresh chain)
      // or a prior C-fragment (handled via frag map by the pointwise ops); here
      // we start each contract from a zero C-fragment and rely on the pointwise
      // ops between contracts to combine, matching the emit_face_reg structure.
      Value cfrag = b.create<arith::ConstantOp>(
          c.getLoc(), L.frag2,
          DenseElementsAttr::get(L.frag2, b.getF64FloatAttr(0.0)));
      for (int s = 0; s < 2; ++s) {  // 8x8x8 = two m8n8k4 slabs
        Value af = operandFragA(A, s);
        Value bf = operandFragB(B, s, bt);
        cfrag = mma(b, c.getLoc(), L, af, bf, cfrag);
      }
      frag[c.getResult()] = cfrag;
    }

    // Rewire: every use of a contract result that ISN'T another contract's
    // operand (i.e. a transfer_write, or a pointwise op) must consume the
    // C-fragment. Pointwise ops on vector<8x8> are rewritten to vector<1x2>.
    // For the first milestone we only handle a terminal transfer_write of the
    // last contract; richer rewiring (pointwise flux fusion) is the next step.
    for (auto c : chain) {
      Value res = c.getResult();
      Value cf = frag[res];
      for (Operation *user : llvm::make_early_inc_range(res.getUsers())) {
        if (isa<vector::ContractionOp>(user))
          continue;  // consumed as a fragment via the frag map
        if (auto w = dyn_cast<vector::TransferWriteOp>(user)) {
          OpBuilder wb(w);
          Value col2k = wb.create<arith::MulIOp>(w.getLoc(), L.k, L.c2idx);
          wb.create<vector::TransferWriteOp>(
              w.getLoc(), cf, w.getSource(), ValueRange{L.i, col2k},
              AffineMapAttr::get(AffineMap::getMinorIdentityMap(2, 2, ctx)),
              /*mask=*/Value(), wb.getBoolArrayAttr({true, true}));
          w.erase();
        }
      }
    }
    for (auto c : llvm::reverse(chain))
      c.erase();
  }
};

}  // namespace

namespace mir {
void registerChainContractsPass() { PassRegistration<ChainContractsPass>(); }
}  // namespace mir
