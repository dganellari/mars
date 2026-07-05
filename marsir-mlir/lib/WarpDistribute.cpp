// --mir-warp-distribute: distribute vector code inside
// vector.warp_execute_on_lane_0 regions across the 32 lanes of a warp, using
// upstream's VectorDistribution pattern library (which ships in the MLIR
// install but is not exposed by any release pass -- same situation as the
// hoisting utilities, same remedy: link it from our own pass).
//
// What distributes today: transfer_write at the region boundary (each lane
// writes its slice), then propagation through elementwise ops, constants,
// broadcasts and transfer_read (each lane computes only its slice; cross-lane
// values become gpu.shuffle). vector.contract has NO upstream propagation
// pattern -- the WarpOpContractToMma pattern below bridges m8n8k4 f64
// contracts to nvgpu.mma.sync on per-lane DMMA fragments.
//
// Distribution map: greedily split the TRAILING dims of a vector so their
// product covers the 32 lanes (8x8 -> lanes (8,4), per-lane vector<1x2>;
// 8x8x8 -> (8,4), per-lane vector<1x2x8> ... trailing-major delinearization).

#include "mir/MirPasses.h"

#include "mlir/Dialect/Arith/IR/Arith.h"
#include "mlir/Dialect/GPU/IR/GPUDialect.h"
#include "mlir/Dialect/MemRef/IR/MemRef.h"
#include "mlir/Dialect/NVGPU/IR/NVGPUDialect.h"
#include "mlir/Dialect/SCF/IR/SCF.h"
#include "mlir/Dialect/Vector/IR/VectorOps.h"
#include "mlir/Dialect/Vector/Transforms/VectorDistribution.h"
#include "mlir/IR/AffineMap.h"
#include "mlir/IR/Builders.h"
#include "mlir/IR/PatternMatch.h"
#include "mlir/Pass/Pass.h"
#include "mlir/Transforms/GreedyPatternRewriteDriver.h"

#include <numeric>

using namespace mlir;

namespace {

// Split warpSize over the leading dims of `shape` (largest dims first is not
// needed: we walk from dim 0) such that each chosen dim divides evenly.
// Returns the affine map whose results are the distributed dims.
static AffineMap distributionMap(Value v) {
  auto ty = dyn_cast<VectorType>(v.getType());
  MLIRContext *ctx = v.getContext();
  if (!ty)
    return AffineMap::get(ctx);
  int64_t lanes = 32;
  SmallVector<AffineExpr> dims;
  for (int64_t d = 0; d < ty.getRank() && lanes > 1; ++d) {
    int64_t sz = ty.getDimSize(d);
    int64_t take = std::gcd(sz, lanes);
    if (take > 1) {
      dims.push_back(getAffineDimExpr(d, ctx));
      lanes /= take;
    }
  }
  return AffineMap::get(ty.getRank(), 0, dims, ctx);
}

static Value shuffleFromIdx(Location loc, OpBuilder &b, Value val, Value srcIdx,
                            int64_t warpSz) {
  Value idx32 = b.create<arith::IndexCastOp>(loc, b.getI32Type(), srcIdx);
  Value width = b.create<arith::ConstantOp>(
      loc, b.getI32IntegerAttr((int32_t)warpSz));
  auto shuf = b.create<gpu::ShuffleOp>(loc, val, idx32, width,
                                       gpu::ShuffleMode::IDX);
  return shuf.getShuffleResult();
}

// Distribute an m8n8k4-shaped f64 vector.contract that is yielded from a warp
// region. A and the C accumulator are re-yielded so upstream read propagation
// distributes them -- the default gcd map already produces the hardware
// fragment layouts A[L/4,L%4] (1x1) and C[L/4,2(L%4)] (1x2). B needs the
// TRANSPOSED layout B[L%4,L/4], which the map callback cannot express (it is
// invoked on the warp result, so per-op tags are invisible); since B is a
// transfer_read of a region-external memref, its per-lane fragment read is
// emitted directly outside the region with explicit lane indexing instead.
// nvgpu.mma.sync consumes the three fragments outside. Upstream has no
// contract propagation -- this is the mir bridge to tensor cores.
struct WarpOpContractToMma
    : public OpRewritePattern<vector::WarpExecuteOnLane0Op> {
  using OpRewritePattern<vector::WarpExecuteOnLane0Op>::OpRewritePattern;

  LogicalResult matchAndRewrite(vector::WarpExecuteOnLane0Op warp,
                                PatternRewriter &rewriter) const override {
    auto yield = cast<vector::YieldOp>(
        warp.getWarpRegion().front().getTerminator());
    MLIRContext *ctx = warp.getContext();

    // Canonical m8n8k4 iteration space: m=d0, n=d1, k=d2 (k reduced). The
    // sum-factorization emits contracts on different tensor axes, so the same
    // mma must be recovered from two operand layouts:
    //   standard : A=(m,k), B=(k,n)  -- the plain matmul (contract axis 0)
    //   transp-B : A=(m,k), B=(n,k)  -- B stored transposed (contract axis 1;
    //              e.g. in @ W^T, where vectorization gives rhs=(m,k) lhs=(n,k))
    // A and C are identical in both (default gcd distribution already yields
    // the hardware A[L/4,L%4] and C[L/4,2(L%4)] fragments); only B's storage
    // and hence its per-lane fragment index differ.
    auto d0 = getAffineDimExpr(0, ctx), d1 = getAffineDimExpr(1, ctx),
         d2 = getAffineDimExpr(2, ctx);
    auto mk = AffineMap::get(3, 0, {d0, d2}, ctx);  // A (m,k)
    auto kn = AffineMap::get(3, 0, {d2, d1}, ctx);  // B standard (k,n)
    auto nk = AffineMap::get(3, 0, {d1, d2}, ctx);  // B transposed (n,k)
    auto mn = AffineMap::get(3, 0, {d0, d1}, ctx);  // out (m,n)

    vector::ContractionOp con;
    unsigned resIdx = 0;
    Value aVal, bVal, accVal;
    bool bTransposed = false;
    for (auto [i, v] : llvm::enumerate(yield.getOperands())) {
      auto c = v.getDefiningOp<vector::ContractionOp>();
      if (!c)
        continue;
      auto maps = c.getIndexingMapsArray();
      if (maps.size() != 3 || maps[2] != mn)
        continue;
      auto iters = c.getIteratorTypesArray();
      if (iters.size() != 3 || iters[0] != vector::IteratorType::parallel ||
          iters[1] != vector::IteratorType::parallel ||
          iters[2] != vector::IteratorType::reduction)
        continue;
      // Assign A/B by matching maps; the two inputs may be in either order.
      Value lhs = c.getLhs(), rhs = c.getRhs();
      Value a, b;
      bool bt;
      if (maps[0] == mk && (maps[1] == kn || maps[1] == nk)) {
        a = lhs; b = rhs; bt = (maps[1] == nk);
      } else if (maps[1] == mk && (maps[0] == kn || maps[0] == nk)) {
        a = rhs; b = lhs; bt = (maps[0] == nk);
      } else {
        continue;
      }
      auto at = dyn_cast<VectorType>(a.getType());
      auto bTy = dyn_cast<VectorType>(b.getType());
      auto ct = dyn_cast<VectorType>(c.getResultType());
      // Transposed B is stored [n,k]=8x4; standard B is [k,n]=4x8. Compare
      // inline -- an ArrayRef bound to a ?:-temporary initializer_list would
      // dangle past the full expression.
      bool bShapeOk = bt ? (bTy && bTy.getShape() == ArrayRef<int64_t>({8, 4}))
                         : (bTy && bTy.getShape() == ArrayRef<int64_t>({4, 8}));
      if (!at || !bTy || !ct || !at.getElementType().isF64() ||
          at.getShape() != ArrayRef<int64_t>({8, 4}) || !bShapeOk ||
          ct.getShape() != ArrayRef<int64_t>({8, 8}))
        continue;
      con = c;
      resIdx = i;
      aVal = a; bVal = b; accVal = c.getAcc();
      bTransposed = bt;
      break;
    }
    if (!con)
      return failure();

    // B must be a plain read of data defined outside the region so the
    // fragment read can be re-created outside with lane indexing.
    auto bread = bVal.getDefiningOp<vector::TransferReadOp>();
    if (!bread || !bread.getPermutationMap().isMinorIdentity())
      return failure();
    Region &wregion = warp.getWarpRegion();
    auto definedOutside = [&](Value v) {
      return !wregion.isAncestor(v.getParentRegion());
    };
    if (!definedOutside(bread.getSource()) ||
        !llvm::all_of(bread.getIndices(), definedOutside))
      return failure();

    Location loc = warp.getLoc();
    auto f64 = Float64Type::get(ctx);
    auto fragA = VectorType::get({1, 1}, f64);
    auto fragB = VectorType::get({1, 1}, f64);
    // C flows through the existing result slot (already fragment-typed by the
    // boundary distribution when the map is the default gcd one).
    auto oldResTy = dyn_cast<VectorType>(warp.getResult(resIdx).getType());
    if (!oldResTy || oldResTy.getShape() != ArrayRef<int64_t>({1, 2}))
      return failure();

    // New warp op: same results, plus fragA appended (B is read outside).
    SmallVector<Type> newTypes(warp.getResultTypes());
    newTypes.push_back(fragA);
    auto newWarp = rewriter.create<vector::WarpExecuteOnLane0Op>(
        loc, newTypes, warp.getLaneid(), 32);
    newWarp.getWarpRegion().takeBody(warp.getWarpRegion());

    // Rewire the yield: slot resIdx now yields the ACC operand; append A.
    auto newYield = cast<vector::YieldOp>(
        newWarp.getWarpRegion().front().getTerminator());
    newYield->setOperand(resIdx, accVal);
    SmallVector<Value> yops(newYield.getOperands());
    yops.push_back(aVal);
    rewriter.modifyOpInPlace(newYield, [&] {
      newYield->setOperands(yops);
    });
    rewriter.eraseOp(con);

    // Outside: B fragment read. Hardware wants B_logical[L%4, L/4] (m8n8k4
    // row.col, B col-major). Standard B stored [k,n] -> index [L%4, L/4];
    // transposed B stored [n,k] -> B_logical[L%4,L/4] = Bstore[L/4, L%4].
    rewriter.setInsertionPointAfter(newWarp);
    Value lane = newWarp.getLaneid();
    Value c4v = rewriter.create<arith::ConstantIndexOp>(loc, 4);
    Value rem = rewriter.create<arith::RemUIOp>(loc, lane, c4v);
    Value div = rewriter.create<arith::DivUIOp>(loc, lane, c4v);
    Value brow = bTransposed ? div : rem;
    Value bcol = bTransposed ? rem : div;
    Value bi0 =
        rewriter.create<arith::AddIOp>(loc, bread.getIndices()[0], brow);
    Value bi1 =
        rewriter.create<arith::AddIOp>(loc, bread.getIndices()[1], bcol);
    Value pad = rewriter.create<arith::ConstantOp>(
        loc, f64, rewriter.getF64FloatAttr(0.0));
    auto fragBRead = rewriter.create<vector::TransferReadOp>(
        loc, fragB, bread.getSource(), ValueRange{bi0, bi1},
        AffineMap::getMinorIdentityMap(2, 2, ctx), pad, /*mask=*/Value(),
        rewriter.getBoolArrayAttr({true, true}));

    unsigned n = warp.getNumResults();
    auto mma = rewriter.create<nvgpu::MmaSyncOp>(
        loc, newWarp.getResult(n), fragBRead, newWarp.getResult(resIdx),
        rewriter.getDenseI64ArrayAttr({8, 8, 4}));

    // Old uses: the contract slot now carries the mma result; others map 1:1.
    for (auto [i, old] : llvm::enumerate(warp.getResults()))
      rewriter.replaceAllUsesWith(old, i == resIdx ? mma.getRes()
                                                   : newWarp.getResult(i));
    rewriter.eraseOp(warp);
    return success();
  }
};

struct WarpDistributePass
    : public PassWrapper<WarpDistributePass, OperationPass<>> {
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(WarpDistributePass)

  StringRef getArgument() const final { return "mir-warp-distribute"; }
  StringRef getDescription() const final {
    return "Distribute vector.warp_execute_on_lane_0 regions across lanes "
           "(upstream VectorDistribution patterns)";
  }
  void getDependentDialects(DialectRegistry &registry) const override {
    registry.insert<gpu::GPUDialect, vector::VectorDialect,
                    nvgpu::NVGPUDialect,
                    memref::MemRefDialect, arith::ArithDialect,
                    scf::SCFDialect>();
  }

  void runOnOperation() override {
    MLIRContext *ctx = &getContext();
    {
      RewritePatternSet patterns(ctx);
      patterns.add<WarpOpContractToMma>(ctx);
      vector::populateDistributeTransferWriteOpPatterns(
          patterns, distributionMap, /*maxNumElementsToExtract=*/1);
      vector::populatePropagateWarpVectorDistributionPatterns(
          patterns, distributionMap, shuffleFromIdx);
      if (failed(applyPatternsAndFoldGreedily(getOperation(),
                                              std::move(patterns))))
        signalPassFailure();
    }
  }
};

}  // namespace

namespace mir {
void registerWarpDistributePass() { PassRegistration<WarpDistributePass>(); }
}  // namespace mir
