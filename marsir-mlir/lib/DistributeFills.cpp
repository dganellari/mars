// --mir-distribute-fills: a splat written into memory the warp SHARES is written
// once per element, by one lane, instead of by all 32.
//
// Outside the fragment chain a kernel body runs redundantly in every lane, so a
// zero fill of a 7x64 scratch buffer is 448 stores per lane: the same 448 values
// 32 times. Memory the warp shares (a kernel argument, or a {mir.workgroup}
// slot) needs one copy, so lane L stores elements L, L+32, L+64, ... of the
// flattened window. Consecutive lanes store consecutive elements, so the stores
// coalesce. Other lanes' elements arrive from other lanes, and
// --mir-warp-barriers orders that like any cross-lane store.
//
// Only in-bounds, unmasked, minor-identity writes of a splat constant qualify;
// per-thread memory keeps the redundant fill (each lane needs its own copy).

#include "mir/MirPasses.h"

#include "mlir/Dialect/Arith/IR/Arith.h"
#include "mlir/Dialect/GPU/IR/GPUDialect.h"
#include "mlir/Dialect/MemRef/IR/MemRef.h"
#include "mlir/Dialect/SCF/IR/SCF.h"
#include "mlir/Dialect/Vector/IR/VectorOps.h"
#include "mlir/IR/Builders.h"
#include "mlir/Interfaces/FunctionInterfaces.h"
#include "mlir/Interfaces/ViewLikeInterface.h"
#include "mlir/Pass/Pass.h"

using namespace mlir;

namespace {

// Shared = rooted, through views, at an argument of the enclosing function.
static bool isShared(Value mem) {
  while (Operation *d = mem.getDefiningOp()) {
    auto v = dyn_cast<ViewLikeOpInterface>(d);
    if (!v)
      return false;
    mem = v.getViewSource();
  }
  auto ba = cast<BlockArgument>(mem);
  Operation *owner = ba.getOwner()->getParentOp();
  return isa<FunctionOpInterface>(owner) && ba.getOwner()->isEntryBlock();
}

struct DistributeFillsPass
    : public PassWrapper<DistributeFillsPass, OperationPass<>> {
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(DistributeFillsPass)

  StringRef getArgument() const final { return "mir-distribute-fills"; }
  StringRef getDescription() const final {
    return "Spread a splat write into warp-shared memory across the 32 lanes";
  }
  void getDependentDialects(DialectRegistry &registry) const override {
    registry.insert<arith::ArithDialect, gpu::GPUDialect, memref::MemRefDialect,
                    scf::SCFDialect>();
  }

  void runOnOperation() override {
    SmallVector<vector::TransferWriteOp> fills;
    getOperation()->walk([&](vector::TransferWriteOp w) {
      auto cst = w.getVector().getDefiningOp<arith::ConstantOp>();
      auto dv = cst ? dyn_cast<DenseElementsAttr>(cst.getValue()) : nullptr;
      if (!dv || !dv.isSplat() || w.getMask() ||
          !w.getPermutationMap().isMinorIdentity() ||
          !isa<MemRefType>(w.getSource().getType()) ||
          llvm::any_of(w.getInBoundsValues(), [](bool b) { return !b; }) ||
          !isShared(w.getSource()))
        return;
      fills.push_back(w);
    });

    for (vector::TransferWriteOp w : fills) {
      OpBuilder b(w);
      Location loc = w.getLoc();
      auto vt = w.getVectorType();
      auto dv = cast<DenseElementsAttr>(
          w.getVector().getDefiningOp<arith::ConstantOp>().getValue());
      Value scalar = b.create<arith::ConstantOp>(
          loc, cast<TypedAttr>(dv.getSplatValue<Attribute>()));
      Value lane = b.create<gpu::ThreadIdOp>(loc, gpu::Dimension::x);
      Value total = b.create<arith::ConstantIndexOp>(loc, vt.getNumElements());
      Value step = b.create<arith::ConstantIndexOp>(loc, 32);
      auto loop = b.create<scf::ForOp>(loc, lane, total, step);
      OpBuilder lb = OpBuilder::atBlockBegin(loop.getBody());

      // Element j of the flattened window -> its position in the window.
      const int64_t r = vt.getRank();
      SmallVector<Value> pos(r);
      Value rest = loop.getInductionVar();
      for (int64_t d = r - 1; d >= 0; --d) {
        if (d == 0) {
          pos[d] = rest;
          break;
        }
        Value dim = lb.create<arith::ConstantIndexOp>(loc, vt.getDimSize(d));
        pos[d] = lb.create<arith::RemUIOp>(loc, rest, dim);
        rest = lb.create<arith::DivUIOp>(loc, rest, dim);
      }
      // Minor identity: the vector covers the LAST r dimensions of the buffer.
      SmallVector<Value> idx(w.getIndices().begin(), w.getIndices().end());
      const int64_t lead = (int64_t)idx.size() - r;
      for (int64_t d = 0; d < r; ++d)
        idx[lead + d] = lb.create<arith::AddIOp>(loc, idx[lead + d], pos[d]);
      lb.create<memref::StoreOp>(loc, scalar, w.getSource(), idx);
      w.erase();
    }
  }
};

}  // namespace

namespace mir {
void registerDistributeFillsPass() { PassRegistration<DistributeFillsPass>(); }
}  // namespace mir
