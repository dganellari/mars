// --mir-forward-transfers: store-to-load forwarding for vector transfers.
//
// Tiling a contraction stages its accumulator in a scratch buffer, so the
// vectorized IR reads back a value it just wrote:
//
//   %a = memref.alloc()
//   vector.transfer_write %zero, %a[0, 0]
//   %acc = vector.transfer_read  %a[0, 0]      <- this IS %zero
//   %r   = vector.contract %x, %y, %acc
//   vector.transfer_write %r, %a[0, 0]         <- and a later read IS %r
//
// Upstream --canonicalize does not forward these. Left alone they cost a real
// round-trip, and worse, they hide the dataflow from --mir-chain-contracts:
// the accumulator looks like an opaque memory read instead of a zero splat, and
// one contraction's result looks like a leaf read instead of the previous
// contraction's C-fragment, so nothing chains in registers.
//
// SAFETY: only memref.alloc buffers whose every user is a vector transfer are
// considered. That rules out aliasing through subview/reshape entirely -- no
// view of the buffer can exist -- and distinct allocations cannot alias each
// other. A write invalidates every recorded entry for its own buffer.

#include "mir/MirPasses.h"

#include "mlir/Dialect/MemRef/IR/MemRef.h"
#include "mlir/Dialect/Vector/IR/VectorOps.h"
#include "mlir/IR/Builders.h"
#include "mlir/Pass/Pass.h"

using namespace mlir;

namespace {

// A buffer whose entire dataflow we can see.
static bool isPrivateAlloc(Value m) {
  if (!m.getDefiningOp<memref::AllocOp>())
    return false;
  for (Operation *u : m.getUsers())
    if (!isa<vector::TransferReadOp, vector::TransferWriteOp>(u))
      return false;
  return true;
}

static bool sameIndices(ValueRange a, ValueRange b) {
  if (a.size() != b.size())
    return false;
  for (auto [x, y] : llvm::zip(a, b))
    if (x != y)
      return false;
  return true;
}

struct ForwardTransfersPass
    : public PassWrapper<ForwardTransfersPass, OperationPass<>> {
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(ForwardTransfersPass)

  StringRef getArgument() const final { return "mir-forward-transfers"; }
  StringRef getDescription() const final {
    return "Forward vector.transfer_write values to matching transfer_reads on "
           "private alloc buffers (store-to-load)";
  }
  void getDependentDialects(DialectRegistry &registry) const override {
    registry.insert<memref::MemRefDialect, vector::VectorDialect>();
  }

  void runOnOperation() override {
    struct Entry { Value mem; SmallVector<Value> idx; Value val; };

    SmallVector<Block *> blocks;
    getOperation()->walk([&](Operation *op) {
      for (Region &r : op->getRegions())
        for (Block &b : r)
          blocks.push_back(&b);
    });

    for (Block *blk : blocks) {
      SmallVector<Entry> live;
      for (Operation &opRef : llvm::make_early_inc_range(*blk)) {
        Operation *op = &opRef;

        if (auto w = dyn_cast<vector::TransferWriteOp>(op)) {
          if (!isPrivateAlloc(w.getSource()))
            continue;
          // A write to this buffer supersedes everything known about it.
          llvm::erase_if(live, [&](const Entry &e) { return e.mem == w.getSource(); });
          live.push_back({w.getSource(),
                          SmallVector<Value>(w.getIndices().begin(),
                                             w.getIndices().end()),
                          w.getVector()});
          continue;
        }

        if (auto r = dyn_cast<vector::TransferReadOp>(op)) {
          if (!isPrivateAlloc(r.getSource()))
            continue;
          for (const Entry &e : live) {
            if (e.mem != r.getSource() || !sameIndices(e.idx, r.getIndices()))
              continue;
            if (e.val.getType() != r.getType())
              break;   // a differently-shaped view of the same window
            r.getResult().replaceAllUsesWith(e.val);
            r.erase();
            break;
          }
        }
      }
    }
  }
};

}  // namespace

namespace mir {
void registerForwardTransfersPass() { PassRegistration<ForwardTransfersPass>(); }
}  // namespace mir
