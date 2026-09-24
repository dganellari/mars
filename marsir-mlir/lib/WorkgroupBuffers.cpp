// --mir-workgroup-buffers: give a warp kernel's scratch buffers to the WARP.
//
// A memref.alloc in a kernel body is per-thread memory (an alloca after
// --promote-buffers-to-stack): 32 private copies of a value the warp computes
// once. A C-fragment can go into such a buffer only after a full gather (the
// coherence rule of --mir-chain-contracts), and every copy is local-memory
// traffic. Workgroup (shared) memory holds ONE copy the whole warp sees: each
// lane stores its own fragment entries, and --mir-warp-barriers orders the
// stages that read them.
//
// Every static, identity-layout alloc at the top level of a {mir.kernel}
// func.func becomes a reinterpret_cast view of a new trailing argument marked
// {mir.workgroup}; --mir-gpu-wrap turns those into gpu.func workgroup
// attributions. Buffers whose live ranges (in the body's top-level order) do not
// overlap share one slot: the scratch of one sweep direction is reused by the
// next instead of adding up. Shared memory per block bounds how many warps an SM
// holds, so the reuse matters for occupancy, not only for footprint.
//
// Runs AFTER --mir-forward-transfers, which only forwards through allocs.

#include "mir/MirPasses.h"

#include "mlir/Dialect/Func/IR/FuncOps.h"
#include "mlir/Dialect/MemRef/IR/MemRef.h"
#include "mlir/IR/Builders.h"
#include "mlir/IR/BuiltinOps.h"
#include "mlir/Interfaces/CallInterfaces.h"
#include "mlir/Pass/Pass.h"

using namespace mlir;

namespace {

constexpr StringLiteral kKernelAttr = "mir.kernel";
constexpr StringLiteral kWorkgroupAttr = "mir.workgroup";

struct WorkgroupBuffersPass
    : public PassWrapper<WorkgroupBuffersPass, OperationPass<ModuleOp>> {
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(WorkgroupBuffersPass)

  StringRef getArgument() const final { return "mir-workgroup-buffers"; }
  StringRef getDescription() const final {
    return "Move a {mir.kernel} function's static scratch allocs into shared "
           "{mir.workgroup} slots, reusing a slot once its buffer is dead";
  }
  void getDependentDialects(DialectRegistry &registry) const override {
    registry.insert<memref::MemRefDialect>();
  }

  void runOnOperation() override {
    for (auto fn : getOperation().getOps<func::FuncOp>())
      if (fn->hasAttr(kKernelAttr) && !fn.isDeclaration())
        runOn(fn);
  }

  void runOn(func::FuncOp fn) {
    Block *entry = &fn.getBody().front();
    DenseMap<Operation *, unsigned> pos;
    unsigned k = 0;
    for (Operation &op : *entry)
      pos[&op] = k++;
    auto topLevel = [&](Operation *op) {
      while (op->getBlock() != entry)
        op = op->getParentOp();
      return op;
    };

    struct Buf {
      memref::AllocOp op;
      unsigned start, end;
      int slot;
    };
    SmallVector<Buf> bufs;
    for (auto a : entry->getOps<memref::AllocOp>()) {
      MemRefType t = a.getType();
      if (!t.hasStaticShape() || !t.getLayout().isIdentity() ||
          t.getMemorySpace() || t.getRank() == 0)
        continue;
      // Live range: to the last top-level op that touches the buffer or any
      // view of it. A buffer that leaves through a call or a terminator has a
      // lifetime this cannot see -- keep it private.
      unsigned end = pos[a];
      bool escapes = false;
      SmallVector<Value> work{a.getResult()};
      while (!work.empty()) {
        Value v = work.pop_back_val();
        for (Operation *u : v.getUsers()) {
          if (isa<CallOpInterface>(u) || u->hasTrait<OpTrait::IsTerminator>())
            escapes = true;
          end = std::max(end, pos[topLevel(u)]);
          for (Value r : u->getResults())
            if (isa<MemRefType>(r.getType()))
              work.push_back(r);
        }
      }
      if (!escapes)
        bufs.push_back({a, pos[a], end, -1});
    }
    if (bufs.empty())
      return;

    // First fit in program order. A slot is sized by the largest buffer it holds.
    struct Slot {
      Type elem;
      int64_t size;
      unsigned freeAfter;
    };
    SmallVector<Slot> slots;
    for (Buf &b : bufs) {
      MemRefType t = b.op.getType();
      for (int s = 0; s < (int)slots.size(); ++s)
        if (slots[s].elem == t.getElementType() && slots[s].freeAfter < b.start) {
          b.slot = s;
          break;
        }
      if (b.slot < 0) {
        b.slot = (int)slots.size();
        slots.push_back({t.getElementType(), 0, 0});
      }
      slots[b.slot].size = std::max(slots[b.slot].size, t.getNumElements());
      slots[b.slot].freeAfter = b.end;
    }

    OpBuilder ab(fn.getContext());
    SmallVector<Value> slotArgs;
    for (Slot &s : slots) {
      const unsigned i = fn.getNumArguments();
      fn.insertArgument(i, MemRefType::get({s.size}, s.elem),
                        ab.getDictionaryAttr(ab.getNamedAttr(kWorkgroupAttr,
                                                             ab.getUnitAttr())),
                        fn.getLoc());
      slotArgs.push_back(fn.getArgument(i));
    }

    for (Buf &b : bufs) {
      MemRefType t = b.op.getType();
      SmallVector<int64_t> strides(t.getRank(), 1);
      for (int64_t d = t.getRank() - 2; d >= 0; --d)
        strides[d] = strides[d + 1] * t.getDimSize(d + 1);
      // A slot is never freed; a dealloc of the view would free the argument.
      for (Operation *u : llvm::make_early_inc_range(b.op->getUsers()))
        if (isa<memref::DeallocOp>(u))
          u->erase();
      OpBuilder rb(b.op);
      Value view = rb.create<memref::ReinterpretCastOp>(
          b.op.getLoc(), t, slotArgs[b.slot], /*offset=*/0, t.getShape(), strides);
      b.op.getResult().replaceAllUsesWith(view);
      b.op.erase();
    }
  }
};

}  // namespace

namespace mir {
void registerWorkgroupBuffersPass() { PassRegistration<WorkgroupBuffersPass>(); }
}  // namespace mir
