// --mir-trace-shared (warp emulation only): before every load and store of a
// {mir.workgroup} buffer, a call __mir_shared_trace(id, offset, bytes, store), with
// offset the element offset into the buffer. The emulator's runtime
// (tools/emu_warp_runtime.h, MIR_EMU_TRACE) groups the k-th execution of an access
// by the lanes of a warp into one warp instruction and counts its shared-memory
// wavefronts as the hardware does (32 banks of 4 bytes, 64-bit accesses per half
// warp, one word read by several lanes is one broadcast). So the bank conflicts of
// a schedule are measured on the CPU, access by access.
//
// Ids number the accesses in walk order. With map-file=<path>, each id's access is
// written there with the kernel it is in.

#include "mir/MirPasses.h"

#include "mlir/Dialect/Arith/IR/Arith.h"
#include "mlir/Dialect/Func/IR/FuncOps.h"
#include "mlir/Dialect/MemRef/IR/MemRef.h"
#include "mlir/IR/BuiltinOps.h"
#include "mlir/Interfaces/ViewLikeInterface.h"
#include "mlir/Pass/Pass.h"

#include "llvm/Support/FileSystem.h"
#include "llvm/Support/raw_ostream.h"

using namespace mlir;

namespace {

constexpr StringLiteral kWorkgroupAttr = "mir.workgroup";
constexpr StringLiteral kHook = "__mir_shared_trace";

Value viewRoot(Value v) {
  while (auto view = v.getDefiningOp<ViewLikeOpInterface>()) v = view.getViewSource();
  return v;
}

struct TraceSharedPass : public PassWrapper<TraceSharedPass, OperationPass<ModuleOp>> {
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(TraceSharedPass)

  TraceSharedPass() = default;
  TraceSharedPass(const TraceSharedPass &other) : PassWrapper(other) {}

  Option<std::string> mapFile{*this, "map-file", llvm::cl::desc("write each trace id's access to this file")};

  StringRef getArgument() const final { return "mir-trace-shared"; }
  StringRef getDescription() const final {
    return "Warp emulation: report every shared-memory access (for bank-conflict counts)";
  }

  void runOnOperation() override {
    ModuleOp mod = getOperation();
    OpBuilder mb = OpBuilder::atBlockBegin(mod.getBody());
    Type i64 = mb.getI64Type();
    auto hook = mod.lookupSymbol<func::FuncOp>(kHook);
    if (!hook) {
      hook = mb.create<func::FuncOp>(mod.getLoc(), kHook, mb.getFunctionType({i64, i64, i64, i64}, {}));
      hook.setPrivate();
    }
    std::string text;
    llvm::raw_string_ostream map(text);
    int64_t id = 0;
    mod.walk([&](func::FuncOp f) {
      if (f.isExternal())
        return;
      DenseSet<Value> workgroup;
      for (unsigned i = 0; i < f.getNumArguments(); ++i)
        if (f.getArgAttr(i, kWorkgroupAttr))
          workgroup.insert(f.getArgument(i));
      if (workgroup.empty())
        return;
      SmallVector<Operation *> accesses;
      f.walk([&](Operation *op) {
        if (auto l = dyn_cast<memref::LoadOp>(op); l && workgroup.contains(viewRoot(l.getMemRef())))
          accesses.push_back(op);
        if (auto s = dyn_cast<memref::StoreOp>(op); s && workgroup.contains(viewRoot(s.getMemRef())))
          accesses.push_back(op);
      });
      for (Operation *op : accesses) {
        const bool store = isa<memref::StoreOp>(op);
        Value m = store ? cast<memref::StoreOp>(op).getMemRef() : cast<memref::LoadOp>(op).getMemRef();
        ValueRange idx = store ? cast<memref::StoreOp>(op).getIndices() : cast<memref::LoadOp>(op).getIndices();
        auto type = cast<MemRefType>(m.getType());
        OpBuilder b(op);
        Location loc = op->getLoc();
        auto meta = b.create<memref::ExtractStridedMetadataOp>(loc, m);
        Value offset = meta.getOffset();
        for (unsigned d = 0; d < idx.size(); ++d)
          offset = b.create<arith::AddIOp>(loc, offset, b.create<arith::MulIOp>(loc, idx[d], meta.getStrides()[d]));
        auto constant = [&](int64_t v) { return b.create<arith::ConstantIntOp>(loc, v, 64).getResult(); };
        b.create<func::CallOp>(loc, hook,
                               ValueRange{constant(id), b.create<arith::IndexCastOp>(loc, i64, offset),
                                          constant(type.getElementTypeBitWidth() / 8), constant(store)});
        map << id << (store ? " st " : " ld ") << f.getName() << ": ";
        op->print(map, OpPrintingFlags().useLocalScope());
        map << "\n";
        ++id;
      }
    });
    if (mapFile.empty())
      return;
    std::error_code ec;
    llvm::raw_fd_ostream out(mapFile, ec, llvm::sys::fs::OF_Text);
    if (ec) {
      mod.emitError("mir-trace-shared: cannot write ") << mapFile << ": " << ec.message();
      return signalPassFailure();
    }
    out << map.str();
  }
};

}  // namespace

namespace mir {
void registerTraceSharedPass() { PassRegistration<TraceSharedPass>(); }
}  // namespace mir
