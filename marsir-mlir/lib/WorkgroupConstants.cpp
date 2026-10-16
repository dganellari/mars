// --mir-workgroup-constants: copy a kernel's small constant tables into shared
// memory once per block.
//
// Lookup tables (memref.get_global of a constant memref.global, e.g. the
// reference-hex tables of --mir-lower-dss) are read with lane-dependent indices,
// so constant memory would serialise them, and from global memory they compete
// in L1 with the streaming data and often come from L2: in the DSS each copy's
// address is a chain table -> table -> value, and every miss in that chain is
// latency. In shared memory a lookup is a few cycles. This is what the
// hand-written MARS kernel does (build_hex_tables, once per thread block).
//
// At kernel entry every thread copies a strided share of each table into a new
// {mir.workgroup} argument (shared memory after --mir-gpu-wrap), then one
// gpu.barrier, then every use of the global reads the copy. Worth it only when a
// block handles many elements (grid-stride, see --mir-warp-per-element), so the
// copy is paid once per block.
//
// Runs on {mir.kernel} func.funcs after --mir-warp-per-element.

#include "mir/MirPasses.h"

#include "mlir/Dialect/Arith/IR/Arith.h"
#include "mlir/Dialect/Func/IR/FuncOps.h"
#include "mlir/Dialect/GPU/IR/GPUDialect.h"
#include "mlir/Dialect/MemRef/IR/MemRef.h"
#include "mlir/Dialect/SCF/IR/SCF.h"
#include "mlir/IR/Builders.h"
#include "mlir/IR/BuiltinOps.h"
#include "mlir/Pass/Pass.h"

using namespace mlir;

namespace {

constexpr int64_t kMaxBytes = 16 * 1024;   // a table, not a data array

struct WorkgroupConstantsPass : public PassWrapper<WorkgroupConstantsPass, OperationPass<ModuleOp>> {
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(WorkgroupConstantsPass)

  StringRef getArgument() const final { return "mir-workgroup-constants"; }
  StringRef getDescription() const final {
    return "Copy the small constant tables a kernel reads into shared memory once per block";
  }
  void getDependentDialects(DialectRegistry &registry) const override {
    registry.insert<arith::ArithDialect, gpu::GPUDialect, memref::MemRefDialect, scf::SCFDialect>();
  }

  void runOnOperation() override {
    ModuleOp mod = getOperation();
    for (auto func : mod.getOps<func::FuncOp>()) {
      if (!func->hasAttr("mir.kernel") || func.isDeclaration())
        continue;
      SmallVector<memref::GetGlobalOp> gets;
      func.walk([&](memref::GetGlobalOp get) {
        auto global = mod.lookupSymbol<memref::GlobalOp>(get.getName());
        auto type = get.getType();
        if (global && global.getConstant() && type.hasStaticShape() && type.getRank() == 1 &&
            type.getNumElements() * ((type.getElementTypeBitWidth() + 7) / 8) <= kMaxBytes)
          gets.push_back(get);
      });
      if (gets.empty())
        continue;

      Block &entry = func.getBody().front();
      OpBuilder b = OpBuilder::atBlockBegin(&entry);
      Location loc = func.getLoc();
      Value tid = b.create<gpu::ThreadIdOp>(loc, gpu::Dimension::x);
      Value threads = b.create<gpu::BlockDimOp>(loc, gpu::Dimension::x);
      DenseMap<StringAttr, Value> copies;   // one copy per table
      for (memref::GetGlobalOp get : gets) {
        StringAttr name = get.getNameAttr().getAttr();
        auto it = copies.find(name);
        if (it == copies.end()) {
          const unsigned at = func.getNumArguments();
          func.insertArgument(at, get.getType(),
                              b.getDictionaryAttr(b.getNamedAttr("mir.workgroup", b.getUnitAttr())), loc);
          Value copy = func.getArgument(at);
          Value table = b.create<memref::GetGlobalOp>(loc, get.getType(), name);
          Value size = b.create<arith::ConstantIndexOp>(loc, get.getType().getNumElements());
          b.create<scf::ForOp>(loc, tid, size, threads, ValueRange{},
                               [&](OpBuilder &lb, Location l, Value i, ValueRange) {
                                 Value v = lb.create<memref::LoadOp>(l, table, ValueRange{i});
                                 lb.create<memref::StoreOp>(l, v, copy, ValueRange{i});
                                 lb.create<scf::YieldOp>(l);
                               });
          it = copies.insert({name, copy}).first;
        }
        get.replaceAllUsesWith(it->second);
        get.erase();
      }
      b.create<gpu::BarrierOp>(loc);
    }
  }
};

}  // namespace

namespace mir {
void registerWorkgroupConstantsPass() { PassRegistration<WorkgroupConstantsPass>(); }
}  // namespace mir
