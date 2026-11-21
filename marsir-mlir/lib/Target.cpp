// --mir-set-target: record the GPU's properties on the module, as an upstream DLTI target
// system spec (device "GPU"), so passes read them instead of assuming NVIDIA.
//
//   module attributes {dlti.target_system_spec = #dlti.target_system_spec<"GPU":
//       #dlti.target_device_spec<#dlti.dl_entry<"mir.vendor", "amd">,
//                                #dlti.dl_entry<"mir.wave_size", 64 : i64>,
//                                #dlti.dl_entry<"mir.max_shared_bytes", 65536 : i64>>>}

#include "mir/Target.h"

#include "mlir/Dialect/DLTI/DLTI.h"
#include "mlir/IR/Builders.h"
#include "mlir/IR/BuiltinOps.h"
#include "mlir/Interfaces/DataLayoutInterfaces.h"
#include "mlir/Pass/Pass.h"
#include "mlir/Pass/PassRegistry.h"

using namespace mlir;
using namespace mir;

mir::GpuTarget mir::getGpuTarget(Operation *op) {
  GpuTarget t;
  auto module = isa<ModuleOp>(op) ? cast<ModuleOp>(op) : op->getParentOfType<ModuleOp>();
  // Kernels sit in a gpu.module inside the host module: take the nearest module with a spec.
  while (module && !module.getTargetSystemSpec())
    module = module->getParentOfType<ModuleOp>();
  if (!module)
    return t;
  DataLayout layout(module);
  auto gpu = StringAttr::get(op->getContext(), "GPU");
  auto get = [&](const char *key) {
    return layout.getDevicePropertyValue(gpu, StringAttr::get(op->getContext(), key));
  };
  // An entry of the wrong type keeps the default rather than asserting.
  if (auto v = get(kVendorKey); v && isa<StringAttr>(*v))
    t.vendor = cast<StringAttr>(*v).getValue();
  if (auto v = get(kWaveSizeKey); v && isa<IntegerAttr>(*v))
    t.waveSize = cast<IntegerAttr>(*v).getInt();
  if (auto v = get(kMaxSharedBytesKey); v && isa<IntegerAttr>(*v))
    t.maxSharedBytes = cast<IntegerAttr>(*v).getInt();
  return t;
}

namespace {

struct SetTargetPass : public PassWrapper<SetTargetPass, OperationPass<ModuleOp>> {
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(SetTargetPass)

  SetTargetPass() = default;
  SetTargetPass(const SetTargetPass &other) : PassWrapper(other) {}

  StringRef getArgument() const final { return "mir-set-target"; }
  StringRef getDescription() const final {
    return "Record the GPU's wave size and shared-memory limit on the module (DLTI)";
  }
  void getDependentDialects(DialectRegistry &registry) const override {
    registry.insert<DLTIDialect>();
  }

  Option<std::string> target{*this, "target", llvm::cl::desc("nvidia (sm_90) or amd (gfx942)"),
                             llvm::cl::init("nvidia")};

  void runOnOperation() override {
    // sm_90: 48 KB of static shared memory per block without the opt-in; gfx942: 64 KB LDS.
    const bool amd = target == "amd";
    if (!amd && target != "nvidia") {
      getOperation().emitError("--mir-set-target: unknown target '") << target << "'";
      return signalPassFailure();
    }
    MLIRContext *ctx = &getContext();
    Builder b(ctx);
    auto entry = [&](const char *key, Attribute value) -> DataLayoutEntryInterface {
      return DataLayoutEntryAttr::get(StringAttr::get(ctx, key), value);
    };
    SmallVector<DataLayoutEntryInterface> entries = {
        entry(kVendorKey, b.getStringAttr(amd ? "amd" : "nvidia")),
        entry(kWaveSizeKey, b.getI64IntegerAttr(amd ? 64 : 32)),
        entry(kMaxSharedBytesKey, b.getI64IntegerAttr(amd ? 65536 : 49152))};
    auto device = TargetDeviceSpecAttr::get(ctx, entries);
    auto system = TargetSystemSpecAttr::get(
        ctx, {DeviceIDTargetDeviceSpecPair(StringAttr::get(ctx, "GPU"), device)});
    getOperation()->setAttr(DLTIDialect::kTargetSystemDescAttrName, system);
  }
};

}  // namespace

void mir::registerSetTargetPass() { PassRegistration<SetTargetPass>(); }
