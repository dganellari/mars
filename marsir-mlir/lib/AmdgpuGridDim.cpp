// --mir-amdgpu-grid-dim: on AMD targets, gpu.grid_dim (blocks in the grid) read from the HSA
// dispatch packet the kernel receives with every launch. --convert-gpu-to-rocdl would turn it
// into a call to __ockl_get_num_groups, a ROCm device-library function, which then has to be
// linked; the packet holds the same numbers (workgroup_size_x/y/z: u16 at bytes 4/6/8,
// grid_size_x/y/z in work-items: u32 at 12/16/20; code object v4 and v5), and the block count
// is ceil(grid size / workgroup size). Runs before --convert-gpu-to-rocdl.

#include "mir/MirPasses.h"
#include "mir/Target.h"

#include "mlir/Dialect/Arith/IR/Arith.h"
#include "mlir/Dialect/GPU/IR/GPUDialect.h"
#include "mlir/Dialect/LLVMIR/LLVMDialect.h"
#include "mlir/IR/Builders.h"
#include "mlir/IR/BuiltinOps.h"
#include "mlir/Pass/Pass.h"

using namespace mlir;

namespace {

struct AmdgpuGridDimPass : public PassWrapper<AmdgpuGridDimPass, OperationPass<ModuleOp>> {
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(AmdgpuGridDimPass)

  StringRef getArgument() const final { return "mir-amdgpu-grid-dim"; }
  StringRef getDescription() const final {
    return "On AMD targets, read gpu.grid_dim from the dispatch packet (no device library)";
  }
  void getDependentDialects(DialectRegistry &registry) const override {
    registry.insert<arith::ArithDialect, LLVM::LLVMDialect>();
  }

  void runOnOperation() override {
    if (!mir::getGpuTarget(getOperation()).isAmd())
      return;
    SmallVector<gpu::GridDimOp> uses;
    getOperation().walk([&](gpu::GridDimOp op) { uses.push_back(op); });
    for (gpu::GridDimOp op : uses) {
      OpBuilder b(op);
      Location loc = op.getLoc();
      const int64_t d = static_cast<int64_t>(op.getDimension());
      auto i8 = b.getI8Type(), i16 = b.getI16Type(), i32 = b.getI32Type();
      auto constant = LLVM::LLVMPointerType::get(b.getContext(), 4);
      Value packet = b.create<LLVM::CallIntrinsicOp>(loc, TypeRange{constant},
                                                     b.getStringAttr("llvm.amdgcn.dispatch.ptr"), ValueRange{})
                         .getResult(0);
      auto field = [&](int64_t offset, Type type, unsigned align) -> Value {
        Value at = b.create<LLVM::GEPOp>(loc, constant, i8, packet, ArrayRef<LLVM::GEPArg>{offset});
        auto load = b.create<LLVM::LoadOp>(loc, type, at, align);
        load.setInvariant(true);
        return load;
      };
      Value grid = field(12 + 4 * d, i32, 4);
      Value width = b.create<arith::ExtUIOp>(loc, i32, field(4 + 2 * d, i16, 2));
      Value one = b.create<arith::ConstantIntOp>(loc, 1, 32);
      Value blocks = b.create<arith::DivUIOp>(
          loc, b.create<arith::SubIOp>(loc, b.create<arith::AddIOp>(loc, grid, width), one), width);
      op.replaceAllUsesWith(b.create<arith::IndexCastUIOp>(loc, b.getIndexType(), blocks).getResult());
      op.erase();
    }
  }
};

}  // namespace

namespace mir {
void registerAmdgpuGridDimPass() { PassRegistration<AmdgpuGridDimPass>(); }
}  // namespace mir
