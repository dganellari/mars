// --mir-gpu-wrap: move a func.func marked {mir.kernel} into a gpu.module as a
// gpu.func kernel. This is the last step of the tet pipeline, where the operator
// arrives as a bufferized, batched func.func rather than as a gpu.func from the
// hex emitter:
//
//   mir  -> --convert-mir-to-linalg
//        -> one-shot-bufferize + buffer-results-to-out-params   (memrefs, DPS)
//        -> --mir-batch-elements                                (warp per element)
//        -> --mir-gpu-wrap                                      (a real kernel)
//
// The marker is an attribute rather than "any function that uses gpu.block_id",
// so that wrapping stays something the caller asks for.

#include "mir/MirPasses.h"
#include "mir/ViewTypes.h"

#include "mlir/Dialect/Func/IR/FuncOps.h"
#include "mlir/Dialect/GPU/IR/GPUDialect.h"
#include "mlir/IR/Builders.h"
#include "mlir/IR/BuiltinOps.h"
#include "mlir/Pass/Pass.h"

using namespace mlir;

namespace {

constexpr StringLiteral kKernelAttr = "mir.kernel";
constexpr StringLiteral kWorkgroupAttr = "mir.workgroup";
constexpr StringLiteral kModuleName = "mir_kernels";

struct GpuWrapPass : public PassWrapper<GpuWrapPass, OperationPass<ModuleOp>> {
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(GpuWrapPass)

  StringRef getArgument() const final { return "mir-gpu-wrap"; }
  StringRef getDescription() const final {
    return "Move each func.func marked {mir.kernel} into a gpu.module as a "
           "gpu.func kernel";
  }
  void getDependentDialects(DialectRegistry &registry) const override {
    registry.insert<gpu::GPUDialect, func::FuncDialect>();
  }

  void runOnOperation() override {
    ModuleOp mod = getOperation();
    SmallVector<func::FuncOp> marked;
    for (auto f : mod.getOps<func::FuncOp>())
      if (f->hasAttr(kKernelAttr))
        marked.push_back(f);
    if (marked.empty())
      return;

    OpBuilder b(mod.getBodyRegion());
    b.setInsertionPointToEnd(mod.getBody());

    // Reuse the gpu.module if a previous run (or the hex path) already made one.
    gpu::GPUModuleOp gmod;
    for (auto m : mod.getOps<gpu::GPUModuleOp>())
      if (m.getName() == kModuleName) {
        gmod = m;
        break;
      }
    if (!gmod)
      gmod = b.create<gpu::GPUModuleOp>(mod.getLoc(), kModuleName);

    for (func::FuncOp f : marked) {
      // gpu.module carries a gpu.module_end terminator; insert before it.
      Block *gbody = gmod.getBody();
      OpBuilder gb = gbody->mightHaveTerminator()
                         ? OpBuilder::atBlockTerminator(gbody)
                         : OpBuilder::atBlockEnd(gbody);
      auto gfunc = gb.create<gpu::GPUFuncOp>(f.getLoc(), f.getName(),
                                             f.getFunctionType());
      gfunc->setAttr(gpu::GPUDialect::getKernelFuncAttrName(), gb.getUnitAttr());
      // takeBody replaces the entry block the builder just made, so the argument
      // values the body already refers to come across unchanged.
      gfunc.getBody().takeBody(f.getBody());
      // A kernel ends in gpu.return, and only a void one can: a result would
      // have to come back through an out-parameter.
      gfunc.walk([&](func::ReturnOp r) {
        OpBuilder rb(r);
        rb.create<gpu::ReturnOp>(r.getLoc(), r.getOperands());
        r.erase();
      });
      // {mir.workgroup} arguments (--mir-workgroup-buffers) are the block's
      // shared memory, not something the host passes: each becomes a workgroup
      // attribution, and every view of it moves to the workgroup address space.
      // Argument attributes live on the func.func, not on its body.
      SmallVector<unsigned> wg;
      for (unsigned i = 0; i < f.getNumArguments(); ++i)
        if (f.getArgAttr(i, kWorkgroupAttr))
          wg.push_back(i);
      // The numeric shared-memory space (3 on NVVM, and LDS on AMDGPU), not
      // #gpu.address_space<workgroup>: the separate vector/memref-to-LLVM passes
      // of the lowering have no mapping for the attribute (only
      // --convert-gpu-to-nvvm does), and would leave every transfer on shared
      // memory unconverted.
      auto space = b.getI64IntegerAttr(3);
      for (unsigned i : wg) {
        auto t = cast<MemRefType>(gfunc.getArgument(i).getType());
        BlockArgument attr = gfunc.addWorkgroupAttribution(
            MemRefType::get(t.getShape(), t.getElementType(), t.getLayout(), space),
            gfunc.getLoc());
        gfunc.getArgument(i).replaceAllUsesWith(attr);
      }
      for (unsigned i : llvm::reverse(wg))
        gfunc.eraseArgument(i);
      mir::refreshViewTypes(gfunc);
      f.erase();
    }
  }
};

}  // namespace

namespace mir {
void registerGpuWrapPass() { PassRegistration<GpuWrapPass>(); }
}  // namespace mir
