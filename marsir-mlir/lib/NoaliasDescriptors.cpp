// --mir-noalias-descriptors: after the NVVM lowering, mark the data pointer of
// every memref argument of a {mir.noalias} kernel `llvm.noalias`.
//
// A kernel marked {mir.noalias} promises that no two of its memref arguments
// overlap (like __restrict__ on every pointer of a CUDA kernel). With that, LLVM
// may keep loads of buffers the kernel never writes on the read-only data path
// (ld.global.nc), as nvcc does for `const T* __restrict__`.
//
// The lowering passes drop argument and function attributes of a gpu.func, so
// --mir-gpu-wrap records the marked kernels on the gpu.module, which survives
// (mir.noalias_kernels). In the lowered function each memref argument arrives as
// separate scalars and is rebuilt into a descriptor struct; the data (aligned)
// pointer is the value inserted at position 1. Only that pointer gets the
// attribute: the allocated pointer at position 0 is the same address and is never
// dereferenced, so marking it too would make two "noalias" names for one buffer.
//
// The backend uses ld.global.nc only for arguments that are noalias AND readonly,
// and it does not infer readonly here, so a pointer whose derived values are only
// ever loaded from (no store through them, no store of them, no other use that
// could write) also gets `llvm.readonly`.

#include "mir/MirPasses.h"

#include "mlir/Dialect/GPU/IR/GPUDialect.h"
#include "mlir/Dialect/LLVMIR/LLVMDialect.h"
#include "mlir/IR/BuiltinOps.h"
#include "mlir/IR/TypeUtilities.h"
#include "mlir/Pass/Pass.h"

using namespace mlir;
using mir::kMirNoaliasKernelsAttr;

namespace {

// True when nothing derived from `ptr` (through descriptor insert/extract, address
// arithmetic and casts) is stored to or escapes. Loads are fine, and so is the
// alignment assumption (memref.assume_alignment lowers to ptrtoint, and, icmp,
// llvm.assume): its integer only feeds the test. Integer fields extracted from a
// descriptor (sizes, strides) are not the pointer.
bool onlyRead(Value ptr) {
  SmallVector<Value> work{ptr};
  DenseSet<Value> seen{ptr};
  auto follow = [&](Value r) {
    if (seen.insert(r).second)
      work.push_back(r);
  };
  while (!work.empty()) {
    Value v = work.pop_back_val();
    const bool integer = isa<IntegerType>(v.getType()) || isa<IntegerType>(getElementTypeOrSelf(v.getType()));
    for (OpOperand &use : v.getUses()) {
      Operation *user = use.getOwner();
      if (integer) {   // the address as an integer: only the alignment test may use it
        if (isa<LLVM::AndOp>(user) || isa<LLVM::ICmpOp>(user))
          follow(user->getResult(0));
        else if (!isa<LLVM::AssumeOp>(user))
          return false;
        continue;
      }
      if (isa<LLVM::LoadOp, LLVM::AssumeOp>(user))
        continue;
      if (auto ext = dyn_cast<LLVM::ExtractValueOp>(user)) {
        if (isa<LLVM::LLVMPointerType, LLVM::LLVMStructType>(ext.getType()))
          follow(ext.getResult());
        continue;
      }
      if (isa<LLVM::InsertValueOp, LLVM::GEPOp, LLVM::AddrSpaceCastOp, LLVM::PtrToIntOp>(user)) {
        follow(user->getResult(0));
        continue;
      }
      return false;   // a store, a call, an escape: not provably read-only
    }
  }
  return true;
}

struct NoaliasDescriptorsPass : public PassWrapper<NoaliasDescriptorsPass, OperationPass<ModuleOp>> {
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(NoaliasDescriptorsPass)

  StringRef getArgument() const final { return "mir-noalias-descriptors"; }
  StringRef getDescription() const final {
    return "Mark the data pointers of the memref arguments of {mir.noalias} kernels "
           "llvm.noalias (and llvm.readonly when never written), after the NVVM lowering";
  }

  void runOnOperation() override {
    getOperation()->walk([&](gpu::GPUModuleOp gmod) {
      auto kernels = gmod->getAttrOfType<ArrayAttr>(kMirNoaliasKernelsAttr);
      if (!kernels)
        return;
      for (Attribute name : kernels) {
        auto fn = gmod.lookupSymbol<LLVM::LLVMFuncOp>(cast<StringAttr>(name));
        if (!fn) {
          gmod.emitError("{mir.noalias} kernel ") << name << " not found after lowering";
          return signalPassFailure();
        }
        fn.walk([&](LLVM::InsertValueOp ins) {
          auto desc = dyn_cast<LLVM::LLVMStructType>(ins.getContainer().getType());
          auto arg = dyn_cast<BlockArgument>(ins.getValue());
          // A memref descriptor: (allocated ptr, aligned ptr, offset, ...).
          if (!desc || !arg || arg.getOwner() != &fn.getBody().front() || ins.getPosition() != ArrayRef<int64_t>{1} ||
              desc.getBody().size() < 3 || !isa<LLVM::LLVMPointerType>(desc.getBody()[0]) ||
              !isa<LLVM::LLVMPointerType>(desc.getBody()[1]))
            return;
          fn.setArgAttr(arg.getArgNumber(), LLVM::LLVMDialect::getNoAliasAttrName(), UnitAttr::get(fn.getContext()));
          if (onlyRead(arg))
            fn.setArgAttr(arg.getArgNumber(), LLVM::LLVMDialect::getReadonlyAttrName(), UnitAttr::get(fn.getContext()));
        });
      }
      gmod->removeAttr(kMirNoaliasKernelsAttr);
    });
  }
};

}  // namespace

namespace mir {
void registerNoaliasDescriptorsPass() { PassRegistration<NoaliasDescriptorsPass>(); }
}  // namespace mir
