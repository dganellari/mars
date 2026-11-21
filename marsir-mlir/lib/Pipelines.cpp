// Registered pipelines: one mir-opt call instead of a driver script chaining several.
//
//   --mir-lower-to-nvvm="maxnreg=83 vector=true"
//       a scheduled module (func.func {mir.kernel}) -> gpu.module with the kernel as PTX
//       (gpu.binary, format=isa): --mir-gpu-wrap, the NVVM target (sm_90, O3), the
//       conversions to LLVM, and gpu-module-to-binary.
//   vector: the tensor-core / vector kernels (vector and nvgpu conversions);
//   math:   math.fma -> llvm.intr.fma (the cell-wise contractions; off for block DSS);
//   noalias: --mir-noalias-descriptors before serialization (the DSS kernels).

#include "mir/MirPasses.h"

#include "mlir/Pass/PassManager.h"
#include "mlir/Pass/PassOptions.h"
#include "mlir/Pass/PassRegistry.h"

using namespace mlir;

namespace {

struct NvvmLoweringOptions : public PassPipelineOptions<NvvmLoweringOptions> {
  Option<int> maxnreg{*this, "maxnreg", llvm::cl::desc("register cap on each kernel (0: none)"),
                      llvm::cl::init(0)};
  Option<std::string> chip{*this, "chip", llvm::cl::desc("NVIDIA architecture"), llvm::cl::init("sm_90")};
  Option<int> optLevel{*this, "O", llvm::cl::desc("NVPTX optimization level"), llvm::cl::init(3)};
  Option<bool> vector{*this, "vector", llvm::cl::desc("vector and nvgpu (tensor-core) conversions"),
                      llvm::cl::init(false)};
  Option<bool> math{*this, "math", llvm::cl::desc("math to LLVM (fma)"), llvm::cl::init(true)};
  Option<bool> noalias{*this, "noalias", llvm::cl::desc("--mir-noalias-descriptors before serialization"),
                       llvm::cl::init(false)};
};

void buildNvvmLowering(OpPassManager &pm, const NvvmLoweringOptions &o) {
  std::string p = "mir-gpu-wrap{maxnreg=" + std::to_string(o.maxnreg) + "},";
  p += "nvvm-attach-target{module=mir_kernels chip=" + o.chip + " O=" + std::to_string(o.optLevel) + "},";
  if (o.vector)
    p += "convert-vector-to-scf,canonicalize,";
  p += "convert-scf-to-cf,";
  if (o.vector)
    p += "convert-nvgpu-to-nvvm,";
  p += "expand-strided-metadata,lower-affine,";
  if (o.vector)
    p += "convert-vector-to-llvm,";
  if (o.math)
    p += "convert-math-to-llvm,";
  // convert-gpu-to-nvvm runs on gpu.module ops (the command line nests it implicitly).
  p += "gpu.module(convert-gpu-to-nvvm),convert-arith-to-llvm,convert-index-to-llvm,convert-cf-to-llvm,"
       "finalize-memref-to-llvm,reconcile-unrealized-casts,";
  if (o.noalias)
    p += "mir-noalias-descriptors,";
  p += "gpu-module-to-binary{format=isa}";
  if (failed(parsePassPipeline(p, pm, llvm::errs())))
    llvm::report_fatal_error("mir-lower-to-nvvm: bad pipeline " + llvm::Twine(p));
}

}  // namespace

namespace mir {
void registerPipelines() {
  PassPipelineRegistration<NvvmLoweringOptions>(
      "mir-lower-to-nvvm", "Scheduled MARSIR kernel -> PTX (gpu-wrap, NVVM target, LLVM, binary)",
      buildNvvmLowering);
}
}  // namespace mir
