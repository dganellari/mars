// Registered pipelines: one mir-opt call instead of a driver script chaining several.
//
//   --mir-lower-to-nvvm="maxnreg=83 vector=true"
//       a scheduled module (func.func {mir.kernel}) -> gpu.module with the kernel as PTX
//       (gpu.binary, format=isa): --mir-gpu-wrap, the NVVM target (sm_90, O3), the
//       conversions to LLVM, and gpu-module-to-binary.
//   vector: the tensor-core / vector kernels (vector and nvgpu conversions);
//   math:   math.fma -> llvm.intr.fma (the cell-wise contractions; off for block DSS);
//   noalias: --mir-noalias-descriptors before serialization (the DSS kernels).
//   --mir-lower-to-rocdl="chip=gfx942": the same for AMD, up to the LLVM dialect.

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

// --mir-bp-schedule="spec=bp3 p=3 target=amd": the cell-wise schedule of a CEED BP operator
// (bp1, bp3, bp1_tri, bp3_tri) with the launch shapes and rewrites measured fastest per target.
// lanes / warps / even-odd / pad / lane-owned override the tables.
struct BpScheduleOptions : public PassPipelineOptions<BpScheduleOptions> {
  Option<std::string> spec{*this, "spec", llvm::cl::desc("bp1, bp3, bp1_tri or bp3_tri"), llvm::cl::init("bp3")};
  Option<std::string> target{*this, "target", llvm::cl::desc("nvidia or amd"), llvm::cl::init("nvidia")};
  Option<int> p{*this, "p", llvm::cl::desc("polynomial degree"), llvm::cl::init(3)};
  Option<int> q{*this, "q", llvm::cl::desc("quadrature points per direction (0: p + 2)"), llvm::cl::init(0)};
  Option<int> lanes{*this, "lanes", llvm::cl::desc("lanes per element (-1: the table)"), llvm::cl::init(-1)};
  Option<int> warps{*this, "warps", llvm::cl::desc("warps per block, used as given (0: the table)"),
                    llvm::cl::init(0)};
  Option<bool> evenOdd{*this, "even-odd", llvm::cl::desc("the even-odd split where the table has it"),
                       llvm::cl::init(true)};
  Option<bool> pad{*this, "pad", llvm::cl::desc("scratch padding against bank conflicts"), llvm::cl::init(true)};
  Option<bool> laneOwned{*this, "lane-owned", llvm::cl::desc("lane-owned scratch in registers"),
                         llvm::cl::init(false)};
  Option<bool> columns{*this, "columns", llvm::cl::desc("3-D phases keep z in the thread (column schedule)"),
                       llvm::cl::init(false)};
};

void buildBpSchedule(OpPassManager &pm, const BpScheduleOptions &o) {
  const std::string s = o.spec;
  const int p = o.p, q = o.q ? (int)o.q : p + 2;
  const bool amd = o.target == "amd";
  const int wave = amd ? 64 : 32;
  auto in = [&](std::initializer_list<std::pair<const char *, int>> table) {
    for (auto [sp, deg] : table)
      if (s == sp && p == deg) return true;
    return false;
  };
  // At p = 7 (q = 9) scratch per element limits the resident warps (BP3: four q^3 buffers,
  // 8 warps per SM), so two warps share an element: GH200 operator -26% BP3, -16% BP1. Three
  // or four warps lose the register blocking and are 2.4x slower. On MI300A 64 lanes are one
  // wavefront (757112: 2 elements per block fastest, which fit-shared gives).
  int lanes = o.lanes >= 0 ? (int)o.lanes : (in({{"bp1", 7}, {"bp3", 7}}) ? 64 : 0);
  // Warps per block: the CVFEM operator's measured default, which fit-shared halves (fewer
  // elements per block) until the static shared memory fits the target's limit; an explicit
  // count is used as given. MI300A, BP3 p = 3: one wave per block, 1.01 ms against 1.18 ms
  // with four (757112: shared memory per element limits the waves per CU either way).
  int warps = o.warps;
  if (!warps && lanes > wave) warps = lanes / wave;
  if (!warps && amd && in({{"bp3", 3}})) warps = 1;
  const bool fit = !warps;
  if (!warps) warps = p >= 3 ? 4 : 8;
  // The even-odd split (--mir-even-odd) is faster at p = 7 (BP3 10%, BP1 32%, BP3 trilinear
  // 12%) and for bp3_tri p = 3 (3%), but bp3 p = 3 is 6% slower with it and p = 1 gains
  // nothing. bp1_tri is not timed; it follows bp1.
  const bool evenOdd = o.evenOdd && (p == 7 || in({{"bp3_tri", 3}}));
  // Less scratch is not always faster (docs/design/cellwise_bp_fcg.md, "Scratch reuse"): bp1
  // p = 7 is 5% slower with largest-first packing (4 instead of 2 elements per block), and
  // bp3 p = 1 keeps the form measured fastest.
  const bool firstFit = in({{"bp3", 1}, {"bp1", 7}});
  // One thread per 1D line (max-block = q), unrolled; the matrices' constant address space
  // steers which dimension a thread keeps.
  std::string pl = "mir-cellwise-schedule{warps=" + std::to_string(warps) + " max-block=" + std::to_string(q) +
                   " unroll-reductions=true const-space=4" + (lanes ? " lanes=" + std::to_string(lanes) : "") +
                   (o.columns ? " columns=true" : "") +
                   "},canonicalize,cse,";
  if (evenOdd)
    pl += "mir-even-odd,canonicalize,cse,";
  pl += "mir-fuse-phases,canonicalize,cse,";
  pl += std::string("mir-warp-per-element{largest-first=") + (firstFit ? "false" : "true") +
        " fit-shared=" + (fit ? "true" : "false") + " pad-scratch=" + (o.pad ? "true" : "false") +
        " lane-owned=" + (o.laneOwned ? "true" : "false") + "},canonicalize,cse,loop-invariant-code-motion";
  if (failed(parsePassPipeline(pl, pm, llvm::errs())))
    llvm::report_fatal_error("mir-bp-schedule: bad pipeline " + llvm::Twine(pl));
}

struct RocdlLoweringOptions : public PassPipelineOptions<RocdlLoweringOptions> {
  Option<std::string> chip{*this, "chip", llvm::cl::desc("AMD architecture"), llvm::cl::init("gfx942")};
  Option<int> optLevel{*this, "O", llvm::cl::desc("AMDGPU optimization level"), llvm::cl::init(3)};
  Option<bool> noalias{*this, "noalias", llvm::cl::desc("--mir-noalias-descriptors at the end (the DSS kernels)"),
                       llvm::cl::init(false)};
};

// The AMD counterpart up to the LLVM dialect; translating and compiling the kernels (llvm-
// translate, opt, llc) is outside MLIR (test/amd_target.py): gpu-module-to-binary would need
// the ROCm device libraries even for kernels that call none.
void buildRocdlLowering(OpPassManager &pm, const RocdlLoweringOptions &o) {
  std::string p = "mir-gpu-wrap,mir-amdgpu-grid-dim,";
  p += "rocdl-attach-target{module=mir_kernels chip=" + o.chip + " O=" + std::to_string(o.optLevel) + "},";
  p += "convert-vector-to-scf,canonicalize,convert-scf-to-cf,expand-strided-metadata,lower-affine,";
  p += "convert-amdgpu-to-rocdl{chipset=" + o.chip + "},convert-vector-to-llvm,convert-math-to-llvm,";
  // convert-gpu-to-rocdl runs on gpu.module ops (the command line nests it implicitly).
  p += "gpu.module(convert-gpu-to-rocdl{chipset=" + o.chip + "}),";
  p += "convert-arith-to-llvm,convert-index-to-llvm,convert-cf-to-llvm,finalize-memref-to-llvm,"
       "reconcile-unrealized-casts";
  if (o.noalias)
    p += ",mir-noalias-descriptors";
  if (failed(parsePassPipeline(p, pm, llvm::errs())))
    llvm::report_fatal_error("mir-lower-to-rocdl: bad pipeline " + llvm::Twine(p));
}

}  // namespace

namespace mir {
void registerPipelines() {
  PassPipelineRegistration<BpScheduleOptions>(
      "mir-bp-schedule", "Cell-wise schedule of a CEED BP operator, tuned per target", buildBpSchedule);
  PassPipelineRegistration<RocdlLoweringOptions>(
      "mir-lower-to-rocdl", "Scheduled MARSIR kernel -> AMDGPU LLVM dialect (gpu-wrap, grid size, ROCDL target)",
      buildRocdlLowering);
  PassPipelineRegistration<NvvmLoweringOptions>(
      "mir-lower-to-nvvm", "Scheduled MARSIR kernel -> PTX (gpu-wrap, NVVM target, LLVM, binary)",
      buildNvvmLowering);
}
}  // namespace mir
