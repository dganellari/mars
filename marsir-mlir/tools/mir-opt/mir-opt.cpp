// mir-opt: the dialect's optimizer/driver, like mlir-opt but with the mir
// dialect (and, later, the mir->linalg lowering passes) registered.
#include "mir/MirDialect.h"
#include "mir/MirPasses.h"
#include "mir/Target.h"

#include "mlir/IR/DialectRegistry.h"
#include "mlir/InitAllDialects.h"
#include "mlir/InitAllExtensions.h"
#include "mlir/Target/LLVMIR/Dialect/All.h"
#include "mlir/InitAllPasses.h"
#include "mlir/Tools/mlir-opt/MlirOptMain.h"

int main(int argc, char **argv) {
  mlir::DialectRegistry registry;
  mlir::registerAllDialects(registry);
  mlir::registerAllExtensions(registry);
  // gpu-module-to-binary translates the kernels to LLVM IR (--mir-lower-to-nvvm).
  mlir::registerAllToLLVMIRTranslations(registry);
  mlir::registerAllPasses();

  registry.insert<mir::MirDialect>();
  mir::registerMirPasses();
  mir::registerHoistTransferPairsPass();
  mir::registerLowerCopiesPass();
  mir::registerWarpDistributePass();
  mir::registerWarpWrapPass();
  mir::registerChainContractsPass();
  mir::registerBatchElementsPass();
  mir::registerGpuWrapPass();
  mir::registerForwardTransfersPass();
  mir::registerEmulateWarpPass();
  mir::registerWorkgroupBuffersPass();
  mir::registerDistributeFillsPass();
  mir::registerWarpBarriersPass();
  mir::registerHoistInvariantReadsPass();
  mir::registerUnrollLoopsPass();
  mir::registerForwardOwnedPass();
  mir::registerLowerCellwisePass();
  mir::registerWarpPerElementPass();
  mir::registerNoaliasDescriptorsPass();
  mir::registerCellwiseSchedulePass();
  mir::registerFusePhasesPass();
  mir::registerTraceSharedPass();
  mir::registerEvenOddPass();
  mir::registerSetTargetPass();
  mir::registerMarkKernelPass();
  mir::registerPipelines();
  mir::registerAmdgpuGridDimPass();

  return mlir::asMainReturnCode(mlir::MlirOptMain(
      argc, argv, "mir optimizer driver (Stage 2 MLIR dialect)\n", registry));
}
