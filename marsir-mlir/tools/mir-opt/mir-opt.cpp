// mir-opt: the dialect's optimizer/driver, like mlir-opt but with the mir
// dialect (and, later, the mir->linalg lowering passes) registered.
#include "mir/MirDialect.h"
#include "mir/MirPasses.h"

#include "mlir/IR/DialectRegistry.h"
#include "mlir/InitAllDialects.h"
#include "mlir/InitAllPasses.h"
#include "mlir/Tools/mlir-opt/MlirOptMain.h"

int main(int argc, char **argv) {
  mlir::DialectRegistry registry;
  mlir::registerAllDialects(registry);
  mlir::registerAllPasses();

  registry.insert<mir::MirDialect>();
  mir::registerMirPasses();
  mir::registerHoistTransferPairsPass();
  mir::registerLowerCopiesPass();
  mir::registerWarpDistributePass();
  mir::registerWarpWrapPass();

  return mlir::asMainReturnCode(mlir::MlirOptMain(
      argc, argv, "mir optimizer driver (Stage 2 MLIR dialect)\n", registry));
}
