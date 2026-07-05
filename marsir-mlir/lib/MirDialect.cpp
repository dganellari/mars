#include "mir/MirDialect.h"
#include "mir/MirOps.h"

#include "llvm/ADT/TypeSwitch.h"
#include "mlir/IR/Builders.h"
#include "mlir/IR/DialectImplementation.h"

using namespace mir;

// Generated dialect definitions.
#include "mir/MirOpsDialect.cpp.inc"

// Generated type definitions (the !mir.block type).
#define GET_TYPEDEF_CLASSES
#include "mir/MirOpsTypes.cpp.inc"

void MirDialect::initialize() {
  addOperations<
#define GET_OP_LIST
#include "mir/MirOps.cpp.inc"
      >();
  addTypes<
#define GET_TYPEDEF_LIST
#include "mir/MirOpsTypes.cpp.inc"
      >();
}
