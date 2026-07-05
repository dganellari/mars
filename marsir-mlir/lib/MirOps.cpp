#include "mir/MirOps.h"

#include "mir/MirDialect.h"
#include "mlir/IR/Builders.h"
#include "mlir/IR/OpImplementation.h"

using namespace mir;

// Generated op definitions.
#define GET_OP_CLASSES
#include "mir/MirOps.cpp.inc"

// Zero-operand terminator builder for mir.flux's implicit terminator.
void mir::YieldOp::build(mlir::OpBuilder &, mlir::OperationState &) {}

// unbind(bind(x)) -> x when the round-trip returns the same tensor type.
mlir::OpFoldResult mir::UnbindOp::fold(FoldAdaptor) {
  if (auto bind = getBlock().getDefiningOp<mir::BindOp>())
    if (bind.getData().getType() == getType())
      return bind.getData();
  return {};
}
