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

mlir::LogicalResult mir::DssOp::verify() {
  const bool epilogue = !getEpilogue().empty();
  const int64_t r = getReductions();
  if (r < 0)
    return emitOpError("reductions must be >= 0");
  if (!getFlags().empty() && getFlags().size() != 3)
    return emitOpError("flags are (counting, edge_dirichlet, vert_dirichlet)");
  if (!epilogue) {
    if (getOuts().size() != 1 || r != 0 || !getCopies().empty() || getPartials())
      return emitOpError("without an epilogue: one outs array, no copies, reductions or partials");
    return mlir::success();
  }
  if (getFlags().size() != 3)
    return emitOpError("an epilogue needs flags(counting, edge_dirichlet, vert_dirichlet)");
  if ((r > 0) != (getPartials() != nullptr))
    return emitOpError("partials are needed exactly when reductions > 0");
  mlir::Block &body = getEpilogue().front();
  if (body.getNumArguments() != 3 + getCopies().size())
    return emitOpError("the epilogue takes (sum, dirichlet, counted) and one value per copies array");
  auto yield = llvm::dyn_cast<mir::YieldOp>(body.getTerminator());
  if (!yield || yield.getNumOperands() != getOuts().size() + r)
    return emitOpError("the epilogue yields one value per outs array, then the reduction terms");
  return mlir::success();
}

mlir::LogicalResult mir::CellwiseOp::verify() {
  const int64_t r = getReductions();
  if (r < 0)
    return emitOpError("reductions must be >= 0");
  if ((r > 0) != (getPartials() != nullptr))
    return emitOpError("partials are needed exactly when reductions > 0");
  mlir::Block &body = getBody().front();
  if (body.getNumArguments() != 1 + getInputs().size())
    return emitOpError("the body takes (counted) and one value per inputs array");
  auto yield = llvm::dyn_cast<mir::YieldOp>(body.getTerminator());
  if (!yield || yield.getNumOperands() != getOuts().size() + r)
    return emitOpError("the body yields one value per outs array, then the reduction terms");
  return mlir::success();
}
