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

mlir::LogicalResult mir::FluxOp::verify() {
  if (getResults().empty())
    return emitOpError("needs a result");
  auto grid = llvm::dyn_cast<mlir::RankedTensorType>(getResult(0).getType());
  if (!grid)
    return emitOpError("results must be ranked tensors");
  if (!getDims())
    return mlir::success();
  if (getDims()->size() != getInputs().size())
    return emitOpError("dims needs one entry per input");
  for (auto [input, attr] : llvm::zip(getInputs(), *getDims())) {
    auto dims = llvm::dyn_cast<mlir::DenseI64ArrayAttr>(attr);
    auto type = llvm::dyn_cast<mlir::RankedTensorType>(input.getType());
    if (!dims || !type || type.getRank() != (int64_t)dims.size())
      return emitOpError("each dims entry is an i64 array as long as its input's rank");
    int64_t last = -1;
    for (auto [k, d] : llvm::enumerate(dims.asArrayRef())) {
      if (d <= last || d >= grid.getRank())
        return emitOpError("dims must be increasing grid dimensions");
      if (type.getDimSize(k) != grid.getDimSize(d))
        return emitOpError("an input's extent must match the grid along its dimension");
      last = d;
    }
  }
  return mlir::success();
}

// Node numbers of a block stay well inside 32 bits up to these sizes.
static mlir::LogicalResult verifyBlock(mlir::Operation *op, int64_t degree, int64_t block) {
  if (degree < 1 || degree > 8 || block < 1 || block > 64)
    return op->emitOpError("degree must be in [1, 8] and block in [1, 64]");
  return mlir::success();
}

mlir::LogicalResult mir::BlockSumOp::verify() { return verifyBlock(*this, getDegree(), getBlock()); }

mlir::LogicalResult mir::BlockSkeletonOp::verify() { return verifyBlock(*this, getDegree(), getBlock()); }

mlir::LogicalResult mir::InterpolateOp::verify() {
  if (getCounting() && getAccumulate())
    return emitOpError("counting writes zeros on the other copies, so it cannot accumulate");
  return mlir::success();
}

mlir::LogicalResult mir::CellwiseOp::verify() {
  const int64_t r = getReductions();
  if (r < 0)
    return emitOpError("reductions must be >= 0");
  if ((r > 0) != (getPartials() != nullptr))
    return emitOpError("partials are needed exactly when reductions > 0");
  if (!getDirichlet().empty() && getDirichlet().size() != 5)
    return emitOpError("dirichlet takes (face_code, edge_of, vert_of, edge_dirichlet, vert_dirichlet)");
  const unsigned flags = getDirichlet().empty() ? 1 : 2;
  mlir::Block &body = getBody().front();
  if (body.getNumArguments() != flags + getInputs().size())
    return emitOpError("the body takes (counted[, dirichlet]) and one value per inputs array");
  auto yield = llvm::dyn_cast<mir::YieldOp>(body.getTerminator());
  if (!yield || yield.getNumOperands() != getOuts().size() + r)
    return emitOpError("the body yields one value per outs array, then the reduction terms");
  return mlir::success();
}
