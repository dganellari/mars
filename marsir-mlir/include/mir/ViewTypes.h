// Re-infer the result types of memref view ops after their source changed type
// in place -- a batched argument gaining a leading dimension, or a buffer moving
// into workgroup memory. The shape, layout and memory space of a view all follow
// its source, and a stale result type would fail verification.
#ifndef MIR_VIEWTYPES_H
#define MIR_VIEWTYPES_H

#include "mlir/Dialect/MemRef/IR/MemRef.h"

namespace mir {

inline void refreshViewTypes(mlir::Operation *func) {
  using namespace mlir;
  bool changed = true;
  while (changed) {
    changed = false;
    func->walk([&](memref::SubViewOp sv) {
      auto src = cast<MemRefType>(sv.getSource().getType());
      auto want = cast<MemRefType>(memref::SubViewOp::inferRankReducedResultType(
          sv.getType().getShape(), src, sv.getMixedOffsets(),
          sv.getMixedSizes(), sv.getMixedStrides()));
      if (want != sv.getType()) {
        sv.getResult().setType(want);
        changed = true;
      }
    });
    func->walk([&](memref::CollapseShapeOp cs) {
      auto src = cast<MemRefType>(cs.getSrc().getType());
      MemRefType want = memref::CollapseShapeOp::computeCollapsedType(
          src, cs.getReassociationIndices());
      if (want != cs.getType()) {
        cs.getResult().setType(want);
        changed = true;
      }
    });
    func->walk([&](memref::ExpandShapeOp es) {
      auto src = cast<MemRefType>(es.getSrc().getType());
      FailureOr<MemRefType> want = memref::ExpandShapeOp::computeExpandedType(
          src, es.getType().getShape(), es.getReassociationIndices());
      if (succeeded(want) && *want != es.getType()) {
        es.getResult().setType(*want);
        changed = true;
      }
    });
    // A reinterpret_cast states its own shape and layout; only the memory space
    // follows the source.
    func->walk([&](memref::ReinterpretCastOp rc) {
      auto src = cast<MemRefType>(rc.getSource().getType());
      MemRefType cur = rc.getType();
      if (cur.getMemorySpace() == src.getMemorySpace())
        return;
      rc.getResult().setType(MemRefType::get(cur.getShape(), cur.getElementType(),
                                             cur.getLayout(),
                                             src.getMemorySpace()));
      changed = true;
    });
  }
}

}  // namespace mir

#endif  // MIR_VIEWTYPES_H
