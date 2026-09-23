// Static reasoning about the memory windows vector transfers access.
#ifndef MIR_WINDOWS_H
#define MIR_WINDOWS_H

#include "mlir/Dialect/MemRef/IR/MemRef.h"
#include "mlir/IR/BuiltinTypes.h"

#include <algorithm>

namespace mir {

// A memref value with any memref.cast peeled off (a cast changes only the type).
inline mlir::Value stripCasts(mlir::Value v) {
  while (auto c = v.getDefiningOp<mlir::memref::CastOp>())
    v = c.getSource();
  return v;
}

// Two subviews of the same memref whose static index ranges do not meet in
// some dimension. A transfer that is in bounds (or masked) never touches memory
// outside its window, so accesses through two such windows never overlap.
inline bool disjointWindows(mlir::Value a, mlir::Value b) {
  using namespace mlir;
  auto sa = stripCasts(a).getDefiningOp<memref::SubViewOp>();
  auto sb = stripCasts(b).getDefiningOp<memref::SubViewOp>();
  if (!sa || !sb || stripCasts(sa.getSource()) != stripCasts(sb.getSource()))
    return false;
  ArrayRef<int64_t> oa = sa.getStaticOffsets(), ob = sb.getStaticOffsets();
  ArrayRef<int64_t> za = sa.getStaticSizes(), zb = sb.getStaticSizes();
  ArrayRef<int64_t> ta = sa.getStaticStrides(), tb = sb.getStaticStrides();
  for (size_t d = 0; d < oa.size(); ++d) {
    if (ShapedType::isDynamic(oa[d]) || ShapedType::isDynamic(ob[d]) ||
        ShapedType::isDynamic(za[d]) || ShapedType::isDynamic(zb[d]) ||
        ShapedType::isDynamic(ta[d]) || ShapedType::isDynamic(tb[d]) ||
        za[d] < 1 || zb[d] < 1)
      continue;
    const int64_t endA = oa[d] + (za[d] - 1) * ta[d];
    const int64_t endB = ob[d] + (zb[d] - 1) * tb[d];
    const int64_t loA = std::min(oa[d], endA), hiA = std::max(oa[d], endA);
    const int64_t loB = std::min(ob[d], endB), hiB = std::max(ob[d], endB);
    if (hiA < loB || hiB < loA)
      return true;
  }
  return false;
}

}  // namespace mir

#endif  // MIR_WINDOWS_H
