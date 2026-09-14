// --mir-batch-elements: turn a single-element kernel into the warp-per-element
// batched kernel. This is the pass form of what mlir_warp.py hand-wrote as
// build_full_batched_kernel.
//
// Every gpu.func argument marked {mir.element} gains a leading dynamic batch
// dimension, and the body is rebased onto a per-element subview indexed by
// gpu.block_id x:
//
//   gpu.func @k(%u: memref<8x8x8xf64> {mir.element}, %B: memref<8x8xf64>)
//     ->
//   gpu.func @k(%u: memref<?x8x8x8xf64>, %B: memref<8x8xf64>) {
//     %e  = gpu.block_id x
//     %el = memref.subview %u[%e, 0, 0, 0] [1, 8, 8, 8] [1, 1, 1, 1]
//     ... body, with %u replaced by %el ...
//
// The marker is explicit rather than inferred: the operator matrices (Btil, W)
// are memrefs too, and a shape heuristic would silently batch them.
//
// Launch contract afterwards: grid.x = element count, blockDim.x = 32 (the
// warp), which is what --mir-warp-wrap / --mir-warp-distribute assume.

#include "mir/MirPasses.h"

#include "mlir/Dialect/Func/IR/FuncOps.h"
#include "mlir/Dialect/GPU/IR/GPUDialect.h"
#include "mlir/Dialect/MemRef/IR/MemRef.h"
#include "mlir/IR/Builders.h"
#include "mlir/IR/BuiltinTypes.h"
#include "mlir/Pass/Pass.h"

using namespace mlir;

namespace {

constexpr StringLiteral kElementAttr = "mir.element";

// A subview's result layout is derived from its source, so rebasing the body
// onto a dynamically-offset element view invalidates every subview below it
// (offset: 8 becomes offset: ?). Re-infer until nothing moves.
static void refreshSubViewTypes(Operation *func) {
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
  }
}

// Templated over the function op: the hex path arrives as a gpu.func from the
// emitter, the tet path as a bufferized func.func. Both expose the same
// function_type / arg-attribute / body API.
template <typename FuncT>
static void batchOne(FuncT func) {
  FunctionType fnTy = func.getFunctionType();
  SmallVector<unsigned> marked;
  for (unsigned i = 0, n = fnTy.getNumInputs(); i < n; ++i)
    if (func.getArgAttr(i, kElementAttr))
      marked.push_back(i);
  if (marked.empty())
    return;

  MLIRContext *ctx = func.getContext();
  Block &body = func.getBody().front();
  OpBuilder b = OpBuilder::atBlockBegin(&body);
  Location loc = func.getLoc();

  // Prepend the dynamic batch dimension to every marked argument.
  SmallVector<Type> inTypes(fnTy.getInputs());
  for (unsigned i : marked) {
    auto mt = cast<MemRefType>(inTypes[i]);
    SmallVector<int64_t> shape{ShapedType::kDynamic};
    llvm::append_range(shape, mt.getShape());
    inTypes[i] = MemRefType::get(shape, mt.getElementType(),
                                 MemRefLayoutAttrInterface{}, mt.getMemorySpace());
    body.getArgument(i).setType(inTypes[i]);
  }
  func.setFunctionType(FunctionType::get(ctx, inTypes, fnTy.getResults()));

  Value e = b.create<gpu::BlockIdOp>(loc, gpu::Dimension::x);

  for (unsigned i : marked) {
    Value arg = body.getArgument(i);
    auto batched = cast<MemRefType>(arg.getType());
    ArrayRef<int64_t> elemShape = batched.getShape().drop_front();

    // [%e, 0..] [1, shape..] [1, 1..] -- rank-reduced back to the element shape,
    // so the body sees exactly the rank it had before.
    SmallVector<OpFoldResult> offs{OpFoldResult(e)};
    SmallVector<OpFoldResult> sizes{b.getIndexAttr(1)};
    SmallVector<OpFoldResult> strides{b.getIndexAttr(1)};
    for (int64_t d : elemShape) {
      offs.push_back(b.getIndexAttr(0));
      sizes.push_back(b.getIndexAttr(d));
      strides.push_back(b.getIndexAttr(1));
    }
    auto elemTy = cast<MemRefType>(memref::SubViewOp::inferRankReducedResultType(
        elemShape, batched, offs, sizes, strides));
    Value el = b.create<memref::SubViewOp>(loc, elemTy, arg, offs, sizes, strides);
    arg.replaceAllUsesExcept(el, el.getDefiningOp());

    // Drop the marker so a second run is a no-op rather than a double batch.
    func.removeArgAttr(i, kElementAttr);
  }

  refreshSubViewTypes(func.getOperation());
}

struct BatchElementsPass
    : public PassWrapper<BatchElementsPass, OperationPass<>> {
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(BatchElementsPass)

  StringRef getArgument() const final { return "mir-batch-elements"; }
  StringRef getDescription() const final {
    return "Batch a single-element gpu.func over elements: {mir.element} args "
           "gain a leading dynamic dim and the body is rebased onto a "
           "gpu.block_id x subview";
  }
  void getDependentDialects(DialectRegistry &registry) const override {
    registry.insert<gpu::GPUDialect, memref::MemRefDialect>();
  }

  void runOnOperation() override {
    getOperation()->walk([](gpu::GPUFuncOp f) { batchOne(f); });
    getOperation()->walk([](func::FuncOp f) { batchOne(f); });
  }
};

}  // namespace

namespace mir {
void registerBatchElementsPass() { PassRegistration<BatchElementsPass>(); }
}  // namespace mir
