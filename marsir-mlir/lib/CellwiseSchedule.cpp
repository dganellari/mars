// --mir-cellwise-schedule: the schedule for operators whose elements are small (low
// polynomial degree): many elements per block, one group of lanes per element,
// thread per output point (the shape of MARS's hand-written cvfem_apply_kernel),
// instead of the tensor-core schedule of a p = 7 element. It takes the operator as
// the operator pipeline has it after bufferization (one element; linalg on memrefs)
// and produces the per-element phase form that --mir-warp-per-element maps:
//
//   1. {mir.element} arguments gain a leading batch dimension and the body moves
//      into an element loop scf.parallel {mir.elements}, on per-element views.
//   2. {mir.stage = "in"} arguments are copied into shared scratch first, {mir.stage
//      = "out"} ones are built in scratch and copied out last: the element's field
//      and result are touched many times, the hand kernel stages them too.
//   3. Short loops with constant bounds (the loop over the faces) are unrolled.
//   4. Every linalg op and memref.copy becomes one phase, an scf.parallel {mir.lanes}
//      over its parallel loops; reduction loops run inside, accumulating in a
//      register. A linalg.fill whose buffer the next writer overwrites entirely as
//      a reduction output is not run at all: the accumulator starts at its value.
//   5. Temporaries (memref.alloc) become shared scratch {mir.scratch}.
//
// Options: lanes per element (default: n^2 of the element's n^3 nodes, at most 32)
// and warps per block (default 8: 256 threads), set on the element loop for
// --mir-warp-per-element.

#include "mir/MirPasses.h"
#include "mir/ViewTypes.h"

#include "mlir/Dialect/Arith/IR/Arith.h"
#include "mlir/Dialect/Func/IR/FuncOps.h"
#include "mlir/Dialect/Linalg/IR/Linalg.h"
#include "mlir/Dialect/Linalg/Transforms/Transforms.h"
#include "mlir/Dialect/MemRef/IR/MemRef.h"
#include "mlir/Dialect/SCF/IR/SCF.h"
#include "mlir/Dialect/SCF/Utils/Utils.h"
#include "mlir/Dialect/Utils/StaticValueUtils.h"
#include "mlir/IR/Builders.h"
#include "mlir/IR/BuiltinOps.h"
#include "mlir/IR/IRMapping.h"
#include "mlir/IR/PatternMatch.h"
#include "mlir/Interfaces/ViewLikeInterface.h"
#include "mlir/Pass/Pass.h"

using namespace mlir;

namespace {

constexpr StringLiteral kElementAttr = "mir.element", kStageAttr = "mir.stage";

Value rootOf(Value m) {
  while (Operation *def = m.getDefiningOp()) {
    auto view = dyn_cast<ViewLikeOpInterface>(def);
    if (!view)
      break;
    m = view.getViewSource();
  }
  return m;
}

Value cidx(OpBuilder &b, Location loc, int64_t v) { return b.create<arith::ConstantIndexOp>(loc, v); }

// One phase over `shape`: scf.parallel {mir.lanes}, body(ivs).
void phase(OpBuilder &b, Location loc, ArrayRef<int64_t> shape,
           function_ref<void(OpBuilder &, Location, ValueRange)> body) {
  SmallVector<Value> lbs, ubs, steps;
  for (int64_t s : shape) {
    lbs.push_back(cidx(b, loc, 0));
    ubs.push_back(cidx(b, loc, s));
    steps.push_back(cidx(b, loc, 1));
  }
  auto par = b.create<scf::ParallelOp>(loc, lbs, ubs, steps,
                                       [&](OpBuilder &pb, Location pl, ValueRange ivs) { body(pb, pl, ivs); });
  par->setAttr("mir.lanes", b.getUnitAttr());
}

// Indices of an operand under a projected-permutation indexing map.
SmallVector<Value> indices(AffineMap map, ValueRange loopIvs) {
  SmallVector<Value> idx;
  for (AffineExpr e : map.getResults()) idx.push_back(loopIvs[cast<AffineDimExpr>(e).getPosition()]);
  return idx;
}

// A linalg op as one phase; `init` (if set) is the value of a fused fill of its
// single output, which the op overwrites entirely.
//
// Register blocking: when an input of a reduction does not depend on one of the
// parallel dimensions (u in u . Btil does not depend on the face l), that
// dimension (if short, and if the others still give every one of the element's
// `lanes` lanes a point) runs inside the thread: the input is loaded once per
// reduction step and the outputs along the dimension accumulate in registers. The
// summation order of each output is unchanged. This is how a hand kernel's thread
// keeps a column in registers and applies it for all faces.
constexpr int64_t kMaxBlock = 8;

LogicalResult lowerLinalg(linalg::LinalgOp op, Value init, int64_t lanes) {
  for (AffineMap m : op.getIndexingMapsArray())
    if (!m.isProjectedPermutation())
      return op->emitOpError("cellwise schedule: indexing maps must be projected permutations");
  OpBuilder b(op);
  Location loc = op.getLoc();
  SmallVector<int64_t> ranges = op.getStaticLoopRanges();
  SmallVector<utils::IteratorType> iters = op.getIteratorTypesArray();
  SmallVector<unsigned> par, red;
  for (unsigned d = 0; d < iters.size(); ++d) (iters[d] == utils::IteratorType::parallel ? par : red).push_back(d);
  SmallVector<OpOperand *> ins = op.getDpsInputOperands();
  SmallVector<OpOperand *> outs;
  for (OpOperand &o : op.getDpsInitsMutable()) outs.push_back(&o);
  if (outs.size() != 1) {
    // Several outputs: a pointwise op (a mir.flux with several results, e.g. the
    // three gradient outputs of a Galerkin Laplacian). One phase over its points
    // computes them all from one load of the inputs.
    if (!red.empty() || init)
      return op->emitOpError("cellwise schedule: several outputs only for an elementwise op");
    Block *region = op.getBlock();
    phase(b, loc, ranges, [&](OpBuilder &pb, Location pl, ValueRange ivs) {
      IRMapping map;
      for (OpOperand *in : ins) {
        Value v = in->get();
        map.map(region->getArgument(in->getOperandNumber()),
                isa<MemRefType>(v.getType())
                    ? pb.create<memref::LoadOp>(pl, v, indices(op.getMatchingIndexingMap(in), ivs)).getResult()
                    : v);
      }
      for (Operation &rop : region->without_terminator()) pb.clone(rop, map);
      for (unsigned k = 0; k < outs.size(); ++k)
        pb.create<memref::StoreOp>(pl, map.lookup(region->getTerminator()->getOperand(k)), outs[k]->get(),
                                   indices(op.getMatchingIndexingMap(outs[k]), ivs));
    });
    op->erase();
    return success();
  }

  int blocked = -1;   // the parallel loop dimension run inside the thread
  if (!red.empty()) {
    int best = 0;
    for (unsigned d : par) {
      if (ranges[d] > kMaxBlock || ranges[d] < 2)
        continue;
      int64_t rest = 1;
      for (unsigned e : par)
        if (e != d)
          rest *= ranges[e];
      if (rest < lanes)
        continue;
      int reuse = 0;
      for (OpOperand *in : ins)
        if (isa<MemRefType>(in->get().getType()) && !op.getMatchingIndexingMap(in).isFunctionOfDim(d))
          ++reuse;
      if (reuse > best) {
        best = reuse;
        blocked = (int)d;
      }
    }
  }
  SmallVector<unsigned> dims;
  SmallVector<int64_t> parShape;
  for (unsigned d : par)
    if ((int)d != blocked) {
      dims.push_back(d);
      parShape.push_back(ranges[d]);
    }
  if (parShape.empty())
    parShape.push_back(1);
  const int64_t B = blocked >= 0 ? ranges[blocked] : 1;
  Block *region = op.getBlock();
  BlockArgument outArg = region->getArgument(outs[0]->getOperandNumber());
  AffineMap outMap = op.getMatchingIndexingMap(outs[0]);

  phase(b, loc, parShape, [&](OpBuilder &pb, Location pl, ValueRange pivs) {
    SmallVector<Value> ivs(iters.size());
    for (unsigned k = 0; k < dims.size(); ++k) ivs[dims[k]] = pivs[k];
    auto at = [&](OpBuilder &ob, Location ol, int64_t kb) {
      if (blocked >= 0)
        ivs[blocked] = cidx(ob, ol, kb);
    };
    // One reduction step for the B outputs; inputs that do not depend on the
    // blocked dimension are loaded once.
    auto step = [&](OpBuilder &bb, Location bl, ValueRange accs) {
      SmallVector<Value> once(ins.size());
      for (unsigned i = 0; i < ins.size(); ++i) {
        Value v = ins[i]->get();
        AffineMap m = op.getMatchingIndexingMap(ins[i]);
        if (isa<MemRefType>(v.getType()) && (blocked < 0 || !m.isFunctionOfDim(blocked)))
          once[i] = bb.create<memref::LoadOp>(bl, v, indices(m, ivs));
      }
      SmallVector<Value> next;
      for (int64_t kb = 0; kb < B; ++kb) {
        at(bb, bl, kb);
        IRMapping map;
        for (unsigned i = 0; i < ins.size(); ++i) {
          Value v = ins[i]->get();
          Value x = once[i] ? once[i]
                            : (isa<MemRefType>(v.getType())
                                   ? bb.create<memref::LoadOp>(bl, v, indices(op.getMatchingIndexingMap(ins[i]), ivs))
                                         .getResult()
                                   : v);
          map.map(region->getArgument(ins[i]->getOperandNumber()), x);
        }
        map.map(outArg, accs[kb]);
        for (Operation &rop : region->without_terminator()) bb.clone(rop, map);
        next.push_back(map.lookup(region->getTerminator()->getOperand(0)));
      }
      return next;
    };
    // The output's current value is read only when the region uses it.
    SmallVector<Value> starts;
    for (int64_t kb = 0; kb < B; ++kb) {
      at(pb, pl, kb);
      starts.push_back(init ? init
                            : !outArg.use_empty()
                                ? pb.create<memref::LoadOp>(pl, outs[0]->get(), indices(outMap, ivs)).getResult()
                                : Value(pb.create<arith::ConstantOp>(pl, pb.getZeroAttr(outArg.getType()))));
    }
    std::function<SmallVector<Value>(OpBuilder &, Location, unsigned, ValueRange)> nest =
        [&](OpBuilder &nb, Location nl, unsigned k, ValueRange accs) -> SmallVector<Value> {
      if (k == red.size())
        return step(nb, nl, accs);
      auto loop = nb.create<scf::ForOp>(nl, cidx(nb, nl, 0), cidx(nb, nl, ranges[red[k]]), cidx(nb, nl, 1), accs,
                                        [&](OpBuilder &lb, Location ll, Value iv, ValueRange a) {
                                          ivs[red[k]] = iv;
                                          lb.create<scf::YieldOp>(ll, nest(lb, ll, k + 1, a));
                                        });
      return SmallVector<Value>(loop.getResults());
    };
    SmallVector<Value> results = nest(pb, pl, 0, starts);
    for (int64_t kb = 0; kb < B; ++kb) {
      at(pb, pl, kb);
      pb.create<memref::StoreOp>(pl, results[kb], outs[0]->get(), indices(outMap, ivs));
    }
  });
  op->erase();
  return success();
}

LogicalResult lowerFill(linalg::FillOp fill) {
  OpBuilder b(fill);
  auto type = cast<MemRefType>(fill.getOutputs()[0].getType());
  Value value = fill.getInputs()[0], buffer = fill.getOutputs()[0];
  phase(b, fill.getLoc(), type.getShape(),
        [&](OpBuilder &pb, Location pl, ValueRange ivs) { pb.create<memref::StoreOp>(pl, value, buffer, ivs); });
  fill->erase();
  return success();
}

void lowerCopy(OpBuilder &b, Location loc, Value src, Value dst) {
  phase(b, loc, cast<MemRefType>(dst.getType()).getShape(), [&](OpBuilder &pb, Location pl, ValueRange ivs) {
    pb.create<memref::StoreOp>(pl, pb.create<memref::LoadOp>(pl, src, ivs), dst, ivs);
  });
}

struct CellwiseSchedulePass : public PassWrapper<CellwiseSchedulePass, OperationPass<ModuleOp>> {
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(CellwiseSchedulePass)

  CellwiseSchedulePass() = default;
  CellwiseSchedulePass(const CellwiseSchedulePass &other) : PassWrapper(other) {}
  Option<int64_t> lanes{*this, "lanes", llvm::cl::desc("lanes per element (0: n^2, at most 32)"), llvm::cl::init(0)};
  Option<int64_t> warps{*this, "warps", llvm::cl::desc("warps per block"), llvm::cl::init(8)};

  StringRef getArgument() const final { return "mir-cellwise-schedule"; }
  StringRef getDescription() const final {
    return "Schedule a small-element operator as per-element phases (element loop, staged "
           "field, one phase per linalg op with fused fills, scratch temporaries)";
  }
  void getDependentDialects(DialectRegistry &registry) const override {
    registry.insert<arith::ArithDialect, linalg::LinalgDialect, memref::MemRefDialect, scf::SCFDialect>();
  }

  LogicalResult schedule(func::FuncOp func) {
    MLIRContext *ctx = func.getContext();
    Block &entry = func.getBody().front();
    Location loc = func.getLoc();
    SmallVector<unsigned> marked;
    for (unsigned i = 0; i < func.getNumArguments(); ++i)
      if (func.getArgAttr(i, kElementAttr))
        marked.push_back(i);
    if (marked.empty())
      return success();

    // 1. Batch dimension and element loop.
    SmallVector<Type> types(func.getArgumentTypes());
    SmallVector<MemRefType> elementTypes(func.getNumArguments());
    for (unsigned i : marked) {
      auto t = cast<MemRefType>(types[i]);
      elementTypes[i] = t;
      SmallVector<int64_t> shape{ShapedType::kDynamic};
      llvm::append_range(shape, t.getShape());
      types[i] = MemRefType::get(shape, t.getElementType());
      entry.getArgument(i).setType(types[i]);
    }
    func.setFunctionType(FunctionType::get(ctx, types, func.getResultTypes()));
    OpBuilder b = OpBuilder::atBlockBegin(&entry);
    Value count = b.create<memref::DimOp>(loc, entry.getArgument(marked[0]), 0);
    auto loop = b.create<scf::ParallelOp>(loc, ValueRange{cidx(b, loc, 0)}, ValueRange{count},
                                          ValueRange{cidx(b, loc, 1)});
    loop->setAttr("mir.elements", b.getUnitAttr());
    Block *body = loop.getBody();
    // Everything after the loop but the return moves into it.
    SmallVector<Operation *> move;
    for (Operation &op : llvm::make_range(std::next(Block::iterator(loop)), Block::iterator(entry.getTerminator())))
      move.push_back(&op);
    for (Operation *op : move) op->moveBefore(body->getTerminator());
    OpBuilder vb = OpBuilder::atBlockBegin(body);
    Value w = loop.getInductionVars()[0];
    DenseMap<unsigned, Value> views;
    for (unsigned i : marked) {
      Value arg = entry.getArgument(i);
      auto batched = cast<MemRefType>(arg.getType());
      SmallVector<OpFoldResult> offs{OpFoldResult(w)}, sizes{vb.getIndexAttr(1)}, strides{vb.getIndexAttr(1)};
      for (int64_t d : elementTypes[i].getShape()) {
        offs.push_back(vb.getIndexAttr(0));
        sizes.push_back(vb.getIndexAttr(d));
        strides.push_back(vb.getIndexAttr(1));
      }
      auto viewType = cast<MemRefType>(memref::SubViewOp::inferRankReducedResultType(
          elementTypes[i].getShape(), batched, offs, sizes, strides));
      Value view = vb.create<memref::SubViewOp>(loc, viewType, arg, offs, sizes, strides);
      arg.replaceUsesWithIf(view, [&](OpOperand &u) {
        return u.getOwner() != view.getDefiningOp() && loop->isProperAncestor(u.getOwner());
      });
      views[i] = view;
      func.removeArgAttr(i, kElementAttr);
    }
    mir::refreshViewTypes(func.getOperation());

    // 2. Staging in shared scratch.
    for (unsigned i : marked) {
      auto stage = func.getArgAttrOfType<StringAttr>(i, kStageAttr);
      if (!stage)
        continue;
      Value view = views[i];
      OpBuilder sb(view.getContext());
      sb.setInsertionPointAfter(view.getDefiningOp());
      auto scratch = sb.create<memref::AllocaOp>(loc, MemRefType::get(elementTypes[i].getShape(),
                                                                       elementTypes[i].getElementType()));
      scratch->setAttr("mir.scratch", sb.getUnitAttr());
      view.replaceUsesWithIf(scratch, [&](OpOperand &u) { return u.getOwner() != scratch.getOperation(); });
      if (stage.getValue() == "in") {
        sb.setInsertionPointAfter(scratch);
        lowerCopy(sb, loc, view, scratch);
      } else {
        OpBuilder eb(body->getTerminator());
        lowerCopy(eb, loc, scratch, view);
      }
      func.removeArgAttr(i, kStageAttr);
    }
    mir::refreshViewTypes(func.getOperation());   // views of a staged buffer lose the dynamic offset

    // 3. Unroll the short constant loops.
    IRRewriter rewriter(ctx);
    SmallVector<scf::ForOp> fors;
    body->walk([&](scf::ForOp f) { fors.push_back(f); });
    for (scf::ForOp f : fors) {
      auto lb = getConstantIntValue(f.getLowerBound()), ub = getConstantIntValue(f.getUpperBound()),
           st = getConstantIntValue(f.getStep());
      if (!lb || !ub || !st || f->getParentOfType<scf::ParallelOp>() != loop)
        continue;
      const int64_t trip = (*ub - *lb + *st - 1) / *st;
      if (trip <= 1) {
        if (failed(f.promoteIfSingleIteration(rewriter)))
          return f.emitOpError("cellwise schedule: cannot inline a one-trip loop");
      } else if (failed(loopUnrollByFactor(f, (uint64_t)trip))) {
        return f.emitOpError("cellwise schedule: cannot unroll");
      }
    }

    // Shape: n^2 lanes per element for an element of n^3 nodes.
    int64_t L = lanes;
    if (L == 0) {
      int64_t nodes = 1;
      for (int64_t d : elementTypes[marked[0]].getShape()) nodes *= d;
      const int64_t n = (int64_t)std::llround(std::cbrt((double)nodes));
      L = 1;
      while (L * 2 <= std::min<int64_t>(n * n, 32)) L *= 2;
    }

    // 4. Linalg ops, fills and copies into phases.
    SmallVector<Operation *> ops;
    for (Operation &op : *body) ops.push_back(&op);
    DenseMap<Value, linalg::FillOp> pendingFill;   // root -> fill not yet run
    auto materialize = [&](Value root) -> LogicalResult {
      auto it = pendingFill.find(root);
      if (it == pendingFill.end())
        return success();
      linalg::FillOp fill = it->second;
      pendingFill.erase(it);
      return lowerFill(fill);
    };
    auto touches = [&](Operation *op, function_ref<LogicalResult(Value)> f) -> LogicalResult {
      for (Value v : op->getOperands())
        if (isa<MemRefType>(v.getType()) && failed(f(rootOf(v))))
          return failure();
      return success();
    };
    for (Operation *op : ops) {
      if (auto fill = dyn_cast<linalg::FillOp>(op)) {
        Value out = fill.getOutputs()[0];
        if (failed(materialize(rootOf(out))))
          return failure();
        if (rootOf(out) == out && out.getDefiningOp<memref::AllocOp>()) {   // a whole temporary
          pendingFill[out] = fill;
          continue;
        }
        if (failed(lowerFill(fill)))
          return failure();
        continue;
      }
      if (auto lop = dyn_cast<linalg::LinalgOp>(op)) {
        Value init;
        Value out = lop.getDpsInits()[0];
        const bool reduction = lop.getNumReductionLoops() > 0;
        auto it = pendingFill.find(out);
        if (reduction && it != pendingFill.end()) {   // the fill's value starts the accumulator
          init = it->second.getInputs()[0];
          linalg::FillOp fill = it->second;
          pendingFill.erase(it);
          fill->erase();
        }
        if (failed(touches(op, materialize)))
          return failure();
        if (!isa<linalg::GenericOp>(op)) {
          rewriter.setInsertionPoint(op);
          FailureOr<linalg::GenericOp> g = linalg::generalizeNamedOp(rewriter, lop);
          if (failed(g))
            return op->emitOpError("cellwise schedule: cannot generalize");
          lop = cast<linalg::LinalgOp>(g->getOperation());
        }
        if (failed(lowerLinalg(lop, init, L)))
          return failure();
        continue;
      }
      if (auto copy = dyn_cast<memref::CopyOp>(op)) {
        if (failed(touches(op, materialize)))
          return failure();
        OpBuilder cb(copy);
        lowerCopy(cb, copy.getLoc(), copy.getSource(), copy.getTarget());
        copy->erase();
        continue;
      }
      if (op->getNumRegions() > 0 && !isa<scf::ParallelOp>(op))
        return op->emitOpError("cellwise schedule: unexpected op with regions in the element body");
      if (failed(touches(op, materialize)))
        return failure();
    }
    for (auto &[root, fill] : llvm::make_early_inc_range(pendingFill)) {
      (void)root;
      if (failed(lowerFill(fill)))
        return failure();
    }

    // 5. Temporaries into shared scratch.
    SmallVector<memref::AllocOp> allocs;
    body->walk([&](memref::AllocOp a) { allocs.push_back(a); });
    for (memref::AllocOp a : allocs) {
      OpBuilder ab(a);
      auto s = ab.create<memref::AllocaOp>(a.getLoc(), a.getType());
      s->setAttr("mir.scratch", ab.getUnitAttr());
      a.replaceAllUsesWith(s.getResult());
      a->erase();
    }
    SmallVector<memref::DeallocOp> deallocs;
    body->walk([&](memref::DeallocOp d) { deallocs.push_back(d); });
    for (memref::DeallocOp d : deallocs) d->erase();

    loop->setAttr("mir.lanes_per_element", b.getI64IntegerAttr(L));
    loop->setAttr("mir.warps_per_block", b.getI64IntegerAttr(warps));
    return success();
  }

  void runOnOperation() override {
    for (auto func : getOperation().getOps<func::FuncOp>())
      if (func->hasAttr("mir.kernel") && failed(schedule(func)))
        return signalPassFailure();
  }
};

}  // namespace

namespace mir {
void registerCellwiseSchedulePass() { PassRegistration<CellwiseSchedulePass>(); }
}  // namespace mir
