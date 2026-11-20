// --mir-mark-kernel: say which func.func is the kernel and what each argument is, by position,
// before any pass that keys on it; optionally rename the kernel. Bufferization renames the
// arguments, so roles are given by index (negative: from the end).
//
//   --mir-mark-kernel="func=bp3_apply elements=0,-2,-1 constants=1,2,3,4 reflect=1,1,-1,-1"
//       u, G, y per element ({mir.element}); B, Bt, D, Dt shared, in constant memory, mirror-
//       symmetric (+1) or antisymmetric (-1) for --mir-even-odd; the func gets {mir.kernel}.
//   stage-in=0 stage-out=-1    the element's field staged into / its result out of shared
//                              memory ({mir.stage = "in"/"out"}, the cell-wise schedule).
//   --mir-mark-kernel="func=bp3_apply rename=bp3_apply_p3"
//       the kernel and the tables named after it (bp3_apply__m0 -> bp3_apply_p3__m0).

#include "mir/MirPasses.h"

#include "mlir/Dialect/Func/IR/FuncOps.h"
#include "mlir/IR/Builders.h"
#include "mlir/IR/BuiltinOps.h"
#include "mlir/IR/SymbolTable.h"
#include "mlir/Pass/Pass.h"

using namespace mlir;

namespace {

struct MarkKernelPass : public PassWrapper<MarkKernelPass, OperationPass<ModuleOp>> {
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(MarkKernelPass)

  MarkKernelPass() = default;
  MarkKernelPass(const MarkKernelPass &other) : PassWrapper(other) {}

  StringRef getArgument() const final { return "mir-mark-kernel"; }
  StringRef getDescription() const final {
    return "Mark the kernel function and its argument roles by position; optionally rename it";
  }

  Option<std::string> func{*this, "func", llvm::cl::desc("the kernel func.func")};
  Option<std::string> rename{*this, "rename", llvm::cl::desc("new name for the kernel and its tables")};
  ListOption<int64_t> elements{*this, "elements", llvm::cl::desc("per-element arguments ({mir.element})")};
  ListOption<int64_t> constants{*this, "constants", llvm::cl::desc("shared constant-memory arguments")};
  ListOption<int64_t> reflect{*this, "reflect", llvm::cl::desc("mir.reflect of each constant (+1 / -1)")};
  Option<int64_t> stageIn{*this, "stage-in", llvm::cl::desc("argument staged into shared memory"),
                          llvm::cl::init(kNone)};
  Option<int64_t> stageOut{*this, "stage-out", llvm::cl::desc("argument staged out of shared memory"),
                           llvm::cl::init(kNone)};

  static constexpr int64_t kNone = INT64_MIN;

  void runOnOperation() override {
    ModuleOp module = getOperation();
    auto f = module.lookupSymbol<func::FuncOp>(func);
    if (!f) {
      module.emitError("--mir-mark-kernel: no func.func @") << func;
      return signalPassFailure();
    }
    const int64_t n = f.getNumArguments();
    auto at = [&](int64_t i) -> int64_t { return i < 0 ? n + i : i; };
    for (int64_t i : llvm::concat<int64_t>(elements, constants))
      if (at(i) < 0 || at(i) >= n) {
        f.emitError("--mir-mark-kernel: argument ") << i << " out of range";
        return signalPassFailure();
      }
    if (!reflect.empty() && reflect.size() != constants.size()) {
      f.emitError("--mir-mark-kernel: reflect needs one value per constant");
      return signalPassFailure();
    }
    Builder b(&getContext());
    if (!elements.empty() || !constants.empty()) {
      f->setAttr("mir.kernel", b.getUnitAttr());
      for (int64_t i : elements)
        f.setArgAttr(at(i), "mir.element", b.getUnitAttr());
      for (auto [j, i] : llvm::enumerate(constants)) {
        f.setArgAttr(at(i), "mir.constant", b.getUnitAttr());
        if (!reflect.empty())
          f.setArgAttr(at(i), "mir.reflect", b.getI64IntegerAttr(reflect[j]));
      }
      if (stageIn != kNone)
        f.setArgAttr(at(stageIn), "mir.stage", b.getStringAttr("in"));
      if (stageOut != kNone)
        f.setArgAttr(at(stageOut), "mir.stage", b.getStringAttr("out"));
    }
    if (!rename.empty()) {
      // The kernel's tables are named <kernel>__<suffix>; they follow the kernel.
      const std::string prefix = func + "__";
      for (Operation &op : llvm::make_early_inc_range(module.getBody()->getOperations())) {
        auto sym = dyn_cast<SymbolOpInterface>(op);
        if (!sym || !sym.getName().starts_with(prefix))
          continue;
        auto to = b.getStringAttr(rename + "__" + sym.getName().drop_front(prefix.size()).str());
        if (failed(SymbolTable::replaceAllSymbolUses(&op, to, module)))
          return signalPassFailure();
        SymbolTable::setSymbolName(&op, to);
      }
      if (failed(SymbolTable::replaceAllSymbolUses(f, b.getStringAttr(rename), module)))
        return signalPassFailure();
      SymbolTable::setSymbolName(f, b.getStringAttr(rename));
    }
  }
};

}  // namespace

namespace mir {
void registerMarkKernelPass() { PassRegistration<MarkKernelPass>(); }
}  // namespace mir
