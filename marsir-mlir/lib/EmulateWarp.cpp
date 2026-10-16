// --mir-emulate-warp: TEST ONLY. Rewrite the warp primitives of a kernel body
// (still a func.func, before --mir-gpu-wrap) into calls to a host runtime, so the
// SAME body the GPU gets runs on the CPU with one thread per lane
// (tools/emu_hl_mma.cpp). Per-thread memory stays per thread -- allocas live on
// each thread's own stack -- so a fragment stored piecewise into private memory
// fails here the way it failed on the GPU, and a missing barrier can show up as
// a wrong or unstable result.
//
//   gpu.thread_id x       -> call @mir_emu_lane()        (y, z -> 0: 32x1x1 block)
//   gpu.block_id x        -> call @mir_emu_block()       (y, z -> 0)
//   gpu.grid_dim x        -> call @mir_emu_grid()        (y, z -> 1)
//   gpu.block_dim x       -> 32                          (y, z -> 1)
//   gpu.shuffle f64/i32   -> call @mir_emu_shfl(v, src) / @mir_emu_shfl_i32(v, src),
//                            src = offset (idx) or lane ^ offset (xor)
//   gpu.barrier           -> call @mir_emu_barrier()
//   nvgpu.mma.sync m8n8k4 -> call @mir_emu_mma(a, b, c0, c1), @mir_emu_mma_hi()
//   {mir.workgroup} arg   -> memref.get_global of one process-wide buffer. The
//                            harness runs one block at a time, so one copy is
//                            the block's. It is filled with NaN at entry: GPU
//                            shared memory starts undefined, and a read before
//                            any write must not pass by reading a lucky zero.
//                            (Integer buffers get all bits set.)

#include "mir/MirPasses.h"

#include "mlir/Dialect/Arith/IR/Arith.h"
#include "mlir/Dialect/Func/IR/FuncOps.h"
#include "mlir/Dialect/GPU/IR/GPUDialect.h"
#include "mlir/Dialect/Linalg/IR/Linalg.h"
#include "mlir/Dialect/MemRef/IR/MemRef.h"
#include "mlir/Dialect/NVGPU/IR/NVGPUDialect.h"
#include "mlir/Dialect/SCF/IR/SCF.h"
#include "mlir/Dialect/Vector/IR/VectorOps.h"
#include "mlir/IR/Builders.h"
#include "mlir/IR/BuiltinOps.h"
#include "mlir/Pass/Pass.h"

#include <limits>

using namespace mlir;

namespace {

constexpr StringLiteral kWorkgroupAttr = "mir.workgroup";

struct EmulateWarpPass
    : public PassWrapper<EmulateWarpPass, OperationPass<ModuleOp>> {
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(EmulateWarpPass)

  StringRef getArgument() const final { return "mir-emulate-warp"; }
  StringRef getDescription() const final {
    return "TEST ONLY: rewrite warp primitives into calls to a host runtime that "
           "runs one thread per lane";
  }
  void getDependentDialects(DialectRegistry &registry) const override {
    registry.insert<arith::ArithDialect, func::FuncDialect, linalg::LinalgDialect,
                    memref::MemRefDialect, scf::SCFDialect,
                    vector::VectorDialect>();
  }

  void runOnOperation() override {
    ModuleOp mod = getOperation();
    MLIRContext *ctx = &getContext();
    OpBuilder mb(ctx);
    mb.setInsertionPointToStart(mod.getBody());
    Type f64 = mb.getF64Type(), i32 = mb.getI32Type(), idx = mb.getIndexType();

    auto declare = [&](StringRef name, TypeRange in, TypeRange out) {
      if (mod.lookupSymbol<func::FuncOp>(name))
        return;
      auto f = mb.create<func::FuncOp>(mod.getLoc(), name,
                                       mb.getFunctionType(in, out));
      f.setPrivate();
    };
    declare("mir_emu_lane", {}, {idx});
    declare("mir_emu_block", {}, {idx});
    declare("mir_emu_grid", {}, {idx});
    declare("mir_emu_barrier", {}, {});
    declare("mir_emu_shfl", {f64, i32}, {f64});
    declare("mir_emu_shfl_i32", {i32, i32}, {i32});
    declare("mir_emu_mma", {f64, f64, f64, f64}, {f64});
    declare("mir_emu_mma_hi", {}, {f64});

    SmallVector<func::FuncOp> kernels;
    for (auto f : mod.getOps<func::FuncOp>())
      if (!f.isDeclaration())
        kernels.push_back(f);

    for (func::FuncOp fn : kernels) {
      SmallVector<Operation *> ops;
      fn.walk([&](Operation *op) {
        if (isa<gpu::ThreadIdOp, gpu::BlockIdOp, gpu::GridDimOp, gpu::BlockDimOp, gpu::ShuffleOp, gpu::BarrierOp,
                nvgpu::MmaSyncOp>(op))
          ops.push_back(op);
      });
      for (Operation *op : ops) {
        OpBuilder b(op);
        Location loc = op->getLoc();
        auto call = [&](StringRef name, TypeRange res, ValueRange args) {
          return b.create<func::CallOp>(loc, name, res, args);
        };
        if (auto t = dyn_cast<gpu::ThreadIdOp>(op)) {
          Value v = t.getDimension() == gpu::Dimension::x
                        ? call("mir_emu_lane", {idx}, {}).getResult(0)
                        : b.create<arith::ConstantIndexOp>(loc, 0).getResult();
          t.replaceAllUsesWith(v);
        } else if (auto t = dyn_cast<gpu::BlockIdOp>(op)) {
          Value v = t.getDimension() == gpu::Dimension::x
                        ? call("mir_emu_block", {idx}, {}).getResult(0)
                        : b.create<arith::ConstantIndexOp>(loc, 0).getResult();
          t.replaceAllUsesWith(v);
        } else if (auto t = dyn_cast<gpu::GridDimOp>(op)) {
          Value v = t.getDimension() == gpu::Dimension::x
                        ? call("mir_emu_grid", {idx}, {}).getResult(0)
                        : b.create<arith::ConstantIndexOp>(loc, 1).getResult();
          t.replaceAllUsesWith(v);
        } else if (auto t = dyn_cast<gpu::BlockDimOp>(op)) {
          t.replaceAllUsesWith(b.create<arith::ConstantIndexOp>(
                                    loc, t.getDimension() == gpu::Dimension::x ? 32 : 1)
                                   .getResult());
        } else if (auto s = dyn_cast<gpu::ShuffleOp>(op)) {
          const bool isF64 = s.getValue().getType().isF64();
          if ((s.getMode() != gpu::ShuffleMode::IDX && s.getMode() != gpu::ShuffleMode::XOR) ||
              (!isF64 && !s.getValue().getType().isInteger(32))) {
            s.emitOpError("mir-emulate-warp: only idx/xor shuffles of f64/i32 are emulated");
            return signalPassFailure();
          }
          Value src = s.getOffset();
          if (s.getMode() == gpu::ShuffleMode::XOR) {
            Value lane = b.create<arith::IndexCastOp>(loc, i32, call("mir_emu_lane", {idx}, {}).getResult(0));
            src = b.create<arith::XOrIOp>(loc, lane, src);
          }
          Value r = isF64 ? call("mir_emu_shfl", {f64}, {s.getValue(), src}).getResult(0)
                          : call("mir_emu_shfl_i32", {i32}, {s.getValue(), src}).getResult(0);
          s.getShuffleResult().replaceAllUsesWith(r);
          s.getValid().replaceAllUsesWith(
              b.create<arith::ConstantIntOp>(loc, 1, 1).getResult());
        } else if (isa<gpu::BarrierOp>(op)) {
          call("mir_emu_barrier", {}, {});
        } else if (auto m = dyn_cast<nvgpu::MmaSyncOp>(op)) {
          auto shape = m.getMmaShapeAsArray();
          auto cTy = cast<VectorType>(m.getMatrixC().getType());
          if (shape[0] != 8 || shape[1] != 8 || shape[2] != 4 ||
              !cTy.getElementType().isF64()) {
            m.emitOpError("mir-emulate-warp: only f64 m8n8k4 is emulated");
            return signalPassFailure();
          }
          auto at = [&](Value v, int64_t c) -> Value {
            return b.create<vector::ExtractOp>(loc, v, ArrayRef<int64_t>{0, c});
          };
          Value d0 = call("mir_emu_mma", {f64},
                          {at(m.getMatrixA(), 0), at(m.getMatrixB(), 0),
                           at(m.getMatrixC(), 0), at(m.getMatrixC(), 1)})
                         .getResult(0);
          Value d1 = call("mir_emu_mma_hi", {f64}, {}).getResult(0);
          Value d = b.create<arith::ConstantOp>(
              loc, cTy, DenseElementsAttr::get(cTy, b.getF64FloatAttr(0.0)));
          d = b.create<vector::InsertOp>(loc, d0, d, ArrayRef<int64_t>{0, 0});
          d = b.create<vector::InsertOp>(loc, d1, d, ArrayRef<int64_t>{0, 1});
          m.getResult().replaceAllUsesWith(d);
        }
        op->erase();
      }

      // Workgroup buffers become process-wide globals, poisoned by lane 0.
      Block &entry = fn.getBody().front();
      SmallVector<unsigned> wg;
      for (unsigned i = 0; i < fn.getNumArguments(); ++i)
        if (fn.getArgAttr(i, kWorkgroupAttr))
          wg.push_back(i);
      if (!wg.empty()) {
        OpBuilder b = OpBuilder::atBlockBegin(&entry);
        Location loc = fn.getLoc();
        Value lane = b.create<func::CallOp>(loc, "mir_emu_lane", TypeRange{idx},
                                            ValueRange{})
                         .getResult(0);
        Value isZero = b.create<arith::CmpIOp>(
            loc, arith::CmpIPredicate::eq, lane,
            b.create<arith::ConstantIndexOp>(loc, 0));
        auto ifOp = b.create<scf::IfOp>(loc, isZero, /*withElseRegion=*/false);
        OpBuilder tb = ifOp.getThenBodyBuilder();
        Value nan = tb.create<arith::ConstantOp>(
            loc, tb.getF64FloatAttr(std::numeric_limits<double>::quiet_NaN()));
        for (unsigned i : wg) {
          auto ty = dyn_cast<MemRefType>(fn.getArgument(i).getType());
          if (!ty || !ty.hasStaticShape()) {
            fn.emitOpError("mir-emulate-warp: a workgroup buffer must be a "
                           "static memref");
            return signalPassFailure();
          }
          std::string name =
              ("mir_emu_wg_" + fn.getName() + "_" + Twine(i)).str();
          OpBuilder gb(ctx);
          gb.setInsertionPoint(fn);
          gb.create<memref::GlobalOp>(loc, name, gb.getStringAttr("private"), ty,
                                      /*initial_value=*/gb.getUnitAttr(),
                                      /*constant=*/false,
                                      /*alignment=*/gb.getI64IntegerAttr(64));
          OpBuilder eb(ctx);
          eb.setInsertionPoint(ifOp);
          Value g = eb.create<memref::GetGlobalOp>(loc, ty, name);
          fn.getArgument(i).replaceAllUsesWith(g);
          Type et = ty.getElementType();
          Value poison = nan;
          if (isa<IntegerType>(et))
            poison = tb.create<arith::ConstantIntOp>(loc, -1, et.getIntOrFloatBitWidth());
          else if (!et.isF64())
            poison = tb.create<arith::ConstantOp>(
                loc, tb.getFloatAttr(et, std::numeric_limits<double>::quiet_NaN()));
          tb.create<linalg::FillOp>(loc, ValueRange{poison}, ValueRange{g});
        }
        b.setInsertionPointAfter(ifOp);
        b.create<func::CallOp>(loc, "mir_emu_barrier", TypeRange{}, ValueRange{});
        for (unsigned i : llvm::reverse(wg))
          fn.eraseArgument(i);
      }
      fn->setAttr("llvm.emit_c_interface", UnitAttr::get(ctx));
    }
  }
};

}  // namespace

namespace mir {
void registerEmulateWarpPass() { PassRegistration<EmulateWarpPass>(); }
}  // namespace mir
