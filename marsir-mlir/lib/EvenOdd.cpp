// --mir-even-odd: halve the multiply-adds of a 1D sweep with a mirror-symmetric matrix.
//
// After --mir-cellwise-schedule a sweep is unrolled per line: y[i] = init[i] +
// sum_j M[i][j] x[j], each term a math.fma whose matrix operand is a load at constant
// indices from a {mir.constant} global. GLL nodes and Gauss points are symmetric about
// the element centre, so M[r-1-i][c-1-j] = s M[i][j] (s = +1 for B and Bt, -1 for the
// derivative D and Dt; the global carries s as mir.reflect). With e[j] = x[j] + x[c-1-j]
// and o[j] = x[j] - x[c-1-j], j < c/2, and A, S the half sums and differences of
// M[i][j] and M[i][c-1-j], i < r/2:
//
//   ye = A e (+ M[i][c/2] x[c/2] if c is odd),  yo = S o
//   y[i] = ye + yo,  y[r-1-i] = s (ye - yo)
//   y[r/2] (r odd) = sum_j M[r/2][j] (s > 0 ? e[j] : o[j]) (+ M[r/2][c/2] x[c/2])
//
// about half the fma, plus c adds for e, o and r for the recombination. The global
// then holds the folded coefficients [A, S, column c/2, row r/2, centre], which the
// host forms from M (mars_marsir_fold.hpp); <global>__fold = {r, c, s} tells it to.
// A global is folded only when every read of it is part of a complete sweep.

#include "mir/MirPasses.h"

#include "mlir/Dialect/Arith/IR/Arith.h"
#include "mlir/Dialect/Math/IR/Math.h"
#include "mlir/Dialect/MemRef/IR/MemRef.h"
#include "mlir/Dialect/Utils/StaticValueUtils.h"
#include "mlir/IR/BuiltinOps.h"
#include "mlir/IR/SymbolTable.h"
#include "mlir/Pass/Pass.h"

#include "llvm/ADT/MapVector.h"
#include "llvm/ADT/SetVector.h"

using namespace mlir;

namespace {

struct Term {
  math::FmaOp fma;
  memref::LoadOp load;
  int64_t i, j;
  Value x;
};

// One output of a sweep: init + sum_j M[i][j] x[j], each j once, in fma order.
struct Chain {
  int64_t i;
  Value init;
  SmallVector<Term> terms;   // in chain order
  SmallVector<Value> x;      // x[j]
};

std::optional<Term> asTerm(math::FmaOp fma, StringRef global) {
  for (int k = 0; k < 2; ++k) {
    auto load = fma->getOperand(k).getDefiningOp<memref::LoadOp>();
    if (!load || load.getIndices().size() != 2) continue;
    auto get = load.getMemRef().getDefiningOp<memref::GetGlobalOp>();
    if (!get || get.getName() != global) continue;
    auto i = getConstantIntValue(load.getIndices()[0]), j = getConstantIntValue(load.getIndices()[1]);
    if (!i || !j) continue;
    return Term{fma, load, *i, *j, fma->getOperand(1 - k)};
  }
  return std::nullopt;
}

bool isZero(Value v) {
  FloatAttr f;
  auto c = v.getDefiningOp<arith::ConstantOp>();
  return c && (f = dyn_cast<FloatAttr>(c.getValue())) && f.getValue().isZero();
}

struct EvenOddPass : public PassWrapper<EvenOddPass, OperationPass<ModuleOp>> {
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(EvenOddPass)

  StringRef getArgument() const final { return "mir-even-odd"; }
  StringRef getDescription() const final {
    return "Even-odd split of unrolled 1D sweeps over mirror-symmetric constant matrices";
  }
  void getDependentDialects(DialectRegistry &registry) const override {
    registry.insert<arith::ArithDialect, math::MathDialect, memref::MemRefDialect>();
  }

  void runOnOperation() override {
    ModuleOp module = getOperation();
    for (auto global : llvm::to_vector(module.getOps<memref::GlobalOp>())) {
      auto reflect = global->getAttrOfType<IntegerAttr>("mir.reflect");
      auto type = global.getType();
      if (!reflect || type.getRank() != 2 || !type.hasStaticShape()) continue;
      fold(module, global, type.getDimSize(0), type.getDimSize(1), reflect.getInt());
    }
  }

  // Chains over `global`, or nothing when a read of it is not part of one.
  std::optional<SmallVector<Chain>> chains(ModuleOp module, StringRef name, int64_t r, int64_t c) {
    llvm::DenseMap<Operation *, Term> terms;
    bool other = false;
    // The table changes layout, so every access must be one this pass rewrites: the symbol
    // only reached through get_global, each of those only read by loads (an alias such as
    // reinterpret_cast would read folded coefficients as M), each load used.
    auto uses = SymbolTable::getSymbolUses(StringAttr::get(module.getContext(), name), module);
    if (!uses) return std::nullopt;
    for (const SymbolTable::SymbolUse &use : *uses) {
      auto get = dyn_cast<memref::GetGlobalOp>(use.getUser());
      if (!get) return std::nullopt;
      for (Operation *user : get->getUsers())
        if (!isa<memref::LoadOp>(user)) return std::nullopt;
    }
    module.walk([&](memref::LoadOp load) {
      auto get = load.getMemRef().getDefiningOp<memref::GetGlobalOp>();
      if (!get || get.getName() != name) return;
      if (load->use_empty()) { other = true; return; }
      for (Operation *user : load->getUsers()) {
        auto fma = dyn_cast<math::FmaOp>(user);
        std::optional<Term> t = fma ? asTerm(fma, name) : std::nullopt;
        if (!t || t->load != load) { other = true; continue; }
        terms[fma] = *t;
      }
    });
    if (other) return std::nullopt;
    SmallVector<Chain> out;
    for (auto &[op, t] : terms) {
      Operation *prev = t.fma.getOperand(2).getDefiningOp();
      if (prev && terms.count(prev)) continue;   // not the first term of its chain
      Chain ch{t.i, t.fma.getOperand(2), {}, SmallVector<Value>(c)};
      Term cur = t;
      while (true) {
        if (cur.i != ch.i || cur.j < 0 || cur.j >= c || ch.x[cur.j]) return std::nullopt;
        ch.x[cur.j] = cur.x;
        ch.terms.push_back(cur);
        Value res = cur.fma.getResult();
        if (!res.hasOneUse()) break;
        Operation *next = *res.getUsers().begin();
        auto it = terms.find(next);
        if (it == terms.end() || it->second.fma.getOperand(2) != res) break;
        cur = it->second;
      }
      if ((int64_t)ch.terms.size() != c) return std::nullopt;
      out.push_back(std::move(ch));
    }
    size_t covered = 0;
    for (Chain &ch : out) covered += ch.terms.size();
    if (covered != terms.size() || out.empty()) return std::nullopt;
    for (Chain &ch : out)
      if (ch.i < 0 || ch.i >= r) return std::nullopt;
    return out;
  }

  void fold(ModuleOp module, memref::GlobalOp global, int64_t r, int64_t c, int64_t s) {
    auto found = chains(module, global.getSymName(), r, c);
    if (!found) return;
    // Sweeps: chains in one block reading the same x, one per row.
    llvm::MapVector<std::pair<Block *, Value>, SmallVector<Chain *>> sweeps;
    for (Chain &ch : *found) sweeps[{ch.terms.front().fma->getBlock(), ch.x[0]}].push_back(&ch);
    for (auto &[key, list] : sweeps) {
      if ((int64_t)list.size() != r) return;
      SmallVector<Chain *> byRow(r, nullptr);
      for (Chain *ch : list) {
        if (byRow[ch->i] || ch->x != list.front()->x) return;
        byRow[ch->i] = ch;
      }
      list = byRow;
    }

    const int64_t hr = r / 2, hc = c / 2, mr = r % 2, mc = c % 2;
    const int64_t offS = hr * hc, offCol = 2 * hr * hc, offRow = offCol + hr * mc, offMid = offRow + mr * hc;
    const int64_t size = offMid + mr * mc;
    auto oldType = global.getType();
    auto newType = MemRefType::get({size}, oldType.getElementType(), MemRefLayoutAttrInterface(), oldType.getMemorySpace());
    global.setType(newType);
    OpBuilder gb(global);
    gb.setInsertionPointAfter(global);
    auto recType = MemRefType::get({3}, gb.getI64Type());
    gb.create<memref::GlobalOp>(global.getLoc(), gb.getStringAttr((global.getSymName() + "__fold").str()),
                                StringAttr(), TypeAttr::get(recType),
                                DenseElementsAttr::get(RankedTensorType::get({3}, gb.getI64Type()), ArrayRef<int64_t>{r, c, s}),
                                gb.getUnitAttr(), IntegerAttr());

    llvm::SetVector<Operation *> dead;
    module.walk([&](memref::GetGlobalOp get) {
      if (get.getName() != global.getSymName()) return;
      get.getResult().setType(newType);
    });
    for (auto &[key, rows] : sweeps) {
      // After the last term of the sweep, where every x and init is defined.
      Operation *last = rows.front()->terms.back().fma;
      for (Chain *ch : rows)
        for (Term &t : ch->terms)
          if (last->isBeforeInBlock(t.fma)) last = t.fma;
      for (Chain *ch : rows)
        for (Operation *user : ch->terms.back().fma->getUsers())
          if (user->getBlock() == last->getBlock() && user->isBeforeInBlock(last)) {
            global.emitError("even-odd: a sweep result is used before the sweep ends");
            return signalPassFailure();
          }
      memref::LoadOp anyLoad = rows.front()->terms.front().load;
      Value table = anyLoad.getMemRef();
      Location loc = last->getLoc();
      OpBuilder b(last->getContext());
      b.setInsertionPointAfter(last);
      ArrayRef<Value> x = rows.front()->x;
      auto coef = [&](int64_t k) {
        Value idx = b.create<arith::ConstantIndexOp>(loc, k);
        return b.create<memref::LoadOp>(loc, table, ValueRange{idx}).getResult();
      };
      // sum_k coef[k] v[k], k over `idx`
      auto dot = [&](ArrayRef<std::pair<int64_t, Value>> terms) {
        Value acc = b.create<arith::MulFOp>(loc, coef(terms[0].first), terms[0].second);
        for (auto [k, v] : terms.drop_front()) acc = b.create<math::FmaOp>(loc, coef(k), v, acc);
        return acc;
      };
      auto plusInit = [&](Value init, Value v) { return isZero(init) ? v : b.create<arith::AddFOp>(loc, init, v).getResult(); };
      SmallVector<Value> e, o;
      for (int64_t j = 0; j < hc; ++j) {
        e.push_back(b.create<arith::AddFOp>(loc, x[j], x[c - 1 - j]));
        o.push_back(b.create<arith::SubFOp>(loc, x[j], x[c - 1 - j]));
      }
      SmallVector<Value> y(r);
      for (int64_t i = 0; i < hr; ++i) {
        SmallVector<std::pair<int64_t, Value>> te, to;
        for (int64_t j = 0; j < hc; ++j) {
          te.push_back({i * hc + j, e[j]});
          to.push_back({offS + i * hc + j, o[j]});
        }
        if (mc) te.push_back({offCol + i, x[hc]});
        Value ye = hc || mc ? dot(te) : Value(), yo = hc ? dot(to) : Value();
        Value top = yo ? b.create<arith::AddFOp>(loc, ye, yo).getResult() : ye;
        Value bot = yo ? (s > 0 ? b.create<arith::SubFOp>(loc, ye, yo).getResult() : b.create<arith::SubFOp>(loc, yo, ye).getResult())
                       : (s > 0 ? ye : b.create<arith::NegFOp>(loc, ye).getResult());
        y[i] = plusInit(rows[i]->init, top);
        y[r - 1 - i] = plusInit(rows[r - 1 - i]->init, bot);
      }
      if (mr) {
        SmallVector<std::pair<int64_t, Value>> tm;
        for (int64_t j = 0; j < hc; ++j) tm.push_back({offRow + j, s > 0 ? e[j] : o[j]});
        if (mc) tm.push_back({offMid, x[hc]});
        y[hr] = plusInit(rows[hr]->init, dot(tm));
      }
      for (int64_t i = 0; i < r; ++i) {
        rows[i]->terms.back().fma.getResult().replaceAllUsesWith(y[i]);
        for (Term &t : llvm::reverse(rows[i]->terms)) {
          dead.insert(t.fma);
          dead.insert(t.load);
        }
      }
    }
    // A load CSE shared between sweeps dies only with its last fma.
    SmallVector<Operation *> pending(dead.begin(), dead.end());
    for (bool progress = true; progress;) {
      progress = false;
      for (Operation *&op : pending)
        if (op && op->use_empty()) {
          op->erase();
          op = nullptr;
          progress = true;
        }
    }
  }
};

}  // namespace

namespace mir {
void registerEvenOddPass() { PassRegistration<EvenOddPass>(); }
}  // namespace mir
