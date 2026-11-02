// --mir-lower-cellwise: lower mir.cellwise (one per-copy phase per element, see
// lowerOne(CellwiseOp)) and mir.dss: the latter to the staged gather of the hand-written MARS
// kernel (unstructured_dss_kernel in
// backend/distributed/unstructured/solvers/mars_cellwise_topology.hpp), as plain
// loops that say nothing about the GPU:
//
//   scf.parallel (w) over the listed elements          {mir.elements}
//     e = elements[w], stage = per-element scratch     {mir.scratch}
//     scf.parallel over the n^3 nodes:  stage <- own copies      {mir.lanes}
//     scf.parallel over the face nodes: stage += the other copy  {mir.lanes}
//     scf.parallel over the edge and corner nodes:
//                                       stage <- sum of the star {mir.lanes}
//     scf.parallel over the n^3 nodes:  out <- stage             {mir.lanes}
//
// Iterations of one phase touch disjoint stage slots, so each phase is a true
// parallel loop; the phases themselves run in order. --mir-warp-per-element maps
// this onto one warp per element. Before that mapping the loops also run as they
// are on the CPU.
//
// The reference-hex tables (face nodes, face orientation maps, edge and corner
// nodes) come from MARS's own mars_cellwise_hex.hpp, the header every hand
// kernel uses, so the node numbering has one definition. They become constant
// globals; in the hand kernel each thread block builds them in shared memory.
//
// Summation order is the hand kernel's, so the results are bit-identical: a
// face node is own + other, a star sum starts from 0 and adds the entries in
// order, four loads in flight at a time.
//
// With an epilogue (the hand kernels' `visit`), two more phases:
//   scf.parallel over the 32 entity bits, reducing with OR: the element's
//     Dirichlet mask (mapped onto lanes this is the hand kernel's ballot)
//   the write-back runs the epilogue on every copy, stores its per-copy results
//     and reduces its terms over the element with scf.reduce; the element's
//     sums go to partials[w * reductions + j].

#include "mir/MirOps.h"
#include "mir/MirPasses.h"

#include "backend/distributed/unstructured/solvers/mars_cellwise_hex.hpp"

#include "mlir/Dialect/Arith/IR/Arith.h"
#include "mlir/Dialect/Func/IR/FuncOps.h"
#include "mlir/Dialect/MemRef/IR/MemRef.h"
#include "mlir/Dialect/SCF/IR/SCF.h"
#include "mlir/IR/Builders.h"
#include "mlir/IR/IRMapping.h"
#include "mlir/IR/BuiltinOps.h"
#include "mlir/Pass/Pass.h"

#include <functional>
#include <string>
#include <vector>

using namespace mlir;

namespace {

namespace cw = mars::cellwise;

constexpr int kMaxDegree = 8;   // node numbers stay below 2^15 (i16 tables)

struct HexTables {
  std::vector<int16_t> face_l;    // [6][F]: local node of face f at position q
  std::vector<int8_t> face_map;   // [64][F]: frame codes (own | other << 3), own q -> other's q
  std::vector<int16_t> edge_l;    // [12][n]: local node of edge k at position t
  std::vector<int16_t> corner_l;  // [8]
  std::vector<int8_t> bit;        // [n^3]: entity bit of every local node
};

// Same loops as build_hex_tables<P> in mars_cellwise_topology.hpp.
template <int P>
HexTables buildTables() {
  constexpr int n = P + 1, F = (P - 1) * (P - 1);
  HexTables t;
  if constexpr (P > 1) {
    for (int i = 0; i < 64 * F; ++i) {
      const int q = i % F, code = i / F;
      int I, J, i2, j2;
      cw::face_to_canonical<P>(code & 7, 1 + q / (P - 1), 1 + q % (P - 1), I, J);
      cw::canonical_to_face<P>(code >> 3, I, J, i2, j2);
      t.face_map.push_back((int8_t)((i2 - 1) * (P - 1) + j2 - 1));
    }
    for (int i = 0; i < 6 * F; ++i)
      t.face_l.push_back((int16_t)cw::face_node<P>(i / F, 1 + (i % F) / (P - 1), 1 + (i % F) % (P - 1)));
  }
  for (int i = 0; i < 12 * n; ++i)
    t.edge_l.push_back((int16_t)cw::edge_node<P>(i / n, i % n));
  for (int c = 0; c < 8; ++c)
    t.corner_l.push_back((int16_t)cw::corner_node<P>(c));
  for (int l = 0; l < n * n * n; ++l)
    t.bit.push_back((int8_t)cw::entity_bit(cw::local_entity<P>(l / (n * n), (l / n) % n, l % n)));
  return t;
}

template <int... Ps>
HexTables tablesFor(int p, std::integer_sequence<int, Ps...>) {
  HexTables t;
  ((p == Ps + 1 ? (t = buildTables<Ps + 1>(), 0) : 0), ...);
  return t;
}

// Small arithmetic vocabulary, so the loop bodies below read like the hand kernel.
// Values tagged i32 stay 32-bit: divisions by the small constants are cheap there,
// and only element offsets need 64 bits (index).
struct Emit {
  OpBuilder &b;
  Location loc;

  Value idx(int64_t v) { return b.create<arith::ConstantIndexOp>(loc, v); }
  Value i32(int64_t v) { return b.create<arith::ConstantIntOp>(loc, v, 32); }
  // Zero of a float type (f32 or f64: the precision of the arrays a kernel works on).
  Value zero(Type t) {
    auto ft = cast<FloatType>(t);
    return b.create<arith::ConstantFloatOp>(loc, APFloat::getZero(ft.getFloatSemantics()), ft);
  }
  Value add(Value x, Value y) { return b.create<arith::AddIOp>(loc, x, y); }
  Value sub(Value x, Value y) { return b.create<arith::SubIOp>(loc, x, y); }
  Value mul(Value x, Value y) { return b.create<arith::MulIOp>(loc, x, y); }
  Value div(Value x, Value y) { return b.create<arith::DivUIOp>(loc, x, y); }   // operands >= 0
  Value rem(Value x, Value y) { return b.create<arith::RemUIOp>(loc, x, y); }
  Value band(Value x, Value y) { return b.create<arith::AndIOp>(loc, x, y); }
  Value lt(Value x, Value y) { return b.create<arith::CmpIOp>(loc, arith::CmpIPredicate::slt, x, y); }
  Value eq(Value x, Value y) { return b.create<arith::CmpIOp>(loc, arith::CmpIPredicate::eq, x, y); }
  Value select(Value c, Value x, Value y) { return b.create<arith::SelectOp>(loc, c, x, y); }
  Value addf(Value x, Value y) { return b.create<arith::AddFOp>(loc, x, y); }
  Value load(Value m, Value i) { return b.create<memref::LoadOp>(loc, m, ValueRange{i}); }
  void store(Value v, Value m, Value i) { b.create<memref::StoreOp>(loc, v, m, ValueRange{i}); }
  // Non-negative i32 <-> index.
  Value toIndex(Value v) { return b.create<arith::IndexCastUIOp>(loc, b.getIndexType(), v); }
  Value toI32(Value v) { return b.create<arith::IndexCastUIOp>(loc, b.getI32Type(), v); }
  // A narrow table entry, widened to i32.
  Value entry(Value table, Value i32index) {
    return b.create<arith::ExtUIOp>(loc, b.getI32Type(), load(table, toIndex(i32index)));
  }
};

// The operands of one mir.dss, plus the tables and sizes its lowering derives.
struct Dss {
  int P, n, n3, F, Ed;
  Value values, ghost, minus, faceNbr, faceCode, edgeOf, edgeOff, edgeEnt, vertOf, vertOff, vertEnt;
  SmallVector<Value> outs, copies;
  Value counting, edgeDir, vertDir, partials;   // with an epilogue
  Value faceL, faceMap, edgeL, cornerL, bit;   // constant globals
  Value localValues;                      // size(values): copies past it are ghosts
  Type real;                              // the copies' element type (f32 or f64)

  // Local copy i: values[i], or values[i] - minus[i] for the DSS of a residual.
  Value local(OpBuilder &b, Location loc, Value i) const {
    Emit x{b, loc};
    const Value v = x.load(values, i);
    return minus ? b.create<arith::SubFOp>(loc, v, x.load(minus, i)).getResult() : v;
  }

  // Copy i (index): a local copy or a ghost copy.
  Value value(OpBuilder &b, Location loc, Value i) const {
    Emit x{b, loc};
    auto branch = b.create<scf::IfOp>(
        loc, x.lt(i, localValues),
        [&](OpBuilder &tb, Location l) { tb.create<scf::YieldOp>(l, local(tb, l, i)); },
        [&](OpBuilder &eb, Location l) {
          Emit y{eb, l};
          eb.create<scf::YieldOp>(l, y.load(ghost, y.sub(i, localValues)));
        });
    return branch.getResult(0);
  }

  // Stage slot of local node l: one pad slot per row of n spreads a face's strided
  // nodes over the shared-memory banks (stage_slot<P> of the hand kernel).
  Value slot(Emit &x, Value l) const { return x.add(l, x.div(l, x.idx(n))); }

  // Sum over star entries [begin, end) of copy(entry), in order, four entries at a
  // time so that their loads overlap (star_sum of the hand kernel).
  Value starSum(OpBuilder &b, Location loc, Value ent, Value begin, Value end,
                function_ref<Value(Emit &, Value)> copyOf) const {
    Emit x{b, loc};
    auto loop = b.create<scf::ForOp>(
        loc, x.toIndex(begin), x.toIndex(end), x.idx(4), ValueRange{x.zero(real)},
        [&](OpBuilder &lb, Location l, Value m, ValueRange acc) {
          Emit y{lb, l};
          const Value stop = y.toIndex(end);
          Value valid[4], copy[4], v[4];
          for (int j = 0; j < 4; ++j) {
            const Value mj = y.add(m, y.idx(j));
            valid[j] = y.lt(mj, stop);
            copy[j] = lb.create<scf::IfOp>(
                            l, valid[j],
                            [&](OpBuilder &tb, Location tl) {
                              Emit z{tb, tl};
                              tb.create<scf::YieldOp>(tl, copyOf(z, z.load(ent, mj)));
                            },
                            [&](OpBuilder &eb, Location el) { eb.create<scf::YieldOp>(el, Emit{eb, el}.idx(0)); })
                          .getResult(0);
          }
          for (int j = 0; j < 4; ++j)
            v[j] = lb.create<scf::IfOp>(
                         l, valid[j],
                         [&](OpBuilder &tb, Location tl) { tb.create<scf::YieldOp>(tl, value(tb, tl, copy[j])); },
                         [&](OpBuilder &eb, Location el) { eb.create<scf::YieldOp>(el, Emit{eb, el}.zero(real)); })
                       .getResult(0);
          // A skipped entry leaves the sum untouched (adding 0.0 would turn -0.0 into +0.0).
          Value s = acc[0];
          for (int j = 0; j < 4; ++j) s = y.select(valid[j], y.addf(s, v[j]), s);
          lb.create<scf::YieldOp>(l, s);
        });
    return loop.getResult(0);
  }
};

Value constantTable(OpBuilder &mb, ModuleOp mod, Location loc, const std::string &name,
                    Type elem, ArrayRef<int64_t> data) {
  auto type = MemRefType::get({(int64_t)data.size()}, elem);
  if (!mod.lookupSymbol(name)) {
    OpBuilder gb = OpBuilder::atBlockBegin(mod.getBody());
    auto tensor = RankedTensorType::get({(int64_t)data.size()}, elem);
    SmallVector<APInt> bits;
    for (int64_t v : data)
      bits.push_back(APInt(elem.getIntOrFloatBitWidth(), v, /*isSigned=*/true));
    gb.create<memref::GlobalOp>(loc, gb.getStringAttr(name), gb.getStringAttr("private"),
                                TypeAttr::get(type), DenseElementsAttr::get(tensor, bits),
                                gb.getUnitAttr(), IntegerAttr());
  }
  return mb.create<memref::GetGlobalOp>(loc, type, name);
}

template <typename T>
SmallVector<int64_t> widen(const std::vector<T> &v) {
  return SmallVector<int64_t>(v.begin(), v.end());
}

// One phase: scf.parallel over [0, trips) marked {mir.lanes}.
void lanes(OpBuilder &ob, Location l, int64_t trips, function_ref<void(OpBuilder &, Location, Value)> body) {
  Emit y{ob, l};
  auto par = ob.create<scf::ParallelOp>(l, ValueRange{y.idx(0)}, ValueRange{y.idx(trips)}, ValueRange{y.idx(1)},
                                        [&](OpBuilder &pb, Location pl, ValueRange iv) { body(pb, pl, iv[0]); });
  par->setAttr("mir.lanes", ob.getUnitAttr());
}

// The same with reductions: body returns one term per init value, combined with
// `combine`; returns the reduced values.
SmallVector<Value> reduceLanes(OpBuilder &ob, Location l, int64_t trips, ValueRange inits,
                               function_ref<Value(OpBuilder &, Location, Value, Value)> combine,
                               function_ref<SmallVector<Value>(OpBuilder &, Location, Value)> body) {
  Emit y{ob, l};
  auto par = ob.create<scf::ParallelOp>(
      l, ValueRange{y.idx(0)}, ValueRange{y.idx(trips)}, ValueRange{y.idx(1)}, inits,
      [&](OpBuilder &pb, Location pl, ValueRange iv, ValueRange) {
        const SmallVector<Value> terms = body(pb, pl, iv[0]);
        auto red = pb.create<scf::ReduceOp>(pl, terms);
        for (unsigned i = 0; i < terms.size(); ++i) {
          Block &blk = red.getReductions()[i].front();
          OpBuilder rb = OpBuilder::atBlockEnd(&blk);
          rb.create<scf::ReduceReturnOp>(pl, combine(rb, pl, blk.getArgument(0), blk.getArgument(1)));
        }
      });
  par->setAttr("mir.lanes", ob.getUnitAttr());
  return SmallVector<Value>(par.getResults());
}

Value addf(OpBuilder &rb, Location rl, Value a, Value c) { return rb.create<arith::AddFOp>(rl, a, c); }

// Marks a per-element loop whose phases have `copies` iterations: an element with
// fewer than 32 copies gets that many lanes (rounded down to a power of two), so
// that --mir-warp-per-element puts several elements in a warp instead of leaving
// most lanes idle (p = 1: 8 copies, 8 lanes).
void markElements(OpBuilder &b, scf::ParallelOp loop, int64_t copies) {
  loop->setAttr("mir.elements", b.getUnitAttr());
  if (copies >= 32)
    return;
  int64_t lanes = 1;
  while (lanes * 2 <= copies) lanes *= 2;
  loop->setAttr("mir.lanes_per_element", b.getI64IntegerAttr(lanes));
}

// Bit `bit` of the 32-bit mask `word`, as i1.
Value maskBit(Emit &z, Value word, Value bit) {
  return z.b.create<arith::CmpIOp>(z.loc, arith::CmpIPredicate::ne,
                                   z.band(z.b.create<arith::ShRUIOp>(z.loc, word, bit), z.i32(1)), z.i32(0));
}

// Inlines a per-copy region (mir.dss epilogue, mir.cellwise body) with `args`, stores
// its first outs.size() results to outs at t and returns the rest (reduction terms).
SmallVector<Value> inlineBody(Emit &z, Block &body, ArrayRef<Value> args, ArrayRef<Value> outs, Value t) {
  IRMapping map;
  for (unsigned j = 0; j < args.size(); ++j) map.map(body.getArgument(j), args[j]);
  for (Operation &bop : body.without_terminator()) z.b.clone(bop, map);
  auto yield = cast<mir::YieldOp>(body.getTerminator());
  for (unsigned j = 0; j < outs.size(); ++j) z.store(map.lookup(yield.getOperand(j)), outs[j], t);
  SmallVector<Value> terms;
  for (unsigned j = outs.size(); j < yield.getNumOperands(); ++j) terms.push_back(map.lookup(yield.getOperand(j)));
  return terms;
}

// The tables that say which entities of an element are Dirichlet entities.
struct DirichletTables {
  Value faceCode, edgeOf, vertOf, edgeDir, vertDir;
};

// Element e's Dirichlet mask: bit b set when entity b (faces 0-5, edges 8-19, corners
// 20-27) is a Dirichlet entity. A phase over the 32 bits reducing with OR; on a warp
// this is the hand kernel's ballot.
Value dirichletMask(OpBuilder &eb, Location l, Value e, const DirichletTables &dt) {
  Emit y{eb, l};
  return reduceLanes(
      eb, l, 32, ValueRange{y.i32(0)},
      [](OpBuilder &rb, Location rl, Value a, Value c) -> Value { return rb.create<arith::OrIOp>(rl, a, c); },
      [&](OpBuilder &pb, Location pl, Value slotIv) -> SmallVector<Value> {
        Emit z{pb, pl};
        const Value bit = z.toI32(slotIv), mask = z.i32(0x7fffffff), zero = z.i32(0);
        auto nonzero01 = [&](Emit &t, Value v8) {
          const Value v = t.b.create<arith::ExtUIOp>(t.loc, t.b.getI32Type(), v8);
          return t.select(t.b.create<arith::CmpIOp>(t.loc, arith::CmpIPredicate::ne, v, t.i32(0)), t.i32(1),
                          t.i32(0));
        };
        auto corner = [&](OpBuilder &cb, Location cl) {
          Emit t{cb, cl};
          const Value c = t.sub(bit, t.i32(cw::kVertexBit));
          auto inRange = cb.create<arith::AndIOp>(
              cl, cb.create<arith::CmpIOp>(cl, arith::CmpIPredicate::sge, c, zero), t.lt(c, t.i32(8)));
          auto r = cb.create<scf::IfOp>(
              cl, inRange,
              [&](OpBuilder &tb, Location tl) {
                Emit u{tb, tl};
                const Value g = u.toIndex(u.load(dt.vertOf, u.add(u.mul(e, u.idx(8)), u.toIndex(c))));
                tb.create<scf::YieldOp>(tl, nonzero01(u, u.load(dt.vertDir, g)));
              },
              [&](OpBuilder &fb, Location fl) { fb.create<scf::YieldOp>(fl, Emit{fb, fl}.i32(0)); });
          cb.create<scf::YieldOp>(cl, r.getResult(0));
        };
        auto edge = [&](OpBuilder &cb, Location cl) {
          Emit t{cb, cl};
          const Value k = t.sub(bit, t.i32(cw::kEdgeBit));
          auto inRange = cb.create<arith::AndIOp>(
              cl, cb.create<arith::CmpIOp>(cl, arith::CmpIPredicate::sge, k, zero), t.lt(k, t.i32(12)));
          auto r = cb.create<scf::IfOp>(
              cl, inRange,
              [&](OpBuilder &tb, Location tl) {
                Emit u{tb, tl};
                const Value own = u.load(dt.edgeOf, u.add(u.mul(e, u.idx(12)), u.toIndex(k)));
                tb.create<scf::YieldOp>(tl, nonzero01(u, u.load(dt.edgeDir, u.toIndex(u.band(own, mask)))));
              },
              corner);
          cb.create<scf::YieldOp>(cl, r.getResult(0));
        };
        auto isDirichlet = pb.create<scf::IfOp>(
            pl, z.lt(bit, z.i32(6)),
            [&](OpBuilder &tb, Location tl) {
              Emit u{tb, tl};
              const Value code = tb.create<arith::ExtUIOp>(
                  tl, tb.getI32Type(), u.load(dt.faceCode, u.add(u.mul(e, u.idx(6)), u.toIndex(bit))));
              const Value high = tb.create<arith::ShRUIOp>(tl, code, u.i32(6));
              tb.create<scf::YieldOp>(tl, u.select(tb.create<arith::CmpIOp>(tl, arith::CmpIPredicate::ne,
                                                                            high, u.i32(0)),
                                                   u.i32(1), u.i32(0)));
            },
            edge);
        return {pb.create<arith::ShLIOp>(pl, isDirichlet.getResult(0), bit).getResult()};
      })[0];
}

LogicalResult lowerOne(mir::DssOp op) {
  const int P = (int)op.getDegree();
  if (P < 1 || P > kMaxDegree)
    return op.emitOpError("degree must be in [1, ") << kMaxDegree << "]";
  const HexTables tables = tablesFor(P, std::make_integer_sequence<int, kMaxDegree>{});

  OpBuilder b(op);
  Location loc = op.getLoc();
  ModuleOp mod = op->getParentOfType<ModuleOp>();
  Emit x{b, loc};

  Dss d;
  d.P = P;
  d.n = P + 1;
  d.n3 = d.n * d.n * d.n;
  d.F = (P - 1) * (P - 1);
  d.Ed = P - 1;
  d.values = op.getValues();
  d.ghost = op.getGhost();
  d.real = cast<MemRefType>(d.values.getType()).getElementType();
  d.minus = op.getMinus();
  d.faceNbr = op.getFaceNbr();
  d.faceCode = op.getFaceCode();
  d.edgeOf = op.getEdgeOf();
  d.edgeOff = op.getEdgeOff();
  d.edgeEnt = op.getEdgeEnt();
  d.vertOf = op.getVertOf();
  d.vertOff = op.getVertOff();
  d.vertEnt = op.getVertEnt();
  d.outs.assign(op.getOuts().begin(), op.getOuts().end());
  d.copies.assign(op.getCopies().begin(), op.getCopies().end());
  Region *epilogue = op.getEpilogue().empty() ? nullptr : &op.getEpilogue();
  const int64_t reductions = op.getReductions();
  if (epilogue) {
    d.counting = op.getFlags()[0];
    d.edgeDir = op.getFlags()[1];
    d.vertDir = op.getFlags()[2];
    d.partials = op.getPartials();
  }
  const std::string tag = "__mir_dss_p" + std::to_string(P);
  if (d.F > 0) {
    d.faceL = constantTable(b, mod, loc, tag + "_face_l", b.getI16Type(), widen(tables.face_l));
    d.faceMap = constantTable(b, mod, loc, tag + "_face_map", b.getI8Type(), widen(tables.face_map));
  }
  d.edgeL = constantTable(b, mod, loc, tag + "_edge_l", b.getI16Type(), widen(tables.edge_l));
  d.cornerL = constantTable(b, mod, loc, tag + "_corner_l", b.getI16Type(), widen(tables.corner_l));
  if (epilogue)
    d.bit = constantTable(b, mod, loc, tag + "_bit", b.getI8Type(), widen(tables.bit));
  d.localValues = b.create<memref::DimOp>(loc, d.values, 0);
  Value count = b.create<memref::DimOp>(loc, op.getElements(), 0);
  Value elements = op.getElements();

  auto elementLoop = b.create<scf::ParallelOp>(
      loc, ValueRange{x.idx(0)}, ValueRange{count}, ValueRange{x.idx(1)},
      [&](OpBuilder &eb, Location l, ValueRange w) {
        Emit y{eb, l};
        const Value e32 = y.load(elements, w[0]);
        const Value e = y.toIndex(e32);
        const Value base = y.mul(e, y.idx(d.n3));
        auto stage = eb.create<memref::AllocaOp>(
            l, MemRefType::get({(int64_t)d.n * d.n * (d.n + 1)}, d.real));
        stage->setAttr("mir.scratch", eb.getUnitAttr());

        // The element's Dirichlet mask: bit b set when entity b is a Dirichlet entity.
        Value dmask;
        if (epilogue)
          dmask = dirichletMask(eb, l, e, {d.faceCode, d.edgeOf, d.vertOf, d.edgeDir, d.vertDir});

        lanes(eb, l, d.n3, [&](OpBuilder &pb, Location pl, Value node) {
          Emit z{pb, pl};
          z.store(d.local(pb, pl, z.add(base, node)), stage, d.slot(z, node));
        });

        if (d.F > 0)
          lanes(eb, l, 6LL * d.F, [&](OpBuilder &pb, Location pl, Value pos) {
            Emit z{pb, pl};
            const Value p32 = z.toI32(pos), F = z.i32(d.F);
            const Value f = z.div(p32, F);
            const Value at = z.add(z.mul(e, z.idx(6)), z.toIndex(f));
            const Value nbr = z.load(d.faceNbr, at);
            const Value code = pb.create<arith::ExtUIOp>(pl, pb.getI32Type(), z.load(d.faceCode, at));
            pb.create<scf::IfOp>(pl, pb.create<arith::CmpIOp>(pl, arith::CmpIPredicate::sge, nbr, z.i32(0)),
                                 [&](OpBuilder &tb, Location tl) {
                                   Emit t{tb, tl};
                                   const Value q = t.entry(
                                       d.faceMap, t.add(t.mul(t.band(code, t.i32(63)), F), t.sub(p32, t.mul(f, F))));
                                   const Value other = t.entry(
                                       d.faceL, t.add(t.mul(t.rem(nbr, t.i32(6)), F), q));
                                   const Value copy = t.add(t.mul(t.toIndex(t.div(nbr, t.i32(6))), t.idx(d.n3)),
                                                            t.toIndex(other));
                                   const Value s = d.slot(t, t.toIndex(t.entry(d.faceL, p32)));
                                   t.store(t.addf(t.load(stage, s), d.value(tb, tl, copy)), stage, s);
                                   tb.create<scf::YieldOp>(tl);
                                 });
          });

        const int edgeSlots = 12 * d.Ed;
        lanes(eb, l, edgeSlots + 8, [&](OpBuilder &pb, Location pl, Value pos) {
          Emit z{pb, pl};
          const Value p32 = z.toI32(pos);
          auto corner = [&](OpBuilder &cb, Location cl) {
            Emit t{cb, cl};
            const Value c = t.sub(p32, t.i32(edgeSlots));
            const Value g = t.toIndex(t.load(d.vertOf, t.add(t.mul(e, t.idx(8)), t.toIndex(c))));
            const Value s = d.starSum(cb, cl, d.vertEnt, t.load(d.vertOff, g), t.load(d.vertOff, t.add(g, t.idx(1))),
                                      [&](Emit &u, Value entry) {
                                        return u.add(u.mul(u.toIndex(u.div(entry, u.i32(8))), u.idx(d.n3)),
                                                     u.toIndex(u.entry(d.cornerL, u.rem(entry, u.i32(8)))));
                                      });
            t.store(s, stage, d.slot(t, t.toIndex(t.entry(d.cornerL, c))));
          };
          if (edgeSlots == 0) {
            corner(pb, pl);
            return;
          }
          pb.create<scf::IfOp>(
              pl, z.lt(p32, z.i32(edgeSlots)),
              [&](OpBuilder &tb, Location tl) {
                Emit t{tb, tl};
                const Value k = t.div(p32, t.i32(d.Ed));
                const Value pos1 = t.sub(t.add(p32, t.i32(1)), t.mul(k, t.i32(d.Ed)));   // position t on the edge
                const Value own = t.load(d.edgeOf, t.add(t.mul(e, t.idx(12)), t.toIndex(k)));
                const Value mask = t.i32(0x7fffffff), zero = t.i32(0), Pc = t.i32(P);
                const Value g = t.toIndex(t.band(own, mask));
                // Canonical position on the edge: the star's entries are in the edge's
                // canonical direction, each with its own reversal bit.
                const Value tc = t.select(t.lt(own, zero), t.sub(Pc, pos1), pos1);
                const Value s = d.starSum(
                    tb, tl, d.edgeEnt, t.load(d.edgeOff, g), t.load(d.edgeOff, t.add(g, t.idx(1))),
                    [&](Emit &u, Value entry) {
                      const Value m = u.band(entry, mask);
                      const Value at = u.select(u.lt(entry, zero), u.sub(Pc, tc), tc);
                      return u.add(u.mul(u.toIndex(u.div(m, u.i32(12))), u.idx(d.n3)),
                                   u.toIndex(u.entry(d.edgeL, u.add(u.mul(u.rem(m, u.i32(12)), u.i32(d.n)), at))));
                    });
                t.store(s, stage, d.slot(t, t.toIndex(t.entry(d.edgeL, t.add(t.mul(k, t.i32(d.n)), pos1)))));
                tb.create<scf::YieldOp>(tl);
              },
              [&](OpBuilder &cb, Location cl) {
                corner(cb, cl);
                cb.create<scf::YieldOp>(cl);
              });
        });

        if (!epilogue) {
          lanes(eb, l, d.n3, [&](OpBuilder &pb, Location pl, Value node) {
            Emit z{pb, pl};
            z.store(z.load(stage, d.slot(z, node)), d.outs[0], z.add(base, node));
          });
          return;
        }
        const Value countingBits = y.load(d.counting, e);
        SmallVector<Value> zeros(reductions, reductions ? y.zero(cast<MemRefType>(d.partials.getType()).getElementType())
                                                        : Value());
        const SmallVector<Value> sums = reduceLanes(
            eb, l, d.n3, zeros,
            addf, [&](OpBuilder &pb, Location pl, Value node) -> SmallVector<Value> {
              Emit z{pb, pl};
              const Value t = z.add(base, node);
              const Value bit = z.entry(d.bit, z.toI32(node));
              SmallVector<Value> args{z.load(stage, d.slot(z, node)), maskBit(z, dmask, bit),
                                      maskBit(z, countingBits, bit)};
              for (Value c : d.copies) args.push_back(z.load(c, t));
              return inlineBody(z, epilogue->front(), args, d.outs, t);
            });
        for (int64_t j = 0; j < reductions; ++j)
          y.store(sums[j], d.partials, y.add(y.mul(w[0], y.idx(reductions)), y.idx(j)));
      });
  markElements(b, elementLoop, d.n3);
  op.erase();
  return success();
}

// mir.cellwise: one phase per element over its n^3 copies, the body inlined per
// copy, its reduction terms summed over the element (scf.reduce).
LogicalResult lowerOne(mir::CellwiseOp op) {
  const int P = (int)op.getDegree();
  if (P < 1 || P > kMaxDegree)
    return op.emitOpError("degree must be in [1, ") << kMaxDegree << "]";
  const int n = P + 1, n3 = n * n * n;
  OpBuilder b(op);
  Location loc = op.getLoc();
  ModuleOp mod = op->getParentOfType<ModuleOp>();
  Emit x{b, loc};
  Value bitTable;
  const bool dirichlet = !op.getDirichlet().empty();
  if (op.getCounting() || dirichlet) {
    const HexTables tables = tablesFor(P, std::make_integer_sequence<int, kMaxDegree>{});
    bitTable = constantTable(b, mod, loc, "__mir_dss_p" + std::to_string(P) + "_bit", b.getI8Type(),
                             widen(tables.bit));
  }
  const int64_t reductions = op.getReductions();
  const SmallVector<Value> inputs(op.getInputs()), outs(op.getOuts());
  Value elements = op.getElements(), counting = op.getCounting(), partials = op.getPartials();
  Value count = b.create<memref::DimOp>(loc, elements, 0);
  Block &body = op.getBody().front();

  auto elementLoop = b.create<scf::ParallelOp>(
      loc, ValueRange{x.idx(0)}, ValueRange{count}, ValueRange{x.idx(1)},
      [&](OpBuilder &eb, Location l, ValueRange w) {
        Emit y{eb, l};
        const Value e = y.toIndex(y.load(elements, w[0]));
        const Value base = y.mul(e, y.idx(n3));
        const Value bits = counting ? y.load(counting, e) : Value();
        const Value dmask =
            dirichlet ? dirichletMask(eb, l, e, {op.getDirichlet()[0], op.getDirichlet()[1], op.getDirichlet()[2],
                                                 op.getDirichlet()[3], op.getDirichlet()[4]})
                      : Value();
        auto copy = [&](OpBuilder &pb, Location pl, Value node) {
          Emit z{pb, pl};
          const Value t = z.add(base, node);
          const Value bit = bitTable ? z.entry(bitTable, z.toI32(node)) : Value();
          SmallVector<Value> args{counting ? maskBit(z, bits, bit) : Value(pb.create<arith::ConstantIntOp>(pl, 1, 1))};
          if (dirichlet) args.push_back(maskBit(z, dmask, bit));
          for (Value in : inputs) args.push_back(z.load(in, t));
          return inlineBody(z, body, args, outs, t);
        };
        if (reductions == 0) {
          lanes(eb, l, n3, [&](OpBuilder &pb, Location pl, Value node) { copy(pb, pl, node); });
          return;
        }
        SmallVector<Value> zeros(reductions, reductions ? y.zero(cast<MemRefType>(partials.getType()).getElementType())
                                                        : Value());
        const SmallVector<Value> sums = reduceLanes(eb, l, n3, zeros, addf, copy);
        for (int64_t j = 0; j < reductions; ++j)
          y.store(sums[j], partials, y.add(y.mul(w[0], y.idx(reductions)), y.idx(j)));
      });
  markElements(b, elementLoop, n3);
  op.erase();
  return success();
}

// mir.interpolate: out (+)= M (x) M (x) M in per element, by sum factorization:
// three phases, one axis each (the last first), through two scratch buffers.
LogicalResult lowerOne(mir::InterpolateOp op) {
  const int64_t R = op.getRows(), C = op.getCols();
  OpBuilder b(op);
  Location loc = op.getLoc();
  Emit x{b, loc};
  Value elements = op.getElements(), counting = op.getCounting(), M = op.getMatrix(), in = op.getInput(),
        out = op.getOutput();
  Value bitTable;
  if (counting) {
    const int P = (int)R - 1;
    if (P < 1 || P > kMaxDegree)
      return op.emitOpError("counting needs rows - 1 in [1, ") << kMaxDegree << "]";
    const HexTables tables = tablesFor(P, std::make_integer_sequence<int, kMaxDegree>{});
    bitTable = constantTable(b, op->getParentOfType<ModuleOp>(), loc, "__mir_dss_p" + std::to_string(P) + "_bit",
                             b.getI8Type(), widen(tables.bit));
  }
  const bool accumulate = op.getAccumulate();
  const Type real = cast<MemRefType>(in.getType()).getElementType();   // f32 or f64
  Value count = b.create<memref::DimOp>(loc, elements, 0);
  // sum_{q < C} M[row][q] * value(q), in order of q.
  auto contract = [&](Emit &z, Value row, function_ref<Value(Emit &, Value)> value) {
    auto loop = z.b.create<scf::ForOp>(
        z.loc, z.idx(0), z.idx(C), z.idx(1), ValueRange{z.zero(real)},
        [&](OpBuilder &lb, Location ll, Value q, ValueRange acc) {
          Emit u{lb, ll};
          const Value m = u.load(M, u.add(u.mul(row, u.idx(C)), q));
          lb.create<scf::YieldOp>(ll, u.addf(acc[0], lb.create<arith::MulFOp>(ll, m, value(u, q)).getResult()));
        });
    return loop.getResult(0);
  };
  auto elementLoop = b.create<scf::ParallelOp>(
      loc, ValueRange{x.idx(0)}, ValueRange{count}, ValueRange{x.idx(1)},
      [&](OpBuilder &eb, Location l, ValueRange w) {
        Emit y{eb, l};
        const Value e = y.toIndex(y.load(elements, w[0]));
        const Value baseIn = y.mul(e, y.idx(C * C * C)), baseOut = y.mul(e, y.idx(R * R * R));
        const Value bits = counting ? y.load(counting, e) : Value();
        auto scratch = [&](int64_t size) {
          auto a = eb.create<memref::AllocaOp>(l, MemRefType::get({size}, real));
          a->setAttr("mir.scratch", eb.getUnitAttr());
          return a.getResult();
        };
        const Value s1 = scratch(C * C * R), s2 = scratch(C * R * R);
        lanes(eb, l, C * C * R, [&](OpBuilder &pb, Location pl, Value q) {   // s1[j][k][c]
          Emit z{pb, pl};
          const Value jk = z.div(q, z.idx(R)), c = z.rem(q, z.idx(R));
          z.store(contract(z, c, [&](Emit &u, Value m) { return u.load(in, u.add(baseIn, u.add(u.mul(jk, u.idx(C)), m))); }),
                  s1, q);
        });
        lanes(eb, l, C * R * R, [&](OpBuilder &pb, Location pl, Value q) {   // s2[j][b][c]
          Emit z{pb, pl};
          const Value j = z.div(q, z.idx(R * R)), bb = z.rem(z.div(q, z.idx(R)), z.idx(R)), c = z.rem(q, z.idx(R));
          z.store(contract(z, bb,
                           [&](Emit &u, Value k) { return u.load(s1, u.add(u.mul(u.add(u.mul(j, u.idx(C)), k), u.idx(R)), c)); }),
                  s2, q);
        });
        lanes(eb, l, R * R * R, [&](OpBuilder &pb, Location pl, Value q) {   // out[a][b][c]
          Emit z{pb, pl};
          const Value a = z.div(q, z.idx(R * R)), bc = z.rem(q, z.idx(R * R));
          Value v = contract(z, a, [&](Emit &u, Value j) { return u.load(s2, u.add(u.mul(j, u.idx(R * R)), bc)); });
          const Value t = z.add(baseOut, q);
          if (accumulate) v = z.addf(z.load(out, t), v);
          if (counting)
            v = pb.create<arith::SelectOp>(pl, maskBit(z, bits, z.entry(bitTable, z.toI32(q))), v, z.zero(real));
          z.store(v, out, t);
        });
      });
  elementLoop->setAttr("mir.elements", b.getUnitAttr());
  op.erase();
  return success();
}

// ---- block-structured DSS: mir.block_sum, mir.block_skeleton ---------------------------
// The bodies of block_node_sum and block_skeleton_node in MARS's
// mars_cellwise_block_dss.hpp, in the same summation order, so the same bits. Both are
// flat: one work item per block node or skeleton node, 32 items to an "element" with
// one lane each, 8 elements to a thread block, so --mir-warp-per-element maps item t to
// consecutive threads as the hand kernels' grid-stride loops do.

struct BlockGeom {
  int64_t P, B, n, n3, B3, Pm;
};

BlockGeom blockGeom(int64_t P, int64_t B) {
  const int64_t n = P + 1;
  return {P, B, n, n * n * n, B * B * B, B * P};
}

// The copies of block node coordinate X (index) along one axis (axis_copies): how many
// (1 or 2, index), and element index l and local coordinate a (i32) of the first (the
// representative) and the second. Everything inside a block is 32-bit, as in the hand
// kernel: 64-bit division by a constant costs several times the instructions.
struct AxisCopies {
  Value count, l0, a0, l1, a1;
};

AxisCopies axisCopies(Emit &z, const BlockGeom &g, Value Xindex) {
  const Value X = z.toI32(Xindex);
  const Value last = z.eq(X, z.i32(g.Pm)), q = z.div(X, z.i32(g.P)), r = z.rem(X, z.i32(g.P));
  AxisCopies c;
  c.l0 = z.select(last, z.i32(g.B - 1), q);
  c.a0 = z.select(last, z.i32(g.P), r);
  const Value two = z.band(z.band(z.eq(r, z.i32(0)), z.lt(z.i32(0), X)), z.lt(X, z.i32(g.Pm)));
  c.count = z.select(two, z.idx(2), z.idx(1));
  c.l1 = z.sub(c.l0, z.i32(1));
  c.a1 = z.i32(g.P);
  return c;
}

// Copy index of element (l0, l1, l2) of block blk (index), node (a0, a1, a2) (i32).
Value blockCopy(Emit &z, const BlockGeom &g, Value blk, const Value (&l)[3], const Value (&a)[3]) {
  const Value le = z.add(z.mul(z.add(z.mul(l[0], z.i32(g.B)), l[1]), z.i32(g.B)), l[2]);
  const Value node = z.add(z.mul(z.add(z.mul(a[0], z.i32(g.n)), a[1]), z.i32(g.n)), a[2]);
  const Value e = z.add(z.mul(blk, z.idx(g.B3)), z.toIndex(le));
  return z.add(z.mul(e, z.idx(g.n3)), z.toIndex(node));
}

Value blockRepresentative(Emit &z, const BlockGeom &g, Value blk, const Value (&x)[3]) {
  const AxisCopies c[3] = {axisCopies(z, g, x[0]), axisCopies(z, g, x[1]), axisCopies(z, g, x[2])};
  return blockCopy(z, g, blk, {c[0].l0, c[1].l0, c[2].l0}, {c[0].a0, c[1].a0, c[2].a0});
}

// The copies of block node x of block blk in lexicographic order (for_block_copies):
// with `acc`, folds body(copy, acc) and returns the result; without, runs body(copy).
Value forBlockCopies(Emit &z, const BlockGeom &g, Value blk, const Value (&x)[3], Value acc,
                     function_ref<Value(Emit &, Value, Value)> body) {
  const AxisCopies c[3] = {axisCopies(z, g, x[0]), axisCopies(z, g, x[1]), axisCopies(z, g, x[2])};
  const SmallVector<Value> inits = acc ? SmallVector<Value>{acc} : SmallVector<Value>{};
  std::function<Value(Emit &, int, Value (&)[3], Value (&)[3], ValueRange)> level =
      [&](Emit &y, int d, Value (&l)[3], Value (&a)[3], ValueRange carried) -> Value {
    auto loop = y.b.create<scf::ForOp>(
        y.loc, y.idx(0), c[d].count, y.idx(1), carried,
        [&](OpBuilder &lb, Location ll, Value i, ValueRange in) {
          Emit u{lb, ll};
          const Value first = u.eq(i, u.idx(0));
          l[d] = u.select(first, c[d].l0, c[d].l1);
          a[d] = u.select(first, c[d].a0, c[d].a1);
          Value r = d < 2 ? level(u, d + 1, l, a, in) : body(u, blockCopy(u, g, blk, l, a), in.empty() ? Value() : in[0]);
          if (in.empty())
            lb.create<scf::YieldOp>(ll);
          else
            lb.create<scf::YieldOp>(ll, r);
        });
    return loop.getNumResults() ? loop.getResult(0) : Value();
  };
  Value l[3], a[3];
  return level(z, 0, l, a, inits);
}

// Flat loop over `total` items: elements of 32 items, one lane each, item(t) per item.
void flatItems(OpBuilder &b, Location loc, Value total, function_ref<void(Emit &, Value)> item) {
  Emit x{b, loc};
  const Value chunks = x.div(x.add(total, x.idx(31)), x.idx(32));
  auto loop = b.create<scf::ParallelOp>(
      loc, ValueRange{x.idx(0)}, ValueRange{chunks}, ValueRange{x.idx(1)},
      [&](OpBuilder &eb, Location l, ValueRange w) {
        lanes(eb, l, 32, [&](OpBuilder &pb, Location pl, Value lane) {
          Emit z{pb, pl};
          const Value t = z.add(z.mul(w[0], z.idx(32)), lane);
          pb.create<scf::IfOp>(pl, z.lt(t, total), [&](OpBuilder &tb, Location tl) {
            Emit u{tb, tl};
            item(u, t);
            tb.create<scf::YieldOp>(tl);
          });
        });
      });
  loop->setAttr("mir.elements", b.getUnitAttr());
  loop->setAttr("mir.warps_per_block", b.getI64IntegerAttr(8));
}

LogicalResult lowerOne(mir::BlockSumOp op) {
  const BlockGeom g = blockGeom(op.getDegree(), op.getBlock());
  OpBuilder b(op);
  Location loc = op.getLoc();
  Emit x{b, loc};
  const Value values = op.getValues(), out = op.getOut();
  const int64_t m = g.Pm + 1;
  const Value blocks = x.div(b.create<memref::DimOp>(loc, values, 0), x.idx(g.B3 * g.n3));
  const Value total = x.mul(blocks, x.idx(m * m * m));
  const Type real = cast<MemRefType>(values.getType()).getElementType();
  flatItems(b, loc, total, [&](Emit &z, Value t) {
    // One 64-bit division for the block; the node inside it in 32 bits.
    const Value blk = z.div(t, z.idx(m * m * m));
    const Value r = z.toI32(z.sub(t, z.mul(blk, z.idx(m * m * m))));
    const Value xyz[3] = {z.toIndex(z.div(r, z.i32(m * m))), z.toIndex(z.rem(z.div(r, z.i32(m)), z.i32(m))),
                          z.toIndex(z.rem(r, z.i32(m)))};
    const Value s = forBlockCopies(z, g, blk, xyz, z.zero(real),
                                   [&](Emit &u, Value k, Value acc) { return u.addf(acc, u.load(values, k)); });
    forBlockCopies(z, g, blk, xyz, Value(), [&](Emit &u, Value k, Value) {
      u.store(s, out, k);
      return Value();
    });
  });
  op.erase();
  return success();
}

LogicalResult lowerOne(mir::BlockSkeletonOp op) {
  const BlockGeom g = blockGeom(op.getDegree(), op.getBlock());
  OpBuilder b(op);
  Location loc = op.getLoc();
  Emit x{b, loc};
  const Value v = op.getValues(), faceNbr = op.getFaceNbr(), faceCode = op.getFaceCode(), edgeOff = op.getEdgeOff(),
              edgeEnt = op.getEdgeEnt(), vertOff = op.getVertOff(), vertEnt = op.getVertEnt();
  const int64_t Pm = g.Pm, m = Pm - 1;
  const Type real = cast<MemRefType>(v.getType()).getElementType();
  const Value blocks = x.div(b.create<memref::DimOp>(loc, faceNbr, 0), x.idx(6));
  const Value edges = x.sub(b.create<memref::DimOp>(loc, edgeOff, 0), x.idx(1));
  const Value vertices = x.sub(b.create<memref::DimOp>(loc, vertOff, 0), x.idx(1));
  const Value faces = x.mul(blocks, x.idx(6 * m * m)), edgeItems = x.mul(edges, x.idx(m));
  const Value total = x.add(x.add(faces, edgeItems), vertices);
  // Coordinates on face f at (i, j) and on edge k at position t (block_face_xyz,
  // block_edge_xyz): o0 / o1, the other axes in increasing order.
  auto faceXyz = [&](Emit &z, Value f, Value i, Value j, Value (&xyz)[3]) {
    const Value axis = z.div(f, z.idx(2)), side = z.mul(z.rem(f, z.idx(2)), z.idx(Pm));
    const Value a0 = z.eq(axis, z.idx(0)), a1 = z.eq(axis, z.idx(1)), a2 = z.eq(axis, z.idx(2));
    xyz[0] = z.select(a0, side, i);
    xyz[1] = z.select(a0, i, z.select(a1, side, j));
    xyz[2] = z.select(a2, side, j);
  };
  auto edgeXyz = [&](Emit &z, Value k, Value t, Value (&xyz)[3]) {
    const Value axis = z.div(k, z.idx(4));
    const Value side0 = z.mul(z.rem(z.div(k, z.idx(2)), z.idx(2)), z.idx(Pm)), side1 = z.mul(z.rem(k, z.idx(2)), z.idx(Pm));
    const Value a0 = z.eq(axis, z.idx(0)), a1 = z.eq(axis, z.idx(1)), a2 = z.eq(axis, z.idx(2));
    xyz[0] = z.select(a0, t, side0);
    xyz[1] = z.select(a0, side0, z.select(a1, t, side1));
    xyz[2] = z.select(a2, t, side1);
  };
  auto cornerXyz = [&](Emit &z, Value c, Value (&xyz)[3]) {
    const Value q = z.rem(c, z.idx(4));
    const Value bx = z.b.create<arith::OrIOp>(z.loc, z.eq(q, z.idx(1)), z.eq(q, z.idx(2)));
    const Value by = z.lt(z.idx(1), q);
    xyz[0] = z.select(bx, z.idx(Pm), z.idx(0));
    xyz[1] = z.select(by, z.idx(Pm), z.idx(0));
    xyz[2] = z.mul(z.div(c, z.idx(4)), z.idx(Pm));
  };
  auto scatter = [&](Emit &z, Value blk, const Value (&xyz)[3], Value s) {
    forBlockCopies(z, g, blk, xyz, Value(), [&](Emit &u, Value k, Value) {
      u.store(s, v, k);
      return Value();
    });
  };
  // A star [begin, end) of 2 or more entries: sum one copy per block in order, then
  // write the sum to every copy. at(entry) gives the entry's block and node.
  auto star = [&](Emit &z, Value off, Value ent, Value g0,
                  function_ref<Value(Emit &, Value, Value (&)[3])> at) {
    const Value begin = z.toIndex(z.load(off, g0)), end = z.toIndex(z.load(off, z.add(g0, z.idx(1))));
    z.b.create<scf::IfOp>(z.loc, z.lt(z.idx(1), z.sub(end, begin)), [&](OpBuilder &tb, Location tl) {
      Emit u{tb, tl};
      auto sum = tb.create<scf::ForOp>(tl, begin, end, u.idx(1), ValueRange{u.zero(real)},
                                       [&](OpBuilder &lb, Location ll, Value e, ValueRange acc) {
                                         Emit w{lb, ll};
                                         Value xyz[3];
                                         const Value blk = at(w, w.load(ent, e), xyz);
                                         lb.create<scf::YieldOp>(
                                             ll, w.addf(acc[0], w.load(v, blockRepresentative(w, g, blk, xyz))));
                                       });
      tb.create<scf::ForOp>(tl, begin, end, u.idx(1), ValueRange{},
                            [&](OpBuilder &lb, Location ll, Value e, ValueRange) {
                              Emit w{lb, ll};
                              Value xyz[3];
                              const Value blk = at(w, w.load(ent, e), xyz);
                              scatter(w, blk, xyz, sum.getResult(0));
                              lb.create<scf::YieldOp>(ll);
                            });
      tb.create<scf::YieldOp>(tl);
    });
  };
  auto vertexItem = [&](Emit &z, Value t) {
    star(z, vertOff, vertEnt, t, [&](Emit &w, Value ent, Value (&xyz)[3]) {
      const Value c = w.toIndex(ent);
      cornerXyz(w, w.rem(c, w.idx(8)), xyz);
      return w.div(c, w.idx(8));
    });
  };
  auto edgeItem = [&](Emit &z, Value t) {
    const Value T = z.add(z.idx(1), z.rem(t, z.idx(m)));
    star(z, edgeOff, edgeEnt, z.div(t, z.idx(m)), [&](Emit &w, Value ent, Value (&xyz)[3]) {
      const Value k = w.toIndex(w.band(ent, w.i32(0x7fffffff)));
      const Value pos = w.select(w.lt(ent, w.i32(0)), w.sub(w.idx(Pm), T), T);
      edgeXyz(w, w.rem(k, w.idx(12)), pos, xyz);
      return w.div(k, w.idx(12));
    });
  };
  auto faceItem = [&](Emit &z, Value t) {
    const Value c = z.div(t, z.idx(m * m)), q = z.rem(t, z.idx(m * m));
    const Value nbr = z.load(faceNbr, c);
    z.b.create<scf::IfOp>(z.loc, z.lt(z.toI32(c), nbr), [&](OpBuilder &tb, Location tl) {   // the lower copy only
      Emit u{tb, tl};
      const Value code = u.toIndex(u.entry(faceCode, u.toI32(c)));
      const Value own = u.rem(code, u.idx(8)), oth = u.rem(u.div(code, u.idx(8)), u.idx(8));
      auto bit = [&](Value word, int64_t k) { return u.eq(u.rem(u.div(word, u.idx(k)), u.idx(2)), u.idx(1)); };
      const Value i = u.add(u.idx(1), u.div(q, u.idx(m))), j = u.add(u.idx(1), u.rem(q, u.idx(m)));
      // face_to_canonical with the own frame, then canonical_to_face with the other's.
      const Value di = u.select(bit(own, 1), u.sub(u.idx(Pm), i), i), dj = u.select(bit(own, 2), u.sub(u.idx(Pm), j), j);
      const Value I = u.select(bit(own, 4), di, dj), J = u.select(bit(own, 4), dj, di);
      const Value ei = u.select(bit(oth, 4), I, J), ej = u.select(bit(oth, 4), J, I);
      const Value i2 = u.select(bit(oth, 1), u.sub(u.idx(Pm), ei), ei), j2 = u.select(bit(oth, 2), u.sub(u.idx(Pm), ej), ej);
      const Value nb = u.toIndex(nbr);
      Value xa[3], xb[3];
      faceXyz(u, u.rem(c, u.idx(6)), i, j, xa);
      faceXyz(u, u.rem(nb, u.idx(6)), i2, j2, xb);
      const Value blkA = u.div(c, u.idx(6)), blkB = u.div(nb, u.idx(6));
      const Value s = u.addf(u.load(v, blockRepresentative(u, g, blkA, xa)), u.load(v, blockRepresentative(u, g, blkB, xb)));
      scatter(u, blkA, xa, s);
      scatter(u, blkB, xb, s);
      tb.create<scf::YieldOp>(tl);
    });
  };
  flatItems(b, loc, total, [&](Emit &z, Value t) {
    if (m == 0) {   // no interior face or edge nodes: corners only
      vertexItem(z, t);
      return;
    }
    z.b.create<scf::IfOp>(
        z.loc, z.lt(t, faces), [&](OpBuilder &tb, Location tl) {
          Emit u{tb, tl};
          faceItem(u, t);
          tb.create<scf::YieldOp>(tl);
        },
        [&](OpBuilder &eb, Location el) {
          Emit u{eb, el};
          const Value te = u.sub(t, faces);
          eb.create<scf::IfOp>(
              el, u.lt(te, edgeItems), [&](OpBuilder &tb, Location tl) {
                Emit w{tb, tl};
                edgeItem(w, te);
                tb.create<scf::YieldOp>(tl);
              },
              [&](OpBuilder &vb, Location vl) {
                Emit w{vb, vl};
                vertexItem(w, w.sub(te, edgeItems));
                vb.create<scf::YieldOp>(vl);
              });
          eb.create<scf::YieldOp>(el);
        });
  });
  op.erase();
  return success();
}

struct LowerCellwisePass : public PassWrapper<LowerCellwisePass, OperationPass<ModuleOp>> {
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(LowerCellwisePass)

  StringRef getArgument() const final { return "mir-lower-cellwise"; }
  StringRef getDescription() const final {
    return "Lower mir.dss (the staged gather: per-element phases over own copies, face "
           "pairs, edge/corner stars, write-back), mir.cellwise (a per-copy phase) and the "
           "block-structured DSS (mir.block_sum, mir.block_skeleton: flat item loops) to "
           "per-element parallel loops";
  }
  void getDependentDialects(DialectRegistry &registry) const override {
    registry.insert<arith::ArithDialect, memref::MemRefDialect, scf::SCFDialect>();
  }

  void runOnOperation() override {
    SmallVector<Operation *> ops;
    getOperation()->walk([&](Operation *op) {
      if (isa<mir::DssOp, mir::CellwiseOp, mir::InterpolateOp, mir::BlockSumOp, mir::BlockSkeletonOp>(op))
        ops.push_back(op);
    });
    for (Operation *op : ops)
      if (failed(isa<mir::DssOp>(op)           ? lowerOne(cast<mir::DssOp>(op))
                 : isa<mir::CellwiseOp>(op)    ? lowerOne(cast<mir::CellwiseOp>(op))
                 : isa<mir::BlockSumOp>(op)    ? lowerOne(cast<mir::BlockSumOp>(op))
                 : isa<mir::BlockSkeletonOp>(op) ? lowerOne(cast<mir::BlockSkeletonOp>(op))
                                               : lowerOne(cast<mir::InterpolateOp>(op))))
        return signalPassFailure();
  }
};

}  // namespace

namespace mir {
void registerLowerCellwisePass() { PassRegistration<LowerCellwisePass>(); }
}  // namespace mir
