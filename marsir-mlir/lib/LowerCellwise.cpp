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
  Value f64(double v) { return b.create<arith::ConstantFloatOp>(loc, APFloat(v), b.getF64Type()); }
  Value add(Value x, Value y) { return b.create<arith::AddIOp>(loc, x, y); }
  Value sub(Value x, Value y) { return b.create<arith::SubIOp>(loc, x, y); }
  Value mul(Value x, Value y) { return b.create<arith::MulIOp>(loc, x, y); }
  Value div(Value x, Value y) { return b.create<arith::DivUIOp>(loc, x, y); }   // operands >= 0
  Value rem(Value x, Value y) { return b.create<arith::RemUIOp>(loc, x, y); }
  Value band(Value x, Value y) { return b.create<arith::AndIOp>(loc, x, y); }
  Value lt(Value x, Value y) { return b.create<arith::CmpIOp>(loc, arith::CmpIPredicate::slt, x, y); }
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
  Value values, ghost, faceNbr, faceCode, edgeOf, edgeOff, edgeEnt, vertOf, vertOff, vertEnt;
  SmallVector<Value> outs, copies;
  Value counting, edgeDir, vertDir, partials;   // with an epilogue
  Value faceL, faceMap, edgeL, cornerL, bit;   // constant globals
  Value localValues;                      // size(values): copies past it are ghosts

  // Copy i (index): a local copy or a ghost copy.
  Value value(OpBuilder &b, Location loc, Value i) const {
    Emit x{b, loc};
    auto branch = b.create<scf::IfOp>(
        loc, x.lt(i, localValues),
        [&](OpBuilder &tb, Location l) { tb.create<scf::YieldOp>(l, Emit{tb, l}.load(values, i)); },
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
        loc, x.toIndex(begin), x.toIndex(end), x.idx(4), ValueRange{x.f64(0.0)},
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
                         [&](OpBuilder &eb, Location el) { eb.create<scf::YieldOp>(el, Emit{eb, el}.f64(0.0)); })
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
            l, MemRefType::get({(int64_t)d.n * d.n * (d.n + 1)}, eb.getF64Type()));
        stage->setAttr("mir.scratch", eb.getUnitAttr());

        // The element's Dirichlet mask: bit b set when entity b is a Dirichlet entity.
        Value dmask;
        if (epilogue)
          dmask = reduceLanes(
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
                        const Value g = u.toIndex(u.load(d.vertOf, u.add(u.mul(e, u.idx(8)), u.toIndex(c))));
                        tb.create<scf::YieldOp>(tl, nonzero01(u, u.load(d.vertDir, g)));
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
                        const Value own = u.load(d.edgeOf, u.add(u.mul(e, u.idx(12)), u.toIndex(k)));
                        tb.create<scf::YieldOp>(tl, nonzero01(u, u.load(d.edgeDir, u.toIndex(u.band(own, mask)))));
                      },
                      corner);
                  cb.create<scf::YieldOp>(cl, r.getResult(0));
                };
                auto isDirichlet = pb.create<scf::IfOp>(
                    pl, z.lt(bit, z.i32(6)),
                    [&](OpBuilder &tb, Location tl) {
                      Emit u{tb, tl};
                      const Value code = tb.create<arith::ExtUIOp>(
                          tl, tb.getI32Type(), u.load(d.faceCode, u.add(u.mul(e, u.idx(6)), u.toIndex(bit))));
                      const Value high = tb.create<arith::ShRUIOp>(tl, code, u.i32(6));
                      tb.create<scf::YieldOp>(tl, u.select(tb.create<arith::CmpIOp>(tl, arith::CmpIPredicate::ne,
                                                                                    high, u.i32(0)),
                                                           u.i32(1), u.i32(0)));
                    },
                    edge);
                return {pb.create<arith::ShLIOp>(pl, isDirichlet.getResult(0), bit).getResult()};
              })[0];

        lanes(eb, l, d.n3, [&](OpBuilder &pb, Location pl, Value node) {
          Emit z{pb, pl};
          z.store(z.load(d.values, z.add(base, node)), stage, d.slot(z, node));
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
        SmallVector<Value> zeros(reductions, y.f64(0.0));
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
  elementLoop->setAttr("mir.elements", b.getUnitAttr());
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
  if (op.getCounting()) {
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
        auto copy = [&](OpBuilder &pb, Location pl, Value node) {
          Emit z{pb, pl};
          const Value t = z.add(base, node);
          SmallVector<Value> args{counting ? maskBit(z, bits, z.entry(bitTable, z.toI32(node)))
                                           : Value(pb.create<arith::ConstantIntOp>(pl, 1, 1))};
          for (Value in : inputs) args.push_back(z.load(in, t));
          return inlineBody(z, body, args, outs, t);
        };
        if (reductions == 0) {
          lanes(eb, l, n3, [&](OpBuilder &pb, Location pl, Value node) { copy(pb, pl, node); });
          return;
        }
        SmallVector<Value> zeros(reductions, y.f64(0.0));
        const SmallVector<Value> sums = reduceLanes(eb, l, n3, zeros, addf, copy);
        for (int64_t j = 0; j < reductions; ++j)
          y.store(sums[j], partials, y.add(y.mul(w[0], y.idx(reductions)), y.idx(j)));
      });
  elementLoop->setAttr("mir.elements", b.getUnitAttr());
  op.erase();
  return success();
}

struct LowerCellwisePass : public PassWrapper<LowerCellwisePass, OperationPass<ModuleOp>> {
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(LowerCellwisePass)

  StringRef getArgument() const final { return "mir-lower-cellwise"; }
  StringRef getDescription() const final {
    return "Lower mir.dss (the staged gather: per-element phases over own copies, face "
           "pairs, edge/corner stars, write-back) and mir.cellwise (a per-copy phase) to "
           "per-element parallel loops";
  }
  void getDependentDialects(DialectRegistry &registry) const override {
    registry.insert<arith::ArithDialect, memref::MemRefDialect, scf::SCFDialect>();
  }

  void runOnOperation() override {
    SmallVector<Operation *> ops;
    getOperation()->walk([&](Operation *op) {
      if (isa<mir::DssOp, mir::CellwiseOp>(op))
        ops.push_back(op);
    });
    for (Operation *op : ops)
      if (failed(isa<mir::DssOp>(op) ? lowerOne(cast<mir::DssOp>(op)) : lowerOne(cast<mir::CellwiseOp>(op))))
        return signalPassFailure();
  }
};

}  // namespace

namespace mir {
void registerLowerCellwisePass() { PassRegistration<LowerCellwisePass>(); }
}  // namespace mir
