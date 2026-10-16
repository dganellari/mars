// --mir-lower-dss: lower mir.dss to the staged gather of the hand-written MARS
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
// prefetch-descriptors (default on): the first phase also loads the element's
// 26 entity descriptors (face neighbours and frame codes, edge entries and star
// bounds, corner star bounds) into a second scratch buffer, so their latency
// overlaps the staging loads and the later phases find them in shared memory
// instead of starting a chain of dependent global loads. The hand kernel does
// the same with one descriptor per lane and shuffles.

#include "mir/MirOps.h"
#include "mir/MirPasses.h"

#include "backend/distributed/unstructured/solvers/mars_cellwise_hex.hpp"

#include "mlir/Dialect/Arith/IR/Arith.h"
#include "mlir/Dialect/Func/IR/FuncOps.h"
#include "mlir/Dialect/MemRef/IR/MemRef.h"
#include "mlir/Dialect/SCF/IR/SCF.h"
#include "mlir/IR/Builders.h"
#include "mlir/IR/BuiltinOps.h"
#include "mlir/Pass/Pass.h"

#include <algorithm>
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
  Value values, ghost, faceNbr, faceCode, edgeOf, edgeOff, edgeEnt, vertOf, vertOff, vertEnt, out;
  Value faceL, faceMap, edgeL, cornerL;   // constant globals
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

// Descriptor scratch layout (i32): face neighbours, face codes, edge entries,
// edge star bounds (begin, end), corner star bounds (begin, end).
constexpr int kDescFaceNbr = 0, kDescFaceCode = 6, kDescEdge = 12, kDescEdgeOff = 24, kDescVertOff = 48,
              kDescSize = 64, kEntities = 26;

LogicalResult lowerOne(mir::DssOp op, bool prefetch) {
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
  d.out = op.getOut();
  const std::string tag = "__mir_dss_p" + std::to_string(P);
  if (d.F > 0) {
    d.faceL = constantTable(b, mod, loc, tag + "_face_l", b.getI16Type(), widen(tables.face_l));
    d.faceMap = constantTable(b, mod, loc, tag + "_face_map", b.getI8Type(), widen(tables.face_map));
  }
  d.edgeL = constantTable(b, mod, loc, tag + "_edge_l", b.getI16Type(), widen(tables.edge_l));
  d.cornerL = constantTable(b, mod, loc, tag + "_corner_l", b.getI16Type(), widen(tables.corner_l));
  d.localValues = b.create<memref::DimOp>(loc, d.values, 0);
  Value count = b.create<memref::DimOp>(loc, op.getElements(), 0);
  Value elements = op.getElements();

  auto lanes = [&](OpBuilder &ob, Location l, int64_t trips,
                   function_ref<void(OpBuilder &, Location, Value)> body) {
    Emit y{ob, l};
    auto par = ob.create<scf::ParallelOp>(
        l, ValueRange{y.idx(0)}, ValueRange{y.idx(trips)}, ValueRange{y.idx(1)},
        [&](OpBuilder &pb, Location pl, ValueRange iv) { body(pb, pl, iv[0]); });
    par->setAttr("mir.lanes", ob.getUnitAttr());
  };

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
        Value desc;
        if (prefetch) {
          auto a = eb.create<memref::AllocaOp>(l, MemRefType::get({kDescSize}, eb.getI32Type()));
          a->setAttr("mir.scratch", eb.getUnitAttr());
          desc = a;
        }
        const Value mask = y.i32(0x7fffffff);

        // The element's entity descriptors, from the tables or from the prefetched copy.
        auto faceNbr = [&](Emit &z, Value f) {
          return prefetch ? z.load(desc, z.toIndex(z.add(f, z.i32(kDescFaceNbr))))
                          : z.load(d.faceNbr, z.add(z.mul(e, z.idx(6)), z.toIndex(f)));
        };
        auto faceCode = [&](Emit &z, Value f) {
          return prefetch ? z.load(desc, z.toIndex(z.add(f, z.i32(kDescFaceCode))))
                          : Value(z.b.create<arith::ExtUIOp>(z.loc, z.b.getI32Type(),
                                                             z.load(d.faceCode, z.add(z.mul(e, z.idx(6)), z.toIndex(f)))));
        };
        auto edgeEntry = [&](Emit &z, Value k) {
          return prefetch ? z.load(desc, z.toIndex(z.add(k, z.i32(kDescEdge))))
                          : z.load(d.edgeOf, z.add(z.mul(e, z.idx(12)), z.toIndex(k)));
        };
        // [begin, end) of the star of local edge k (own = its edge entry) or corner c.
        auto edgeStar = [&](Emit &z, Value k, Value own, Value &begin, Value &end) {
          if (prefetch) {
            const Value at = z.toIndex(z.add(z.mul(k, z.i32(2)), z.i32(kDescEdgeOff)));
            begin = z.load(desc, at);
            end = z.load(desc, z.add(at, z.idx(1)));
            return;
          }
          const Value g = z.toIndex(z.band(own, mask));
          begin = z.load(d.edgeOff, g);
          end = z.load(d.edgeOff, z.add(g, z.idx(1)));
        };
        auto vertStar = [&](Emit &z, Value c, Value &begin, Value &end) {
          if (prefetch) {
            const Value at = z.toIndex(z.add(z.mul(c, z.i32(2)), z.i32(kDescVertOff)));
            begin = z.load(desc, at);
            end = z.load(desc, z.add(at, z.idx(1)));
            return;
          }
          const Value g = z.toIndex(z.load(d.vertOf, z.add(z.mul(e, z.idx(8)), z.toIndex(c))));
          begin = z.load(d.vertOff, g);
          end = z.load(d.vertOff, z.add(g, z.idx(1)));
        };

        // Descriptor slot i (i < 26): face i, edge i - 6, corner i - 18, from the tables.
        auto loadDescriptor = [&](OpBuilder &ob, Location ol, Value i) {
          Emit z{ob, ol};
          auto ifFace = ob.create<scf::IfOp>(ol, z.lt(i, z.i32(6)), /*withElseRegion=*/true);
          {
            OpBuilder tb = ifFace.getThenBodyBuilder();
            Emit t{tb, ol};
            const Value at = t.add(t.mul(e, t.idx(6)), t.toIndex(i));
            t.store(t.load(d.faceNbr, at), desc, t.toIndex(t.add(i, t.i32(kDescFaceNbr))));
            t.store(tb.create<arith::ExtUIOp>(ol, tb.getI32Type(), t.load(d.faceCode, at)), desc,
                    t.toIndex(t.add(i, t.i32(kDescFaceCode))));
          }
          OpBuilder eb2 = ifFace.getElseBodyBuilder();
          Emit z2{eb2, ol};
          auto ifEdge = eb2.create<scf::IfOp>(ol, z2.lt(i, z2.i32(18)), /*withElseRegion=*/true);
          {
            OpBuilder tb = ifEdge.getThenBodyBuilder();
            Emit t{tb, ol};
            const Value k = t.sub(i, t.i32(6));
            const Value own = t.load(d.edgeOf, t.add(t.mul(e, t.idx(12)), t.toIndex(k)));
            const Value g = t.toIndex(t.band(own, mask));
            const Value at = t.toIndex(t.add(t.mul(k, t.i32(2)), t.i32(kDescEdgeOff)));
            t.store(own, desc, t.toIndex(t.add(k, t.i32(kDescEdge))));
            t.store(t.load(d.edgeOff, g), desc, at);
            t.store(t.load(d.edgeOff, t.add(g, t.idx(1))), desc, t.add(at, t.idx(1)));
          }
          {
            OpBuilder tb = ifEdge.getElseBodyBuilder();
            Emit t{tb, ol};
            const Value c = t.sub(i, t.i32(18));
            const Value g = t.toIndex(t.load(d.vertOf, t.add(t.mul(e, t.idx(8)), t.toIndex(c))));
            const Value at = t.toIndex(t.add(t.mul(c, t.i32(2)), t.i32(kDescVertOff)));
            t.store(t.load(d.vertOff, g), desc, at);
            t.store(t.load(d.vertOff, t.add(g, t.idx(1))), desc, t.add(at, t.idx(1)));
          }
        };

        const int64_t first = prefetch ? std::max<int64_t>(d.n3, kEntities) : d.n3;
        lanes(eb, l, first, [&](OpBuilder &pb, Location pl, Value node) {
          Emit z{pb, pl};
          auto stageOwn = [&](OpBuilder &sb, Location sl) {
            Emit t{sb, sl};
            t.store(t.load(d.values, t.add(base, node)), stage, d.slot(t, node));
          };
          if (first == d.n3) {
            stageOwn(pb, pl);
          } else {
            auto guard = pb.create<scf::IfOp>(pl, z.lt(node, z.idx(d.n3)), /*withElseRegion=*/false);
            OpBuilder gb = guard.getThenBodyBuilder();
            stageOwn(gb, pl);
          }
          if (prefetch) {
            const Value i = z.toI32(node);
            auto guard = pb.create<scf::IfOp>(pl, z.lt(i, z.i32(kEntities)), /*withElseRegion=*/false);
            OpBuilder gb = guard.getThenBodyBuilder();
            loadDescriptor(gb, pl, i);
          }
        });

        if (d.F > 0)
          lanes(eb, l, 6LL * d.F, [&](OpBuilder &pb, Location pl, Value pos) {
            Emit z{pb, pl};
            const Value p32 = z.toI32(pos), F = z.i32(d.F);
            const Value f = z.div(p32, F);
            const Value nbr = faceNbr(z, f);
            const Value code = faceCode(z, f);
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
            Value begin, end;
            vertStar(t, c, begin, end);
            const Value s = d.starSum(cb, cl, d.vertEnt, begin, end,
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
                const Value own = edgeEntry(t, k);
                const Value zero = t.i32(0), Pc = t.i32(P);
                Value begin, end;
                edgeStar(t, k, own, begin, end);
                // Canonical position on the edge: the star's entries are in the edge's
                // canonical direction, each with its own reversal bit.
                const Value tc = t.select(t.lt(own, zero), t.sub(Pc, pos1), pos1);
                const Value s = d.starSum(
                    tb, tl, d.edgeEnt, begin, end,
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

        lanes(eb, l, d.n3, [&](OpBuilder &pb, Location pl, Value node) {
          Emit z{pb, pl};
          z.store(z.load(stage, d.slot(z, node)), d.out, z.add(base, node));
        });
      });
  elementLoop->setAttr("mir.elements", b.getUnitAttr());
  op.erase();
  return success();
}

struct LowerDssPass : public PassWrapper<LowerDssPass, OperationPass<ModuleOp>> {
  MLIR_DEFINE_EXPLICIT_INTERNAL_INLINE_TYPE_ID(LowerDssPass)

  LowerDssPass() = default;
  LowerDssPass(const LowerDssPass &other) : PassWrapper(other) {}
  Option<bool> prefetch{*this, "prefetch-descriptors",
                        llvm::cl::desc("Load each element's descriptors in the first phase"), llvm::cl::init(true)};

  StringRef getArgument() const final { return "mir-lower-dss"; }
  StringRef getDescription() const final {
    return "Lower mir.dss to the staged gather: per-element parallel phases over "
           "own copies, face pairs, edge/corner stars and the write-back";
  }
  void getDependentDialects(DialectRegistry &registry) const override {
    registry.insert<arith::ArithDialect, memref::MemRefDialect, scf::SCFDialect>();
  }

  void runOnOperation() override {
    SmallVector<mir::DssOp> ops;
    getOperation()->walk([&](mir::DssOp op) { ops.push_back(op); });
    for (mir::DssOp op : ops)
      if (failed(lowerOne(op, prefetch)))
        return signalPassFailure();
  }
};

}  // namespace

namespace mir {
void registerLowerDssPass() { PassRegistration<LowerDssPass>(); }
}  // namespace mir
