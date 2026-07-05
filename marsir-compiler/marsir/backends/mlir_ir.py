"""MLIR backend: emit the operator as a mir-dialect partial-assembly chain.

Stage 2 twin of the source emitters. The SAME parsed operator spec that feeds
host_cpp / cuda here becomes textual MLIR in the `mir` dialect -- the
B^T . D . B sum-factorization skeleton:

  B   : mir.contract   -- sweep the field to the quadrature grid (per axis)
  D   : mir.flux        -- the authored pointwise flux, walked from the AST
  B^T : mir.contract   -- integrate back to modes

Only the flux region differs between operators; the contractions are the fixed
skeleton. The emitted text round-trips through the real dialect parser
(marsir-mlir/.../mir-opt), so it is the on-ramp to the mir -> linalg -> gpu
lowering path, not just a pretty dump.
"""

from ..ir import ElementApply  # noqa: F401  (kept for callers/type context)
from ..expr import Num, Var, Unary, Binary
from .common import AUTOGEN_BANNER

# AST binary op -> arith float op. Unary '-' maps to arith.negf separately.
_ARITH = {
    "+": "arith.addf",
    "-": "arith.subf",
    "*": "arith.mulf",
    "/": "arith.divf",
}


def _fmt_f64(text):
    # MLIR float literals need a decimal point (1e-6 does not lex, 1.0e-6 does).
    v = float(text)
    s = repr(v)
    if "e" in s or "E" in s:
        mant, _, exp = s.partition("e")
        if "." not in mant:
            mant += ".0"
        return mant + "e" + exp
    if "." not in s:
        s += ".0"
    return s


def _fresh(counter):
    name = "%%t%d" % counter[0]
    counter[0] += 1
    return name


def to_mlir_ssa(expr, valnames, lines, counter):
    """Recursively lower a flux AST node to SSA arith ops over f64 block args.

    valnames maps a flux-input name (deriv, g2, ...) to its region block-arg SSA
    name. Emits one `%tN = arith.<op> ... : f64` line per interior node and
    returns the SSA name holding the node's value.
    """
    if isinstance(expr, Num):
        name = _fresh(counter)
        lines.append("    %s = arith.constant %s : f64" % (name, _fmt_f64(expr.value)))
        return name
    if isinstance(expr, Var):
        # Flux inputs are block args; anything else is an unknown constant name,
        # which the mir dialect has no slot for -- reject early with a clear msg.
        if expr.name not in valnames:
            raise ValueError("flux references %r which is not a block input; "
                             "the mlir backend only supports flux-input names"
                             % expr.name)
        return valnames[expr.name]
    if isinstance(expr, Unary):
        operand = to_mlir_ssa(expr.operand, valnames, lines, counter)
        name = _fresh(counter)
        lines.append("    %s = arith.negf %s : f64" % (name, operand))
        return name
    if isinstance(expr, Binary):
        lhs = to_mlir_ssa(expr.lhs, valnames, lines, counter)
        rhs = to_mlir_ssa(expr.rhs, valnames, lines, counter)
        name = _fresh(counter)
        lines.append("    %s = %s %s, %s : f64" % (name, _ARITH[expr.op], lhs, rhs))
        return name
    raise TypeError("unknown expr node %r" % (expr,))


def emit(ea, p=7):
    o = ea.op
    n = p + 1
    t3 = "tensor<%dx%dx%dxf64>" % (n, n, n)   # the element field on the grid
    t2 = "tensor<%dx%dxf64>" % (n, n)         # a 1D reference operator matrix

    # Block args of the flux region, one per flux input, sorted for a stable ABI.
    inputs = sorted(ea.free_vars)
    metric_used = [g for g in ("g0", "g1", "g2") if g in ea.free_vars]

    # func arguments: field, the two always-present operator matrices, the
    # tangential-derivative matrix (only if the flux has cross terms), then one
    # 3D metric tensor per metric component the flux reads.
    args = ["%%u: %s" % t3, "%%Dtil: %s" % t2, "%%W: %s" % t2]
    if ea.needs_tangential:
        args.append("%%D: %s" % t2)
    for g in metric_used:
        args.append("%%%s: %s" % (g, t3))

    banner = AUTOGEN_BANNER.format(backend="mlir", name=o.name, flux=o.flux_src)

    lines = [banner.rstrip("\n"),
             "// RUN: mir-opt %s | mir-opt",
             "func.func @%s_pa(%s) -> %s {" % (o.name, ", ".join(args), t3)]

    # B sweeps: normal derivative (Dtil, axis 0) and the tangential derivatives
    # (D, axes 1/2). Each flux input name maps to the tensor that carries it.
    input_tensor = {}
    if ea.needs_deriv:
        lines.append("  %%deriv = mir.contract %%u, %%Dtil {axis = 0 : i64} : "
                     "(%s, %s) -> %s" % (t3, t2, t3))
        input_tensor["deriv"] = "%deriv"
    if ea.needs_dt1:
        lines.append("  %%dt1 = mir.contract %%u, %%D {axis = 1 : i64} : "
                     "(%s, %s) -> %s" % (t3, t2, t3))
        input_tensor["dt1"] = "%dt1"
    if ea.needs_dt2:
        lines.append("  %%dt2 = mir.contract %%u, %%D {axis = 2 : i64} : "
                     "(%s, %s) -> %s" % (t3, t2, t3))
        input_tensor["dt2"] = "%dt2"
    for g in metric_used:
        input_tensor[g] = "%" + g

    # D: the authored flux as a mir.flux region over per-point f64 block args.
    ins_tensors = ", ".join(input_tensor[name] for name in inputs)
    ins_types = ", ".join(t3 for _ in inputs)
    blk = ", ".join("%%in_%s: f64" % name for name in inputs)
    lines.append("  %%flux = mir.flux ins(%s) : (%s) -> %s {"
                 % (ins_tensors, ins_types, t3))
    lines.append("  ^bb0(%s):" % blk)

    valnames = dict((name, "%in_" + name) for name in inputs)
    counter = [0]
    result = to_mlir_ssa(o.flux, valnames, lines, counter)
    lines.append("    mir.yield %s : f64" % result)
    lines.append("  }")

    # B^T: integrate the flux back to modal coefficients.
    lines.append("  %%y = mir.contract %%flux, %%W {axis = 0 : i64} : "
                 "(%s, %s) -> %s" % (t3, t2, t3))
    lines.append("  return %%y : %s" % t3)
    lines.append("}")
    return "\n".join(lines) + "\n"


def _apply_body(L, ea, p, ctr, v, uval, gval, indent="  ", y_init=None):
    """Emit the per-element Knaus Alg-2 apply into L; returns the final y SSA
    name. The face loop l is a REAL scf.for with iter_args(y) (not unrolled):
    temporaries then exist once per dir and bufferize to a few reusable
    buffers instead of ~70 live allocas (the p=7 97KB/thread local-spill
    pathology of the unrolled form). dir stays unrolled: contraction `axis`
    is an attribute and cannot be dynamic."""
    o = ea.op
    n, P = p + 1, p
    I = indent
    t3 = "tensor<%dx%dx%dxf64>" % (n, n, n)
    t2 = "tensor<%dx%dxf64>" % (n, n)
    tPn = "tensor<%dx%dxf64>" % (P, n)
    tG = "tensor<3x%dx%dx%dx3xf64>" % (P, n, n)
    inputs = sorted(ea.free_vars)
    metric_used = [g for g in ("g0", "g1", "g2") if g in ea.free_vars]

    def fresh(tag):
        v[0] += 1
        return "%%%s%d" % (tag, v[0])

    def all_faces_type(d):
        dims = [n, n, n]
        dims[d] = P
        return "tensor<%dx%dx%dxf64>" % tuple(dims)

    def face_slice(J, src, src_ty, d, lvar):
        off = ["0", "0", "0"]; off[d] = lvar
        sz = [str(n)] * 3; sz[d] = "1"
        r = fresh("f")
        L.append("%s%s = tensor.extract_slice %s[%s] [%s] [1, 1, 1] : %s to %s"
                 % (J, r, src, ", ".join(off), ", ".join(sz), src_ty, t2))
        return r

    def g_slice(J, d, lvar, c):
        r = fresh("g")
        L.append("%s%s = tensor.extract_slice %s[%d, %s, 0, 0, %d] "
                 "[1, 1, %d, %d, 1] [1, 1, 1, 1, 1] : %s to %s"
                 % (J, r, gval, d, lvar, c, n, n, tG, t2))
        return r

    def contract2(J, src, mat, axis):
        r = fresh("c")
        L.append("%s%s = mir.contract %s, %s {axis = %d : i64} : "
                 "(%s, %s) -> %s" % (J, r, src, mat, axis, t2, t2, t2))
        return r

    def flux2(J, tensors, body_fn, names):
        r = fresh("x")
        tys = ", ".join(t2 for _ in tensors)
        blk = ", ".join("%%a%d: f64" % i for i in range(len(tensors)))
        L.append("%s%s = mir.flux ins(%s) : (%s) -> %s {"
                 % (J, r, ", ".join(tensors), tys, t2))
        L.append("%s^bb0(%s):" % (J, blk))
        valnames = dict(zip(names, ["%%a%d" % i for i in range(len(tensors))]))
        body_fn(valnames)
        L.append("%s}" % J)
        return r

    ze = fresh("ze")
    L.append("%s%s = arith.constant 0.0 : f64" % (I, ze))
    if y_init is None:
        y_init = fresh("ye")
        L.append("%s%s = tensor.empty() : %s" % (I, y_init, t3))
    y = fresh("y")
    L.append("%s%s = linalg.fill ins(%s : f64) outs(%s : %s) -> %s"
             % (I, y, ze, y_init, t3, t3))
    c0 = fresh("i")
    L.append("%s%s = arith.constant 0 : index" % (I, c0))
    c1 = fresh("i")
    L.append("%s%s = arith.constant 1 : index" % (I, c1))
    cP = fresh("i")
    L.append("%s%s = arith.constant %d : index" % (I, cP, P))

    for d in range(3):
        fa_ty = all_faces_type(d)
        interp_all = deriv_all = None
        if ea.needs_tangential:
            interp_all = fresh("ia")
            L.append("%s%s = mir.contract %s, %%Btil {axis = %d : i64} : "
                     "(%s, %s) -> %s" % (I, interp_all, uval, d, t3, tPn, fa_ty))
        if ea.needs_deriv:
            deriv_all = fresh("da")
            L.append("%s%s = mir.contract %s, %%Dtil {axis = %d : i64} : "
                     "(%s, %s) -> %s" % (I, deriv_all, uval, d, t3, tPn, fa_ty))

        lv = fresh("l")
        yloop = fresh("y")
        yiter = fresh("ya")
        L.append("%s%s = scf.for %s = %s to %s step %s iter_args(%s = %s) -> (%s) {"
                 % (I, yloop, lv, c0, cP, c1, yiter, y, t3))
        J = I + "  "

        tensors, names = [], []
        if ea.needs_deriv:
            tensors.append(face_slice(J, deriv_all, fa_ty, d, lv))
            names.append("deriv")
        if ea.needs_tangential:
            iface = face_slice(J, interp_all, fa_ty, d, lv)
            if ea.needs_dt1:
                tensors.append(contract2(J, iface, "%D", 0)); names.append("dt1")
            if ea.needs_dt2:
                tensors.append(contract2(J, iface, "%D", 1)); names.append("dt2")
        for g in metric_used:
            tensors.append(g_slice(J, d, lv, int(g[1])))
            names.append(g)
        order = [names.index(x) for x in inputs]
        tensors = [tensors[i] for i in order]

        def flux_body(valnames):
            res = to_mlir_ssa(o.flux, valnames, L, ctr)
            L.append("    mir.yield %s : f64" % res)
        fl = flux2(J, tensors, flux_body, inputs)

        tmp = contract2(J, fl, "%W", 1)
        intf = contract2(J, tmp, "%W", 0)

        lp1 = fresh("l")
        L.append("%s%s = arith.addi %s, %s : index" % (J, lp1, lv, c1))
        ycur = yiter
        for plane, opname in ((lv, "arith.subf"), (lp1, "arith.addf")):
            off = ["0", "0", "0"]; off[d] = plane
            sz = [str(n)] * 3; sz[d] = "1"
            cur = fresh("pl")
            L.append("%s%s = tensor.extract_slice %s[%s] [%s] [1, 1, 1] : "
                     "%s to %s" % (J, cur, ycur, ", ".join(off), ", ".join(sz),
                                   t3, t2))
            def upd_body(valnames, _op=opname):
                r = "%%t%d" % ctr[0]; ctr[0] += 1
                L.append("    %s = %s %%a0, %%a1 : f64" % (r, _op))
                L.append("    mir.yield %s : f64" % r)
            upd = flux2(J, [cur, intf], upd_body, ["cur", "intf"])
            ynew = fresh("y")
            L.append("%s%s = tensor.insert_slice %s into %s[%s] [%s] "
                     "[1, 1, 1] : %s into %s"
                     % (J, ynew, upd, ycur, ", ".join(off), ", ".join(sz),
                        t2, t3))
            ycur = ynew
        L.append("%sscf.yield %s : %s" % (J, ycur, t3))
        L.append("%s}" % I)
        y = yloop
    return y


def emit_full(ea, p=7):
    """Single-element full Knaus apply (tensor-semantic func)."""
    o = ea.op
    n, P = p + 1, p
    t3 = "tensor<%dx%dx%dxf64>" % (n, n, n)
    t2 = "tensor<%dx%dxf64>" % (n, n)
    tPn = "tensor<%dx%dxf64>" % (P, n)
    tG = "tensor<3x%dx%dx%dx3xf64>" % (P, n, n)
    metric_used = [g for g in ("g0", "g1", "g2") if g in ea.free_vars]

    args = ["%%u: %s" % t3, "%%Btil: %s" % tPn, "%%Dtil: %s" % tPn,
            "%%W: %s" % t2]
    if ea.needs_tangential:
        args.append("%%D: %s" % t2)
    if metric_used:
        args.append("%%G: %s" % tG)

    banner = AUTOGEN_BANNER.format(backend="mlir-full", name=o.name,
                                   flux=o.flux_src)
    L = [banner.rstrip("\n"),
         "// RUN: mir-opt %s | mir-opt",
         "func.func @%s_apply(%s) -> %s {" % (o.name, ", ".join(args), t3)]
    ctr, v = [0], [0]
    y = _apply_body(L, ea, p, ctr, v, "%u", "%G")
    L.append("  return %s : %s" % (y, t3))
    L.append("}")
    return "\n".join(L) + "\n"


def emit_full_batched(ea, p=7, tpb=128):
    """The FUSED batched form: one memref-level scf.parallel (blocks x threads)
    over elements -- thread t of block b applies the WHOLE operator to element
    e = b*tpb + t. mir stays tensor-semantic inside via the bufferization
    boundary ops: bufferization.to_tensor(restrict) on the element subviews in,
    materialize_in_destination out. Lowers to ONE gpu kernel (thread-per-
    element), unlike the naive path's one-launch-per-linalg-op."""
    o = ea.op
    n, P = p + 1, p
    n2, n3 = n * n, n * n * n
    t3 = "tensor<%dx%dx%dxf64>" % (n, n, n)
    t2 = "tensor<%dx%dxf64>" % (n, n)
    tPn = "tensor<%dx%dxf64>" % (P, n)
    tG = "tensor<3x%dx%dx%dx3xf64>" % (P, n, n)
    metric_used = [g for g in ("g0", "g1", "g2") if g in ea.free_vars]

    mU = "memref<?x%dx%dx%dxf64>" % (n, n, n)
    mUe = "memref<%dx%dx%dxf64, strided<[%d, %d, 1], offset: ?>>" % (n, n, n, n2, n)
    gdims = (P, n, n)
    gstr = (P * n * n * 3, n * n * 3, n * 3, 3)
    mG = "memref<?x3x%dx%dx%dx3xf64>" % gdims
    mGe = ("memref<3x%dx%dx%dx3xf64, strided<[%d, %d, %d, %d, 1], offset: ?>>"
           % (gdims + gstr))
    m2 = "memref<%dx%dxf64>" % (n, n)
    mPn = "memref<%dx%dxf64>" % (P, n)

    args = ["%%U: %s" % mU, "%%Y: %s" % mU, "%%Btilm: %s" % mPn,
            "%%Dtilm: %s" % mPn, "%%Wm: %s" % m2]
    if ea.needs_tangential:
        args.append("%%Dm: %s" % m2)
    if metric_used:
        args.append("%%Gm: %s" % mG)

    banner = AUTOGEN_BANNER.format(backend="mlir-full-batched", name=o.name,
                                   flux=o.flux_src)
    dm = "(d0) -> (d0)"
    L = [banner.rstrip("\n"),
         "func.func @%s_apply_batched(%s) {" % (o.name, ", ".join(args)),
         "  %c0 = arith.constant 0 : index",
         "  %c1 = arith.constant 1 : index",
         "  %%ctpb = arith.constant %d : index" % tpb,
         "  %%E = memref.dim %%U, %%c0 : %s" % mU,
         "  %B = arith.divui %E, %ctpb : index",
         "  scf.parallel (%b) = (%c0) to (%B) step (%c1) {",
         "    %be = arith.muli %b, %ctpb : index",
         "    scf.parallel (%t) = (%c0) to (%ctpb) step (%c1) {",
         "      %e = arith.addi %be, %t : index",
         "      %%us = memref.subview %%U[%%e, 0, 0, 0] [1, %d, %d, %d] "
         "[1, 1, 1, 1] : %s to %s" % (n, n, n, mU, mUe),
         "      %%u = bufferization.to_tensor %%us restrict : %s" % mUe,
         "      %%Btil = bufferization.to_tensor %%Btilm restrict : %s" % mPn,
         "      %%Dtil = bufferization.to_tensor %%Dtilm restrict : %s" % mPn,
         "      %%W = bufferization.to_tensor %%Wm restrict : %s" % m2]
    if ea.needs_tangential:
        L.append("      %%D = bufferization.to_tensor %%Dm restrict : %s" % m2)
    if metric_used:
        L.append("      %%gs = memref.subview %%Gm[%%e, 0, 0, 0, 0, 0] "
                 "[1, 3, %d, %d, %d, 3] [1, 1, 1, 1, 1, 1] : %s to %s"
                 % (P, n, n, mG, mGe))
        L.append("      %%G = bufferization.to_tensor %%gs restrict : %s" % mGe)

    L += ["      %%ys = memref.subview %%Y[%%e, 0, 0, 0] [1, %d, %d, %d] "
          "[1, 1, 1, 1] : %s to %s" % (n, n, n, mU, mUe),
          "      %%yseed = bufferization.to_tensor %%ys restrict writable : %s"
          % mUe]
    ctr, v = [0], [0]
    y = _apply_body(L, ea, p, ctr, v, "%u", "%G", indent="      ",
                    y_init="%yseed")
    L += ["      bufferization.materialize_in_destination %s in writable %%ys "
          ": (%s, %s) -> ()" % (y, t3, mUe),
          "      scf.reduce",
          "    } {mapping = [#gpu.loop_dim_map<processor = thread_x, "
          "map = %s, bound = %s>]}" % (dm, dm),
          "    scf.reduce",
          "  } {mapping = [#gpu.loop_dim_map<processor = block_x, "
          "map = %s, bound = %s>]}" % (dm, dm),
          "  return",
          "}"]
    return "\n".join(L) + "\n"
