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


def emit_full(ea, p=7):
    """The FULL Knaus Alg-2 apply as mir IR: per direction, ALL SCS faces at
    once via one rectangular contraction (Btil/Dtil are Pxn: the contracted
    axis changes size n->P), then per face l: tangential contractions, the
    authored flux with the per-(dir,l) metric slices, two W integrations, and
    the +/- scatter onto the two face-bounding planes of y. The (dir, l) loops
    are emitted UNROLLED -- generating that text is the front-end's job.

    G layout matches the host metric: tensor<3xPxnxnx3xf64>, g[dir,l,s,r,c]
    with c = (0:tang2/r, 1:tang1/s, 2:normal)."""
    o = ea.op
    n, P = p + 1, p
    t3 = "tensor<%dx%dx%dxf64>" % (n, n, n)
    t2 = "tensor<%dx%dxf64>" % (n, n)
    tPn = "tensor<%dx%dxf64>" % (P, n)
    tG = "tensor<3x%dx%dx%dx3xf64>" % (P, n, n)

    inputs = sorted(ea.free_vars)
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
         "func.func @%s_apply(%s) -> %s {" % (o.name, ", ".join(args), t3),
         "  %%y0 = arith.constant dense<0.0> : %s" % t3]

    ctr = [0]           # SSA temp counter for flux bodies
    v = [0]             # value counter for tensors

    def fresh(tag):
        v[0] += 1
        return "%%%s%d" % (tag, v[0])

    def all_faces_type(d):
        dims = [n, n, n]
        dims[d] = P
        return "tensor<%dx%dx%dxf64>" % tuple(dims)

    def face_slice(src, src_ty, d, l):
        """extract the l-th face (rank-reducing) along axis d."""
        off = ["0", "0", "0"]; off[d] = str(l)
        sz = [str(n)] * 3; sz[d] = "1"
        r = fresh("f")
        L.append("  %s = tensor.extract_slice %s[%s] [%s] [1, 1, 1] : %s to %s"
                 % (r, src, ", ".join(off), ", ".join(sz), src_ty, t2))
        return r

    def g_slice(d, l, c):
        r = fresh("g")
        L.append("  %s = tensor.extract_slice %%G[%d, %d, 0, 0, %d] "
                 "[1, 1, %d, %d, 1] [1, 1, 1, 1, 1] : %s to %s"
                 % (r, d, l, c, n, n, tG, t2))
        return r

    def contract2(src, mat, axis):
        r = fresh("c")
        L.append("  %s = mir.contract %s, %s {axis = %d : i64} : "
                 "(%s, %s) -> %s" % (r, src, mat, axis, t2, t2, t2))
        return r

    def flux2(tensors, body_fn, names):
        """rank-2 mir.flux over `tensors`; body built by body_fn(valnames)."""
        r = fresh("x")
        tys = ", ".join(t2 for _ in tensors)
        blk = ", ".join("%%a%d: f64" % i for i in range(len(tensors)))
        L.append("  %s = mir.flux ins(%s) : (%s) -> %s {"
                 % (r, ", ".join(tensors), tys, t2))
        L.append("  ^bb0(%s):" % blk)
        valnames = dict(zip(names, ["%%a%d" % i for i in range(len(tensors))]))
        body_fn(valnames)
        L.append("  }")
        return r

    y = "%y0"
    for d in range(3):
        fa_ty = all_faces_type(d)
        interp_all = deriv_all = None
        if ea.needs_tangential:
            interp_all = fresh("ia")
            L.append("  %s = mir.contract %%u, %%Btil {axis = %d : i64} : "
                     "(%s, %s) -> %s" % (interp_all, d, t3, tPn, fa_ty))
        if ea.needs_deriv:
            deriv_all = fresh("da")
            L.append("  %s = mir.contract %%u, %%Dtil {axis = %d : i64} : "
                     "(%s, %s) -> %s" % (deriv_all, d, t3, tPn, fa_ty))
        for l in range(P):
            tensors, names = [], []
            if ea.needs_deriv:
                tensors.append(face_slice(deriv_all, fa_ty, d, l))
                names.append("deriv")
            if ea.needs_tangential:
                iface = face_slice(interp_all, fa_ty, d, l)
                if ea.needs_dt1:
                    tensors.append(contract2(iface, "%D", 0)); names.append("dt1")
                if ea.needs_dt2:
                    tensors.append(contract2(iface, "%D", 1)); names.append("dt2")
            for g in metric_used:
                tensors.append(g_slice(d, l, int(g[1])))
                names.append(g)
            # reorder to the sorted flux-input order
            order = [names.index(x) for x in inputs]
            tensors = [tensors[i] for i in order]

            def flux_body(valnames):
                res = to_mlir_ssa(o.flux, valnames, L, ctr)
                L.append("    mir.yield %s : f64" % res)
            fl = flux2(tensors, flux_body, inputs)

            tmp = contract2(fl, "%W", 1)
            intf = contract2(tmp, "%W", 0)

            # y[plane l] -= intf ; y[plane l+1] += intf
            for plane, opname in ((l, "arith.subf"), (l + 1, "arith.addf")):
                off = ["0", "0", "0"]; off[d] = str(plane)
                sz = [str(n)] * 3; sz[d] = "1"
                cur = fresh("p")
                L.append("  %s = tensor.extract_slice %s[%s] [%s] [1, 1, 1] : "
                         "%s to %s" % (cur, y, ", ".join(off), ", ".join(sz),
                                       t3, t2))
                def upd_body(valnames, _op=opname):
                    r = "%%t%d" % ctr[0]; ctr[0] += 1
                    L.append("    %s = %s %%a0, %%a1 : f64" % (r, _op))
                    L.append("    mir.yield %s : f64" % r)
                upd = flux2([cur, intf], upd_body, ["cur", "intf"])
                ynew = fresh("y")
                L.append("  %s = tensor.insert_slice %s into %s[%s] [%s] "
                         "[1, 1, 1] : %s into %s"
                         % (ynew, upd, y, ", ".join(off), ", ".join(sz),
                            t2, t3))
                y = ynew

    L.append("  return %s : %s" % (y, t3))
    L.append("}")
    return "\n".join(L) + "\n"
