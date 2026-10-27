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
from ..expr import Num, Var, Unary, Binary, parse_expr
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


def _apply_body(L, ea, p, ctr, v, uval, gval, indent="  ", y_init=None,
                pad_faces=False):
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
    # pad_faces: Btil/Dtil arrive with a zero (P+1)-th row, so every sweep tile
    # is a full n-row tile. A P-row tile needs out-of-bounds guards on its last
    # row, and LLVM sinks a tensor-core mma whose result only that guarded store
    # uses into the guard -- undefined behavior, since every lane must execute
    # mma.sync. The extra plane is computed (in tiles already n rows tall) and
    # never read: the face loop still runs over the P real faces.
    nf = n if pad_faces else P
    tPn = "tensor<%dx%dxf64>" % (nf, n)
    tG = "tensor<3x%dx3x%dx%dxf64>" % (P, n, n)
    inputs = sorted(ea.free_vars)
    metric_used = [g for g in ("g0", "g1", "g2") if g in ea.free_vars]

    def fresh(tag):
        v[0] += 1
        return "%%%s%d" % (tag, v[0])

    def all_faces_type(d):
        dims = [n, n, n]
        dims[d] = nf
        return "tensor<%dx%dx%dxf64>" % tuple(dims)

    def face_slice(J, src, src_ty, d, lvar):
        off = ["0", "0", "0"]; off[d] = lvar
        sz = [str(n)] * 3; sz[d] = "1"
        r = fresh("f")
        L.append("%s%s = tensor.extract_slice %s[%s] [%s] [1, 1, 1] : %s to %s"
                 % (J, r, src, ", ".join(off), ", ".join(sz), src_ty, t2))
        return r

    # The metric is COMPONENT-MAJOR, [dir][face][component][row][col]: each
    # component plane is contiguous, so a warp reads it as one 16-byte load per
    # lane. Interleaving the components (component last) made every component
    # read two scalar loads spread over three times the cache lines, and the
    # tensor-core kernel is bound by the load/store unit, not by DRAM.
    def g_slice(J, d, lvar, c):
        r = fresh("g")
        L.append("%s%s = tensor.extract_slice %s[%d, %s, %d, 0, 0] "
                 "[1, 1, 1, %d, %d] [1, 1, 1, 1, 1] : %s to %s"
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

    # Direction 2 first: its Y planes are element-strided (each access touches
    # every cache line of the element), and the first direction only writes Y
    # -- its reads see the zero fill -- so the strided reads disappear. The
    # three directions are independent sums; only the rounding order changes.
    for d in (2, 1, 0):
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


def emit_full(ea, p=7, pad_faces=False, real="f64"):
    """Single-element full Knaus apply (tensor-semantic func). pad_faces: see
    _apply_body; Btil/Dtil are then n x n with a zero last row. real: the scalar
    type of every tensor and constant (f64, or f32 for a reduced-precision
    preconditioner)."""
    o = ea.op
    n, P = p + 1, p
    t3 = "tensor<%dx%dx%dxf64>" % (n, n, n)
    t2 = "tensor<%dx%dxf64>" % (n, n)
    tPn = "tensor<%dx%dxf64>" % (n if pad_faces else P, n)
    tG = "tensor<3x%dx3x%dx%dxf64>" % (P, n, n)
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
    y = _apply_body(L, ea, p, ctr, v, "%u", "%G", pad_faces=pad_faces)
    L.append("  return %s : %s" % (y, t3))
    L.append("}")
    text = "\n".join(L) + "\n"
    return text if real == "f64" else text.replace("f64", real)


def emit_galerkin(ga, p=7, q=None, real="f64"):
    """A form-galerkin operator (marsir.ir.GalerkinApply) on one element as mir: u
    (n^3 GLL values, n = p + 1) -> B along each axis to the q^3 Gauss points (q =
    p + 2 by default, the CEED benchmark choice) -> the reference gradient, D along
    each axis, for the components the outputs read -> the pointwise outputs, one
    mir.flux -> back through the transposed sweeps of the outputs given:
    w = v + sum_d Dt along d of dv_d, then y = Bt along each axis of w.
    Arguments: u, B (q x n), Bt (n x q), then D and Dt (q x q) when a gradient is
    read, then G (K x q x q x q: the K geometry components read). real: the scalar
    type (f64, or f32).

    geometry = trilinear: instead of G, gq and wq (q: the Gauss points and weights on
    [-1, 1]) and X (3 x 2 x 2 x 2: the element's corners by component, then corner bit
    along each axis). A scalar flux first turns the corners into the 21 coefficients
    the Jacobian needs: x_a = sum over masks m of c[a][m] prod_{d in m} xi_d, and J reads
    the 7 with m != 0 (the mixed ones are shared by the three derivatives). The point
    flux takes those and the coordinate and weight along each axis (dims broadcasts
    them) and forms, in registers, J[a][b] (3 FMAs each), the adjugate and det J, and
    from them the inputs the spec reads: wdetj = w det J and m_de = w adj adj^T / det
    (= w det J J^-1 J^-T). Nothing per Gauss point is stored or staged."""
    o = ga.op
    n = p + 1
    q = q or p + 2

    def ty(dims):
        return "tensor<%sxf64>" % "x".join(str(d) for d in dims)

    tq = ty([q, q, q])
    grads = [d for d in range(3) if "du%d" % d in ga.free_vars]
    args = ["%%u: %s" % ty([n, n, n]), "%%B: %s" % ty([q, n]), "%%Bt: %s" % ty([n, q])]
    if grads:
        args += ["%%D: %s" % ty([q, q]), "%%Dt: %s" % ty([q, q])]
    if ga.geometry:
        args.append("%%G: %s" % ty([ga.geometry, q, q, q]))
    if ga.trilinear:
        args += ["%%gq: %s" % ty([q]), "%%wq: %s" % ty([q]), "%%X: %s" % ty([3, 2, 2, 2])]
    banner = AUTOGEN_BANNER.format(backend="mlir-galerkin", name=o.name, flux=o.flux_src)
    L = [banner.rstrip("\n"), "// RUN: mir-opt %s | mir-opt",
         "func.func @%s_apply(%s) -> %s {" % (o.name, ", ".join(args), ty([n, n, n]))]
    v = [0]

    def fresh(tag):
        v[0] += 1
        return "%%%s%d" % (tag, v[0])

    def contract(x, dims, mat, rows, axis):
        out = list(dims)
        out[axis] = rows
        r = fresh("c")
        L.append("  %s = mir.contract %s, %s {axis = %d : i64} : (%s, %s) -> %s"
                 % (r, x, mat, axis, ty(dims), ty([rows, dims[axis]]), ty(out)))
        return r, out

    def sweep(x, dims, mat, rows):
        for axis in range(3):
            x, dims = contract(x, dims, mat, rows, axis)
        return x

    uq = sweep("%u", [n, n, n], "%B", q)
    inputs = []   # (name, SSA value)
    if ga.needs_u:
        inputs.append(("u", uq))
    for d in grads:
        inputs.append(("du%d" % d, contract(uq, [q, q, q], "%D", q, d)[0]))
    for k in range(ga.geometry):
        g = fresh("g")
        L.append("  %s = tensor.extract_slice %%G[%d, 0, 0, 0] [1, %d, %d, %d] [1, 1, 1, 1] : %s to %s"
                 % (g, k, q, q, q, ty([ga.geometry, q, q, q]), tq))
        inputs.append(("g%d" % k, g))
    geo_names = sorted(ga.free_vars & ({"wdetj"} | {"m%d%d" % (d, e) for d in range(3) for e in range(d, 3)}))
    dims = [[0, 1, 2] for _ in inputs]   # grid dimensions each flux input varies along
    if ga.trilinear and geo_names:
        for d in range(3):
            inputs += [("xi%d" % d, "%gq"), ("w%d" % d, "%wq")]
            dims += [[d], [d]]
        # The element's scalars come first, ahead of every sweep, so they do not split
        # the point flux from the contractions it fuses with.
        corners, pre = [], []
        for a in range(3):
            for c in range(8):
                x = fresh("x")
                pre.append("  %s = tensor.extract_slice %%X[%d, %d, %d, %d] [1, 1, 1, 1] [1, 1, 1, 1] : %s to tensor<f64>"
                           % (x, a, c >> 2, (c >> 1) & 1, c & 1, ty([3, 2, 2, 2])))
                corners.append(x)
        # c[a][m] = 1/8 sum over corners of X[a][corner] prod_{d in m} s_d, s_d = -1 or +1
        # for the corner's bit along axis d (bit d of m: xi_d is a factor).
        coef = [(a, m) for a in range(3) for m in range(1, 8)]
        cbody, cnames, cctr = [], {"X%d%d" % (a, c): "%%a%d" % (8 * a + c) for a in range(3) for c in range(8)}, [0]
        cyields = []
        for a, m in coef:
            terms = []
            for c in range(8):
                bits = (c >> 2, (c >> 1) & 1, c & 1)
                sign = 1
                for d in range(3):
                    if m >> d & 1 and not bits[d]:
                        sign = -sign
                terms.append(("+ " if sign > 0 else "- ") + "X%d%d" % (a, c))
            cyields.append(to_mlir_ssa(parse_expr("0.125*(%s)" % " ".join(terms).lstrip("+ ")), cnames, cbody, cctr))
        cx = fresh("x")
        pre.append("  %s:%d = mir.flux {dims = [%s]} ins(%s) : (%s) -> (%s) {"
                   % (cx, len(coef), ", ".join(["array<i64>"] * 24), ", ".join(corners), ", ".join(["tensor<f64>"] * 24),
                      ", ".join(["tensor<f64>"] * len(coef))))
        pre.append("  ^bb0(%s):" % ", ".join("%%a%d: f64" % i for i in range(24)))
        pre.extend(cbody)
        pre.append("    mir.yield %s : %s" % (", ".join(cyields), ", ".join(["f64"] * len(coef))))
        pre.append("  }")
        L[3:3] = pre
        for i, (a, m) in enumerate(coef):
            inputs.append(("c%d%d" % (a, m), "%s#%d" % (cx, i)))
            dims.append([])
    names = [nm for nm, _ in inputs]
    outs = [nm for nm, _ in ga.outputs]
    x = fresh("x")
    head = "%s:%d" % (x, len(outs)) if len(outs) > 1 else x
    res_ty = "(%s)" % ", ".join([tq] * len(outs)) if len(outs) > 1 else tq
    in_tys = [ty([q] * len(dm)) if dm else "tensor<f64>" for dm in dims]
    attr = ("{dims = [%s]} " % ", ".join("array<i64%s>" % ((": " + ", ".join(map(str, dm))) if dm else "")
                                          for dm in dims)) if any(dm != [0, 1, 2] for dm in dims) else ""
    L.append("  %s = mir.flux %sins(%s) : (%s) -> %s {"
             % (head, attr, ", ".join(val for _, val in inputs), ", ".join(in_tys), res_ty))
    L.append("  ^bb0(%s):" % ", ".join("%%a%d: f64" % i for i in range(len(inputs))))
    valnames = {nm: "%%a%d" % i for i, nm in enumerate(names)}
    ctr = [0]
    body = []
    if ga.trilinear and geo_names:
        # J[a][b] = d x_a / d xi_b in Horner form, o1 < o2 the other axes:
        # (c[b] + c[b,o1] xi_o1) + xi_o2 (c[b,o2] + c[b,o1,o2] xi_o1).
        valnames["w"] = to_mlir_ssa(parse_expr("w0*w1*w2"), valnames, body, ctr)
        for a in range(3):
            for b in range(3):
                o1, o2 = [ax for ax in range(3) if ax != b]
                mb, m1, m2 = 1 << b, 1 << o1, 1 << o2
                text = "(c%d%d + c%d%d*xi%d) + xi%d*(c%d%d + c%d7*xi%d)" % (
                    a, mb, a, mb | m1, o1, o2, a, mb | m2, a, o1)
                valnames["J%d%d" % (a, b)] = to_mlir_ssa(parse_expr(text), valnames, body, ctr)
        # The adjugate (adj = det J^-1), det J, then the geometry the spec reads.
        adj = {(0, 0): "J11*J22 - J12*J21", (0, 1): "J02*J21 - J01*J22", (0, 2): "J01*J12 - J02*J11",
               (1, 0): "J12*J20 - J10*J22", (1, 1): "J00*J22 - J02*J20", (1, 2): "J02*J10 - J00*J12",
               (2, 0): "J10*J21 - J11*J20", (2, 1): "J01*J20 - J00*J21", (2, 2): "J00*J11 - J01*J10"}
        for (d, r), text in adj.items():
            valnames["adj%d%d" % (d, r)] = to_mlir_ssa(parse_expr(text), valnames, body, ctr)
        valnames["det"] = to_mlir_ssa(parse_expr("J00*adj00 + J01*adj10 + J02*adj20"), valnames, body, ctr)
        for name in geo_names:
            if name == "wdetj":
                text = "w*det"
            else:
                d, e = int(name[1]), int(name[2])
                text = "w*(%s)/det" % " + ".join("adj%d%d*adj%d%d" % (d, r, e, r) for r in range(3))
            valnames[name] = to_mlir_ssa(parse_expr(text), valnames, body, ctr)
    yields = [to_mlir_ssa(expr, valnames, body, ctr) for _, expr in ga.outputs]
    L.extend(body)
    L.append("    mir.yield %s : %s" % (", ".join(yields), ", ".join(["f64"] * len(yields))))
    L.append("  }")
    result = {nm: (x if len(outs) == 1 else "%s#%d" % (x, i)) for i, nm in enumerate(outs)}

    w = result.get("v")
    for d in range(3):
        if "dv%d" % d not in result:
            continue
        t, _ = contract(result["dv%d" % d], [q, q, q], "%Dt", q, d)
        if w is None:
            w = t
            continue
        s = fresh("w")
        L.append("  %s = mir.flux ins(%s, %s) : (%s, %s) -> %s {" % (s, w, t, tq, tq, tq))
        L.append("  ^bb0(%a0: f64, %a1: f64):")
        L.append("    %t = arith.addf %a0, %a1 : f64")
        L.append("    mir.yield %t : f64")
        L.append("  }")
        w = s
    y = sweep(w, [q, q, q], "%Bt", n)
    L.append("  return %s : %s" % (y, ty([n, n, n])))
    L.append("}")
    text = "\n".join(L) + "\n"
    return text if real == "f64" else text.replace("f64", real)


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
    tG = "tensor<3x%dx3x%dx%dxf64>" % (P, n, n)
    metric_used = [g for g in ("g0", "g1", "g2") if g in ea.free_vars]

    mU = "memref<?x%dx%dx%dxf64>" % (n, n, n)
    mUe = "memref<%dx%dx%dxf64, strided<[%d, %d, 1], offset: ?>>" % (n, n, n, n2, n)
    gdims = (P, n, n)
    gstr = (P * 3 * n * n, 3 * n * n, n * n, n)
    mG = "memref<?x3x%dx3x%dx%dxf64>" % gdims
    mGe = ("memref<3x%dx3x%dx%dxf64, strided<[%d, %d, %d, %d, 1], offset: ?>>"
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
         # Ceil-divide: a floor here silently drops the last E mod tpb elements.
         "  %ctpbm1 = arith.subi %ctpb, %c1 : index",
         "  %Epad = arith.addi %E, %ctpbm1 : index",
         "  %B = arith.divui %Epad, %ctpb : index",
         "  scf.parallel (%b) = (%c0) to (%B) step (%c1) {",
         "    %be = arith.muli %b, %ctpb : index",
         "    scf.parallel (%t) = (%c0) to (%ctpb) step (%c1) {",
         "      %e = arith.addi %be, %t : index",
         # ...and the last block is then partly out of range, so guard it.
         "      %inb = arith.cmpi ult, %e, %E : index",
         "      scf.if %inb {",
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
                 "[1, 3, %d, 3, %d, %d] [1, 1, 1, 1, 1, 1] : %s to %s"
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
          "      }",
          "      scf.reduce",
          "    } {mapping = [#gpu.loop_dim_map<processor = thread_x, "
          "map = %s, bound = %s>]}" % (dm, dm),
          "    scf.reduce",
          "  } {mapping = [#gpu.loop_dim_map<processor = block_x, "
          "map = %s, bound = %s>]}" % (dm, dm),
          "  return",
          "}"]
    return "\n".join(L) + "\n"


# The cell-wise DSS kernel: direct stiffness summation on element-local vectors
# (marsir-mlir: mir.dss, lowered by --mir-lower-dss). The arguments are the flat
# device arrays of MARS's UnstructuredTopology, in the order of mir.dss.
DSS_ARGS = [("values", "f64"), ("ghost", "f64"), ("elements", "i32"),
            ("face_nbr", "i32"), ("face_code", "i8"),
            ("edge_of", "i32"), ("edge_off", "i32"), ("edge_ent", "i32"),
            ("vert_of", "i32"), ("vert_off", "i32"), ("vert_ent", "i32"),
            ("out", "f64")]


def emit_dss(p=7, real="f64"):
    a = {name: f"%{name}" for name, _ in DSS_ARGS}
    args = [(name, real if t == "f64" else t) for name, t in DSS_ARGS]
    sig = ", ".join(f"%{name}: memref<?x{t}>" for name, t in args)
    types = ", ".join(f"memref<?x{t}>" for _, t in args)
    name = f"dss_p{p}" + ("" if real == "f64" else "_" + real)
    return "\n".join([
        # mir.noalias: no two arguments overlap. The gather reads neighbours'
        # copies while it writes `out`, so it cannot run in place anyway.
        f"func.func @{name}({sig}) attributes {{mir.kernel, mir.noalias}} {{",
        f"  mir.dss ins({a['values']}, {a['ghost']}) elements({a['elements']})",
        f"          faces({a['face_nbr']}, {a['face_code']})",
        f"          edges({a['edge_of']}, {a['edge_off']}, {a['edge_ent']})",
        f"          vertices({a['vert_of']}, {a['vert_off']}, {a['vert_ent']})",
        f"          outs({a['out']}) {{degree = {p} : i64}} : {types}",
        "  return",
        "}",
        ""])


def emit_block_dss(p, block):
    """The block-structured DSS (MARS's mars_cellwise_block_dss.hpp, the block-node form)
    as two kernels, run in this order on one stream:
      block_sum_p{p}_b{B}(values, out)          every block node: out at its copies in
                                                the block = the sum of values over them
      block_skeleton_p{p}_b{B}(out, face_nbr, face_code, edge_off, edge_ent, vert_off,
                               vert_ent)        the nodes shared between blocks, in place,
                                                from one copy per block
    The tables are MARS's UnstructuredTopology built over the blocks as hexes."""
    tag = f"p{p}_b{block}"
    attrs = f"{{degree = {p} : i64, block = {block} : i64}}"
    return "\n".join([
        # mir.noalias: values and out are different arrays; the skeleton tables are
        # read-only and never overlap the copies.
        f"func.func @block_sum_{tag}(%values: memref<?xf64>, %out: memref<?xf64>) "
        "attributes {mir.kernel, mir.noalias} {",
        f"  mir.block_sum ins(%values) outs(%out) {attrs} : memref<?xf64>, memref<?xf64>",
        "  return",
        "}",
        f"func.func @block_skeleton_{tag}(%values: memref<?xf64>, %face_nbr: memref<?xi32>, "
        "%face_code: memref<?xi8>, %edge_off: memref<?xi32>, %edge_ent: memref<?xi32>, "
        "%vert_off: memref<?xi32>, %vert_ent: memref<?xi32>) attributes {mir.kernel, mir.noalias} {",
        "  mir.block_skeleton ins(%values) faces(%face_nbr, %face_code) edges(%edge_off, %edge_ent)",
        f"          vertices(%vert_off, %vert_ent) {attrs}",
        "          : memref<?xf64>, memref<?xi32>, memref<?xi8>, memref<?xi32>, memref<?xi32>, memref<?xi32>, "
        "memref<?xi32>",
        "  return",
        "}",
        ""])


def emit_dss_precond(p=7, az=True, zz=True):
    """The Jacobi preconditioner of the cell-wise solver as one mir.dss with an
    epilogue: z = 0 on Dirichlet copies, DSS(q) / diag elsewhere, and with az / zz
    the per-element partials of <a, z> and <z, z>, each node counted once (the
    epilogue of unstructured_precondition_kernel in mars_cellwise_topology.hpp)."""
    base = [n for n, _ in DSS_ARGS if n != "out"]
    flags = [("counting", "i32"), ("edge_dirichlet", "i8"), ("vert_dirichlet", "i8")]
    copies = [("diag", "f64")] + ([("a", "f64")] if az else [])
    nred = int(az) + int(zz)
    outs = [("z", "f64")] + ([("partials", "f64")] if nred else [])
    args = [(n, t) for n, t in DSS_ARGS if n != "out"] + flags + copies + outs
    sig = ", ".join(f"%{n}: memref<?x{t}>" for n, t in args)
    types = ", ".join(f"memref<?x{t}>" for _, t in args)
    name = f"precond_p{p}" + ("_az" if az else "") + ("_zz" if zz else "")
    body = ["  %zero = arith.constant 0.0 : f64",
            "  %one = arith.constant 1.0 : f64",
            "  %q = arith.divf %sum, %diag_t : f64",
            "  %v = arith.select %dir, %zero, %q : f64",
            "  %w = arith.select %cnt, %one, %zero : f64"]
    ys = ["%v"]
    if az:
        body += ["  %av = arith.mulf %a_t, %v : f64", "  %avw = arith.mulf %av, %w : f64"]
        ys.append("%avw")
    if zz:
        body += ["  %vv = arith.mulf %v, %v : f64", "  %vvw = arith.mulf %vv, %w : f64"]
        ys.append("%vvw")
    bargs = "%sum: f64, %dir: i1, %cnt: i1, %diag_t: f64" + (", %a_t: f64" if az else "")
    partials = " partials(%partials)" if nred else ""
    return "\n".join([
        f"func.func @{name}({sig}) attributes {{mir.kernel, mir.noalias}} {{",
        f"  mir.dss ins(%values, %ghost) elements(%elements) faces(%face_nbr, %face_code)",
        f"          edges(%edge_of, %edge_off, %edge_ent) vertices(%vert_of, %vert_off, %vert_ent)",
        f"          flags(%counting, %edge_dirichlet, %vert_dirichlet)",
        f"          copies({', '.join('%' + n for n, _ in copies)}) outs(%z){partials}",
        f"          {{degree = {p} : i64, reductions = {nred} : i64}} : {types} epilogue {{",
        f"  ^bb0({bargs}):",
        *["  " + line for line in body],
        f"    mir.yield {', '.join(ys)} : {', '.join(['f64'] * len(ys))}",
        "  }",
        "  return",
        "}",
        ""])


def emit_dss_precond_all(p=7):
    """All four variants the Krylov solvers call (no dots, <a, z>, <z, z>, both), in
    one module: precond_p{p}, precond_p{p}_az, precond_p{p}_zz, precond_p{p}_az_zz."""
    return "".join(emit_dss_precond(p, az, zz) for az in (False, True) for zz in (False, True))


# ---- BiCGStab on element-local vectors (bicgstab of mars_cellwise_krylov.hpp) -----
# Krylov scalars live in a device array, in the order of its enum Scalar.
KRYLOV = {"rho": 0, "alpha": 1, "omega": 2, "beta": 3, "rr": 4, "t0": 5, "t1": 6}


def _scalar(name, sc="%sc"):
    return [f"  %{name}_i = arith.constant {KRYLOV[name]} : index",
            f"  %{name} = memref.load {sc}[%{name}_i] : memref<?xf64>"]


def _cellwise(name, p, inputs, outs, body, yields, counting=False, reductions=0, scalars=()):
    """A kernel of one mir.cellwise op. inputs / outs name f64 arrays (an array may be
    both); body lines see %counted and %<name>_t for every input."""
    args = ["elements"] + (["counting"] if counting else []) + list(dict.fromkeys(inputs + outs)) + ["sc"] \
        + (["partials"] if reductions else [])
    types = {"elements": "i32", "counting": "i32"}
    sig = ", ".join(f"%{a}: memref<?x{types.get(a, 'f64')}>" for a in args)
    opnds = ["%elements"] + (["%counting"] if counting else []) + [f"%{a}" for a in inputs + outs] \
        + (["%partials"] if reductions else [])
    otypes = ["memref<?xi32>"] + (["memref<?xi32>"] if counting else []) + ["memref<?xf64>"] * (
        len(inputs) + len(outs) + (1 if reductions else 0))
    lines = [f"func.func @{name}({sig}) attributes {{mir.kernel}} {{"]
    for s_ in scalars:
        lines += _scalar(s_)
    lines += [f"  mir.cellwise elements(%elements)" + (" counting(%counting)" if counting else "")
              + (f" inputs({', '.join('%' + a for a in inputs)})" if inputs else "")
              + (f" outs({', '.join('%' + a for a in outs)})" if outs else "")
              + (" partials(%partials)" if reductions else "")
              + f" {{degree = {p} : i64, reductions = {reductions} : i64}} : {', '.join(otypes)} body {{",
              "  ^bb0(%counted: i1" + "".join(f", %{a}_t: f64" for a in inputs) + "):"]
    lines += ["    " + l for l in body]
    lines += [f"    mir.yield {', '.join(yields)} : {', '.join(['f64'] * len(yields))}", "  }", "  return", "}", ""]
    return "\n".join(lines)


def _reduce_block(n, op="arith.addf"):
    return ", ".join("{\n      ^bb0(%l: f64, %r: f64):\n        %x = " + op
                     + " %l, %r : f64\n        scf.reduce.return %x : f64\n      }" for _ in range(n))


def _sum_range(nv, lo, hi, src="%src"):
    """Lines that sum src[i * nv + j] over i in [lo, hi) as one {mir.lanes} phase into
    %s#j (one lane per i, a fixed-order tree across lanes)."""
    res = "%s" if nv == 1 else f"%s:{nv}"
    loads = []
    for j in range(nv):
        loads += [f"      %k{j} = arith.addi %k, %c{j} : index",
                  f"      %v{j} = memref.load {src}[%k{j}] : memref<?xf64>"]
    zeros = ", ".join(["%zero"] * nv)
    return [f"    {res} = scf.parallel (%i) = ({lo}) to ({hi}) step (%c1) init ({zeros}) -> ({', '.join(['f64'] * nv)}) {{",
            f"      %k = arith.muli %i, %nv : index", *loads,
            f"      scf.reduce({', '.join(f'%v{j}' for j in range(nv))} : {', '.join(['f64'] * nv)}) {_reduce_block(nv)}",
            "    } {mir.lanes}"]


def _sums(nv):
    return ["%s"] if nv == 1 else [f"%s#{j}" for j in range(nv)]


def _consts(nv):
    return ["  %c0 = arith.constant 0 : index", "  %c1 = arith.constant 1 : index",
            f"  %nv = arith.constant {nv} : index", "  %zero = arith.constant 0.0 : f64"] + \
        [f"  %c{j} = arith.constant {j} : index" for j in range(2, nv)]


# What a finished reduction does with its totals v (apply_step of the krylov header).
def _step_lines(step, v):
    st = lambda name, val: [f"    %{name}_si = arith.constant {KRYLOV[name]} : index",
                             f"    memref.store {val}, %sc[%{name}_si] : memref<?xf64>"]
    ld = lambda name: [f"    %{name}_li = arith.constant {KRYLOV[name]} : index",
                       f"    %{name}_old = memref.load %sc[%{name}_li] : memref<?xf64>"]
    if step == "start":
        return (["    %one = arith.constant 1.0 : f64", "    %zf = arith.constant 0.0 : f64"]
                + st("rr", v[0]) + st("rho", v[0]) + st("alpha", "%one") + st("omega", "%one") + st("beta", "%zf"))
    if step == "alpha":
        return ld("rho") + [f"    %a = arith.divf %rho_old, {v[0]} : f64"] + st("alpha", "%a")
    if step == "omega":
        return [f"    %w = arith.divf {v[0]}, {v[1]} : f64"] + st("omega", "%w")
    if step == "rho":   # beta = (v1 / rho) * (alpha / omega), then rho = v1
        return (ld("rho") + ld("alpha") + ld("omega")
                + [f"    %q1 = arith.divf {v[1]}, %rho_old : f64", "    %q2 = arith.divf %alpha_old, %omega_old : f64",
                   "    %beta_new = arith.mulf %q1, %q2 : f64"]
                + st("rr", v[0]) + st("beta", "%beta_new") + st("rho", v[1]))
    raise ValueError(step)


STEP_NV = {"start": 1, "alpha": 1, "omega": 2, "rho": 2}


def emit_bicgstab(p=7):
    """The kernels of a cell-wise BiCGStab besides the operator and the preconditioner:
      bicg_update_p_p{p}   p = r + beta (p - omega v)
      bicg_update_s_p{p}   s = r - alpha v
      bicg_update_xr_p{p}  x += alpha p + omega s, r = s - omega t, per-element partials
                           of <r, r> and <rh, r> (each node counted once)
      fold_nv{1,2}         per-element partials -> one partial per block (a contiguous
                           range each, fixed order); grid = size(dst) / nv
      finish_<step>        block partials -> totals -> the Krylov step (one block)
      totals_nv{1,2}       block partials -> totals in sc[t0], sc[t1] (for MPI)
      step_<step>          sc[t0], sc[t1] -> the Krylov step (after MPI)
    and for the right-preconditioned form (A M y = b, x = M y; the residual is the
    assembled b - A x, so an FP32 V-cycle can be M, see cellwise_mixed_ref.py):
      rbicg_update_p_p{p}  p = r + beta (p - omega v), and pw = p on the counting copy
                           of each node, 0 elsewhere (an unassembled form for M)
      rbicg_update_s_p{p}  s = r - alpha v, and sw likewise
      rbicg_update_xr_p{p} x += alpha ph + omega sh (ph = M pw, sh = M sw), r = s - omega t,
                           per-element partials of <r, r> and <rh, r>"""
    out = []
    out.append(_cellwise(f"bicg_update_p_p{p}", p, ["r", "v", "p"], ["p"],
                         ["%ov = arith.mulf %omega, %v_t : f64", "%d = arith.subf %p_t, %ov : f64",
                          "%bd = arith.mulf %beta, %d : f64", "%pn = arith.addf %r_t, %bd : f64"],
                         ["%pn"], scalars=("beta", "omega")))
    out.append(_cellwise(f"bicg_update_s_p{p}", p, ["r", "v"], ["s"],
                         ["%av = arith.mulf %alpha, %v_t : f64", "%sn = arith.subf %r_t, %av : f64"],
                         ["%sn"], scalars=("alpha",)))
    out.append(_cellwise(f"bicg_update_xr_p{p}", p, ["x", "p", "s", "t", "rh"], ["x", "r"],
                         ["%ap = arith.mulf %alpha, %p_t : f64", "%os = arith.mulf %omega, %s_t : f64",
                          "%inc = arith.addf %ap, %os : f64", "%xn = arith.addf %x_t, %inc : f64",
                          "%ot = arith.mulf %omega, %t_t : f64", "%rn = arith.subf %s_t, %ot : f64",
                          "%one = arith.constant 1.0 : f64", "%zw = arith.constant 0.0 : f64",
                          "%w = arith.select %counted, %one, %zw : f64",
                          "%rr = arith.mulf %rn, %rn : f64", "%rrw = arith.mulf %rr, %w : f64",
                          "%hr = arith.mulf %rh_t, %rn : f64", "%hrw = arith.mulf %hr, %w : f64"],
                         ["%xn", "%rn", "%rrw", "%hrw"], counting=True, reductions=2, scalars=("alpha", "omega")))
    weigh = ["%one = arith.constant 1.0 : f64", "%zw = arith.constant 0.0 : f64",
             "%w = arith.select %counted, %one, %zw : f64"]
    out.append(_cellwise(f"rbicg_update_p_p{p}", p, ["r", "v", "p"], ["p", "pw"],
                         ["%ov = arith.mulf %omega, %v_t : f64", "%d = arith.subf %p_t, %ov : f64",
                          "%bd = arith.mulf %beta, %d : f64", "%pn = arith.addf %r_t, %bd : f64"] + weigh +
                         ["%pwn = arith.mulf %pn, %w : f64"],
                         ["%pn", "%pwn"], counting=True, scalars=("beta", "omega")))
    out.append(_cellwise(f"rbicg_update_s_p{p}", p, ["r", "v"], ["s", "sw"],
                         ["%av = arith.mulf %alpha, %v_t : f64", "%sn = arith.subf %r_t, %av : f64"] + weigh +
                         ["%swn = arith.mulf %sn, %w : f64"],
                         ["%sn", "%swn"], counting=True, scalars=("alpha",)))
    out.append(_cellwise(f"rbicg_update_xr_p{p}", p, ["x", "ph", "sh", "s", "t", "rh"], ["x", "r"],
                         ["%ap = arith.mulf %alpha, %ph_t : f64", "%os = arith.mulf %omega, %sh_t : f64",
                          "%inc = arith.addf %ap, %os : f64", "%xn = arith.addf %x_t, %inc : f64",
                          "%ot = arith.mulf %omega, %t_t : f64", "%rn = arith.subf %s_t, %ot : f64"] + weigh +
                         ["%rr = arith.mulf %rn, %rn : f64", "%rrw = arith.mulf %rr, %w : f64",
                          "%hr = arith.mulf %rh_t, %rn : f64", "%hrw = arith.mulf %hr, %w : f64"],
                         ["%xn", "%rn", "%rrw", "%hrw"], counting=True, reductions=2, scalars=("alpha", "omega")))
    for nv in (1, 2):
        out.append("\n".join(
            [f"func.func @fold_nv{nv}(%src: memref<?xf64>, %dst: memref<?xf64>, %count: index) "
             "attributes {mir.kernel, mir.noalias} {"] + _consts(nv) +
            ["  %dn = memref.dim %dst, %c0 : memref<?xf64>", "  %blocks = arith.divui %dn, %nv : index",
             "  %bm1 = arith.subi %blocks, %c1 : index", "  %cb = arith.addi %count, %bm1 : index",
             "  %per = arith.divui %cb, %blocks : index",
             "  scf.parallel (%b) = (%c0) to (%blocks) step (%c1) {",
             "    %lo = arith.muli %b, %per : index", "    %lp = arith.addi %lo, %per : index",
             "    %hi = arith.minui %lp, %count : index"] + _sum_range(nv, "%lo", "%hi") +
            ["    %o = arith.muli %b, %nv : index"] +
            sum(([f"    %o{j} = arith.addi %o, %c{j} : index", f"    memref.store {s_}, %dst[%o{j}] : memref<?xf64>"]
                 for j, s_ in enumerate(_sums(nv))), []) +
            ["  } {mir.elements}", "  return", "}", ""]))
        out.append("\n".join(
            [f"func.func @totals_nv{nv}(%src: memref<?xf64>, %sc: memref<?xf64>, %count: index) "
             "attributes {mir.kernel, mir.noalias} {"] + _consts(nv) +
            ["  scf.parallel (%b) = (%c0) to (%c1) step (%c1) {"] + _sum_range(nv, "%c0", "%count") +
            sum(([f"    %ti{j} = arith.constant {KRYLOV['t0'] + j} : index",
                  f"    memref.store {s_}, %sc[%ti{j}] : memref<?xf64>"] for j, s_ in enumerate(_sums(nv))), []) +
            ["  } {mir.elements}", "  return", "}", ""]))
    for step, nv in STEP_NV.items():
        out.append("\n".join(
            [f"func.func @finish_{step}(%src: memref<?xf64>, %sc: memref<?xf64>, %count: index) "
             "attributes {mir.kernel, mir.noalias} {"] + _consts(nv) +
            ["  scf.parallel (%b) = (%c0) to (%c1) step (%c1) {"] + _sum_range(nv, "%c0", "%count") +
            _step_lines(step, _sums(nv)) + ["  } {mir.elements}", "  return", "}", ""]))
        out.append("\n".join(
            [f"func.func @step_{step}(%sc: memref<?xf64>) attributes {{mir.kernel}} {{",
             "  %c0 = arith.constant 0 : index", "  %c1 = arith.constant 1 : index",
             "  scf.parallel (%b) = (%c0) to (%c1) step (%c1) {"] +
            sum(([f"    %tl{j} = arith.constant {KRYLOV['t0'] + j} : index",
                  f"    %tv{j} = memref.load %sc[%tl{j}] : memref<?xf64>"] for j in range(nv)), []) +
            _step_lines(step, [f"%tv{j}" for j in range(nv)]) + ["  } {mir.elements}", "  return", "}", ""]))
    return "".join(out)


# ---- p-multigrid V-cycle (PMultigrid of mars_cellwise_pmg.hpp) ------------------
_TABLES = [("face_nbr", "i32"), ("face_code", "i8"), ("edge_of", "i32"), ("edge_off", "i32"), ("edge_ent", "i32"),
           ("vert_of", "i32"), ("vert_off", "i32"), ("vert_ent", "i32")]
_FLAGS = [("counting", "i32"), ("edge_dirichlet", "i8"), ("vert_dirichlet", "i8")]


def _chebyshev(p, kind, real="f64"):
    """One Chebyshev-Jacobi step at degree p (unstructured_chebyshev_kernel): z = 0 on
    Dirichlet copies, DSS(b - y) / diag elsewhere; first0: d = x = c2 z (x was 0, no
    y); first: d = c2 z, x += d; next: d = c1 d + c2 z, x += d. One argument list
    for the three: (b, ghost, y, elements, tables, flags, diag, d, x, c1, c2)."""
    args = [("b", real), ("ghost", real), ("y", real), ("elements", "i32")] + _TABLES + _FLAGS + \
        [("diag", real), ("d", real), ("x", real)]
    sig = ", ".join(f"%{n}: memref<?x{t}>" for n, t in args) + f", %c1: {real}, %c2: {real}"
    copies = {"first0": ["diag"], "first": ["diag", "x"], "next": ["diag", "d", "x"]}[kind]
    minus = "" if kind == "first0" else " minus(%y)"
    mr = f"memref<?x{real}>"
    opnd_types = [mr, mr] + ([] if kind == "first0" else [mr]) + \
        ["memref<?xi32>"] + [f"memref<?x{t}>" for _, t in _TABLES + _FLAGS] + [mr] * (len(copies) + 2)
    body = [f"  %zero = arith.constant 0.0 : {real}",
            f"  %q = arith.divf %sum, %diag_t : {real}",
            f"  %z = arith.select %dir, %zero, %q : {real}",
            f"  %cz = arith.mulf %c2, %z : {real}"]
    if kind == "first0":
        ys = ["%cz", "%cz"]
    elif kind == "first":
        body += [f"  %xn = arith.addf %x_t, %cz : {real}"]
        ys = ["%cz", "%xn"]
    else:
        body += [f"  %cd = arith.mulf %c1, %d_t : {real}", f"  %dn = arith.addf %cd, %cz : {real}",
                 f"  %xn = arith.addf %x_t, %dn : {real}"]
        ys = ["%dn", "%xn"]
    bargs = f"%sum: {real}, %dir: i1, %cnt: i1" + "".join(f", %{c}_t: {real}" for c in copies)
    name = f"cheb_{kind}_p{p}" + ("" if real == "f64" else "_" + real)
    return "\n".join([
        f"func.func @{name}({sig}) attributes {{mir.kernel, mir.noalias}} {{",
        f"  mir.dss ins(%b, %ghost){minus} elements(%elements) faces(%face_nbr, %face_code)",
        "          edges(%edge_of, %edge_off, %edge_ent) vertices(%vert_of, %vert_off, %vert_ent)",
        "          flags(%counting, %edge_dirichlet, %vert_dirichlet)",
        f"          copies({', '.join('%' + c for c in copies)}) outs(%d, %x)",
        f"          {{degree = {p} : i64, reductions = 0 : i64}} : {', '.join(opnd_types)} epilogue {{",
        f"  ^bb0({bargs}):", *["  " + l for l in body],
        f"    mir.yield {', '.join(ys)} : {real}, {real}", "  }", "  return", "}", ""])


def _masked_copy(p, real="f64"):
    """r = b, zero on Dirichlet copies (the V-cycle's residual with pre = 0, then
    zero_dirichlet)."""
    dt = [("face_code", "i8"), ("edge_of", "i32"), ("vert_of", "i32"), ("edge_dirichlet", "i8"),
          ("vert_dirichlet", "i8")]
    args = [("elements", "i32")] + dt + [("b", real), ("r", real)]
    sig = ", ".join(f"%{n}: memref<?x{t}>" for n, t in args)
    types = ", ".join(f"memref<?x{t}>" for _, t in args)
    name = f"masked_copy_p{p}" + ("" if real == "f64" else "_" + real)
    return "\n".join([
        f"func.func @{name}({sig}) attributes {{mir.kernel, mir.noalias}} {{",
        "  mir.cellwise elements(%elements) dirichlet(%face_code, %edge_of, %vert_of, %edge_dirichlet, %vert_dirichlet)"
        f" inputs(%b) outs(%r) {{degree = {p} : i64}} : {types} body {{",
        f"  ^bb0(%cnt: i1, %dir: i1, %b_t: {real}):",
        f"    %zero = arith.constant 0.0 : {real}",
        f"    %rv = arith.select %dir, %zero, %b_t : {real}",
        f"    mir.yield %rv : {real}", "  }", "  return", "}", ""])


def _interpolate(rows, cols, add, real="f64", counted=False):
    name = f"interp_{rows}x{cols}" + ("_add" if add else "") + ("_counted" if counted else "") + \
        ("" if real == "f64" else "_" + real)
    mr = f"memref<?x{real}>"
    cnt_arg, cnt_op, cnt_type = ((", %counting: memref<?xi32>", " counting(%counting)", ", memref<?xi32>")
                                 if counted else ("", "", ""))
    return "\n".join([
        f"func.func @{name}(%elements: memref<?xi32>{cnt_arg}, %matrix: {mr}, %in: {mr}, "
        f"%out: {mr}) attributes {{mir.kernel, mir.noalias}} {{",
        f"  mir.interpolate elements(%elements){cnt_op} matrix(%matrix) ins(%in) outs(%out) {{rows = {rows} : i64, "
        f"cols = {cols} : i64{', accumulate' if add else ''}}} : memref<?xi32>{cnt_type}, {mr}, {mr}, {mr}",
        "  return", "}", ""])


def _convert(p, src, dst):
    """dst = src in another precision, on the listed elements' copies (the V-cycle's
    input and output when it runs in single precision)."""
    op = "arith.truncf" if src == "f64" else "arith.extf"
    return "\n".join([
        f"func.func @to_{dst}_p{p}(%elements: memref<?xi32>, %src: memref<?x{src}>, %dst: memref<?x{dst}>) "
        "attributes {mir.kernel, mir.noalias} {",
        f"  mir.cellwise elements(%elements) inputs(%src) outs(%dst) {{degree = {p} : i64}} "
        f": memref<?xi32>, memref<?x{src}>, memref<?x{dst}> body {{",
        f"  ^bb0(%cnt: i1, %v: {src}):",
        f"    %w = {op} %v : {src} to {dst}",
        f"    mir.yield %w : {dst}", "  }", "  return", "}", ""])


def _dots(p, az, zz):
    """Per-element partials of <a, z> / <z, z> after the V-cycle (wdots_kernel)."""
    body, ys = ["%one = arith.constant 1.0 : f64", "%zw = arith.constant 0.0 : f64",
                "%w = arith.select %counted, %one, %zw : f64"], []
    if az:
        body += ["%av = arith.mulf %a_t, %z_t : f64", "%avw = arith.mulf %av, %w : f64"]
        ys.append("%avw")
    if zz:
        body += ["%zz2 = arith.mulf %z_t, %z_t : f64", "%zzw = arith.mulf %zz2, %w : f64"]
        ys.append("%zzw")
    name = f"pmg_dots_p{p}" + ("_az" if az else "") + ("_zz" if zz else "")
    return _cellwise(name, p, ["z"] + (["a"] if az else []), [], body, ys, counting=True, reductions=len(ys))


def emit_pmg(real="f64"):
    """The per-iteration kernels of the p-multigrid V-cycle 7 -> 3 -> 1 besides the
    operators (the operator pipeline generates those at p = 7, 3, 1). real = f32: the
    same kernels in single precision (names end in _f32), the dot products left out
    (they belong to the FP64 Krylov solver) and added: the conversions to_f32_p* /
    to_f64_p* (at p = 7 the V-cycle's input and output, at p = 1 those of an FP64
    coarse solve, BoomerAMG), and the continuous prolongation 3 -> 7: interp_8x4_counted
    writes the counting copies only, dss_p7 spreads them. In FP32 the copies of a node
    that the plain prolongation gives differ by rounding (each element sums in its own
    frame's order), and the right-preconditioned solver cannot see such jumps in its
    residual: they stay in x as an error near 1e-8."""
    out = [_chebyshev(p, k, real) for p in (7, 3, 1) for k in ("first0", "first", "next")]
    out += [_masked_copy(p, real) for p in (7, 3)]
    out += [_interpolate(4, 8, False, real), _interpolate(2, 4, False, real)]          # restrict 7->3, 3->1
    out += [_interpolate(8, 4, add, real) for add in (False, True)]                     # prolong 3->7
    out += [_interpolate(4, 2, add, real) for add in (False, True)]                     # prolong 1->3
    if real == "f64":
        out += [_dots(7, az, zz) for az, zz in ((True, False), (False, True), (True, True))]
    else:
        out += [_convert(p, a, b) for p in (7, 1) for a, b in (("f64", real), (real, "f64"))]
        out += [_interpolate(8, 4, False, real, counted=True), emit_dss(7, real)]
    return "".join(out)
