"""MLIR backend: emit the operator as a mir-dialect partial-assembly chain.

Stage 2 twin of the source emitters. The SAME parsed operator spec that feeds
host_cpp / cuda here becomes textual MLIR in the `mir` dialect -- the
B^T . D . B sum-factorization skeleton:

  B   : mir.contract   -- sweep the field to the quadrature grid (per axis)
  D   : mir.flux        -- the authored pointwise flux, walked from the AST
  B^T : mir.contract   -- integrate back to modes

Only the flux region differs between operators; the contractions are the fixed
skeleton. The emitted text round-trips through the real dialect parser
(mirage-mlir/.../mir-opt), so it is the on-ramp to the mir -> linalg -> gpu
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
