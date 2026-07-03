"""mirage IR: the operator spec and the synthesized element-apply schedule.

The front-end produces an `Operator` (the math the user authored). The
synthesizer expands it into an `ElementApply` -- the fixed gather -> contract
-> flux -> integrate -> scatter Knaus Alg-2 skeleton -- recording only what the
authored flux actually needs (which tangential derivatives, whether the metric
has cross terms). Backends walk `ElementApply` and emit code; they never see the
raw spec text.

This mirrors the libCEED split: the user writes the pointwise action (our
`flux`), the library owns the E/B/D gather-contract-scatter structure.
"""

from .expr import Var, Num, Binary, simplify, _is_zero

# Names the flux expression may reference. These are the quantities the Alg-2
# skeleton makes available at each subcontrol-surface quadrature slot (s,r):
#   deriv : normal derivative of u at the SCS face  (Dtil . u)
#   dt1   : tangential derivative along the s-axis   (D_s . interp)
#   dt2   : tangential derivative along the r-axis   (D_r . interp)
#   g0,g1,g2 : metric components (tang2 r-axis, tang1 s-axis, normal)
FLUX_INPUTS = {"deriv", "dt1", "dt2", "g0", "g1", "g2"}

SUPPORTED_ELEMENTS = {"hex"}
SUPPORTED_FORMS = {"diffusion"}


class Operator:
    """What the user authored in the .op spec."""
    def __init__(self, name, element, form, flux, flux_src):
        self.name = name
        self.element = element
        self.form = form
        self.flux = flux
        self.flux_src = flux_src  # original text, for the IR dump / provenance


class ElementApply:
    """Synthesized element-local schedule for one operator.

    The `diffusion` form is a per-direction (dir=0,1,2) sweep over the p
    subcontrol-surface faces l; at each face the 4-step contraction runs over
    the (p+1)^2 tangential slots. Backends emit that fixed skeleton and inline
    `flux` at step 2. The flags below let a backend drop work the authored flux
    does not use (e.g. an orthogonal-cube operator with no cross terms skips the
    two tangential D-contractions entirely).
    """
    def __init__(self, op, needs_deriv, needs_dt1, needs_dt2, needs_metric,
                 free_vars=None, block=None):
        self.op = op
        self.needs_deriv = needs_deriv
        self.needs_dt1 = needs_dt1
        self.needs_dt2 = needs_dt2
        self.needs_metric = needs_metric
        self.free_vars = free_vars if free_vars is not None else set()
        self.block = block          # the typed batch this operator acts on

    @property
    def needs_tangential(self):
        return self.needs_dt1 or self.needs_dt2


class Block:
    """A typed batch of alike elements -- Mirage's analogue of NektarIR's
    `!nir.block`. It formalizes AS A TYPE the metadata that was otherwise
    scattered as loose flags (which element shape, basis, precision, batch), and
    is the connective tissue toward a real MLIR dialect (Stage 2).

    The concrete static fields (fields, shape, basis, scheme) come from the
    authored operator. The remaining axes -- Order (P), Deformed (curved vs
    affine = our PerPoint vs Affine metric mode), Precision, Batch (E) -- are the
    codegen configuration, carried as TEMPLATE PARAMETERS in the emitted kernel,
    so the one generated header covers all of them.
    """
    def __init__(self, fields, shape, basis, scheme):
        self.fields = fields        # e.g. ["u"]
        self.shape = shape          # "hex" | "tet"
        self.basis = basis          # "gll" (GLL-Lagrange nodal)
        self.scheme = scheme        # "knaus-scs" (GLL solution nodes + Gauss SCS pts)

    def render(self):
        f = ", ".join(self.fields)
        return (f"!mir.block<Fields: [{f}], Shape: {self.shape}, Basis: {self.basis}, "
                f"Scheme: {self.scheme}, Order: P, Deformed: {{affine|curved}}, "
                f"Precision: (sol=f64, metric=f64), Batch: E=MirageLaunchDefault<P>>")


def block_for(op):
    # The operator reads a single solution field through the flux; the batch is a
    # collection of `op.element` elements with the HO-CVFEM (Knaus) GLL scheme.
    return Block(fields=["u"], shape=op.element, basis="gll", scheme="knaus-scs")


def synthesize(op: Operator) -> ElementApply:
    if op.element not in SUPPORTED_ELEMENTS:
        raise ValueError(f"element {op.element!r} not supported "
                         f"(have {sorted(SUPPORTED_ELEMENTS)})")
    if op.form not in SUPPORTED_FORMS:
        raise ValueError(f"form {op.form!r} not supported "
                         f"(have {sorted(SUPPORTED_FORMS)})")

    fv = op.flux.free_vars()
    unknown = fv - FLUX_INPUTS
    if unknown:
        raise ValueError(
            f"flux references unknown name(s) {sorted(unknown)}; "
            f"available at a quadrature slot: {sorted(FLUX_INPUTS)}")

    return ElementApply(
        op=op,
        needs_deriv="deriv" in fv,
        needs_dt1="dt1" in fv,
        needs_dt2="dt2" in fv,
        needs_metric=bool(fv & {"g0", "g1", "g2"}),
        free_vars=fv,
        block=block_for(op),
    )


# The flux inputs that are LINEAR functionals of the solution u: deriv = Dtil.u,
# dt1/dt2 = D.(Btil.u). The metric g0,g1,g2 is geometry -- constant w.r.t. u. So
# the operator's u-dependence enters ONLY through these, and the Jacobian-action
# (JVP) in a direction du is, by the chain rule,
#     J*du flux = sum_x (d flux / d x)|_u  *  x_dot
# where x_dot is the SAME contraction applied to du (deriv_dot, dt1_dot, dt2_dot).
LINEAR_INPUTS = ("deriv", "dt1", "dt2")


def jacobian_action_flux(op: Operator):
    """Symbolic form differentiation: build the JVP flux + the per-input partials.

    Returns (jvp_flux, partials) where partials[x] = simplify(d flux / d x) and
    jvp_flux = sum over x with nonzero partial of  partial * x_dot. For a flux
    that is LINEAR in u (e.g. the Laplacian) the partials are u-independent and
    jvp_flux reduces to the primal flux with x -> x_dot, so J*du == A*du. For a
    nonlinear flux the partials reference the primal state, so the emitted JVP
    kernel gathers both u (for the partials) and du (for the x_dot).
    """
    partials = {}
    terms = []
    for x in LINEAR_INPUTS:
        p = simplify(op.flux.diff(x))
        partials[x] = p
        if _is_zero(p):
            continue
        terms.append(Binary("*", p, Var(x + "_dot")))
    if not terms:
        jvp = Num("0")
    else:
        jvp = terms[0]
        for t in terms[1:]:
            jvp = Binary("+", jvp, t)
    return simplify(jvp), partials


def dump(ea: ElementApply) -> str:
    """Human-readable IR dump (Stage-1 deliverable: parse -> IR -> text)."""
    o = ea.op
    lines = [
        f"operator {o.name}",
        f"  element : {o.element}",
        f"  form    : {o.form}  (Knaus Alg-2 sum-factorized diffusion apply)",
        f"  flux    : {o.flux_src}",
        f"  uses    : {', '.join(sorted(ea.free_vars)) or '(none)'}",
        f"  block   : {ea.block.render()}",
        "  schedule (per element):",
        "    gather   u_local <- u_global[edof]        ; y_local <- 0",
        "    for dir in 0,1,2:",
        "      for l in 0..p-1:                          # subcontrol-surface faces",
        "        step1 contract-normal   : interp = Btil.u ; deriv = Dtil.u",
    ]
    if ea.needs_tangential:
        td = []
        if ea.needs_dt2:
            td.append("dt2 = D_r.interp")
        if ea.needs_dt1:
            td.append("dt1 = D_s.interp")
        lines.append(f"        step2 contract-tangential: {' ; '.join(td)}")
    else:
        lines.append("        step2 contract-tangential: (skipped -- flux has no cross terms)")
    lines += [
        f"        step2 flux              : flux = {o.flux_src}",
        "        step3 integrate         : tmp = W_r.flux ; intf = W_s.tmp",
        "        step4 scatter-face      : y[l] -= intf ; y[l+1] += intf",
        "    scatter  y_global[edof] += y_local",
    ]
    # Symbolic form differentiation: the Jacobian-action operator, derived from
    # the same flux (not hand-written). Self-adjoint/linear flux -> J*du flux
    # equals the primal flux with x -> x_dot.
    jvp, partials = jacobian_action_flux(o)
    lines.append("  jacobian (symbolic form differentiation):")
    for x in LINEAR_INPUTS:
        lines.append(f"    d flux / d {x:<5} = {partials[x].render({})}")
    lines.append(f"    J*du flux           = {jvp.render({})}")
    return "\n".join(lines)
