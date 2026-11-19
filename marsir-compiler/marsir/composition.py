"""Add scalar, stored-geometry Galerkin forms before element scheduling."""
import re

from .expr import Binary, parse_expr
from .ir import GALERKIN_INPUTS, GALERKIN_OUTPUTS, TRILINEAR_INPUTS, Operator, synthesize

_IDENTIFIER = r"[A-Za-z_][A-Za-z0-9_]*"
_TERM = re.compile(r"(" + _IDENTIFIER + r")(?:\s*\(([^()]*)\))?")
_BINDING = re.compile(r"(g[0-5])\s*=\s*(" + _IDENTIFIER + r")")


def coefficient_name(name):
    return (re.fullmatch(_IDENTIFIER, name) is not None
            and not re.fullmatch(r"g[0-9]+", name)
            and name not in GALERKIN_INPUTS | TRILINEAR_INPUTS | set(GALERKIN_OUTPUTS))


def _bindings(text, form):
    aliases = {}
    if text is None:
        return aliases
    for item in text.split(","):
        match = _BINDING.fullmatch(item.strip())
        if not match:
            raise ValueError(f"expected a stored slot binding such as g0=sigma_w, got {item!r}")
        slot, field = match.groups()
        if slot not in form.free_vars or slot in aliases:
            raise ValueError(f"unused or repeated slot {slot!r} in {form.op.name}")
        if not coefficient_name(field):
            raise ValueError(f"invalid or reserved quadrature field {field!r}")
        aliases[slot] = field
    return aliases


def compose(name, expression, resolve):
    terms = expression.split("+")
    if len(terms) < 2:
        raise ValueError("a sum needs at least two operator references")
    outputs = {}
    fields = []
    used_slots = set()
    for term in terms:
        match = _TERM.fullmatch(term.strip())
        if not match:
            raise ValueError(f"invalid operator reference {term.strip()!r}")
        source, bindings = match.groups()
        op = resolve(source)
        if op.element != "hex" or op.form != "galerkin" or op.geometry != "stored":
            raise ValueError(f"{source}: sums require stored-geometry Galerkin hex operators")
        # Validate each operand before merging names. A field from the next
        # operand must never disguise an undeclared input in this one.
        form = synthesize(op)
        aliases = _bindings(bindings, form)
        slots = {v for v in form.free_vars if re.fullmatch(r"g[0-5]", v)} - aliases.keys()
        overlap = sorted(used_slots & slots)
        if overlap:
            raise ValueError(f"ambiguous stored slots {overlap} in {source}; bind to named fields")
        used_slots.update(slots)
        for field in list(op.coefficients) + list(aliases.values()):
            if field not in fields:
                fields.append(field)
        for output, expr in op.outputs:
            if aliases:
                expr = parse_expr(expr.render(aliases))
            outputs[output] = Binary("+", outputs[output], expr) if output in outputs else expr

    ordered = [(k, outputs[k]) for k in GALERKIN_OUTPUTS if k in outputs]
    return Operator(name, "hex", "galerkin", None, expression,
                    outputs=ordered, coefficients=tuple(fields))
