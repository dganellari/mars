"""Front-end: parse an operator spec (.op) into an `Operator`.

Grammar (deliberately tiny, line-oriented):

    operator <name> {
        element = hex
        form    = diffusion
        flux    = g2*deriv + g0*dt2 + g1*dt1
    }

'#' starts a comment. Only one operator per file. `flux` is an arithmetic
expression over the quadrature-slot inputs (see mir.ir.FLUX_INPUTS).
"""

import re

from .expr import parse_expr
from .ir import Operator

_HEADER = re.compile(r"^\s*operator\s+([A-Za-z_]\w*)\s*\{\s*$")
_KV = re.compile(r"^\s*([A-Za-z_]\w*)\s*=\s*(.+?)\s*$")


def _strip_comment(line: str) -> str:
    # No strings in the grammar, so a bare '#' always starts a comment.
    h = line.find("#")
    return line if h < 0 else line[:h]


def parse_spec(text: str) -> Operator:
    name = None
    fields = {}
    in_body = False
    closed = False

    for raw in text.splitlines():
        line = _strip_comment(raw)
        if not line.strip():
            continue
        if not in_body:
            m = _HEADER.match(line)
            if not m:
                raise ValueError(f"expected 'operator <name> {{', got: {raw!r}")
            name = m.group(1)
            in_body = True
            continue
        if line.strip() == "}":
            closed = True
            in_body = False
            continue
        m = _KV.match(line)
        if not m:
            raise ValueError(f"expected 'key = value' or '}}', got: {raw!r}")
        key, val = m.group(1), m.group(2)
        if key in fields:
            raise ValueError(f"duplicate key {key!r}")
        fields[key] = val

    if name is None:
        raise ValueError("no operator block found")
    if not closed:
        raise ValueError("operator block not closed with '}'")

    for req in ("element", "form", "flux"):
        if req not in fields:
            raise ValueError(f"missing required key {req!r}")

    return Operator(
        name=name,
        element=fields["element"],
        form=fields["form"],
        flux=parse_expr(fields["flux"]),
        flux_src=fields["flux"],
    )


def parse_spec_file(path: str) -> Operator:
    with open(path) as f:
        return parse_spec(f.read())
