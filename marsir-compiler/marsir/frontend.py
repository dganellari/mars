"""Front-end: parse an operator spec (.op) into an `Operator`.

Grammar (deliberately tiny, line-oriented):

    operator <name> {
        element = hex
        form    = diffusion
        flux    = g2*deriv + g0*dt2 + g1*dt1
    }

'#' starts a comment. Only one operator per file. `flux` is an arithmetic
expression over the quadrature-slot inputs (see mir.ir.FLUX_INPUTS).

`form = galerkin` (a standard Galerkin operator on Gauss points) has no `flux`;
instead it gives one or more of the pointwise outputs `v` (tested against the
basis values) and `dv0`, `dv1`, `dv2` (tested against the reference gradient), as
expressions over the inputs of mir.ir.GALERKIN_INPUTS:

    operator bp1 {
        element = hex
        form    = galerkin
        v       = g0*u
    }

`geometry = trilinear` computes the geometry inside the kernel from each element's 8
corners; the outputs then read wdetj and m00..m22 instead of g0..g5 (mir.ir.TRILINEAR_INPUTS).
"""

import re
from pathlib import Path

from .composition import coefficient_name, compose
from .expr import parse_expr
from .ir import GALERKIN_OUTPUTS, Operator

_HEADER = re.compile(r"^\s*operator\s+([A-Za-z_]\w*)\s*\{\s*$")
_KV = re.compile(r"^\s*([A-Za-z_]\w*)\s*=\s*(.+?)\s*$")


def _strip_comment(line: str) -> str:
    # No strings in the grammar, so a bare '#' always starts a comment.
    h = line.find("#")
    return line if h < 0 else line[:h]


def parse_spec(text: str, resolve=None) -> Operator:
    lines = [_strip_comment(raw).strip() for raw in text.splitlines()]
    lines = [line for line in lines if line]
    if len(lines) == 1:
        match = re.fullmatch(r"operator\s+([A-Za-z_]\w*)\s*=\s*(.+)", lines[0])
        if match:
            if resolve is None:
                raise ValueError("operator references require parse_spec_file or a resolver")
            return compose(match.group(1), match.group(2), resolve)

    name = None
    fields = {}
    in_body = False
    closed = False

    for line in lines:
        if not in_body:
            if closed:
                raise ValueError("expected only one operator per file")
            m = _HEADER.match(line)
            if not m:
                raise ValueError(f"expected 'operator <name> {{', got: {line!r}")
            name = m.group(1)
            in_body = True
            continue
        if line.strip() == "}":
            closed = True
            in_body = False
            continue
        m = _KV.match(line)
        if not m:
            raise ValueError(f"expected 'key = value' or '}}', got: {line!r}")
        key, val = m.group(1), m.group(2)
        if key in fields:
            raise ValueError(f"duplicate key {key!r}")
        fields[key] = val

    if name is None:
        raise ValueError("no operator block found")
    if not closed:
        raise ValueError("operator block not closed with '}'")

    for req in ("element", "form"):
        if req not in fields:
            raise ValueError(f"missing required key {req!r}")

    if fields["form"] == "galerkin":
        outputs = [(k, fields[k]) for k in GALERKIN_OUTPUTS if k in fields]
        if not outputs:
            raise ValueError(f"form galerkin needs at least one of {GALERKIN_OUTPUTS}")
        coefficients = ()
        if "fields" in fields:
            coefficients = tuple(name.strip() for name in fields["fields"].split(","))
            for name_field in coefficients:
                if not coefficient_name(name_field):
                    raise ValueError(f"invalid or reserved quadrature field {name_field!r}")
            if len(set(coefficients)) != len(coefficients):
                raise ValueError("duplicate quadrature field")
            if fields.get("geometry", "stored") != "stored":
                raise ValueError("named quadrature fields require stored geometry")
        extra = set(fields) - {"element", "form", "geometry", "fields"} - set(GALERKIN_OUTPUTS)
        if extra:
            raise ValueError(f"unknown key(s) {sorted(extra)} for form galerkin")
        return Operator(
            name=name,
            element=fields["element"],
            form="galerkin",
            flux=None,
            flux_src="; ".join(f"{k} = {v}" for k, v in outputs),
            coefficients=coefficients,
            outputs=[(k, parse_expr(v)) for k, v in outputs],
            geometry=fields.get("geometry", "stored"),
        )

    if "fields" in fields:
        raise ValueError("named quadrature fields require a Galerkin form")
    if "flux" not in fields:
        raise ValueError("missing required key 'flux'")
    return Operator(
        name=name,
        element=fields["element"],
        form=fields["form"],
        flux=parse_expr(fields["flux"]),
        flux_src=fields["flux"],
    )


def parse_spec_file(path: str) -> Operator:
    def load(file, parents):
        if file in parents:
            chain = " -> ".join(p.stem for p in parents + (file,))
            raise ValueError(f"cyclic operator reference: {chain}")

        def resolve(source):
            target = file.with_name(source + ".op")
            if not target.is_file():
                raise ValueError(f"{file.name}: missing operator file {target.name}")
            result = load(target, parents + (file,))
            if result.name != source:
                raise ValueError(f"{target.name} declares {result.name!r}, expected {source!r}")
            return result

        return parse_spec(file.read_text(), resolve)

    return load(Path(path), ())
