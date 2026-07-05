#!/usr/bin/env python3
"""Drop the (already-applied) transform-dialect schedule module from stdin IR,
so downstream tools need not parse transform ops."""
import sys
lines = sys.stdin.read().splitlines(True)
out, depth, skipping = [], 0, False
for ln in lines:
    if not skipping and ln.lstrip().startswith("module attributes {transform.with_named_sequence}"):
        skipping = True
        depth = ln.count("{") - ln.count("}")
        continue
    if skipping:
        depth += ln.count("{") - ln.count("}")
        if depth <= 0:
            skipping = False
        continue
    out.append(ln)
sys.stdout.write("".join(out))
