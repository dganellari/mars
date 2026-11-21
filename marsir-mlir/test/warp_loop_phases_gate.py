#!/usr/bin/env python3
# Loops of phases in --mir-warp-per-element (the planes of a z-march schedule): phases in an
# scf.for of the element body get a barrier between them inside the loop and one at the end of
# each trip (the next trip rewrites the plane another lane may still read); a loop whose trip
# count is not a constant is rejected, since every lane must reach the same barriers.
# Run from marsir-mlir/: python3 test/warp_loop_phases_gate.py

import os
import re
import subprocess
import sys

MIROPT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "build", "tools", "mir-opt", "mir-opt")

KERNEL = """func.func @k(%x: memref<?x4x4xf64>, %y: memref<?x4x4xf64>) {{
  %c0 = arith.constant 0 : index
  %c1 = arith.constant 1 : index
  %c4 = arith.constant 4 : index
  %n = memref.dim %x, %c0 : memref<?x4x4xf64>
  scf.parallel (%e) = (%c0) to (%n) step (%c1) {{
    %plane = memref.alloca() {{mir.scratch}} : memref<4xf64>
    scf.for %c = %c0 to {UB} step %c1 {{
      scf.parallel (%a) = (%c0) to (%c4) step (%c1) {{
        %v = memref.load %x[%e, %a, %c] : memref<?x4x4xf64>
        memref.store %v, %plane[%a] : memref<4xf64>
        scf.reduce
      }} {{mir.lanes}}
      scf.parallel (%a) = (%c0) to (%c4) step (%c1) {{
        %b = arith.subi %c4, %c1 : index
        %i = arith.subi %b, %a : index
        %v = memref.load %plane[%i] : memref<4xf64>
        memref.store %v, %y[%e, %a, %c] : memref<?x4x4xf64>
        scf.reduce
      }} {{mir.lanes}}
    }}
    scf.reduce
  }} {{mir.elements, mir.lanes_per_element = 4 : i64, mir.warps_per_block = 1 : i64}}
  return
}}
"""


def run(ub):
    r = subprocess.run([MIROPT, "-", "--mir-warp-per-element"], input=KERNEL.format(UB=ub), capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


def plane_loop(ir):
    """Lines of the plane loop's body: from `scf.for ... to %c4 step %c1 {` to its closing brace."""
    lines = ir.splitlines()
    start = next(i for i, l in enumerate(lines) if re.search(r"scf\.for %\w+ = %c0 to %c4 step %c1 \{", l))
    depth, out = 0, []
    for l in lines[start:]:
        depth += l.count("{") - l.count("}")
        out.append(l)
        if depth == 0:
            return out
    return out


ok = True
rc, ir = run("%c4")
body = plane_loop(ir) if rc == 0 else []
barriers = [i for i, l in enumerate(body) if "gpu.barrier" in l]
stores = [i for i, l in enumerate(body) if "memref.store" in l and "reinterpret_cast" in l]
loads = [i for i, l in enumerate(body) if "memref.load %reinterpret_cast" in l]
between = rc == 0 and stores and loads and any(stores[0] < b < loads[0] for b in barriers)
at_end = rc == 0 and barriers and barriers[-1] > loads[-1] and body[barriers[-1] + 1].strip() == "}"
print(f"  barrier between the phases inside the loop: {'ok' if between else 'MISSING'}")
print(f"  barrier at the end of each trip:            {'ok' if at_end else 'MISSING'}")
ok &= bool(between and at_end)
rc, out = run("%n")
rejected = rc != 0 and "needs constant bounds" in out
print(f"  non-constant trip count rejected:           {'ok' if rejected else 'NOT REJECTED'}")
ok &= rejected
print("WARP LOOP-OF-PHASES GATE:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
