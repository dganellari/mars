#!/usr/bin/env python3
# Which tables --mir-even-odd may fold. Folding changes the table layout, so a table with
# any access the pass does not rewrite must stay as it is. The rejected cases must come
# out identical to the input; an alias read of a folded table once returned 4.5 instead
# of 4 (review internal-notes/bp_even_odd_review_2026-11-08.md). Numerical correctness
# of folded sweeps is the BP emulator gate (bp_gpu_pipeline.py --emulate).
# Run from marsir-mlir/: python3 test/even_odd_gate.py

import os
import subprocess
import sys

MIROPT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "build", "tools", "mir-opt", "mir-opt")

# y0 = M00 x0 + M01 x1, y1 = M10 x0 + M11 x1, returns y0 + y1 + EXTRA
SWEEP = """module {{
 memref.global @M : memref<2x2xf64> = uninitialized {{mir.reflect = 1 : i64}}
 func.func @test(%x0: f64, %x1: f64) -> f64 {{
  %g = memref.get_global @M : memref<2x2xf64>
  %c0 = arith.constant 0 : index
  %c1 = arith.constant 1 : index
  %zero = arith.constant 0.0 : f64
  %a = memref.load %g[%c0, %c0] : memref<2x2xf64>
  %b = memref.load %g[%c0, %c1] : memref<2x2xf64>
  %c = memref.load %g[%c1, %c0] : memref<2x2xf64>
  %d = memref.load %g[%c1, %c1] : memref<2x2xf64>
  %y00 = math.fma %a, %x0, %zero : f64
  %y0 = math.fma %b, %x1, %y00 : f64
  %y10 = math.fma %c, %x0, %zero : f64
  %y1 = math.fma %d, %x1, %y10 : f64
  %sum = arith.addf %y0, %y1 : f64
{extra}
 }}
}}
"""

CASES = {
    "complete sweep": (True, "  return %sum : f64"),
    "alias read": (False, """  %alias = memref.reinterpret_cast %g to offset: [0], sizes: [2, 2], strides: [2, 1] : memref<2x2xf64> to memref<2x2xf64, strided<[2, 1]>>
  %e = memref.load %alias[%c0, %c0] : memref<2x2xf64, strided<[2, 1]>>
  %r = arith.addf %sum, %e : f64
  return %r : f64"""),
    "direct read outside a sweep": (False, """  %e = memref.load %g[%c0, %c0] : memref<2x2xf64>
  %r = arith.addf %sum, %e : f64
  return %r : f64"""),
    "unused load": (False, """  %e = memref.load %g[%c1, %c0] : memref<2x2xf64>
  return %sum : f64"""),
}


def opt(ir, *flags):
    p = subprocess.run([MIROPT, "-", *flags], input=ir, capture_output=True, text=True)
    if p.returncode:
        sys.exit(f"mir-opt failed:\n{p.stderr}")
    return p.stdout


failed = 0
for name, (folds, extra) in CASES.items():
    ir = SWEEP.format(extra=extra)
    out = opt(ir, "--mir-even-odd")
    if folds:
        ok = "@M__fold" in out and "memref<2xf64>" in out and "memref<2x2xf64>" not in out
    else:
        ok = out == opt(ir)
    print(f"  {name:28} {'folded' if '@M__fold' in out else 'kept':6}  {'PASS' if ok else 'FAIL'}")
    failed += not ok
sys.exit(1 if failed else 0)
