#!/usr/bin/env python3
# Lane-owned scratch in --mir-warp-per-element=lane-owned: a buffer every phase touches only at
# its own points (same iteration space, the phase's induction variables in the same positions,
# constants for dimensions kept in the thread) becomes a private array; one that a phase reads at
# another lane's point stays in shared memory. Numerical checks: bp_gpu_pipeline.py --lane-owned
# --emulate (with --emu-tsan) on a z-march mir.
# Run from marsir-mlir/: python3 test/lane_owned_gate.py

import os
import subprocess
import sys

MIROPT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "build", "tools", "mir-opt", "mir-opt")

# w[a][k]: phase 1 writes w[a][k] for its point a (k kept in the thread), phase 2 reads w[{READ}][k].
KERNEL = """func.func @k(%x: memref<?x4xf64>, %y: memref<?x4xf64>) {{
  %c0 = arith.constant 0 : index
  %c1 = arith.constant 1 : index
  %c3 = arith.constant 3 : index
  %c4 = arith.constant 4 : index
  %n = memref.dim %x, %c0 : memref<?x4xf64>
  scf.parallel (%e) = (%c0) to (%n) step (%c1) {{
    %w = memref.alloca() {{mir.scratch}} : memref<4x2xf64>
    scf.parallel (%a) = (%c0) to (%c4) step (%c1) {{
      %v = memref.load %x[%e, %a] : memref<?x4xf64>
      memref.store %v, %w[%a, %c0] : memref<4x2xf64>
      memref.store %v, %w[%a, %c1] : memref<4x2xf64>
      scf.reduce
    }} {{mir.lanes}}
    scf.parallel (%a) = (%c0) to (%c4) step (%c1) {{
      %r = arith.subi %c3, %a : index
      %u = memref.load %w[{READ}, %c1] : memref<4x2xf64>
      memref.store %u, %y[%e, %a] : memref<?x4xf64>
      scf.reduce
    }} {{mir.lanes}}
    scf.reduce
  }} {{mir.elements, mir.lanes_per_element = 4 : i64, mir.warps_per_block = 1 : i64}}
  return
}}
"""


def run(read):
    r = subprocess.run([MIROPT, "-", "--mir-warp-per-element=lane-owned=true"], input=KERNEL.format(READ=read),
                       capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


ok = True
rc, ir = run("%a")   # own point: lane-owned
private = rc == 0 and "mir.workgroup" not in ir and "memref<1x2xf64>" in ir
print(f"  own points -> private array, no shared buffer: {'ok' if private else 'NO'}")
ok &= private
rc, ir = run("%r")   # another lane's point: must stay shared
shared = rc == 0 and "mir.workgroup" in ir and "gpu.barrier" in ir
print(f"  another lane's point -> shared, with a barrier: {'ok' if shared else 'NO'}")
ok &= shared
print("LANE-OWNED GATE:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
