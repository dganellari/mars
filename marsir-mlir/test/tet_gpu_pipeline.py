#!/usr/bin/env python3
"""P4: the collapsed-tet operator from high-level mir to a GPU kernel, by passes
only -- no Python low-level emit.

    test/tet_laplacian.mlir                 high-level mir (hand-authored)
      --convert-mir-to-linalg               ragged scf nests + linalg
      one-shot-bufferize                    tensors -> memrefs
      --buffer-results-to-out-params        destination-passing (a kernel cannot
                                            return a buffer)
      --convert-linalg-to-loops             the NVVM pipeline does not take linalg
      --mir-batch-elements                  warp-per-element (gpu.block_id x)
      --mir-gpu-wrap                        gpu.module { gpu.func ... kernel }
      --gpu-lower-to-nvvm-pipeline          PTX

Numerics for the operator are gated separately by exec_gate.py gate 8 (against
tet_galerkin.py) and the four ragged stages by gate 7.

Run from marsir-mlir/:  python3 test/tet_gpu_pipeline.py
"""
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
LLVM = "/opt/homebrew/opt/llvm"
MIROPT = os.path.join(ROOT, "build", "tools", "mir-opt", "mir-opt")
MLIROPT = os.path.join(LLVM, "bin", "mlir-opt")

# Which arguments of tet_laplacian carry one element's data. A/Ad/AT/AdT/B/Bd/
# C/Cd (1..8) are the PKD reference tables -- identical for every element, so
# batching them would be wrong. They share a type with the metrics, which is why
# --mir-batch-elements takes an explicit marker instead of guessing from shapes.
# In a real front-end these markers come from the operator spec.
PER_ELEMENT = {0} | set(range(9, 18))     # u, then g00..g22; the out-param is appended


def run(cmd, data):
    p = subprocess.run(cmd, input=data, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"FAILED: {' '.join(cmd)}\n{p.stderr[:2000]}")
    return p.stdout


def main():
    src = open(os.path.join(HERE, "tet_laplacian.mlir")).read()
    ir = run([MIROPT, "-", "--convert-mir-to-linalg"], src)
    ir = run([MLIROPT, "-",
              "--one-shot-bufferize=bufferize-function-boundaries "
              "function-boundary-type-conversion=identity-layout-map",
              "--buffer-results-to-out-params", "--convert-linalg-to-loops"], ir)

    m = re.search(r"(func\.func @tet_laplacian\()(.*?)(\) \{)", ir, re.S)
    args = m.group(2).split(", ")
    marks = PER_ELEMENT | {len(args) - 1}          # + the appended out-parameter
    args = [a + " {mir.element}" if i in marks else a for i, a in enumerate(args)]
    ir = (ir[:m.start(2)] + ", ".join(args) + ") attributes {mir.kernel} {"
          + ir[m.end(3):])

    ir = run([MIROPT, "-", "--mir-batch-elements", "--mir-gpu-wrap"], ir)
    ir = ir.replace("gpu.module @mir_kernels {",
                    'gpu.module @mir_kernels [#nvvm.target<chip = "sm_90", O = 3>] {')
    ir = run([MLIROPT, "-",
              "--gpu-lower-to-nvvm-pipeline=cubin-chip=sm_90 cubin-format=isa"], ir)

    out = os.path.join(ROOT, "generated", "tet_laplacian_sm90.ptx")
    ptx = run([sys.executable, os.path.join(HERE, "extract_ptx.py"), out, "0"], ir)
    print(ptx.strip())

    p = open(out).read()
    checks = [
        ("kernel entry point",  p.count(".visible .entry tet_laplacian") == 1),
        ("element index (ctaid)", "ctaid.x" in p),
        ("fp64 multiply-add",   p.count("mul.rn.f64") > 0 and p.count("add.rn.f64") > 0),
        ("no tensor cores (tets do not fit m8n8k4)", p.count("mma.sync") == 0),
    ]
    for name, good in checks:
        print(f"  {'ok  ' if good else 'FAIL'} {name}")
    ok = all(g for _, g in checks)
    print("TET GPU PIPELINE PASS" if ok else "TET GPU PIPELINE FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
