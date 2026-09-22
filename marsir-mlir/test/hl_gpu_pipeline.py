#!/usr/bin/env python3
"""The Knaus operator from HIGH-LEVEL mir to a tensor-core GPU kernel, by passes
only -- no Python low-level emit anywhere in the path.

    laplacian.op -> mlir_ir.emit_full        high-level mir: mir.contract + mir.flux
      --convert-mir-to-linalg                linalg
      one-shot-bufferize                     memrefs
      --buffer-results-to-out-params         DPS (a kernel cannot return a buffer)
      dmma_schedule.mlir                     tile m8n8k4 + vectorize
      --convert-linalg-to-loops              whatever tiling left
      --mir-forward-transfers                store-to-load, so the chain is visible
      --mir-chain-contracts                  per-lane nvgpu.mma.sync + gpu.shuffle
      LICM + --mir-hoist-transfer-pairs      accumulator into a register
      --promote-buffers-to-stack             no device-side malloc in the kernel
      --mir-batch-elements                   one warp per element (grid = E)
      --mir-gpu-wrap                         gpu.module { gpu.func ... kernel }
      explicit NVVM lowering                 -> PTX

WHY THE LOWERING IS SPELLED OUT rather than --gpu-lower-to-nvvm-pipeline: that
pipeline leaves memref.subview / collapse_shape / expand_shape and nvgpu.mma.sync
unconverted for this kernel, and the index-typed loop arguments then fail LLVM
legalisation with "'llvm.cond_br' op operand #1 ... got 'index'". Adding
--expand-strided-metadata and --convert-nvgpu-to-nvvm in the right order fixes it.

Numerics are gated on the GPU by tools/run_hl_mma.cpp against the Knaus oracle.

Run from marsir-mlir/:  python3 test/hl_gpu_pipeline.py [p]
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

# U (element field), the metric, and the output are per element; Btil/Dtil/Dm/W
# are the reference operators and are shared by every element.
PER_ELEMENT = {0, 5}          # the out-param is appended and marked separately


def run(cmd, data):
    p = subprocess.run(cmd, input=data, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"FAILED: {' '.join(cmd[:3])} ...\n{p.stderr[:2000]}")
    return p.stdout


def main():
    p = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    src = run([sys.executable, "-c", f"""
import sys; sys.path.insert(0, {os.path.join(ROOT, '..', 'marsir-compiler')!r})
from marsir import parse_spec_file, synthesize
from marsir.backends import mlir_ir
ea = synthesize(parse_spec_file({os.path.join(ROOT, '..', 'marsir-compiler', 'specs', 'laplacian.op')!r}))
sys.stdout.write(mlir_ir.emit_full(ea, p={p}))
"""], "")

    ir = run([MIROPT, "-", "--convert-mir-to-linalg",
              "--one-shot-bufferize=bufferize-function-boundaries=true "
              "function-boundary-type-conversion=identity-layout-map"], src)
    ir = run([MLIROPT, "-", "--buffer-results-to-out-params"], ir)
    ir = run([MLIROPT, "-",
              f"--transform-preload-library=transform-library-paths={HERE}/dmma_schedule.mlir",
              "--transform-interpreter"], ir)
    ir = run([MLIROPT, "-", "--convert-linalg-to-loops", "--canonicalize", "--cse"], ir)
    ir = run([MIROPT, "-", "--mir-forward-transfers"], ir)
    ir = run([MIROPT, "-", "--mir-chain-contracts", "--canonicalize", "--cse"], ir)
    mma_ir = ir.count("nvgpu.mma.sync")
    shfl_ir = ir.count("gpu.shuffle")
    ir = run([MLIROPT, "-", "--loop-invariant-code-motion"], ir)
    ir = run([MIROPT, "-", "--mir-hoist-transfer-pairs"], ir)
    # Bufferization leaves memref.alloc for every temporary. Inside a kernel those
    # become DEVICE-SIDE malloc calls: with grid = E each block allocates, the 8 MB
    # device heap is gone almost immediately, malloc returns null and the kernel
    # faults with an illegal memory access. Promote them to stack allocations.
    # This has to run AFTER --mir-forward-transfers, which keys on memref.alloc to
    # decide what it may safely forward.
    ir = run([MLIROPT, "-",
              "--promote-buffers-to-stack=max-alloc-size-in-bytes=65536",
              "--canonicalize", "--cse"], ir)

    m = re.search(r"(func\.func @laplacian_apply\()(.*?)(\) \{)", ir, re.S)
    args = m.group(2).split(", ")
    marks = PER_ELEMENT | {len(args) - 1}
    args = [a + " {mir.element}" if i in marks else a for i, a in enumerate(args)]
    ir = (ir[:m.start(2)] + ", ".join(args) + ") attributes {mir.kernel} {"
          + ir[m.end(3):])

    ir = run([MIROPT, "-", "--mir-batch-elements", "--mir-gpu-wrap"], ir)
    ir = ir.replace("gpu.module @mir_kernels {",
                    'gpu.module @mir_kernels [#nvvm.target<chip = "sm_90", O = 3>] {')
    ir = run([MLIROPT, "-",
              "--convert-vector-to-scf", "--canonicalize",
              "--convert-scf-to-cf",
              "--convert-nvgpu-to-nvvm",
              "--expand-strided-metadata",
              "--lower-affine",
              "--convert-vector-to-llvm",
              "--convert-gpu-to-nvvm",
              "--convert-arith-to-llvm", "--convert-index-to-llvm",
              "--convert-cf-to-llvm",
              "--finalize-memref-to-llvm", "--reconcile-unrealized-casts"], ir)
    ir = run([MLIROPT, "-", "--gpu-module-to-binary=format=isa"], ir)

    out = os.path.join(ROOT, "generated", f"hl_full_p{p}_sm90.ptx")
    print(run([sys.executable, os.path.join(HERE, "extract_ptx.py"), out, "0"], ir).strip())

    ptx = open(out).read()
    checks = [
        ("single kernel entry", ptx.count(".visible .entry") == 1),
        ("element index (ctaid.x)", "ctaid.x" in ptx),
        ("lane index (tid.x)", "tid.x" in ptx),
        ("fp64 tensor-core mma", ptx.count("mma.sync.aligned.m8n8k4") > 0),
        ("register relayout (shfl)", ptx.count("shfl.sync") > 0),
        ("no shared memory", ptx.count(".shared") == 0),
    ]
    for name, good in checks:
        print(f"  {'ok  ' if good else 'FAIL'} {name}")
    print(f"  IR had {mma_ir} nvgpu.mma.sync and {shfl_ir} gpu.shuffle; "
          f"PTX has {ptx.count('mma.sync.aligned.m8n8k4')} mma, "
          f"{ptx.count('shfl.sync')} shfl")
    ok = all(g for _, g in checks)
    print("HIGH-LEVEL mir -> TENSOR-CORE PTX: PASS" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
