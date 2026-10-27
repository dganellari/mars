#!/usr/bin/env python3
"""The block-structured DSS from HIGH-LEVEL mir to a GPU module, by passes only.

    mlir_ir.emit_block_dss       mir.block_sum + mir.block_skeleton, one kernel each
      --mir-lower-cellwise       flat item loops (block nodes; skeleton nodes), the
                                 bodies of MARS's mars_cellwise_block_dss.hpp in the
                                 same summation order
      --mir-warp-per-element     32 items to an element, 8 elements to a block,
                                 grid-stride (mapped before canonicalizing, which would
                                 merge the element and lane loops)
      --mir-gpu-wrap, NVVM, --mir-noalias-descriptors   as in dss_gpu_pipeline.py
      -> generated/block_dss_p{p}_b{B}_sm90.ptx

Kernel contract (every argument a 1-D memref descriptor, 256 threads per block, any
grid; run in this order on one stream):
  block_sum_p{p}_b{B}(values, out)
  block_skeleton_p{p}_b{B}(out, face_nbr, face_code, edge_off, edge_ent, vert_off, vert_ent)

--emulate: the same kernels with --mir-emulate-warp, run on the CPU by
tools/emu_block_dss.cpp against the hand bodies (bit for bit) and the assembled sums of a
rotated-block mesh; --emu-tsan adds ThreadSanitizer.

Run from marsir-mlir/:  python3 test/block_dss_gpu_pipeline.py p B [--emulate [--emu-tsan]]
"""
import os
import re
import subprocess
import sys

import dss_gpu_pipeline as d


def main():
    argv = [a for a in sys.argv[1:] if not a.startswith("--")]
    p, block = int(argv[0]), int(argv[1])
    src = d.run([sys.executable, "-c", f"""
import sys; sys.path.insert(0, {os.path.join(d.MARS, 'marsir-compiler')!r})
from marsir.backends import mlir_ir
sys.stdout.write(mlir_ir.emit_block_dss({p}, {block}))
"""], "")
    ir = d.run([d.MIROPT, "-", "--mir-lower-cellwise", "--mir-warp-per-element", "--canonicalize", "--cse",
                "--loop-invariant-code-motion"], src)
    if "--emulate" in sys.argv:
        return emulate(ir, p, block)

    ir = d.run([d.MIROPT, "-", "--mir-gpu-wrap"], ir)
    ir = ir.replace("gpu.module @mir_kernels ",
                    'gpu.module @mir_kernels [#nvvm.target<chip = "sm_90", O = 3>] ', 1)
    ir = d.run([d.MLIROPT, "-",
                "--convert-scf-to-cf",
                "--expand-strided-metadata",
                "--lower-affine",
                "--convert-gpu-to-nvvm",
                "--convert-arith-to-llvm", "--convert-index-to-llvm",
                "--convert-cf-to-llvm",
                "--finalize-memref-to-llvm", "--reconcile-unrealized-casts"], ir)
    ir = d.run([d.MIROPT, "-", "--mir-noalias-descriptors"], ir)
    ir = d.run([d.MLIROPT, "-", "--gpu-module-to-binary=format=isa"], ir)
    out = os.path.join(d.ROOT, "generated", f"block_dss_p{p}_b{block}_sm90.ptx")
    print(d.run([sys.executable, os.path.join(d.HERE, "extract_ptx.py"), out, "0"], ir).strip())
    ptx = open(out).read()
    checks = [
        ("kernel entries", ptx.count(".visible .entry") == 2),
        ("element index (ctaid.x)", "ctaid.x" in ptx),
        ("lane index (tid.x)", "tid.x" in ptx),
        ("no local memory", ptx.count(".local") == 0),
        ("no shared memory, no barriers", ".shared" not in ptx and "bar.sync" not in ptx and "barrier.sync" not in ptx),
        ("values on the read-only path (block_sum)", "ld.global.nc.f64" in ptx),
        ("warp collectives convergent", subprocess.run(
            [sys.executable, os.path.join(d.HERE, "ptx_convergence.py"), out], capture_output=True).returncode == 0),
    ]
    for name, good in checks:
        print(f"  {'ok  ' if good else 'FAIL'} {name}")
    ok = all(g for _, g in checks)
    print("HIGH-LEVEL mir.block_sum + mir.block_skeleton -> PTX: PASS" if ok else "FAILED")
    return 0 if ok else 1


def emulate(ir, p, block):
    ir = d.run([d.MIROPT, "-", "--mir-emulate-warp"], ir)
    work = os.path.join(d.ROOT, "build", "emu_block_dss")
    os.makedirs(work, exist_ok=True)
    tsan = "--emu-tsan" in sys.argv
    ir = d.run([d.MLIROPT, "-",
                "--convert-scf-to-cf",
                "--expand-strided-metadata",
                "--lower-affine",
                "--convert-arith-to-llvm", "--convert-index-to-llvm",
                "--convert-cf-to-llvm",
                "--finalize-memref-to-llvm",
                "--convert-func-to-llvm",
                "--reconcile-unrealized-casts"]
               + (["--ensure-debug-info-scope-on-llvm-func"] if tsan else []), ir)
    ll = d.run([os.path.join(d.LLVM, "bin", "mlir-translate"), "--mlir-to-llvmir"], ir)
    if tsan:   # TSan instruments only functions marked sanitize_thread
        ll = re.sub(r"^(define .*\))( #\d+)? \{$", r"\1 sanitize_thread\2 {", ll, flags=re.M)
    with open(os.path.join(work, "kernel.ll"), "w") as f:
        f.write(ll)
    san = ["-fsanitize=thread", "-g", "-O1"] if tsan else ["-O2"]
    cc = "/usr/bin/clang" if tsan else os.path.join(d.LLVM, "bin", "clang")
    d.run([cc, *san, "-c", os.path.join(work, "kernel.ll"), "-o", os.path.join(work, "kernel.o")], "")
    exe = os.path.join(work, f"emu_block_dss_p{p}_b{block}")
    blocks_per_axis = 3 if p * block <= 6 else 2
    d.run(["/usr/bin/clang++", "-std=c++20", "-ffp-contract=off", *san, f"-DMIR_BP={p}", f"-DMIR_BB={block}",
           f"-DMIR_BM={blocks_per_axis}", "-I", d.MARS, "-I", os.path.join(d.ROOT, "tools"),
           os.path.join(d.ROOT, "tools", "emu_block_dss.cpp"), os.path.join(work, "kernel.o"), "-o", exe], "")
    env = dict(os.environ, TSAN_OPTIONS="halt_on_error=1 exitcode=66")
    r = subprocess.run([exe], capture_output=True, text=True, env=env)
    if tsan and "ThreadSanitizer: data race" in r.stderr:
        print("\n".join(r.stderr.splitlines()[:24]))
        print("BLOCK DSS WARP EMULATION (TSan): DATA RACE")
        return 1
    print(r.stdout.strip() or r.stderr.strip())
    return r.returncode


if __name__ == "__main__":
    sys.exit(main())
