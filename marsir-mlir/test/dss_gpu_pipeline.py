#!/usr/bin/env python3
"""The cell-wise DSS from HIGH-LEVEL mir to a GPU kernel, by passes only.

    mlir_ir.emit_dss                one mir.dss op on the topology tables
      --mir-lower-cellwise              per-element parallel phases (own copies, face
                                    pairs, edge/corner stars, write-back); the
                                    reference-hex tables become constant globals
      --mir-warp-per-element        elements -> blocks (grid-stride), phases ->
                                    lane-cyclic loops, a barrier between phases,
                                    scratch -> shared
      --mir-gpu-wrap                gpu.func kernel (the tables move along)
      explicit NVVM lowering        as in hl_gpu_pipeline.py
      --mir-noalias-descriptors     the arguments' data pointers noalias, so
                                    read-only buffers load as ld.global.nc
      -> PTX

Kernel contract: dss_p{p}(values, ghost, elements, face_nbr, face_code, edge_of,
edge_off, edge_ent, vert_of, vert_off, vert_ent, out), every argument a 1-D
memref descriptor; block = 32, any grid (grid-stride over the listed elements).
Measured on GH200 at 64^3 (mars_cellwise_unstructured_dss --dss-ptx): bit-identical
to the hand kernel, 1.25-1.27 ms plain / 1.36-1.37 ms rotated against its
1.05 / 1.19 ms. Tried and measured slower: tables in shared memory (2.2-3.3 ms),
fewer blocks striding over the elements (1.47 ms), descriptors prefetched in the
staging phase (1.31 ms).

--emulate: the same kernel body with --mir-emulate-warp instead of
--mir-gpu-wrap, run on the CPU by tools/emu_dss.cpp (32 threads per element) and
compared bit for bit with the hand kernel's formula; --emu-tsan adds
ThreadSanitizer, which reports any shared access no barrier orders.

--precond: the Jacobi preconditioner instead, in the four variants the Krylov
solvers call (no dots, <a, z>, <z, z>, both; mlir_ir.emit_dss_precond_all: mir.dss
with an epilogue), one module, precond_p{p}_sm90.ptx. The CPU gate checks the
variant with both dots.

--bicgstab: the other BiCGStab kernels (vector updates as mir.cellwise; partial
folding, finish and step kernels as structured parallel loops; mlir_ir.emit_bicgstab),
bicgstab_p{p}_sm90.ptx, CPU gate tools/emu_bicg.cpp. Mapped before canonicalizing:
the one-iteration element loop of a finish kernel would otherwise fold away.

--pmg: the per-iteration kernels of the p-multigrid V-cycle 7 -> 3 -> 1 besides the
operators (Chebyshev steps, masked residual copies, transfers as mir.interpolate,
dot products; mlir_ir.emit_pmg), pmg_sm90.ptx, CPU gate tools/emu_pmg.cpp. The
operators at p = 3 and 1 come from hl_gpu_pipeline.py (p = 3 with --no-chain).

Run from marsir-mlir/:  python3 test/dss_gpu_pipeline.py [p] [--precond | --bicgstab | --pmg] [--emulate [--emu-tsan]]
"""
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
MARS = os.path.dirname(ROOT)
LLVM = "/opt/homebrew/opt/llvm"
MIROPT = os.path.join(ROOT, "build", "tools", "mir-opt", "mir-opt")
MLIROPT = os.path.join(LLVM, "bin", "mlir-opt")


def run(cmd, data):
    p = subprocess.run(cmd, input=data, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"FAILED: {' '.join(cmd[:3])} ...\n{p.stderr[:2000]}")
    return p.stdout


def main():
    argv = [a for a in sys.argv[1:] if not a.startswith("--")]
    p = int(argv[0]) if argv else 7
    precond = "--precond" in sys.argv
    bicg = "--bicgstab" in sys.argv
    pmg = "--pmg" in sys.argv
    emitter = "emit_pmg" if pmg else "emit_bicgstab" if bicg else ("emit_dss_precond_all" if precond else "emit_dss")
    src = run([sys.executable, "-c", f"""
import sys; sys.path.insert(0, {os.path.join(MARS, 'marsir-compiler')!r})
from marsir.backends import mlir_ir
sys.stdout.write(mlir_ir.{emitter}({'' if pmg else p}))
"""], "")
    early_cleanup = [] if bicg or pmg else ["--canonicalize", "--cse"]
    ir = run([MIROPT, "-", "--mir-lower-cellwise", *early_cleanup,
              "--mir-warp-per-element", "--canonicalize", "--cse",
              "--loop-invariant-code-motion"], src)
    if "--emulate" in sys.argv:
        return emulate(ir, p, precond, bicg, pmg)

    ir = run([MIROPT, "-", "--mir-gpu-wrap"], ir)
    ir = ir.replace("gpu.module @mir_kernels ",
                    'gpu.module @mir_kernels [#nvvm.target<chip = "sm_90", O = 3>] ', 1)
    ir = run([MLIROPT, "-",
              "--convert-scf-to-cf",
              "--expand-strided-metadata",
              "--lower-affine",
              "--convert-gpu-to-nvvm",
              "--convert-arith-to-llvm", "--convert-index-to-llvm",
              "--convert-cf-to-llvm",
              "--finalize-memref-to-llvm", "--reconcile-unrealized-casts"], ir)
    ir = run([MIROPT, "-", "--mir-noalias-descriptors"], ir)
    ir = run([MLIROPT, "-", "--gpu-module-to-binary=format=isa"], ir)

    name = "pmg" if pmg else f"bicgstab_p{p}" if bicg else (f"precond_p{p}" if precond else f"dss_p{p}")
    out = os.path.join(ROOT, "generated", f"{name}_sm90.ptx")
    print(run([sys.executable, os.path.join(HERE, "extract_ptx.py"), out, "0"], ir).strip())
    ptx = open(out).read()
    checks = [
        ("kernel entries", ptx.count(".visible .entry") == (20 if pmg else 15 if bicg else 4 if precond else 1)),
        ("element index (ctaid.x)", "ctaid.x" in ptx),
        ("lane index (tid.x)", "tid.x" in ptx),
        ("no local memory", ptx.count(".local") == 0),
        ("values on the read-only path", "ld.global.nc.f64" in ptx),
        ("phase barriers", bicg or pmg or ptx.count("bar.sync") + ptx.count("barrier.sync") >= (4 if p > 1 else 3)),
        ("warp collectives convergent", subprocess.run(
            [sys.executable, os.path.join(HERE, "ptx_convergence.py"), out],
            capture_output=True).returncode == 0),
    ]
    for name, good in checks:
        print(f"  {'ok  ' if good else 'FAIL'} {name}")
    shared = sum(int(x) for x in re.findall(r"\.shared \.align \d+ \.b8 \S+\[(\d+)\]", ptx))
    print(f"  shared memory {shared} B/block")
    ok = all(g for _, g in checks)
    print("HIGH-LEVEL mir.dss -> PTX: PASS" if ok else "FAILED")
    return 0 if ok else 1


def emulate(ir, p, precond, bicg, pmg=False):
    ir = run([MIROPT, "-", "--mir-emulate-warp"], ir)
    work = os.path.join(ROOT, "build", "emu_dss")
    os.makedirs(work, exist_ok=True)
    tsan = "--emu-tsan" in sys.argv
    ir = run([MLIROPT, "-",
              "--convert-linalg-to-loops",   # the emulator's NaN fill of the scratch
              "--convert-scf-to-cf",
              "--expand-strided-metadata",
              "--lower-affine",
              "--convert-arith-to-llvm", "--convert-index-to-llvm",
              "--convert-cf-to-llvm",
              "--finalize-memref-to-llvm",
              "--convert-func-to-llvm",
              "--reconcile-unrealized-casts"]
             + (["--ensure-debug-info-scope-on-llvm-func"] if tsan else []), ir)
    ll = run([os.path.join(LLVM, "bin", "mlir-translate"), "--mlir-to-llvmir"], ir)
    if tsan:   # TSan instruments only functions marked sanitize_thread
        ll = re.sub(r"^(define .*\))( #\d+)? \{$", r"\1 sanitize_thread\2 {", ll, flags=re.M)
    with open(os.path.join(work, "kernel.ll"), "w") as f:
        f.write(ll)
    san = ["-fsanitize=thread", "-g", "-O1"] if tsan else ["-O2"]
    cc = "/usr/bin/clang" if tsan else os.path.join(LLVM, "bin", "clang")
    run([cc, *san, "-c", os.path.join(work, "kernel.ll"),
         "-o", os.path.join(work, "kernel.o")], "")
    exe = os.path.join(work, f"emu_{'pmg' if pmg else 'bicg' if bicg else 'precond' if precond else 'dss'}_p{p}")
    run(["/usr/bin/clang++", "-std=c++20", "-ffp-contract=off", *san, f"-DMIR_DSS_P={p}"]
        + (["-DMIR_DSS_PRECOND"] if precond else [])
        + (["-DMIR_EMU_TSAN"] if tsan else []) +
        ["-I", MARS, os.path.join(ROOT, "tools", "emu_pmg.cpp" if pmg else "emu_bicg.cpp" if bicg else "emu_dss.cpp"),
         os.path.join(work, "kernel.o"), "-o", exe], "")
    env = dict(os.environ, TSAN_OPTIONS="halt_on_error=1 exitcode=66")
    r = subprocess.run([exe], capture_output=True, text=True, env=env)
    if tsan and "ThreadSanitizer: data race" in r.stderr:
        print("\n".join(r.stderr.splitlines()[:24]))
        print("WARP EMULATION (TSan): DATA RACE -- a shared access no gpu.barrier orders")
        return 1
    print(r.stdout.strip() or r.stderr.strip())
    return r.returncode


if __name__ == "__main__":
    sys.exit(main())
