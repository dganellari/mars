#!/usr/bin/env python3
"""The cell-wise DSS from HIGH-LEVEL mir to a GPU kernel, by passes only.

    mlir_ir.emit_dss                one mir.dss op on the topology tables
      --mir-lower-dss               per-element parallel phases (own copies, face
                                    pairs, edge/corner stars, write-back); the
                                    reference-hex tables become constant globals
      --mir-warp-per-element        elements -> blocks (grid-stride), phases ->
                                    lane-cyclic loops, a barrier between phases,
                                    scratch -> shared
      [--mir-workgroup-constants]   with --shared-tables: the reference-hex
                                    tables copied to shared memory per block
      --mir-gpu-wrap                gpu.func kernel (the tables move along)
      explicit NVVM lowering        as in hl_gpu_pipeline.py
      --mir-noalias-descriptors     the arguments' data pointers noalias, so
                                    read-only buffers load as ld.global.nc
      -> PTX

Kernel contract: dss_p{p}(values, ghost, elements, face_nbr, face_code, edge_of,
edge_off, edge_ent, vert_of, vert_off, vert_ent, out), every argument a 1-D
memref descriptor; block = 32, any grid (grid-stride over the listed elements).
--shared-tables writes dss_p{p}_shtab_sm90.ptx.

--emulate: the same kernel body with --mir-emulate-warp instead of
--mir-gpu-wrap, run on the CPU by tools/emu_dss.cpp (32 threads per element) and
compared bit for bit with the hand kernel's formula; --emu-tsan adds
ThreadSanitizer, which reports any shared access no barrier orders.

Run from marsir-mlir/:  python3 test/dss_gpu_pipeline.py [p] [--emulate [--emu-tsan]]
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
    src = run([sys.executable, "-c", f"""
import sys; sys.path.insert(0, {os.path.join(MARS, 'marsir-compiler')!r})
from marsir.backends import mlir_ir
sys.stdout.write(mlir_ir.emit_dss({p}))
"""], "")
    shared_tables = "--shared-tables" in sys.argv
    ir = run([MIROPT, "-", "--mir-lower-dss", "--canonicalize", "--cse",
              "--mir-warp-per-element", "--canonicalize", "--cse",
              "--loop-invariant-code-motion"]
             + (["--mir-workgroup-constants"] if shared_tables else []), src)
    if "--emulate" in sys.argv:
        return emulate(ir, p)

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

    out = os.path.join(ROOT, "generated", f"dss_p{p}{'_shtab' if shared_tables else ''}_sm90.ptx")
    print(run([sys.executable, os.path.join(HERE, "extract_ptx.py"), out, "0"], ir).strip())
    ptx = open(out).read()
    checks = [
        ("single kernel entry", ptx.count(".visible .entry") == 1),
        ("element index (ctaid.x)", "ctaid.x" in ptx),
        ("lane index (tid.x)", "tid.x" in ptx),
        ("no local memory", ptx.count(".local") == 0),
        ("values on the read-only path", "ld.global.nc.f64" in ptx),
        ("phase barriers", ptx.count("bar.sync") + ptx.count("barrier.sync") >= (4 if p > 1 else 3)),
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


def emulate(ir, p):
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
    exe = os.path.join(work, f"emu_dss_p{p}")
    run(["/usr/bin/clang++", "-std=c++20", *san, f"-DMIR_DSS_P={p}"]
        + (["-DMIR_EMU_TSAN"] if tsan else []) +
        ["-I", MARS, os.path.join(ROOT, "tools", "emu_dss.cpp"),
         os.path.join(work, "kernel.o"), "-o", exe], "")
    env = dict(os.environ, TSAN_OPTIONS="halt_on_error=1 exitcode=66")
    r = subprocess.run([exe], capture_output=True, text=True, env=env)
    if tsan and "ThreadSanitizer: data race" in r.stderr:
        print("\n".join(r.stderr.splitlines()[:24]))
        print("DSS WARP EMULATION (TSan): DATA RACE -- a shared access no gpu.barrier orders")
        return 1
    print(r.stdout.strip() or r.stderr.strip())
    return r.returncode


if __name__ == "__main__":
    sys.exit(main())
