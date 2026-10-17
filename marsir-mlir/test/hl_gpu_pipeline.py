#!/usr/bin/env python3
"""The Knaus operator from HIGH-LEVEL mir to a tensor-core GPU kernel, by passes
only -- no Python low-level emit anywhere in the path.

    laplacian.op -> mlir_ir.emit_full        high-level mir: mir.contract + mir.flux
      --convert-mir-to-linalg                linalg
      one-shot-bufferize                     memrefs
      --buffer-results-to-out-params         DPS, Y accumulated in the out-param
      dmma_schedule.mlir                     tile m8n8k4 + vectorize
      --convert-linalg-to-loops              whatever tiling left
      --mir-forward-transfers                store-to-load, so the chain is visible;
                                             write-only scratch is dropped
      --mir-workgroup-buffers                remaining scratch -> shared slots, reused
      --mir-chain-contracts                  per-lane nvgpu.mma.sync + gpu.shuffle
      LICM + --mir-hoist-invariant-reads     operator fragments read once, not per face
      --mir-hoist-transfer-pairs             accumulator into a register
      --mir-unroll-loops + --mir-forward-owned  Y planes stay in registers across faces
      --mir-distribute-fills                 a shared fill is 1/32 per lane
      --mir-warp-barriers                    gpu.barrier where lanes can conflict
      --promote-buffers-to-stack             no device-side malloc in the kernel
      --mir-batch-elements                   one warp per element (grid = E)
      --mir-gpu-wrap                         gpu.func kernel, slots -> workgroup
      explicit NVVM lowering                 -> PTX

WHY THE LOWERING IS SPELLED OUT rather than --gpu-lower-to-nvvm-pipeline: that
pipeline leaves memref.subview / collapse_shape / expand_shape and nvgpu.mma.sync
unconverted for this kernel, and the index-typed loop arguments then fail LLVM
legalisation with "'llvm.cond_br' op operand #1 ... got 'index'". Adding
--expand-strided-metadata and --convert-nvgpu-to-nvvm in the right order fixes it.

Numerics are gated on the GPU by tools/run_hl_mma.cpp against the Knaus oracle.

--cellwise: the schedule for small elements (low degree, the p-multigrid's coarse
levels) instead of the tensor-core one. After the same lowering to bufferized linalg:
      --mir-cellwise-schedule      element loop, U staged in and Y out of shared
                                   memory, face loops unrolled, one phase per linalg
                                   op (fills fused into reductions), temporaries in
                                   shared scratch; n^2 lanes per element, 8 warps
      --mir-fuse-phases            neighbouring phases with lane-local dependences
                                   merged, handed-over values kept in registers
      --mir-warp-per-element       several elements per block, barriers where phases
                                   share memory
  -> hl_cellwise_p{p}_sm90.ptx: kernel laplacian_apply_p{p} (named by degree, so a
     profile tells the levels apart), launch shape published as
     laplacian_apply_p{p}__launch. --warps=W sets the warps per block (default 4
     for p >= 3, else 8: measured), written to hl_cellwise_p{p}_w{W}_sm90.ptx when
     not the default.

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
    argv = [a for a in sys.argv[1:] if not a.startswith("--")]
    p = int(argv[0]) if argv else 7
    # --no-chain: identical pipeline minus --mir-chain-contracts, so the leftover
    # vector.contract ops lower to ordinary FMA code. Same batching, wrapping,
    # lowering and kernel signature -- a control for bisecting a GPU failure
    # between the chain pass and everything around it.
    no_chain = "--no-chain" in sys.argv
    # --chain-opts=a=false,b=false  forwards switches to --mir-chain-contracts;
    # --no-hoist skips --mir-hoist-transfer-pairs; --tag=x names the output.
    # Together they bisect a numerical failure to one feature of the chain pass.
    chain_opts = next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--chain-opts=")), "")
    no_hoist = "--no-hoist" in sys.argv
    cellwise = "--cellwise" in sys.argv
    tag_arg = next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--tag=")), "")
    src = run([sys.executable, "-c", f"""
import sys; sys.path.insert(0, {os.path.join(ROOT, '..', 'marsir-compiler')!r})
from marsir import parse_spec_file, synthesize
from marsir.backends import mlir_ir
ea = synthesize(parse_spec_file({os.path.join(ROOT, '..', 'marsir-compiler', 'specs', 'laplacian.op')!r}))
sys.stdout.write(mlir_ir.emit_full(ea, p={p}, pad_faces=True))
"""], "")

    ir = run([MIROPT, "-", "--convert-mir-to-linalg",
              "--one-shot-bufferize=bufferize-function-boundaries=true "
              "function-boundary-type-conversion=identity-layout-map"], src)
    # Canonicalize first so the returned buffer is the alloc itself, not a chain
    # of memref iter_args; hoist-static-allocs then accumulates Y IN the out-
    # param instead of a private copy that is copied out at the end.
    ir = run([MLIROPT, "-", "--canonicalize",
              "--buffer-results-to-out-params=hoist-static-allocs=true"], ir)

    # The kernel and its per-element arguments, before any pass that keys on them.
    m = re.search(r"(func\.func @laplacian_apply\()(.*?)(\) \{)", ir, re.S)
    args = m.group(2).split(", ")
    marks = PER_ELEMENT | {len(args) - 1}
    args = [a + " {mir.element}" if i in marks else a for i, a in enumerate(args)]
    if cellwise:   # the element's field is staged into shared memory, its result out of it
        args[0] = args[0].replace("{mir.element}", '{mir.element, mir.stage = "in"}')
        args[-1] = args[-1].replace("{mir.element}", '{mir.element, mir.stage = "out"}')
    ir = (ir[:m.start(2)] + ", ".join(args) + ") attributes {mir.kernel} {"
          + ir[m.end(3):])
    if cellwise:
        # Measured on GH200 at 64^3 (p-multigrid V-cycle): p = 3 runs 11% faster with 4
        # warps per block (30 KB of shared memory) than with 8 (61 KB, half the blocks per SM).
        default_warps = 4 if p >= 3 else 8
        warps = next((int(a.split("=", 1)[1]) for a in sys.argv if a.startswith("--warps=")), default_warps)
        ir = run([MIROPT, "-", f"--mir-cellwise-schedule=warps={warps}", "--canonicalize", "--cse",
                  "--mir-fuse-phases", "--canonicalize", "--cse",
                  "--mir-warp-per-element", "--canonicalize", "--cse", "--loop-invariant-code-motion"], ir)
        if "--emulate" in sys.argv:
            return emulate(ir, p)
        ir = ir.replace("func.func @laplacian_apply(", f"func.func @laplacian_apply_p{p}(")
        return finish(ir, p, "" if warps == default_warps else f"_w{warps}", no_chain=True, prefix="hl_cellwise")
    ir = run([MLIROPT, "-",
              f"--transform-preload-library=transform-library-paths={HERE}/dmma_schedule.mlir",
              "--transform-interpreter"], ir)
    ir = run([MLIROPT, "-", "--convert-linalg-to-loops", "--canonicalize", "--cse"], ir)
    # Forwarding first (it keys on allocs), then what is left of the scratch
    # moves to workgroup memory, so fragments are stored there piecewise.
    ir = run([MIROPT, "-", "--mir-forward-transfers", "--mir-workgroup-buffers"], ir)
    if not no_chain:
        flag = "--mir-chain-contracts" + (
            "=" + " ".join(chain_opts.split(",")) if chain_opts else "")
        ir = run([MIROPT, "-", flag, "--canonicalize", "--cse"], ir)
    mma_ir = ir.count("nvgpu.mma.sync")
    shfl_ir = ir.count("gpu.shuffle")
    # LICM hoists pure ops only; the invariant operator reads it leaves behind are
    # hoisted next, and a second LICM takes the extracts that fed off them.
    ir = run([MLIROPT, "-", "--loop-invariant-code-motion"], ir)
    ir = run([MIROPT, "-", "--mir-hoist-invariant-reads"], ir)
    ir = run([MLIROPT, "-", "--loop-invariant-code-motion"], ir)
    if not no_hoist:
        ir = run([MIROPT, "-", "--mir-hoist-transfer-pairs"], ir)
    # Unrolled, every plane and tile has constant offsets; a lane's own
    # fragment of a Y plane then stays in registers between the faces that
    # update it, and the first direction's planes start from the zero fill.
    ir = run([MIROPT, "-", "--mir-unroll-loops", "--canonicalize", "--cse",
              "--mir-hoist-invariant-reads",   # now merges the unrolled reloads
              "--mir-forward-owned", "--canonicalize", "--cse"], ir)
    # Barriers last: they must see the final access pattern, fills included.
    ir = run([MIROPT, "-", "--mir-distribute-fills", "--mir-warp-barriers"], ir)
    barriers = ir.count("gpu.barrier")
    # Bufferization leaves memref.alloc for every temporary. Inside a kernel those
    # become DEVICE-SIDE malloc calls: with grid = E each block allocates, the 8 MB
    # device heap is gone almost immediately, malloc returns null and the kernel
    # faults with an illegal memory access. Promote them to stack allocations.
    # This has to run AFTER --mir-forward-transfers, which keys on memref.alloc to
    # decide what it may safely forward.
    ir = run([MLIROPT, "-",
              "--promote-buffers-to-stack=max-alloc-size-in-bytes=65536",
              "--canonicalize", "--cse"], ir)

    if "--emulate" in sys.argv:
        return emulate(run([MIROPT, "-", "--mir-batch-elements"], ir), p)

    tag = "_nochain" if no_chain else (("_" + tag_arg) if tag_arg else "")
    return finish(run([MIROPT, "-", "--mir-batch-elements"], ir), p, tag, no_chain,
                  stats=(mma_ir, shfl_ir, barriers))


def finish(ir, p, tag, no_chain, prefix="hl_full", stats=None):
    """The batched kernel (before --mir-gpu-wrap) -> PTX in generated/, and its checks."""
    ir = run([MIROPT, "-", "--mir-gpu-wrap"], ir)
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

    out = os.path.join(ROOT, "generated", f"{prefix}_p{p}{tag}_sm90.ptx")
    print(run([sys.executable, os.path.join(HERE, "extract_ptx.py"), out, "0"], ir).strip())

    ptx = open(out).read()
    checks = [
        ("single kernel entry", ptx.count(".visible .entry") == 1),
        ("element index (ctaid.x)", "ctaid.x" in ptx),
        ("lane index (tid.x)", no_chain or "tid.x" in ptx),
        ("fp64 tensor-core mma", no_chain or ptx.count("mma.sync.aligned.m8n8k4") > 0),
        ("register relayout (shfl)", no_chain or ptx.count("shfl.sync") > 0),
        ("no local memory", no_chain or ptx.count(".local") == 0),
        # every mma/shfl/barrier reached by all 32 lanes (see ptx_convergence.py)
        ("warp collectives convergent", subprocess.run(
            [sys.executable, os.path.join(HERE, "ptx_convergence.py"), out],
            capture_output=True).returncode == 0),
    ]
    for name, good in checks:
        print(f"  {'ok  ' if good else 'FAIL'} {name}")
    shared = sum(int(x) for x in re.findall(r"\.shared \.align \d+ \.b8 \S+\[(\d+)\]", ptx))
    print(f"  shared memory {shared} B/block, {ptx.count('bar.sync') + ptx.count('barrier.sync')} barriers in PTX")
    if stats:
        print(f"  IR had {stats[0]} nvgpu.mma.sync and {stats[1]} gpu.shuffle, {stats[2]} gpu.barrier; "
              f"PTX has {ptx.count('mma.sync.aligned.m8n8k4')} mma, {ptx.count('shfl.sync')} shfl")
    ok = all(g for _, g in checks)
    print(("HIGH-LEVEL mir -> TENSOR-CORE PTX: PASS" if prefix == "hl_full" else "HIGH-LEVEL mir -> PTX: PASS")
          if ok else "FAILED")
    return 0 if ok else 1


def emulate(ir, p):
    """--emulate: the batched kernel body -- everything the GPU build runs up to
    --mir-gpu-wrap -- with its warp primitives turned into host runtime calls,
    compiled for the CPU and run with 32 threads per element by
    tools/emu_hl_mma.cpp against the Knaus oracle."""
    E = next((int(a.split("=", 1)[1]) for a in sys.argv if a.startswith("--emu-e=")), 3)
    ir = run([MIROPT, "-", "--mir-emulate-warp"], ir)
    work = os.path.join(ROOT, "build", "emu")
    os.makedirs(work, exist_ok=True)
    # TSan reports point at lines of THIS text (the debug scopes below), so keep it.
    with open(os.path.join(work, "kernel_in.mlir"), "w") as f:
        f.write(ir)
    tsan = "--emu-tsan" in sys.argv
    ir = run([MLIROPT, "-",
              "--convert-linalg-to-loops",
              "--convert-vector-to-scf", "--canonicalize",
              "--convert-scf-to-cf",
              "--expand-strided-metadata",
              "--lower-affine",
              "--convert-vector-to-llvm",
              "--convert-arith-to-llvm", "--convert-index-to-llvm",
              "--convert-cf-to-llvm",
              "--finalize-memref-to-llvm",
              "--convert-func-to-llvm",
              "--reconcile-unrealized-casts"]
             + (["--ensure-debug-info-scope-on-llvm-func"] if tsan else []), ir)
    ll = run([os.path.join(LLVM, "bin", "mlir-translate"), "--mlir-to-llvmir"], ir)
    # --emu-tsan: ThreadSanitizer, with the lane exchanges hidden from it (see
    # emu_hl_mma.cpp), so it reports shared accesses no gpu.barrier orders. TSan
    # only instruments functions marked sanitize_thread, which clang adds for C
    # sources but mlir-translate does not.
    if tsan:
        # A short tile's out-of-bounds access lowers to llvm.masked.load/store,
        # which TSan does not instrument; scalarized, it is ordinary loads/stores.
        ll = run([os.path.join(LLVM, "bin", "opt"), "-S",
                  "-passes=scalarize-masked-mem-intrin"], ll)
        ll = re.sub(r"^(define .*\))( #\d+)? \{$", r"\1 sanitize_thread\2 {", ll,
                    flags=re.M)
    with open(os.path.join(work, "kernel.ll"), "w") as f:
        f.write(ll)
    san = ["-fsanitize=thread", "-g", "-O1"] if tsan else ["-O2"]
    cc = "/usr/bin/clang" if tsan else os.path.join(LLVM, "bin", "clang")
    run([cc, *san, "-c", os.path.join(work, "kernel.ll"),
         "-o", os.path.join(work, "kernel.o")], "")
    exe = os.path.join(work, "emu_hl_mma")
    run(["/usr/bin/clang++", "-std=c++20", *san] + (["-DMIR_EMU_TSAN"] if tsan else []) +
        [os.path.join(ROOT, "tools", "emu_hl_mma.cpp"),
         os.path.join(work, "kernel.o"), "-o", exe], "")
    env = dict(os.environ, TSAN_OPTIONS="halt_on_error=1 exitcode=66")
    r = subprocess.run([exe, str(E), str(p)], capture_output=True, text=True, env=env)
    if tsan and "ThreadSanitizer: data race" in r.stderr:
        print("\n".join(r.stderr.splitlines()[:24]))
        print("WARP EMULATION (TSan): DATA RACE -- a shared access no gpu.barrier orders")
        return 1
    print(r.stdout.strip() or r.stderr.strip())
    return r.returncode


if __name__ == "__main__":
    sys.exit(main())
