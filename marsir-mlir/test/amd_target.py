"""MARSIR kernels for AMD MI300A (gfx942): runs a pipeline script (hl_gpu_pipeline.py, or a
module that calls its finish(), e.g. bp_gpu_pipeline) with the same MARSIR passes, but
lowers through ROCDL instead of NVVM and writes generated/<kernel>_gfx942[_covN].s.

  python3 test/amd_target.py [bp_gpu_pipeline] <the pipeline's own arguments>
  AMD_COV=4 ...   code object version 4, for an older ROCm runtime (default: LLVM's, 5)

Then: llvm-mc -triple=amdgcn-amd-amdhsa -mcpu=gfx942 -filetype=obj x.s -o x.o, and on the
target clang++ -target amdgcn-amd-amdhsa -mcpu=gfx942 x.o -o x.hsaco (the ROCm assembler may
not know LLVM 19's directives, so assemble with this LLVM and only link there).
The module carries the AMD target (--mir-set-target=target=amd: 64-lane wavefronts, 64 KB
LDS, MFMA fragments), which the passes read; --mir-amdgpu-grid-dim reads the grid size from
the dispatch packet. At the LLVM level, each explained where it is applied: private-space
allocas, and opt -O3 for gfx942 before llc (as clang runs it).
gpu-module-to-binary is not used: it requires the ROCm device bitcode."""
import os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__))
TEST = HERE
sys.path.insert(0, TEST)
import hl_gpu_pipeline as hl
import importlib

def finish_amd(ir, p, tag, no_chain, prefix="hl_full", stats=None, direct=False, maxnreg=0):
    ir = hl.run([hl.MIROPT, "-", "--mir-gpu-wrap", "--mir-amdgpu-grid-dim"], ir)
    ir = hl.run([hl.MLIROPT, "-", "--rocdl-attach-target=module=mir_kernels chip=gfx942 O=3",
                 "--convert-vector-to-scf", "--canonicalize", "--convert-scf-to-cf",
                 "--expand-strided-metadata", "--lower-affine",
                 "--convert-amdgpu-to-rocdl=chipset=gfx942",
                 "--convert-vector-to-llvm", "--convert-math-to-llvm",
                 "--convert-gpu-to-rocdl=chipset=gfx942",
                 "--convert-arith-to-llvm", "--convert-index-to-llvm", "--convert-cf-to-llvm",
                 "--finalize-memref-to-llvm", "--reconcile-unrealized-casts"], ir)
    # Skip gpu-module-to-binary (it wants ROCm device bitcode): lift the gpu.module body into
    # a plain LLVM-dialect module, translate it, and let llc emit gfx942 assembly.
    start = ir.index("gpu.module @mir_kernels")
    brace = ir.index("{", ir.index("]", start)) if "attributes" not in ir[start:ir.index("{", start)] else ir.index("{", ir.index("attributes", start) + 10)
    # the body opens at the first "{" after the attribute dictionary closes
    head = ir[start:]
    depth = 0; body_open = None
    i = head.index("attributes")
    j = head.index("{", i); depth = 1; k = j + 1
    while depth:                      # skip the attribute dictionary
        depth += {"{": 1, "}": -1}.get(head[k], 0); k += 1
    body_open = head.index("{", k)
    depth = 1; k = body_open + 1
    while depth:
        depth += {"{": 1, "}": -1}.get(head[k], 0); k += 1
    body = head[body_open + 1:k - 1]
    dl = re.search(r'llvm\.data_layout = "([^"]*)"', head).group(1)
    mod = ('module attributes {llvm.data_layout = "%s", llvm.target_triple = "amdgcn-amd-amdhsa"} {\n' % dl) + body + "\n}\n"
    ll = hl.run([os.path.join(hl.LLVM, "bin", "mlir-translate"), "--mlir-to-llvmir"], mod)
    # AMDGPU keeps stack memory in address space 5 (data layout A5); MLIR's memref lowering
    # emits generic-space allocas, which instruction selection rejects. Allocate in space 5
    # and cast to a generic pointer for the uses.
    def private(m):
        ind, name, rest = m.group(1), m.group(2), m.group(3)
        tmp = "%p5." + name[1:]   # numbered values cannot take a suffix
        return f"{ind}{tmp} = alloca {rest}, addrspace(5)\n{ind}{name} = addrspacecast ptr addrspace(5) {tmp} to ptr"
    ll, nalloca = re.subn(r"^(\s*)(%[\w.]+) = alloca (.*?)$", private, ll, flags=re.M)
    if nalloca: print(f"  {nalloca} allocas moved to address space 5")
    # --mir-amdgpu-grid-dim reads the grid size from the dispatch packet, so no ROCm
    # device-library function may remain.
    assert "@__ockl_" not in ll, "a ROCm device-library call survived the lowering"
    out = os.path.join(os.path.dirname(HERE), "generated", f"{prefix}_p{p}{tag}_gfx942.s")
    cov = os.environ.get("AMD_COV")   # code object version, for an older ROCm runtime
    llc = [os.path.join(hl.LLVM, "bin", "llc"), "-march=amdgcn", "-mcpu=gfx942", "-O3"]
    if cov:
        llc.append(f"--amdhsa-code-object-version={cov}")
        out = out.replace("_gfx942.s", f"_gfx942_cov{cov}.s")
    # Kernel arrays arrive as generic pointers inside MLIR's memref descriptors, and llc alone
    # keeps loads through them as flat_load (address space resolved per access, and each one
    # also waits like an LDS access). opt's AMDGPU pipeline, as clang runs it for HIP, proves
    # them global (global_load). Measured on MI300A (Beverin 757106): BP3 p=3 1.63 -> 1.18 ms,
    # p=7 2.07 -> 1.56 ms, 15 fewer registers.
    ll = hl.run([os.path.join(hl.LLVM, "bin", "opt"), "-mtriple=amdgcn-amd-amdhsa", "-mcpu=gfx942",
                 "-passes=default<O3>", "-S", "-o", "-"], ll)
    if os.environ.get("AMD_DUMP_LL"):   # the LLVM IR llc receives (for reading)
        open(os.environ["AMD_DUMP_LL"], "w").write(ll)
    text = hl.run(llc + ["-o", "-"], ll)
    open(out, "w").write(text)
    print(f"wrote {out}: {text.count(chr(10))} lines, "
          f"{text.count('v_fma_f64')} v_fma_f64, {text.count('ds_bpermute')} ds_bpermute, "
          f"{text.count('v_mfma_f64_4x4x4_4b_f64')} v_mfma_f64_4x4x4_4b_f64, {text.count('s_barrier')} s_barrier, {text.count('ds_read')+text.count('ds_write')} LDS ops; "
          f"nvvm left in IR: {'nvvm' in mod}")
    for key in ('.vgpr_count', '.sgpr_count', '.group_segment_fixed_size', '.wavefront_size', '.private_segment_fixed_size'):
        m = re.search(re.escape(key) + r':\s*(\d+)', text)
        if m: print(f"  {key[1:]} = {m.group(1)}")

_run = hl.run
def run_amd(cmd, ir):
    """Every mir-opt call first records the AMD target on the module (--mir-set-target), so
    the passes read a 64-lane wavefront and 64 KB of LDS and the tensor-core pass picks the
    AMD fragments. Direct hand-off is NVIDIA-only, so it is dropped."""
    cmd = list(cmd)
    if cmd and cmd[0] == hl.MIROPT and "-" in cmd:
        cmd.insert(cmd.index("-") + 1, "--mir-set-target=target=amd")
    for n, a in enumerate(cmd):
        if isinstance(a, str) and a.startswith("--mir-chain-contracts="):
            opts = [o for o in a.split("=", 1)[1].split() if not o.startswith("direct-handoff")]
            cmd[n] = "--mir-chain-contracts" + ("=" + " ".join(opts) if opts else "")
    return _run(cmd, ir)
hl.run = run_amd
hl.finish = finish_amd
# argv[1] may name a pipeline module that calls hl.finish (e.g. bp_gpu_pipeline); default hl.
if len(sys.argv) > 1 and sys.argv[1].endswith("_pipeline"):
    mod = importlib.import_module(sys.argv.pop(1))
    if hasattr(mod, "finish"):
        mod.finish = finish_amd
    if hasattr(mod, "TARGET"):   # the pipeline's per-target tuning tables
        mod.TARGET = "amd"
    if hasattr(mod, "run"):
        mod.run = run_amd
    mod.main()
else:
    hl.main()
