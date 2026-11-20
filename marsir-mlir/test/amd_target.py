"""MARSIR kernels for AMD MI300A (gfx942): runs a pipeline script (hl_gpu_pipeline.py, or a
module that calls its finish(), e.g. bp_gpu_pipeline) with the same MARSIR passes, but
lowers through ROCDL instead of NVVM and writes generated/<kernel>_gfx942[_covN].s.

  python3 test/amd_target.py [bp_gpu_pipeline] <the pipeline's own arguments>
  AMD_COV=4 ...   code object version 4, for an older ROCm runtime (default: LLVM's, 5)

Then: llvm-mc -triple=amdgcn-amd-amdhsa -mcpu=gfx942 -filetype=obj x.s -o x.o, and on the
target clang++ -target amdgcn-amd-amdhsa -mcpu=gfx942 x.o -o x.hsaco (the ROCm assembler may
not know LLVM 19's directives, so assemble with this LLVM and only link there).
Three AMD-specific fixes, each explained where it is applied: private-space allocas,
__ockl_get_num_groups from the dispatch packet, and the tensor-core pass's target=amd.
gpu-module-to-binary is not used: it requires the ROCm device bitcode."""
import os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__))
TEST = HERE
sys.path.insert(0, TEST)
import hl_gpu_pipeline as hl
import importlib

def finish_amd(ir, p, tag, no_chain, prefix="hl_full", stats=None, direct=False, maxnreg=0):
    ir = hl.run([hl.MIROPT, "-", "--mir-gpu-wrap"], ir)
    ir = ir.replace("gpu.module @mir_kernels {",
                    'gpu.module @mir_kernels [#rocdl.target<chip = "gfx942", O = 3>] {')
    ir = hl.run([hl.MLIROPT, "-",
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
    # MLIR lowers gpu.grid_dim to __ockl_get_num_groups, a ROCm device-library function; define
    # it here from the HSA dispatch packet (workgroup_size_x/y/z: u16 at 4/6/8, grid_size_x/y/z:
    # u32 at 12/16/20), valid for code object v4 and v5, so no device library is needed.
    decl = re.search(r"^declare (\S+) @__ockl_get_num_groups\((\S+)\)(.*)$", ll, re.M)
    if decl:
        ret, arg = decl.group(1), decl.group(2)
        assert (ret, arg) == ("i64", "i32"), (ret, arg)
        body = """define internal i64 @__ockl_get_num_groups(i32 %dim) alwaysinline {
  %dp = call ptr addrspace(4) @llvm.amdgcn.dispatch.ptr()
  %d = zext i32 %dim to i64
  %goff = shl i64 %d, 2
  %gaddr = add i64 %goff, 12
  %woff = shl i64 %d, 1
  %waddr = add i64 %woff, 4
  %gp = getelementptr i8, ptr addrspace(4) %dp, i64 %gaddr
  %wp = getelementptr i8, ptr addrspace(4) %dp, i64 %waddr
  %grid = load i32, ptr addrspace(4) %gp, align 4
  %w16 = load i16, ptr addrspace(4) %wp, align 2
  %w = zext i16 %w16 to i32
  %sum = add i32 %grid, %w
  %sum1 = sub i32 %sum, 1
  %n = udiv i32 %sum1, %w
  %n64 = zext i32 %n to i64
  ret i64 %n64
}
declare ptr addrspace(4) @llvm.amdgcn.dispatch.ptr()"""
        ll = ll.replace(decl.group(0), body)
        print("  __ockl_get_num_groups defined from the dispatch packet")
    out = os.path.join(os.path.dirname(HERE), "generated", f"{prefix}_p{p}{tag}_gfx942.s")
    cov = os.environ.get("AMD_COV")   # code object version, for an older ROCm runtime
    llc = [os.path.join(hl.LLVM, "bin", "llc"), "-march=amdgcn", "-mcpu=gfx942", "-O3"]
    if cov:
        llc.append(f"--amdhsa-code-object-version={cov}")
        out = out.replace("_gfx942.s", f"_gfx942_cov{cov}.s")
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
    """The tensor-core pass with the AMD fragments; direct hand-off is NVIDIA-only."""
    cmd = list(cmd)
    for n, a in enumerate(cmd):
        if isinstance(a, str) and a.startswith("--mir-chain-contracts"):
            opts = [o for o in (a.split("=", 1)[1].split() if "=" in a else []) if not o.startswith("direct-handoff")]
            cmd[n] = "--mir-chain-contracts=" + " ".join(["target=amd"] + opts)
    return _run(cmd, ir)
hl.run = run_amd
hl.finish = finish_amd
# argv[1] may name a pipeline module that calls hl.finish (e.g. bp_gpu_pipeline); default hl.
if len(sys.argv) > 1 and sys.argv[1].endswith("_pipeline"):
    mod = importlib.import_module(sys.argv.pop(1))
    if hasattr(mod, "finish"):
        mod.finish = finish_amd
    if hasattr(mod, "run"):
        mod.run = run_amd
    mod.main()
else:
    hl.main()
