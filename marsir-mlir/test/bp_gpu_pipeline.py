#!/usr/bin/env python3
"""Galerkin operators from their specs (marsir-compiler/specs/bp1.op, bp3.op: the CEED
benchmark problems BP1, mass, and BP3, Laplacian), through the cell-wise schedule to
PTX, or to the CPU warp emulation.

  marsir-compiler mlir_ir.emit_galerkin    the spec as mir: B along each axis to the q^3
                                           Gauss points, the pointwise outputs, back
  --convert-mir-to-linalg, bufferize       as for the CVFEM operator (hl_gpu_pipeline.py)
  --mir-cellwise-schedule, --mir-fuse-phases, --mir-warp-per-element
  -> generated/{spec}_p{p}[_f32]_sm90.ptx  kernel {spec}_apply_p{p}[_f32]; arguments u,
     [G (E x K x q^3)], y; the 1D matrices B (q x n), Bt (n x q) [, D, Dt (q x q)] in
     constant memory, {kernel}__m0 .. __m3, filled by the host; launch shape published as
     {kernel}__launch

--emulate: tools/emu_bp.cpp runs the emulated kernel on random u, matrices and geometry
against direct six-index sums. --q=Q sets the Gauss points (default p + 2), --f32 the
precision, --warps=W the warps per block.

Specs bp1_tri / bp3_tri (geometry = trilinear) take gq, wq (q: Gauss points and weights) and the
corners X (E x 3 x 2 x 2 x 2) instead of G.

Run from marsir-mlir/:  python3 test/bp_gpu_pipeline.py {bp1|bp3|bp1_tri|bp3_tri} [p] [--emulate] [--lanes=L] [--maxnreg=N]
                        [--warps=W] [--no-pad]   (--no-pad: no scratch padding against bank conflicts)
"""
import os
import re
import sys

import hl_gpu_pipeline as hl

# The schedule and its per-target tables (lanes, warps, even-odd, packing) are the registered
# pipeline --mir-bp-schedule (lib/Pipelines.cpp). Two tables stay here: the front end's buffer
# reuse and the NVIDIA register cap, which belong to steps outside that pipeline.
# Less scratch is not always faster (docs/design/cellwise_bp_fcg.md, "Scratch reuse"): bp3
# p = 1 is 26% slower with any two of the three rewrites at the same launch shape.
NO_REUSE = {("bp3", 1)}
# BP1 p = 7 runs a further 25% faster on GH200 with a 72-register cap (28 B spilled; 40.6%
# occupancy, the shared-memory limit); BP3 is slower with a cap (it spills more).
MAXNREG = {("bp1", 7): 72}
# Per target: amd_target.py sets TARGET = "amd" (--mir-bp-schedule picks that target's tables).
TARGET = "nvidia"


def main():
    spec = sys.argv[1]
    p = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2].isdigit() else 3
    q = next((int(a.split("=", 1)[1]) for a in sys.argv if a.startswith("--q=")), p + 2)
    real = "f32" if "--f32" in sys.argv else "f64"
    # --mir-file=F: the operator's mir from F instead of the emitter (schedule experiments).
    mir_file = next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--mir-file=")), None)
    src = open(mir_file).read() if mir_file else hl.run([sys.executable, "-c", f"""
import sys; sys.path.insert(0, {os.path.join(hl.ROOT, '..', 'marsir-compiler')!r})
from marsir import parse_spec_file, synthesize
from marsir.backends import mlir_ir
ga = synthesize(parse_spec_file({os.path.join(hl.ROOT, '..', 'marsir-compiler', 'specs', spec + '.op')!r}))
sys.stdout.write(mlir_ir.emit_galerkin(ga, p={p}, q={q}, real={real!r}))
"""], "")
    reuse = "false" if (spec, p) in NO_REUSE else "true"
    ir = hl.run([hl.MIROPT, "-", f"--convert-mir-to-linalg=reuse-buffers={reuse}",
                 "--one-shot-bufferize=bufferize-function-boundaries=true "
                 "function-boundary-type-conversion=identity-layout-map"], src)
    ir = hl.run([hl.MLIROPT, "-", "--canonicalize", "--buffer-results-to-out-params=hoist-static-allocs=true"], ir)

    # Per element: u (first), the geometry G (the last input, when the spec reads any)
    # and the out-param y (last); the 1D matrices are shared. Bufferization renames the
    # arguments, so they are found by position. u and y are not staged in shared memory:
    # with one thread per 1D line the first sweep reads u once and the last writes y once.
    elements = [0, -1] + ([-2] if "%G:" in src or "%X:" in src else [])   # stored geometry, or the corners
    # The 1D matrices B, Bt [, D, Dt] follow u: constant memory, filled by the host. GLL
    # nodes and Gauss points are symmetric, so B and Bt are mirror-symmetric and D, Dt
    # antisymmetric (mir.reflect, for --mir-even-odd).
    constants = [1, 2, 3, 4] if "%D:" in src else [1, 2]
    reflect = [1, 1, -1, -1][:len(constants)]
    ir = hl.run([hl.MIROPT, "-", f"--mir-mark-kernel=func={spec}_apply elements={','.join(map(str, elements))} "
                 f"constants={','.join(map(str, constants))} reflect={','.join(map(str, reflect))}"], ir)

    opts = [f"spec={spec}", f"p={p}", f"q={q}", f"target={TARGET}"]
    for flag, key in (("--lanes=", "lanes"), ("--warps=", "warps")):
        value = next((a.split("=", 1)[1] for a in sys.argv if a.startswith(flag)), None)
        if value is not None:
            opts.append(f"{key}={value}")
    opts += ["even-odd=false"] * ("--no-even-odd" in sys.argv) + ["pad=false"] * ("--no-pad" in sys.argv)
    opts += ["lane-owned=true"] * ("--lane-owned" in sys.argv) + ["columns=true"] * ("--columns" in sys.argv)
    ir = hl.run([hl.MIROPT, "-", "--mir-bp-schedule=" + " ".join(opts)], ir)
    if "--emulate" in sys.argv:
        grad = "%D:" in src
        tri = "%X:" in src
        k = re.search(r"%G: tensor<(\d+)x", src)
        # The pointwise action of bp1_tri / bp3_tri is that of bp1 / bp3 on the geometry
        # the harness forms from the corners.
        ir = re.sub(r"(memref<[0-9x]*f(?:32|64)(?:, strided<\[[^\]]*\](?:, offset: [?\d]+)?>)?), 4>", r"\1>", ir)
        return hl.emulate(ir, p, harness="emu_bp.cpp",
                          defines=[f"-DBP_NAME={spec}_apply", f"-DBP_SPEC_{spec.split('_')[0].upper()}",
                                   f"-DBP_Q={q}", f"-DBP_K={k.group(1) if k else (6 if grad else 1)}"]
                          + (["-DBP_GRAD"] if grad else []) + (["-DBP_TRILINEAR"] if tri else [])
                          + (["-DBP_EVEN_ODD"] if "__fold" in ir else [])
                          + [f"-DBP_FOLD_M{k}=1" for k in range(4) if f"__m{k}__fold" in ir])
    suffix = "_f32" if real == "f32" else ""
    named = hl.run([hl.MIROPT, "-", f"--mir-mark-kernel=func={spec}_apply rename={spec}_apply_p{p}{suffix}"], ir)
    maxnreg = next((int(a.split("=", 1)[1]) for a in sys.argv if a.startswith("--maxnreg=")), MAXNREG.get((spec, p), 0))
    return hl.finish(named, p, suffix, no_chain=True, prefix=spec, maxnreg=maxnreg)

if __name__ == "__main__":
    sys.exit(main())
