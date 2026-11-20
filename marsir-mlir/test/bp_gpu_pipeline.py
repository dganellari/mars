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

# Less scratch is not always faster. Measured on GH200 (docs/design/cellwise_bp_fcg.md,
# "Scratch reuse"): bp3 p = 1 is 26% slower with any two of the three rewrites at the same
# launch shape, and bp1 p = 7 5% slower with largest-first packing (4 instead of 2
# elements per block). Those two keep the forms measured fastest.
NO_REUSE = {("bp3", 1)}
FIRST_FIT = {("bp3", 1), ("bp1", 7)}
# The even-odd split (--mir-even-odd) is faster at p = 7 (BP3 10%, BP1 32%, BP3 trilinear
# 12%) and for bp3_tri p = 3 (3%), but bp3 p = 3 is 6% slower with it and p = 1 gains
# nothing (same doc). bp1_tri is not timed; it follows bp1.
EVEN_ODD = {(s, 7) for s in ("bp1", "bp3", "bp1_tri", "bp3_tri")} | {("bp3_tri", 3)}
# At p = 7 (q = 9) scratch per element limits the resident warps (BP3: four q^3 buffers,
# 8 warps per SM), so two warps share an element: GH200 operator -26% BP3, -16% BP1, and
# BP1 a further -25% with a 72-register cap (28 B spilled; 40.6% occupancy, the shared-
# memory limit). BP3 is slower with a cap (it spills more). Three or four warps lose the
# register blocking and are 2.4x slower. --lanes=L / --maxnreg=N override.
LANES = {("bp1", 7): 64, ("bp3", 7): 64}
MAXNREG = {("bp1", 7): 72}
# Per target: amd_target.py sets TARGET = "amd". Measured on MI300A (Beverin 757112): BP3 p = 3
# runs 1.01 ms with one wave per block against 1.18 ms with four (the NVIDIA default). Shared
# memory per element limits the waves per CU either way; smaller blocks pack it better.
TARGET = "nvidia"
WARPS_AMD = {("bp3", 3): 1}


def main():
    spec = sys.argv[1]
    p = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2].isdigit() else 3
    q = next((int(a.split("=", 1)[1]) for a in sys.argv if a.startswith("--q=")), p + 2)
    real = "f32" if "--f32" in sys.argv else "f64"
    src = hl.run([sys.executable, "-c", f"""
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

    # Warps per block: the CVFEM operator's measured default, which
    # --mir-warp-per-element=fit-shared halves (fewer elements per block) until the static
    # shared memory fits the target's limit; an explicit --warps=W is used as given.
    given = next((int(a.split("=", 1)[1]) for a in sys.argv if a.startswith("--warps=")), None)
    # --lanes=L: lanes per element; a multiple of 32 puts L / 32 warps on each element.
    lanes = next((int(a.split("=", 1)[1]) for a in sys.argv if a.startswith("--lanes=")), LANES.get((spec, p), 0))
    if lanes > 32 and given is None:
        given = lanes // 32
    if given is None and TARGET == "amd":
        given = WARPS_AMD.get((spec, p))
    warps = given or (4 if p >= 3 else 8)
    # One thread per 1D line (max-block = q), unrolled. The matrices' constant address
    # space steers which dimension a thread keeps, so the CPU gate schedules with it too
    # and only drops it before the CPU build.
    sched = (f"--mir-cellwise-schedule=warps={warps} max-block={q} unroll-reductions=true const-space=4"
             + (f" lanes={lanes}" if lanes else ""))
    even_odd = ["--mir-even-odd", "--canonicalize", "--cse"] if (spec, p) in EVEN_ODD and "--no-even-odd" not in sys.argv else []
    ir = hl.run([hl.MIROPT, "-", sched, "--canonicalize", "--cse", *even_odd,
                 "--mir-fuse-phases", "--canonicalize", "--cse",
                 f"--mir-warp-per-element=largest-first={'false' if (spec, p) in FIRST_FIT else 'true'} "
                 f"fit-shared={'false' if given is not None else 'true'} "
                 f"pad-scratch={'false' if '--no-pad' in sys.argv else 'true'}",
                 "--canonicalize", "--cse", "--loop-invariant-code-motion"], ir)
    if "--emulate" in sys.argv:
        grad = "%D:" in src
        tri = "%X:" in src
        k = re.search(r"%G: tensor<(\d+)x", src)
        # The pointwise action of bp1_tri / bp3_tri is that of bp1 / bp3 on the geometry
        # the harness forms from the corners.
        ir = re.sub(r"(memref<[0-9x]*f(?:32|64)), 4>", r"\1>", ir)
        return hl.emulate(ir, p, harness="emu_bp.cpp",
                          defines=[f"-DBP_NAME={spec}_apply", f"-DBP_SPEC_{spec.split('_')[0].upper()}",
                                   f"-DBP_Q={q}", f"-DBP_K={k.group(1) if k else (6 if grad else 1)}"]
                          + (["-DBP_GRAD"] if grad else []) + (["-DBP_TRILINEAR"] if tri else [])
                          + (["-DBP_EVEN_ODD"] if "__fold" in ir else []))
    suffix = "_f32" if real == "f32" else ""
    named = hl.run([hl.MIROPT, "-", f"--mir-mark-kernel=func={spec}_apply rename={spec}_apply_p{p}{suffix}"], ir)
    maxnreg = next((int(a.split("=", 1)[1]) for a in sys.argv if a.startswith("--maxnreg=")), MAXNREG.get((spec, p), 0))
    return hl.finish(named, p, suffix, no_chain=True, prefix=spec, maxnreg=maxnreg)

if __name__ == "__main__":
    sys.exit(main())
