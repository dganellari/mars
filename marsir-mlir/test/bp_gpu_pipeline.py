#!/usr/bin/env python3
"""Galerkin operators from their specs (marsir-compiler/specs/bp1.op, bp3.op: the CEED
benchmark problems BP1, mass, and BP3, Laplacian), through the cell-wise schedule to
PTX, or to the CPU warp emulation.

  marsir-compiler mlir_ir.emit_galerkin    the spec as mir: B along each axis to the q^3
                                           Gauss points, the pointwise outputs, back
  --convert-mir-to-linalg, bufferize       as for the CVFEM operator (hl_gpu_pipeline.py)
  --mir-cellwise-schedule, --mir-fuse-phases, --mir-warp-per-element
  -> generated/{spec}_p{p}[_f32]_sm90.ptx  kernel {spec}_apply_p{p}[_f32]; arguments u,
     B (q x n), Bt (n x q), [D, Dt (q x q)], [G (E x K x q^3)], y; launch shape
     published as {kernel}__launch

--emulate: tools/emu_bp.cpp runs the emulated kernel on random u, matrices and geometry
against direct six-index sums. --q=Q sets the Gauss points (default p + 2), --f32 the
precision, --warps=W the warps per block.

Specs bp1_tri / bp3_tri (geometry = trilinear) take gq, wq (q: Gauss points and weights) and the
corners X (E x 3 x 2 x 2 x 2) instead of G.

Run from marsir-mlir/:  python3 test/bp_gpu_pipeline.py {bp1|bp3|bp1_tri|bp3_tri} [p] [--emulate]
"""
import os
import re
import sys

import hl_gpu_pipeline as hl


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
    ir = hl.run([hl.MIROPT, "-", "--convert-mir-to-linalg",
                 "--one-shot-bufferize=bufferize-function-boundaries=true "
                 "function-boundary-type-conversion=identity-layout-map"], src)
    ir = hl.run([hl.MLIROPT, "-", "--canonicalize", "--buffer-results-to-out-params=hoist-static-allocs=true"], ir)

    # Per element: u (first), the geometry G (the last input, when the spec reads any)
    # and the out-param y (last); the 1D matrices are shared. Bufferization renames the
    # arguments, so they are found by position. The field is staged into shared
    # memory, the result out of it.
    m = re.search(rf"(func\.func @{spec}_apply\()(.*?)(\) \{{)", ir, re.S)
    args = m.group(2).split(", ")
    args[0] += ' {mir.element, mir.stage = "in"}'
    args[-1] += ' {mir.element, mir.stage = "out"}'
    if "%G:" in src or "%X:" in src:   # stored geometry, or the corners (geometry = trilinear)
        args[-2] += " {mir.element}"
    ir = ir[:m.start(2)] + ", ".join(args) + ") attributes {mir.kernel} {" + ir[m.end(3):]

    # Warps per block: the CVFEM operator's measured default, halved (fewer elements
    # per block) until the static shared memory fits in 48 KB.
    given = next((int(a.split("=", 1)[1]) for a in sys.argv if a.startswith("--warps=")), None)
    warps = given or (4 if p >= 3 else 8)
    bufferized = ir
    while True:
        ir = hl.run([hl.MIROPT, "-", f"--mir-cellwise-schedule=warps={warps}", "--canonicalize", "--cse",
                     "--mir-fuse-phases", "--canonicalize", "--cse",
                     "--mir-warp-per-element", "--canonicalize", "--cse", "--loop-invariant-code-motion"],
                    bufferized)
        if "--emulate" in sys.argv:
            grad = "%D:" in src
            tri = "%X:" in src
            k = re.search(r"%G: tensor<(\d+)x", src)
            # The pointwise action of bp1_tri / bp3_tri is that of bp1 / bp3 on the geometry
            # the harness forms from the corners.
            return hl.emulate(ir, p, harness="emu_bp.cpp",
                              defines=[f"-DBP_NAME={spec}_apply", f"-DBP_SPEC_{spec.split('_')[0].upper()}",
                                       f"-DBP_Q={q}", f"-DBP_K={k.group(1) if k else (6 if grad else 1)}"]
                              + (["-DBP_GRAD"] if grad else []) + (["-DBP_TRILINEAR"] if tri else []))
        suffix = "_f32" if real == "f32" else ""
        named = ir.replace(f"func.func @{spec}_apply(", f"func.func @{spec}_apply_p{p}{suffix}(")
        status = hl.finish(named, p, suffix, no_chain=True, prefix=spec)
        ptx = open(os.path.join(hl.ROOT, "generated", f"{spec}_p{p}{suffix}_sm90.ptx")).read()
        shared = sum(int(x) for x in re.findall(r"\.shared \.align \d+ \.b8 \S+\[(\d+)\]", ptx))
        if shared <= 48 * 1024:
            return status
        if given is not None or warps == 1:
            print(f"  {shared} B of shared memory per block > 48 KB even with {warps} warp(s): FAILED")
            return 1
        warps //= 2
        print(f"  {shared} B of shared memory > 48 KB: retrying with {warps} warp(s) per block")


if __name__ == "__main__":
    sys.exit(main())
