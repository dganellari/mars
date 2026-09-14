#!/usr/bin/env python3
"""Numerical EXECUTION gate for the mir dialect lowering.

Lowers mir IR through the CPU pipeline (same --convert-mir-to-linalg used by the
tensor-core path, incl. the collapse/matmul unfolding), EXECUTES it with
mlir-cpu-runner, and compares against a NumPy oracle:

  gate 1: mir.contract           == einsum('ip,pjk->ijk', D, u)
  gate 2: PA chain with the REAL CVFEM flux
          y = W^T . ( g2*deriv + g0*dt2 + g1*dt1 ) contracted back
          where deriv/dt1/dt2 are mir.contract along axes 0/1/2.

Run from marsir-mlir/:  python3 test/exec_gate.py       (needs build/tools/mir-opt)
"""

import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
LLVM = "/opt/homebrew/opt/llvm"
MIROPT = os.path.join(ROOT, "build", "tools", "mir-opt", "mir-opt")
MLIROPT = os.path.join(LLVM, "bin", "mlir-opt")
RUNNER = os.path.join(LLVM, "bin", "mlir-cpu-runner")
LIBS = [os.path.join(LLVM, "lib", "libmlir_runner_utils.dylib"),
        os.path.join(LLVM, "lib", "libmlir_c_runner_utils.dylib")]

N = 4  # order p=3 -> n=4; small enough for readable dense constants


def fmt_tensor(a):
    """numpy array -> MLIR dense<...> literal (nested brackets)."""
    if a.ndim == 1:
        return "[" + ", ".join(f"{v:.17e}" for v in a) + "]"
    return "[" + ", ".join(fmt_tensor(s) for s in a) + "]"


def cst(name, a):
    shape = "x".join(str(d) for d in a.shape) + "xf64"
    return f"    %{name} = arith.constant dense<{fmt_tensor(a)}> : tensor<{shape}>"


def run_pipeline(payload):
    """mir -> linalg -> loops -> LLVM -> execute; return printed floats."""
    env = dict(os.environ, PATH=f"{LLVM}/bin:" + os.environ.get("PATH", ""))
    p1 = subprocess.run(
        [MIROPT, "-", "--convert-mir-to-linalg",
         "--one-shot-bufferize=bufferize-function-boundaries=true "
         "function-boundary-type-conversion=identity-layout-map"],
        input=payload, capture_output=True, text=True, env=env)
    if p1.returncode:
        sys.exit(f"mir-opt failed:\n{p1.stderr}")
    p2 = subprocess.run(
        [MLIROPT, "-", "--convert-linalg-to-loops", "--convert-scf-to-cf",
         "--expand-strided-metadata", "--lower-affine", "--finalize-memref-to-llvm",
         "--convert-arith-to-llvm", "--convert-math-to-llvm",
         "--convert-func-to-llvm", "--convert-cf-to-llvm",
         "--reconcile-unrealized-casts"],
        input=p1.stdout, capture_output=True, text=True, env=env)
    if p2.returncode:
        sys.exit(f"mlir-opt (CPU lowering) failed:\n{p2.stderr}")
    p3 = subprocess.run(
        [RUNNER, "-e", "main", "-entry-point-result=void"]
        + [f"-shared-libs={l}" for l in LIBS],
        input=p2.stdout, capture_output=True, text=True, env=env)
    if p3.returncode:
        sys.exit(f"mlir-cpu-runner failed:\n{p3.stderr}")
    vals = []
    for line in p3.stdout.splitlines():
        if "base@" in line or "sizes" in line or "rank" in line:
            continue  # printMemrefF64 metadata header, not data
        if "[" in line:
            for tok in line.replace("[", " ").replace("]", " ").replace(",", " ").split():
                try:
                    vals.append(float(tok))
                except ValueError:
                    pass
    return np.array(vals)


def diff_block(result_var, expected_arr, n=N):
    """IR that prints |result - expected| * 1e12 (a second mir.flux -- we dogfood
    the dialect for the verification math). Printed values < 1.0 <=> err < 1e-12,
    immune to printMemrefF64's ~7-digit output precision."""
    return f"""{cst("exp", expected_arr)}
    %scale = arith.constant 1.0e12 : f64
    %diff = mir.flux ins(%{result_var}, %exp)
         : (tensor<{n}x{n}x{n}xf64>, tensor<{n}x{n}x{n}xf64>) -> tensor<{n}x{n}x{n}xf64> {{
    ^bb0(%a: f64, %b: f64):
      %d = arith.subf %a, %b : f64
      %ad = math.absf %d : f64
      %sd = arith.mulf %ad, %scale : f64
      mir.yield %sd : f64
    }}
    %p = tensor.cast %diff : tensor<{n}x{n}x{n}xf64> to tensor<*xf64>
    call @printMemrefF64(%p) : (tensor<*xf64>) -> ()"""


def gate1():
    rng = np.random.RandomState(3)
    D = rng.uniform(-1, 1, (N, N))
    u = rng.uniform(-1, 1, (N, N, N))
    expected = np.einsum("ip,pjk->ijk", D, u)

    payload = f"""
func.func private @printMemrefF64(tensor<*xf64>)
func.func @main() {{
{cst("D", D)}
{cst("u", u)}
    %y = mir.contract %u, %D {{axis = 0 : i64}}
         : (tensor<{N}x{N}x{N}xf64>, tensor<{N}x{N}xf64>) -> tensor<{N}x{N}x{N}xf64>
{diff_block("y", expected)}
    return
}}
"""
    got = run_pipeline(payload)
    scaled = np.max(np.abs(got))
    ok = len(got) == N ** 3 and scaled < 1.0
    print(f"gate 1 (mir.contract == einsum)         : max|err| = {scaled:.3e}e-12  {'PASS' if ok else 'FAIL'}")
    return ok


def gate2():
    rng = np.random.RandomState(7)
    Dt = rng.uniform(-1, 1, (N, N))   # normal-derivative operator (Dtil)
    Dm = rng.uniform(-1, 1, (N, N))   # tangential derivative operator (D)
    W = rng.uniform(-1, 1, (N, N))    # integration operator
    u = rng.uniform(-1, 1, (N, N, N))
    g0 = rng.uniform(-1, 1, (N, N, N))
    g1 = rng.uniform(-1, 1, (N, N, N))
    g2 = rng.uniform(-1, 1, (N, N, N))

    deriv = np.einsum("ip,pjk->ijk", Dt, u)
    dt1 = np.einsum("jp,ipk->ijk", Dm, u)
    dt2 = np.einsum("kp,ijp->ijk", Dm, u)
    flux = g2 * deriv + g0 * dt2 + g1 * dt1          # the AUTHORED CVFEM flux
    expected = np.einsum("ip,pjk->ijk", W, flux)     # transpose/integrate sweep

    payload = f"""
func.func private @printMemrefF64(tensor<*xf64>)
func.func @main() {{
{cst("Dt", Dt)}
{cst("Dm", Dm)}
{cst("W", W)}
{cst("u", u)}
{cst("g0", g0)}
{cst("g1", g1)}
{cst("g2", g2)}
    %deriv = mir.contract %u, %Dt {{axis = 0 : i64}}
         : (tensor<{N}x{N}x{N}xf64>, tensor<{N}x{N}xf64>) -> tensor<{N}x{N}x{N}xf64>
    %dt1 = mir.contract %u, %Dm {{axis = 1 : i64}}
         : (tensor<{N}x{N}x{N}xf64>, tensor<{N}x{N}xf64>) -> tensor<{N}x{N}x{N}xf64>
    %dt2 = mir.contract %u, %Dm {{axis = 2 : i64}}
         : (tensor<{N}x{N}x{N}xf64>, tensor<{N}x{N}xf64>) -> tensor<{N}x{N}x{N}xf64>
    %flux = mir.flux ins(%deriv, %dt1, %dt2, %g0, %g1, %g2)
         : (tensor<{N}x{N}x{N}xf64>, tensor<{N}x{N}x{N}xf64>, tensor<{N}x{N}x{N}xf64>,
            tensor<{N}x{N}x{N}xf64>, tensor<{N}x{N}x{N}xf64>, tensor<{N}x{N}x{N}xf64>)
           -> tensor<{N}x{N}x{N}xf64> {{
    ^bb0(%d: f64, %t1: f64, %t2: f64, %a0: f64, %a1: f64, %a2: f64):
      %m2 = arith.mulf %a2, %d : f64
      %m0 = arith.mulf %a0, %t2 : f64
      %m1 = arith.mulf %a1, %t1 : f64
      %s0 = arith.addf %m2, %m0 : f64
      %s1 = arith.addf %s0, %m1 : f64
      mir.yield %s1 : f64
    }}
    %y = mir.contract %flux, %W {{axis = 0 : i64}}
         : (tensor<{N}x{N}x{N}xf64>, tensor<{N}x{N}xf64>) -> tensor<{N}x{N}x{N}xf64>
{diff_block("y", expected)}
    return
}}
"""
    got = run_pipeline(payload)
    scaled = np.max(np.abs(got))
    ok = len(got) == N ** 3 and scaled < 1.0
    print(f"gate 2 (PA chain, real CVFEM flux)      : max|err| = {scaled:.3e}e-12  {'PASS' if ok else 'FAIL'}")
    return ok


def gate3():
    """THE END-TO-END GATE: specs/laplacian.op -> Python front-end -> mir dialect
    (--backend mlir) -> CPU execution == NumPy oracle. Validates the whole
    two-stage architecture: one .op feeds Stage 1 (source) AND Stage 2 (MLIR)."""
    compiler = os.path.abspath(os.path.join(ROOT, "..", "marsir-compiler"))
    emitted = subprocess.run(
        [sys.executable, "marsir-emit.py", "specs/laplacian.op", "--backend", "mlir"],
        cwd=compiler, capture_output=True, text=True)
    if emitted.returncode:
        sys.exit(f"marsir-emit failed:\n{emitted.stderr}")
    pa_func = emitted.stdout  # func.func @laplacian_pa(%u,%Dtil,%W,%D,%g0,%g1,%g2), n=8

    n = 8
    rng = np.random.RandomState(11)
    u = rng.uniform(-1, 1, (n, n, n))
    Dt = rng.uniform(-1, 1, (n, n))
    W = rng.uniform(-1, 1, (n, n))
    Dm = rng.uniform(-1, 1, (n, n))
    g0 = rng.uniform(-1, 1, (n, n, n))
    g1 = rng.uniform(-1, 1, (n, n, n))
    g2 = rng.uniform(-1, 1, (n, n, n))

    deriv = np.einsum("ip,pjk->ijk", Dt, u)
    dt1 = np.einsum("jp,ipk->ijk", Dm, u)
    dt2 = np.einsum("kp,ijp->ijk", Dm, u)
    flux = g2 * deriv + g0 * dt2 + g1 * dt1
    expected = np.einsum("ip,pjk->ijk", W, flux)

    t3, t2 = f"tensor<{n}x{n}x{n}xf64>", f"tensor<{n}x{n}xf64>"
    payload = f"""{pa_func}
func.func private @printMemrefF64(tensor<*xf64>)
func.func @main() {{
{cst("u", u)}
{cst("Dt", Dt)}
{cst("W", W)}
{cst("Dm", Dm)}
{cst("g0", g0)}
{cst("g1", g1)}
{cst("g2", g2)}
    %y = call @laplacian_pa(%u, %Dt, %W, %Dm, %g0, %g1, %g2)
         : ({t3}, {t2}, {t2}, {t2}, {t3}, {t3}, {t3}) -> {t3}
{diff_block("y", expected, n)}
    return
}}
"""
    got = run_pipeline(payload)
    scaled = np.max(np.abs(got))
    ok = len(got) == n ** 3 and scaled < 1.0
    print(f"gate 3 (.op -> mir dialect -> executed)  : max|err| = {scaled:.3e}e-12  {'PASS' if ok else 'FAIL'}")
    return ok


def gate4():
    """Ragged (collapsed-tet) sweep: mir.simplex_contract executed vs a NumPy
    triangular reference -- out[p,q,k] = sum_{r<=D-p-q} u[p,q,r]*C[p,q,r,k],
    zero outside the simplex set. The bounds depend on OUTPUT indices, which is
    exactly what linalg cannot express (hence the dedicated op + scf lowering)."""
    D = 3
    W, nq = D + 1, D + 1
    rng = np.random.RandomState(5)
    u = rng.uniform(-1, 1, (W, W, W))
    C = rng.uniform(-1, 1, (W, W, W, nq))
    expected = np.zeros((W, W, nq))
    for p in range(W):
        for q in range(W - p):
            for k in range(nq):
                expected[p, q, k] = sum(u[p, q, r] * C[p, q, r, k]
                                        for r in range(W - p - q))

    payload = f"""
func.func private @printMemrefF64(tensor<*xf64>)
func.func @main() {{
{cst("u", u)}
{cst("C", C)}
    %y = mir.simplex_contract %u, %C {{degree = {D} : i64}}
         : (tensor<{W}x{W}x{W}xf64>, tensor<{W}x{W}x{W}x{nq}xf64>) -> tensor<{W}x{W}x{nq}xf64>
{cst("exp", expected)}
    %scale = arith.constant 1.0e12 : f64
    %diff = mir.flux ins(%y, %exp)
         : (tensor<{W}x{W}x{nq}xf64>, tensor<{W}x{W}x{nq}xf64>) -> tensor<{W}x{W}x{nq}xf64> {{
    ^bb0(%a: f64, %b: f64):
      %d = arith.subf %a, %b : f64
      %ad = math.absf %d : f64
      %sd = arith.mulf %ad, %scale : f64
      mir.yield %sd : f64
    }}
    %pr = tensor.cast %diff : tensor<{W}x{W}x{nq}xf64> to tensor<*xf64>
    call @printMemrefF64(%pr) : (tensor<*xf64>) -> ()
    return
}}
"""
    got = run_pipeline(payload)
    scaled = np.max(np.abs(got))
    ok = len(got) == W * W * nq and scaled < 1.0
    print(f"gate 4 (ragged tet simplex_contract)    : max|err| = {scaled:.3e}e-12  {'PASS' if ok else 'FAIL'}")
    return ok





def gate5():
    """THE FULL OPERATOR: laplacian.op -> emit_full (whole Knaus Alg-2: all 3
    directions, per-face metric flux, W integrations, +/- plane scatter) ->
    CPU-executed == a direct NumPy port of applyHoCvfemElement."""
    compiler = os.path.abspath(os.path.join(ROOT, "..", "marsir-compiler"))
    sys.path.insert(0, compiler)
    from marsir import parse_spec_file, synthesize
    from marsir.backends import mlir_ir

    p = 3
    n, P = p + 1, p
    ea = synthesize(parse_spec_file(os.path.join(compiler, "specs", "laplacian.op")))
    apply_func = mlir_ir.emit_full(ea, p=p)

    rng = np.random.RandomState(21)
    u = rng.uniform(-1, 1, (n, n, n))
    Btil = rng.uniform(-1, 1, (P, n))
    Dtil = rng.uniform(-1, 1, (P, n))
    Dm = rng.uniform(-1, 1, (n, n))
    W = rng.uniform(-1, 1, (n, n))
    G = rng.uniform(-1, 1, (3, P, n, n, 3))

    # NumPy port of the Knaus Alg-2 host reference (applyHoCvfemElement).
    y = np.zeros((n, n, n))
    for d in range(3):
        U = np.moveaxis(u, d, 0)
        Y = np.moveaxis(y, d, 0)           # view: += updates y
        for l in range(P):
            interp = np.einsum("q,qsr->sr", Btil[l], U)
            deriv = np.einsum("q,qsr->sr", Dtil[l], U)
            dt2 = np.einsum("rq,sq->sr", Dm, interp)
            dt1 = np.einsum("sq,qr->sr", Dm, interp)
            g = G[d, l]
            flux = g[..., 2] * deriv + g[..., 0] * dt2 + g[..., 1] * dt1
            tmp = np.einsum("rq,sq->sr", W, flux)
            intf = np.einsum("sq,qr->sr", W, tmp)
            Y[l] -= intf
            Y[l + 1] += intf
    expected = y

    t3 = f"tensor<{n}x{n}x{n}xf64>"
    t2 = f"tensor<{n}x{n}xf64>"
    tPn = f"tensor<{P}x{n}xf64>"
    tG = f"tensor<3x{P}x{n}x{n}x3xf64>"
    payload = f"""{apply_func}
func.func private @printMemrefF64(tensor<*xf64>)
func.func @main() {{
{cst("u", u)}
{cst("Bt", Btil)}
{cst("Dt", Dtil)}
{cst("W", W)}
{cst("Dm", Dm)}
{cst("G", G)}
    %y = call @laplacian_apply(%u, %Bt, %Dt, %W, %Dm, %G)
         : ({t3}, {tPn}, {tPn}, {t2}, {t2}, {tG}) -> {t3}
{diff_block("y", expected, n)}
    return
}}
"""
    got = run_pipeline(payload)
    scaled = np.max(np.abs(got))
    ok = len(got) == n ** 3 and scaled < 1.0
    print(f"gate 5 (FULL operator == Knaus oracle)   : max|err| = {scaled:.3e}e-12  {'PASS' if ok else 'FAIL'}")
    return ok




def gate6():
    """FUSED BATCHED operator: emit_full_batched (one scf.parallel over
    elements, mir inside via bufferization boundary ops) executed on CPU for
    E=4 elements == per-element NumPy oracle. Validates the exact form that
    lowers to the single GPU kernel."""
    compiler = os.path.abspath(os.path.join(ROOT, "..", "marsir-compiler"))
    sys.path.insert(0, compiler)
    from marsir import parse_spec_file, synthesize
    from marsir.backends import mlir_ir

    p, E = 3, 4
    n, P = p + 1, p
    ea = synthesize(parse_spec_file(os.path.join(compiler, "specs", "laplacian.op")))
    fn = mlir_ir.emit_full_batched(ea, p=p, tpb=2)   # 2 blocks x 2 threads

    rng = np.random.RandomState(33)
    U = rng.uniform(-1, 1, (E, n, n, n))
    Btil = rng.uniform(-1, 1, (P, n))
    Dtil = rng.uniform(-1, 1, (P, n))
    Dm = rng.uniform(-1, 1, (n, n))
    W = rng.uniform(-1, 1, (n, n))
    G = rng.uniform(-1, 1, (E, 3, P, n, n, 3))

    def oracle(u, g):
        y = np.zeros((n, n, n))
        for d in range(3):
            Uv, Yv = np.moveaxis(u, d, 0), np.moveaxis(y, d, 0)
            for l in range(P):
                interp = np.einsum("q,qsr->sr", Btil[l], Uv)
                deriv = np.einsum("q,qsr->sr", Dtil[l], Uv)
                dt2 = np.einsum("rq,sq->sr", Dm, interp)
                dt1 = np.einsum("sq,qr->sr", Dm, interp)
                gg = g[d, l]
                flux = gg[..., 2] * deriv + gg[..., 0] * dt2 + gg[..., 1] * dt1
                tmp = np.einsum("rq,sq->sr", W, flux)
                intf = np.einsum("sq,qr->sr", W, tmp)
                Yv[l] -= intf
                Yv[l + 1] += intf
        return y
    expected = np.stack([oracle(U[e], G[e]) for e in range(E)])

    tU = f"tensor<{E}x{n}x{n}x{n}xf64>"
    mUs = f"memref<{E}x{n}x{n}x{n}xf64>"
    mUd = f"memref<?x{n}x{n}x{n}xf64>"
    tGt = f"tensor<{E}x3x{P}x{n}x{n}x3xf64>"
    mGs = f"memref<{E}x3x{P}x{n}x{n}x3xf64>"
    mGd = f"memref<?x3x{P}x{n}x{n}x3xf64>"
    t2m, tPnm = f"memref<{n}x{n}xf64>", f"memref<{P}x{n}xf64>"

    def buf(name, arr, tty, mty):
        return (f"{cst(name + 'c', arr)}\n"
                f"    %{name}m = memref.alloc() : {mty}\n"
                f"    bufferization.materialize_in_destination %{name}c in "
                f"writable %{name}m : (tensor<{'x'.join(str(d) for d in arr.shape)}xf64>, {mty}) -> ()")

    payload = f"""{fn}
func.func private @printMemrefF64(tensor<*xf64>)
func.func @main() {{
{buf("U", U, tU, mUs)}
{buf("Y", np.zeros_like(U), tU, mUs)}
{buf("Bt", Btil, "", tPnm)}
{buf("Dt", Dtil, "", tPnm)}
{buf("W", W, "", t2m)}
{buf("Dm", Dm, "", t2m)}
{buf("G", G, tGt, mGs)}
    %Ud = memref.cast %Um : {mUs} to {mUd}
    %Yd = memref.cast %Ym : {mUs} to {mUd}
    %Gd = memref.cast %Gm : {mGs} to {mGd}
    call @laplacian_apply_batched(%Ud, %Yd, %Btm, %Dtm, %Wm, %Dmm, %Gd)
        : ({mUd}, {mUd}, {tPnm}, {tPnm}, {t2m}, {t2m}, {mGd}) -> ()
    %Yt = bufferization.to_tensor %Ym restrict : {mUs}
{cst("exp4", expected)}
    %scale = arith.constant 1.0e12 : f64
    %diff = mir.flux ins(%Yt, %exp4) : ({tU}, {tU}) -> {tU} {{
    ^bb0(%a: f64, %b: f64):
      %d = arith.subf %a, %b : f64
      %ad = math.absf %d : f64
      %sd = arith.mulf %ad, %scale : f64
      mir.yield %sd : f64
    }}
    %pr = tensor.cast %diff : {tU} to tensor<*xf64>
    call @printMemrefF64(%pr) : (tensor<*xf64>) -> ()
    return
}}
"""
    got = run_pipeline(payload)
    scaled = np.max(np.abs(got))
    ok = len(got) == E * n ** 3 and scaled < 1.0
    print(f"gate 6 (FUSED batched operator, E=4)     : max|err| = {scaled:.3e}e-12  {'PASS' if ok else 'FAIL'}")
    return ok


def gate7():
    """The four RAGGED stages of the collapsed tet operator, each executed and
    compared against a NumPy reference. gate4 covers only axis=2 forward; the
    full operator also needs the q-axis sweep and both integrate-back
    transposes, and an index slip in any of them is silent otherwise."""
    D = 3
    W = nq = D + 1
    rng = np.random.RandomState(11)
    ok = True

    def shape_of(axis, tr):
        if axis == 2 and not tr:  return (W, W, W), (W, W, W, nq), (W, W, nq)
        if axis == 2 and tr:      return (W, W, nq), (W, W, W, nq), (W, W, W)
        if axis == 1 and not tr:  return (W, W, nq), (W, W, nq),    (W, nq, nq)
        return (W, nq, nq), (W, W, nq), (W, W, nq)

    def reference(axis, tr, x, T, out):
        res = np.zeros(out)
        for p in range(W):
            for o1 in range(nq if (axis == 1 and not tr) else W - p):
                hi2 = (W - p - o1) if (axis == 2 and tr) else nq
                for o2 in range(hi2):
                    if axis == 2 and not tr:
                        res[p, o1, o2] = sum(x[p, o1, r] * T[p, o1, r, o2]
                                             for r in range(W - p - o1))
                    elif axis == 2:
                        res[p, o1, o2] = sum(x[p, o1, k] * T[p, o1, o2, k]
                                             for k in range(nq))
                    elif not tr:
                        res[p, o1, o2] = sum(T[p, q, o1] * x[p, q, o2]
                                             for q in range(W - p))
                    else:
                        res[p, o1, o2] = sum(T[p, o1, j] * x[p, j, o2]
                                             for j in range(nq))
        return res

    def ty(sh):
        return "tensor<" + "x".join(str(d) for d in sh) + "xf64>"

    for axis, tr in ((2, False), (2, True), (1, False), (1, True)):
        xs, Ts, os_ = shape_of(axis, tr)
        x = rng.uniform(-1, 1, xs)
        T = rng.uniform(-1, 1, Ts)
        expected = reference(axis, tr, x, T, os_)
        payload = f"""
func.func private @printMemrefF64(tensor<*xf64>)
func.func @main() {{
{cst("x", x)}
{cst("T", T)}
    %y = mir.simplex_contract %x, %T {{degree = {D} : i64, axis = {axis} : i64, transposed = {str(tr).lower()}}}
         : ({ty(xs)}, {ty(Ts)}) -> {ty(os_)}
{cst("exp", expected)}
    %scale = arith.constant 1.0e12 : f64
    %diff = mir.flux ins(%y, %exp) : ({ty(os_)}, {ty(os_)}) -> {ty(os_)} {{
    ^bb0(%a: f64, %b: f64):
      %d = arith.subf %a, %b : f64
      %ad = math.absf %d : f64
      %sd = arith.mulf %ad, %scale : f64
      mir.yield %sd : f64
    }}
    %pr = tensor.cast %diff : {ty(os_)} to tensor<*xf64>
    call @printMemrefF64(%pr) : (tensor<*xf64>) -> ()
    return
}}
"""
        got = run_pipeline(payload)
        scaled = np.max(np.abs(got)) if len(got) else 9e9
        good = len(got) == int(np.prod(os_)) and scaled < 1.0
        ok &= good
        label = f"axis={axis}{' transposed' if tr else '          '}"
        print(f"gate 7 ({label})           : max|err| = {scaled:.3e}e-12  {'PASS' if good else 'FAIL'}")
    return ok


def gate8():
    """THE FULL COLLAPSED-TET OPERATOR in high-level mir: y = B^T G_c B u on the
    Duffy tetrahedron, sum-factorized, executed and compared against a NumPy
    transcription of tet_galerkin.py (the validated reference).

    Six stages per gradient component. The p sweep is full-range (mir.contract);
    the q and r sweeps are ragged (mir.simplex_contract, axis 1 / 2, forward and
    transposed). This is the tet analogue of gate 5."""
    D = 3
    W = n = D + 1
    rng = np.random.RandomState(23)
    u  = rng.uniform(-1, 1, (W, W, W))
    A  = rng.uniform(-1, 1, (W, n));  Ad = rng.uniform(-1, 1, (W, n))
    B  = rng.uniform(-1, 1, (W, W, n)); Bd = rng.uniform(-1, 1, (W, W, n))
    C  = rng.uniform(-1, 1, (W, W, W, n)); Cd = rng.uniform(-1, 1, (W, W, W, n))
    G  = rng.uniform(-1, 1, (3, 3, n, n, n))

    # --- NumPy oracle: tet_galerkin.py, transcribed ---
    def sweep(Af, Bf, Cf):
        f1 = np.zeros((W, W, n)); f2 = np.zeros((W, n, n)); g = np.zeros((n, n, n))
        for p in range(W):
            for q in range(W - p):
                for k in range(n):
                    f1[p, q, k] = sum(u[p, q, r] * Cf[p, q, r, k] for r in range(W - p - q))
        for p in range(W):
            for j in range(n):
                for k in range(n):
                    f2[p, j, k] = sum(Bf[p, q, j] * f1[p, q, k] for q in range(W - p))
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    g[i, j, k] = sum(Af[p, i] * f2[p, j, k] for p in range(W))
        return g
    grad = [sweep(Ad, B, C), sweep(A, Bd, C), sweep(A, B, Cd)]
    flux = [sum(G[a, b] * grad[b] for b in range(3)) for a in range(3)]

    def tsweep(Af, Bf, Cf, fl):
        h1 = np.zeros((W, n, n)); h2 = np.zeros((W, W, n)); yy = np.zeros((W, W, W))
        for p in range(W):
            for j in range(n):
                for k in range(n):
                    h1[p, j, k] = sum(Af[p, i] * fl[i, j, k] for i in range(n))
        for p in range(W):
            for q in range(W - p):
                for k in range(n):
                    h2[p, q, k] = sum(Bf[p, q, j] * h1[p, j, k] for j in range(n))
        for p in range(W):
            for q in range(W - p):
                for r in range(W - p - q):
                    yy[p, q, r] = sum(Cf[p, q, r, k] * h2[p, q, k] for k in range(n))
        return yy
    expected = (tsweep(Ad, B, C, flux[0]) + tsweep(A, Bd, C, flux[1])
                + tsweep(A, B, Cd, flux[2]))

    T3 = f"tensor<{W}x{W}x{W}xf64>"
    TQ = f"tensor<{n}x{n}x{n}xf64>"
    TF1 = f"tensor<{W}x{W}x{n}xf64>"
    TF2 = f"tensor<{W}x{n}x{n}xf64>"
    TB = f"tensor<{W}x{W}x{n}xf64>"
    TC = f"tensor<{W}x{W}x{W}x{n}xf64>"
    TA = f"tensor<{W}x{n}xf64>"
    TAT = f"tensor<{n}x{W}xf64>"

    decls = [cst("u", u), cst("Am", A), cst("Adm", Ad), cst("AmT", A.T),
             cst("AdmT", Ad.T), cst("Bm", B), cst("Bdm", Bd),
             cst("Cm", C), cst("Cdm", Cd)]
    for a in range(3):
        for b in range(3):
            decls.append(cst(f"g{a}{b}", G[a, b]))

    body = []
    # forward: per component, ragged r sweep -> ragged q sweep -> full p sweep
    for c, (af, bf, cf) in enumerate([("AdmT", "Bm", "Cm"), ("AmT", "Bdm", "Cm"),
                                      ("AmT", "Bm", "Cdm")]):
        body.append(f"""    %f1_{c} = mir.simplex_contract %u, %{cf} {{degree = {D} : i64, axis = 2 : i64}}
         : ({T3}, {TC}) -> {TF1}
    %f2_{c} = mir.simplex_contract %f1_{c}, %{bf} {{degree = {D} : i64, axis = 1 : i64}}
         : ({TF1}, {TB}) -> {TF2}
    %gr{c} = mir.contract %f2_{c}, %{af} {{axis = 0 : i64}} : ({TF2}, {TAT}) -> {TQ}""")
    # pointwise flux: the 3x3 collapse metric maps gradient -> flux
    for a in range(3):
        body.append(f"""    %fx{a} = mir.flux ins(%gr0, %gr1, %gr2, %g{a}0, %g{a}1, %g{a}2)
         : ({TQ}, {TQ}, {TQ}, {TQ}, {TQ}, {TQ}) -> {TQ} {{
    ^bb0(%d0: f64, %d1: f64, %d2: f64, %m0: f64, %m1: f64, %m2: f64):
      %p0 = arith.mulf %m0, %d0 : f64
      %p1 = arith.mulf %m1, %d1 : f64
      %p2 = arith.mulf %m2, %d2 : f64
      %s0 = arith.addf %p0, %p1 : f64
      %s1 = arith.addf %s0, %p2 : f64
      mir.yield %s1 : f64
    }}""")
    # transpose: integrate the flux back to modal coefficients
    for c, (af, bf, cf) in enumerate([("Adm", "Bm", "Cm"), ("Am", "Bdm", "Cm"),
                                      ("Am", "Bm", "Cdm")]):
        body.append(f"""    %h1_{c} = mir.contract %fx{c}, %{af} {{axis = 0 : i64}} : ({TQ}, {TA}) -> {TF2}
    %h2_{c} = mir.simplex_contract %h1_{c}, %{bf} {{degree = {D} : i64, axis = 1 : i64, transposed = true}}
         : ({TF2}, {TB}) -> {TF1}
    %y{c} = mir.simplex_contract %h2_{c}, %{cf} {{degree = {D} : i64, axis = 2 : i64, transposed = true}}
         : ({TF1}, {TC}) -> {T3}""")
    body.append(f"""    %y = mir.flux ins(%y0, %y1, %y2) : ({T3}, {T3}, {T3}) -> {T3} {{
    ^bb0(%a0: f64, %a1: f64, %a2: f64):
      %t0 = arith.addf %a0, %a1 : f64
      %t1 = arith.addf %t0, %a2 : f64
      mir.yield %t1 : f64
    }}""")

    payload = f"""
func.func private @printMemrefF64(tensor<*xf64>)
func.func @main() {{
{chr(10).join(decls)}
{chr(10).join(body)}
{cst("exp", expected)}
    %scale = arith.constant 1.0e12 : f64
    %diff = mir.flux ins(%y, %exp) : ({T3}, {T3}) -> {T3} {{
    ^bb0(%a: f64, %b: f64):
      %d = arith.subf %a, %b : f64
      %ad = math.absf %d : f64
      %sd = arith.mulf %ad, %scale : f64
      mir.yield %sd : f64
    }}
    %pr = tensor.cast %diff : {T3} to tensor<*xf64>
    call @printMemrefF64(%pr) : (tensor<*xf64>) -> ()
    return
}}
"""
    got = run_pipeline(payload)
    scaled = np.max(np.abs(got)) if len(got) else 9e9
    ok = len(got) == W * W * W and scaled < 1.0
    print(f"gate 8 (FULL TET operator vs tet_galerkin): max|err| = {scaled:.3e}e-12  {'PASS' if ok else 'FAIL'}")
    return ok


if __name__ == "__main__":
    ok = gate1() & gate2() & gate3() & gate4() & gate5() & gate6() & gate7() & gate8()
    print("ALL PASS (mir lowering executes correctly)" if ok else "FAILED")
    sys.exit(0 if ok else 1)
