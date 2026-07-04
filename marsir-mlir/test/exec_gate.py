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
         "--expand-strided-metadata", "--finalize-memref-to-llvm",
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


if __name__ == "__main__":
    ok = gate1() & gate2() & gate3() & gate4() & gate5()
    print("ALL PASS (mir lowering executes correctly)" if ok else "FAILED")
    sys.exit(0 if ok else 1)
