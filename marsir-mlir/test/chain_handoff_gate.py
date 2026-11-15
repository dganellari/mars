#!/usr/bin/env python3
# Gate for --mir-chain-contracts=direct-handoff=true: small kernels on 8x8 matrices,
# lowered with and without the option, run on the CPU warp emulator
# (tools/emu_chain_gate.cpp) and compared against numpy. With the option, the
# eligible cases must also be free of shuffles and the others must keep the shuffle
# path. F = X * G stands for geometry applied between contractions.
# Run from marsir-mlir/: python3 test/chain_handoff_gate.py

import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
LLVM = "/opt/homebrew/opt/llvm"
MIROPT = os.path.join(ROOT, "build", "tools", "mir-opt", "mir-opt")
MLIROPT = os.path.join(LLVM, "bin", "mlir-opt")
WORK = os.path.join(ROOT, "build", "emu_chain")

HEAD = """#mk = affine_map<(d0, d1, d2) -> (d0, d2)>
#kn = affine_map<(d0, d1, d2) -> (d2, d1)>
#nk = affine_map<(d0, d1, d2) -> (d1, d2)>
#mn = affine_map<(d0, d1, d2) -> (d0, d1)>
module {
 func.func @k(%X: memref<8x8xf64>, %G: memref<8x8xf64>, %L: memref<8x8xf64>, %R: memref<8x8xf64>,
              %C: memref<8x8xf64>, %Y: memref<8x8xf64>, %Z: memref<8x8xf64>) {
  %c0 = arith.constant 0 : index
  %f0 = arith.constant 0.0 : f64
  %z = arith.constant dense<0.0> : vector<8x8xf64>
"""


def read(v, mem, mt="memref<8x8xf64>", t="8x8"):
    return (f"  %{v} = vector.transfer_read %{mem}[%c0, %c0], %f0 {{in_bounds = [true, true]}} "
            f": {mt}, vector<{t}xf64>\n")


def write(v, mem, mt="memref<8x8xf64>", t="8x8"):
    return (f"  vector.transfer_write %{v}, %{mem}[%c0, %c0] {{in_bounds = [true, true]}} "
            f": vector<{t}xf64>, {mt}\n")


# A short tile pads its missing rows by reading and writing past the END OF THE
# MEMREF, so a 7-row value lives in a 7-row view.
SUB7 = "memref<7x8xf64, strided<[8, 1]>>"


def view7(v, mem):
    return f"  %{v} = memref.subview %{mem}[0, 0] [7, 8] [1, 1] : memref<8x8xf64> to {SUB7}\n"


def contract(res, lhs, lmap, rhs, rmap, acc, lt="8x8", rt="8x8", ot="8x8"):
    return (f"  %{res} = vector.contract {{indexing_maps = [#{lmap}, #{rmap}, #mn], "
            f"iterator_types = [\"parallel\", \"parallel\", \"reduction\"], kind = #vector.kind<add>}} "
            f"%{lhs}, %{rhs}, %{acc} : vector<{lt}xf64>, vector<{rt}xf64> into vector<{ot}xf64>\n")


def mul(res, a, b, t="8x8"):
    return f"  %{res} = arith.mulf %{a}, %{b} : vector<{t}xf64>\n"


LEAVES = "".join(read(v, m) for v, m in [("x", "X"), ("g", "G"), ("l", "L"), ("r", "R"), ("acc", "C")])
F = LEAVES + mul("f", "x", "g")
NAN = np.full((8, 8), np.nan)

# name -> (body, reference(X, G, L, R, C) -> (Y, Z), shuffle-free with the option)
CASES = {
    "pair, transposed leaf, memory accumulator": (
        F + contract("t", "r", "nk", "f", "mk", "z") + contract("y", "l", "mk", "t", "kn", "acc") + write("y", "Y"),
        lambda X, G, L, R, C: (L @ ((X * G) @ R.T) + C, NAN), True),
    "pair, standard leaf, geometry after": (
        F + contract("t", "f", "mk", "r", "kn", "z") + contract("y", "l", "mk", "t", "kn", "z")
        + mul("yg", "y", "g") + write("yg", "Y"),
        lambda X, G, L, R, C: ((L @ ((X * G) @ R)) * G, NAN), True),
    "pair from two leaves": (
        LEAVES + contract("t", "x", "mk", "r", "kn", "z") + contract("y", "l", "mk", "t", "kn", "acc") + write("y", "Y"),
        lambda X, G, L, R, C: (L @ (X @ R) + C, NAN), True),
    "A role, both leaf layouts, stored": (
        F + contract("t", "r", "nk", "f", "mk", "z") + write("t", "Y")
        + contract("u", "f", "mk", "r", "kn", "z") + write("u", "Z"),
        lambda X, G, L, R, C: ((X * G) @ R.T, (X * G) @ R), True),
    "transposed B role": (
        F + contract("y", "l", "mk", "f", "nk", "z") + write("y", "Y"),
        lambda X, G, L, R, C: (L @ (X * G).T, NAN), True),
    "two fragment operands": (
        F + mul("h", "l", "g") + contract("y", "h", "nk", "f", "mk", "z") + write("y", "Y"),
        lambda X, G, L, R, C: ((X * G) @ (L * G).T, NAN), True),
    "two consumers (shuffle fallback)": (
        F + contract("t", "r", "nk", "f", "mk", "z") + contract("y", "l", "mk", "t", "kn", "z")
        + write("y", "Y") + write("t", "Z"),
        lambda X, G, L, R, C: (L @ ((X * G) @ R.T), (X * G) @ R.T), False),
    "producer accumulator (shuffle fallback)": (
        F + contract("t", "r", "nk", "f", "mk", "acc") + contract("y", "l", "mk", "t", "kn", "z") + write("y", "Y"),
        lambda X, G, L, R, C: (L @ ((X * G) @ R.T + C), NAN), False),
    "short m (declined)": (
        LEAVES + view7("X7", "X") + view7("G7", "G") + view7("Y7", "Y")
        + read("x7", "X7", SUB7, "7x8") + read("g7", "G7", SUB7, "7x8") + mul("f7", "x7", "g7", "7x8")
        + "  %z7 = arith.constant dense<0.0> : vector<7x8xf64>\n"
        + contract("t", "r", "nk", "f7", "mk", "z7", rt="7x8", ot="7x8") + write("t", "Y7", SUB7, "7x8"),
        lambda X, G, L, R, C: (np.vstack([((X * G) @ R.T)[:7], np.full((1, 8), np.nan)]), NAN), False),
}


def run(cmd, data=""):
    p = subprocess.run(cmd, input=data, capture_output=True, text=True)
    if p.returncode:
        sys.exit(f"FAILED: {' '.join(cmd[:3])}\n{p.stderr[:2000]}")
    return p.stdout


def build(body, direct, tag):
    """The kernel through the chain pass and the warp emulation, linked with the harness."""
    ir = HEAD + body + "  return\n }\n}\n"
    ir = run([MIROPT, "-", "--mir-chain-contracts" + ("=direct-handoff=true" if direct else ""),
              "--canonicalize", "--cse"], ir)
    shuffles = ir.count("gpu.shuffle")
    ir = run([MIROPT, "-", "--mir-emulate-warp"], ir)
    ir = run([MLIROPT, "-", "--convert-vector-to-scf", "--canonicalize", "--convert-scf-to-cf",
              "--expand-strided-metadata", "--lower-affine", "--convert-vector-to-llvm",
              "--convert-arith-to-llvm", "--convert-index-to-llvm", "--convert-cf-to-llvm",
              "--finalize-memref-to-llvm", "--convert-func-to-llvm", "--reconcile-unrealized-casts"], ir)
    ll = run([os.path.join(LLVM, "bin", "mlir-translate"), "--mlir-to-llvmir"], ir)
    base = os.path.join(WORK, tag)
    with open(base + ".ll", "w") as f:
        f.write(ll)
    run([os.path.join(LLVM, "bin", "clang"), "-O1", "-c", base + ".ll", "-o", base + ".o"])
    run(["/usr/bin/clang++", "-std=c++20", "-O1", os.path.join(ROOT, "tools", "emu_chain_gate.cpp"),
         base + ".o", "-o", base])
    return base, shuffles


def main():
    os.makedirs(WORK, exist_ok=True)
    rng = np.random.default_rng(7)
    ins = [rng.uniform(-1, 1, (8, 8)) for _ in range(5)]
    data = "\n".join(f"{v:.17g}" for m in ins for v in m.ravel())
    failed = 0
    for i, (name, (body, ref, free)) in enumerate(CASES.items()):
        want = ref(*ins)
        for direct in (False, True):
            exe, shuffles = build(body, direct, f"case{i}_{int(direct)}")
            out = np.array([float(v) for v in run([exe], data).split()]).reshape(2, 8, 8)
            same_nan = all((np.isnan(o) == np.isnan(w)).all() for o, w in zip(out, want))
            err = max(np.nanmax(np.abs(o - w)) if not np.isnan(w).all() else 0.0 for o, w in zip(out, want))
            ok = same_nan and err < 1e-12
            if direct:
                ok = ok and (shuffles == 0) == free
            print(f"  {name:44} {'direct ' if direct else 'shuffle'}  {shuffles:2} shuffles  "
                  f"err {err:.1e}  {'PASS' if ok else 'FAIL'}")
            failed += not ok
    print("CHAIN HAND-OFF GATE: " + ("PASS" if not failed else f"{failed} FAILED"))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
