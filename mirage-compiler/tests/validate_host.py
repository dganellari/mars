#!/usr/bin/env python3
"""Local numeric parity: generated host apply vs the hand-written oracle.

No GPU. Emits the host-C++ twin from specs/laplacian.op, compiles it with a
driver that also includes mars::fem::applyHoCvfemElement, and diffs the two on a
distorted hex for p=1..8. If they agree to machine precision, the IR has
captured the operator correctly -- the evidence that de-risks the Alps CUDA gate.

Run from the mirage-compiler/ directory:
    python3 tests/validate_host.py
"""

import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FEM = os.path.abspath(os.path.join(
    ROOT, "..", "backend", "distributed", "unstructured", "fem"))
GEN = os.path.join(ROOT, "generated")

sys.path.insert(0, ROOT)
from mirage import parse_spec_file, synthesize          # noqa: E402
from mirage.backends import host_cpp                     # noqa: E402

DRIVER = r"""
#include "mars_cvfem_ho_apply.hpp"
#include "laplacian_apply_host.hpp"
#include <random>
#include <vector>
#include <cstdio>
#include <cmath>

using namespace mars::fem;

int main() {
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> pert(-0.15, 0.15), uni(-1.0, 1.0);
    int fails = 0;
    const double tol = 1e-12;

    for (int P = 1; P <= 8; ++P) {
        auto op = buildHoCvfemOperators(P);
        // Distorted (non-orthogonal) hex so the metric has real cross terms.
        double corners[8][3];
        const auto& S = hexCornerRef();
        for (int c = 0; c < 8; ++c)
            for (int a = 0; a < 3; ++a)
                corners[c][a] = double(S[c][a]) + pert(rng);
        auto G = computeElementMetric(op, corners);

        const int n = P + 1, n3 = n * n * n;
        std::vector<double> u(n3), y_ref, y_gen;
        for (auto& x : u) x = uni(rng);

        applyHoCvfemElement(op, G, u, y_ref);
        mars::mirage::laplacian_apply_host(op, G, u, y_gen);

        double md = 0, refn = 0;
        for (int i = 0; i < n3; ++i) {
            md = std::max(md, std::fabs(y_gen[i] - y_ref[i]));
            refn = std::max(refn, std::fabs(y_ref[i]));
        }
        double rel = md / (refn > 0 ? refn : 1.0);

        // A.1 constant-nullspace invariant: u = const -> K u = 0.
        std::vector<double> uc(n3, 1.0), yc;
        mars::mirage::laplacian_apply_host(op, G, uc, yc);
        double null = 0;
        for (double v : yc) null = std::max(null, std::fabs(v));

        bool ok = (rel < tol) && (null < 1e-10);
        printf("p=%d  n3=%4d  parity_rel=%.3e  nullspace=%.3e  %s\n",
               P, n3, rel, null, ok ? "OK" : "FAIL");
        if (!ok) ++fails;
    }
    printf(fails ? "\n%d order(s) FAILED\n" : "\nALL PASS (generated == oracle)\n", fails);
    return fails ? 1 : 0;
}
"""


def main():
    os.makedirs(GEN, exist_ok=True)
    # (Re)generate the host twin from the spec.
    op = synthesize(parse_spec_file(os.path.join(ROOT, "specs", "laplacian.op")))
    with open(os.path.join(GEN, "laplacian_apply_host.hpp"), "w") as f:
        f.write(host_cpp.emit(op))

    drv = os.path.join(GEN, "validate_driver.cpp")
    with open(drv, "w") as f:
        f.write(DRIVER)

    exe = os.path.join(GEN, "validate_driver")
    cc = os.environ.get("CXX", "c++")
    compile_cmd = [cc, "-std=c++17", "-O2", f"-I{FEM}", f"-I{GEN}", drv, "-o", exe]
    print("compile:", " ".join(compile_cmd))
    subprocess.run(compile_cmd, check=True)
    print("run:", exe)
    r = subprocess.run([exe])
    sys.exit(r.returncode)


if __name__ == "__main__":
    main()
