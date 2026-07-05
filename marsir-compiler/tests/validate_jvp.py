#!/usr/bin/env python3
"""Local numeric proof of SYMBOLIC FORM DIFFERENTIATION (the Jacobian-action).

No GPU. Two checks, p=1..8, on a distorted hex:

  A. Linear self-adjoint Laplacian: the generated JVP  J(u).du  must equal the
     primal operator applied to du  (A.du)  to reduction order -- because a
     linear operator IS its own Jacobian. Confirms the JVP machinery reduces
     correctly in the linear case.

  B. Nonlinear demo flux (quadratic): the generated JVP must match a central
     finite-difference of the primal residual,
        J(u).du  ~=  (F(u+eps.du) - F(u-eps.du)) / (2 eps).
     This is the independent numerical proof that d flux / d x is correct for a
     genuinely nonlinear operator.

Run from marsir-compiler/:  python3 tests/validate_jvp.py
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
from marsir import parse_spec_file, synthesize          # noqa: E402
from marsir.backends import host_cpp                     # noqa: E402

DRIVER = r"""
#include "mars_cvfem_ho_apply.hpp"
#include "laplacian_apply_host.hpp"
#include "laplacian_jvp_host.hpp"
#include "nldemo_apply_host.hpp"
#include "nldemo_jvp_host.hpp"
#include <random>
#include <vector>
#include <cstdio>
#include <cmath>

using namespace mars::fem;
using namespace mars::marsir;

int main() {
    std::mt19937 rng(2024);
    std::uniform_real_distribution<double> pert(-0.15, 0.15), uni(-1.0, 1.0);
    int fails = 0;

    for (int P = 1; P <= 8; ++P) {
        auto op = buildHoCvfemOperators(P);
        double corners[8][3];
        const auto& S = hexCornerRef();
        for (int c = 0; c < 8; ++c)
            for (int a = 0; a < 3; ++a)
                corners[c][a] = double(S[c][a]) + pert(rng);
        auto G = computeElementMetric(op, corners);

        const int n = P + 1, n3 = n * n * n;
        std::vector<double> u(n3), du(n3);
        for (auto& x : u) x = uni(rng);
        for (auto& x : du) x = uni(rng);

        // --- A: linear Laplacian, J(u).du == A.du (reduction order) ---
        std::vector<double> Adu, Jdu_lin;
        laplacian_apply_host(op, G, du, Adu);
        laplacian_jvp_host(op, G, u, du, Jdu_lin);
        double mA = 0, nA = 0;
        for (int i = 0; i < n3; ++i) {
            mA = std::max(mA, std::fabs(Jdu_lin[i] - Adu[i]));
            nA = std::max(nA, std::fabs(Adu[i]));
        }
        double relA = mA / (nA > 0 ? nA : 1.0);

        // --- B: nonlinear demo, JVP vs central finite difference ---
        const double eps = 1e-6;
        std::vector<double> up(n3), um(n3), Fp, Fm, Jdu;
        for (int i = 0; i < n3; ++i) { up[i] = u[i] + eps * du[i]; um[i] = u[i] - eps * du[i]; }
        nldemo_apply_host(op, G, up, Fp);
        nldemo_apply_host(op, G, um, Fm);
        nldemo_jvp_host(op, G, u, du, Jdu);
        double mB = 0, nB = 0;
        for (int i = 0; i < n3; ++i) {
            double fd = (Fp[i] - Fm[i]) / (2.0 * eps);
            mB = std::max(mB, std::fabs(Jdu[i] - fd));
            nB = std::max(nB, std::fabs(fd));
        }
        double relB = mB / (nB > 0 ? nB : 1.0);

        bool ok = (relA < 1e-12) && (relB < 1e-6);
        printf("p=%d  n3=%4d  A[linear J==A]=%.3e  B[nonlinear J vs FD]=%.3e  %s\n",
               P, n3, relA, relB, ok ? "OK" : "FAIL");
        if (!ok) ++fails;
    }
    printf(fails ? "\n%d order(s) FAILED\n" : "\nALL PASS (symbolic Jacobian verified)\n", fails);
    return fails ? 1 : 0;
}
"""


def _emit(spec, primal_name, jvp_name):
    ea = synthesize(parse_spec_file(os.path.join(ROOT, "specs", spec)))
    with open(os.path.join(GEN, primal_name), "w") as f:
        f.write(host_cpp.emit(ea))
    with open(os.path.join(GEN, jvp_name), "w") as f:
        f.write(host_cpp.emit_jvp(ea))


def main():
    os.makedirs(GEN, exist_ok=True)
    _emit("laplacian.op", "laplacian_apply_host.hpp", "laplacian_jvp_host.hpp")
    _emit("nldemo.op", "nldemo_apply_host.hpp", "nldemo_jvp_host.hpp")

    drv = os.path.join(GEN, "validate_jvp_driver.cpp")
    with open(drv, "w") as f:
        f.write(DRIVER)

    exe = os.path.join(GEN, "validate_jvp_driver")
    cc = os.environ.get("CXX", "c++")
    cmd = [cc, "-std=c++17", "-O2", f"-I{FEM}", f"-I{GEN}", drv, "-o", exe]
    print("compile:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    print("run:", exe)
    sys.exit(subprocess.run([exe]).returncode)


if __name__ == "__main__":
    main()
