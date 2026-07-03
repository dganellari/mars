"""Tet collapsed-coordinate (Duffy/PKD) Galerkin backend -- host emitter.

Emits the sum-factorized tet Laplacian apply  y = B^T G_c B u  on the collapsed
tetrahedron (Karniadakis-Sherwin PKD construction), the device twin target being
rung 9 of the sibling reference (see marsir-compiler/TET_SUMFAC_HANDOFF.md):

  fwd  : modal u -> 3 collapsed-gradient components at the Gauss-Jacobi grid, via
         three ragged 3-stage sweeps (contract r->q->p), one factor differentiated
         per component: (Ad,B,C)->d/da, (A,Bd,C)->d/db, (A,B,Cd)->d/dc.
  flux : per-quad 3x3 metric contraction  flux = G_c . grad  (the Galerkin action).
  bwd  : transpose sweeps (contract i->j->k) integrate the flux back to modal y.

Ragged bounds p+q+r<=D are the simplex analogue of the hex rectangular i,j,k in
[0,p]. A/B/C are the 1D PKD factor tables at the quad nodes; G_c the collapse
metric. Both are built by the reference library (jacobi/gauss_jacobi/pkd/tet_quad)
-- geometry/basis infrastructure, like buildHoCvfemOperators for hex. This emitter
generates only the OPERATOR APPLY.

Structurally a faithful port of the sibling's test_sumfac_op.cpp Asf/fwd/bwd, so
the generated apply == the dense per-mode oracle to ~1e-14 (validate_tet.py).
"""

from .common import AUTOGEN_BANNER


def emit(name="tet_laplacian"):
    banner = AUTOGEN_BANNER.format(backend="tet-host", name=name,
                                   flux="B^T G_c B (collapsed Galerkin Laplacian)")
    fn = f"{name}_apply_host"
    return f"""#pragma once
{banner}
#include <vector>
#include <array>
#include <cstddef>

namespace mars {{
namespace marsir {{

// y = B^T G_c B u on the collapsed (Duffy) tetrahedron, sum-factorized.
// Inputs (built by the PKD/Gauss-Jacobi reference library, order D):
//   A[p*n+i], Ad, B[(p*W+q)*n+j], Bd, C[((p*W+q)*W+r)*n+k], Cd  (W=n=D+1)
//   Gc[(i*n+j)*n+k][0..8]  : 3x3 collapse metric per quad point
//   u : modal coeffs, layout ((p*W+q)*W+r), ragged p+q+r<=D (others 0)
inline void {fn}(int D,
                 const std::vector<double>& A, const std::vector<double>& Ad,
                 const std::vector<double>& B, const std::vector<double>& Bd,
                 const std::vector<double>& C, const std::vector<double>& Cd,
                 const std::vector<std::array<double, 9>>& Gc,
                 const std::vector<double>& u, std::vector<double>& y)
{{
    const int n = D + 1, W = D + 1, NQ = n * n * n;
    std::vector<std::array<double, 3>> grad(NQ, {{0, 0, 0}}), flux(NQ);

    // --- forward: modal u -> collapsed gradient (3 components) at the quad grid ---
    auto sweep = [&](const double* Af, const double* Bf, const double* Cf, int comp) {{
        std::vector<double> f1((size_t)W * W * n, 0.0), f2((size_t)W * n * n, 0.0);
        for (int p = 0; p <= D; ++p) for (int q = 0; q <= D - p; ++q) for (int k = 0; k < n; ++k) {{
            double v = 0;
            for (int r = 0; r <= D - p - q; ++r) v += u[((size_t)p * W + q) * W + r] * Cf[(((size_t)p * W + q) * W + r) * n + k];
            f1[((size_t)p * W + q) * n + k] = v;
        }}
        for (int p = 0; p <= D; ++p) for (int j = 0; j < n; ++j) for (int k = 0; k < n; ++k) {{
            double v = 0;
            for (int q = 0; q <= D - p; ++q) v += Bf[((size_t)p * W + q) * n + j] * f1[((size_t)p * W + q) * n + k];
            f2[((size_t)p * n + j) * n + k] = v;
        }}
        for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) for (int k = 0; k < n; ++k) {{
            double v = 0;
            for (int p = 0; p <= D; ++p) v += Af[(size_t)p * n + i] * f2[((size_t)p * n + j) * n + k];
            grad[((size_t)i * n + j) * n + k][comp] = v;
        }}
    }};
    sweep(Ad.data(), B.data(), C.data(), 0);
    sweep(A.data(), Bd.data(), C.data(), 1);
    sweep(A.data(), B.data(), Cd.data(), 2);

    // --- pointwise flux: collapsed metric maps gradient -> flux (Galerkin action) ---
    for (int x = 0; x < NQ; ++x) {{
        const auto& G = Gc[x];
        const auto& g = grad[x];
        for (int a = 0; a < 3; ++a) flux[x][a] = G[a * 3 + 0] * g[0] + G[a * 3 + 1] * g[1] + G[a * 3 + 2] * g[2];
    }}

    // --- transpose: integrate collapsed flux back to modal y ---
    y.assign((size_t)W * W * W, 0.0);
    auto tsweep = [&](const double* Af, const double* Bf, const double* Cf, int comp) {{
        std::vector<double> h1((size_t)W * n * n, 0.0), h2((size_t)W * W * n, 0.0);
        for (int p = 0; p <= D; ++p) for (int j = 0; j < n; ++j) for (int k = 0; k < n; ++k) {{
            double v = 0;
            for (int i = 0; i < n; ++i) v += Af[(size_t)p * n + i] * flux[((size_t)i * n + j) * n + k][comp];
            h1[((size_t)p * n + j) * n + k] = v;
        }}
        for (int p = 0; p <= D; ++p) for (int q = 0; q <= D - p; ++q) for (int k = 0; k < n; ++k) {{
            double v = 0;
            for (int j = 0; j < n; ++j) v += Bf[((size_t)p * W + q) * n + j] * h1[((size_t)p * n + j) * n + k];
            h2[((size_t)p * W + q) * n + k] = v;
        }}
        for (int p = 0; p <= D; ++p) for (int q = 0; q <= D - p; ++q) for (int r = 0; r <= D - p - q; ++r) {{
            double v = 0;
            for (int k = 0; k < n; ++k) v += Cf[(((size_t)p * W + q) * W + r) * n + k] * h2[((size_t)p * W + q) * n + k];
            y[((size_t)p * W + q) * W + r] += v;
        }}
    }};
    tsweep(Ad.data(), B.data(), C.data(), 0);
    tsweep(A.data(), Bd.data(), C.data(), 1);
    tsweep(A.data(), B.data(), Cd.data(), 2);
}}

}}  // namespace marsir
}}  // namespace mars
"""
