// Memory crossover: assembled-matrix (CSR) footprint vs matrix-free, swept over
// order p, for the high-order Galerkin operator. Host-only (uses HODofHandler,
// no CUDA) so it runs anywhere. This is the quantitative core of the trillion-DOF
// argument: at high p the assembled matrix is enormous (NNZ ~ (2p+1)^3 per row)
// while matrix-free stores only the DOF vectors + a per-element metric.
//
// We count the EXACT global CSR nonzeros from the topological DOF map (two DOFs
// couple iff they share an element), then compare bytes/DOF and extrapolate to
// 1e9 and 1e12 DOFs.

#include "backend/distributed/unstructured/fem/mars_ho_dof_handler.hpp"

#include <algorithm>
#include <array>
#include <cstdio>
#include <unordered_set>
#include <vector>

using namespace mars::fem;

template <int P>
void sweep(int E)
{
    constexpr int n  = P + 1;
    constexpr int N3 = n * n * n;

    auto cgid = [&](int cx, int cy, int cz) { return (cx * (E + 1) + cy) * (E + 1) + cz; };
    std::vector<std::array<int, 8>> ec;
    ec.reserve(size_t(E) * E * E);
    for (int ex = 0; ex < E; ++ex)
        for (int ey = 0; ey < E; ++ey)
            for (int ez = 0; ez < E; ++ez)
                ec.push_back({cgid(ex, ey, ez),     cgid(ex + 1, ey, ez),     cgid(ex + 1, ey + 1, ez),     cgid(ex, ey + 1, ez),
                              cgid(ex, ey, ez + 1), cgid(ex + 1, ey, ez + 1), cgid(ex + 1, ey + 1, ez + 1), cgid(ex, ey + 1, ez + 1)});

    HODofHandler dh;
    dh.build(ec, long(E + 1) * (E + 1) * (E + 1), P);
    const long   nDof = dh.numDof;
    const size_t nEl  = ec.size();

    // Exact NNZ: per-row set of coupled column DOFs (share an element).
    std::vector<std::unordered_set<int>> cols(nDof);
    for (size_t e = 0; e < nEl; ++e)
        for (int I = 0; I < N3; ++I)
        {
            auto& row = cols[dh.elemDof[e * N3 + I]];
            for (int J = 0; J < N3; ++J) row.insert(dh.elemDof[e * N3 + J]);
        }
    long long nnz = 0;
    for (long r = 0; r < nDof; ++r) nnz += (long long)cols[r].size();

    // Footprints. CSR: value(8B) + colIndex(4B) per nonzero (+ rowPtr, negligible).
    const double matBytesPerDof = double(nnz) * 12.0 / double(nDof);
    // Matrix-free: ~5 CG vectors (x,r,p,Ap,b) * 8B/DOF + per-element metric (6
    // doubles/elem, affine) amortized over the DOFs.
    const double mfBytesPerDof = 5.0 * 8.0 + 6.0 * 8.0 * double(nEl) / double(nDof);
    const double nnzPerRow     = double(nnz) / double(nDof);

    auto tb = [](double bytesPerDof, double dofs) { return bytesPerDof * dofs / 1.0e12; };  // → TB

    printf("p=%d E=%2d | nDof=%-9ld nnz=%-12lld nnz/row=%6.0f | "
           "matrix %8.1f B/DOF  matfree %5.1f B/DOF  (%6.0fx) | "
           "@1e12 DOF: matrix %8.1f TB  matfree %6.2f TB\n",
           P, E, nDof, nnz, nnzPerRow,
           matBytesPerDof, mfBytesPerDof, matBytesPerDof / mfBytesPerDof,
           tb(matBytesPerDof, 1.0e12), tb(mfBytesPerDof, 1.0e12));
}

int main()
{
    printf("High-order operator: assembled-CSR vs matrix-free memory footprint\n");
    printf("(structured cube, Galerkin Laplacian; nnz counted exactly from the DOF map)\n\n");
    sweep<1>(40);
    sweep<2>(24);
    sweep<4>(12);
    sweep<7>(6);
    printf("\nNote: matrix B/DOF is scale-invariant (interior-dominated); the @1e12 columns\n"
           "show why matrix-free is the enabling technology for trillion-DOF high order.\n");
    return 0;
}
