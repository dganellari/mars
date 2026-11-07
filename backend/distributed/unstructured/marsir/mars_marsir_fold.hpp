#pragma once
// The folded form of a mirror-symmetric 1D matrix M (r x c, row-major, M[r-1-i][c-1-j] =
// s M[i][j]) that kernels from marsir-mlir --mir-even-odd read in place of M:
//   [A (h_r x h_c), S (h_r x h_c), M[i][c/2] (i < h_r, c odd), M[r/2][j] (j < h_c, r odd),
//    M[r/2][c/2] (both odd)],  h_r = r / 2, h_c = c / 2,
// with A and S the half sums and differences of M[i][j] and M[i][c-1-j]. The kernel
// publishes {r, c, s} as <matrix symbol>__fold. Used by the launcher (GalerkinPtx) and by
// the CPU gate (marsir-mlir/tools/emu_bp.cpp).

#if defined(__CUDACC__) || defined(__HIPCC__)
#define MARS_MARSIR_HD __host__ __device__
#else
#define MARS_MARSIR_HD
#endif

namespace mars {
namespace marsir {

MARS_MARSIR_HD inline long long folded_size(long long r, long long c)
{
    const long long hr = r / 2, hc = c / 2, mr = r % 2, mc = c % 2;
    return 2 * hr * hc + hr * mc + mr * hc + mr * mc;
}

template <typename Real>
MARS_MARSIR_HD inline Real folded_entry(const Real* M, long long r, long long c, long long k)
{
    const long long hr = r / 2, hc = c / 2, mr = r % 2, mc = c % 2;
    const long long block = hr * hc, col = 2 * block, row = col + hr * mc, mid = row + mr * hc;
    if (k < 2 * block) {
        const long long i = (k % block) / hc, j = (k % block) % hc;
        const Real a = M[i * c + j], b = M[i * c + c - 1 - j];
        return k < block ? (a + b) / 2 : (a - b) / 2;
    }
    if (k < row) return M[(k - col) * c + hc];
    if (k < mid) return M[hr * c + (k - row)];
    return M[hr * c + hc];
}

}  // namespace marsir
}  // namespace mars
