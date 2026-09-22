#pragma once
#include <limits>

#if defined(__CUDACC__) || defined(__HIPCC__)
#define MARS_NODE_HD __host__ __device__
#else
#define MARS_NODE_HD
#endif

namespace mars::segregated {

MARS_NODE_HD inline bool finite_coefficient(double value)
{
    constexpr double bound = std::numeric_limits<double>::max();
    return value >= -bound && value <= bound;
}

struct SteadyMomentumNode {
    double density, volume, pseudo_dt, mass_divergence;
    double velocity[3], pressure_gradient[3], force[3], source[3];
};

// Fixed frame: physical pressure and forces per volume; mass_divergence is kg/s.
MARS_NODE_HD inline void steady_momentum_node(const SteadyMomentumNode& in,
                                             double* lhs, double* rhs)
{
    for (int i = 0; i < 9; ++i) lhs[i] = 0;
    for (int i = 0; i < 3; ++i) {
        lhs[3*i+i] = (in.mass_divergence < 0 ? -in.mass_divergence : 0)
                     + in.density*in.volume/in.pseudo_dt;
        rhs[i] = in.mass_divergence*in.velocity[i] - in.pressure_gradient[i]*in.volume
                 + (in.force[i]+in.source[i])*in.volume;
    }
}

// Increment form: change only the scalar diagonal, never add an absolute-u RHS.
MARS_NODE_HD inline void relax_momentum_diagonal(double* block, double alpha)
{
    const double inverse = 1/alpha;
    for (int i = 0; i < 3; ++i) block[3*i+i] *= inverse;
}

// Caller selects the union of eligible boundary nodes, once per owned node.
MARS_NODE_HD inline void relax_boundary_rhs(double* rhs, double factor)
{
    for (int i = 0; i < 3; ++i) rhs[i] *= factor;
}

// Row blocks are node-major 3x3. The entire diagonal-node block is excluded
// from SIMPLEC's neighbor sum; cross-component entries never contribute.
MARS_NODE_HD inline bool momentum_influence(double volume, const double* row_blocks,
                                            int blocks, int diagonal_block, bool consistent,
                                            double* d, double* d_tilde)
{
    constexpr double small = 2.2204460492503131e-16; // pinned reference SMALL (FP64 epsilon)
    if (!(volume > 0) || !finite_coefficient(volume) || blocks <= 0
        || diagonal_block < 0 || diagonal_block >= blocks || !row_blocks) return false;
    for (int i = 0; i < 3; ++i) {
        const double diagonal = row_blocks[9*diagonal_block+3*i+i];
        const double denominator = diagonal+small;
        if (denominator == 0 || !finite_coefficient(denominator)) return false;
        d[i] = volume/denominator;
        if (!finite_coefficient(d[i])) return false;
        d_tilde[i] = 0;
        if (consistent) {
            double neighbors = 0;
            for (int j = 0; j < blocks; ++j)
                if (j != diagonal_block) neighbors += row_blocks[9*j+3*i+i];
            const double denominator_c = diagonal+neighbors+small;
            if (denominator_c == 0 || !finite_coefficient(denominator_c)) return false;
            d_tilde[i] = volume/denominator_c;
            if (!finite_coefficient(d_tilde[i])) return false;
        }
    }
    return true;
}
} // namespace mars::segregated
#undef MARS_NODE_HD
