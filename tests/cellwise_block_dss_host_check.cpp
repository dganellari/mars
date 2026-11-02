// Host check of the block-structured DSS bodies (backend/distributed/unstructured/
// solvers/mars_cellwise_block_dss.hpp) on the rotated-block test mesh of
// tests/cellwise_block_mesh.hpp. Each block stage (the passes, or one thread per block
// node) and the skeleton, run thread by thread in launch order, must give every copy
// the sum over all copies of its node. A non-finite output fails, and a corrupted copy
// of the output must fail the same comparison.
//
// Build and run from the MARS root:
//   c++ -std=c++20 -O2 -I. tests/cellwise_block_dss_host_check.cpp -o /tmp/block_dss_check && /tmp/block_dss_check
#include <cstdio>

#include "tests/cellwise_block_mesh.hpp"

namespace cw = mars::cellwise;

namespace {

template <int P, int B>
int check(int M, unsigned seed)
{
    const auto m = block_test::BlockMesh<P, B>::build(M, seed);
    const cw::BlockSkeleton sk = m.skeleton();
    std::mt19937 rng(seed + 1000);
    int fails = 0;
    for (const bool nodes : {false, true}) {
        std::vector<double> out(m.N, std::nan(""));
        if (nodes) {
            const long long items = m.blocks * cw::block_nodes_per_block<P, B>;
            for (long long t = 0; t < items; ++t) cw::block_node_sum<P, B>(m.in.data(), out.data(), t);
        } else {
            for (long long t = 0; t < m.N; ++t) cw::block_pass_z<P, B>(m.in.data(), out.data(), t);
            if constexpr (B > 1) {
                const long long pairs = m.blocks * cw::block_pairs_per_block<P, B>;
                for (long long t = 0; t < pairs; ++t) cw::block_pass_xy<P, B, 0>(out.data(), t);
                for (long long t = 0; t < pairs; ++t) cw::block_pass_xy<P, B, 1>(out.data(), t);
            }
        }
        const long long threads = cw::block_skeleton_threads<P, B>(sk);
        for (long long t = 0; t < threads; ++t) cw::block_skeleton_node<P, B>(out.data(), sk, t);
        const double err = m.error(out);
        std::vector<double> bad = out;
        bad[rng() % m.N] = std::nan("");
        const bool rejected = !(m.error(bad) <= 1e-13);
        const bool ok = err <= 1e-13 && rejected;
        printf("block DSS p=%d B=%d (%s), %d^3 rotated blocks (%lld copies): max|err| = %.2e, corrupted copy "
               "rejected: %s  %s\n",
               P, B, nodes ? "block nodes" : "passes", M, m.N, err, rejected ? "yes" : "no", ok ? "PASS" : "FAIL");
        fails += !ok;
    }
    return fails;
}

}  // namespace

int main()
{
    int fails = 0;
    fails += check<1, 1>(3, 1);
    fails += check<1, 2>(3, 2);
    fails += check<1, 4>(2, 3);
    fails += check<2, 1>(3, 4);
    fails += check<2, 3>(3, 5);
    fails += check<3, 2>(3, 6);
    fails += check<3, 4>(2, 7);
    fails += check<7, 2>(2, 8);
    return fails;
}
