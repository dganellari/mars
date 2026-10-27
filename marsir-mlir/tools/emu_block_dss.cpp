// CPU gate for the pass-generated block-structured DSS (block_dss_gpu_pipeline.py
// --emulate): the two kernels the GPU gets (block_sum, then block_skeleton), run with
// one thread per lane on the rotated-block test mesh of tests/cellwise_block_mesh.hpp,
// compared BIT FOR BIT with the hand bodies of MARS (block_node_sum and
// block_skeleton_node in backend/distributed/unstructured/solvers/
// mars_cellwise_block_dss.hpp, run on the host in launch order) and, to rounding, with
// the assembled sums of the mesh. Every copy must be written (finite), and a
// corrupted copy of the output must fail the comparison. One block per 8 elements of
// 32 items and, separately, 3 blocks striding over all items.
//
// Defines: MIR_BP (degree), MIR_BB (block size), MIR_BM (blocks per axis).
// Build + run: python3 test/block_dss_gpu_pipeline.py p B --emulate [--emu-tsan]
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <thread>
#include <vector>

#include "emu_warp_runtime.h"
#include "tests/cellwise_block_mesh.hpp"

#define MIR_CAT2(a, b) a##b
#define MIR_CAT(a, b) MIR_CAT2(a, b)
#define MIR_TAG MIR_CAT(MIR_CAT(MIR_CAT(p, MIR_BP), _b), MIR_BB)
#define MIR_SUM MIR_CAT(_mlir_ciface_block_sum_, MIR_TAG)
#define MIR_SKEL MIR_CAT(_mlir_ciface_block_skeleton_, MIR_TAG)
#define MIR_SUM_THREADS MIR_CAT(MIR_CAT(block_sum_, MIR_TAG), __threads)
#define MIR_SKEL_THREADS MIR_CAT(MIR_CAT(block_skeleton_, MIR_TAG), __threads)

using D64 = Desc<double>;
using D32 = Desc<int32_t>;
using D8 = Desc<int8_t>;
extern "C" void MIR_SUM(D64*, D64*);
extern "C" void MIR_SKEL(D64*, D32*, D8*, D32*, D32*, D32*, D32*);
extern "C" int64_t MIR_SUM_THREADS();
extern "C" int64_t MIR_SKEL_THREADS();

namespace cw = mars::cellwise;
constexpr int P = MIR_BP, B = MIR_BB;

namespace {

template <typename T, typename U>
Desc<T> desc(std::vector<U>& v)
{
    T* p = reinterpret_cast<T*>(v.data());
    return {p, p, 0, {(int64_t)v.size()}, {1}};
}

template <typename F>
void launch(int64_t threads, long long grid, F&& kernel)
{
    Warp w;
    w.bar.threads = (int)threads;
    warp = &w;
    for (long long b = 0; b < grid; ++b) {
        std::vector<std::thread> lanes;
        for (int L = 0; L < threads; ++L)
            lanes.emplace_back([&, L, b] {
                tlLane = L;
                tlBlock = b;
                tlGrid = grid;
                kernel();
            });
        for (auto& t : lanes) t.join();
    }
}

}  // namespace

int main()
{
    auto m = block_test::BlockMesh<P, B>::build(MIR_BM, 2026);
    const cw::BlockSkeleton sk = m.skeleton();

    // The hand bodies in launch order.
    std::vector<double> hand(m.N, std::nan(""));
    const long long items = m.blocks * cw::block_nodes_per_block<P, B>, skel = cw::block_skeleton_threads<P, B>(sk);
    for (long long t = 0; t < items; ++t) cw::block_node_sum<P, B>(m.in.data(), hand.data(), t);
    for (long long t = 0; t < skel; ++t) cw::block_skeleton_node<P, B>(hand.data(), sk, t);

    int fails = 0;
    const int64_t per_block = MIR_SUM_THREADS() / 32;   // elements of 32 items per block
    for (const bool stride : {false, true}) {
        std::vector<double> out(m.N, std::nan(""));
        D64 din = desc<double>(m.in), dout = desc<double>(out);
        D32 dfn = desc<int32_t>(m.face_nbr), deo = desc<int32_t>(m.edge_off), dee = desc<int32_t>(m.edge_ent),
            dvo = desc<int32_t>(m.vert_off), dve = desc<int32_t>(m.vert_ent);
        D8 dfc = desc<int8_t>(m.face_code);
        auto grid = [&](long long work) { return stride ? 3 : (work + 32 * per_block - 1) / (32 * per_block); };
        launch(MIR_SUM_THREADS(), grid(items), [&] { MIR_SUM(&din, &dout); });
        launch(MIR_SKEL_THREADS(), grid(skel), [&] { MIR_SKEL(&dout, &dfn, &dfc, &deo, &dee, &dvo, &dve); });

        long long differ = 0;
        for (long long t = 0; t < m.N; ++t) differ += std::memcmp(&out[t], &hand[t], sizeof(double)) != 0;
        const double err = m.error(out);
        std::vector<double> bad = out;
        bad[m.N / 3] = std::nan("");
        const bool rejected = !(m.error(bad) <= 1e-13);
        const bool ok = differ == 0 && err <= 1e-13 && rejected;
        printf("BLOCK DSS WARP EMULATION (p=%d B=%d, %d^3 rotated blocks, %lld copies, %s): %lld copies differ from "
               "the hand kernels' bits, max|err| vs assembled sums %.2e, corrupted copy rejected: %s  %s\n",
               P, B, MIR_BM, m.N, stride ? "3 blocks striding" : "one pass of blocks", differ, err,
               rejected ? "yes" : "no", ok ? "PASS" : "FAIL");
        fails += !ok;
    }
    return fails;
}
