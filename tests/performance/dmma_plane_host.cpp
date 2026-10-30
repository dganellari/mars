// Included after the extracted device functions by test_cuda_dmma_plane.py.
#include "ho_laplacian_host_check.hpp"
constexpr int warpsPerBlock = 4;
template<int Mode>
void run(const std::vector<double>& u, std::vector<double>& y, const std::vector<int>& dof, const std::vector<double>& D,
         const std::vector<double>& G, const std::vector<double>& weights, int elements)
{
    for(int block=0;block<(elements+warpsPerBlock-1)/warpsPerBlock;++block) {
        blockIdx.x=block;
        // The kernel synchronizes only within a warp, so the warps of a block can run one after another.
        for(int w=0;w<warpsPerBlock;++w) {
            std::vector<std::thread> threads;
            for(int lane=0;lane<32;++lane) threads.emplace_back([&,w,lane]{
                threadIdx.x=w*32+lane;
                mars::fem::ho_laplacian_apply_dmma_plane<Mode,warpsPerBlock>(u.data(),y.data(),dof.data(),D.data(),G.data(),weights.data(),elements);
            });
            for(auto& t:threads) t.join();
        }
    }
}
int main() {
    bool ok=true;
    for(int n:{1,2,3,5,9}) for(bool constant:{false,true}) {
        ok=check<0>("DMMA plane",n,constant,run<0>)&&ok;
        ok=check<1>("DMMA plane",n,constant,run<1>)&&ok;
    }
    return ok?0:1;
}
