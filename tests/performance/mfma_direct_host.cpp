// Included after the extracted device functions by test_amd_mfma_direct.py.
// The emulated MFMA uses the measured gfx942 lane layout.
#include "ho_laplacian_host_check.hpp"
template<int Mode>
void run(const std::vector<double>& u, std::vector<double>& y, const std::vector<int>& dof, const std::vector<double>& D,
         const std::vector<double>& G, const std::vector<double>& weights, int elements)
{
    for(int block=0;block<(elements+1)/2;++block) {
        blockIdx.x=block;
        std::vector<std::thread> threads;
        for(int lane=0;lane<64;++lane) threads.emplace_back([&,lane]{
            threadIdx.x=lane;
            mars::fem::ho_laplacian_apply_mfma_plane<Mode>(u.data(),y.data(),dof.data(),D.data(),G.data(),weights.data(),elements);
        });
        for(auto& t:threads) t.join();
    }
}
int main() {
    bool ok=true;
    for(int n:{1,2,3,5}) for(bool constant:{false,true}) {
        ok=check<0>("MFMA plane",n,constant,run<0>)&&ok;
        ok=check<1>("MFMA plane",n,constant,run<1>)&&ok;
    }
    return ok?0:1;
}
