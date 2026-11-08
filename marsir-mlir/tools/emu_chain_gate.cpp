// CPU harness of test/chain_handoff_gate.py: one emulated warp runs a small kernel
// k(X, G, L, R, C, Y, Z) on 8x8 row-major matrices. The five inputs come from stdin,
// the outputs Y and Z go to stdout, NaN where the kernel writes nothing.
#include <cmath>
#include <cstdio>
#include <thread>
#include <vector>

#include "emu_warp_runtime.h"

using D2 = Desc<double, 2>;
extern "C" void _mlir_ciface_k(D2*, D2*, D2*, D2*, D2*, D2*, D2*);

int main()
{
    std::vector<double> m[7];
    for (int a = 0; a < 7; ++a) {
        m[a].assign(64, std::nan(""));
        if (a < 5)
            for (double& x : m[a])
                if (scanf("%lf", &x) != 1) return 2;
    }
    D2 d[7];
    for (int a = 0; a < 7; ++a) d[a] = D2{m[a].data(), m[a].data(), 0, {8, 8}, {8, 1}};
    Warp w;
    warp = &w;
    std::vector<std::thread> lanes;
    for (int L = 0; L < 32; ++L)
        lanes.emplace_back([&, L] {
            tlLane = L;
            _mlir_ciface_k(&d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6]);
        });
    for (auto& t : lanes) t.join();
    for (int a = 5; a < 7; ++a)
        for (double x : m[a]) printf("%.17g\n", x);
    return 0;
}
