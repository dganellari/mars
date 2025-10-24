#ifdef MARS_ENABLE_CUDA
#include <cstone/traversal/collisions_gpu.cu>  // For markMacsGpu
// Include the header/file for findHalosGpu (likely <cstone/halos/exchange_halos_gpu.cuh> or similar)
#include <cstone/halos/exchange_halos_gpu.cuh>  // Adjust based on cstone structure

namespace cstone {
    // Instantiate markMacsGpu<double, unsigned int> (from previous issue)
    template void markMacsGpu<double, unsigned int>(
        const unsigned int*, const int*, const int*, 
        const util::array<double, 4>*, const Box<double>&, 
        const unsigned int*, int, bool, unsigned char*
    );

    // Instantiate findHalosGpu<unsigned int, double> (new issue)
    template void findHalosGpu<unsigned int, double>(
        const unsigned int*,     // keys
        const int*,              // nodeCounts
        const int*,              // leafCounts
        const util::array<double, 3>*,  // centers
        const util::array<double, 3>*,  // sizes
        const unsigned int*,     // haloKeys
        const util::array<double, 3>*,  // haloCenters
        const util::array<double, 3>*,  // haloSizes
        const Box<double>&,      // box
        int,                     // numHalos
        int,                     // numLeaves
        unsigned char*           // halos
    );
    // Add other missing variants if needed
}
#endif