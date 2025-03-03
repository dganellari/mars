/* #include "cstone/domain/domain.hpp"

using Real = double;
using KeyType = unsigned;

int main() {
    int rank = 0, numRanks = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks); */

    // fill x,y,z,h with some initial values on each rank
    //std::vector<Real> x{/*...*/}, y{/*...*/}, z{/*...*/}, h{/*...*/};

    /* int bucketSize = 10;
    cstone::Domain<KeyType, Real> domain(rank, numRanks, bucketSize);

    int nIterations = 10;
    for (int iter = 0; iter < nIterations; ++iter) {
        domain.sync(x, y, z, h);

        // x,y,z,h now contain all particles of a part of the global octree,
        // including their halos.

        std::vector<Real> density(x.size());

        // compute physical quantities, e.g. densities for particles in the assigned ranges:
        // computeDensity(density,x,y,z,h,domain.startIndex(),domain.endIndex());

        // repeat the halo exchange for densities
        domain.exchangeHalos(density);

        // compute more quantities and finally update particle positions in x,y,z and h,
        // then start over
    }

    return;
} */
