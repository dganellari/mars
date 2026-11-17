#pragma once
#include "examples/distributed/unstructured/mars_nonconformal_scalar_runtime.hpp"
// Procedural input only: independent node numbering and different triangulations
// on x=0. Production extraction consumes element-face side-set entries.
static Mesh mesh(int n, double rightK, int streamwise = 2, int rightResolution = 0)
{
    Mesh m;
    for (int side = 0; side < 2; ++side)
    {
        int ny = side ? (rightResolution > 0 ? rightResolution : n + 1) : n;
        int nz = ny, nx = streamwise, base = int(m.x.size());
        auto node = [=](int i, int j, int l) { return base + (i * (ny + 1) + j) * (nz + 1) + l; };
        for (int i = 0; i <= nx; ++i)
            for (int j = 0; j <= ny; ++j)
                for (int l = 0; l <= nz; ++l)
                {
                    m.x.push_back({double(side - 1) + double(i) / nx, double(j) / ny, double(l) / nz});
                    m.side.push_back(side);
                }
        for (int i = 0; i < nx; ++i)
            for (int j = 0; j < ny; ++j)
                for (int l = 0; l < nz; ++l)
                {
                    std::array<int, 3> perm{0, 1, 2};
                    do
                    {
                        int e = int(m.elements.size());
                        m.elements.push_back(e);
                        m.owners.push_back(e % rankCount);
                        m.k.push_back(side ? rightK : 1.);
                        int p[]{i, j, l};
                        m.cells.push_back(node(p[0], p[1], p[2]));
                        for (int d : perm)
                        {
                            ++p[d];
                            m.cells.push_back(node(p[0], p[1], p[2]));
                        }
                        for (int opp = 0; opp < 4; ++opp)
                        {
                            bool seam = true;
                            for (int a = 0; a < 4; ++a)
                                if (a != opp) seam = seam && m.x[m.cells[4 * e + a]].x == 0;
                            if (seam) m.tags.push_back({e, opp, side});
                        }
                    } while (std::next_permutation(perm.begin(), perm.end()));
                }
    }
    m.dofs.resize(m.x.size());
    std::iota(m.dofs.begin(), m.dofs.end(), 0);
    return m;
}
