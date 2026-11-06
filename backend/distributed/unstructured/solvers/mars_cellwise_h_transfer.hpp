#pragma once

#ifdef __CUDACC__
#define BP3_H_HD __host__ __device__
#else
#define BP3_H_HD
#endif

namespace mars::cellwise::h_coarse {

struct Grid {
    int nx, ny, nz;

    BP3_H_HD long long elements() const { return (long long)nx * ny * nz; }
    BP3_H_HD Grid coarse() const { return {nx / 2, ny / 2, nz / 2}; }

    BP3_H_HD long long coarse_element(long long fine) const
    {
        const long long x = fine / ((long long)ny * nz), y = fine / nz % ny, z = fine % nz;
        return ((x / 2) * (ny / 2) + y / 2) * (nz / 2) + z / 2;
    }

    BP3_H_HD int child(long long fine) const
    {
        return int((fine / ((long long)ny * nz) % 2) * 4 + (fine / nz % ny % 2) * 2 + fine % nz % 2);
    }

    BP3_H_HD long long fine_element(long long coarse, int child) const
    {
        const Grid c = this->coarse();
        const long long x = coarse / ((long long)c.ny * c.nz), y = coarse / c.nz % c.ny, z = coarse % c.nz;
        return ((2 * x + (child >> 2)) * ny + 2 * y + ((child >> 1) & 1)) * nz + 2 * z + (child & 1);
    }
};

// Q1 coarse basis at a fine corner. Midpoints are shared by two child cells.
BP3_H_HD inline double weight(int child, int fine_node, int coarse_node)
{
    double w = 1.0;
    for (int bit = 0; bit < 3; ++bit) {
        const double x = 0.5 * (((child >> bit) & 1) + ((fine_node >> bit) & 1));
        w *= (coarse_node >> bit) & 1 ? x : 1.0 - x;
    }
    return w;
}

BP3_H_HD inline double prolong_value(const double* coarse, Grid fine, long long copy)
{
    const long long element = copy / 8, parent = fine.coarse_element(element);
    const int child = fine.child(element), node = int(copy % 8);
    double sum = 0.0;
    for (int j = 0; j < 8; ++j) {
        const double w = weight(child, node, j);
        if (w) sum += w * coarse[parent * 8 + j];
    }
    return sum;
}

BP3_H_HD inline double restrict_value(const double* raw, Grid fine, long long copy)
{
    const long long parent = copy / 8;
    const int node = int(copy % 8);
    double sum = 0.0;
    for (int child = 0; child < 8; ++child) {
        const long long element = fine.fine_element(parent, child);
        for (int i = 0; i < 8; ++i) {
            const double w = weight(child, i, node);
            if (w) sum += w * raw[element * 8 + i];
        }
    }
    return sum;
}

} // namespace mars::cellwise::h_coarse

#undef BP3_H_HD
