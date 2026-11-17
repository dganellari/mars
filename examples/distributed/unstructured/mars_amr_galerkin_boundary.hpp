#pragma once
#include <cmath>
#include "backend/distributed/unstructured/amr/mars_amr_lepp_boundary.hpp"

#if defined(__CUDACC__) || defined(__HIPCC__)
#define MARS_PDE_HD __host__ __device__
#else
#define MARS_PDE_HD
#endif

namespace mars::amr_pde
{
template<class KeyType>
struct BoundaryFace
{
    amr::LeppFaceKey<KeyType> key;
    int labelled;
};

// Count owned-cell faces at SFC owners. A Dirichlet surface must cover every exterior face and no interior face.
template<class Execution, class KeyType, class Owner>
uint64_t check_boundary(const KeyType* keys,
                        std::array<const KeyType*, 4> conn,
                        size_t cells,
                        const typename Execution::template Vector<amr::LeppTaggedFace>& labels,
                        Owner owner,
                        MPI_Comm comm,
                        bool staged)
{
    using Packet   = BoundaryFace<KeyType>;
    const size_t n = 4 * cells + labels.size();
    typename Execution::template Vector<Packet> packets(n);
    typename Execution::template Vector<int> destinations(n);
    auto p    = Execution::data(packets);
    auto d    = Execution::data(destinations);
    auto tags = Execution::data(labels);
    Execution::for_each(n,
                        [=] MARS_PDE_HD(size_t i)
                        {
                            size_t e = i < 4 * cells ? i / 4 : tags[i - 4 * cells].element;
                            int f    = i < 4 * cells ? int(i % 4) : tags[i - 4 * cells].face;
                            KeyType k[3];
                            for (int j = 0, v = 0; v < 4; ++v)
                                if (v != f) k[j++] = keys[conn[v][e]];
                            p[i] = {amr::LeppFaceKey<KeyType>(k[0], k[1], k[2]), i >= 4 * cells};
                            d[i] = owner(p[i].key.keys[0]);
                        });
    amr::AmrSparsePlan<Execution, Packet> plan(destinations, comm, staged);
    auto received   = plan.forward(packets);
    auto r          = Execution::data(received);
    const size_t nr = received.size();
    Execution::sort(r, nr, [] MARS_PDE_HD(const Packet& a, const Packet& b) { return a.key < b.key; });
    typename Execution::template Vector<uint64_t> invalid(nr), boundary(nr);
    auto bad = Execution::data(invalid);
    auto out = Execution::data(boundary);
    Execution::for_each(nr,
                        [=] MARS_PDE_HD(size_t i)
                        {
                            bad[i] = out[i] = 0;
                            if (i && r[i - 1].key == r[i].key) return;
                            size_t cells_on_face = 0, labels_on_face = 0;
                            for (size_t j = i; j < nr && r[j].key == r[i].key; ++j)
                                if (r[j].labelled)
                                    ++labels_on_face;
                                else
                                    ++cells_on_face;
                            out[i] = cells_on_face == 1;
                            bad[i] = !((cells_on_face == 1 && labels_on_face > 0) ||
                                       (cells_on_face == 2 && labels_on_face == 0));
                        });
    amr::require_collectively(Execution::sum(bad, nr) == 0,
                              "Dirichlet side sets must cover every exterior face and no interior face", comm);
    uint64_t local = Execution::sum(out, nr), global = 0;
    MPI_Allreduce(&local, &global, 1, MPI_UINT64_T, MPI_SUM, comm);
    amr::require_collectively(global > 0, "Dirichlet surface is empty", comm);
    return global;
}

template<class KeyType>
struct BoundaryArea
{
    const amr::LeppTaggedFace* labels;
    std::array<const KeyType*, 4> conn;
    const double *x, *y, *z;
    int tag;
    MARS_PDE_HD double operator()(size_t i) const
    {
        if (labels[i].tag != tag) return 0;
        KeyType v[3];
        for (int c = 0, j = 0; c < 4; ++c)
            if (c != labels[i].face) v[j++] = conn[c][labels[i].element];
        double a[3] = {x[v[1]] - x[v[0]], y[v[1]] - y[v[0]], z[v[1]] - z[v[0]]};
        double b[3] = {x[v[2]] - x[v[0]], y[v[2]] - y[v[0]], z[v[2]] - z[v[0]]};
        double c[3] = {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
        return .5 * sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
    }
};

// A four-point tensor Gauss rule after the Duffy map integrates squared quadratic error exactly.
template<class Problem>
MARS_PDE_HD double tet_error(const double (&x)[4][3], const double (&u)[4], const Problem& problem)
{
    constexpr double q[4] = {.0694318442029737124, .330009478207571868, .669990521792428132, .930568155797026288};
    constexpr double w[4] = {.173927422568726929, .326072577431273071, .326072577431273071, .173927422568726929};
    double a[3], b[3], c[3];
    for (int d = 0; d < 3; ++d)
    {
        a[d] = x[1][d] - x[0][d];
        b[d] = x[2][d] - x[0][d];
        c[d] = x[3][d] - x[0][d];
    }
    const double det =
        a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
    double error = 0;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            for (int k = 0; k < 4; ++k)
            {
                double phi[4]   = {0, q[i], (1 - q[i]) * q[j], (1 - q[i]) * (1 - q[j]) * q[k]};
                phi[0]          = 1 - phi[1] - phi[2] - phi[3];
                double point[3] = {}, value = 0;
                for (int v = 0; v < 4; ++v)
                {
                    value += phi[v] * u[v];
                    for (int d = 0; d < 3; ++d)
                        point[d] += phi[v] * x[v][d];
                }
                const double delta = value - problem.u(point[0], point[1], point[2]);
                error += delta * delta * w[i] * w[j] * w[k] * (1 - q[i]) * (1 - q[i]) * (1 - q[j]) * det;
            }
    return error;
}
} // namespace mars::amr_pde
#undef MARS_PDE_HD
