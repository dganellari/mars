#pragma once
#include "backend/distributed/unstructured/fem/mars_nonconformal_scalar.hpp"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <limits>
#if defined(__CUDACC__) || defined(__HIPCC__)
#include "backend/distributed/unstructured/fem/mars_nonconformal_scalar_device.hpp"
#include <thrust/copy.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#endif

using namespace mars::fem;
using Face        = InterfaceTetFace<double>;
using Point       = InterfacePoint<double>;
using Stencil     = InterfaceStencil<double>;
inline int checks = 0;
inline void require(bool ok, const char* name)
{
    ++checks;
    if (!ok)
    {
        std::fprintf(stderr, "FAIL: %s\n", name);
        std::exit(1);
    }
}
inline void near(double a, double b, const char* name) { require(std::abs(a - b) < 2e-12 * (1 + std::abs(b)), name); }

inline std::vector<Face> fixture()
{
    std::vector<Face> f{{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, -1}}, {0, 1, 2, 3}, 10, 0, 3, 2.}};
    Vec3<double> x[]   = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {.5, 0, 0}, {.5, .5, 0}, {0, .5, 0}, {0, 0, 1}};
    int triangles[][3] = {{0, 3, 5}, {3, 1, 4}, {5, 4, 2}, {3, 4, 5}};
    for (int i = 0; i < 4; ++i)
    {
        Face b{};
        b.id          = 20 + i;
        b.side        = 1;
        b.opposite    = 3;
        b.diffusivity = 2.;
        for (int j = 0; j < 3; ++j)
        {
            b.x[j]   = x[triangles[i][j]];
            b.dof[j] = 4 + triangles[i][j];
        }
        b.x[3]   = x[6];
        b.dof[3] = 10;
        f.push_back(b);
    }
    return f;
}
// Two triangle-centroid rules for each median-dual face patch: six IPs per face.
inline std::vector<Point> points(const std::vector<Face>& faces)
{
    std::vector<Point> p;
    for (int f = 0; f < int(faces.size()); ++f)
    {
        double area = f ? .125 : .5;
        for (int row = 0; row < 3; ++row)
            for (int edge = 1; edge <= 2; ++edge)
            {
                Point q{};
                q.face                       = f;
                q.rowCorner                  = row;
                q.area                       = area / 6;
                q.bary[row]                  = 11. / 18;
                q.bary[(row + edge) % 3]     = 5. / 18;
                q.bary[(row + 3 - edge) % 3] = 1. / 9;
                p.push_back(q);
            }
    }
    return p;
}
inline double affine(Vec3<double> x) { return 1 + 2 * x.x - 3 * x.y + 4 * x.z; }
inline double applyStencil(const Stencil& s, const std::vector<Face>& f)
{
    double values[11]{};
    for (auto& e : f)
        for (int j = 0; j < 4; ++j)
            values[e.dof[j]] = affine(e.x[j]);
    double y = 0;
    for (int j = 0; j < 8; ++j)
        y += s.value[j] * values[s.column[j]];
    return y;
}
inline void deviceCheck(const std::vector<Point>& p,
                        const std::vector<Face>& faces,
                        const std::vector<int>& offsets,
                        const std::vector<int>& candidates,
                        InterfacePairingOptions<double> opts)
{
#if defined(__CUDACC__) || defined(__HIPCC__)
    std::vector<Stencil> expected(p.size());
    buildNonconformalScalarStencilsHost(int(p.size()), p.data(), faces.data(), int(faces.size()), offsets.data(),
                                        candidates.data(), int(candidates.size()), opts, expected.data());
    thrust::device_vector<Point> dp(p.begin(), p.end());
    thrust::device_vector<Face> df(faces.begin(), faces.end());
    thrust::device_vector<int> offsets_d(offsets.begin(), offsets.end()),
        candidates_d(candidates.begin(), candidates.end());
    thrust::device_vector<Stencil> output(p.size());
    cudaStream_t stream;
    require(cudaStreamCreate(&stream) == cudaSuccess, "create GPU stream");
    require(buildNonconformalScalarStencils(int(p.size()), thrust::raw_pointer_cast(dp.data()),
                                            thrust::raw_pointer_cast(df.data()), int(faces.size()),
                                            thrust::raw_pointer_cast(offsets_d.data()),
                                            thrust::raw_pointer_cast(candidates_d.data()), int(candidates.size()), opts,
                                            thrust::raw_pointer_cast(output.data()), stream) == cudaSuccess,
            "launch device pairing/stencils");
    require(cudaStreamSynchronize(stream) == cudaSuccess, "device completion");
    thrust::host_vector<Stencil> result = output;
    require(cudaStreamDestroy(stream) == cudaSuccess, "destroy stream");
    for (size_t i = 0; i < p.size(); ++i)
    {
        require(result[i].status == expected[i].status, "device status");
        if (expected[i].status != InterfaceStatus::valid) continue;
        require(result[i].row == expected[i].row && result[i].donor == expected[i].donor, "device pairing");
        for (int j = 0; j < 8; ++j)
        {
            require(result[i].column[j] == expected[i].column[j], "device column");
            near(result[i].value[j], expected[i].value[j], "device coefficient");
        }
    }
#else
    (void)p;
    (void)faces;
    (void)offsets;
    (void)candidates;
    (void)opts;
#endif
}
