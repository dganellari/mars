#pragma once
#include "fixture.hpp"
#include "backend/distributed/unstructured/fem/mars_nonconformal_overlap.hpp"

inline std::vector<Stencil> conservativeFixture(const std::vector<Face>& faces,
                                                const std::vector<InterfaceFacePair>& pairs,
                                                double tolerance,
                                                bool device = false)
{
    std::vector<Stencil> output;
#if defined(__CUDACC__) || defined(__HIPCC__)
    if (device)
    {
        thrust::device_vector<Face> df(faces);
        thrust::device_vector<InterfaceFacePair> dp(pairs);
        thrust::device_vector<Stencil> ds(pairs.size() * nonconformal_pair_capacity);
        thrust::device_vector<int> dc(pairs.size());
        thrust::device_vector<InterfaceStatus> dst(pairs.size());
        cudaStream_t stream;
        require(cudaStreamCreate(&stream) == cudaSuccess, "overlap stream");
        require(buildNonconformalConservativePairs(
                    int(pairs.size()), thrust::raw_pointer_cast(dp.data()), thrust::raw_pointer_cast(df.data()),
                    int(faces.size()), InterfacePairingOptions<double>{tolerance, 1e-12},
                    thrust::raw_pointer_cast(ds.data()), thrust::raw_pointer_cast(dc.data()),
                    thrust::raw_pointer_cast(dst.data()), stream) == cudaSuccess,
                "overlap device launch");
        require(cudaStreamSynchronize(stream) == cudaSuccess, "overlap device completion");
        thrust::host_vector<Stencil> hs          = ds;
        thrust::host_vector<int> hc              = dc;
        thrust::host_vector<InterfaceStatus> hst = dst;
        require(cudaStreamDestroy(stream) == cudaSuccess, "overlap stream destroy");
        for (size_t i = 0; i < pairs.size(); ++i)
        {
            require(hst[i] == InterfaceStatus::valid, "device overlap geometry");
            require(hc[i] >= 0 && hc[i] <= nonconformal_pair_capacity && hc[i] % 2 == 0, "device overlap count");
            output.insert(output.end(), hs.begin() + i * nonconformal_pair_capacity,
                          hs.begin() + i * nonconformal_pair_capacity + hc[i]);
        }
        return output;
    }
#else
    (void)device;
#endif
    for (auto p : pairs)
    {
        Stencil out[nonconformal_pair_capacity];
        int count = 0;
        require(nonconformalConservativePair(faces[p.first], faces[p.second],
                                             InterfacePairingOptions<double>{tolerance, 1e-12}, out,
                                             nonconformal_pair_capacity, count) == InterfaceStatus::valid,
                "host overlap geometry");
        output.insert(output.end(), out, out + count);
    }
    return output;
}
inline void overlapChecks()
{
    using namespace nonconformal_overlap_detail;
    Point2<double> a[]{{0, 0}, {1, 0}, {0, 1}}, b[]{{0, 0}, {1, 0}, {0, 1}}, c;
    double area;
    require(overlap(a, b, area, c), "triangle intersection");
    near(area, .5, "triangle area");
    near(c.x, 1. / 3, "triangle centroid x");
    near(c.y, 1. / 3, "triangle centroid y");
    b[0] = {.5, 0};
    b[1] = {1.5, 0};
    b[2] = {.5, 1};
    require(overlap(a, b, area, c), "partial triangle intersection");
    near(area, .125, "partial overlap area");
    near(c.x, 2. / 3, "partial centroid x");
    near(c.y, 1. / 6, "partial centroid y");
    b[0] = {1, 0};
    b[1] = {2, 0};
    b[2] = {1, 1};
    require(overlap(a, b, area, c), "point contact");
    near(area, 0, "zero contact area");
    std::vector<InterfaceFacePair> pairs{{0, 1}, {0, 2}, {0, 3}, {0, 4}};
    for (double scale : {0.001, 1., 1000.})
        for (bool rotate : {false, true})
            for (bool varying : {false, true})
                for (bool uneven : {false, true})
                {
                    auto f = fixture();
                    for (auto& face : f)
                    {
                        if (varying) face.diffusivity = 1. + .1 * (face.id - 10);
                        for (int corner = 0; corner < 4; ++corner)
                        {
                            auto& x = face.x[corner];
                            if (uneven)
                            {
                                if (face.dof[corner] == 7) x = {.23, 0, 0};
                                if (face.dof[corner] == 8) x = {.63, .37, 0};
                                if (face.dof[corner] == 9) x = {0, .41, 0};
                            }
                            x = {scale * x.x, scale * x.y, scale * x.z};
                            if (rotate)
                                x = {.6 * x.x - .8 * x.z + 2 * scale, x.y - 3 * scale, .8 * x.x + .6 * x.z + 4 * scale};
                        }
                    }
                    auto h = conservativeFixture(f, pairs, 1e-11 * scale);
                    require(!h.empty(), "overlap stencils present");
                    auto verify = [&](const std::vector<Stencil>& s)
                    {
                        double column[11]{}, row[11]{}, affineFlux[2]{}, jumpFlux[2]{};
                        for (size_t i = 0; i < s.size(); ++i)
                        {
                            double sum = 0;
                            for (int j = 0; j < 8; ++j)
                            {
                                sum += s[i].value[j];
                                column[s[i].column[j]] += s[i].value[j];
                                if (s[i].column[j] < 4) jumpFlux[i % 2] += s[i].value[j];
                            }
                            row[s[i].row] += sum;
                            affineFlux[i % 2] += applyStencil(s[i], f);
                            if (i % 2)
                                for (int j = 0; j < 8; ++j)
                                    require(s[i].column[j] == s[i - 1].column[j] && s[i].value[j] == -s[i - 1].value[j],
                                            "equal opposite interface coefficient");
                        }
                        for (int i = 0; i < 11; ++i)
                        {
                            near(column[i] / scale, 0, "all basis fields conserve interface flux");
                            near(row[i] / scale, 0, "constant nullspace on every row");
                        }
                        if (!varying)
                        {
                            double gradientNormal = rotate ? -.8 * 2. + .6 * 4. : 4.;
                            near(affineFlux[0] / (scale * scale), -gradientNormal, "analytic affine flux");
                            near(affineFlux[1] / (scale * scale), gradientNormal, "opposite affine flux");
                            near(jumpFlux[0] / scale, 1., "constant jump integrated area");
                            near(jumpFlux[1] / scale, -1., "opposite jump integrated area");
                        }
                    };
                    verify(h);
                    auto reversed = pairs;
                    for (auto& pair : reversed)
                        std::swap(pair.first, pair.second);
                    auto reverse = conservativeFixture(f, reversed, 1e-11 * scale);
                    double forwardMatrix[11][11]{}, reverseMatrix[11][11]{};
                    for (const auto& s : h)
                        for (int j = 0; j < 8; ++j)
                            forwardMatrix[s.row][s.column[j]] += s.value[j];
                    for (const auto& s : reverse)
                        for (int j = 0; j < 8; ++j)
                            reverseMatrix[s.row][s.column[j]] += s.value[j];
                    for (int i = 0; i < 11; ++i)
                        for (int j = 0; j < 11; ++j)
                            near(forwardMatrix[i][j] / scale, reverseMatrix[i][j] / scale,
                                 "pair reversal matrix invariance");
#if defined(__CUDACC__) || defined(__HIPCC__)
                    auto d = conservativeFixture(f, pairs, 1e-11 * scale, true);
                    verify(d);
                    // Compare assembled coefficients: clipping at a touching edge may
                    // split a zero-area patch differently under device contraction.
                    double hm[11][11]{}, dm[11][11]{};
                    for (const auto& s : h)
                        for (int j = 0; j < 8; ++j)
                            hm[s.row][s.column[j]] += s.value[j];
                    for (const auto& s : d)
                        for (int j = 0; j < 8; ++j)
                            dm[s.row][s.column[j]] += s.value[j];
                    for (int i = 0; i < 11; ++i)
                        for (int j = 0; j < 11; ++j)
                            near(dm[i][j] / scale, hm[i][j] / scale, "device overlap matrix");
#endif
                }
    auto f = fixture();
    Stencil out[72];
    int count    = 0;
    auto options = InterfacePairingOptions<double>{1e-12};
    require(nonconformalConservativePair(f[0], f[1], options, out, 71, count) == InterfaceStatus::invalid_input &&
                count == 0,
            "overlap capacity refusal");
    auto bad   = f[1];
    bad.dof[0] = f[0].dof[0];
    require(nonconformalConservativePair(f[0], bad, options, out, 72, count) == InterfaceStatus::invalid_input,
            "shared dof refusal");
    bad = f[1];
    for (auto& x : bad.x)
        x.z += .01;
    require(nonconformalConservativePair(f[0], bad, options, out, 72, count) == InterfaceStatus::invalid_geometry,
            "noncoplanar refusal");
    bad = f[1];
    for (auto& x : bad.x)
        x.x += 2;
    require(nonconformalConservativePair(f[0], bad, options, out, 72, count) == InterfaceStatus::valid && count == 0,
            "disjoint faces");
}
