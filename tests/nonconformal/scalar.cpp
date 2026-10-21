#include "fixture.hpp"
#include "overlap.hpp"
int main()
{
    overlapChecks();
    auto faces = fixture();
    auto p     = points(faces);
    std::vector<int> offsets{0}, candidates;
    for (auto q : p)
    {
        if (q.face == 0)
            candidates.insert(candidates.end(), {1, 2, 3, 4});
        else
            candidates.push_back(0);
        offsets.push_back(candidates.size());
    }
    std::vector<Stencil> out(p.size());
    InterfacePairingOptions<double> opts{1e-12, 1e-12};
    require(buildNonconformalScalarStencilsHost(int(p.size()), p.data(), faces.data(), int(faces.size()),
                                                offsets.data(), candidates.data(), int(candidates.size()), opts,
                                                out.data()),
            "batch valid");
    deviceCheck(p, faces, offsets, candidates, opts);
    double side[2]{};
    for (size_t i = 0; i < p.size(); ++i)
    {
        auto& s = out[i];
        require(s.status == InterfaceStatus::valid, "nonmatching face pairing");
        double sum = 0;
        for (double v : s.value)
            sum += v;
        near(sum, 0, "constant field null flux");
        double expected = p[i].area * (p[i].face == 0 ? -8. : 8.);
        near(applyStencil(s, faces), expected, "affine normal flux");
        side[faces[p[i].face].side] += applyStencil(s, faces);
        require(s.row == faces[p[i].face].dof[p[i].rowCorner], "current CVFEM row");
        for (int j = 0; j < 4; ++j)
            require(s.column[j] != s.column[4 + j], "distinct-side DOFs");
    }
    near(side[0], -4, "integrated current flux");
    near(side[1], 4, "integrated donor flux");
    near(side[0] + side[1], 0, "affine dual-pass balance");
    // Unequal diffusivities: continuous trace and continuous physical diffusive flux.
    auto heterogeneous = faces;
    for (auto& f : heterogeneous)
        if (f.side == 1) f.diffusivity = 4.;
    std::vector<Stencil> hs(p.size());
    buildNonconformalScalarStencilsHost(int(p.size()), p.data(), heterogeneous.data(), 5, offsets.data(),
                                        candidates.data(), int(candidates.size()), opts, hs.data());
    double hv[11]{};
    for (auto& f : heterogeneous)
        for (int j = 0; j < 4; ++j)
            hv[f.dof[j]] = 1 + 2 * f.x[j].x - 3 * f.x[j].y + (f.side == 0 ? 4. : 2.) * f.x[j].z;
    for (size_t i = 0; i < p.size(); ++i)
    {
        require(hs[i].status == InterfaceStatus::valid, "heterogeneous pairing");
        double value = 0;
        for (int j = 0; j < 8; ++j)
            value += hs[i].value[j] * hv[hs[i].column[j]];
        near(value, p[i].area * (p[i].face == 0 ? -8. : 8.), "heterogeneous continuous flux");
        double jumpFlux = 0;
        for (int j = 0; j < 4; ++j)
            jumpFlux += hs[i].value[j];
        near(jumpFlux, 3 * p[i].area, "positive penalty for current-side jump");
    }
    deviceCheck(p, heterogeneous, offsets, candidates, opts);
    int c[]        = {1, 2, 3, 4};
    auto q         = p[0];
    auto reference = nonconformalScalarStencil(q, faces.data(), 5, c, 4, opts);
    std::reverse(c, c + 4);
    auto reordered = nonconformalScalarStencil(q, faces.data(), 5, c, 4, opts);
    require(reference.donor == reordered.donor, "candidate order independence");
    // Independent analytic derivative of the interface flux, including opposing-element columns.
    const Face* b = nullptr;
    for (auto& e : faces)
        if (e.id == reference.donor) b = &e;
    require(b != nullptr, "donor lookup");
    double ua[4] = {.7, -.2, .6, 1.1}, ub[4] = {-.3, .8, 1.2, .4};
    // Both fixture heights are one. Tangential derivatives do not enter the normal flux.
    auto flux = [&]()
    {
        double x = q.bary[1], y = q.bary[2];
        double x0 = b->x[0].x, y0 = b->x[0].y;
        double dx1 = b->x[1].x - x0, dy1 = b->x[1].y - y0, dx2 = b->x[2].x - x0, dy2 = b->x[2].y - y0;
        double d  = dx1 * dy2 - dx2 * dy1;
        double w1 = ((x - x0) * dy2 - (y - y0) * dx2) / d, w2 = (dx1 * (y - y0) - dy1 * (x - x0)) / d;
        double w0 = 1 - w1 - w2, pa = q.bary[0] * ua[0] + q.bary[1] * ua[1] + q.bary[2] * ua[2];
        // Extrapolate each face's affine trace to the apex projection (0,0).
        double v1 = (-x0 * dy2 + y0 * dx2) / d, v2 = (-dx1 * y0 + dy1 * x0) / d;
        double atApex = (1 - v1 - v2) * ub[0] + v1 * ub[1] + v2 * ub[2];
        double qa = -2 * (ua[0] - ua[3]), qb = 2 * (ub[3] - atApex);
        return q.area * (.5 * (qa - qb) + 2 * (pa - (w0 * ub[0] + w1 * ub[1] + w2 * ub[2])));
    };
    double base = flux(), linear = 0;
    for (int j = 0; j < 4; ++j)
        linear += reference.value[j] * ua[j] + reference.value[4 + j] * ub[j];
    near(linear, base, "independent flux matches stencil");
    for (int j = 0; j < 8; ++j)
    {
        double& v = j < 4 ? ua[j] : ub[j - 4];
        v += 1;
        near(flux() - base, reference.value[j], "eight-column flux Jacobian");
        v -= 1;
    }
    // Rotate, translate and rescale both meshes; integrated diffusion scales as length.
    for (double length : {1e-4, 2., 1e4})
    {
        auto moved = faces;
        for (auto& e : moved)
            for (auto& x : e.x)
                x = {length * (x.z + 3), length * (x.x - 2), length * (x.y + 1)};
        auto qp = q;
        qp.area *= length * length;
        auto ms = nonconformalScalarStencil(qp, moved.data(), 5, c, 4,
                                            InterfacePairingOptions<double>{1e-12 * length, 1e-12});
        require(ms.status == InterfaceStatus::valid, "scaled rotated geometry");
        for (int j = 0; j < 8; ++j)
            near(ms.value[j], reference.value[j] * length, "scaled stencil");
    }
    auto permuted = faces;
    for (auto& f : permuted)
    {
        std::swap(f.x[0], f.x[1]);
        std::swap(f.dof[0], f.dof[1]);
    }
    auto permutedPoint = q;
    std::swap(permutedPoint.bary[0], permutedPoint.bary[1]);
    permutedPoint.rowCorner = 1;
    auto ps                 = nonconformalScalarStencil(permutedPoint, permuted.data(), 5, c, 4, opts);
    require(ps.status == InterfaceStatus::valid, "reversed connectivity accepted");
    near(applyStencil(ps, permuted), applyStencil(reference, faces), "connectivity orientation invariant flux");
    // Face-edge ties select the same stable face ID regardless of candidate order.
    auto edge    = q;
    edge.bary[0] = .5;
    edge.bary[1] = .5;
    edge.bary[2] = 0;
    auto e1      = nonconformalScalarStencil(edge, faces.data(), 5, c, 4, opts);
    std::reverse(c, c + 4);
    auto e2 = nonconformalScalarStencil(edge, faces.data(), 5, c, 4, opts);
    require(e1.status == InterfaceStatus::valid && e1.donor == e2.donor, "edge donor tie");
    auto bad = faces;
    for (auto& e : bad)
        if (e.side == 1)
            for (auto& x : e.x)
                x.z += .01;
    require(nonconformalScalarStencil(q, bad.data(), 5, c, 4, opts).status == InterfaceStatus::no_donor,
            "gap rejected");
    deviceCheck(p, bad, offsets, candidates, opts);
    q.bary[0] = -.1;
    require(nonconformalScalarStencil(q, faces.data(), 5, c, 4, opts).status == InterfaceStatus::invalid_input,
            "invalid point rejected");
    q           = p[0];
    bad         = faces;
    bad[1].x[3] = bad[1].x[0];
    int one     = 1;
    require(nonconformalScalarStencil(q, bad.data(), 5, &one, 1, opts).status == InterfaceStatus::invalid_geometry,
            "degenerate donor rejected");
    bad           = faces;
    bad[1].x[3].z = -1;
    require(nonconformalScalarStencil(q, bad.data(), 5, &one, 1, opts).status == InterfaceStatus::no_donor,
            "same-facing donor rejected");
    bad           = faces;
    bad[1].dof[0] = 0;
    require(nonconformalScalarStencil(q, bad.data(), 5, &one, 1, opts).status == InterfaceStatus::invalid_input,
            "merged side DOFs rejected");
    bad                = faces;
    bad[0].diffusivity = std::numeric_limits<double>::quiet_NaN();
    require(nonconformalScalarStencil(q, bad.data(), 5, c, 4, opts).status == InterfaceStatus::invalid_input,
            "nonfinite coefficient rejected");
    one = 99;
    require(nonconformalScalarStencil(q, faces.data(), 5, &one, 1, opts).status == InterfaceStatus::invalid_input,
            "invalid candidate index rejected");
    offsets[1] = int(candidates.size()) + 1;
    buildNonconformalScalarStencilsHost(int(p.size()), p.data(), faces.data(), 5, offsets.data(), candidates.data(),
                                        int(candidates.size()), opts, out.data());
    require(out[0].status == InterfaceStatus::invalid_input, "invalid candidate span rejected");
    deviceCheck(p, faces, offsets, candidates, opts);
    require(buildNonconformalScalarStencilsHost<double>(0, nullptr, nullptr, 0, nullptr, nullptr, 0, opts, nullptr),
            "empty batch");
    require(!buildNonconformalScalarStencilsHost<double>(1, nullptr, nullptr, 0, nullptr, nullptr, 0, opts, nullptr),
            "null batch rejected");
    require(nonconformalScalarStencil(q, faces.data(), 5, nullptr, 0, opts).status == InterfaceStatus::no_donor,
            "empty candidates rejected");
    q.rowCorner = 3;
    require(nonconformalScalarStencil(q, faces.data(), 5, c, 4, opts).status == InterfaceStatus::invalid_input,
            "off-face row rejected");
    std::printf("PASS: %d stationary non-conformal scalar checks\n", checks);
}
