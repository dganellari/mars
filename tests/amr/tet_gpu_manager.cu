#include "tet_gpu_checks.hpp"
#include "public_tet_mesh.hpp"
#include "backend/distributed/unstructured/amr/mars_amr_manager.hpp"

namespace amr_test
{
using Manager = mars::amr::AmrManager<mars::TetTag, KeyType, RealType>;

void label_boundary(Domain& domain, int rank)
{
    auto cells = snapshot(domain, rank);
    std::array<std::vector<KeyType>, 7> keys;
    for (const auto& cell : cells)
        for (int f = 0; f < 4; ++f)
        {
            int tag = boundary_plane(tet_face(cell, f));
            if (tag < 0) continue;
            for (int c = 0; c < 4; ++c)
                if (c != f)
                {
                    keys[tag].push_back(cell.keys[c]);
                    if (tag == 0) keys[6].push_back(cell.keys[c]);
                }
        }
    for (int tag = 0; tag < 7; ++tag)
        domain.set_side_set_from_face_keys("boundary" + std::to_string(tag), upload(keys[tag]));
}

void check_boundary(Domain& domain, const std::vector<Tet>& global, bool imported)
{
    auto keys = read_device(domain.getLocalToGlobalSfcMap().data(), domain.getNodeCount());
    auto x    = read_device(domain.getNodeX().data(), keys.size());
    auto y    = read_device(domain.getNodeY().data(), keys.size());
    auto z    = read_device(domain.getNodeZ().data(), keys.size());
    std::map<KeyType, Point> coordinates;
    for (size_t i = 0; i < keys.size(); ++i)
        coordinates[keys[i]] = {x[i], y[i], z[i]};
    require(domain.sideSetNames().size() == 7, "Boundary set names changed");
    for (int tag = 0; tag < 7; ++tag)
    {
        auto name = "boundary" + std::to_string(tag);
        require(domain.has_side_set_face_keys(name), "Boundary face ancestry was lost");
        const auto& device_faces = domain.side_set_face_keys(name);
        auto face_keys           = read_device(device_faces.data(), device_faces.size());
        std::vector<Face> local;
        std::set<KeyType> expected_nodes;
        for (size_t i = 0; i < face_keys.size(); i += 3)
        {
            Face face;
            for (int c = 0; c < 3; ++c)
            {
                face[c] = coordinates.at(face_keys[i + c]);
                expected_nodes.insert(face_keys[i + c]);
            }
            std::sort(face.begin(), face.end());
            local.push_back(face);
        }
        auto node_ids = read_device(domain.sideSetNodes(name).data(), domain.sideSetNodes(name).size());
        if (imported)
        {
            // A file-loaded rank can hold a boundary node without holding a complete incident boundary face.
            expected_nodes.clear();
            int plane = tag == 6 ? 0 : tag;
            for (const auto& [key, point] : coordinates)
                if (std::abs(point[plane / 2] - plane % 2) < 2.e-6) expected_nodes.insert(key);
        }
        std::set<KeyType> actual_nodes;
        for (int i : node_ids)
            actual_nodes.insert(keys.at(i));
        require(actual_nodes == expected_nodes, "Boundary node list disagrees with its tagged faces");
        auto faces = gather(local);
        std::set<Face> actual(faces.begin(), faces.end()), expected;
        for (const auto& tet : global)
            for (int f = 0; f < 4; ++f)
            {
                auto face = tet_face(tet, f);
                if (boundary_plane(face) == (tag == 6 ? 0 : tag)) expected.insert(face);
            }
        require(!expected.empty() && actual == expected, "Boundary triangles lost or relabeled during AMR");
    }
}

Manager::Config manager_config(bool staged)
{
    Manager::Config config;
    config.tet_refinement            = mars::amr::TetRefinement::ConformingBisection;
    config.stage_refinement_messages = staged;
    config.strategy                  = mars::amr::MarkingStrategy::Doerfler;
    // Allow the additional bisections required by conformity closure.
    config.maxLevels       = 32;
    config.bucketSize      = 8;
    config.refineFraction  = 0.5;
    config.coarsenFraction = 0;
    return config;
}

void check_manager_rounds(Manager& manager, int rank)
{
    auto fields = initial_fields(manager.domain());
    (void)check_identity(manager.domain(), rank);

    RealType affine_tolerance = 1.e-12;
    for (int round = 0; round < 3; ++round)
    {
        auto& current = manager.domain();
        auto local    = snapshot(current, rank);
        auto before   = gather(local);
        check_mesh(before);
        // Needed ghosts are refreshed by the manager; deeper halo copies must never become field sources.
        poison_ghost_fields(current, fields);
        std::vector<RealType> errors(current.getElementCount(), 0);
        for (size_t e = 0; e < local.size(); ++e)
        {
            RealType center = 0;
            for (const auto& point : local[e].points)
                center += point[round % 3] / 4;
            errors[current.startIndex() + e] = center < 0.4 ? 1 : 0;
        }
        auto device_errors = upload(errors);
        std::array<Values, 2> next_fields;
        auto stats    = manager.adaptMeshMultiField(device_errors.data(), {fields[0].data(), fields[1].data()},
                                                    {&next_fields[0], &next_fields[1]});
        auto& updated = manager.domain();
        auto after    = gather(snapshot(updated, rank));
        require(after.size() > before.size() && stats.elementsRefined > 0, "AMR did not refine the marked region");
        check_mesh(after);
        check_boundary(updated, after);
        // Midpoint averaging cannot increase an existing field error. Re-encoding adds at most
        // one coordinate grid interval per axis, weighted by this affine field's gradient.
        affine_tolerance += affine_quantization_bound(updated);
        auto valid_nodes = current_nodes(updated);
        check_fields(updated.getNodeX().data(), updated.getNodeY().data(), updated.getNodeZ().data(),
                     updated.getNodeCount(), next_fields, affine_tolerance, &valid_nodes);
        (void)check_identity(updated, rank);
        fields = std::move(next_fields);
        if (rank == 0)
            std::cout << "Manager round " << round + 1 << ": " << before.size() << " -> " << after.size()
                      << " tets, affine bound=" << affine_tolerance << ", conformity/volume/identity/fields PASS\n";
    }
}

void check_manager(std::unique_ptr<Domain> domain, int rank, int ranks, bool staged)
{
    Manager manager(manager_config(staged));
    manager.initialize(std::move(domain), rank, ranks);
    label_boundary(manager.domain(), rank);
    check_boundary(manager.domain(), gather(snapshot(manager.domain(), rank)));
    check_manager_rounds(manager, rank);
}

void check_exodus_manager(const std::string& path, int rank, int ranks, bool staged)
{
#ifndef MARS_HAVE_NETCDF
    throw std::runtime_error("The Exodus AMR gate requires netCDF in the installed MARS package");
#else
    auto config          = manager_config(staged);
    config.storeSideSets = true;
    Manager manager(config);
    manager.initialize(path, rank, ranks);
    auto global = gather(snapshot(manager.domain(), rank));
    require(global.size() == 384, "Imported mesh has the wrong tetrahedron count");
    check_mesh(global);
    check_boundary(manager.domain(), global, true);
    check_manager_rounds(manager, rank);
    cudaCheckCall(cudaDeviceSynchronize());
#endif
}

using KeyFace = std::array<KeyType, 3>;

std::pair<double, std::map<std::string, double>> check_public_state(Domain& domain, int rank)
{
    auto global = gather(snapshot(domain, rank));
    std::map<KeyFace, int> incidence;
    std::map<KeyType, Point> points;
    double volume = 0;
    for (const auto& tet : global)
    {
        double det = determinant(tet.points);
        require(det > 0, "Public mesh contains an inverted tet");
        volume += det / 6;
        for (int c = 0; c < 4; ++c)
            points[tet.keys[c]] = tet.points[c];
        for (int f = 0; f < 4; ++f)
        {
            KeyFace face;
            int j = 0;
            for (int c = 0; c < 4; ++c)
                if (c != f) face[j++] = tet.keys[c];
            std::sort(face.begin(), face.end());
            ++incidence[face];
        }
    }
    std::set<KeyFace> tagged;
    std::map<std::string, double> area;
    for (const auto& name : domain.sideSetNames())
    {
        const auto& device = domain.side_set_face_keys(name);
        auto local         = read_device(device.data(), device.size());
        std::vector<KeyFace> faces;
        for (size_t i = 0; i < local.size(); i += 3)
        {
            KeyFace face{local[i], local[i + 1], local[i + 2]};
            std::sort(face.begin(), face.end());
            faces.push_back(face);
        }
        auto all = gather(faces);
        std::set<KeyFace> unique(all.begin(), all.end());
        require(!unique.empty(), "Lost a public boundary label");
        double total = 0;
        for (auto face : unique)
        {
            require(incidence.at(face) == 1 && tagged.insert(face).second,
                    "Duplicate or interior public boundary label");
            const auto& p = points.at(face[0]);
            const auto& q = points.at(face[1]);
            const auto& r = points.at(face[2]);
            Point a, b, c;
            for (int d = 0; d < 3; ++d)
            {
                a[d] = q[d] - p[d];
                b[d] = r[d] - p[d];
            }
            for (int d = 0; d < 3; ++d)
                c[d] = a[(d + 1) % 3] * b[(d + 2) % 3] - a[(d + 2) % 3] * b[(d + 1) % 3];
            total += .5 * std::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
        }
        area[name] = total;
    }
    for (const auto& [face, count] : incidence)
        require(count == (tagged.count(face) ? 1 : 2), "Hanging or unlabelled public mesh face");
    (void)check_identity(domain, rank);
    return {volume, area};
}

void check_public_manager(const std::string& path, int rank, int ranks, bool staged)
{
    auto input       = read_public_tet(path);
    const auto& mesh = input.mesh;
    Domain::DeviceCoordsTuple xyz{upload(mesh.xyz[0]), upload(mesh.xyz[1]), upload(mesh.xyz[2])};
    Domain::DeviceConnectivityTuple conn;
    std::array<std::vector<KeyType>, 4> columns;
    for (size_t e = 0; e < mesh.conn[0].size(); ++e)
        if (int(e % ranks) == rank)
            for (int c = 0; c < 4; ++c)
                columns[c].push_back(mesh.conn[c][e]);
    [&]<size_t... I>(std::index_sequence<I...>)
    { ((std::get<I>(conn) = upload(columns[I])), ...); }(std::make_index_sequence<4>{});
    Manager manager(manager_config(staged));
    manager.initialize(std::make_unique<Domain>(std::move(xyz), std::move(conn), rank, ranks, 8), rank, ranks);
    auto& domain = manager.domain();
    std::set<KeyFace> local_faces;
    for (const auto& tet : snapshot(domain, rank))
        for (int f = 0; f < 4; ++f)
        {
            KeyFace face;
            int j = 0;
            for (int c = 0; c < 4; ++c)
                if (c != f) face[j++] = tet.keys[c];
            std::sort(face.begin(), face.end());
            local_faces.insert(face);
        }
    for (const auto& [tag, faces] : input.boundary)
    {
        std::vector<KeyType> keys;
        for (const auto& face : faces)
        {
            KeyFace key;
            for (int c = 0; c < 3; ++c)
                key[c] = cstone::sfc3D<cstone::SfcKind<KeyType>>(mesh.xyz[0][face[c]], mesh.xyz[1][face[c]],
                                                                 mesh.xyz[2][face[c]], domain.getBoundingBox());
            std::sort(key.begin(), key.end());
            if (local_faces.count(key)) keys.insert(keys.end(), key.begin(), key.end());
        }
        domain.set_side_set_from_face_keys("mfem-" + std::to_string(tag), upload(keys));
    }
    const auto reference = check_public_state(domain, rank);
    auto fields          = initial_fields(domain);
    double bound         = 1.e-12;
    for (int round = 0; round < 3; ++round)
    {
        auto& before = manager.domain();
        auto local   = snapshot(before, rank);
        std::vector<double> error(before.getElementCount(), 0);
        for (size_t e = 0; e < local.size(); ++e)
        {
            Point center{};
            for (const auto& point : local[e].points)
                for (int d = 0; d < 3; ++d)
                    center[d] += point[d] / 4;
            error[before.startIndex() + e] =
                1. / (1 + center[0] * center[0] + center[1] * center[1] + center[2] * center[2]);
        }
        auto device_error = upload(error);
        std::array<Values, 2> next;
        poison_ghost_fields(before, fields);
        auto stats = manager.adaptMeshMultiField(device_error.data(), {fields[0].data(), fields[1].data()},
                                                 {&next[0], &next[1]});
        require(stats.elementsRefined > 0, "Public LEPP case did not refine");
        auto& after   = manager.domain();
        auto measured = check_public_state(after, rank);
        require(std::abs(measured.first - reference.first) <= 1.e-4 * reference.first, "Public mesh volume changed");
        require(measured.second.size() == reference.second.size(), "Public labels changed");
        for (const auto& [name, area] : reference.second)
            require(std::abs(measured.second.at(name) - area) <= 1.e-4 * area, "Public boundary label area changed");
        bound += affine_quantization_bound(after);
        auto current = current_nodes(after);
        check_fields(after.getNodeX().data(), after.getNodeY().data(), after.getNodeZ().data(), after.getNodeCount(),
                     next, bound, &current);
        fields = std::move(next);
        if (rank == 0)
            std::cout << "[public-lepp] round=" << round << " cells=" << stats.elementsAfterGlobal
                      << " volume=" << measured.first << " affine-bound=" << bound << " PASS\n";
    }
}

} // namespace amr_test
