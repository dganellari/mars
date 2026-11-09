#include "tet_gpu_checks.hpp"
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

void check_boundary(Domain& domain, const std::vector<Tet>& global)
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

void check_manager(std::unique_ptr<Domain> domain, int rank, int ranks, bool staged)
{
    Manager::Config config;
    config.tet_refinement            = mars::amr::TetRefinement::ConformingBisection;
    config.stage_refinement_messages = staged;
    config.strategy                  = mars::amr::MarkingStrategy::Doerfler;
    // The existing marker uses focus-tree depth, which need not equal the number of bisections.
    config.maxLevels       = 32;
    config.bucketSize      = 8;
    config.refineFraction  = 0.5;
    config.coarsenFraction = 0;
    Manager manager(config);
    manager.initialize(std::move(domain), rank, ranks);
    auto fields = initial_fields(manager.domain());
    (void)check_identity(manager.domain(), rank);

    label_boundary(manager.domain(), rank);
    check_boundary(manager.domain(), gather(snapshot(manager.domain(), rank)));

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

} // namespace amr_test
