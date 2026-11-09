#include "tet_reference.hpp"
#include "tet_gpu_checks.hpp"
#include "backend/distributed/unstructured/amr/mars_amr_tet_bisection.hpp"

namespace amr_test
{
void check_dependent_device_lepp(int rank, int ranks, bool staged)
{
    using Execution = mars::amr::LeppDeviceExecution;
    std::vector<RealType> x{0, 2, 1, 1, 1}, y{0, 0, .5, .2, .2}, z{0, 0, 0, .3, -3};
    mars::amr::TetNodeKey<KeyType, RealType> encode{cstone::Box<RealType>(-4, 4)};
    std::vector<KeyType> keys(5);
    for (int i = 0; i < 5; ++i)
        keys[i] = encode(x[i], y[i], z[i]);
    std::array<std::vector<KeyType>, 4> full_conn{{{0, 0}, {1, 2}, {2, 1}, {3, 4}}}, local_conn;
    std::vector<uint8_t> marks;
    std::vector<mars::amr::LeppTaggedFace> full_faces, local_faces;
    for (int e = 0; e < 2; ++e)
    {
        for (int f = 0; f < 3; ++f)
            full_faces.push_back({uint64_t(e), f, 3 * e + f});
        if (e % ranks != rank) continue;
        size_t local_e = marks.size();
        for (int c = 0; c < 4; ++c)
            local_conn[c].push_back(full_conn[c][e]);
        for (int f = 0; f < 3; ++f)
            local_faces.push_back({local_e, f, 3 * e + f});
        marks.push_back(e == 0);
    }
    auto d_x = upload(x), d_y = upload(y), d_z = upload(z);
    auto d_keys  = upload(keys);
    auto d_marks = upload(marks);
    std::array<cstone::DeviceVector<KeyType>, 4> conn;
    for (int c = 0; c < 4; ++c)
        conn[c] = upload(local_conn[c]);
    Execution::Vector<mars::amr::LeppTaggedFace> faces(local_faces.size());
    Execution::upload(local_faces.data(), local_faces.size(), Execution::data(faces));
    mars::amr::TetLeppView<KeyType, RealType> view{{conn[0].data(), conn[1].data(), conn[2].data(), conn[3].data()},
                                                   d_keys.data(),
                                                   d_x.data(),
                                                   d_y.data(),
                                                   d_z.data(),
                                                   d_marks.data(),
                                                   marks.size(),
                                                   keys.size()};
    auto result = mars::amr::refine_tet_lepp<Execution>(
        view, encode, [ranks] __host__ __device__(KeyType key) { return int(((key >> 32) * uint64_t(ranks)) >> 31); },
        MPI_COMM_WORLD, staged, faces);
    std::vector<uint8_t> all_marks{1, 0};
    mars::amr::TetLeppView<KeyType, RealType> host_view{
        {full_conn[0].data(), full_conn[1].data(), full_conn[2].data(), full_conn[3].data()},
        keys.data(),
        x.data(),
        y.data(),
        z.data(),
        all_marks.data(),
        2,
        5};
    auto reference = refine_reference(host_view, encode.box, full_faces);
    std::array<Values, 2> transferred;
    std::vector<RealType> initial(5);
    for (int f = 0; f < 2; ++f)
    {
        for (int i = 0; i < 5; ++i)
            initial[i] = f ? affine({x[i], y[i], z[i]}) : 2.5;
        auto old = upload(initial);
        transferred[f] =
            mars::amr::transfer_lepp_field<Execution>(old.data(), 5, result.parent_edges.data(),
                                                      result.parent_edges.size(), result.birth_offsets, MPI_COMM_WORLD);
    }
    check_fields(result.x.data(), result.y.data(), result.z.data(), result.keys.size(), transferred, 1.e-12);
    int dependent = 0;
    for (uint64_t edge : read_device(result.parent_edges.data(), result.parent_edges.size()))
        dependent |= (edge & UINT64_C(0xffffffff)) >= 5 || (edge >> 32) >= 5;
    MPI_Allreduce(MPI_IN_PLACE, &dependent, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    require(dependent, "Device fixture did not exercise successive midpoint generations");
    auto local  = snapshot({result.conn[0].data(), result.conn[1].data(), result.conn[2].data(), result.conn[3].data()},
                           0, result.conn[0].size(), result.x.data(), result.y.data(), result.z.data(),
                           result.keys.size(), nullptr, rank);
    auto global = gather(local);
    std::vector<std::array<Point, 4>> expected, actual;
    for (size_t e = 0; e < reference.conn[0].size(); ++e)
    {
        std::array<Point, 4> points;
        for (int c = 0; c < 4; ++c)
        {
            auto i    = reference.conn[c][e];
            points[c] = {reference.x[i], reference.y[i], reference.z[i]};
        }
        require(determinant(points) > 0, "Reference LEPP inverted a tetrahedron");
        std::sort(points.begin(), points.end());
        expected.push_back(points);
    }
    for (auto tet : global)
    {
        require(determinant(tet.points) > 0, "Device LEPP inverted a tetrahedron");
        std::sort(tet.points.begin(), tet.points.end());
        actual.push_back(tet.points);
    }
    std::sort(expected.begin(), expected.end());
    std::sort(actual.begin(), actual.end());
    require(actual == expected, "Successive device LEPP generations differ from the CPU core");
    struct TaggedFace
    {
        Face points;
        int tag;
    };
    std::vector<TaggedFace> local_tags;
    for (auto face : read_device(Execution::data(result.tagged_faces), result.tagged_faces.size()))
        local_tags.push_back({tet_face(local.at(face.element), face.face), face.tag});
    std::set<std::pair<Face, int>> actual_tags, expected_tags;
    for (const auto& tagged : gather(local_tags))
        require(actual_tags.emplace(tagged.points, tagged.tag).second, "Duplicate tagged output face");
    for (auto face : reference.tagged_faces)
    {
        Tet tet;
        for (int c = 0; c < 4; ++c)
        {
            auto i        = reference.conn[c][face.element];
            tet.points[c] = {reference.x[i], reference.y[i], reference.z[i]};
        }
        expected_tags.emplace(tet_face(tet, face.face), face.tag);
    }
    require(actual_tags == expected_tags, "Dependent refinement lost boundary ancestry");
    if (rank == 0) std::cout << "Dependent device LEPP: CPU parity, empty ranks, fields and boundary ancestry PASS\n";
}

} // namespace amr_test
