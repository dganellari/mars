#include "backend/distributed/unstructured/amr/mars_amr_lepp.hpp"
#include "backend/distributed/unstructured/amr/mars_amr_lepp_boundary.hpp"
#include "backend/distributed/unstructured/amr/mars_amr_field_migration.hpp"

#include <iostream>
#include "public_tet_mesh.hpp"
#include <map>
#include <set>
#include <string>

namespace traffic
{
uint64_t sends = 0;
}

extern "C" int
MPI_Isend(const void* buffer, int count, MPI_Datatype type, int rank, int tag, MPI_Comm comm, MPI_Request* request)
{
    ++traffic::sends;
    return PMPI_Isend(buffer, count, type, rank, tag, comm, request);
}

extern "C" int
MPI_Issend(const void* buffer, int count, MPI_Datatype type, int rank, int tag, MPI_Comm comm, MPI_Request* request)
{
    ++traffic::sends;
    return PMPI_Issend(buffer, count, type, rank, tag, comm, request);
}

namespace
{
using KeyType   = uint64_t;
using Point     = std::array<double, 3>;
using Tet       = std::array<Point, 4>;
using Face      = std::array<Point, 3>;
using Execution = mars::amr::LeppHostExecution;
using Result    = mars::amr::TetLeppResult<Execution, KeyType, double>;

struct NodeKey
{
    KeyType operator()(double x, double y, double z) const
    {
        auto fixed = [](double v) { return KeyType((v + 8) * 65536 + .5); };
        return fixed(x) | (fixed(y) << 21) | (fixed(z) << 42);
    }
};

auto node_owner(MPI_Comm comm)
{
    int ranks;
    MPI_Comm_size(comm, &ranks);
    return [ranks](KeyType key) { return int(((key >> 32) * uint64_t(ranks)) >> 31); };
}

void require(bool valid, const char* message)
{
    if (!valid) throw std::runtime_error(message);
}

double volume(const Tet& tet)
{
    Point a{}, b{}, c{};
    for (int d = 0; d < 3; ++d)
    {
        a[d] = tet[1][d] - tet[0][d];
        b[d] = tet[2][d] - tet[0][d];
        c[d] = tet[3][d] - tet[0][d];
    }
    return (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) +
            a[2] * (b[0] * c[1] - b[1] * c[0])) /
           6;
}

double area(const Face& face)
{
    Point a{}, b{}, cross{};
    for (int d = 0; d < 3; ++d)
    {
        a[d] = face[1][d] - face[0][d];
        b[d] = face[2][d] - face[0][d];
    }
    cross = {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
    return .5 * std::sqrt(cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]);
}

Face face_of(const Tet& tet, int opposite)
{
    Face face;
    for (int i = 0, j = 0; i < 4; ++i)
        if (i != opposite) face[j++] = tet[i];
    std::sort(face.begin(), face.end());
    return face;
}

struct TaggedTet
{
    Tet points;
    int tags[4]    = {-1, -1, -1, -1};
    uint8_t marked = 0;
};

template<class T>
std::vector<T> gather(const std::vector<T>& local, MPI_Comm comm)
{
    int ranks;
    MPI_Comm_size(comm, &ranks);
    int bytes = int(local.size() * sizeof(T));
    std::vector<int> counts(ranks), offsets(ranks);
    MPI_Allgather(&bytes, 1, MPI_INT, counts.data(), 1, MPI_INT, comm);
    std::exclusive_scan(counts.begin(), counts.end(), offsets.begin(), 0);
    size_t total = std::accumulate(counts.begin(), counts.end(), size_t(0));
    std::vector<T> result(total / sizeof(T));
    MPI_Allgatherv(local.data(), bytes, MPI_BYTE, result.data(), counts.data(), offsets.data(), MPI_BYTE, comm);
    return result;
}

std::map<Face, int> surface(const std::vector<TaggedTet>& cells)
{
    std::map<Face, int> counts;
    for (const auto& cell : cells)
        for (int f = 0; f < 4; ++f)
            ++counts[face_of(cell.points, f)];
    return counts;
}

void label_boundary(std::vector<TaggedTet>& cells)
{
    auto faces = surface(cells);
    std::map<Face, int> labels;
    int tag = 0;
    for (const auto& [face, count] : faces)
        if (count == 1) labels[face] = tag++;
    for (auto& cell : cells)
        for (int f = 0; f < 4; ++f)
        {
            auto found   = labels.find(face_of(cell.points, f));
            cell.tags[f] = found == labels.end() ? -1 : found->second;
        }
}

void check_geometry(const std::vector<TaggedTet>& before, const std::vector<TaggedTet>& after)
{
    double old_volume = 0, new_volume = 0;
    std::map<int, double> old_area, new_area;
    for (const auto& cell : before)
    {
        old_volume += volume(cell.points);
        for (int f = 0; f < 4; ++f)
            if (cell.tags[f] >= 0) old_area[cell.tags[f]] += area(face_of(cell.points, f));
    }
    auto counts = surface(after);
    std::set<Face> tagged;
    for (const auto& cell : after)
    {
        require(volume(cell.points) > 0, "LEPP produced a nonpositive tet");
        new_volume += volume(cell.points);
        for (int f = 0; f < 4; ++f)
            if (cell.tags[f] >= 0)
            {
                auto face = face_of(cell.points, f);
                require(counts.at(face) == 1 && tagged.insert(face).second,
                        "Boundary tag on an interior/duplicate face");
                new_area[cell.tags[f]] += area(face);
            }
    }
    require(std::abs(old_volume - new_volume) < 1.e-11 * old_volume, "LEPP volume changed");
    require(old_area.size() == new_area.size(), "Lost boundary label");
    for (const auto& [tag, a] : old_area)
        require(std::abs(a - new_area.at(tag)) < 1.e-11 * a, "Boundary area/label changed");
    for (const auto& [face, count] : counts)
        require(count == (tagged.count(face) ? 1 : 2), "Hanging or nonmanifold interior face");
}

std::vector<Tet> canonical(const std::vector<TaggedTet>& mesh)
{
    std::vector<Tet> tets;
    for (const auto& cell : mesh)
    {
        auto tet = cell.points;
        std::sort(tet.begin(), tet.end());
        tets.push_back(tet);
    }
    std::sort(tets.begin(), tets.end());
    return tets;
}

struct Input
{
    std::array<std::vector<KeyType>, 4> conn;
    std::vector<KeyType> keys;
    std::vector<double> x, y, z;
    std::vector<uint8_t> marks;
    std::vector<mars::amr::LeppTaggedFace> faces;
    mars::amr::TetLeppView<KeyType, double> view() const
    {
        return {{conn[0].data(), conn[1].data(), conn[2].data(), conn[3].data()},
                keys.data(),
                x.data(),
                y.data(),
                z.data(),
                marks.data(),
                marks.size(),
                keys.size()};
    }
};

Input partition(const std::vector<TaggedTet>& mesh, int rank, int ranks, int layout)
{
    std::vector<TaggedTet> owned;
    for (size_t i = 0; i < mesh.size(); ++i)
    {
        int owner = layout == 1 ? 0 : int((i + size_t(layout)) % ranks);
        if (owner == rank) owned.push_back(mesh[i]);
    }
    std::set<Point> points;
    for (const auto& cell : owned)
        for (const auto& point : cell.points)
            points.insert(point);
    std::vector<Point> ordered(points.begin(), points.end());
    if (rank % 2 == 0) std::reverse(ordered.begin(), ordered.end());
    std::map<Point, KeyType> ids;
    Input input;
    for (const auto& point : ordered)
    {
        ids[point] = input.keys.size();
        input.keys.push_back(NodeKey{}(point[0], point[1], point[2]));
        input.x.push_back(point[0]);
        input.y.push_back(point[1]);
        input.z.push_back(point[2]);
    }
    for (const auto& cell : owned)
    {
        size_t e = input.marks.size();
        for (int c = 0; c < 4; ++c)
        {
            input.conn[c].push_back(ids.at(cell.points[c]));
            if (cell.tags[c] >= 0) input.faces.push_back({e, c, cell.tags[c]});
        }
        input.marks.push_back(cell.marked);
    }
    return input;
}

std::vector<TaggedTet>
run(const std::vector<TaggedTet>& mesh, MPI_Comm comm, bool staged, int layout, bool* dependent_midpoint = nullptr)
{
    int rank, ranks;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &ranks);
    Input input = partition(mesh, rank, ranks, layout);
    auto result =
        mars::amr::refine_tet_lepp<Execution>(input.view(), NodeKey{}, node_owner(comm), comm, staged, input.faces);
    for (bool constant : {true, false})
    {
        std::vector<double> values(input.keys.size());
        for (size_t i = 0; i < values.size(); ++i)
            values[i] = constant ? 3. : 1 + 2 * input.x[i] - 3 * input.y[i] + .5 * input.z[i];
        auto transferred =
            mars::amr::transfer_lepp_field<Execution>(values.data(), input.keys.size(), result.parent_edges.data(),
                                                      result.parent_edges.size(), result.birth_offsets, comm);
        for (size_t i = 0; i < transferred.size(); ++i)
        {
            double expected = constant ? 3. : 1 + 2 * result.x[i] - 3 * result.y[i] + .5 * result.z[i];
            require(std::abs(transferred[i] - expected) < 1.e-12, "Dependent P1 transfer failed");
        }
    }
    int dependent = 0;
    for (uint64_t edge : result.parent_edges)
        dependent |= (edge & UINT64_C(0xffffffff)) >= input.keys.size() || (edge >> 32) >= input.keys.size();
    MPI_Allreduce(MPI_IN_PLACE, &dependent, 1, MPI_INT, MPI_MAX, comm);
    if (dependent_midpoint) *dependent_midpoint |= dependent;
    std::vector<TaggedTet> output(result.conn[0].size());
    for (size_t e = 0; e < output.size(); ++e)
        for (int c = 0; c < 4; ++c)
        {
            size_t i            = result.conn[c][e];
            output[e].points[c] = {result.x[i], result.y[i], result.z[i]};
        }
    for (auto face : result.tagged_faces)
        output[face.element].tags[face.face] = face.tag;
    return gather(output, comm);
}

std::vector<TaggedTet> witness()
{
    Point a{0, 0, 0}, b{2, 0, 0}, c{1, .5, 0}, d{1, .2, .3}, e{-.5, 0, -1};
    std::vector<TaggedTet> mesh(2);
    mesh[0].points = {a, b, c, d};
    mesh[0].marked = 1;
    mesh[1].points = {a, c, b, e};
    label_boundary(mesh);
    return mesh;
}

std::vector<TaggedTet> edge_star(int valence)
{
    std::vector<TaggedTet> mesh(valence);
    for (int i = 0; i < valence; ++i)
    {
        double a = 2 * std::acos(-1.) * i / valence, b = 2 * std::acos(-1.) * (i + 1) / valence;
        Point p{std::cos(a), std::sin(a), 0}, q{std::cos(b), std::sin(b), 0};
        if (i + 1 == valence) q = {1, 0, 0};
        mesh[i].points = {Point{0, 0, -2}, Point{0, 0, 2}, p, q};
        if (volume(mesh[i].points) < 0) std::swap(mesh[i].points[0], mesh[i].points[1]);
    }
    mesh[0].marked = 1;
    label_boundary(mesh);
    return mesh;
}

std::vector<TaggedTet> cube()
{
    const Point points[8] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    const int tets[6][4]  = {{0, 1, 2, 6}, {0, 2, 3, 6}, {0, 3, 7, 6}, {0, 7, 4, 6}, {0, 4, 5, 6}, {0, 5, 1, 6}};
    std::vector<TaggedTet> mesh(6);
    for (int e = 0; e < 6; ++e)
    {
        for (int c = 0; c < 4; ++c)
            mesh[e].points[c] = points[tets[e][c]];
        mesh[e].marked = e == 0;
    }
    label_boundary(mesh);
    return mesh;
}

void check_case(std::vector<TaggedTet> mesh, bool staged, int cycles, bool require_dependency)
{
    bool dependent = false;
    for (int cycle = 0; cycle < cycles; ++cycle)
    {
        auto reference = run(mesh, MPI_COMM_SELF, staged, 0, &dependent);
        for (int layout : {0, 1, 2})
        {
            auto distributed = run(mesh, MPI_COMM_WORLD, staged, layout, &dependent);
            require(canonical(distributed) == canonical(reference), "Distributed LEPP topology differs from one rank");
            check_geometry(mesh, distributed);
        }
        mesh = std::move(reference);
        for (size_t e = 0; e < mesh.size(); ++e)
            mesh[e].marked = ((e + cycle) % 3) == 0;
    }
    if (require_dependency) require(dependent, "Fixture did not exercise a new midpoint as a parent");
}
template<class Function>
void reject_collectively(Function function, MPI_Comm comm)
{
    int caught = 0, ranks;
    try
    {
        function();
    }
    catch (const std::runtime_error&)
    {
        caught = 1;
    }
    MPI_Comm_size(comm, &ranks);
    MPI_Allreduce(MPI_IN_PLACE, &caught, 1, MPI_INT, MPI_SUM, comm);
    require(caught == ranks, "Malformed input was not rejected on every rank");
}

void check_face_migration(MPI_Comm comm, bool staged)
{
    int rank, ranks;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &ranks);
    auto mesh        = witness();
    auto input       = partition(mesh, rank, ranks, 2);
    using FaceRecord = mars::amr::LeppFaceRecord<KeyType>;
    using FaceKey    = mars::amr::LeppFaceKey<KeyType>;
    std::vector<FaceRecord> sources;
    if (rank == 0)
        for (const auto& cell : mesh)
            for (int f = 0; f < 4; ++f)
                if (cell.tags[f] >= 0)
                {
                    auto face = face_of(cell.points, f);
                    KeyType keys[3];
                    for (int c = 0; c < 3; ++c)
                        keys[c] = NodeKey{}(face[c][0], face[c][1], face[c][2]);
                    FaceKey key(keys[0], keys[1], keys[2]);
                    sources.push_back({key, cell.tags[f]});
                    sources.push_back({key, cell.tags[f]}); // Shared copies must not duplicate labels.
                    sources.push_back({key, 1000});         // A face can belong to more than one named set.
                }
    std::array<const KeyType*, 4> conn{input.conn[0].data(), input.conn[1].data(), input.conn[2].data(),
                                       input.conn[3].data()};
    auto matched = mars::amr::match_lepp_faces<Execution, KeyType>(sources, input.keys.data(), conn, input.marks.size(),
                                                                   node_owner(comm), comm, staged);
    using Identity = std::array<uint64_t, 3>;
    std::set<Identity> expected, actual;
    for (auto face : input.faces)
    {
        expected.insert({face.element, uint64_t(face.face), uint64_t(face.tag)});
        expected.insert({face.element, uint64_t(face.face), 1000});
    }
    for (auto face : matched)
        require(actual.insert({face.element, uint64_t(face.face), uint64_t(face.tag)}).second,
                "Duplicate face label after migration");
    require(actual == expected, "Face labels failed to follow distributed cells");
    auto refined =
        mars::amr::refine_tet_lepp<Execution>(input.view(), NodeKey{}, node_owner(comm), comm, staged, matched);
    std::map<std::pair<uint64_t, int>, std::set<int>> labels;
    for (auto face : refined.tagged_faces)
        labels[{face.element, face.face}].insert(face.tag);
    for (const auto& [face, tags] : labels)
        require(tags.size() == 2 && tags.count(1000), "Refinement lost overlapping face labels");

    if (rank == 0) sources.push_back({FaceKey(1, 2, 3), 0});
    reject_collectively(
        [&]
        {
            mars::amr::match_lepp_faces<Execution, KeyType>(sources, input.keys.data(), conn, input.marks.size(),
                                                            node_owner(comm), comm, staged);
        },
        comm);
}

void check_field_migration(MPI_Comm comm, bool staged)
{
    int rank, ranks;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &ranks);
    size_t n = size_t(4 * ranks);
    std::vector<KeyType> keys(n), sources(n);
    std::iota(keys.begin(), keys.end(), KeyType(0));
    sources = keys;
    std::reverse(keys.begin(), keys.end());
    std::vector<uint8_t> owned(n), valid(n, rank == 0);
    for (size_t i = 0; i < n; ++i)
        owned[i] = int(keys[i] % ranks) == rank;
    auto owner      = [=](KeyType key) { return int(key % ranks); };
    using Migration = mars::amr::AmrFieldMigration<Execution, KeyType, double>;
    Migration plan(sources.data(), n, valid.data(), keys.data(), owned.data(), n, owner, comm, staged);
    for (int field = 0; field < 3; ++field)
    {
        std::vector<double> values(n, -12345), result;
        if (rank == 0)
            for (size_t i = 0; i < n; ++i)
                values[i] = 3 + field * double(sources[i]);
        plan.apply(values.data(), n, result);
        for (size_t i = 0; i < n; ++i)
            require(result[i] == (owned[i] ? 3 + field * double(keys[i]) : 0),
                    "Reused field plan returned a wrong value");
    }
    if (rank == 0) keys[0] = n + ranks; // Owned by rank zero; no source supplies it.
    if (rank == 0) owned[0] = 1;
    reject_collectively(
        [&] { Migration missing(sources.data(), n, valid.data(), keys.data(), owned.data(), n, owner, comm, staged); },
        comm);
    keys = sources;
    for (size_t i = 0; i < n; ++i)
        owned[i] = int(keys[i] % ranks) == rank;
    sources.push_back(0);
    valid.push_back(rank == 0);
    Migration conflicting(sources.data(), sources.size(), valid.data(), keys.data(), owned.data(), n, owner, comm,
                          staged);
    std::vector<double> values(sources.size(), 1), result;
    values.back() = 2;
    reject_collectively([&] { conflicting.apply(values.data(), values.size(), result); }, comm);
}

void check_local_communication(MPI_Comm comm, bool staged)
{
    int rank, ranks;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &ranks);
    // Each disconnected tetrahedron lies strictly inside its coordinator's key interval.
    double z = -8 + 32 * (rank + .25) / ranks;
    double h = .01 / ranks;
    TaggedTet cell{{Point{0, 0, z}, Point{h, 0, z}, Point{0, h, z}, Point{0, 0, z + h}}};
    cell.marked = 1;
    std::vector<TaggedTet> mesh{cell};
    label_boundary(mesh);
    auto input = partition(mesh, 0, 1, 0);
    auto owner = node_owner(comm);
    for (auto key : input.keys)
        require(owner(key) == rank, "Local routing fixture escaped its interval");
    using FaceRecord = mars::amr::LeppFaceRecord<KeyType>;
    std::vector<FaceRecord> sources;
    for (int f = 0; f < 4; ++f)
    {
        KeyType keys[3];
        for (int c = 0, j = 0; c < 4; ++c)
            if (c != f) keys[j++] = input.keys[input.conn[c][0]];
        sources.push_back({{keys[0], keys[1], keys[2]}, f});
    }
    traffic::sends = 0;
    std::array<const KeyType*, 4> conn{input.conn[0].data(), input.conn[1].data(), input.conn[2].data(),
                                       input.conn[3].data()};
    auto faces =
        mars::amr::match_lepp_faces<Execution, KeyType>(sources, input.keys.data(), conn, 1, owner, comm, staged);
    auto result = mars::amr::refine_tet_lepp<Execution>(input.view(), NodeKey{}, owner, comm, staged, faces);
    mars::amr::validate_lepp_node_keys<Execution>(result.keys.data(), result.x.data(), result.y.data(), result.z.data(),
                                                  result.keys.size(), owner, comm, staged);
    require(result.conn[0].size() == 2 && result.tagged_faces.size() == 6, "Local LEPP fixture failed");
    require(traffic::sends == 0, "Local edge/node/face records entered MPI");
}

void check_sparse_plan(MPI_Comm comm, bool staged)
{
    int rank, ranks;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &ranks);
    mars::amr::AmrCommunicator private_comm(comm);
    for (bool local : {true, false})
    {
        int peer = local ? rank : (rank + 1) % ranks;
        std::vector<int> destinations{rank, peer, rank, peer};
        traffic::sends = 0;
        mars::amr::AmrSparsePlan<Execution, uint64_t> plan(destinations, private_comm, staged);
        require(plan.outgoing_peers() == size_t(peer != rank), "Sparse plan has an unexpected peer");
        require(plan.outgoing_records() == (peer == rank ? 0u : 2u), "Local records entered the remote buffer");
        require(plan.local_records() == (peer == rank ? 4u : 2u), "Local records were not retained");
        for (int round = 0; round < 3; ++round)
        {
            std::vector<uint64_t> values{uint64_t(4 * rank), uint64_t(4 * rank + 1), uint64_t(4 * rank + 2),
                                         uint64_t(4 * rank + 3)};
            auto received = plan.forward(values);
            for (auto& value : received)
                value += uint64_t(round + 7);
            auto returned = plan.reverse(received);
            for (size_t i = 0; i < values.size(); ++i)
                require(returned[i] == values[i] + uint64_t(round + 7), "Sparse plan reverse order changed");
        }
        if (peer == rank)
            require(traffic::sends == 0, "Self messages reached MPI");
        else
            require(traffic::sends > 0, "MPI send counter missed remote traffic");
    }
}

void check_public_mesh(MPI_Comm comm, bool staged, const std::string& path)
{
    auto fixture = amr_test::read_public_tet(path);
    std::vector<TaggedTet> mesh(fixture.mesh.conn[0].size());
    std::map<std::array<uint64_t, 3>, int> tags;
    for (const auto& [tag, faces] : fixture.boundary)
        for (auto face : faces)
        {
            std::sort(face.begin(), face.end());
            tags[face] = tag;
        }
    for (size_t e = 0; e < mesh.size(); ++e)
    {
        for (int c = 0; c < 4; ++c)
            for (int d = 0; d < 3; ++d)
                mesh[e].points[c][d] = fixture.mesh.xyz[d][fixture.mesh.conn[c][e]];
        for (int f = 0; f < 4; ++f)
        {
            std::array<uint64_t, 3> face;
            int at = 0;
            for (int c = 0; c < 4; ++c)
                if (c != f) face[at++] = fixture.mesh.conn[c][e];
            std::sort(face.begin(), face.end());
            auto it = tags.find(face);
            if (it != tags.end()) mesh[e].tags[f] = it->second;
        }
    }
    for (int round = 0; round < 3; ++round)
    {
        for (size_t e = 0; e < mesh.size(); ++e)
            mesh[e].marked = (e + round) % 3 == 0;
        auto after = run(mesh, comm, staged, round % 2);
        require(after.size() > mesh.size(), "Public mesh did not refine");
        check_geometry(mesh, after);
        mesh = std::move(after);
    }
}

void check_generations(MPI_Comm comm, bool staged)
{
    int rank, ranks;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &ranks);
    for (KeyType initial : {KeyType(0), KeyType(5)})
    {
        auto input = partition(witness(), rank, ranks, 1);
        std::vector<KeyType> generations(input.marks.size(), initial);
        auto view        = input.view();
        view.generations = generations.data();
        auto result      = mars::amr::refine_tet_lepp<Execution>(view, NodeKey{}, node_owner(comm), comm, staged,
                                                                 input.faces, initial + 2);
        auto all         = gather(result.generations, comm);
        std::sort(all.begin(), all.end());
        require(all == std::vector<KeyType>({initial + 1, initial + 1, initial + 1, initial + 2, initial + 2}),
                "LEPP generation does not count actual parent bisections");
        reject_collectively(
            [&]
            {
                mars::amr::refine_tet_lepp<Execution>(view, NodeKey{}, node_owner(comm), comm, staged, input.faces,
                                                      initial + 1);
            },
            comm);
    }
}

void check_rejections(MPI_Comm comm, bool staged)
{
    int rank, ranks;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &ranks);
    auto input = partition(witness(), rank, ranks, 1);
    if (rank == 0) input.conn[0][0] = input.keys.size();
    reject_collectively(
        [&] { mars::amr::refine_tet_lepp<Execution>(input.view(), NodeKey{}, node_owner(comm), comm, staged); }, comm);
    input = partition(witness(), rank, ranks, 1);
    if (rank == 0) input.faces[0].face = 4;
    reject_collectively(
        [&] {
            mars::amr::refine_tet_lepp<Execution>(input.view(), NodeKey{}, node_owner(comm), comm, staged, input.faces);
        },
        comm);
    input = partition(witness(), rank, ranks, 1);
    if (rank == 0) input.z.assign(input.z.size(), 0);
    reject_collectively(
        [&] { mars::amr::refine_tet_lepp<Execution>(input.view(), NodeKey{}, node_owner(comm), comm, staged); }, comm);
    input = partition(witness(), rank, ranks, 1);
    reject_collectively(
        [&]
        {
            mars::amr::refine_tet_lepp<Execution>(
                input.view(), [](double, double, double) { return KeyType(0); }, node_owner(comm), comm, staged);
        },
        comm);
    KeyType key = 7;
    double x = rank, y = 0, z = 0;
    if (ranks > 1)
        reject_collectively(
            [&] { mars::amr::validate_lepp_node_keys<Execution>(&key, &x, &y, &z, 1, node_owner(comm), comm, staged); },
            comm);
    std::vector<double> old{1, 2};
    std::vector<uint64_t> parents{uint64_t(2) << 32};
    std::vector<size_t> batches{2, 3};
    reject_collectively(
        [&] { mars::amr::transfer_lepp_field<Execution>(old.data(), 2, parents.data(), 1, batches, comm); }, comm);
    if (ranks > 1)
    {
        parents[0] = uint64_t(1) << 32;
        if (rank) batches.insert(batches.begin() + 1, 2);
        reject_collectively(
            [&] { mars::amr::transfer_lepp_field<Execution>(old.data(), 2, parents.data(), 1, batches, comm); }, comm);
    }
}
} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank, ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    try
    {
        char processor[MPI_MAX_PROCESSOR_NAME];
        int processor_length;
        MPI_Get_processor_name(processor, &processor_length);
        std::cout << "CPU rank=" << rank << " host=" << std::string(processor, processor_length) << '\n';
        bool staged = argc == 2 && std::string(argv[1]) == "--staged";
        require(argc == 1 || staged, "Unknown argument");
        check_sparse_plan(MPI_COMM_WORLD, staged);
        check_local_communication(MPI_COMM_WORLD, staged);
        auto initial = witness();
        auto refined = run(initial, MPI_COMM_WORLD, staged, 0);
        require(refined.size() == 5, "LEPP witness did not split the neighbour first");
        check_geometry(initial, refined);
        check_case(witness(), staged, 3, false);
        check_case(edge_star(48), staged, 2, false);
        check_case(cube(), staged, 5, false);
        auto propagation         = witness();
        propagation[1].points[3] = {1, .2, -3};
        label_boundary(propagation);
        check_case(propagation, staged, 2, true);
        check_face_migration(MPI_COMM_WORLD, staged);
        check_field_migration(MPI_COMM_WORLD, staged);
        check_public_mesh(MPI_COMM_WORLD, staged, MARS_PUBLIC_TET_MESH);
        check_generations(MPI_COMM_WORLD, staged);
        check_rejections(MPI_COMM_WORLD, staged);
        MPI_Comm split;
        MPI_Comm_split(MPI_COMM_WORLD, rank % 2, rank, &split);
        auto split_output = run(witness(), split, staged, 0);
        require(split_output.size() == 5, "LEPP used the wrong communicator");
        check_face_migration(split, staged);
        check_field_migration(split, staged);
        MPI_Comm_free(&split);
        if (rank == 0)
            std::cout << "PASS LEPP: ranks=" << ranks << " staged=" << staged
                      << " legacy witness, empty ranks, permuted identities, 48-cell star, tagged faces, dependent "
                         "transfer, repeated repartition\n";
    }
    catch (const std::exception& error)
    {
        std::cerr << "rank " << rank << ": " << error.what() << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
    MPI_Finalize();
}
