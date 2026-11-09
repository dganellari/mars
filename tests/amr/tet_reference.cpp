#include "tet_reference.hpp"

namespace amr_test
{
void require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

RealType determinant(const std::array<Point, 4>& p)
{
    Point a{}, b{}, c{};
    for (int d = 0; d < 3; ++d)
    {
        a[d] = p[1][d] - p[0][d];
        b[d] = p[2][d] - p[0][d];
        c[d] = p[3][d] - p[0][d];
    }
    return a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
}

RealType check_mesh(const std::vector<Tet>& global)
{
    using Face = std::array<Point, 3>;
    std::map<Face, int> faces;
    std::set<std::array<Point, 4>> cells;
    RealType volume = 0;
    for (const auto& tet : global)
    {
        for (const auto& point : tet.points)
            for (RealType value : point)
                require(std::isfinite(value), "Non-finite mesh coordinate");
        RealType det = determinant(tet.points);
        require(det > 1.e-14, "Inverted or degenerate tetrahedron");
        volume += det / 6;
        auto canonical = tet.points;
        std::sort(canonical.begin(), canonical.end());
        require(cells.insert(canonical).second, "An element has more than one owner");
        for (int opposite = 0; opposite < 4; ++opposite)
        {
            Face face;
            int j = 0;
            for (int c = 0; c < 4; ++c)
                if (c != opposite) face[j++] = tet.points[c];
            std::sort(face.begin(), face.end());
            ++faces[face];
        }
    }
    for (const auto& [face, count] : faces)
    {
        bool boundary = false;
        for (int d = 0; d < 3; ++d)
            for (RealType value : {0.0, 1.0})
                boundary |= std::abs(face[0][d] - value) < 2.e-6 && std::abs(face[1][d] - value) < 2.e-6 &&
                            std::abs(face[2][d] - value) < 2.e-6;
        require(count == (boundary ? 1 : 2), "Unmatched interior face, overlapping cells or nonmanifold face");
    }
    require(std::abs(volume - 1) < 1.e-5, "Cube volume changed");
    return volume;
}

RealType affine(const Point& p) { return 1 + 2 * p[0] - 3 * p[1] + 0.5 * p[2]; }

Edge longest_edge(const Tet& tet)
{
    Edge best{};
    RealType length = -1;
    for (int a = 0; a < 4; ++a)
        for (int b = a + 1; b < 4; ++b)
        {
            RealType dx      = tet.points[a][0] - tet.points[b][0];
            RealType dy      = tet.points[a][1] - tet.points[b][1];
            RealType dz      = tet.points[a][2] - tet.points[b][2];
            RealType squared = dx * dx + dy * dy + dz * dz;
            Edge edge        = {tet.keys[a], tet.keys[b]};
            std::sort(edge.begin(), edge.end());
            if (squared > length || (squared == length && edge < best))
            {
                best   = edge;
                length = squared;
            }
        }
    return best;
}

Face tet_face(const Tet& tet, int opposite)
{
    Face face;
    for (int c = 0, j = 0; c < 4; ++c)
        if (c != opposite) face[j++] = tet.points[c];
    std::sort(face.begin(), face.end());
    return face;
}

int boundary_plane(const Face& face)
{
    for (int d = 0; d < 3; ++d)
        for (int side = 0; side < 2; ++side)
            if (std::abs(face[0][d] - side) < 2.e-6 && std::abs(face[1][d] - side) < 2.e-6 &&
                std::abs(face[2][d] - side) < 2.e-6)
                return 2 * d + side;
    return -1;
}

HostResult refine_reference(mars::amr::TetLeppView<KeyType, RealType> view,
                            const cstone::Box<RealType>& box,
                            const std::vector<mars::amr::LeppTaggedFace>& faces)
{
    auto key = [box](RealType x, RealType y, RealType z)
    { return KeyType(cstone::sfc3D<cstone::SfcKind<KeyType>>(x, y, z, box)); };
    return mars::amr::refine_tet_lepp<mars::amr::LeppHostExecution>(
        view, key, [](KeyType) { return 0; }, MPI_COMM_SELF, false, faces);
}
} // namespace amr_test
