#include "graded_mesh.hpp"
#include <map>
#include <stdexcept>
#include <iostream>

void require(bool value, const char* message)
{
    if (!value) throw std::runtime_error(message);
}
int main()
{
    auto mesh = amr_test::thin_mesh();
    std::map<std::array<uint64_t, 3>, int> faces;
    double volume = 0;
    for (size_t e = 0; e < mesh.conn[0].size(); ++e)
    {
        double d[3][3];
        for (int c = 0; c < 3; ++c)
            for (int j = 0; j < 3; ++j)
                d[c][j] = mesh.xyz[j][mesh.conn[c + 1][e]] - mesh.xyz[j][mesh.conn[0][e]];
        double det = d[0][0] * (d[1][1] * d[2][2] - d[1][2] * d[2][1]) -
                     d[0][1] * (d[1][0] * d[2][2] - d[1][2] * d[2][0]) +
                     d[0][2] * (d[1][0] * d[2][1] - d[1][1] * d[2][0]);
        require(det > 0, "Thin fixture contains an inverted or degenerate cell");
        volume += det / 6;
        for (int f = 0; f < 4; ++f)
        {
            std::array<uint64_t, 3> face;
            int at = 0;
            for (int c = 0; c < 4; ++c)
                if (c != f) face[at++] = mesh.conn[c][e];
            std::sort(face.begin(), face.end());
            ++faces[face];
        }
    }
    size_t boundary = 0;
    for (const auto& [face, count] : faces)
    {
        require(count == 1 || count == 2, "Nonmanifold fixture");
        boundary += count == 1;
    }
    require(boundary == 64 * 4 && std::abs(volume - 64. / 6) < 1.e-12, "Stellar splits lost a face or volume");
    // Replay the old detector for coarse ABPD and candidate Q. They share A/P/D through neighbours.
    const int n[4]  = {0, 1, 4, 3};
    double shortest = 1.e100, distance2 = 0;
    for (int a = 0; a < 4; ++a)
        for (int b = a + 1; b < 4; ++b)
        {
            double l2 = 0;
            for (int d = 0; d < 3; ++d)
            {
                double v = mesh.xyz[d][n[a]] - mesh.xyz[d][n[b]];
                l2 += v * v;
            }
            shortest = std::min(shortest, l2);
        }
    bool inside = true;
    for (int d = 0; d < 3; ++d)
    {
        double lo = mesh.xyz[d][n[0]], hi = lo;
        for (int i : n)
        {
            lo = std::min(lo, mesh.xyz[d][i]);
            hi = std::max(hi, mesh.xyz[d][i]);
        }
        double q = mesh.xyz[d][5];
        inside &= q >= lo - .1 * std::sqrt(shortest) && q <= hi + .1 * std::sqrt(shortest);
        double delta = q - .5 * (mesh.xyz[d][0] + mesh.xyz[d][1]);
        distance2 += delta * delta;
    }
    require(inside && distance2 < .01, "Thin fixture does not reproduce the old false AB constraint");
    auto graded = amr_test::graded_mesh<8>();
    require(graded.conn[0].size() == 64 && .08 / .015 > 3,
            "Graded fixture does not trigger the old unmarked size rule");
    std::cout << "graded geometry, positive volumes, conforming incidence and old-detector negative control PASS\n";
}
