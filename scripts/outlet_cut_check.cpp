// GPT/Codex, 2026-09-11. Public Kuhn cube, independent of any application mesh.
#include "../backend/distributed/unstructured/fem/mars_outlet_flux.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

int main()
{
    int checks = 0;
    auto near = [&](double actual, double expected) {
        ++checks;
        if (std::abs(actual-expected) > 1e-13) throw std::runtime_error("signed cut gate failed");
    };
    for (double left : {-1., 0., 1.})
        for (double right : {-1., 0., 1.})
            for (double cut : {-1., 0., 1.})
            {
                const int orientation = outlet_cut_weight(left, right, cut);
                near(orientation, -outlet_cut_weight(right, left, cut));
                // Independent row scatter, including endpoint ties and same-side cancellation.
                const double rows = (left <= cut ? 2.7 : 0.) + (right <= cut ? -2.7 : 0.);
                near(orientation*2.7, rows);
            }
    const double cube[8][3] = {{0,0,0},{1,0,0},{1,1,0},{0,1,0},
                               {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    const int tets[6][4] = {{0,1,2,6},{0,2,3,6},{0,3,7,6},{0,7,4,6},{0,4,5,6},{0,5,1,6}};
    const int edges[6][2] = {{0,1},{1,2},{0,2},{0,3},{1,3},{2,3}};
    for (int axis = 0; axis < 3; ++axis)
        for (double cut : {.1,.5,.9})
        {
            double signed_flux = 0, old_flux = 0;
            for (const auto& tet : tets)
            {
                double coordinates[4][3], gradient[4][3], det;
                for (int i = 0; i < 4; ++i)
                    for (int d = 0; d < 3; ++d) coordinates[i][d] = cube[tet[i]][d];
                outlet_tet_gradient(coordinates, det, gradient);
                for (const auto& edge : edges)
                {
                    const int l = edge[0], r = edge[1];
                    const double area = det*(gradient[r][axis]-gradient[l][axis])/24;
                    const double flux = .1*area;
                    const int weight = outlet_cut_weight(coordinates[l][axis], coordinates[r][axis], cut);
                    signed_flux += weight*flux;
                    if (weight != 0) old_flux += flux;
                }
            }
            near(signed_flux, .1);
            if (axis == 0) near(old_flux, .1*5/6);
        }
    std::cout << "PASS: " << checks << " signed cut checks; old x-flux=1/12, corrected=1/10\n";
}
