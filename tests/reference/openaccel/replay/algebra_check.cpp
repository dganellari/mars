// GPT/Codex, 2026-09-21. Independent invariants beyond the captured channel states.
#include "mars_segregated_tet_interior.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace mars::segregated;

int checks = 0;
void close(double a, double b)
{
    ++checks;
    if (!std::isfinite(a) || !std::isfinite(b)
        || std::abs(a-b) > 2e-9*std::max({1., std::abs(a), std::abs(b)}))
        throw std::runtime_error("algebra check failed");
}

int main()
{
    TetInteriorInput x{};
    const double coordinates[] = {0,0,0, 1,0,0, 0,1,0, 0,0,1};
    const double grad[] = {-1,-1,-1, 1,0,0, 0,1,0, 0,0,1};
    const int edges[] = {0,1, 0,2, 0,3, 1,2, 1,3, 2,3};
    std::copy(coordinates, coordinates+12, x.coordinates);
    std::copy(edges, edges+12, x.edges);
    for (int n = 0; n < 4; ++n) {
        x.density[n] = 2;
        x.viscosity[n] = 0.2+0.1*n;
        x.pressure[n] = 1+2*coordinates[3*n]-3*coordinates[3*n+1]+5*coordinates[3*n+2];
        for (int j = 0; j < 3; ++j) {
            x.velocity[3*n+j] = 0.1*(1+n+2*j);
            x.pressure_gradient[3*n+j] = j == 0 ? 2 : (j == 1 ? -3 : 5);
            x.influence_lhs[3*n+j] = x.influence_rhs[3*n+j] = 0.3+0.04*n+0.2*j;
            x.velocity_blend[3*n+j] = 0.4;
            for (int k = 0; k < 3; ++k) x.velocity_gradient[9*n+3*j+k] = 0.2*(n+j-k);
        }
    }
    for (int s = 0; s < 6; ++s) {
        x.stored_flux[s] = s%2 ? -0.2*(s+1) : 0.3*(s+1);
        for (int n = 0; n < 4; ++n) {
            x.velocity_shape[4*s+n] = 0.25;
            x.coordinate_shape[4*s+n] = 0.25;
        }
        for (int j = 0; j < 12; ++j) x.shape_gradient[12*s+j] = grad[j];
        for (int j = 0; j < 3; ++j) x.area[3*s+j] = 0.01*(s+1)*(j-1.5);
    }
    TetInteriorOutput base{}, changed{};
    tet_interior(x, base);
    // Affine p with its exact reconstructed gradient has no Rhie-Chow defect.
    for (int s = 0; s < 6; ++s) {
        double advective = 0;
        for (int n = 0; n < 4; ++n)
            for (int j = 0; j < 3; ++j) advective += 2*0.25*x.velocity[3*n+j]*x.area[3*s+j];
        close(base.flux[s], advective);
    }
    // Freeze reconstruction and coefficients: d(rhs)/dp = -A for SIMPLE dl=dr.
    const double epsilon = 1e-6;
    for (int n = 0; n < 4; ++n) {
        auto y = x; y.pressure[n] += epsilon;
        tet_interior(y, changed);
        for (int r = 0; r < 4; ++r) close((changed.rhs[r]-base.rhs[r])/epsilon, -base.lhs[4*r+n]);
    }
    // SIMPLEC's LHS coefficient must not leak into the RHS flux.
    auto y = x;
    for (double& d : y.influence_lhs) d *= 3;
    tet_interior(y, changed);
    for (int i = 0; i < 16; ++i) close(changed.lhs[i], 3*base.lhs[i]);
    for (int s = 0; s < 6; ++s) close(changed.flux[s], base.flux[s]);
    for (int stage = 0; stage < 2; ++stage) {
        x.stage = stage; tet_interior(x, base);
        const int components = stage ? 3 : 1, rows = 4*components;
        for (int c = 0; c < components; ++c) {
            double sum = 0;
            for (int n = 0; n < 4; ++n) sum += base.rhs[n*components+c];
            close(sum, 0);
            for (int j = 0; j < rows; ++j) {
                sum = 0;
                for (int n = 0; n < 4; ++n) sum += base.lhs[(n*components+c)*rows+j];
                close(sum, 0);
            }
        }
        if (stage) for (int j = 0; j < 12; ++j) {
            y = x; y.velocity[j] += epsilon;
            tet_interior(y, changed);
            for (int r = 0; r < 12; ++r) close((changed.rhs[r]-base.rhs[r])/epsilon, -base.lhs[12*r+j]);
        }
    }
    // Single unit-area face: transpose stress creates an off-component entry -mu.
    x = TetInteriorInput{}; x.stage = 1;
    std::copy(edges, edges+12, x.edges);
    for (int n = 0; n < 4; ++n) { x.viscosity[n] = 2; x.velocity_shape[n] = 0.25; }
    std::copy(grad, grad+12, x.shape_gradient);
    x.area[0] = 1;
    tet_interior(x, base);
    close(base.lhs[12*1+3*2], -2);
    close(base.lhs[12*4+3*2], 2);
    std::cout << "PASS: " << checks << " interior algebra checks\n";
}
