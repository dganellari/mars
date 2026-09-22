#pragma once
#include <cmath>

#if defined(__CUDACC__) || defined(__HIPCC__)
#define MARS_BOUNDARY_HD __host__ __device__
#else
#define MARS_BOUNDARY_HD
#endif

namespace mars::segregated {
// Stages: pressure inlet/outlet/wall, momentum inlet/outlet/wall.
// Tet-local indexing except momentum wall, whose block uses three face nodes.
struct BoundaryInput {
    int stage = 0;
    int face_nodes[3]{}, nearest[3]{}, opposing[3]{}, reversal[3]{};
    double area[9]{}, shape[9]{}, gradient[36]{}, velocity[12]{};
    double boundary_velocity[9]{}, viscosity[3]{}, density[3]{}, pressure[4]{};
    double pressure_gradient[12]{}, influence_lhs[9]{}, influence_rhs[9]{};
    double bc_multiplier[4]{}, stored_flux[3]{}, wall_coefficient[3]{};
};
struct BoundaryOutput {
    double lhs[144]{}, rhs[12]{}, flux[3]{};
};

// Inputs are reference workspaces after side-value substitution. No row pinning,
// trace refresh, flux relaxation, or changes to reversal flags occur here.
MARS_BOUNDARY_HD inline void boundary_block(const BoundaryInput& x, BoundaryOutput& out)
{
    out = BoundaryOutput{};
    const bool momentum = x.stage >= 3;
    const int width = momentum ? (x.stage == 5 ? 9 : 12) : 4;
    for (int sample = 0; sample < 3; ++sample) {
        const int row = x.nearest[sample];
        const double* a = x.area+3*sample;
        const double* shape = x.shape+3*sample;
        const double* grad = x.gradient+12*sample;
        // The velocity-specified momentum inlet has no reversal branch.
        if ((x.stage == 0 || x.stage == 1 || x.stage == 4) && x.reversal[sample]) continue;
        if (!momentum) {
            double rho = 0;
            for (int f = 0; f < 3; ++f) rho += shape[f]*x.density[f];
            if (x.stage != 1) {
                for (int i = 0; i < 3; ++i) out.flux[sample] += rho*x.boundary_velocity[3*sample+i]*a[i];
            } else {
                for (int node = 0; node < 4; ++node) {
                    double entry = 0;
                    for (int i = 0; i < 3; ++i) {
                        double d = 0;
                        for (int f = 0; f < 3; ++f) d += shape[f]*x.influence_lhs[3*f+i];
                        entry += -rho*d*grad[3*node+i]*a[i];
                    }
                    out.lhs[4*row+node] += entry*x.bc_multiplier[node];
                }
                for (int i = 0; i < 3; ++i) {
                    double u = 0, d = 0, dp = 0;
                    for (int f = 0; f < 3; ++f) {
                        u += shape[f]*x.velocity[3*x.face_nodes[f]+i];
                        d += shape[f]*x.influence_rhs[3*f+i];
                    }
                    for (int node = 0; node < 4; ++node) dp += grad[3*node+i]*x.pressure[node];
                    const double g = .5*(x.pressure_gradient[3*row+i]+x.pressure_gradient[3*x.opposing[sample]+i]);
                    out.flux[sample] += (rho*u-rho*d*(dp-g))*a[i];
                }
            }
            out.rhs[row] -= out.flux[sample];
            continue;
        }
        double area_squared = 0;
        for (int i = 0; i < 3; ++i) area_squared += a[i]*a[i];
        const double area_magnitude = sqrt(area_squared);
        double n[3];
        for (int i = 0; i < 3; ++i) n[i] = a[i]/area_magnitude;
        if (x.stage == 5) {
            double u[3]{};
            for (int f = 0; f < 3; ++f)
                for (int i = 0; i < 3; ++i) u[i] += shape[f]*x.velocity[3*f+i];
            for (int i = 0; i < 3; ++i) {
                double tangential = 0, prescribed = 0;
                for (int j = 0; j < 3; ++j) {
                    const double projection = (i == j ? 1.0 : 0.0)-n[i]*n[j];
                    tangential += projection*u[j];
                    prescribed += projection*x.boundary_velocity[3*sample+j];
                    for (int f = 0; f < 3; ++f)
                        out.lhs[(3*row+i)*width+3*f+j] += x.wall_coefficient[sample]*projection*shape[f];
                }
                out.rhs[3*row+i] -= x.wall_coefficient[sample]*(tangential-prescribed);
            }
            continue;
        }
        double mu = 0;
        for (int f = 0; f < 3; ++f) mu += shape[f]*x.viscosity[f];
        for (int i = 0; i < 3; ++i) {
            const int r = 3*row+i;
            const double advected = x.stage == 3 ? x.boundary_velocity[3*sample+i] : x.velocity[r];
            out.rhs[r] -= x.stored_flux[sample]*advected;
            if (x.stage == 4) out.lhs[r*width+r] += x.stored_flux[sample];
            for (int node = 0; node < 4; ++node) {
                // -mu (grad u + grad u^T) A, projected tangentially at the outlet.
                for (int component = 0; component < 3; ++component) {
                    double entry = 0;
                    for (int k = 0; k < 3; ++k) {
                        const double projection = x.stage == 3 ? (i == k ? 1.0 : 0.0)
                                                              : (i == k ? 1.0 : 0.0)-n[i]*n[k];
                        for (int j = 0; j < 3; ++j)
                            entry += -mu*projection*((k == component ? grad[3*node+j] : 0.0)
                                                 +(j == component ? grad[3*node+k] : 0.0))*a[j];
                    }
                    const int c = 3*node+component;
                    out.rhs[r] -= entry*x.velocity[c];
                    out.lhs[r*width+c] += entry*(x.stage == 3 ? x.bc_multiplier[node] : 1.0);
                }
            }
        }
    }
}
} // namespace mars::segregated
#undef MARS_BOUNDARY_HD
