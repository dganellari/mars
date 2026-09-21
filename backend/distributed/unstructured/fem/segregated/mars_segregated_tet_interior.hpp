// GPT/Codex, 2026-09-21. Interior algebra for the pinned OpenAccel contract.
// Geometry, reconstructed gradients and influence coefficients are inputs.
#pragma once

#if defined(__CUDACC__) || defined(__HIPCC__)
#define MARS_SEGREGATED_HD __host__ __device__
#else
#define MARS_SEGREGATED_HD
#endif

namespace mars::segregated {

struct TetInteriorInput {
    int stage; // 0: pressure, 1: momentum
    int edges[12];
    double coordinates[12], velocity[12], density[4];
    double velocity_shape[24], coordinate_shape[24], shape_gradient[72], area[18];
    double viscosity[4], velocity_blend[12], velocity_gradient[36], stored_flux[6];
    double pressure[4], pressure_gradient[12], influence_lhs[12], influence_rhs[12];
    double density_blend[4], density_gradient[12];
};

struct TetInteriorOutput {
    double lhs[144], rhs[12], flux[6];
};

// Incompressible, fixed mesh/frame, no body force or NSO. No boundary/time terms.
MARS_SEGREGATED_HD inline void tet_interior(const TetInteriorInput& in, TetInteriorOutput& out)
{
    for (int i = 0; i < 144; ++i) out.lhs[i] = 0;
    for (int i = 0; i < 12; ++i) out.rhs[i] = 0;
    for (int s = 0; s < 6; ++s) {
        const int left = in.edges[2*s], right = in.edges[2*s+1];
        const double* area = in.area + 3*s;
        const double* shape = in.velocity_shape + 4*s;
        const double* coord_shape = in.coordinate_shape + 4*s;
        const double* grad = in.shape_gradient + 12*s;
        double point[3] = {};
        for (int n = 0; n < 4; ++n)
            for (int j = 0; j < 3; ++j)
                point[j] += coord_shape[n]*in.coordinates[3*n+j];

        if (in.stage == 0) {
            double u[3] = {}, dp[3] = {}, dl[3] = {}, dr[3] = {};
            for (int n = 0; n < 4; ++n)
                for (int j = 0; j < 3; ++j) {
                    u[j] += shape[n]*in.velocity[3*n+j];
                    dp[j] += grad[3*n+j]*in.pressure[n];
                    dl[j] += shape[n]*in.influence_lhs[3*n+j];
                    dr[j] += shape[n]*in.influence_rhs[3*n+j];
                }
            double volume_flux = 0;
            for (int j = 0; j < 3; ++j) {
                const double reconstructed = 0.5*(in.pressure_gradient[3*left+j]
                                                   + in.pressure_gradient[3*right+j]);
                volume_flux += u[j]*area[j];
                volume_flux -= dr[j]*(dp[j]-reconstructed)*area[j];
            }
            const int up = volume_flux > 0 ? left : right;
            double correction = 0;
            for (int j = 0; j < 3; ++j)
                correction += in.density_blend[up]*(point[j]-in.coordinates[3*up+j])
                              *in.density_gradient[3*up+j];
            const double rho = in.density[up]+correction;
            out.flux[s] = rho*volume_flux;
            for (int n = 0; n < 4; ++n) {
                double value = 0;
                for (int j = 0; j < 3; ++j) value -= rho*dl[j]*grad[3*n+j]*area[j];
                out.lhs[4*left+n] += value;
                out.lhs[4*right+n] -= value;
            }
            out.rhs[left] -= out.flux[s];
            out.rhs[right] += out.flux[s];
        } else {
            const double mass_flux = in.stored_flux[s];
            out.flux[s] = mass_flux;
            const int up = mass_flux > 0 ? left : right;
            const double absolute_flux = mass_flux >= 0 ? mass_flux : -mass_flux;
            double mu = 0;
            for (int n = 0; n < 4; ++n) mu += shape[n]*in.viscosity[n];
            for (int i = 0; i < 3; ++i) {
                const int l = 3*left+i, r = 3*right+i;
                double correction = 0;
                for (int j = 0; j < 3; ++j)
                    correction += in.velocity_blend[3*up+i]*(point[j]-in.coordinates[3*up+j])
                                  *in.velocity_gradient[9*up+3*i+j];
                const double flux = mass_flux*(in.velocity[3*up+i]+correction);
                out.rhs[l] -= flux;
                out.rhs[r] += flux;
                const double positive = 0.5*(mass_flux+absolute_flux);
                const double negative = 0.5*(mass_flux-absolute_flux);
                out.lhs[12*l+3*left+i] += positive;
                out.lhs[12*r+3*left+i] -= positive;
                out.lhs[12*l+3*right+i] += negative;
                out.lhs[12*r+3*right+i] -= negative;
            }
            // Both grad(u) and its transpose belong to the implicit stress block.
            for (int n = 0; n < 4; ++n)
                for (int i = 0; i < 3; ++i) {
                    const int l = 3*left+i, r = 3*right+i;
                    double diagonal = 0;
                    for (int j = 0; j < 3; ++j) {
                        diagonal -= mu*grad[3*n+j]*area[j];
                        const double value = -mu*grad[3*n+i]*area[j];
                        out.lhs[12*l+3*n+j] += value;
                        out.lhs[12*r+3*n+j] -= value;
                        out.rhs[l] -= value*in.velocity[3*n+j];
                        out.rhs[r] += value*in.velocity[3*n+j];
                    }
                    out.lhs[12*l+3*n+i] += diagonal;
                    out.lhs[12*r+3*n+i] -= diagonal;
                    out.rhs[l] -= diagonal*in.velocity[3*n+i];
                    out.rhs[r] += diagonal*in.velocity[3*n+i];
                }
        }
    }
}
} // namespace mars::segregated
#undef MARS_SEGREGATED_HD
