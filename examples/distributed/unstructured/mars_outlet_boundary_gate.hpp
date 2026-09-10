#pragma once

// GPT/Codex, 2026-09-09. Synthetic gates; no case geometry or field data.
#include "../../../backend/distributed/unstructured/fem/mars_outlet_flux.hpp"
#include "../../../backend/distributed/unstructured/fem/mars_outlet_iteration.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#if defined(__CUDACC__) || defined(__HIPCC__)
#define MARS_OUTLET_GATE_HD __host__ __device__
#else
#define MARS_OUTLET_GATE_HD
#endif

namespace mars::outlet_gate {

struct BoundaryInput {
    double coords[4][3]{};
    int face[3]{1, 2, 3};
    int opposite = 0;
    double area[3]{};
    double velocity[3][3]{};
    double pressure = 0;
    double trace[3]{};
    double reconstructed[4][3]{};
    double coefficient[3]{2, 2, 2};
};

struct BoundaryOutput {
    double samples[3]{};
    double derivative = 0;
    double determinant = 0;
    double gradients[4][3]{};
};

MARS_OUTLET_GATE_HD inline BoundaryOutput evaluate(const BoundaryInput& input)
{
    BoundaryOutput output;
    outlet_tet_gradient(input.coords, output.determinant, output.gradients);
    outlet_facet_sample_flux(input.coords, input.face, input.opposite, input.area,
                             input.velocity, input.pressure, input.trace,
                             input.reconstructed, input.coefficient, output.samples,
                             &output.derivative);
    return output;
}

struct ExpectedBoundary {
    std::string name;
    double samples[3]{};
    double derivative = 0;
    double analytic_gradient[3]{};
    bool affine = false;
};

struct Checks {
    int count = 0;
    int failures = 0;

    void require(bool condition, const std::string& name)
    {
        ++count;
        if (!condition) {
            ++failures;
            std::fprintf(stderr, "FAIL %s\n", name.c_str());
        }
    }

    void near(double got, double expected, double tolerance, const std::string& name)
    {
        ++count;
        if (!std::isfinite(got) || std::abs(got - expected) > tolerance) {
            ++failures;
            std::fprintf(stderr, "FAIL %s: got %.17g expected %.17g (tol %.3g)\n",
                         name.c_str(), got, expected, tolerance);
        }
    }
};

inline double dot(const double a[3], const double b[3])
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline void cross(const double a[3], const double b[3], double result[3])
{
    result[0] = a[1]*b[2] - a[2]*b[1];
    result[1] = a[2]*b[0] - a[0]*b[2];
    result[2] = a[0]*b[1] - a[1]*b[0];
}

inline void set_face_geometry(BoundaryInput& input, int opposite)
{
    input.opposite = opposite;
    int face_index = 0;
    for (int node = 0; node < 4; ++node)
        if (node != opposite) input.face[face_index++] = node;
    double edge1[3], edge2[3], inward[3];
    for (int d = 0; d < 3; ++d) {
        edge1[d] = input.coords[input.face[1]][d] - input.coords[input.face[0]][d];
        edge2[d] = input.coords[input.face[2]][d] - input.coords[input.face[0]][d];
        inward[d] = input.coords[opposite][d] - input.coords[input.face[0]][d];
    }
    cross(edge1, edge2, input.area);
    const double factor = dot(input.area, inward) > 0 ? -0.5 : 0.5;
    for (int d = 0; d < 3; ++d) input.area[d] *= factor;
}

inline BoundaryInput unit_tet()
{
    BoundaryInput input;
    for (int d = 0; d < 3; ++d) input.coords[d+1][d] = 1;
    set_face_geometry(input, 0);
    return input;
}

inline double geometric_derivative(const BoundaryInput& input)
{
    // N_opposite falls from 1 to 0 over the normal altitude; no matrix inverse.
    double inward[3];
    for (int d = 0; d < 3; ++d)
        inward[d] = input.coords[input.opposite][d] - input.coords[input.face[0]][d];
    const double area_squared = dot(input.area, input.area);
    const double coefficient = (input.coefficient[0] + input.coefficient[1] +
                                input.coefficient[2])/3;
    return -coefficient*area_squared/(3*dot(input.area, inward));
}

inline void boundary_cases(std::vector<BoundaryInput>& inputs,
                           std::vector<ExpectedBoundary>& expected)
{
    const double linear_gradient[3] = {2, -3, 5};
    const double transform[4][3][3] = {
        {{1,0,0}, {0,1,0}, {0,0,1}},
        {{1,0,0}, {0,1,0}, {0,0,1}},
        {{2,0,0}, {0,2,0}, {0,0,2}},
        {{1,.3,.2}, {.1,1.2,.4}, {.2,-.1,.8}}
    };
    for (int geometry = 0; geometry < 4; ++geometry) {
        for (int opposite = 0; opposite < 4; ++opposite) {
            BoundaryInput input;
            for (int node = 0; node < 4; ++node)
                for (int d = 0; d < 3; ++d)
                    input.coords[node][d] = (geometry == 1 ? (d == 0 ? 3 : d == 1 ? -4 : 2) : 0)
                        + (node == 0 ? 0 : transform[geometry][d][node-1]);
            set_face_geometry(input, opposite);
            ExpectedBoundary affine;
            affine.name = "affine geometry " + std::to_string(geometry) +
                          " opposite " + std::to_string(opposite);
            affine.affine = true;
            for (int d = 0; d < 3; ++d) affine.analytic_gradient[d] = linear_gradient[d];
            input.pressure = 7 + dot(linear_gradient, input.coords[opposite]);
            for (int node = 0; node < 4; ++node)
                for (int d = 0; d < 3; ++d) input.reconstructed[node][d] = linear_gradient[d];
            for (int r = 0; r < 3; ++r) {
                input.trace[r] = 7 + dot(linear_gradient, input.coords[input.face[r]]);
                input.velocity[r][r] = r+1; // Unequal samples catch triangle-mean scattering.
                affine.samples[r] = (r+1)*input.area[r]/3;
            }
            affine.derivative = geometric_derivative(input);
            inputs.push_back(input);
            expected.push_back(affine);

            BoundaryInput constant = input;
            constant.pressure = -6;
            ExpectedBoundary zero;
            zero.name = "constant geometry " + std::to_string(geometry) +
                        " opposite " + std::to_string(opposite);
            zero.derivative = affine.derivative;
            for (int node = 0; node < 4; ++node)
                for (int d = 0; d < 3; ++d) constant.reconstructed[node][d] = 0;
            for (int r = 0; r < 3; ++r) {
                constant.trace[r] = -6;
                for (int d = 0; d < 3; ++d) constant.velocity[r][d] = 0;
            }
            inputs.push_back(constant);
            expected.push_back(zero);
        }
    }

    BoundaryInput input = unit_tet();
    input.pressure = 1;
    input.coefficient[0] = 2;
    input.coefficient[1] = 4;
    input.coefficient[2] = 6;
    ExpectedBoundary face_only;
    face_only.name = "face-only diffusivity gives total flux 6";
    face_only.derivative = 2;
    for (double& sample : face_only.samples) sample = 2;
    inputs.push_back(input);
    expected.push_back(face_only);

    input = unit_tet();
    for (int d = 0; d < 3; ++d) {
        input.reconstructed[0][d] = d+3;
        input.reconstructed[d+1][d] = 3*(d+1);
    }
    ExpectedBoundary reconstructed;
    reconstructed.name = "half face-mean/opposite reconstructed gradient";
    reconstructed.derivative = 1;
    for (double& sample : reconstructed.samples) sample = 3;
    inputs.push_back(input);
    expected.push_back(reconstructed);
}

inline void check_boundary_outputs(const std::vector<BoundaryInput>& inputs,
                                   const std::vector<ExpectedBoundary>& expected,
                                   const std::vector<BoundaryOutput>& outputs, Checks& checks)
{
    checks.require(inputs.size() == outputs.size(), "boundary result count");
    if (inputs.size() != outputs.size()) return;
    for (std::size_t index = 0; index < inputs.size(); ++index) {
        const auto& input = inputs[index];
        const auto& output = outputs[index];
        const auto& want = expected[index];
        for (int r = 0; r < 3; ++r)
            checks.near(output.samples[r], want.samples[r], 1e-11, want.name + " sample " + std::to_string(r));
        checks.near(output.derivative, want.derivative, 1e-12, want.name + " derivative");
        checks.require(std::abs(output.determinant) > .1, want.name + " nondegenerate");
        for (int d = 0; d < 3; ++d) {
            double gradient_sum = 0;
            for (int node = 0; node < 4; ++node) gradient_sum += output.gradients[node][d];
            checks.near(gradient_sum, 0, 1e-12, want.name + " partition of unity");
            if (want.affine) {
                double gradient = input.pressure*output.gradients[input.opposite][d];
                for (int r = 0; r < 3; ++r)
                    gradient += input.trace[r]*output.gradients[input.face[r]][d];
                checks.near(gradient, want.analytic_gradient[d], 1e-11, want.name + " analytic gradient");
            }
        }
    }
}

using Vector = std::array<double, 4>;
using Velocity = std::array<double, 12>;
using Matrix = std::array<Vector, 4>;
using Divergence = std::array<Velocity, 4>;

inline double norm(const Vector& x)
{
    double result = 0;
    for (double value : x) result += value*value;
    return std::sqrt(result);
}

inline Vector multiply(const Matrix& matrix, const Vector& x)
{
    Vector result{};
    for (int row = 0; row < 4; ++row)
        for (int col = 0; col < 4; ++col) result[row] += matrix[row][col]*x[col];
    return result;
}

inline Vector solve(Matrix matrix, Vector rhs, Checks& checks)
{
    for (int col = 0; col < 4; ++col) {
        int pivot = col;
        for (int row = col+1; row < 4; ++row)
            if (std::abs(matrix[row][col]) > std::abs(matrix[pivot][col])) pivot = row;
        checks.require(std::abs(matrix[pivot][col]) > 1e-12, "small operator nonsingular");
        if (std::abs(matrix[pivot][col]) <= 1e-12) return {};
        std::swap(matrix[col], matrix[pivot]);
        std::swap(rhs[col], rhs[pivot]);
        const double diagonal = matrix[col][col];
        for (int j = col; j < 4; ++j) matrix[col][j] /= diagonal;
        rhs[col] /= diagonal;
        for (int row = 0; row < 4; ++row) {
            if (row == col) continue;
            const double factor = matrix[row][col];
            for (int j = col; j < 4; ++j) matrix[row][j] -= factor*matrix[col][j];
            rhs[row] -= factor*rhs[col];
        }
    }
    return rhs;
}

inline constexpr int edges[6][2] = {{0,1}, {1,2}, {0,2}, {0,3}, {1,3}, {2,3}};
inline constexpr double shape_gradient[4][3] = {{-1,-1,-1}, {1,0,0}, {0,1,0}, {0,0,1}};
inline constexpr Matrix stiffness = {{{.5,-1./6,-1./6,-1./6},
                                       {-1./6,1./6,0,0},
                                       {-1./6,0,1./6,0},
                                       {-1./6,0,0,1./6}}};

inline Divergence divergence()
{
    Divergence result{};
    // Median-dual tet SCS area is V*(grad N_R-grad N_L)/4, V=1/6.
    for (const auto& edge : edges) {
        const int left = edge[0], right = edge[1];
        for (int d = 0; d < 3; ++d) {
            const double half_area = (shape_gradient[right][d] - shape_gradient[left][d])/48;
            result[left][3*left+d] += half_area;
            result[left][3*right+d] += half_area;
            result[right][3*left+d] -= half_area;
            result[right][3*right+d] -= half_area;
        }
    }
    for (int row = 1; row < 4; ++row)
        for (int d = 0; d < 3; ++d) result[row][3*row+d] += 1./6;
    return result;
}

inline Velocity gradient(const Divergence& b, const Vector& pressure, const Vector& trace)
{
    Velocity result{};
    for (int velocity = 0; velocity < 12; ++velocity)
        for (int node = 0; node < 4; ++node) result[velocity] -= 24*b[node][velocity]*pressure[node];
    for (int row = 1; row < 4; ++row)
        for (int d = 0; d < 3; ++d) result[3*row+d] += 4*trace[row];
    return result;
}

inline Matrix jacobian(const Divergence& b, const Velocity& mask, double h, double coefficient)
{
    Matrix result{};
    for (int row = 0; row < 4; ++row) {
        for (int col = 0; col < 4; ++col) {
            for (int velocity = 0; velocity < 12; ++velocity)
                result[row][col] += 24*b[row][velocity]*mask[velocity]*b[col][velocity];
            result[row][col] += coefficient*stiffness[row][col]/h;
        }
        if (row > 0) result[row][0] += coefficient/(2*h);
    }
    return result;
}

inline Vector residual(const Velocity& velocity, const Vector& pressure,
                       const Vector& trace, double coefficient, double* boundary_sum = nullptr)
{
    Vector result{};
    double compact_gradient[3]{};
    const double reconstructed[3] = {.2, -.3, .7};
    for (int node = 0; node < 4; ++node)
        for (int d = 0; d < 3; ++d) compact_gradient[d] += pressure[node]*shape_gradient[node][d];
    for (const auto& edge : edges) {
        const int left = edge[0], right = edge[1];
        double flux = 0;
        for (int d = 0; d < 3; ++d) {
            const double area = (shape_gradient[right][d] - shape_gradient[left][d])/24;
            flux += (.5*(velocity[3*left+d] + velocity[3*right+d]) +
                     coefficient*(reconstructed[d] - compact_gradient[d]))*area;
        }
        result[left] += flux;
        result[right] -= flux;
    }
    BoundaryInput input = unit_tet();
    input.pressure = pressure[0];
    for (int r = 0; r < 3; ++r) {
        input.trace[r] = trace[r+1];
        input.coefficient[r] = coefficient;
        for (int d = 0; d < 3; ++d) input.velocity[r][d] = velocity[3*(r+1)+d];
    }
    for (int node = 0; node < 4; ++node)
        for (int d = 0; d < 3; ++d) input.reconstructed[node][d] = reconstructed[d];
    const auto output = evaluate(input);
    double total = -.17; // Fixed exterior source belongs in both sides of the balance.
    result[0] += total;
    for (int r = 0; r < 3; ++r) {
        result[r+1] += output.samples[r];
        total += output.samples[r];
    }
    if (boundary_sum) *boundary_sum = total;
    return result;
}

inline void algebra_checks(Checks& checks)
{
    const auto b = divergence();
    Velocity mask;
    mask.fill(1);
    const Matrix h_expected = {{{33./48,5./48,5./48,5./48},
                                {5./48,77./48,-1./48,-1./48},
                                {5./48,-1./48,77./48,-1./48},
                                {5./48,-1./48,-1./48,77./48}}};
    const Matrix h_computed = jacobian(b, mask, 1, 0);
    for (int row = 0; row < 4; ++row)
        for (int col = 0; col < 4; ++col)
            checks.near(h_computed[row][col], h_expected[row][col], 1e-13, "independent boundary Gram fixture");

    const Vector pressure{2, -1, .5, 3}, trace{0, .7, -.3, 1.1}, direction{.3, -.7, 1.2, .4};
    Velocity velocity;
    for (int i = 0; i < 12; ++i) velocity[i] = (i%5-2)*.13;
    const auto g = gradient(b, pressure, trace);
    double work = 0;
    for (int i = 0; i < 12; ++i) work += velocity[i]*g[i]/24;
    for (int row = 0; row < 4; ++row)
        for (int i = 0; i < 12; ++i) work += pressure[row]*b[row][i]*velocity[i];
    for (int row = 1; row < 4; ++row)
        for (int d = 0; d < 3; ++d) work -= trace[row]*velocity[3*row+d]/6;
    checks.near(work, 0, 1e-13, "pressure work with boundary adjoint");
    const Vector constant{3, 3, 3, 3};
    for (double value : gradient(b, constant, constant))
        checks.near(value, 0, 1e-13, "matching constant pressure and trace gradient");
    checks.require(norm(multiply(h_computed, constant)) > 1, "free outlet temporal term anchors constant mode");

    constexpr double coefficient = .3;
    const auto base = residual(velocity, pressure, trace, coefficient);
    for (int startup = 0; startup < 2; ++startup) {
        const double rho = 2.7, dt = .14;
        const double h = (startup == 0 ? dt : 2*dt/3)/rho;
        for (int prescribed = 0; prescribed < 2; ++prescribed) {
            mask.fill(1);
            if (prescribed) for (int d = 0; d < 3; ++d) mask[d] = 0;
            const Matrix j = jacobian(b, mask, h, coefficient);
            const auto action = multiply(j, direction);
            const auto delta_gradient = gradient(b, direction, {});
            for (double epsilon : {1e-3, 1e-4, 1e-5}) {
                auto perturbed_pressure = pressure;
                auto perturbed_velocity = velocity;
                for (int node = 0; node < 4; ++node) perturbed_pressure[node] += epsilon*direction[node];
                for (int i = 0; i < 12; ++i)
                    perturbed_velocity[i] -= epsilon*h*mask[i]*delta_gradient[i];
                const auto perturbed = residual(perturbed_velocity, perturbed_pressure, trace, coefficient);
                Vector error{};
                for (int node = 0; node < 4; ++node)
                    error[node] = (perturbed[node]-base[node])/epsilon-h*action[node];
                checks.require(norm(error)/std::max(1e-12, h*norm(action)) < 1e-8,
                               "full velocity-pressure JVP, BDF" + std::to_string(startup+1) +
                               " prescribed=" + std::to_string(prescribed));
            }
            Vector rhs;
            for (int node = 0; node < 4; ++node) rhs[node] = -base[node]/h;
            const auto correction = solve(j, rhs, checks);
            auto corrected_pressure = pressure;
            auto corrected_velocity = velocity;
            const auto correction_gradient = gradient(b, correction, {});
            for (int node = 0; node < 4; ++node) corrected_pressure[node] += correction[node];
            for (int i = 0; i < 12; ++i) corrected_velocity[i] -= h*mask[i]*correction_gradient[i];
            checks.require(norm(residual(corrected_velocity, corrected_pressure, trace, coefficient)) < 1e-12,
                           "exact small correction closes all continuity rows");
        }
    }

    double exterior = 0;
    const auto all_rows = residual(velocity, pressure, trace, coefficient, &exterior);
    double sum = 0;
    for (double value : all_rows) sum += value;
    checks.near(sum, exterior, 1e-13, "all rows equal unique exterior flux plus prescribed source");
    checks.require(std::abs(all_rows[0]-exterior) > 1e-3, "interior-only residual cannot certify boundary balance");

    mask.fill(1);
    const auto j = jacobian(b, mask, 1, coefficient);
    Matrix approximate{};
    for (int row = 0; row < 4; ++row) {
        for (int col = 0; col < 4; ++col) approximate[row][col] = 1.3*stiffness[row][col];
        if (row > 0) approximate[row][0] += .15;
    }
    // Symmetric outlet-node subspace: solve the exact 2x2 generalized eigenproblem.
    const double aa = approximate[0][0], ab = 3*approximate[0][1];
    const double ac = approximate[1][0], ad = approximate[1][1];
    const double ja = j[0][0], jb = 3*j[0][1];
    const double jc = j[1][0], jd = j[1][1]+j[1][2]+j[1][3];
    const double quadratic = aa*ad-ab*ac;
    const double linear = -(ja*ad+jd*aa-jb*ac-jc*ab);
    const double constant_term = ja*jd-jb*jc;
    const double lambda = (-linear+std::sqrt(linear*linear-4*quadratic*constant_term))/(2*quadratic);
    checks.near(lambda, 13.0478622748, 1e-9, "compact correction dominant eigenvalue");
    const double y = -(ja-lambda*aa)/(jb-lambda*ab);
    const Vector eigenvector{1, y, y, y};
    const auto initial = multiply(j, eigenvector);
    const auto approximate_step = solve(approximate, initial, checks);
    const auto residual_step = multiply(j, approximate_step);
    for (double omega : {1., .25, .1}) {
        Vector next;
        for (int node = 0; node < 4; ++node) next[node] = initial[node]-omega*residual_step[node];
        checks.near(norm(next)/norm(initial), std::abs(1-omega*lambda), 1e-11,
                    "measured residual contraction at damping " + std::to_string(omega));
    }
    checks.require(std::abs(1-.25*lambda) > 1 && std::abs(1-.1*lambda) < 1,
                   "quarter damping fails, tenth contracts this mode only");

    Vector remaining = multiply(j, {1., -.2, .3, .7});
    const double initial_norm = norm(remaining);
    int iterations = 0;
    while (norm(remaining) > 1e-10*initial_norm && iterations < 100) {
        Vector rhs;
        for (int node = 0; node < 4; ++node) rhs[node] = -remaining[node];
        const auto step = solve(approximate, rhs, checks);
        const auto delta = multiply(j, step);
        double product = 0, square = 0;
        for (int node = 0; node < 4; ++node) {
            product += 24*remaining[node]*delta[node];
            square += 24*delta[node]*delta[node];
        }
        const double omega = outlet_correction_damping(product, square, 1);
        checks.require(omega > 0 && omega <= 1, "production damping accepts residual descent");
        Vector next;
        for (int node = 0; node < 4; ++node) next[node] = remaining[node]+omega*delta[node];
        const bool contracts = outlet_correction_contracts(norm(remaining), norm(next), omega);
        checks.require(contracts, "production Armijo check contracts actual mixed-mode residual");
        if (!contracts) break;
        remaining = next;
        ++iterations;
    }
    checks.require(iterations > 1 && iterations < 100 && norm(remaining) <= 1e-10*initial_norm,
                   "residual-optimal correction converges the compact counterexample");
}

inline void iteration_checks(Checks& checks)
{
    // R=(1,2), delta R=(-3,-1), V=(1,4) gives omega=14/37.
    checks.near(outlet_correction_damping(-3.5, 9.25, 1), 14./37, 1e-15,
                "production volume-weighted optimum");
    checks.near(outlet_correction_damping(-2, 1, .25), .25, 1e-15, "production maximum damping clamp");
    const double infinity = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (const auto& values : std::vector<std::array<double,3>>{
            {0,1,1}, {1,1,1}, {-1,0,1}, {-1,-1,1}, {-1,1,0}, {-1,1,-1}, {-1,1,1.1},
            {infinity,1,1}, {-infinity,1,1}, {-1,infinity,1}, {-1,1,infinity},
            {nan,1,1}, {-1,nan,1}, {-1,1,nan}})
        checks.near(outlet_correction_damping(values[0], values[1], values[2]), 0, 0,
                    "production damping rejects invalid or non-descent direction");
    for (const auto& values : std::vector<std::array<double,3>>{
            {0,0,1}, {1,1,1}, {1,2,.1}, {1,-.1,.1}, {1,.5,0}, {1,.5,-1}, {1,.5,1.1},
            {1,1-1e-5,1}, {infinity,.5,1}, {1,infinity,1}, {nan,.5,1},
            {1,nan,1}, {1,.5,nan}, {1,.5,infinity}})
        checks.require(!outlet_correction_contracts(values[0], values[1], values[2]),
                       "production contraction rejects stagnation, growth, nonfinite or weak decrease");
    checks.require(outlet_correction_contracts(1,.9998,1), "production Armijo accepts sufficient decrease");
}

inline void derivative_checks(Checks& checks)
{
    std::vector<BoundaryInput> inputs;
    std::vector<ExpectedBoundary> expected;
    boundary_cases(inputs, expected);
    for (const auto& input : inputs) {
        const auto base = evaluate(input);
        for (int startup = 0; startup < 2; ++startup) {
            const double rho = 2.7, dt = .14;
            const double dt_effective = startup == 0 ? dt : 2*dt/3;
            checks.near(base.derivative*rho/dt_effective,
                        geometric_derivative(input)*rho/dt_effective, 1e-11,
                        "compact matrix partial includes rho/dt_eff");
        }
        for (double epsilon : {1e-3, 1e-4, 1e-5}) {
            auto perturbed = input;
            perturbed.pressure += epsilon;
            const auto result = evaluate(perturbed);
            for (int sample = 0; sample < 3; ++sample)
                checks.near((result.samples[sample]-base.samples[sample])/epsilon,
                            base.derivative, 1e-8*std::max(1., std::abs(base.derivative)),
                            "frozen-trace opposite pressure finite difference");
            perturbed = input;
            for (double& trace : perturbed.trace) trace += epsilon;
            const auto shifted_trace = evaluate(perturbed);
            for (int sample = 0; sample < 3; ++sample)
                checks.near((shifted_trace.samples[sample]-base.samples[sample])/epsilon,
                            -base.derivative, 1e-8*std::max(1., std::abs(base.derivative)),
                            "common trace pressure finite difference");
        }
    }
}

inline void geometry_ownership_checks(Checks& checks)
{
    // Two independently constructed triangles share one vertex. Their normals
    // are either perpendicular or opposite; scalar area must survive both.
    for (bool opposed : {false, true}) {
        BoundaryInput facets[2];
        facets[0].coords[0][2] = -1;
        facets[0].coords[2][0] = 1;
        facets[0].coords[3][1] = 2;
        if (opposed) {
            facets[1].coords[0][2] = 1;
            facets[1].coords[2][0] = 1;
            facets[1].coords[3][1] = -2;
        } else {
            facets[1].coords[0][0] = -1;
            facets[1].coords[2][1] = 1;
            facets[1].coords[3][2] = 2;
        }
        const int global_nodes[2][3] = {{0,1,2}, {0,3,4}};
        std::array<double, 5> scalar_area{}, serial_residual{};
        std::array<std::array<double, 3>, 5> vector_area{};
        double exterior = 0;
        for (int facet = 0; facet < 2; ++facet) {
            auto& input = facets[facet];
            set_face_geometry(input, 0);
            input.pressure = 3;
            const double area = std::sqrt(dot(input.area, input.area));
            checks.near(area, 1, 1e-13, "triangle-coordinate area");
            for (int r = 0; r < 3; ++r) {
                const int node = global_nodes[facet][r];
                input.trace[r] = 3;
                scalar_area[node] += area/3;
                for (int d = 0; d < 3; ++d) {
                    vector_area[node][d] += input.area[d]/3;
                    input.velocity[r][d] = .1*(node+1)*(d+1);
                }
            }
            const auto output = evaluate(input);
            for (int r = 0; r < 3; ++r) {
                serial_residual[global_nodes[facet][r]] += output.samples[r];
                exterior += output.samples[r];
            }
        }
        checks.near(scalar_area[0], 2./3, 1e-13, "shared vertex physical area survives bent/opposed normals");
        const double vector_weight = std::sqrt(dot(vector_area[0].data(), vector_area[0].data()));
        checks.near(vector_weight, opposed ? 0 : std::sqrt(2.)/3, 1e-13,
                    "norm of summed vector is not physical area");
        double area_sum = 0, moment = 0;
        for (int node = 0; node < 5; ++node) {
            area_sum += scalar_area[node];
            if (node != 0) moment += 3*scalar_area[node];
        }
        checks.near(area_sum, 2, 1e-13, "unique triangle scalar area sum");
        checks.near(moment/area_sum, 2, 1e-13, "coordinate-derived physical pressure mean");
        const double wrong_mean = moment/(4./3+vector_weight);
        checks.require(std::abs(wrong_mean-2) > .2, "vector-norm pressure mean mutation is detected");

        for (int ranks : {1,2,4}) {
            std::array<double, 5> reduced{};
            int empty_ranks = 0;
            for (int rank = 0; rank < ranks; ++rank) {
                std::array<double, 5> local{};
                int owned_facets = 0;
                for (int facet = 0; facet < 2; ++facet) {
                    if (facet%ranks != rank) continue;
                    ++owned_facets;
                    const auto output = evaluate(facets[facet]);
                    for (int r = 0; r < 3; ++r) local[global_nodes[facet][r]] += output.samples[r];
                }
                if (!owned_facets) ++empty_ranks;
                for (int node = 0; node < 5; ++node) reduced[node] += local[node];
            }
            double global_sum = 0;
            for (int node = 0; node < 5; ++node) {
                checks.near(reduced[node], serial_residual[node], 1e-13,
                            "owned-facet scatter partition model " + std::to_string(ranks));
                global_sum += reduced[node];
            }
            checks.near(global_sum, exterior, 1e-13, "partition model global boundary conservation");
            checks.require(ranks != 4 || empty_ranks == 2, "four-partition fixture has empty-facet ranks");
        }
    }
}

inline int run_host_checks()
{
    Checks checks;
    std::vector<BoundaryInput> inputs;
    std::vector<ExpectedBoundary> expected;
    boundary_cases(inputs, expected);
    std::vector<BoundaryOutput> outputs;
    for (const auto& input : inputs) outputs.push_back(evaluate(input));
    check_boundary_outputs(inputs, expected, outputs, checks);
    derivative_checks(checks);
    algebra_checks(checks);
    iteration_checks(checks);
    geometry_ownership_checks(checks);
    std::printf("%s: %d host checks of the production outlet evaluator and synthetic correction algebra\n",
                checks.failures ? "FAIL" : "PASS", checks.count);
    std::printf("The small correction fixture is not a flow, CUDA scatter, or MPI validation.\n");
    return checks.failures ? 1 : 0;
}

} // namespace mars::outlet_gate

#undef MARS_OUTLET_GATE_HD
