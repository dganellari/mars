// The pump solver's hex (HEX8) pieces on one element, against closed forms, and an explicit
// instantiation of the pump stepper for hexes (it must compile).
//
//   mars_pump_hex_gate        prints PASS/FAIL per check, exit code 1 on any failure
//
// The element is a frustum: a 2 x 2 square at z = 0 under a 1 x 1 square at z = 1, shifted by
// (0.2, 0.1), so the hex is not a parallelepiped and its faces are planar trapezoids. Its volume
// is h (A0 + A1 + sqrt(A0 A1)) / 3 = 7/3 (a bounding box gives 2 x 2 x 1 = 4).
//
// Checks:
//   lumped mass        sum of the eight sub-control volumes = 7/3, each one > 0
//   convection         C_ij = int N_i (a . grad N_j): rows sum to 0 (sum_j N_j = 1), and
//                      sum_ij C_ij x_j = a_x V (the field x has a . grad x = a_x)
//   SUPG               S = C(supg) - C: symmetric, rows sum to 0, positive semi-definite on
//                      random vectors
#include "backend/distributed/unstructured/fem/mars_ns_pump_solver.hpp"

#include <cstdio>
#include <random>

// The hex stepper's whole pump path (setup, step, penalised projection, inlet ramp).
template void setupNSStepper<uint64_t, double, HexTag>(NSStepper<uint64_t, double, HexTag>&, double, double,
                                                       CvfemKernelVariant);
template void runNsStep<uint64_t, double, HexTag>(NSStepper<uint64_t, double, HexTag>&, double, double, double);
template void setPenalisedProjection<uint64_t, double, HexTag>(NSStepper<uint64_t, double, HexTag>&, const double*,
                                                               double);
template void rescaleInletVelocityTarget<uint64_t, double, HexTag>(NSStepper<uint64_t, double, HexTag>&, double);

namespace
{

int failures = 0;

void check(bool ok, const char* what, double value, double expected)
{
    std::printf("%s %-40s %.15g (expected %.15g)\n", ok ? "PASS" : "FAIL", what, value, expected);
    if (!ok) ++failures;
}

void cuda_check(cudaError_t error, const char* operation)
{
    if (error == cudaSuccess) return;
    std::fprintf(stderr, "FAIL %s: %s\n", operation, cudaGetErrorString(error));
    MPI_Abort(MPI_COMM_WORLD, 1);
}

template<class T>
T* to_device(const std::vector<T>& h)
{
    T* d = nullptr;
    cuda_check(cudaMalloc(&d, h.size() * sizeof(T)), "cudaMalloc");
    cuda_check(cudaMemcpy(d, h.data(), h.size() * sizeof(T), cudaMemcpyHostToDevice), "copy to device");
    return d;
}

template<class T>
std::vector<T> to_host(const T* d, size_t n)
{
    std::vector<T> h(n);
    cuda_check(cudaMemcpy(h.data(), d, n * sizeof(T), cudaMemcpyDeviceToHost), "copy to host");
    return h;
}

// One element, nodes 0..7, a dense 8 x 8 CSR (one DOF per node).
struct OneHex
{
    std::vector<double> x, y, z;
    uint64_t* conn[8];
    double *dx, *dy, *dz;
    int *rowPtr, *colInd, *nodeToDof;
    uint8_t *ownership, *isBdry;

    OneHex()
    {
        // Exodus order: (-,-,-) (+,-,-) (+,+,-) (-,+,-) at z = 0, then the same at z = 1.
        const double b[4][2] = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
        for (int k = 0; k < 4; ++k)
        {
            x.push_back(b[k][0]);
            y.push_back(b[k][1]);
            z.push_back(0);
        }
        for (int k = 0; k < 4; ++k)
        {
            x.push_back(0.2 + 0.5 * b[k][0]);
            y.push_back(0.1 + 0.5 * b[k][1]);
            z.push_back(1);
        }
        for (int c = 0; c < 8; ++c) conn[c] = to_device(std::vector<uint64_t>{uint64_t(c)});
        dx = to_device(x);
        dy = to_device(y);
        dz = to_device(z);
        std::vector<int> rp(9), ci(64), n2d(8);
        for (int i = 0; i <= 8; ++i) rp[i] = 8 * i;
        for (int k = 0; k < 64; ++k) ci[k] = k % 8;
        for (int i = 0; i < 8; ++i) n2d[i] = i;
        rowPtr    = to_device(rp);
        colInd    = to_device(ci);
        nodeToDof = to_device(n2d);
        ownership = to_device(std::vector<uint8_t>(8, 1));
        isBdry    = to_device(std::vector<uint8_t>(8, 0));
    }
};

std::vector<double> convection(const OneHex& h, const double a[3], int supg)
{
    double* u    = to_device(std::vector<double>(8, a[0]));
    double* v    = to_device(std::vector<double>(8, a[1]));
    double* w    = to_device(std::vector<double>(8, a[2]));
    double* vals = to_device(std::vector<double>(64, 0.0));
    assembleConvectionHexKernel<uint64_t, double><<<1, 32>>>(
        h.conn[0], h.conn[1], h.conn[2], h.conn[3], h.conn[4], h.conn[5], h.conn[6], h.conn[7], h.dx, h.dy, h.dz,
        u, v, w, h.nodeToDof, h.ownership, h.isBdry, h.rowPtr, h.colInd, 8, 1e-6, 1e-3, supg, vals, 1);
    cuda_check(cudaGetLastError(), "assembleConvectionHexKernel");
    cuda_check(cudaDeviceSynchronize(), "assembleConvectionHexKernel");
    std::vector<double> c = to_host(vals, 64);
    cudaFree(u);
    cudaFree(v);
    cudaFree(w);
    cudaFree(vals);
    return c;
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    const double V = 7.0 / 3.0;
    OneHex h;

    // Lumped mass.
    double* mass = to_device(std::vector<double>(8, 0.0));
    computeLumpedMassPerNodeKernel<uint64_t, double, HexTag><<<1, 32>>>(
        h.conn[0], h.conn[1], h.conn[2], h.conn[3], h.conn[4], h.conn[5], h.conn[6], h.conn[7], h.dx, h.dy, h.dz,
        mass, 0, 1);
    cuda_check(cudaGetLastError(), "computeLumpedMassPerNodeKernel");
    std::vector<double> m = to_host(mass, 8);
    double sum = 0, mmin = m[0];
    for (double mi : m)
    {
        sum += mi;
        mmin = std::min(mmin, mi);
    }
    check(std::fabs(sum - V) < 1e-12, "lumped mass: sum = volume", sum, V);
    check(mmin > 0, "lumped mass: every sub-volume > 0", mmin, 0);

    // Convection, plain and with SUPG.
    const double a[3] = {0.7, -0.3, 0.5};
    std::vector<double> c = convection(h, a, 0);
    std::vector<double> cs = convection(h, a, 1);
    double rowMax = 0, total = 0, supgRowMax = 0, asym = 0, scale = 0;
    for (int i = 0; i < 8; ++i)
    {
        double row = 0, srow = 0;
        for (int j = 0; j < 8; ++j)
        {
            row += c[8 * i + j];
            srow += cs[8 * i + j] - c[8 * i + j];
            total += c[8 * i + j] * h.x[j];
            asym = std::max(asym, std::fabs((cs[8 * i + j] - c[8 * i + j]) - (cs[8 * j + i] - c[8 * j + i])));
            scale = std::max(scale, std::fabs(c[8 * i + j]));
        }
        rowMax = std::max(rowMax, std::fabs(row));
        supgRowMax = std::max(supgRowMax, std::fabs(srow));
    }
    check(rowMax < 1e-12 * scale, "convection: rows sum to 0", rowMax, 0);
    check(std::fabs(total - a[0] * V) < 1e-12, "convection: sum C x = a_x V", total, a[0] * V);
    check(supgRowMax < 1e-12 * scale, "SUPG: rows sum to 0", supgRowMax, 0);
    check(asym < 1e-14 * scale, "SUPG: symmetric", asym, 0);
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> uni(-1, 1);
    double qmin = 1;
    for (int t = 0; t < 100; ++t)
    {
        double r[8], q = 0;
        for (int i = 0; i < 8; ++i) r[i] = uni(rng);
        for (int i = 0; i < 8; ++i)
            for (int j = 0; j < 8; ++j) q += r[i] * (cs[8 * i + j] - c[8 * i + j]) * r[j];
        qmin = std::min(qmin, q);
    }
    check(qmin > -1e-14 * scale, "SUPG: positive semi-definite", qmin, 0);
    double diag = 0;
    for (int i = 0; i < 8; ++i) diag += cs[8 * i + i] - c[8 * i + i];
    check(diag > 0, "SUPG: active (trace > 0)", diag, 0);

    std::printf("%s\n", failures ? "FAILED" : "ALL PASSED");
    MPI_Finalize();
    return failures ? 1 : 0;
}
