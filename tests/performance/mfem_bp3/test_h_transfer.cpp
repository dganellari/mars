#include "h_transfer.hpp"
#include <cmath>
#include <iostream>
#include <limits>

using namespace mfem;

void check(bool ok)
{
    int local = ok ? 1 : 0, all = 0;
    MPI_Allreduce(&local, &all, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (!all) MPI_Abort(MPI_COMM_WORLD, 1);
}

double relative_difference(const Vector& a, const Vector& b)
{
    double local[2] = {0., 1.}, global[2];
    const double* av = a.HostRead();
    const double* bv = b.HostRead();
    for (int i = 0; i < a.Size(); ++i) {
        local[0] = (!std::isfinite(av[i]) || !std::isfinite(bv[i]))
                       ? std::numeric_limits<double>::infinity()
                       : std::max(local[0], std::abs(av[i] - bv[i]));
        local[1] = std::max(local[1], std::abs(bv[i]));
    }
    MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return std::isfinite(global[0]) && std::isfinite(global[1])
               ? global[0] / global[1] : std::numeric_limits<double>::infinity();
}

int run(const char* config)
{
    Hypre::Init();
    HYPRE_SetGpuAwareMPI(1);
    HYPRE_SetSpMVUseVendor(0);
    HYPRE_SetSpGemmUseVendor(0);
    int rank = Mpi::WorldRank(), size = Mpi::WorldSize(), local_rank = 0;
    MPI_Comm node;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &node);
    MPI_Comm_rank(node, &local_rank);
    MPI_Comm_free(&node);
    int gpus = 1;
#ifdef MFEM_USE_CUDA
    gpus = std::max(1, CuGetDeviceCount());
#endif
    Device device(config, local_rank % gpus);
    Device::SetGPUAwareMPI(true);
    Vector invalid(1), finite(1);
    finite = 0.;
    invalid = std::numeric_limits<double>::quiet_NaN();
    check(std::isinf(relative_difference(finite, invalid)));
    invalid = std::numeric_limits<double>::infinity();
    check(std::isinf(relative_difference(invalid, finite)));
    int dims[3] = {0, 0, 0};
    MPI_Dims_create(size, 3, dims);
    Mesh serial = Mesh::MakeCartesian3D(2, 2, 2, Element::HEXAHEDRON, 1., 1., 1., false);
    std::unique_ptr<int[]> partition(serial.CartesianPartitioning(dims));
    auto* mesh = new ParMesh(MPI_COMM_WORLD, serial, partition.get());
    H1_FECollection collection(1, 3);
    auto* fes = new ParFiniteElementSpace(mesh, &collection);
    ParFiniteElementSpaceHierarchy hierarchy(mesh, fes, true, true);
    double worst = 0.;
    for (int level = 0; level < 2; ++level) {
        bp3_transfer::add_assembled_level(hierarchy);
        auto& coarse = hierarchy.GetFESpaceAtLevel(level);
        auto& fine = hierarchy.GetFESpaceAtLevel(level + 1);
        TrueTransferOperator reference(coarse, fine);
        auto* assembled = hierarchy.GetProlongationAtLevel(level);
        check(dynamic_cast<HypreParMatrix*>(assembled) != nullptr);
        Vector x(coarse.GetTrueVSize()), y(fine.GetTrueVSize());
        Vector px(y.Size()), ref_px(y.Size()), pty(x.Size()), ref_pty(x.Size());
        x.UseDevice(true); y.UseDevice(true); px.UseDevice(true); ref_px.UseDevice(true);
        pty.UseDevice(true); ref_pty.UseDevice(true);
        double* hx = x.HostWrite();
        double* hy = y.HostWrite();
        for (int i = 0; i < x.Size(); ++i) hx[i] = std::sin(.17 * (i + 1 + 31 * rank));
        for (int i = 0; i < y.Size(); ++i) hy[i] = std::cos(.23 * (i + 1 + 19 * rank));
        assembled->Mult(x, px); reference.Mult(x, ref_px);
        assembled->MultTranspose(y, pty); reference.MultTranspose(y, ref_pty);
        double error = std::max(relative_difference(px, ref_px), relative_difference(pty, ref_pty));
        double sums[2] = {px * y, x * pty}, global[2];
        MPI_Allreduce(sums, global, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        check(std::isfinite(global[0]) && std::isfinite(global[1]));
        error = std::max(error, std::abs(global[0] - global[1]) / (1. + std::abs(global[0])));
        x = 1.; assembled->Mult(x, px); ref_px = 1.;
        error = std::max(error, relative_difference(px, ref_px));
        worst = std::max(worst, error);
        check(std::isfinite(error) && error <= 1e-12);
    }
    if (!rank) std::cout << "PASS assembled/native transfer, transpose, adjoint and constants; error " << worst << '\n';
    return 0;
}
int main(int argc, char** argv)
{
    Mpi::Init(argc, argv);
    int status = run(argc > 1 ? argv[1] : "cpu");
    Hypre::Finalize();
    Mpi::Finalize();
    return status;
}
