// Public communication fixture for validating the optional metadata interposer.
#include <HYPRE.h>
#include <HYPRE_IJ_mv.h>
#include <HYPRE_parcsr_ls.h>
#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <dlfcn.h>

// Exercise the hook with an invalid public index without sending it to MPI.
static HYPRE_Int partition_setup(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector) {
    using Lookup = HYPRE_Int(*)(MPI_Comm,HYPRE_BigInt,HYPRE_BigInt,HYPRE_BigInt,HYPRE_Int*);
    auto lookup = reinterpret_cast<Lookup>(dlsym(RTLD_DEFAULT, "hypre_GetAssumedPartitionProcFromRow"));
    int rank, ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    HYPRE_Int owner = -1;
    const HYPRE_BigInt total = 40 * ranks;
    if (!lookup || lookup(MPI_COMM_WORLD, 40 * rank, 0, total, &owner) || owner != rank) return 1;
    if (lookup(MPI_COMM_WORLD, total, 0, total, &owner) || owner != ranks) return 1;
    return 0;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    setenv("SLURM_PROCID", std::to_string(rank).c_str(), 1);
    HYPRE_Initialize();
    constexpr HYPRE_Int local = 40;
    const HYPRE_BigInt begin = rank * local, end = begin + local - 1, total = ranks * local;
    HYPRE_IJMatrix matrix;
    HYPRE_IJMatrixCreate(MPI_COMM_WORLD, begin, end, begin, end, &matrix);
    HYPRE_IJMatrixSetObjectType(matrix, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(matrix);
    for (HYPRE_BigInt row = begin; row <= end; ++row) {
        HYPRE_Int n = 0;
        HYPRE_BigInt columns[3]; HYPRE_Complex values[3];
        columns[n] = row; values[n++] = 2;
        if (row) { columns[n] = row - 1; values[n++] = -1; }
        if (row + 1 < total) { columns[n] = row + 1; values[n++] = -1; }
        HYPRE_IJMatrixSetValues(matrix, 1, &n, &row, columns, values);
    }
    HYPRE_IJMatrixAssemble(matrix);
    HYPRE_IJVector rhs, solution;
    HYPRE_IJVectorCreate(MPI_COMM_WORLD, begin, end, &rhs);
    HYPRE_IJVectorCreate(MPI_COMM_WORLD, begin, end, &solution);
    for (auto vector : {rhs, solution}) {
        HYPRE_IJVectorSetObjectType(vector, HYPRE_PARCSR);
        HYPRE_IJVectorInitialize(vector);
        HYPRE_IJVectorAssemble(vector);
    }
    HYPRE_ParCSRMatrix a; HYPRE_ParVector b, x;
    HYPRE_IJMatrixGetObject(matrix, reinterpret_cast<void**>(&a));
    HYPRE_IJVectorGetObject(rhs, reinterpret_cast<void**>(&b));
    HYPRE_IJVectorGetObject(solution, reinterpret_cast<void**>(&x));
    HYPRE_Solver solver, amg;
    HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, &solver);
    HYPRE_BoomerAMGCreate(&amg);
    HYPRE_BoomerAMGSetCoarsenType(amg, 10);
    HYPRE_BoomerAMGSetInterpType(amg, 6);
    HYPRE_BoomerAMGSetMaxCoarseSize(amg, 9);
    const bool partition_test = argc == 2 && std::string(argv[1]) == "--partition-check";
    HYPRE_ParCSRGMRESSetPrecond(solver, HYPRE_BoomerAMGSolve,
                              partition_test ? partition_setup : HYPRE_BoomerAMGSetup, amg);
    int code = HYPRE_ParCSRGMRESSetup(solver, a, b, x);
    HYPRE_ParCSRGMRESDestroy(solver);
    HYPRE_BoomerAMGDestroy(amg);
    HYPRE_IJVectorDestroy(rhs); HYPRE_IJVectorDestroy(solution); HYPRE_IJMatrixDestroy(matrix);
    HYPRE_Finalize();
    int failed = code != 0, any;
    MPI_Allreduce(&failed, &any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (!rank) std::puts(any ? "FAIL: public setup trace fixture" : "PASS: public setup trace fixture");
    MPI_Finalize();
    return any;
}
