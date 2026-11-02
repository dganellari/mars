#include <HYPRE.h>
#include <HYPRE_IJ_mv.h>
#include <HYPRE_parcsr_ls.h>
#include <HYPRE_utilities.h>
#include <mpi.h>

#include <cmath>
#include <cstdio>
#include <vector>

// A small SPD Dirichlet stencil exercises real Hypre cycling with both hook builds.
int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    HYPRE_Initialize();
    constexpr int side = 8, n = side*side*side;
    HYPRE_IJMatrix a;
    HYPRE_IJMatrixCreate(MPI_COMM_WORLD, 0, n-1, 0, n-1, &a);
    HYPRE_IJMatrixSetObjectType(a, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(a);
    for (HYPRE_BigInt row = 0; row < n; ++row) {
        HYPRE_BigInt cols[7] = {row};
        HYPRE_Complex vals[7] = {6};
        HYPRE_Int count = 1;
        const int xyz[3] = {int(row)/64, (int(row)/8)%8, int(row)%8};
        const int strides[3] = {64,8,1};
        for (int d = 0; d < 3; ++d) for (int s : {-1,1}) {
            if (xyz[d]+s < 0 || xyz[d]+s >= side) continue;
            cols[count] = row+s*strides[d]; vals[count++] = -1;
        }
        HYPRE_IJMatrixSetValues(a, 1, &count, &row, cols, vals);
    }
    HYPRE_IJMatrixAssemble(a);
    HYPRE_ParCSRMatrix matrix;
    HYPRE_IJMatrixGetObject(a, reinterpret_cast<void**>(&matrix));
    std::vector<HYPRE_BigInt> indices(n);
    std::vector<HYPRE_Complex> ones(n,1), zeros(n,0), solution(n);
    for (int i = 0; i < n; ++i) indices[i] = i;
    HYPRE_IJVector b, x;
    for (auto vector : {&b,&x}) {
        HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, n-1, vector);
        HYPRE_IJVectorSetObjectType(*vector, HYPRE_PARCSR);
        HYPRE_IJVectorInitialize(*vector);
        HYPRE_IJVectorSetValues(*vector, n, indices.data(), vector == &b ? ones.data() : zeros.data());
        HYPRE_IJVectorAssemble(*vector);
    }
    HYPRE_ParVector rhs, sol;
    HYPRE_IJVectorGetObject(b, reinterpret_cast<void**>(&rhs));
    HYPRE_IJVectorGetObject(x, reinterpret_cast<void**>(&sol));
    HYPRE_Solver amg;
    HYPRE_BoomerAMGCreate(&amg);
    HYPRE_BoomerAMGSetCoarsenType(amg,8);
    HYPRE_BoomerAMGSetInterpType(amg,6);
    HYPRE_BoomerAMGSetPMaxElmts(amg,4);
    HYPRE_BoomerAMGSetRelaxType(amg,18);
    HYPRE_BoomerAMGSetRelaxOrder(amg,0);
    HYPRE_BoomerAMGSetCycleRelaxType(amg,18,3);
    HYPRE_BoomerAMGSetCycleNumSweeps(amg,4,3);
    HYPRE_BoomerAMGSetMaxCoarseSize(amg,9);
    HYPRE_BoomerAMGSetTol(amg,0);
    HYPRE_BoomerAMGSetMaxIter(amg,3);
    HYPRE_BoomerAMGSetup(amg,matrix,rhs,sol);
    for (int i = 0; i < 3; ++i) {
        HYPRE_ParVectorSetConstantValues(sol,0);
        HYPRE_BoomerAMGSolve(amg,matrix,rhs,sol);
    }
    HYPRE_IJVectorGetValues(x,n,indices.data(),solution.data());
    int status = HYPRE_GetError();
    for (double v : solution) { if (!std::isfinite(v)) status = 1; std::printf("%a\n",v); }
    HYPRE_BoomerAMGDestroy(amg);
    HYPRE_IJVectorDestroy(x); HYPRE_IJVectorDestroy(b); HYPRE_IJMatrixDestroy(a);
    HYPRE_Finalize();
    MPI_Finalize();
    return status != 0;
}
