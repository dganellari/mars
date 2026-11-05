#pragma once
#include <HYPRE_parcsr_ls.h>

// Explicit lifetime hooks isolate the optional adapter to the BP3 coarse solver.
extern "C" HYPRE_Int mars_bp3_amg_setup(HYPRE_Solver, HYPRE_ParCSRMatrix,
                                        HYPRE_ParVector, HYPRE_ParVector, int);
extern "C" HYPRE_Int mars_bp3_amg_destroy(HYPRE_Solver);
