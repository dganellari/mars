#include <HYPRE.h>
#include <HYPRE_krylov.h>
#include <HYPRE_parcsr_ls.h>
#include <_hypre_parcsr_ls.h>

#include <cmath>
#include <dlfcn.h>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

static void checked(HYPRE_Int code) {
    if (code != 0) throw std::runtime_error("Hypre call failed");
}

static std::string quoted(const char* text) {
    std::string result = "\"";
    for (const unsigned char c : std::string(text)) {
        if (c < 32) throw std::runtime_error("Invalid metadata");
        if (c == '"' || c == '\\') result += '\\';
        result += c;
    }
    return result + '"';
}

int main(int argc, char** argv) {
    try {
        checked(hypre_MPI_Init(&argc, &argv));
        HYPRE_Int ranks = 0;
        checked(hypre_MPI_Comm_size(hypre_MPI_COMM_WORLD, &ranks));
        if (ranks != 1) throw std::runtime_error("One rank required");
        HYPRE_Int major, minor, patch, release;
        checked(HYPRE_VersionNumber(&major, &minor, &patch, &release));
        if (release != HYPRE_RELEASE_NUMBER) throw std::runtime_error("Header version mismatch");
        checked(HYPRE_Init());
        HYPRE_Solver amg = nullptr, gmres = nullptr, flex = nullptr;
        checked(HYPRE_BoomerAMGCreate(&amg));
        checked(HYPRE_ParCSRGMRESCreate(hypre_MPI_COMM_WORLD, &gmres));
        checked(HYPRE_ParCSRFlexGMRESCreate(hypre_MPI_COMM_WORLD, &flex));
        auto* data = reinterpret_cast<hypre_ParAMGData*>(amg);
        std::map<std::string, double> values;

        // Verify available getters before using matching installed headers for the rest.
        auto integer = [&](const char* name, HYPRE_Int (*get)(void*, HYPRE_Int*), HYPRE_Int expected) {
            HYPRE_Int value = 0;
            checked(get(amg, &value));
            if (value != expected) throw std::runtime_error("Header layout mismatch");
            values[name] = static_cast<double>(value);
        };
        auto real = [&](const char* name, HYPRE_Int (*get)(void*, HYPRE_Real*), HYPRE_Real expected) {
            HYPRE_Real value = 0;
            checked(get(amg, &value));
            if (value != expected) throw std::runtime_error("Header layout mismatch");
            values[name] = static_cast<double>(value);
        };
#define READ_INT(label, field) integer(label, hypre_BoomerAMGGet##field, hypre_ParAMGData##field(data))
#define READ_REAL(label, field) real(label, hypre_BoomerAMGGet##field, hypre_ParAMGData##field(data))
        READ_INT("amg_coarsen_type", CoarsenType);
        READ_INT("amg_interp_type", InterpType);
        READ_INT("amg_relax_order", RelaxOrder);
        READ_INT("amg_p_max_elmts", PMaxElmts);
        READ_INT("amg_max_levels", MaxLevels);
        READ_INT("amg_min_coarse_size", MinCoarseSize);
        READ_INT("amg_max_coarse_size", MaxCoarseSize);
        READ_INT("amg_num_functions", NumFunctions);
        READ_INT("amg_coarsen_cut_factor", CoarsenCutFactor);
        READ_INT("amg_cycle_type", CycleType);
        READ_INT("amg_fcycle", FCycle);
        READ_REAL("amg_strong_threshold", StrongThreshold);
        READ_REAL("amg_trunc_factor", TruncFactor);
        READ_REAL("amg_jacobi_trunc_threshold", JacobiTruncThreshold);
        READ_REAL("amg_max_row_sum", MaxRowSum);
#undef READ_INT
#undef READ_REAL
        for (int cycle = 1; cycle <= 3; ++cycle) {
            const std::string suffix = cycle == 1 ? "down" : cycle == 2 ? "up" : "coarse";
            HYPRE_Int relax = 0, sweeps = 0;
            checked(hypre_BoomerAMGGetCycleRelaxType(amg, &relax, cycle));
            checked(hypre_BoomerAMGGetCycleNumSweeps(amg, &sweeps, cycle));
            if (relax != hypre_ParAMGDataGridRelaxType(data)[cycle] ||
                sweeps != hypre_ParAMGDataNumGridSweeps(data)[cycle])
                throw std::runtime_error("Header layout mismatch");
            values["amg_relax_" + suffix] = relax;
            values["amg_sweeps_" + suffix] = sweeps;
        }
        values["amg_agg_num_levels"] = hypre_ParAMGDataAggNumLevels(data);
        values["amg_agg_interp_type"] = hypre_ParAMGDataAggInterpType(data);
        values["amg_agg_trunc_factor"] = hypre_ParAMGDataAggTruncFactor(data);
        values["amg_num_paths"] = hypre_ParAMGDataNumPaths(data);
        values["amg_keep_transpose"] = hypre_ParAMGDataKeepTranspose(data);
        values["amg_nodal"] = hypre_ParAMGDataNodal(data);
        values["amg_nodal_diag"] = hypre_ParAMGDataNodalDiag(data);

        auto krylov = [&](const char* prefix, HYPRE_Solver solver,
                          HYPRE_Int (*get_kdim)(HYPRE_Solver, HYPRE_Int*),
                          HYPRE_Int (*get_min)(HYPRE_Solver, HYPRE_Int*)) {
            HYPRE_Int kdim = 0, min_iter = 0;
            checked(get_kdim(solver, &kdim));
            checked(get_min(solver, &min_iter));
            values[std::string(prefix) + "_restart_dimension"] = kdim;
            values[std::string(prefix) + "_minimum_iterations"] = min_iter;
        };
        krylov("gmres", gmres, HYPRE_GMRESGetKDim, HYPRE_GMRESGetMinIter);
        krylov("flexgmres", flex, HYPRE_FlexGMRESGetKDim, HYPRE_FlexGMRESGetMinIter);
        HYPRE_MemoryLocation memory;
        HYPRE_ExecutionPolicy execution;
        checked(HYPRE_GetMemoryLocation(&memory));
        checked(HYPRE_GetExecutionPolicy(&execution));
        Dl_info library{};
        if (!dladdr(reinterpret_cast<void*>(&HYPRE_BoomerAMGCreate), &library) || !library.dli_fname)
            throw std::runtime_error("Library identity unavailable");
        const char* mpi_library = "";
#ifndef HYPRE_SEQUENTIAL
        Dl_info mpi{};
        if (!dladdr(reinterpret_cast<void*>(&MPI_Init), &mpi) || !mpi.dli_fname)
            throw std::runtime_error("MPI identity unavailable");
        mpi_library = mpi.dli_fname;
#endif
        for (const auto& entry : values)
            if (!std::isfinite(entry.second)) throw std::runtime_error("Nonfinite default");
        checked(HYPRE_ParCSRFlexGMRESDestroy(flex));
        checked(HYPRE_ParCSRGMRESDestroy(gmres));
        checked(HYPRE_BoomerAMGDestroy(amg));
        checked(HYPRE_Finalize());
        checked(hypre_MPI_Finalize());

        std::cout << "{\"schema\":\"mars-hypre-defaults-raw-v1\",\"version\":\""
                  << major << '.' << minor << '.' << patch << "\",\"library_path\":"
                  << quoted(library.dli_fname) << ",\"mpi_library_path\":" << quoted(mpi_library)
                  << ",\"build_cuda\":"
#ifdef HYPRE_USING_CUDA
                  << "true"
#else
                  << "false"
#endif
                  << ",\"memory_location\":" << int(memory)
                  << ",\"execution_policy\":" << int(execution)
                  << ",\"header_version_matches\":true,\"getter_layout_checks_passed\":true,\"defaults\":{";
        bool first = true;
        for (const auto& entry : values) {
            std::cout << (first ? "" : ",") << quoted(entry.first.c_str()) << ':'
                      << std::setprecision(17) << entry.second;
            first = false;
        }
        std::cout << "}}\n";
        return 0;
    } catch (...) {
        std::cerr << "Hypre default probe failed; no matrix or solve was requested.\n";
        return 1;
    }
}
