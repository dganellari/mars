#ifndef MAIN_FUNCS_H
#define MAIN_FUNCS_H

#include <cmath>

void addFieldsToMesh(stk::mesh::MetaData& metaData)
{
    add_scalar_field<scalar>("phi", metaData, stk::topology::NODE_RANK);
    add_vector_field<scalar>("phiGrad", metaData, stk::topology::NODE_RANK);
    add_scalar_field<scalar>("beta", metaData, stk::topology::NODE_RANK);
    add_scalar_field<scalar>("gamma", metaData, stk::topology::NODE_RANK);
    add_scalar_N_field<scalar, FEhex<scalar>::numScsIp>(
        "mdot", metaData, stk::topology::ELEMENT_RANK);
    add_scalar_field<scalar>("mdotEdge", metaData, stk::topology::EDGE_RANK);
    add_vector_field<scalar>("areaVEdge", metaData, stk::topology::EDGE_RANK);
    add_scalar_N_field<scalar, FEhex<scalar>::numScsIp * SPATIAL_DIM>(
        "areaV", metaData, stk::topology::ELEMENT_RANK);
    // add_scalar_MxN_field<scalar,
    //                      FEhex<scalar,
    //                      Kokkos::DefaultExecutionSpace>::numScsIp,
    //                      SPATIAL_DIM>(
    //     "areaV", metaData, stk::topology::ELEMENT_RANK); // ip x dim
}

template <typename scalar, size_t BLOCKSIZE>
void getDataSize(size_t& total_bytes,
                 const Coefficients<BLOCKSIZE>& coeffs,
                 const stk::mesh::BulkData& bulkData,
                 const int verbose = 0)
{
    // minimum field bytes accessed (depending on cache/shared memory
    // utilization, this may be more in actual measurement)
    size_t total_field_bytes = 0;

    const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();

    const stk::mesh::Field<scalar>& GammaSTKFieldRef =
        *metaData.get_field<scalar>(stk::topology::NODE_RANK, "gamma");
    const stk::mesh::Field<scalar>& phiSTKFieldRef =
        *metaData.get_field<scalar>(stk::topology::NODE_RANK, "phi");
    const stk::mesh::Field<scalar>& gradPhiSTKFieldRef =
        *metaData.get_field<scalar>(stk::topology::NODE_RANK, "phiGrad");

    // clang-format off
        total_field_bytes += stk::mesh::get_total_ngp_field_allocation_bytes(GammaSTKFieldRef);
        // total_field_bytes += stk::mesh::get_total_ngp_field_allocation_bytes(mdotSTKFieldRef);
        total_field_bytes += stk::mesh::get_total_ngp_field_allocation_bytes(phiSTKFieldRef);
        total_field_bytes += stk::mesh::get_total_ngp_field_allocation_bytes(gradPhiSTKFieldRef);
        // total_field_bytes += stk::mesh::get_total_ngp_field_allocation_bytes(betaSTKFieldRef);
        // total_field_bytes += stk::mesh::get_total_ngp_field_allocation_bytes(coordinatesRef);
    // clang-format on

    const size_t total_buckets =
        bulkData
            .get_buckets(stk::topology::NODE_RANK, metaData.universal_part())
            .size();

    if (verbose > 0)
    {
        std::cout << "Total buckets: " << total_buckets << std::endl;
        std::cout << "Total NGP field bytes (MB): "
                  << total_field_bytes / 1024 / 1024 << std::endl;
    }

    const typename Coefficients<BLOCKSIZE>::Matrix& A = coeffs.getAMatrix();
    const typename Coefficients<BLOCKSIZE>::Vector& b = coeffs.getBVector();

    const size_t total_lin_sys_bytes =
        (A.nnzGlobal() + b.size()) *
        sizeof(typename Coefficients<BLOCKSIZE>::DataType);

    if (verbose > 0)
    {
        std::cout << "Total linear system bytes (MB): "
                  << total_lin_sys_bytes / 1024 / 1024 << std::endl;

        std::cout << "A.nRows() " << A.nRows() << std::endl;
        std::cout << "A.nGlobalRows() " << A.nGlobalRows() << std::endl;
        std::cout << "A.nnzBlocks() " << A.nnzBlocks() << std::endl;
    }

    // NOTE: [fab4100@posteo.net; 2025-05-16] The minimum amount of floating
    // point data bytes required to stream the problem through the hardware.
    // This is the ideal amount of bytes if there was a perfect cache
    // (read/write to DRAM once).
    const size_t n_read = 1;
    const size_t n_write = 1;
    total_bytes = n_read * total_field_bytes + n_write * total_lin_sys_bytes;
}

template <typename Matrix>
typename Matrix::DataType frobeniusNorm(const Matrix& A)
{
    using TReal = typename Matrix::DataType;
    TReal norm = 0;
    const auto values = A.valuesRef();
    for (size_t i = 0; i < values.size(); i++)
    {
        const TReal a_ij = values[i];
        norm += a_ij * a_ij;
    }
    assert(norm > 0);
    return std::sqrt(norm);
}

template <typename scalar,
          typename label,
          size_t BLOCKSIZE,
          size_t SPATIAL_DIM,
          bool includeAdv = true,
          bool isShifted = true>
void runAssembler(
    const std::unique_ptr<
        AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM, includeAdv, isShifted>>&
        assembler,
    Coefficients<BLOCKSIZE>& coeffs,
    const stk::mesh::MetaData& metaData,
    const stk::mesh::BulkData& bulkData,
    const size_t total_bytes,
    const int verbose = 0)
{
    // reset A and b
    coeffs.zeroALL();

    typename Coefficients<BLOCKSIZE>::Matrix& A = coeffs.getAMatrix();
    typename Coefficients<BLOCKSIZE>::Vector& b = coeffs.getBVector();

    /******* FIRST PASS (warmup)  *******/
    auto startTime = std::chrono::high_resolution_clock::now();

    assembler->assemble(A, b, metaData, bulkData);

    Kokkos::fence();

    auto stopTime = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
        stopTime - startTime);

    if (verbose > 0)
    {
        std::cout << "Assembler: " << std::setw(32) << std::left
                  << std::setfill('.') << assembler->name() << ": "
                  << std::setfill(' ') << std::scientific
                  << duration.count() * 1.0e-3 << " milliseconds (warmup)\n";
    }

    A.syncToHost();
#ifdef USE_KOKKOS
    b.syncToHost();
#endif

    const auto matrix_norm = frobeniusNorm(A);

    /******* SECOND PASS  *******/
    double t_accum = 0;
    for (int i = 0; i < n_reps; i++)
    {
        // reset A and b
        coeffs.zeroALL();

        startTime = std::chrono::high_resolution_clock::now();

        assembler->assemble(A, b, metaData, bulkData);

        Kokkos::fence();

        stopTime = std::chrono::high_resolution_clock::now();

        duration = std::chrono::duration_cast<std::chrono::microseconds>(
            stopTime - startTime);
        t_accum += duration.count();
    }
    std::cout << "Assembler: " << std::setw(32) << std::left
              << std::setfill('.') << assembler->name() << ": "
              << std::setfill(' ') << std::scientific
              << t_accum * 1.0e-3 / n_reps << " milliseconds @ "
              << n_reps * total_bytes * 1.0e-9 / (t_accum * 1.0e-6)
              << " GB/s (average of " << n_reps
              << " samples) [matrix norm: " << matrix_norm << "]\n";

    std::cout.unsetf(std::ios::fixed);

    A.syncToHost();
#ifdef USE_KOKKOS
    b.syncToHost();
#endif
    if (verbose > 1)
    {
        A.dump(15, 15, 8, 2);
        for (int i = 0; i < 15; ++i)
        {
            std::cout << i << "\t" << b[i] << std::endl;
        }
    }
}

std::set<std::string> commaListToSet(const std::string& input)
{
    std::set<std::string> result;
    std::stringstream ss(input);
    std::string item;
    while (std::getline(ss, item, ','))
    {
        result.insert(item);
    }
    return result;
}

void parseOptions(int argc,
                  char* argv[],
                  std::string& baseName,
                  std::set<std::string>& assemblerNames,
                  int& verbose,
                  int& n_reps)
{
    // clang-format off
    try
    {
        cxxopts::Options options(
           argv[0], "Description"); //"Example program with optional -a argument");

        // add optionals
        options.add_options()
            ("a,assemblers","Comma-separated list of assembler names",cxxopts::value<std::string>())
            ("v,verbose", "Verbosity level (default 0)",cxxopts::value<int>())
            ("r,repetitions", "Number of repetitions for all kernels (default 10)",cxxopts::value<int>())
            ("h,help", "Print help");

        // add positionals
        options.add_options()
            ("meshBaseName","Mesh Base Name (without .exo)",cxxopts::value<std::string>());

        options.positional_help("<meshBaseName>");

        options.parse_positional({"meshBaseName"});

        auto result = options.parse(argc, argv);

        if (result.count("help"))
        {
            std::cout << options.help() << std::endl;
            exit(0);
        }

        if (result.count("verbose"))
        {
            verbose = result["verbose"].as<int>();
        }
        
        if (result.count("repetitions"))
        {
            n_reps = result["repetitions"].as<int>();
            if (n_reps <= 0) 
            {
                throw cxxopts::exceptions::exception("Repetitions must be > 0");
            }
        }

        if (result.count("assemblers"))
        {
            assemblerNames.clear();
            assemblerNames =
                commaListToSet(result["assemblers"].as<std::string>());
        }

        baseName = result["meshBaseName"].as<std::string>();
    }
    catch (const cxxopts::exceptions::exception& e)
    {
        std::cerr << "Error parsing options: " << e.what() << std::endl;
        exit(1);
    }
    // clang-format on

    // if (argc < 2)
    // {
    //     std::cerr << "Usage: " << argv[0] << " <baseName>\n";
    //     exit(1);
    // }

    // check validity of assembler names
    for (auto& asmName : assemblerNames)
    {
        if (strToAsm(asmName) == AsmType::NONE)
        {
            exit(1);
        }
    }

    if (verbose > 0)
    {
        std::cout << "Running the following assemblers: " << std::endl;
        for (auto& asmName : assemblerNames)
        {
            std::cout << "    " << asmName << std::endl;
        }
        std::cout << std::endl;

        // std::string baseName = argv[1];
        std::cout << "BaseName: " << baseName << "\n";
        std::cout << "Mesh File: " << baseName + ".exo" << "\n";
        std::cout << "Graph File: " << baseName + ".bin" << "\n";
    }
}
#endif // MAIN_FUNCS_H
