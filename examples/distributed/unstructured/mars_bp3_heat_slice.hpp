// Optional post-solve output: keep only owned elements intersecting the requested window.
struct SliceCell
{
    unsigned long long key;
    double corners[24], value[64], lift[64];
};
static_assert(sizeof(SliceCell) == 153 * sizeof(double));

struct SliceWindow
{
    double xlo, xhi, ylo, yhi, z;
    const double* corners;
    __device__ bool operator()(long long e) const
    {
        double lo[3], hi[3];
        for (int a = 0; a < 3; ++a)
        {
            lo[a] = hi[a] = corners[e * 24 + a * 8];
            for (int k = 1; k < 8; ++k)
            {
                const double v = corners[e * 24 + a * 8 + k];
                lo[a]          = fmin(lo[a], v);
                hi[a]          = fmax(hi[a], v);
            }
        }
        return lo[0] <= xhi && hi[0] >= xlo && lo[1] <= yhi && hi[1] >= ylo && lo[2] <= z && hi[2] >= z;
    }
};

__global__ void slice_gather(const long long* ids,
                             long long count,
                             const unsigned long long* keys,
                             const double* corners,
                             const double* solution,
                             const double* lift,
                             SliceCell* out)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= count * 64) return;
    const long long i = t / 64, e = ids[i];
    const int node = t % 64;
    if (node == 0) out[i].key = keys[e];
    if (node < 24) out[i].corners[node] = corners[e * 24 + node];
    out[i].value[node] = solution[e * 64 + node];
    out[i].lift[node]  = lift[e * 64 + node];
}

void export_heat_slice(const std::string& prefix,
                       SliceWindow window,
                       long long elements,
                       const unsigned long long* keys,
                       const double* corners,
                       const double* solution,
                       const double* lift,
                       const std::vector<double>& nodes,
                       double true_residual,
                       double boundary_error,
                       cudaStream_t stream)
{
    int rank = 0, ranks = 0;
    MARS_CELLWISE_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
    MARS_CELLWISE_MPI(MPI_Comm_size(MPI_COMM_WORLD, &ranks));
    window.corners = corners;
    thrust::device_vector<long long> ids(elements);
    auto first            = thrust::make_counting_iterator<long long>(0);
    auto end              = thrust::copy_if(mars::gpuPar.on(stream), first, first + elements, ids.begin(), window);
    const long long count = end - ids.begin();
    // Bound temporary storage and host transfers independently of the full mesh size.
    constexpr long long chunk = 8192;
    thrust::device_vector<SliceCell> d_cells(std::min(count, chunk));
    std::vector<SliceCell> h_cells(std::min(count, chunk));
    BP_CK(cudaStreamSynchronize(0));
    const std::string name = prefix + "-rank" + std::to_string(rank);
    FILE* file             = fopen((name + ".bin").c_str(), "wb");
    if (!file) throw std::runtime_error("Cannot open slice output: " + name);
    bool ok = true;
    for (long long start = 0; start < count && ok; start += chunk)
    {
        const long long n = std::min(chunk, count - start);
        slice_gather<<<blocks_for(n * 64), cellwise::kThreads, 0, stream>>>(
            thrust::raw_pointer_cast(ids.data()) + start, n, keys, corners, solution, lift,
            thrust::raw_pointer_cast(d_cells.data()));
        BP_CK(cudaGetLastError());
        // This opt-in host copy is for file output only, after all solve timers stop.
        BP_CK(cudaMemcpyAsync(h_cells.data(), thrust::raw_pointer_cast(d_cells.data()), n * sizeof(SliceCell),
                              cudaMemcpyDeviceToHost, stream));
        BP_CK(cudaStreamSynchronize(stream));
        ok = fwrite(h_cells.data(), sizeof(SliceCell), n, file) == (size_t)n;
    }
    if (fclose(file) != 0) ok = false;
    if (!ok) throw std::runtime_error("Cannot write slice output: " + name);
    file = fopen((name + ".json").c_str(), "w");
    if (!file) throw std::runtime_error("Cannot open slice metadata: " + name);
    const int wrote  = fprintf(file,
                               "{\"format\":\"mars-heat-slice-v1\",\"p\":3,\"rank\":%d,\"ranks\":%d,\"elements\":%lld,"
                                "\"record_bytes\":%zu,\"window\":[%.17g,%.17g,%.17g,%.17g,%.17g],"
                                "\"gll\":[%.17g,%.17g,%.17g,%.17g],\"true_residual\":%.17g,\"boundary_error\":%.17g}\n",
                               rank, ranks, count, sizeof(SliceCell), window.xlo, window.xhi, window.ylo, window.yhi,
                               window.z, nodes[0], nodes[1], nodes[2], nodes[3], true_residual, boundary_error);
    const int closed = fclose(file);
    if (wrote < 0 || closed != 0) throw std::runtime_error("Cannot write slice metadata: " + name);
    long long total = 0;
    MARS_CELLWISE_MPI(MPI_Allreduce(&count, &total, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD));
    if (total == 0) throw std::runtime_error("Slice window selected no elements");
    if (rank == 0)
        printf("SLICE EXPORT: %lld owned elements on %d ranks; computed dimensionless temperature, boundary lift "
               "stored separately\n",
               total, ranks);
}
