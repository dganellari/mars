#pragma once

namespace heat
{

using cellwise::topo_detail::raw;

constexpr double hot_temperature = 353.15, cold_temperature = 293.15, conductivity = 1.0;

__global__ void tag_nodes(const int* nodes, long long count, unsigned bit, unsigned* tags, long long size, int* bad)
{
    const long long i = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= count) return;
    const int node = nodes[i];
    if (node < 0 || node >= size)
        atomicExch(bad, 1);
    else
        atomicOr(tags + node, bit);
}

__global__ void classify_faces(cellwise::domain_detail::Corners<uint64_t> ids,
                               long long first,
                               long long count,
                               const unsigned char* mirror,
                               const unsigned* tags,
                               const unsigned char* code,
                               unsigned char* labels,
                               int* bad)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= count * 6) return;
    labels[t] = 0;
    if (!(code[t] & 64)) return;
    const long long e = t / 6;
    unsigned bits     = 7;
    for (int j = 0; j < 4; ++j)
    {
        const int k        = cellwise::face_corner(t % 6, j);
        const int original = mirror[e] && (k & 1) ? k ^ 2 : k;
        bits &= tags[ids.key[original][first + e]];
    }
    // A physical face needs one complete label; corner-only intersections do not define a face.
    if (bits != 1 && bits != 2 && bits != 4)
        atomicExch(bad, 1);
    else
        labels[t] = bits == 1 ? 1 : bits == 2 ? 2 : 3;
}

void boundary_faces(MeshDomain& domain,
                    cellwise::domain_detail::Corners<uint64_t> ids,
                    const thrust::device_vector<unsigned char>& mirror,
                    long long first,
                    long long count,
                    const cellwise::UnstructuredTopology& topology,
                    const std::string& model,
                    thrust::device_vector<unsigned char>& labels)
{
    try
    {
        const std::map<std::string, unsigned> prescribed =
            model == "wing" ? std::map<std::string, unsigned>{{"blade", 1}, {"bottom", 2}, {"in", 2}, {"out", 2},
                                                              {"side", 2},  {"top", 2},    {"sym", 4}}
                            : std::map<std::string, unsigned>{{"hot", 1}, {"cold", 2}, {"insulated", 4}};
        const auto names = domain.sideSetNames();
        if (names.size() != prescribed.size()) throw std::runtime_error("Unexpected thermal boundary sets");
        thrust::device_vector<unsigned> tags(domain.getNodeX().size(), 0);
        thrust::device_vector<int> bad(1, 0);
        for (const auto& name : names)
        {
            const auto found = prescribed.find(name);
            if (found == prescribed.end()) throw std::runtime_error("Unknown thermal boundary: " + name);
            const auto& nodes = domain.sideSetNodes(name);
            if (!nodes.empty())
                tag_nodes<<<blocks_for(nodes.size()), cellwise::kThreads>>>(thrust::raw_pointer_cast(nodes.data()),
                                                                            nodes.size(), found->second, raw(tags),
                                                                            tags.size(), raw(bad));
        }
        labels.resize(count * 6);
        if (count)
            classify_faces<<<blocks_for(count * 6), cellwise::kThreads>>>(
                ids, first, count, raw(mirror), raw(tags), raw(topology.face_code), raw(labels), raw(bad));
        BP_CK(cudaGetLastError());
        BP_CK(cudaDeviceSynchronize());
        int invalid = bad[0];
        MARS_CELLWISE_MPI(MPI_Allreduce(MPI_IN_PLACE, &invalid, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD));
        if (invalid) throw std::runtime_error("Missing, overlapping or invalid thermal face labels");
        long long counts[3];
        for (int i = 0; i < 3; ++i)
            counts[i] = thrust::count(labels.begin(), labels.end(), i + 1);
        MARS_CELLWISE_MPI(MPI_Allreduce(MPI_IN_PLACE, counts, 3, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD));
        if (!counts[0] || !counts[1] || !counts[2]) throw std::runtime_error("An essential thermal boundary is empty");
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (!rank) printf("THERMAL FACES: hot=%lld cold=%lld insulated=%lld\n", counts[0], counts[1], counts[2]);
    }
    catch (const std::exception& error)
    {
        fprintf(stderr, "Thermal boundary setup: %s\n", error.what());
        MPI_Abort(MPI_COMM_WORLD, 92);
    }
}

template<int P>
__global__ void node_marks(long long elements, const unsigned char* faces, double* hot, double* cold)
{
    constexpr int n = P + 1, n3 = n * n * n;
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= elements * n3) return;
    const int l = t % n3, xyz[3] = {l / (n * n), (l / n) % n, l % n};
    bool h = false, c = false;
    for (int a = 0; a < 3; ++a)
    {
        if (xyz[a] != 0 && xyz[a] != P) continue;
        const unsigned char label = faces[(t / n3) * 6 + 2 * a + (xyz[a] == P)];
        h |= label == 1;
        c |= label == 2;
    }
    hot[t]  = h;
    cold[t] = c;
}

__global__ void set_face_flags(long long n, const unsigned char* labels, unsigned char* codes)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t < n) codes[t] = (codes[t] & 63) | ((labels[t] == 1 || labels[t] == 2) ? 64 : 0);
}

template<int P>
__global__ void set_star_flags(cellwise::TopologyView T,
                               long long edges,
                               long long vertices,
                               const double* hot,
                               const double* cold,
                               unsigned char* edge_flags,
                               unsigned char* vertex_flags)
{
    constexpr int n3  = (P + 1) * (P + 1) * (P + 1);
    const long long g = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (g < edges)
    {
        bool dir = false;
        for (int i = T.edge_off[g]; i < T.edge_off[g + 1]; ++i)
        {
            const int entry = T.edge_ent[i] & 0x7fffffff, e = entry / 12;
            if (e < T.local)
            {
                const long long node = (long long)e * n3 + cellwise::edge_node<P>(entry % 12, 1);
                dir |= hot[node] > 0 || cold[node] > 0;
            }
        }
        edge_flags[g] = dir;
    }
    if (g < vertices)
    {
        bool dir = false;
        for (int i = T.vert_off[g]; i < T.vert_off[g + 1]; ++i)
        {
            const int entry = T.vert_ent[i], e = entry / 8;
            if (e < T.local)
            {
                const long long node = (long long)e * n3 + cellwise::corner_node<P>(entry % 8);
                dir |= hot[node] > 0 || cold[node] > 0;
            }
        }
        vertex_flags[g] = dir;
    }
}

__global__ void make_lift(long long n, const double* hot, const double* cold, double* lift, int* bad)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= n) return;
    if (hot[t] > 0 && cold[t] > 0) atomicExch(bad, 1);
    lift[t] = hot[t] > 0 ? 1.0 : 0.0;
}

template<int P>
void prepare(cellwise::UnstructuredTopology& T,
             cellwise::UnstructuredHalo* halo,
             const thrust::device_vector<unsigned char>& faces,
             thrust::device_vector<double>& hot,
             thrust::device_vector<double>& cold,
             thrust::device_vector<double>& scratch,
             double* lift,
             cudaStream_t stream)
{
    const long long n = T.local * (P + 1) * (P + 1) * (P + 1);
    if (n) node_marks<P><<<blocks_for(n), cellwise::kThreads, 0, stream>>>(T.local, raw(faces), raw(hot), raw(cold));
    BP_CK(cudaGetLastError());
    // Union boundary marks across every copy before changing the homogeneous correction masks.
    cellwise::dss<P>(raw(hot), raw(scratch), T, halo, stream);
    hot.swap(scratch);
    cellwise::dss<P>(raw(cold), raw(scratch), T, halo, stream);
    cold.swap(scratch);
    if (T.local)
        set_face_flags<<<blocks_for(T.local * 6), cellwise::kThreads, 0, stream>>>(T.local * 6, raw(faces),
                                                                                   raw(T.face_code));
    const long long stars = std::max(T.edges(), T.vertices());
    if (stars)
        set_star_flags<P><<<blocks_for(stars), cellwise::kThreads, 0, stream>>>(
            cellwise::view(T), T.edges(), T.vertices(), raw(hot), raw(cold), raw(T.edge_dirichlet),
            raw(T.vert_dirichlet));
    thrust::device_vector<int> bad(1, 0);
    BP_CK(cudaStreamSynchronize(0));
    if (n) make_lift<<<blocks_for(n), cellwise::kThreads, 0, stream>>>(n, raw(hot), raw(cold), lift, raw(bad));
    BP_CK(cudaGetLastError());
    BP_CK(cudaStreamSynchronize(stream));
    int invalid = bad[0];
    MARS_CELLWISE_MPI(MPI_Allreduce(MPI_IN_PLACE, &invalid, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD));
    if (invalid)
    {
        fprintf(stderr, "Thermal hot and cold boundaries share a node\n");
        MPI_Abort(MPI_COMM_WORLD, 92);
    }
}

template<int P>
__device__ bool dirichlet(long long t, cellwise::TopologyView T)
{
    constexpr int n = P + 1, n3 = n * n * n;
    const long long e = t / n3;
    const int l       = t % n3;
    const auto entity = cellwise::local_entity<P>(l / (n * n), (l / n) % n, l % n);
    if (entity.kind == 1) return T.face_code[e * 6 + entity.index] & 64;
    if (entity.kind == 2) return T.edge_dirichlet[T.edge_of[e * 12 + entity.index] & 0x7fffffff];
    if (entity.kind == 3) return T.vert_dirichlet[T.vert_of[e * 8 + entity.index]];
    return false;
}

__global__ void add_lift(double* x, const double* lift, long long n)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t < n) x[t] += lift[t];
}

struct Sums
{
    double hot, cold, free, energy, bc, low, high, outside, nodes;
};
struct AddSums
{
    __host__ __device__ Sums operator()(Sums a, Sums b) const
    {
        return {a.hot + b.hot,      a.cold + b.cold,      a.free + b.free,       a.energy + b.energy, fmax(a.bc, b.bc),
                fmin(a.low, b.low), fmax(a.high, b.high), a.outside + b.outside, a.nodes + b.nodes};
    }
};
template<int P>
struct Measure
{
    cellwise::TopologyView T;
    const double *x, *lift, *action;
    __device__ Sums operator()(long long t) const
    {
        const bool dir = dirichlet<P>(t, T), hot = lift[t] > 0;
        constexpr int n = P + 1, n3 = n * n * n;
        const int l       = t % n3;
        const auto entity = cellwise::local_entity<P>(l / (n * n), (l / n) % n, l % n);
        const bool owned  = (T.counting[t / n3] >> cellwise::entity_bit(entity)) & 1u;
        // Count each global node once; overshoot is a diagnostic, not a Q3 correctness gate.
        const bool outside = x[t] < -1e-6 || x[t] > 1 + 1e-6;
        return {hot ? action[t] : 0,
                dir && !hot ? action[t] : 0,
                !dir ? action[t] : 0,
                x[t] * action[t],
                dir ? fabs(x[t] - lift[t]) : 0,
                x[t],
                x[t],
                owned && outside ? 1.0 : 0.0,
                owned ? 1.0 : 0.0};
    }
};

template<int P>
Sums measure(const cellwise::UnstructuredTopology& T,
             const double* x,
             const double* lift,
             const double* action,
             cudaStream_t stream)
{
    const long long n = T.local * (P + 1) * (P + 1) * (P + 1);
    Sums local        = thrust::transform_reduce(
        thrust::cuda::par.on(stream), thrust::make_counting_iterator(0LL), thrust::make_counting_iterator(n),
        Measure<P>{cellwise::view(T), x, lift, action}, Sums{0, 0, 0, 0, 0, 1e300, -1e300, 0, 0}, AddSums{});
    double sum[6]    = {local.hot, local.cold, local.free, local.energy, local.outside, local.nodes},
           maxima[3] = {local.bc, -local.low, local.high};
    MARS_CELLWISE_MPI(MPI_Allreduce(MPI_IN_PLACE, sum, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
    MARS_CELLWISE_MPI(MPI_Allreduce(MPI_IN_PLACE, maxima, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD));
    return {sum[0], sum[1], sum[2], sum[3], maxima[0], -maxima[1], maxima[2], sum[4], sum[5]};
}

template<int P>
__global__ void slab_error(const double* x, const double* corners, double* error, long long elements)
{
    constexpr int n = P + 1, n3 = n * n * n;
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= elements * n3) return;
    const int l       = t % n3;
    const double r[3] = {c_gll[l / (n * n)], c_gll[(l / n) % n], c_gll[l % n]};
    double physical_x = 0;
    for (int k = 0; k < 8; ++k)
    {
        double weight = 1;
        for (int a = 0; a < 3; ++a)
            weight *= .5 * (1 + r[a] * (2 * ((k >> (2 - a)) & 1) - 1));
        physical_x += weight * corners[(t / n3) * 24 + k];
    }
    error[t] = fabs(x[t] - (1 - physical_x));
}

} // namespace heat
