#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>
#include <climits>
#include <cmath>
#include <exception>
#include <set>
#include <mpi.h>
#include <thrust/gather.h>
#include <thrust/device_ptr.h>
#include <thrust/execution_policy.h>
#include "mars_gpu_runtime.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_element_traits.hpp"

namespace mars
{
namespace fem
{
namespace detail
{
inline size_t sfemGatherDisplacements(const std::vector<int>& counts, std::vector<int>& displacements)
{
    size_t total = 0;
    displacements.resize(counts.size());
    for (size_t r = 0; r < counts.size(); ++r)
    {
        if (counts[r] < 0 || size_t(counts[r]) > size_t(INT_MAX) - total)
            throw std::runtime_error("SfemWriter: count exceeds MPI gather capacity");
        displacements[r] = int(total);
        total += size_t(counts[r]);
    }
    return total;
}
} // namespace detail

// Writes a distributed mesh and its nodal fields in SFEM's raw format (the layout of smesh's
// mesh_write and sfem::Output), gathered on rank 0 into one serial numbering:
//
//   <dir>/mesh/meta.yaml, x.float32, y.float32, z.float32, i0.int32 .. i<NPE-1>.int32
//   <dir>/out/<name>.<step %09d>.<type>               (scalar)
//   <dir>/out/<name>.<component>.<step %09d>.<type>   (3-vector, one file per component)
//   <dir>/out/time.txt                                (one line per written step)
//
// The mesh is written once, on the first writeFields; each frame then holds only the fields'
// values (no mesh, no ghosts), as float32 unless full precision is asked for. Nodes are numbered
// by their global SFC key (ascending), each written once by its owner (ownership 1); elements are
// the ranks' owned ranges [startIndex, endIndex), in rank order. This bounded gather format requires fewer than 2^31
// nodes and elements, and enough rank-0 memory for the gathered mesh and fields.
// Coincident entities that need distinct identities are unsupported. Read the mesh with
// smesh/SFEM or numpy.fromfile.
//
// Collective: every rank calls writeFields with the same fields. The domain layout must stay fixed.
template<typename KeyType, typename RealType, typename ElementTag = HexTag>
class SfemWriter
{
public:
    using DomainT = ElementDomain<ElementTag, RealType, KeyType, cstone::execution::Gpu>;

    static constexpr int NPE = ElemTraits<ElementTag>::NodesPerElem;

    struct Field
    {
        std::string name;
        const cstone::DeviceVector<RealType>* x; // per local node (nodeCount)
        const cstone::DeviceVector<RealType>* y; // null for a scalar
        const cstone::DeviceVector<RealType>* z; // null for a scalar
    };

    explicit SfemWriter(const std::string& dir, bool full_precision = false)
        : dir_(dir)
        , full_precision_(full_precision)
    {
    }

    // A resumed run (a restart): time.txt is kept and appended to, not started afresh. The mesh is
    // written again on the first writeFields (the same mesh: the same files).
    void appendTimes() { time_started_ = true; }

    // Writes <dir>/out/<field>[.c].<step>.<type> for every field, and appends t to time.txt.
    void writeFields(int step, double t, const DomainT& domain, const std::vector<Field>& fields)
    {
        // A root I/O or allocation failure must stop peers already entering their next gather.
        try
        {
            writeFieldsImpl(step, t, domain, fields);
        }
        catch (const std::exception& e)
        {
            fprintf(stderr, "SfemWriter: collective output failed: %s\n", e.what());
            MPI_Abort(MPI_COMM_WORLD, 1);
            std::terminate();
        }
        catch (...)
        {
            fprintf(stderr, "SfemWriter: collective output failed\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
            std::terminate();
        }
    }

private:
    static void checkMpi(int error)
    {
        if (error != MPI_SUCCESS) throw std::runtime_error("SfemWriter: MPI operation failed");
    }

    static void requireCollectively(bool valid)
    {
        int local = valid ? 0 : 1, global = 0;
        checkMpi(MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD));
        if (global) throw std::runtime_error("SfemWriter: invalid collective output contract");
    }

    void writeFieldsImpl(int step, double t, const DomainT& domain, const std::vector<Field>& fields)
    {
        bool valid = step >= 0 && std::isfinite(t) && !dir_.empty() && dir_.find('\0') == std::string::npos &&
                     domain.getNodeCount() <= size_t(INT_MAX) &&
                     (!mesh_written_ || (&domain == domain_ && domain.getNodeCount() == local_nodes_));
        std::string schema = dir_ + "\n" + std::to_string(step) + "\n" + std::to_string(full_precision_) + "\n" +
                             std::to_string(mesh_written_) + "\n" + std::to_string(time_started_);
        std::set<std::string> names;
        std::set<std::string> output_names;
        for (const auto& f : fields)
        {
            valid = valid && names.insert(f.name).second && !f.name.empty() &&
                    f.name.find_first_of("/\\\n\r") == std::string::npos && f.name.find('\0') == std::string::npos &&
                    f.name != "." && f.name != ".." && f.x && bool(f.y) == bool(f.z);
            const cstone::DeviceVector<RealType>* components[3] = {f.x, f.y, f.z};
            for (auto c : components)
                if (c) valid = valid && c->size() == domain.getNodeCount();
            const int ncomp = f.y ? 3 : 1;
            for (int c = 0; c < ncomp; ++c)
                valid = output_names.insert(f.name + (ncomp == 1 ? "" : "." + std::to_string(c))).second && valid;
            schema += "\n" + f.name + (f.y ? ":3" : ":1");
        }
        valid = valid && schema.size() <= size_t(INT_MAX);
        requireCollectively(valid);
        int rank = 0;
        checkMpi(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
        int length = int(schema.size());
        checkMpi(MPI_Bcast(&length, 1, MPI_INT, 0, MPI_COMM_WORLD));
        std::string expected(size_t(length), '\0');
        if (rank == 0) expected = schema;
        checkMpi(MPI_Bcast(expected.data(), length, MPI_CHAR, 0, MPI_COMM_WORLD));
        double expectedTime = t;
        checkMpi(MPI_Bcast(&expectedTime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD));
        requireCollectively(schema == expected && t == expectedTime);
        if (!mesh_written_) writeMesh(domain);

        std::vector<RealType> h(owned_.size());
        std::vector<double> packed(owned_.size());
        std::vector<double> gathered(rank == 0 ? n_nodes_ : 0);
        for (const auto& f : fields)
        {
            const cstone::DeviceVector<RealType>* comp[3] = {f.x, f.y, f.z};
            const int ncomp                               = f.y ? 3 : 1;
            for (int c = 0; c < ncomp; ++c)
            {
                copyOwned(*comp[c], h);
                std::copy(h.begin(), h.end(), packed.begin());
                gatherToRoot(packed, gathered);
                if (rank != 0) continue;

                char frame[16];
                snprintf(frame, sizeof(frame), "%09d", step);
                const std::string name = dir_ + "/out/" + f.name + (ncomp == 1 ? "" : "." + std::to_string(c)) + "." +
                                         frame + (full_precision_ ? ".float64" : ".float32");
                if (full_precision_)
                {
                    std::vector<double> v(n_nodes_);
                    for (size_t i = 0; i < n_nodes_; ++i)
                        v[order_[i]] = gathered[i];
                    writeRaw(name.c_str(), v.data(), v.size());
                }
                else
                {
                    std::vector<float> v(n_nodes_);
                    for (size_t i = 0; i < n_nodes_; ++i)
                        v[order_[i]] = float(gathered[i]);
                    writeRaw(name.c_str(), v.data(), v.size());
                }
            }
        }

        if (rank == 0)
        {
            const std::string path = dir_ + "/out/time.txt";
            FILE* f                = fopen(path.c_str(), time_started_ ? "a" : "w");
            if (!f) throw std::runtime_error("SfemWriter: cannot open " + path);
            const bool failed = fprintf(f, "%.17g\n", t) < 0;
            if (fclose(f) != 0 || failed) throw std::runtime_error("SfemWriter: time write failed");
        }
        // A successful return means root has closed every file in this frame.
        checkMpi(MPI_Barrier(MPI_COMM_WORLD));
        time_started_ = true;
    }

private:
    std::string dir_;
    bool full_precision_;
    bool mesh_written_     = false;
    bool time_started_     = false;
    const DomainT* domain_ = nullptr;
    size_t local_nodes_    = 0;

    cstone::DeviceVector<int> d_owned_;
    cstone::DeviceVector<RealType> d_field_;
    std::vector<int> owned_;    // this rank's owned local nodes, in the order they are sent
    std::vector<int> counts_;   // rank 0: owned nodes per rank
    std::vector<int> displs_;   // rank 0
    std::vector<size_t> order_; // rank 0: gathered position -> serial node number
    size_t n_nodes_ = 0;        // rank 0: serial node count

    template<typename T>
    static void writeRaw(const char* path, const T* data, size_t n)
    {
        FILE* f = fopen(path, "wb");
        if (!f) throw std::runtime_error(std::string("SfemWriter: cannot open ") + path);
        if (n && fwrite(data, sizeof(T), n, f) != n)
        {
            fclose(f);
            throw std::runtime_error(std::string("SfemWriter: cannot write ") + path);
        }
        if (fclose(f) != 0) throw std::runtime_error("SfemWriter: closing output failed");
    }

    template<typename T>
    static std::vector<T> toHost(const cstone::DeviceVector<T>& d, size_t n)
    {
        std::vector<T> h(n);
        if (n > d.size()) throw std::runtime_error("SfemWriter: short device array");
        if (n && cudaMemcpy(h.data(), d.data(), n * sizeof(T), cudaMemcpyDeviceToHost) != cudaSuccess)
            throw std::runtime_error("SfemWriter: mesh transfer failed");
        return h;
    }

    void copyOwned(const cstone::DeviceVector<RealType>& values, std::vector<RealType>& h)
    {
        h.resize(owned_.size());
        if (h.empty()) return;
        // Only owned values need to cross to the host filesystem API.
        thrust::gather(thrust::device, thrust::device_pointer_cast(d_owned_.data()),
                       thrust::device_pointer_cast(d_owned_.data() + d_owned_.size()),
                       thrust::device_pointer_cast(values.data()), thrust::device_pointer_cast(d_field_.data()));
        if (cudaMemcpy(h.data(), d_field_.data(), h.size() * sizeof(RealType), cudaMemcpyDeviceToHost) != cudaSuccess)
            throw std::runtime_error("SfemWriter: field transfer failed");
    }

    // Owned values (this rank's owned_ order) to rank 0, concatenated in rank order.
    void gatherToRoot(const std::vector<double>& mine, std::vector<double>& all) const
    {
        checkMpi(MPI_Gatherv(mine.data(), int(mine.size()), MPI_DOUBLE, all.data(), counts_.data(), displs_.data(),
                             MPI_DOUBLE, 0, MPI_COMM_WORLD));
    }

    void writeMesh(const DomainT& domain)
    {
        int rank = 0, size = 1;
        checkMpi(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
        checkMpi(MPI_Comm_size(MPI_COMM_WORLD, &size));

        const size_t nodeCount          = domain.getNodeCount();
        const std::vector<uint8_t> own  = toHost(domain.getNodeOwnershipMap(), nodeCount);
        const std::vector<KeyType> keys = toHost(domain.getLocalToGlobalSfcMap(), nodeCount);

        // On one rank every node is owned, whatever the map says.
        owned_.clear();
        for (size_t i = 0; i < nodeCount; ++i)
            if (size == 1 || own[i] == 1) owned_.push_back(int(i));

        d_owned_.resize(owned_.size());
        d_field_.resize(owned_.size());
        if (!owned_.empty() && cudaMemcpy(d_owned_.data(), owned_.data(), owned_.size() * sizeof(int),
                                          cudaMemcpyHostToDevice) != cudaSuccess)
            throw std::runtime_error("SfemWriter: ownership transfer failed");

        requireCollectively(domain.getNodeX().size() == nodeCount && domain.getNodeY().size() == nodeCount &&
                            domain.getNodeZ().size() == nodeCount);
        std::vector<RealType> x, y, z;
        copyOwned(domain.getNodeX(), x);
        copyOwned(domain.getNodeY(), y);
        copyOwned(domain.getNodeZ(), z);

        // Node keys and coordinates to rank 0.
        const int n_owned = int(owned_.size());
        counts_.assign(rank == 0 ? size : 0, 0);
        checkMpi(MPI_Gather(&n_owned, 1, MPI_INT, counts_.data(), 1, MPI_INT, 0, MPI_COMM_WORLD));
        n_nodes_ = rank == 0 ? detail::sfemGatherDisplacements(counts_, displs_) : 0;

        std::vector<double> packed(owned_.size());
        std::vector<double> gx(n_nodes_), gy(n_nodes_), gz(n_nodes_);
        for (int i = 0; i < n_owned; ++i)
            packed[i] = x[i];
        gatherToRoot(packed, gx);
        for (int i = 0; i < n_owned; ++i)
            packed[i] = y[i];
        gatherToRoot(packed, gy);
        for (int i = 0; i < n_owned; ++i)
            packed[i] = z[i];
        gatherToRoot(packed, gz);
        std::vector<uint64_t> mykeys(n_owned), gk(n_nodes_);
        for (int i = 0; i < n_owned; ++i)
            mykeys[i] = uint64_t(keys[owned_[i]]);
        checkMpi(MPI_Gatherv(mykeys.data(), n_owned, MPI_UINT64_T, gk.data(), counts_.data(), displs_.data(),
                             MPI_UINT64_T, 0, MPI_COMM_WORLD));

        // Owned elements' connectivity, as global keys, to rank 0.
        const auto& cs  = domain.getDomain();
        const size_t e0 = size_t(cs.startIndex()), e1 = size_t(cs.endIndex());
        const auto& d_conn = domain.getElementToNodeConnectivity();
        std::vector<KeyType> conn[NPE];
        copyConn(conn, d_conn, std::get<0>(d_conn).size(), std::make_index_sequence<NPE>{});

        requireCollectively(e1 >= e0 && e1 <= std::get<0>(d_conn).size() && e1 - e0 <= size_t(INT_MAX));
        const int n_elem = int(e1 - e0);
        std::vector<int> ecounts(rank == 0 ? size : 0), edispls(rank == 0 ? size : 0);
        checkMpi(MPI_Gather(&n_elem, 1, MPI_INT, ecounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD));
        const size_t n_elements = rank == 0 ? detail::sfemGatherDisplacements(ecounts, edispls) : 0;

        std::vector<uint64_t> mine(n_elem), all(n_elements);
        std::vector<std::vector<uint64_t>> gconn(NPE);
        for (int k = 0; k < NPE; ++k)
        {
            for (int e = 0; e < n_elem; ++e)
            {
                const size_t node = size_t(conn[k][e0 + e]);
                if (node >= keys.size()) throw std::runtime_error("SfemWriter: invalid element node index");
                mine[e] = uint64_t(keys[node]);
            }
            checkMpi(MPI_Gatherv(mine.data(), n_elem, MPI_UINT64_T, all.data(), ecounts.data(), edispls.data(),
                                 MPI_UINT64_T, 0, MPI_COMM_WORLD));
            if (rank == 0) gconn[k] = all;
        }

        if (rank == 0)
        {
            // Serial numbering: nodes by ascending key.
            const std::vector<uint64_t>& k = gk;
            std::vector<size_t> by_key(n_nodes_);
            for (size_t i = 0; i < n_nodes_; ++i)
                by_key[i] = i;
            std::sort(by_key.begin(), by_key.end(), [&k](size_t a, size_t b) { return k[a] < k[b]; });
            order_.assign(n_nodes_, 0);
            std::vector<uint64_t> sorted(n_nodes_);
            for (size_t i = 0; i < n_nodes_; ++i)
            {
                order_[by_key[i]] = i;
                sorted[i]         = k[by_key[i]];
            }
            for (size_t i = 1; i < n_nodes_; ++i)
                if (sorted[i] == sorted[i - 1]) throw std::runtime_error("SfemWriter: a node owned twice");

            std::filesystem::create_directories(dir_ + "/mesh");
            std::filesystem::create_directories(dir_ + "/out");

            std::vector<float> c(n_nodes_);
            const std::vector<double>* g[3] = {&gx, &gy, &gz};
            const char* xyz                 = "xyz";
            for (int d = 0; d < 3; ++d)
            {
                for (size_t i = 0; i < n_nodes_; ++i)
                    c[order_[i]] = float((*g[d])[i]);
                writeRaw((dir_ + "/mesh/" + xyz[d] + ".float32").c_str(), c.data(), n_nodes_);
            }

            const char* idx = "int32";
            for (int kk = 0; kk < NPE; ++kk)
            {
                std::vector<int32_t> ids(n_elements);
                for (size_t e = 0; e < n_elements; ++e)
                {
                    auto it = std::lower_bound(sorted.begin(), sorted.end(), gconn[kk][e]);
                    if (it == sorted.end() || *it != gconn[kk][e])
                        throw std::runtime_error("SfemWriter: an element node owned by no rank");
                    ids[e] = int32_t(it - sorted.begin());
                }
                const std::string path = dir_ + "/mesh/i" + std::to_string(kk) + "." + idx;
                writeRaw(path.c_str(), ids.data(), n_elements);
            }

            const std::string path = dir_ + "/mesh/meta.yaml";
            FILE* f                = fopen(path.c_str(), "w");
            if (!f) throw std::runtime_error("SfemWriter: cannot open " + path);
            fprintf(f, "# smesh mesh meta file (generated by mars_sfem_writer.hpp)\n");
            fprintf(f, "spatial_dimension: 3\n");
            fprintf(f, "elem_num_nodes: %d\n", NPE);
            fprintf(f, "element_type: %s\n", NPE == 4 ? "TET4" : "HEX8");
            fprintf(f, "geom_map: IsoParametric\n");
            fprintf(f, "n_elements: %ld\n", long(n_elements));
            fprintf(f, "n_nodes: %ld\n", long(n_nodes_));
            fprintf(f, "elements:\n");
            for (int kk = 0; kk < NPE; ++kk)
                fprintf(f, "- i%d: i%d.%s\n", kk, kk, idx);
            fprintf(f, "points:\n");
            for (int d = 0; d < 3; ++d)
                fprintf(f, "- %c: %c.float32\n", xyz[d], xyz[d]);
            fprintf(f, "rpath: true\n");
            const bool failed = ferror(f) != 0;
            if (fclose(f) != 0 || failed) throw std::runtime_error("SfemWriter: metadata write failed");
        }
        domain_       = &domain;
        local_nodes_  = nodeCount;
        mesh_written_ = true;
    }

    template<typename ConnT, std::size_t... Is>
    static void copyConn(std::vector<KeyType>* h, const ConnT& d_conn, size_t n, std::index_sequence<Is...>)
    {
        ((h[Is] = toHost(std::get<Is>(d_conn), n)), ...);
    }
};

} // namespace fem
} // namespace mars
