#ifndef GENERATION_MARS_MESH_DISTRIBUTED_GENERATION_HPP_
#define GENERATION_MARS_MESH_DISTRIBUTED_GENERATION_HPP_
#include "mars_base.hpp"
#include "mars_err.hpp"
#include "mars_sfc_code.hpp"
#include "mars_utils_kokkos.hpp"
#include <cstdint>
#include <limits>
#include <stdexcept>
#ifdef MARS_ENABLE_MPI
#ifdef MARS_ENABLE_KOKKOS
#include "mars_distributed_mesh_kokkos.hpp"
#include "mars_sparse_exchange.hpp"

namespace mars {

    template <Integer Dim, Integer ManifoldDim, Integer Type, class KeyType>
    using DMesh =
        Mesh<Dim, ManifoldDim, DistributedImplementation, NonSimplex<Type, DistributedImplementation>, KeyType>;

    template <Integer Type, class SfcKey>
    constexpr std::uint64_t structured_coordinate_limit() {
        using ValueType = typename SfcKey::ValueType;
        constexpr unsigned bits = Type == ElementType::Quad4 && IsMorton<SfcKey>::value
                                      ? std::numeric_limits<ValueType>::digits / 2
                                      : maxTreeLevel<ValueType>::value;
        return std::uint64_t{1} << bits;
    }


    inline Integer partition_offset(Integer total, int rank, int ranks) {
        return (total / ranks) * rank + std::min<Integer>(rank, total % ranks);
    }

    template <class View, class Key>
    MARS_INLINE_FUNCTION Integer sorted_key_bound(const View& keys, Key key, bool inclusive = false) {
        Integer first = 0, last = keys.extent(0);
        while (first < last) {
            const Integer middle = first + (last - first) / 2;
            if (keys(middle) < key || (inclusive && keys(middle) == key)) first = middle + 1;
            else last = middle;
        }
        return first;
    }

    template <Integer Type, class SfcKey>
    typename SfcKey::ValueType calculate_number_of_elements(Unsigned nx, Unsigned ny, Unsigned nz = 0) {
        using Key = typename SfcKey::ValueType;
        static_assert(Type == ElementType::Quad4 || Type == ElementType::Hex8);
        return Key(nx) * Key(ny) * (Type == ElementType::Hex8 ? Key(nz) : Key(1));
    }

    template <Integer Dim, Integer ManifoldDim, Integer Type, class SfcKey>
    auto get_number_of_elements(const DMesh<Dim, ManifoldDim, Type, SfcKey>& mesh) {
        return calculate_number_of_elements<Type, SfcKey>(mesh.get_XDim(), mesh.get_YDim(), mesh.get_ZDim());
    }

    template <Integer Dim, Integer ManifoldDim, Integer Type, class SfcKey>
    void partition_mesh(DMesh<Dim, ManifoldDim, Type, SfcKey>& mesh) {
        using Key = typename SfcKey::ValueType;
        using Policy = Kokkos::RangePolicy<Kokkos::IndexType<Integer>>;
        const auto& context = mesh.get_context();
        const MPI_Comm comm = context->distributed->mpi_comm();
        if (comm == MPI_COMM_NULL) throw std::invalid_argument("Structured partition requires MPI");
        const int rank = mpi::rank(comm), ranks = mpi::size(comm);
        const Integer total = get_number_of_elements(mesh);
        if (!context->distributed->min(static_cast<int>(total >= ranks)))
            throw std::invalid_argument("Structured partition requires at least one cell per rank");
        const Integer first = partition_offset(total, rank, ranks);
        const Integer count = partition_offset(total, rank + 1, ranks) - first;
        const Integer nx = mesh.get_XDim(), ny = mesh.get_YDim();
        ViewVectorType<Key> keys("partition_local_keys", count);
        Kokkos::parallel_for("partition_generate_local_keys", Policy(0, count), KOKKOS_LAMBDA(Integer i) {
            const Integer linear = first + i;
            const Integer x = linear % nx;
            const Integer y = (linear / nx) % ny;
            if constexpr (Type == ElementType::Quad4) keys(i) = encode_sfc_2D<SfcKey>(x, y);
            else keys(i) = encode_sfc_3D<SfcKey>(x, y, linear / nx / ny);
        });
        Kokkos::sort(keys);

        // Select exact global order statistics using O(ranks) metadata. Cell keys stay on device.
        Key minimum = 0, maximum = 0;
        Kokkos::deep_copy(minimum, Kokkos::subview(keys, 0));
        Kokkos::deep_copy(maximum, Kokkos::subview(keys, count - 1));
        minimum = context->distributed->min(minimum);
        maximum = context->distributed->max(maximum);
        std::vector<Key> low(ranks - 1, minimum), high(ranks - 1, maximum);
        ViewVectorType<Key> queries("partition_queries", ranks - 1);
        auto h_queries = Kokkos::create_mirror_view(queries);
        ViewVectorType<Integer> less_equal("partition_counts", ranks - 1);
        auto h_counts = Kokkos::create_mirror_view(less_equal);
        for (;;) {
            bool unfinished = false;
            for (int p = 0; p < ranks - 1; ++p) {
                unfinished |= low[p] < high[p];
                h_queries(p) = low[p] + (high[p] - low[p]) / 2;
            }
            if (!unfinished) break;
            Kokkos::deep_copy(queries, h_queries);
            Kokkos::parallel_for("partition_order_statistics", ranks - 1, KOKKOS_LAMBDA(int p) {
                less_equal(p) = sorted_key_bound(keys, queries(p), true);
            });
            Kokkos::deep_copy(h_counts, less_equal);
            MPI_OR_THROW(MPI_Allreduce, MPI_IN_PLACE, h_counts.data(), ranks - 1,
                         mpi::mpi_traits<Integer>::mpi_type(), MPI_SUM, comm);
            for (int p = 0; p < ranks - 1; ++p) {
                if (h_counts(p) <= partition_offset(total, p + 1, ranks)) low[p] = h_queries(p) + 1;
                else high[p] = h_queries(p);
            }
        }
        ViewVectorType<Key> partition("global_static_partition", 2 * (ranks + 1));
        auto h_partition = Kokkos::create_mirror_view(partition);
        h_partition(0) = minimum;
        for (int p = 1; p < ranks; ++p) h_partition(2 * p) = low[p - 1];
        h_partition(2 * ranks) = maximum + 1;
        for (int p = 0; p <= ranks; ++p) h_partition(2 * p + 1) = partition_offset(total, p, ranks);
        Kokkos::deep_copy(partition, h_partition);

        ViewVectorType<Integer> offsets("partition_send_offsets", ranks + 1);
        Kokkos::parallel_for("partition_send_offsets", ranks + 1, KOKKOS_LAMBDA(int p) {
            offsets(p) = sorted_key_bound(keys, partition(2 * p));
        });
        auto h_offsets = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, offsets);
        std::vector<mpi::PeerSegment> outgoing = mpi::sparse_segments(h_offsets.data(), ranks, 1, comm);
        std::vector<int> destinations, weights;
        for (const auto& peer : outgoing) { destinations.push_back(peer.rank); weights.push_back(peer.count); }
        const int degree = mpi::checked_mpi_count(outgoing.size(), comm);
        MPI_Comm graph = MPI_COMM_NULL;
        struct GraphLifetime {
            MPI_Comm& comm;
            ~GraphLifetime() { if (comm != MPI_COMM_NULL) MPI_Comm_free(&comm); }
        } graph_lifetime{graph};
        // MPI discovers incoming peers from sparse outgoing edges, without a dense count all-to-all.
        MPI_OR_THROW(MPI_Dist_graph_create, comm, 1, &rank, &degree, destinations.data(), weights.data(),
                     MPI_INFO_NULL, 0, &graph);
        int incoming_count = 0, outgoing_count = 0, weighted = 0;
        MPI_OR_THROW(MPI_Dist_graph_neighbors_count, graph, &incoming_count, &outgoing_count, &weighted);
        std::vector<int> sources(incoming_count), receive_counts(incoming_count);
        MPI_OR_THROW(MPI_Dist_graph_neighbors, graph, incoming_count, sources.data(), receive_counts.data(),
                     outgoing_count, destinations.data(), weights.data());
        std::vector<mpi::PeerSegment> incoming;
        for (int i = 0; i < incoming_count; ++i) incoming.push_back({sources[i], receive_counts[i], 0});
        std::sort(incoming.begin(), incoming.end(), [](const auto& a, const auto& b) { return a.rank < b.rank; });
        std::size_t receive_size = 0;
        for (auto& peer : incoming) { peer.offset = receive_size; receive_size += peer.count; }
        if (receive_size != static_cast<std::size_t>(count))
            mpi::abort_invalid_count(graph, "partition size does not match selected splitters");
        {
            mpi::SparseExchange<Key> exchange(graph, 0, outgoing, incoming, count, receive_size);
            Kokkos::deep_copy(exchange.send_buffer(), keys);
            exchange.start();
            exchange.wait();
            const auto local = exchange.recv_buffer();
            Kokkos::sort(local);
            mesh.set_view_sfc(local);
        }
        MPI_OR_THROW(MPI_Comm_free, &graph);
        mesh.set_proc(rank);
        mesh.set_chunk_size(count);
        mesh.set_view_gp(partition);
        mesh.set_partition_local_counts(count, receive_size);
        const int map_valid = mesh.generate_sfc_to_local_map();
        if (!context->distributed->min(map_valid))
            throw std::runtime_error("Structured mesh SFC map construction failed");
    }

    template <class KeyType, Integer Dim, Integer ManifoldDim, Integer Type, bool Meshless = true>
    void generate_distributed_cube(DMesh<Dim, ManifoldDim, Type, KeyType> &mesh,
                                   const Integer xDim,
                                   const Integer yDim,
                                   const Integer zDim) {
        using namespace Kokkos;

        assert(Dim <= 3);
        assert(ManifoldDim <= Dim);

        const context &context = mesh.get_context();

        static_assert(Type == ElementType::Quad4 || Type == ElementType::Hex8);
        const Integer dimensions[] = {xDim, yDim, Type == ElementType::Hex8 ? zDim : 1};
        const auto coordinate_limit = structured_coordinate_limit<Type, KeyType>();
        Integer cells = 1;
        bool valid = Type == ElementType::Hex8 || zDim == 0;
        for (const Integer dimension : dimensions) {
            if (dimension <= 0 || static_cast<std::uint64_t>(dimension) >= coordinate_limit ||
                cells > std::numeric_limits<Integer>::max() / dimension) {
                valid = false;
                break;
            }
            cells *= dimension;
        }
        valid = valid && cells >= num_ranks(context);
        // Reject on every rank before any partition enters communication or allocates keys.
        if (!context->distributed->min(static_cast<int>(valid))) {
            throw std::invalid_argument(
                "Structured dimensions must fit the SFC key and provide at least one cell per rank.");
        }

        bool same_dimensions = true;
        for (const Integer dimension : dimensions) {
            const auto smallest = context->distributed->min(dimension);
            const auto largest = context->distributed->max(dimension);
            same_dimensions = same_dimensions && smallest == largest;
        }
        if (!same_dimensions) throw std::invalid_argument("Structured dimensions must agree across ranks.");

        int proc_num = rank(context);
        mesh.set_XDim(xDim);
        mesh.set_YDim(yDim);
        mesh.set_ZDim(zDim);

        printf("xDim %ld, yDim %ld, zDim %ld\n", xDim, yDim, zDim);
        // partition the mesh and then generate the points and elements.
        partition_mesh(mesh);

        mesh.template create_ghost_layer<Type>();

        if (!Meshless) {
            Kokkos::Timer timer_gen;

            /* the mesh construct depends on template parameters. */
            auto gen_pts = mesh.template generate_points<Type>();

            auto gen_elm = mesh.template generate_elements<Type>();

            double time_gen = timer_gen.seconds();
            std::cout << "Distributed Mesh Generation took: " << time_gen << " seconds. Process: " << proc_num
                      << std::endl;

            if (!gen_pts || !gen_elm) std::cerr << "Not implemented for other dimensions yet" << std::endl;
        }
    }

}  // namespace mars

#endif
#endif
#endif  // GENERATION_MARS_MESH_DISTRIBUTED_GENERATION_HPP_
