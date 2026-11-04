#pragma once

#include <mars_mpi.hpp>

namespace mars::test
{

// Exercise the view collectives through the same staging policy as mesh communication.
inline int check_mpi_view_collectives(MPI_Comm comm)
{
    const int rank = mpi::rank(comm);
    const int ranks = mpi::size(comm);
    int errors = 0;
    int current_errors = 0;

    // Values spanning both bytes expose accidental MPI_CHAR transport of short fields.
    ViewVectorType<short> uniform_global("preflight_uniform_global", 2 * ranks);
    ViewVectorType<short> uniform_local("preflight_uniform_local", 2);
    Kokkos::parallel_for("preflight_scatter_fill", 2 * ranks,
        KOKKOS_LAMBDA(const int i) { uniform_global(i) = rank == 0 ? 513 + i : -1; });
    mpi::scatter(uniform_global, uniform_local, comm);
    Kokkos::parallel_reduce("preflight_scatter_check", 2,
        KOKKOS_LAMBDA(const int i, int& failures) {
            failures += uniform_local(i) != 513 + 2 * rank + i;
        }, current_errors);
    errors += current_errors;

    std::vector<int> counts(ranks);
    for (int peer = 0; peer < ranks; ++peer) counts[peer] = peer + 1;
    const int global_size = ranks * (ranks + 1) / 2;
    const int offset = rank * (rank + 1) / 2;
    ViewVectorType<short> variable_global("preflight_variable_global", global_size);
    ViewVectorType<short> variable_local("preflight_variable_local", rank + 1);
    Kokkos::parallel_for("preflight_scatterv_fill", global_size,
        KOKKOS_LAMBDA(const int i) { variable_global(i) = rank == 0 ? 1025 + i : -1; });
    mpi::scatterv(variable_global, variable_local, counts, comm);
    Kokkos::parallel_reduce("preflight_scatterv_check", rank + 1,
        KOKKOS_LAMBDA(const int i, int& failures) {
            failures += variable_local(i) != 1025 + offset + i;
        }, current_errors);
    errors += current_errors;

    ViewVectorType<short> broadcast_values("preflight_broadcast", 3);
    Kokkos::parallel_for("preflight_broadcast_fill", 3,
        KOKKOS_LAMBDA(const int i) { broadcast_values(i) = rank == 0 ? 769 + i : -1; });
    mpi::broadcast(broadcast_values, comm);
    Kokkos::parallel_reduce("preflight_broadcast_check", 3,
        KOKKOS_LAMBDA(const int i, int& failures) { failures += broadcast_values(i) != 769 + i; }, current_errors);
    errors += current_errors;

    ViewVectorType<int> gathered_scalars("preflight_gathered_scalars", ranks);
    mpi::gather_all_view(rank + 5, gathered_scalars, comm);
    Kokkos::parallel_reduce("preflight_scalar_gather_check", ranks,
        KOKKOS_LAMBDA(const int i, int& failures) { failures += gathered_scalars(i) != i + 5; }, current_errors);
    errors += current_errors;

    Kokkos::parallel_for("preflight_vector_gather_fill", rank + 1,
        KOKKOS_LAMBDA(const int i) { variable_local(i) = 100 * rank + i; });
    ViewVectorType<short> gathered_vectors;
    mpi::gather_all_view(variable_local, gathered_vectors, comm);
    Kokkos::parallel_reduce("preflight_vector_gather_check", ranks,
        KOKKOS_LAMBDA(const int peer, int& failures) {
            const int begin = peer * (peer + 1) / 2;
            for (int i = 0; i <= peer; ++i)
                failures += gathered_vectors(begin + i) != 100 * peer + i;
        }, current_errors);
    errors += current_errors;

    ViewObject<double> reduction_value("preflight_reduction_value");
    Kokkos::parallel_for("preflight_reduction_fill", 1,
        KOKKOS_LAMBDA(const int) { reduction_value(0) = rank + 0.25; });
    auto minimum = mpi::reduce(reduction_value, MPI_MIN, comm);
    auto maximum = mpi::reduce(reduction_value, MPI_MAX, comm);
    Kokkos::parallel_reduce("preflight_reduction_check", 1,
        KOKKOS_LAMBDA(const int, int& failures) {
            failures += minimum(0) != 0.25;
            failures += maximum(0) != ranks - 0.75;
        }, current_errors);
    errors += current_errors;

    ViewObject<short> short_reduction("preflight_short_reduction");
    Kokkos::parallel_for("preflight_short_reduction_fill", 1,
        KOKKOS_LAMBDA(const int) { short_reduction(0) = -1025 + rank; });
    auto short_minimum = mpi::reduce(short_reduction, MPI_MIN, comm);
    auto short_maximum = mpi::reduce(short_reduction, MPI_MAX, comm);
    Kokkos::parallel_reduce("preflight_short_reduction_check", 1,
        KOKKOS_LAMBDA(const int, int& failures) {
            failures += short_minimum(0) != -1025;
            failures += short_maximum(0) != -1026 + ranks;
        }, current_errors);
    errors += current_errors;

    const int source = (rank + ranks - 1) % ranks;
    const int destination = (rank + 1) % ranks;
    std::vector<Integer> send_offsets(ranks + 1, 0), receive_offsets(ranks + 1, 0);
    for (int peer = 0; peer < ranks; ++peer) {
        send_offsets[peer + 1] = send_offsets[peer] + (peer == destination ? 2 : 0);
        receive_offsets[peer + 1] = receive_offsets[peer] + (peer == source ? 2 : 0);
    }
    ViewVectorType<short> outgoing("preflight_short_outgoing", 2), incoming("preflight_short_incoming", 2);
    Kokkos::parallel_for("preflight_short_send_fill", 2,
        KOKKOS_LAMBDA(const int i) { outgoing(i) = 513 + 2 * rank + i; });
    mpi::i_send_recv_view(incoming, receive_offsets.data(), outgoing, send_offsets.data(), comm);
    Kokkos::parallel_reduce("preflight_short_receive_check", 2,
        KOKKOS_LAMBDA(const int i, int& failures) { failures += incoming(i) != 513 + 2 * source + i; }, current_errors);
    return errors + current_errors;
}

} // namespace mars::test
