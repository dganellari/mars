#include <mpi.h>
#include <cstdio>
#include <cstring>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm comm;
    MPI_Comm_dup(MPI_COMM_WORLD, &comm);
    MPI_Comm_set_errhandler(comm, MPI_ERRORS_RETURN);
    int ranks = 0, rank = 0, value = 1, total = 0;
    MPI_Comm_size(comm, &ranks);
    MPI_Comm_rank(comm, &rank);
    MPI_Allreduce(&value, &total, 1, MPI_INT, MPI_SUM, comm);
    int incoming = -1;
    MPI_Request requests[2];
    MPI_Irecv(&incoming, 1, MPI_INT, (rank + ranks - 1) % ranks, 1, comm, &requests[0]);
    MPI_Isend(&rank, 1, MPI_INT, (rank + 1) % ranks, 1, comm, &requests[1]);
    int completed = MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);
    if (argc == 2 && std::strcmp(argv[1], "--abort") == 0) MPI_Abort(comm, MPI_ERR_OTHER);
    MPI_Request request = MPI_REQUEST_NULL;
    // A returned error proves the observer preserves the application's error policy.
    int error = MPI_Isend(&value, 1, MPI_INT, ranks, 0, comm, &request);
    int valid = error != MPI_SUCCESS && total == ranks && completed == MPI_SUCCESS &&
                incoming == (rank + ranks - 1) % ranks &&
                requests[0] == MPI_REQUEST_NULL && requests[1] == MPI_REQUEST_NULL;
    int all_valid = 0;
    MPI_Allreduce(&valid, &all_valid, 1, MPI_INT, MPI_MIN, comm);
    std::printf("MARS_MPI_TRACE_PROBE rank=%d ranks=%d valid=%d\n", rank, ranks, all_valid);
    MPI_Comm_free(&comm);
    MPI_Finalize();
    return all_valid ? 0 : 1;
}
