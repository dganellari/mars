#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <initializer_list>

namespace {
int expected_count = 0;
MPI_Request* expected_requests = nullptr;
MPI_Status* expected_statuses = nullptr;
int returned_error = MPI_SUCCESS;
bool expect_retention = false;
int calls = 0;

void require(bool condition, const char* message) {
    if (!condition) {
        std::fprintf(stderr, "%s\n", message);
        std::abort();
    }
}
}

// Observe the arguments received by MPI without depending on a transport failure.
extern "C" int PMPI_Waitall(int count, MPI_Request* requests, MPI_Status* statuses) {
    ++calls;
    require(count == expected_count && requests == expected_requests, "Request arguments changed");
    if (expect_retention) require(statuses != MPI_STATUSES_IGNORE, "Requested status retention missing");
    else require(statuses == expected_statuses, "Status argument changed without permission");
    if (statuses != MPI_STATUSES_IGNORE)
        for (int i = 0; i < count; ++i) statuses[i].MPI_ERROR = MPI_SUCCESS;
    return returned_error;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    MPI_Request requests[2] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};
    MPI_Status statuses[2];
    expected_requests = requests;
    for (bool retain : {false, true}) {
        setenv("MARS_DIAG_RETAIN_STATUSES", retain ? "1" : "0", 1);
        for (bool ignore : {false, true}) {
            expected_statuses = ignore ? MPI_STATUSES_IGNORE : statuses;
            for (int count : {0, 2}) {
                expected_count = count;
                expect_retention = retain && ignore && count > 0;
                for (int error : {MPI_SUCCESS, MPI_ERR_IN_STATUS}) {
                    returned_error = error;
                    require(MPI_Waitall(count, requests, expected_statuses) == error, "MPI return code changed");
                }
            }
        }
    }
    require(calls == 16, "Missing wait cases");
    std::puts("PASS: original and retained status arguments, empty waits and error returns");
    MPI_Finalize();
}
