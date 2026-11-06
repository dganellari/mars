// CPU-only stand-ins for parsing the real GPU driver's command line.
#define MARS_CELLWISE_MPI(x) x
#define BP_CK(x) x
#define MPI_IN_PLACE nullptr
#define MPI_INT 0
#define MPI_MIN 0
#define MPI_COMM_WORLD 0
struct cudaDeviceProp { int major = 9, minor = 0; };
[[maybe_unused]] void cudaGetDevice(int* device) { *device = 0; }
[[maybe_unused]] void cudaGetDeviceProperties(cudaDeviceProp*, int) {}
[[maybe_unused]] void MPI_Allreduce(void*, int*, int, int, int, int) {}
void MPI_Dims_create(int ranks, int, int* dims)
{
    dims[0] = dims[1] = dims[2] = 1;
    for (int factor = 2; factor <= ranks;) {
        if (ranks % factor) { ++factor; continue; }
        int axis = dims[0] <= dims[1] ? (dims[0] <= dims[2] ? 0 : 2) : (dims[1] <= dims[2] ? 1 : 2);
        dims[axis] *= factor;
        ranks /= factor;
    }
}
