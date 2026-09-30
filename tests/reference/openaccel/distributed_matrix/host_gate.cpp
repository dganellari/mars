// CPU-MPI gate: the adapter's shared functors on host buffers. The CUDA gate runs the
// same run_gates() on device buffers and adds the real Hypre solves.
#include "gate_common.hpp"
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>
using namespace dmatrix_gate;

// Same accessors as mars::fem::SparseMatrix, host storage.
struct HostCsr {
    std::vector<int> offsets, columns; std::vector<double> values;
    void allocate(int rows,int,int nnz) { offsets.assign(std::size_t(rows)+1,0); columns.assign(std::size_t(nnz),0); values.assign(std::size_t(nnz),0); }
    int* rowOffsetsPtr() { return offsets.data(); }
    const int* rowOffsetsPtr() const { return offsets.data(); }
    int* colIndicesPtr() { return columns.data(); }
    const int* colIndicesPtr() const { return columns.data(); }
    double* valuesPtr() { return values.data(); }
    const double* valuesPtr() const { return values.data(); }
};

int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);
    Report report; MPI_Comm_rank(MPI_COMM_WORLD,&report.rank);
    try {
        long long product=0; const long long int_max=std::numeric_limits<int>::max();
        report.result("checked products reject overflow and negative factors",
            checked_product(3,715827882,int_max,product) && product==2147483646 && !checked_product(3,715827883,int_max,product)
            && !checked_product(-1,2,int_max,product) && !checked_product(1LL<<32,1LL<<32,std::numeric_limits<long long>::max(),product)
            && checked_product(0,std::numeric_limits<long long>::max(),int_max,product) && product==0);
        report.result("fault names decode",describe(0)=="none" && describe(capacity|nonfinite)=="capacity; nonfinite owned value or RHS");
        run_gates<1,HostCsr,std::int64_t>(MPI_COMM_WORLD,report);
        run_gates<3,HostCsr,std::int64_t>(MPI_COMM_WORLD,report);
    } catch (const std::exception& e) {
        std::cerr<<"rank "<<report.rank<<" unexpected exception: "<<e.what()<<std::endl; ++report.failures;
    }
    int failures=0; MPI_Allreduce(&report.failures,&failures,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    if (report.rank==0) std::cout<<(failures?"FAIL: ":"PASS: ")<<report.passes<<" passed, "<<report.failures<<" failed, "<<report.skips<<" skipped"<<std::endl;
    MPI_Finalize();
    return failures?1:0;
}
