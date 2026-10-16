#pragma once
// The cell-wise DSS kernel MARSIR generates from mir.dss (a PTX file such as
// marsir-mlir/generated/dss_p7_sm90.ptx, from test/dss_gpu_pipeline.py), on the
// tables of an UnstructuredTopology. Same results, bit for bit, as dss<P> of
// mars_cellwise_topology.hpp, which it can replace pass for pass.
//
// Kernel contract: dss_p{p}(values, ghost, elements, face_nbr, face_code, edge_of,
// edge_off, edge_ent, vert_of, vert_off, vert_ent, out), each a 1-D memref
// descriptor (allocated ptr, aligned ptr, offset, size, stride); grid = listed
// elements, block = 32. The size of `values` is where ghost copies start, so it
// must be the local copy count exactly.

#include "backend/distributed/unstructured/marsir/mars_marsir_ptx_operator.hpp"
#include "backend/distributed/unstructured/solvers/mars_cellwise_topology.hpp"

#include <thrust/device_vector.h>
#include <thrust/sequence.h>

#include <string>

namespace mars {
namespace marsir {

class DssPtx {
public:
    DssPtx(const std::string& ptx_path, int p) : p_(p)
    {
        module_ = load_ptx(ptx_path);
        MARS_MARSIR_CU(cuModuleGetFunction(&fn_, module_, ("dss_p" + std::to_string(p)).c_str()));
    }
    ~DssPtx() { cuModuleUnload(module_); }
    DssPtx(const DssPtx&) = delete;
    DssPtx& operator=(const DssPtx&) = delete;

    // out = DSS(in) on the elements of `set` (all local elements when set.ids is
    // null); `ghost` holds the ghost elements' copies (null on one rank).
    void run(const double* d_in, const double* d_ghost, double* d_out, const cellwise::UnstructuredTopology& T,
             const cellwise::ElementList& set, cudaStream_t stream = 0)
    {
        if (set.count == 0) return;
        const long long n3 = (long long)(p_ + 1) * (p_ + 1) * (p_ + 1);
        const int* ids = set.ids;
        if (!ids) {   // the kernel always reads a list; the whole-mesh pass reads 0 .. local
            if ((long long)all_.size() != T.local) {
                all_.resize(T.local);
                thrust::sequence(all_.begin(), all_.end());
            }
            ids = thrust::raw_pointer_cast(all_.data());
        }
        Arg a[12] = {{d_in, T.local * n3},
                     {d_ghost, (T.elements - T.local) * n3},
                     {ids, set.count},
                     {raw(T.face_nbr), (long long)T.face_nbr.size()},
                     {raw(T.face_code), (long long)T.face_code.size()},
                     {raw(T.edge_of), (long long)T.edge_of.size()},
                     {raw(T.edge_off), (long long)T.edge_off.size()},
                     {raw(T.edge_ent), (long long)T.edge_ent.size()},
                     {raw(T.vert_of), (long long)T.vert_of.size()},
                     {raw(T.vert_off), (long long)T.vert_off.size()},
                     {raw(T.vert_ent), (long long)T.vert_ent.size()},
                     {d_out, T.local * n3}};
        // The driver API takes a pointer to every scalar of every descriptor.
        void* params[12 * 5];
        for (int i = 0; i < 12; ++i) {
            params[5 * i + 0] = &a[i].ptr;
            params[5 * i + 1] = &a[i].ptr;
            params[5 * i + 2] = &zero_;
            params[5 * i + 3] = &a[i].size;
            params[5 * i + 4] = &one_;
        }
        MARS_MARSIR_CU(cuLaunchKernel(fn_, (unsigned)set.count, 1, 1, 32, 1, 1, 0, (CUstream)stream, params, nullptr));
    }

private:
    struct Arg {
        const void* ptr;
        long long size;
    };
    template <typename T>
    static const void* raw(const thrust::device_vector<T>& v)
    {
        return thrust::raw_pointer_cast(v.data());
    }

    CUmodule module_ = nullptr;
    CUfunction fn_ = nullptr;
    int p_;
    long long zero_ = 0, one_ = 1;
    thrust::device_vector<int> all_;
};

}  // namespace marsir
}  // namespace mars
