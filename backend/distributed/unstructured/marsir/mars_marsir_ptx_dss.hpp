#pragma once
// The cell-wise kernels MARSIR generates from mir.dss (PTX files from
// marsir-mlir/test/dss_gpu_pipeline.py), on the tables of an UnstructuredTopology:
//   DssPtx      generated/dss_p7_sm90.ptx: same bits as dss<P> of
//               mars_cellwise_topology.hpp, which it can replace pass for pass;
//   PrecondPtx  generated/precond_p7_sm90.ptx: the Jacobi preconditioner in the
//               four variants of unstructured_precondition_kernel<AZ, ZZ>: same z
//               bits; the dot products as per-element partials;
//   MarsirJacobi  PrecondPtx behind the interface of UnstructuredJacobi, for the
//               Krylov solvers of mars_cellwise_krylov.hpp.
//
// Kernel contract: the arguments of mir.dss in order (values, ghost, elements,
// face_nbr, face_code, edge_of, edge_off, edge_ent, vert_of, vert_off, vert_ent,
// then flags, copies, outs, partials), each a 1-D memref descriptor (allocated
// ptr, aligned ptr, offset, size, stride); block = 32, one block per listed
// element (measured faster than fewer blocks striding over them). The size of
// `values` is where ghost copies start, so it must be the local copy count
// exactly.

#include "backend/distributed/unstructured/marsir/mars_marsir_ptx_operator.hpp"
#include "backend/distributed/unstructured/solvers/mars_cellwise_krylov.hpp"
#include "backend/distributed/unstructured/solvers/mars_cellwise_topology.hpp"

#include <thrust/device_vector.h>
#include <thrust/sequence.h>

#include <string>
#include <vector>

namespace mars {
namespace marsir {

// A loaded mir.dss kernel and the descriptor packing every such kernel shares.
class CellwisePtx {
public:
    CellwisePtx(const std::string& ptx_path, int p) : p_(p) { module_ = load_ptx(ptx_path); }
    ~CellwisePtx() { cuModuleUnload(module_); }
    CellwisePtx(const CellwisePtx&) = delete;
    CellwisePtx& operator=(const CellwisePtx&) = delete;

protected:
    struct Arg {
        const void* ptr;
        long long size;
    };

    long long n3() const { return (long long)(p_ + 1) * (p_ + 1) * (p_ + 1); }

    // The first 11 arguments: copies (local, ghost), the element list and the tables.
    std::vector<Arg> topology_args(const double* d_in, const double* d_ghost, const cellwise::UnstructuredTopology& T,
                                   const cellwise::ElementList& set)
    {
        const int* ids = set.ids;
        if (!ids) {   // the kernel always reads a list; the whole-mesh pass reads 0 .. local
            if ((long long)all_.size() != T.local) {
                all_.resize(T.local);
                thrust::sequence(all_.begin(), all_.end());
            }
            ids = thrust::raw_pointer_cast(all_.data());
        }
        return {{d_in, T.local * n3()},
                {d_ghost, (T.elements - T.local) * n3()},
                {ids, set.count},
                vec(T.face_nbr),
                vec(T.face_code),
                vec(T.edge_of),
                vec(T.edge_off),
                vec(T.edge_ent),
                vec(T.vert_of),
                vec(T.vert_off),
                vec(T.vert_ent)};
    }

    template <typename T>
    static Arg vec(const thrust::device_vector<T>& v)
    {
        return {thrust::raw_pointer_cast(v.data()), (long long)v.size()};
    }

    CUfunction function(const std::string& name) const
    {
        CUfunction fn = nullptr;
        MARS_MARSIR_CU(cuModuleGetFunction(&fn, module_, name.c_str()));
        return fn;
    }

    void launch(CUfunction fn, std::vector<Arg>& args, long long blocks, cudaStream_t stream)
    {
        // The driver API takes a pointer to every scalar of every descriptor.
        std::vector<void*> params(args.size() * 5);
        for (size_t i = 0; i < args.size(); ++i) {
            params[5 * i + 0] = &args[i].ptr;
            params[5 * i + 1] = &args[i].ptr;
            params[5 * i + 2] = &zero_;
            params[5 * i + 3] = &args[i].size;
            params[5 * i + 4] = &one_;
        }
        MARS_MARSIR_CU(
            cuLaunchKernel(fn, (unsigned)blocks, 1, 1, 32, 1, 1, 0, (CUstream)stream, params.data(), nullptr));
    }

    int p_;

private:
    CUmodule module_ = nullptr;
    long long zero_ = 0, one_ = 1;
    thrust::device_vector<int> all_;
};

class DssPtx : public CellwisePtx {
public:
    DssPtx(const std::string& ptx_path, int p) : CellwisePtx(ptx_path, p), fn_(function("dss_p" + std::to_string(p)))
    {
    }

    // out = DSS(in) on the elements of `set` (all local elements when set.ids is
    // null); `ghost` holds the ghost elements' copies (null on one rank).
    void run(const double* d_in, const double* d_ghost, double* d_out, const cellwise::UnstructuredTopology& T,
             const cellwise::ElementList& set, cudaStream_t stream = 0)
    {
        if (set.count == 0) return;
        std::vector<Arg> args = topology_args(d_in, d_ghost, T, set);
        args.push_back({d_out, T.local * n3()});
        launch(fn_, args, set.count, stream);
    }

private:
    CUfunction fn_;
};

class PrecondPtx : public CellwisePtx {
public:
    PrecondPtx(const std::string& ptx_path, int p) : CellwisePtx(ptx_path, p)
    {
        const std::string base = "precond_p" + std::to_string(p);
        fn_[0][0] = function(base);
        fn_[0][1] = function(base + "_zz");
        fn_[1][0] = function(base + "_az");
        fn_[1][1] = function(base + "_az_zz");
    }

    // z = P q on the elements of `set`; with AZ / ZZ also the w-th listed element's
    // share of <a, z> / <z, z> (each node counted once) at partials[NV w ...].
    template <bool AZ, bool ZZ>
    void run(const double* d_q, const double* d_ghost, double* d_z, const double* d_diag, const double* d_a,
             double* d_partials, const cellwise::UnstructuredTopology& T, const cellwise::ElementList& set,
             cudaStream_t stream = 0)
    {
        if (set.count == 0) return;
        constexpr int NV = AZ + ZZ;
        std::vector<Arg> args = topology_args(d_q, d_ghost, T, set);
        args.push_back(vec(T.counting));
        args.push_back(vec(T.edge_dirichlet));
        args.push_back(vec(T.vert_dirichlet));
        args.push_back({d_diag, T.local * n3()});
        if (AZ) args.push_back({d_a, T.local * n3()});
        args.push_back({d_z, T.local * n3()});
        if (NV > 0) args.push_back({d_partials, NV * set.count});
        launch(fn_[AZ][ZZ], args, set.count, stream);
    }

private:
    CUfunction fn_[2][2];
};

// Block partials (NV values each) of the per-element partials src[0 .. count):
// block b sums a contiguous range in a fixed order, so the result is reproducible.
template <int NV>
__global__ void __launch_bounds__(cellwise::kThreads)
sum_partials_kernel(const double* __restrict__ src, long long count, double* __restrict__ dst)
{
    double v[NV] = {};
    const long long per = (count + gridDim.x - 1) / gridDim.x, lo = blockIdx.x * per;
    const long long hi = lo + per < count ? lo + per : count;
    for (long long i = lo + threadIdx.x; i < hi; i += blockDim.x)
#pragma unroll
        for (int j = 0; j < NV; ++j) v[j] += src[i * NV + j];
    cellwise::store_partial<NV>(v, dst);
}

// The generated preconditioner as cellwise::UnstructuredJacobi: z = P q, and the
// block partials of <a, z> / <z, z> in `partial`, returning how many blocks wrote
// them (what the Krylov solvers' Reduction expects). The kernels write one partial
// per element into `scratch` (2 * local elements); sum_partials_kernel folds them
// into a few hundred blocks. `halo` is null on one rank.
struct MarsirJacobi {
    PrecondPtx* ptx;
    const cellwise::UnstructuredTopology* T;
    const double* diag;
    cellwise::UnstructuredHalo* halo;
    double* scratch;

    template <bool AZ, bool ZZ>
    int operator()(const double* q, double* z, const double* a, double* partial, cudaStream_t stream) const
    {
        constexpr int NV = AZ + ZZ;
        long long done = 0;   // elements written so far: the interior pass, then the rest
        cellwise::run_unstructured_pass(halo, T->local, q, nullptr, stream,
                                        [&](const cellwise::ElementList& set, const double* ghost, int) {
                                            ptx->run<AZ, ZZ>(q, ghost, z, diag, a, scratch + done * NV, *T, set,
                                                             stream);
                                            done += set.count;
                                            return 0;
                                        });
        if constexpr (NV == 0) {
            return 0;
        } else {
            // Within Reduction's 4 * 32 * sms doubles for NV = 2; one device per rank.
            static const int blocks = [] {
                int device = 0, sms = 0;
                MARS_CELLWISE_CK(cudaGetDevice(&device));
                MARS_CELLWISE_CK(cudaDeviceGetAttribute(&sms, cudaDevAttrMultiProcessorCount, device));
                return 2 * sms;
            }();
            sum_partials_kernel<NV><<<blocks, cellwise::kThreads, 0, stream>>>(scratch, done, partial);
            MARS_CELLWISE_CK(cudaGetLastError());
            return blocks;
        }
    }
};

}  // namespace marsir
}  // namespace mars
