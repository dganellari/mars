#pragma once
// The cell-wise kernels MARSIR generates from mir.dss (PTX files from
// marsir-mlir/test/dss_gpu_pipeline.py), on the tables of an UnstructuredTopology:
//   DssPtx      generated/dss_p7_sm90.ptx: same bits as dss<P> of
//               mars_cellwise_topology.hpp, which it can replace pass for pass;
//   PrecondPtx  generated/precond_p7_sm90.ptx: the Jacobi preconditioner in the
//               four variants of unstructured_precondition_kernel<AZ, ZZ>: same z
//               bits; the dot products as per-element partials;
//   BicgstabPtx generated/bicgstab_p7_sm90.ptx: the vector updates, partial
//               folding and Krylov scalar steps of BiCGStab (mlir_ir.emit_bicgstab);
//   MarsirJacobi  PrecondPtx behind the interface of UnstructuredJacobi, for the
//               Krylov solvers of mars_cellwise_krylov.hpp.
//
// Kernel contract: the arguments of mir.dss in order (values, ghost, elements,
// face_nbr, face_code, edge_of, edge_off, edge_ent, vert_of, vert_off, vert_ent,
// then flags, copies, outs, partials), each a 1-D memref descriptor (allocated
// ptr, aligned ptr, offset, size, stride), integers (index) last; block = 32, one block per listed
// element (measured faster than fewer blocks striding over them). The size of
// `values` is where ghost copies start, so it must be the local copy count
// exactly.

#include "backend/distributed/unstructured/marsir/mars_marsir_ptx_operator.hpp"
#include "backend/distributed/unstructured/solvers/mars_cellwise_krylov.hpp"
#include "backend/distributed/unstructured/solvers/mars_cellwise_topology.hpp"

#include <thrust/device_vector.h>
#include <thrust/sequence.h>

#include <map>
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

    long long values_per_element() const { return n3(); }

protected:
    struct Arg {
        const void* ptr;
        long long size;
    };

    long long n3() const { return (long long)(p_ + 1) * (p_ + 1) * (p_ + 1); }

    // The element list argument: set.ids, or 0 .. local when it is null.
    Arg element_list(const cellwise::ElementList& set, long long local)
    {
        if (set.ids) return {set.ids, set.count};
        if ((long long)all_.size() != local) {
            all_.resize(local);
            thrust::sequence(all_.begin(), all_.end());
        }
        return {thrust::raw_pointer_cast(all_.data()), local};
    }

    // The first 11 arguments: copies (local, ghost), the element list and the tables.
    std::vector<Arg> topology_args(const double* d_in, const double* d_ghost, const cellwise::UnstructuredTopology& T,
                                   const cellwise::ElementList& set, long long values_per_element = 0)
    {
        const long long m = values_per_element ? values_per_element : n3();
        return {{d_in, T.local * m},
                {d_ghost, (T.elements - T.local) * m},
                element_list(set, T.local),
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

    // The kernel and its launch shape: <name>__launch = {threads per block, elements
    // per block} when the compiler published one, else one warp per element.
    CUfunction function(const std::string& name)
    {
        CUfunction fn = nullptr;
        MARS_MARSIR_CU(cuModuleGetFunction(&fn, module_, name.c_str()));
        Shape shape;
        CUdeviceptr global = 0;
        size_t bytes = 0;
        if (cuModuleGetGlobal(&global, &bytes, module_, (name + "__launch").c_str()) == CUDA_SUCCESS &&
            bytes == sizeof(long long[2])) {
            long long v[2];
            MARS_MARSIR_CU(cuMemcpyDtoH(v, global, sizeof(v)));
            shape = {(unsigned)v[0], v[1]};
        }
        shapes_[fn] = shape;
        return fn;
    }

    // `elements`: how many the kernel visits (its element list, or its loop count).
    // Integer (index) and then real (f64) arguments follow the memrefs.
    void launch(CUfunction fn, std::vector<Arg>& args, long long elements, cudaStream_t stream,
                std::vector<long long> integers = {}, std::vector<double> reals = {})
    {
        const Shape shape = shapes_.at(fn);
        const long long blocks = (elements + shape.per_block - 1) / shape.per_block;
        // The driver API takes a pointer to every scalar of every descriptor.
        std::vector<void*> params(args.size() * 5);
        for (size_t i = 0; i < args.size(); ++i) {
            params[5 * i + 0] = &args[i].ptr;
            params[5 * i + 1] = &args[i].ptr;
            params[5 * i + 2] = &zero_;
            params[5 * i + 3] = &args[i].size;
            params[5 * i + 4] = &one_;
        }
        for (long long& v : integers) params.push_back(&v);
        for (double& v : reals) params.push_back(&v);
        MARS_MARSIR_CU(
            cuLaunchKernel(fn, (unsigned)blocks, 1, 1, shape.threads, 1, 1, 0, (CUstream)stream, params.data(), nullptr));
    }

    int p_;

private:
    struct Shape {
        unsigned threads = 32;
        long long per_block = 1;
    };
    CUmodule module_ = nullptr;
    std::map<CUfunction, Shape> shapes_;
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

// The Krylov scalars (the device array of mars_cellwise_krylov.hpp's Reduction) and
// the Step names follow the krylov header.
class BicgstabPtx : public CellwisePtx {
public:
    BicgstabPtx(const std::string& ptx_path, int p) : CellwisePtx(ptx_path, p)
    {
        const std::string q = std::to_string(p);
        update_p_ = function("bicg_update_p_p" + q);
        update_s_ = function("bicg_update_s_p" + q);
        update_xr_ = function("bicg_update_xr_p" + q);
        fold_[0] = function("fold_nv1");
        fold_[1] = function("fold_nv2");
        totals_[0] = function("totals_nv1");
        totals_[1] = function("totals_nv2");
        const char* steps[4] = {"start", "alpha", "omega", "rho"};
        for (int k = 0; k < 4; ++k) {
            finish_[k] = function(std::string("finish_") + steps[k]);
            step_[k] = function(std::string("step_") + steps[k]);
        }
        int device = 0, sms = 0;
        MARS_CELLWISE_CK(cudaGetDevice(&device));
        MARS_CELLWISE_CK(cudaDeviceGetAttribute(&sms, cudaDevAttrMultiProcessorCount, device));
        fold_blocks_ = 2 * sms;   // within Reduction's 4 * 32 * sms doubles for NV = 2
    }

    // p = r + beta (p - omega v) and s = r - alpha v on the `local` elements.
    void update_p(double* p, const double* r, const double* v, const double* sc, long long local, cudaStream_t s)
    {
        std::vector<Arg> a{all(local), {r, local * n3()}, {v, local * n3()}, {p, local * n3()}, {sc, 7}};
        launch(update_p_, a, local, s);
    }
    void update_s(double* sv, const double* r, const double* v, const double* sc, long long local, cudaStream_t s)
    {
        std::vector<Arg> a{all(local), {r, local * n3()}, {v, local * n3()}, {sv, local * n3()}, {sc, 7}};
        launch(update_s_, a, local, s);
    }
    // x += alpha p + omega s, r = s - omega t, and per-element partials of <r, r>, <rh, r>.
    void update_xr(double* x, double* r, const double* p, const double* sv, const double* t, const double* rh,
                   const double* sc, const cellwise::UnstructuredTopology& T, double* partials, cudaStream_t s)
    {
        const long long L = T.local, m = L * n3();
        std::vector<Arg> a{all(L), vec(T.counting), {x, m}, {p, m}, {sv, m}, {t, m}, {rh, m}, {r, m}, {sc, 7},
                           {partials, 2 * L}};
        launch(update_xr_, a, L, s);
    }

    // Per-element partials (count x NV) -> block partials in `blocks` (fold_blocks()
    // x NV, the layout of the hand kernels' store_partial). Returns fold_blocks().
    int fold(int nv, const double* partials, long long count, double* blocks, cudaStream_t s)
    {
        std::vector<Arg> a{{partials, count * nv}, {blocks, (long long)fold_blocks_ * nv}};
        launch(fold_[nv - 1], a, fold_blocks_, s, {count});
        return fold_blocks_;
    }

    // Block partials -> the Krylov step (step: 0 start, 1 alpha, 2 omega, 3 rho), with
    // the totals summed over the ranks of `comm` on the host between two kernels.
    void finish(int step, int nv, const double* blocks, int count, double* sc, double* host, MPI_Comm comm,
                cudaStream_t s)
    {
        int ranks = 1;
        MARS_CELLWISE_MPI(MPI_Comm_size(comm, &ranks));
        std::vector<Arg> a{{blocks, (long long)count * nv}, {sc, 7}};
        if (ranks == 1) {
            launch(finish_[step], a, 1, s, {count});
            return;
        }
        launch(totals_[nv - 1], a, 1, s, {count});
        MARS_CELLWISE_CK(cudaMemcpyAsync(host, sc + cellwise::kT0, nv * sizeof(double), cudaMemcpyDeviceToHost, s));
        MARS_CELLWISE_CK(cudaStreamSynchronize(s));
        MARS_CELLWISE_MPI(MPI_Allreduce(MPI_IN_PLACE, host, nv, MPI_DOUBLE, MPI_SUM, comm));
        MARS_CELLWISE_CK(cudaMemcpyAsync(sc + cellwise::kT0, host, nv * sizeof(double), cudaMemcpyHostToDevice, s));
        std::vector<Arg> b{{sc, 7}};
        launch(step_[step], b, 1, s);
    }

    int fold_blocks() const { return fold_blocks_; }

private:
    Arg all(long long local)
    {
        const cellwise::ElementList set{nullptr, local};
        return element_list(set, local);
    }

    CUfunction update_p_, update_s_, update_xr_, fold_[2], totals_[2], finish_[4], step_[4];
    int fold_blocks_ = 0;
};

// The generated preconditioner as cellwise::UnstructuredJacobi: z = P q, and the
// block partials of <a, z> / <z, z> in `partial`, returning how many blocks wrote
// them (what the Krylov solvers' Reduction expects). The kernels write one partial
// per element into `scratch` (2 * local elements); the generated fold kernel turns
// them into block partials. `halo` is null on one rank.
struct MarsirJacobi {
    PrecondPtx* ptx;
    BicgstabPtx* vec;
    const cellwise::UnstructuredTopology* T;
    const double* diag;
    cellwise::UnstructuredHalo* halo;
    double* scratch;

    // z = P q and the per-element partials in scratch; returns the element count.
    template <bool AZ, bool ZZ>
    long long per_element(const double* q, double* z, const double* a, cudaStream_t stream) const
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
        return done;
    }

    template <bool AZ, bool ZZ>
    int operator()(const double* q, double* z, const double* a, double* partial, cudaStream_t stream) const
    {
        constexpr int NV = AZ + ZZ;
        const long long count = per_element<AZ, ZZ>(q, z, a, stream);
        if constexpr (NV == 0) return 0;
        else return vec->fold(NV, scratch, count, partial, stream);
    }
};

}  // namespace marsir
}  // namespace mars
