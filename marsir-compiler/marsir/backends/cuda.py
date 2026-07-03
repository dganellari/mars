"""CUDA backend: emit the device apply kernel from the IR.

Structurally a twin of mars::fem::ho_cvfem_apply_kernel (PerPoint, GStore=double
-- the validated path): same __constant__ reference operators, same ho_idx
layout, same one-thread-per-tangential-slot 4-step sweep, same dynamic-smem
carveout launcher. The only operator-specific part is step 2, where the authored
`flux` is inlined and any tangential contraction the flux does not use is pruned.

The generated banks/symbols are namespaced (mars::marsir, marsir_c_*) so the generated
kernel can be compiled alongside the hand kernel in one TU for the parity gate.
It consumes an externally built d_G (geometry is library infrastructure, not part
of the authored operator) -- e.g. mars::fem::ho_cvfem_metric_perpoint_launch.
"""

from ..ir import ElementApply, jacobian_action_flux
from .common import (LAUNCH_DEFAULT, NAMES_CUDA, NAMES_CUDA_JVP,
                     render_flux, AUTOGEN_BANNER)


def _jvp_block(ea: ElementApply) -> str:
    """The Jacobian-action (JVP) device kernel + launcher: y = J(u) . du.

    Device twin of host_cpp.emit_jvp -- two gathers (u for the partials, du for
    the perturbation contractions), the JVP flux at step 2, then the same
    integrate/scatter skeleton. Emitted into the SAME header as the primal so
    they share the one set of __constant__ banks."""
    o = ea.op
    jvp, _ = jacobian_action_flux(o)
    fv = jvp.free_vars()
    K2 = f"{o.name}_jvp_apply_kernel"
    launch2 = f"{o.name}_jvp_apply_launch"
    flux = jvp.render(NAMES_CUDA_JVP)

    p_deriv = "deriv" in fv
    p_dt1, p_dt2 = "dt1" in fv, "dt2" in fv
    p_interp = p_dt1 or p_dt2
    p_u = p_interp or p_deriv
    d_deriv = "deriv_dot" in fv
    d_dt1, d_dt2 = "dt1_dot" in fv, "dt2_dot" in fv
    d_interp = d_dt1 or d_dt2
    metric = bool(fv & {"g0", "g1", "g2"})

    L = [
        "// Jacobian-action (JVP) kernel: y = J(u) . du, from symbolic form",
        "// differentiation of the flux. Two gathers (u for partials, du for the",
        "// perturbation contractions); same integrate/scatter skeleton as primal.",
        "template<typename RealType, int P, int BlockSize, int ElemsPerBlock>",
        "__global__ void __launch_bounds__(BlockSize)",
        f"{K2}(const RealType* __restrict__ d_u,",
        "    const RealType* __restrict__ d_du,",
        "    RealType* __restrict__ d_y,",
        "    const int* __restrict__ d_elemDof,",
        "    const double* __restrict__ d_G,",
        "    size_t numElements,",
        "    const int* __restrict__ d_elemList = nullptr,",
        "    size_t count = 0)",
        "{",
        '    static_assert(std::is_same<RealType, double>::value, "double only");',
        "    constexpr int N = P + 1, NN = N * N, N3 = NN * N, E = ElemsPerBlock;",
        "    constexpr int threadsPerElem = BlockSize / E;",
        "    constexpr int NP = N, NNP = N * NP;",
        "    extern __shared__ char marsir_jvp_smem_raw[];",
        "    RealType* u_sh  = reinterpret_cast<RealType*>(marsir_jvp_smem_raw);",
        "    RealType* du_sh = u_sh + (size_t)E * N3;",
        "    RealType* y_sh  = du_sh + (size_t)E * N3;",
        "    RealType* face  = y_sh + (size_t)E * N3;   // 3 buffers/elem: fA,fB,fC",
        "    const int t = threadIdx.x;",
        "    const int localElem = t / threadsPerElem;",
        "    const int laneInElem = t % threadsPerElem;",
        "    const size_t slot_e = (size_t)blockIdx.x * E + localElem;",
        "    const size_t e = (d_elemList != nullptr) ? (slot_e < count ? (size_t)d_elemList[slot_e] : numElements) : slot_e;",
        "    RealType* my_u  = u_sh  + localElem * N3;",
        "    RealType* my_du = du_sh + localElem * N3;",
        "    RealType* my_y  = y_sh  + localElem * N3;",
        "    RealType* fA = face + localElem * 3 * NNP;",
        "    RealType* fB = fA + NNP;",
        "    RealType* fC = fB + NNP;",
        "    const bool valid = (e < numElements);",
        "    const int* edof = valid ? (d_elemDof + e * N3) : nullptr;",
        "    for (int l = laneInElem; l < N3; l += threadsPerElem) {",
        "        my_u[l]  = (valid && edof[l] >= 0) ? d_u[edof[l]]  : RealType(0);",
        "        my_du[l] = (valid && edof[l] >= 0) ? d_du[edof[l]] : RealType(0);",
        "        my_y[l]  = RealType(0);",
        "    }",
        "    __syncthreads();",
        "    const double* Gbase = d_G + (valid ? e : 0) * (size_t)(3 * P * NN * 3);",
        "    for (int dir = 0; dir < 3; ++dir)",
        "        for (int l = 0; l < P; ++l) {",
        "            constexpr int kIters = (NN + threadsPerElem - 1) / threadsPerElem;",
    ]
    if p_deriv:
        L.append("            RealType deriv_cache[kIters];")
    if d_deriv:
        L.append("            RealType deriv_dot_cache[kIters];")
    L += [
        "            int it = 0;",
        "            for (int slot = laneInElem; slot < NN; slot += threadsPerElem, ++it) {",
        "                const int s = slot / N, r = slot % N;",
    ]
    accs = []
    if p_interp:
        accs.append("bi = 0")
    if p_deriv:
        accs.append("di = 0")
    if d_interp:
        accs.append("bid = 0")
    if d_deriv:
        accs.append("did = 0")
    if accs:
        L.append("                RealType " + ", ".join(accs) + ";")
    L += ["                #pragma unroll", "                for (int q = 0; q < N; ++q) {"]
    if p_u:
        L.append("                    RealType uq  = my_u[marsir_idx<N, NN>(dir, q, s, r)];")
    if d_interp or d_deriv:
        L.append("                    RealType duq = my_du[marsir_idx<N, NN>(dir, q, s, r)];")
    if p_interp:
        L.append("                    bi  += marsir_c_Btil[l * N + q] * uq;")
    if p_deriv:
        L.append("                    di  += marsir_c_Dtil[l * N + q] * uq;")
    if d_interp:
        L.append("                    bid += marsir_c_Btil[l * N + q] * duq;")
    if d_deriv:
        L.append("                    did += marsir_c_Dtil[l * N + q] * duq;")
    L.append("                }")
    if p_interp:
        L.append("                fA[s * NP + r] = bi;")
    if d_interp:
        L.append("                fB[s * NP + r] = bid;")
    if p_deriv:
        L.append("                deriv_cache[it] = di;")
    if d_deriv:
        L.append("                deriv_dot_cache[it] = did;")
    L += ["            }", "            __syncthreads();", "            it = 0;",
          "            for (int slot = laneInElem; slot < NN; slot += threadsPerElem, ++it) {",
          "                const int s = slot / N, r = slot % N;"]
    if p_dt2:
        L += ["                RealType dt2 = 0;", "                #pragma unroll",
              "                for (int q = 0; q < N; ++q) dt2 += marsir_c_D[r * N + q] * fA[s * NP + q];"]
    if p_dt1:
        L += ["                RealType dt1 = 0;", "                #pragma unroll",
              "                for (int q = 0; q < N; ++q) dt1 += marsir_c_D[s * N + q] * fA[q * NP + r];"]
    if d_dt2:
        L += ["                RealType dt2_dot = 0;", "                #pragma unroll",
              "                for (int q = 0; q < N; ++q) dt2_dot += marsir_c_D[r * N + q] * fB[s * NP + q];"]
    if d_dt1:
        L += ["                RealType dt1_dot = 0;", "                #pragma unroll",
              "                for (int q = 0; q < N; ++q) dt1_dot += marsir_c_D[s * N + q] * fB[q * NP + r];"]
    if metric:
        L.append("                const double* g = Gbase + (size_t)(((dir * P + l) * N + s) * N + r) * 3;")
        if "g0" in fv:
            L.append("                const RealType g0 = RealType(g[0]);")
        if "g1" in fv:
            L.append("                const RealType g1 = RealType(g[1]);")
        if "g2" in fv:
            L.append("                const RealType g2 = RealType(g[2]);")
    L += [
        f"                fC[s * NP + r] = {flux};",
        "            }",
        "            __syncthreads();",
        "            for (int slot = laneInElem; slot < NN; slot += threadsPerElem) {",
        "                const int s = slot / N, r = slot % N;",
        "                RealType v = 0;",
        "                #pragma unroll",
        "                for (int q = 0; q < N; ++q) v += marsir_c_W[r * N + q] * fC[s * NP + q];",
        "                fA[s * NP + r] = v;",
        "            }",
        "            __syncthreads();",
        "            for (int slot = laneInElem; slot < NN; slot += threadsPerElem) {",
        "                const int s = slot / N, r = slot % N;",
        "                RealType intf = 0;",
        "                #pragma unroll",
        "                for (int q = 0; q < N; ++q) intf += marsir_c_W[s * N + q] * fA[q * NP + r];",
        "                my_y[marsir_idx<N, NN>(dir, l,     s, r)] -= intf;",
        "                my_y[marsir_idx<N, NN>(dir, l + 1, s, r)] += intf;",
        "            }",
        "            __syncthreads();",
        "        }",
        "    if (e < numElements) {",
        "        for (int l = laneInElem; l < N3; l += threadsPerElem) {",
        "            int dof = edof[l];",
        "            if (dof >= 0) atomicAdd(&d_y[dof], my_y[l]);",
        "        }",
        "    }",
        "}",
        "",
        "// JVP launcher: smem [u_sh|du_sh|y_sh|face] = E*(3*N3 + 3*NN) doubles.",
        "template<typename RealType, int P,",
        "         int BlockSize = MarsirLaunchDefault<P>::Block,",
        "         int ElemsPerBlock = MarsirLaunchDefault<P>::Elems>",
        f"inline MARSIR_gpuError_t {launch2}(const RealType* d_u, const RealType* d_du, RealType* d_y,",
        "                             const int* d_elemDof, const double* d_G,",
        "                             size_t numElements, MARSIR_gpuStream_t stream = 0,",
        "                             const int* d_elemList = nullptr, size_t count = 0)",
        "{",
        '    static_assert(BlockSize % ElemsPerBlock == 0, "Block % E");',
        "    constexpr int N = P + 1, NN = N * N, N3 = NN * N;",
        "    constexpr int smemBytes = (int)(((size_t)ElemsPerBlock * 3 * N3 + (size_t)ElemsPerBlock * 3 * NN) * sizeof(RealType));",
        f"    auto kernel = {K2}<RealType, P, BlockSize, ElemsPerBlock>;",
        "    MARSIR_gpuError_t attrErr = MARSIR_gpuFuncSetAttribute(kernel, MARSIR_gpuMaxDynamicSmem, smemBytes);",
        "    if (attrErr != MARSIR_gpuSuccess) return attrErr;",
        "    const size_t launchElems = (d_elemList != nullptr) ? count : numElements;",
        "    const size_t grid = (launchElems + ElemsPerBlock - 1) / ElemsPerBlock;",
        "    if (grid == 0) return MARSIR_gpuSuccess;",
        "    kernel<<<(unsigned)grid, BlockSize, smemBytes, stream>>>(d_u, d_du, d_y, d_elemDof, d_G, numElements, d_elemList, count);",
        "    return MARSIR_gpuGetLastError();",
        "}",
    ]
    return "\n".join(L)


def _launch_table() -> str:
    lines = ["template<int P> struct MarsirLaunchDefault;"]
    for p, (block, elems) in sorted(LAUNCH_DEFAULT.items()):
        lines.append(
            f"template<> struct MarsirLaunchDefault<{p}> "
            f"{{ static constexpr int Block = {block}; static constexpr int Elems = {elems}; }};"
        )
    return "\n".join(lines)


def emit(ea: ElementApply) -> str:
    o = ea.op
    flux = render_flux(ea, NAMES_CUDA)
    banner = AUTOGEN_BANNER.format(backend="cuda", name=o.name, flux=o.flux_src)
    K = f"{o.name}_apply_kernel"
    up = f"{o.name}_upload_operators"
    launch = f"{o.name}_apply_launch"
    jvp_block = _jvp_block(ea)

    # step 1: normal contraction. bi (interp) only when the flux has tangential
    # terms; di (deriv) only when the flux uses it.
    s1 = []
    s1.append("                for (int slot = laneInElem; slot < NN; slot += threadsPerElem, ++it) {")
    s1.append("                    const int s = slot / N, r = slot % N;")
    s1.append("                    RealType bi = 0, di = 0; (void)bi; (void)di;")
    s1.append("                    #pragma unroll")
    s1.append("                    for (int q = 0; q < N; ++q) {")
    s1.append("                        RealType uq = my_u[marsir_idx<N, NN>(dir, q, s, r)];")
    if ea.needs_tangential:
        s1.append("                        bi += marsir_c_Btil[l * N + q] * uq;")
    if ea.needs_deriv:
        s1.append("                        di += marsir_c_Dtil[l * N + q] * uq;")
    s1.append("                    }")
    if ea.needs_tangential:
        s1.append("                    faceA[s * NP + r] = bi;")
    s1.append("                    deriv_cache[it] = di;")
    s1.append("                }")
    step1 = "\n".join(s1)

    # step 2: tangential contractions (pruned) + metric reads (only used comps) + flux.
    s2 = []
    s2.append("                for (int slot = laneInElem; slot < NN; slot += threadsPerElem, ++it) {")
    s2.append("                    const int s = slot / N, r = slot % N;")
    if ea.needs_dt2:
        s2.append("                    RealType dt2 = 0;")
        s2.append("                    #pragma unroll")
        s2.append("                    for (int q = 0; q < N; ++q) dt2 += marsir_c_D[r * N + q] * faceA[s * NP + q];")
    if ea.needs_dt1:
        s2.append("                    RealType dt1 = 0;")
        s2.append("                    #pragma unroll")
        s2.append("                    for (int q = 0; q < N; ++q) dt1 += marsir_c_D[s * N + q] * faceA[q * NP + r];")
    if ea.needs_metric:
        s2.append("                    const double* g = Gbase + (size_t)(((dir * P + l) * N + s) * N + r) * 3;")
        if "g0" in ea.free_vars:
            s2.append("                    const RealType g0 = RealType(g[0]);")
        if "g1" in ea.free_vars:
            s2.append("                    const RealType g1 = RealType(g[1]);")
        if "g2" in ea.free_vars:
            s2.append("                    const RealType g2 = RealType(g[2]);")
    s2.append(f"                    faceB[s * NP + r] = {flux};")
    s2.append("                }")
    step2 = "\n".join(s2)

    return f"""#pragma once
{banner}
// Portable GPU veneer: this ONE generated header targets NVIDIA (nvcc) and AMD
// (hipcc). The device code below is already vendor-neutral (__global__,
// __constant__, __syncthreads, atomicAdd, __launch_bounds__, <<<>>>); only the
// host runtime API differs, mapped here. Under CUDA these resolve to cuda*
// exactly as before.
#if defined(__HIP_PLATFORM_AMD__) || defined(__HIPCC__)
#include <hip/hip_runtime.h>
#define MARSIR_gpuError_t          hipError_t
#define MARSIR_gpuSuccess          hipSuccess
#define MARSIR_gpuStream_t         hipStream_t
#define MARSIR_gpuMemcpyToSymbol   hipMemcpyToSymbol
#define MARSIR_gpuGetLastError     hipGetLastError
#define MARSIR_gpuFuncSetAttribute hipFuncSetAttribute
#define MARSIR_gpuMaxDynamicSmem   hipFuncAttributeMaxDynamicSharedMemorySize
#else
#include <cuda_runtime.h>
#define MARSIR_gpuError_t          cudaError_t
#define MARSIR_gpuSuccess          cudaSuccess
#define MARSIR_gpuStream_t         cudaStream_t
#define MARSIR_gpuMemcpyToSymbol   cudaMemcpyToSymbol
#define MARSIR_gpuGetLastError     cudaGetLastError
#define MARSIR_gpuFuncSetAttribute cudaFuncSetAttribute
#define MARSIR_gpuMaxDynamicSmem   cudaFuncAttributeMaxDynamicSharedMemorySize
#endif
#include <cstddef>
#include <type_traits>

#ifndef MARSIR_MAX_P
#define MARSIR_MAX_P 8
#endif

namespace mars {{
namespace marsir {{

namespace detail {{
constexpr int kMaxN  = MARSIR_MAX_P + 1;
constexpr int kBDlen = MARSIR_MAX_P * kMaxN;   // Btil / Dtil
constexpr int kDWlen = kMaxN * kMaxN;       // D / W
}}

// Reference operators in constant memory (own banks so this coexists with the
// hand kernel's c_* banks in one TU). This header defines __constant__ symbols,
// so include it from EXACTLY ONE translation unit.
__constant__ double marsir_c_Btil[detail::kBDlen];
__constant__ double marsir_c_Dtil[detail::kBDlen];
__constant__ double marsir_c_D[detail::kDWlen];
__constant__ double marsir_c_W[detail::kDWlen];

// Upload host HoCvfemOperators arrays once per P before launching.
inline MARSIR_gpuError_t {up}(int P,
                        const double* h_Btil, const double* h_Dtil,
                        const double* h_D, const double* h_W)
{{
    const int n = P + 1;
    MARSIR_gpuError_t err;
    err = MARSIR_gpuMemcpyToSymbol(marsir_c_Btil, h_Btil, sizeof(double) * (size_t)P * n); if (err) return err;
    err = MARSIR_gpuMemcpyToSymbol(marsir_c_Dtil, h_Dtil, sizeof(double) * (size_t)P * n); if (err) return err;
    err = MARSIR_gpuMemcpyToSymbol(marsir_c_D,    h_D,    sizeof(double) * (size_t)n * n); if (err) return err;
    err = MARSIR_gpuMemcpyToSymbol(marsir_c_W,    h_W,    sizeof(double) * (size_t)n * n); if (err) return err;
    return MARSIR_gpuSuccess;
}}

template<int N, int NN>
__device__ __forceinline__ int marsir_idx(int dir, int nrm, int t1, int t2)
{{
    if (dir == 0) return nrm * NN + t1 * N + t2;
    if (dir == 1) return t1 * NN + nrm * N + t2;
    return t1 * NN + t2 * N + nrm;
}}

// PerPoint apply, GStore=double. d_y MUST be pre-zeroed; d_G built externally
// (e.g. mars::fem::ho_cvfem_metric_perpoint_launch); operators uploaded via {up}.
template<typename RealType, int P, int BlockSize, int ElemsPerBlock>
__global__ void __launch_bounds__(BlockSize)
{K}(const RealType* __restrict__ d_u,
    RealType* __restrict__ d_y,
    const int* __restrict__ d_elemDof,
    const double* __restrict__ d_G,
    size_t numElements,
    const int* __restrict__ d_elemList = nullptr,
    size_t count = 0)
{{
    static_assert(std::is_same<RealType, double>::value,
        "marsir HO-CVFEM apply is validated for RealType=double only");
    constexpr int N   = P + 1;
    constexpr int NN  = N * N;
    constexpr int N3  = NN * N;
    constexpr int E   = ElemsPerBlock;
    constexpr int threadsPerElem = BlockSize / E;
    constexpr int NP  = N;
    constexpr int NNP = N * NP;

    extern __shared__ char marsir_smem_raw[];
    RealType* u_sh = reinterpret_cast<RealType*>(marsir_smem_raw);
    RealType* y_sh = u_sh + (size_t)E * N3;
    RealType* face = y_sh + (size_t)E * N3;

    const int t          = threadIdx.x;
    const int localElem  = t / threadsPerElem;
    const int laneInElem = t % threadsPerElem;
    const size_t slot_e  = (size_t)blockIdx.x * E + localElem;
    const size_t e       = (d_elemList != nullptr)
                             ? (slot_e < count ? (size_t)d_elemList[slot_e] : numElements)
                             : slot_e;

    RealType* my_u  = u_sh + localElem * N3;
    RealType* my_y  = y_sh + localElem * N3;
    RealType* faceA = face + localElem * 2 * NNP;
    RealType* faceB = faceA + NNP;
    const bool valid = (e < numElements);
    const int* edof  = valid ? (d_elemDof + e * N3) : nullptr;

    for (int l = laneInElem; l < N3; l += threadsPerElem) {{
        my_u[l] = (valid && edof[l] >= 0) ? d_u[edof[l]] : RealType(0);
        my_y[l] = RealType(0);
    }}
    __syncthreads();

    const double* Gbase = d_G + (valid ? e : 0) * (size_t)(3 * P * NN * 3);

    for (int dir = 0; dir < 3; ++dir)
        for (int l = 0; l < P; ++l) {{
            constexpr int kDerivIters = (NN + threadsPerElem - 1) / threadsPerElem;
            RealType deriv_cache[kDerivIters];
            int it = 0;
{step1}
            __syncthreads();

            it = 0;
{step2}
            __syncthreads();

            for (int slot = laneInElem; slot < NN; slot += threadsPerElem) {{
                const int s = slot / N, r = slot % N;
                RealType v = 0;
                #pragma unroll
                for (int q = 0; q < N; ++q) v += marsir_c_W[r * N + q] * faceB[s * NP + q];
                faceA[s * NP + r] = v;
            }}
            __syncthreads();

            for (int slot = laneInElem; slot < NN; slot += threadsPerElem) {{
                const int s = slot / N, r = slot % N;
                RealType intf = 0;
                #pragma unroll
                for (int q = 0; q < N; ++q) intf += marsir_c_W[s * N + q] * faceA[q * NP + r];
                my_y[marsir_idx<N, NN>(dir, l,     s, r)] -= intf;
                my_y[marsir_idx<N, NN>(dir, l + 1, s, r)] += intf;
            }}
            __syncthreads();
        }}

    if (e < numElements) {{
        for (int l = laneInElem; l < N3; l += threadsPerElem) {{
            int dof = edof[l];
            if (dof >= 0) atomicAdd(&d_y[dof], my_y[l]);
        }}
    }}
}}

{_launch_table()}

// Launcher mirrors ho_cvfem_apply_launch_impl: dynamic smem [u_sh|y_sh|face] =
// E*2*(N3+NN) doubles, Hopper carveout opt-in, grid = ceil(elems/E).
template<typename RealType, int P,
         int BlockSize = MarsirLaunchDefault<P>::Block,
         int ElemsPerBlock = MarsirLaunchDefault<P>::Elems>
inline MARSIR_gpuError_t {launch}(const RealType* d_u, RealType* d_y,
                            const int* d_elemDof, const double* d_G,
                            size_t numElements, MARSIR_gpuStream_t stream = 0,
                            const int* d_elemList = nullptr, size_t count = 0)
{{
    static_assert(BlockSize % ElemsPerBlock == 0,
        "BlockSize must be a multiple of ElemsPerBlock");
    constexpr int N = P + 1, NN = N * N, N3 = NN * N;
    constexpr int smemBytes =
        (int)(((size_t)ElemsPerBlock * 2 * N3 + (size_t)ElemsPerBlock * 2 * NN) * sizeof(RealType));

    auto kernel = {K}<RealType, P, BlockSize, ElemsPerBlock>;
    MARSIR_gpuError_t attrErr = MARSIR_gpuFuncSetAttribute(
        kernel, MARSIR_gpuMaxDynamicSmem, smemBytes);
    if (attrErr != MARSIR_gpuSuccess) return attrErr;

    const size_t launchElems = (d_elemList != nullptr) ? count : numElements;
    const size_t grid = (launchElems + ElemsPerBlock - 1) / ElemsPerBlock;
    if (grid == 0) return MARSIR_gpuSuccess;
    kernel<<<(unsigned)grid, BlockSize, smemBytes, stream>>>(
        d_u, d_y, d_elemDof, d_G, numElements, d_elemList, count);
    return MARSIR_gpuGetLastError();
}}

{jvp_block}

}}  // namespace marsir
}}  // namespace mars
"""
