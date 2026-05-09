# MARS CVFEM Kernel Family — Deep Reference

Citation-heavy reference for `backend/distributed/unstructured/fem/`. For
the dispatcher and how it plugs into ElementDomain see
[02_mars_overview.md](02_mars_overview.md).

## 1. Mathematical Foundation

**Control-Volume Finite-Element Method (CVFEM)** — conservative
unstructured discretization combining FE shape functions with flux-based
control volumes (Patankar 1980, Shashkov–Steinberg).

### Sub-Control-Surface (SCS) integration

A hex8 element is split into **12 SCS** — faces connecting element center
to edge midpoints:
- 4 bottom (SCS 0–3, z=−0.5): edges 0-1, 1-2, 2-3, 3-0
- 4 top (SCS 4–7, z=+0.5): edges 4-5, 5-6, 6-7, 7-4
- 4 vertical (SCS 8–11, z=0): edges 0-4, 1-5, 2-6, 3-7

Each SCS = quad with 4 nodes from the 27-point hex subdivision (8 corners
+ 12 edge mids + 6 face centers + 1 center). Integration points are at
face locations (±0.5), aligning element-center flux surfaces to hex faces.

**LRSCV array** — `mars_cvfem_hex_kernel.hpp:88–92`:
```c
__device__ __constant__ int hexLRSCV[24] = {
    0,1, 1,2, 2,3, 0,3,        // bottom 4
    4,5, 5,6, 6,7, 4,7,        // top 4
    0,4, 1,5, 2,6, 3,7         // vertical 4
};
```
nodeL, nodeR define opposing control volumes separated by each SCS face.

### Shape functions

`hexShapeFcn[12][8]` — `mars_cvfem_hex_kernel.hpp:60–85`. At each SCS,
**only 2 of 8 nodes** have non-zero shape functions (0.5 each). Parametric
derivatives at IPs in `hexDerivConst[12][8][3]` (L96–121).

### Flux formulation

Per SCS integration point `ip`:

**Advection:**
```
phi_upwind = upwind(mdot, phi[nodeL], phi[nodeR])
dcorr      = beta_upwind * grad_phi · (coords_ip - coords_upwind)
adv_flux   = mdot * (phi_upwind + dcorr)
```

**Diffusion:**
```
diff_coeff[n] = -gamma_ip * (dndx[n] · areaVec)    ∀ n ∈ [0,8]
```

### Sparsity patterns
- **Full** (27 NNZ/row): all 8×8 element-local pairs → 64 entries/element
- **Graph** (7 NNZ/row, diagonal lumping): only 12 SCS edges (nodeL,
  nodeR) contribute. Missing 45 entries lumped to diagonal — preserves row
  sums (conservation). See `mars_sparsity_builder.hpp:104–127`.

## 2. Reference Element & Shape Functions

### Jacobian — `mars_cvfem_hex_kernel.hpp:124–138`
```cuda
__device__ inline void invert3x3(double J[3][3], double invJ[3][3]) {
    double det = J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1])
               - J[0][1]*(J[1][0]*J[2][2] - J[1][2]*J[2][0])
               + J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
    double invDet = 1.0 / det;
    invJ[0][0] = (J[1][1]*J[2][2] - J[1][2]*J[2][1]) * invDet;
    // ... 8 more lines
}
```

### Physical shape derivatives — `mars_cvfem_hex_kernel.hpp:141–170`
```cuda
__device__ inline void computeShapeDerivatives(
    int ip, const double coords[8][3], double dndx[8][3])
{
    double J[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    for (int node = 0; node < 8; ++node)
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                J[i][j] += hexDerivConst[ip][node][j] * coords[node][i];
    double invJ[3][3];
    invert3x3(J, invJ);
    for (int node = 0; node < 8; ++node)
        for (int i = 0; i < 3; ++i) {
            dndx[node][i] = 0.0;
            for (int j = 0; j < 3; ++j)
                dndx[node][i] += invJ[i][j] * hexDerivConst[ip][node][j];
        }
}
```
Implements ∇_x N = J^{-T} · ∇_ξ N.

### Area vectors — `mars_cvfem_hex_kernel.hpp:265–286`
27-point hex subdivision (L173–212). Each SCS quad's area computed via
**triangulation from quad centroid** (L215–256):
```cuda
__device__ inline void quad_area_by_triangulation(
    const double areacoords[4][3], double areaVec[3])
{
    double xmid[3] = { 0.25*sum_x, 0.25*sum_y, 0.25*sum_z };
    for (int itri = 0; itri < 4; ++itri) {
        int t = (itri + 1) % 4;
        double r2[3] = areacoords[t] - xmid;
        areaVec += 0.5 * cross(r1, r2);
        r1 = r2;
    }
}
```

## 3. NodeData Layout (AoS)

`mars_cvfem_node_data.hpp:14–24`
```cuda
struct alignas(32) NodeData {
    double x, y, z;           // 24 B
    double phi, gamma, beta;  // 24 B
    double gx, gy, gz;        // 24 B
};  // 72 B = 3 L1 sectors (32 B each)
```
Rationale: 32 threads × scattered nodes hit
- SoA (9 arrays): 9 × 32 = 288 sectors
- AoS (1 struct): ~96 sectors
**3× L2 traffic reduction**.

Packing kernel (L39–63): one thread per node, copies 9 SoA arrays into AoS.

## 4. Sparse Matrix CSR

`mars_cvfem_hex_kernel.hpp:9–38`
```cuda
template<typename RealType>
struct CSRMatrix {
    int* rowPtr;          // [numRows+1]
    int* colInd;          // [nnz], sorted within row
    RealType* values;     // [nnz]
    int* diagPtr;         // [numRows] direct diagonal index
    int numRows, nnz;

    __device__ void addValue(int row, int col, RealType val) {
        // Linear search (baseline only) — slow
        for (int i = rowPtr[row]; i < rowPtr[row+1]; ++i)
            if (colInd[i] == col) { atomicAdd(&values[i], val); return; }
    }
    __device__ __forceinline__ RealType* diag(int row) {
        return &values[diagPtr[row]];
    }
};
```
Columns sorted within each row → binary search (`tensor.hpp:231–240`).

## 5. Sparsity Builder — `mars_sparsity_builder.hpp`

### Graph sparsity — L130–273
1. **Edge-list kernel** (L77–127): each thread emits 12×2=24 edge tuples per
   element, plus invalid for ghosts
2. **Thrust pipeline** (L149–213):
   - Allocate edge row/col vectors
   - Append diag entries (`thrust::sequence`)
   - Remove `row < 0` (ghosts)
   - `thrust::sort` by (row, col)
   - `thrust::unique` → dedup
   - atomicAdd per row → counts
   - `thrust::exclusive_scan` → rowPtr
   - Copy colInd
3. **Diagonal lookup** (L247–270): per row, binary search colInd for
   `colInd[k] == row`

Complexity: O(E log E), E = numElements × 24. At 10^9 elements,
E ≈ 2.4×10^10. Thrust GPU sort ~10× faster than CPU.

## 6. Twelve Kernel Variants — Cross-Comparison

| Variant | File | Block | Threads/Elem | Reg | Best When |
|---------|------|-------|--------------|-----|-----------|
| Original (Graph) | mars_cvfem_hex_kernel_graph.hpp | 256 | 1 | ~255 | never |
| Optimized | mars_cvfem_hex_kernel_optimized.hpp | 256 | 1 | ~240 | baseline+10% |
| Shmem | mars_cvfem_hex_kernel_shmem.hpp | 256 | 1 | ~80 | memory-bound |
| Team | mars_cvfem_hex_kernel_optimized.hpp:268 | 256 (8 warps) | 32 | ~160 | extreme contention |
| **Tensor** | mars_cvfem_hex_kernel_tensor.hpp | 256 | 1 | ~210–230 | **CURRENT BEST: 4.60 ms / 10^8 hex on GH200** |
| TensorColored | mars_cvfem_hex_kernel_tensor_colored.hpp | 256 | 1 | ~210–230 | atomics costly |
| TensorAoS | mars_cvfem_hex_kernel_tensor_aos.hpp | 256 | 1 | ~210–230 | nodes dominate |
| TensorPerip | mars_cvfem_hex_kernel_tensor_perip.hpp | 256 | 1 | ~170–200 | regs constrained |
| TensorPeripLb2 | mars_cvfem_hex_kernel_tensor_perip.hpp:277 | 128 | 1 | ~140–170 | NOT GH200 |
| SmemCache | mars_cvfem_hex_kernel_smem_cache.hpp | 256 | 1 | ~180–210 | L2 latency dominates |
| WmmaTensor | mars_cvfem_hex_kernel_tensor_wmma.hpp | 256 (1 warp/elem) | 32 | ~120 | SM80+ FP64 WMMA |
| WgmmaTensor | (incomplete) | 128 | 16 (2 elem/warp) | ~100–130 | SM90+ Hopper async |

### 6.1 Tensor Kernel (current best) — `mars_cvfem_hex_kernel_tensor.hpp:30–291`

**Why it wins**:
1. Full local matrix `lhs[8×8]` in registers
2. Optimized loads (L87–136): unroll-by-4 for ILP and load reordering
3. Direct scatter via atomic; sorted CSR cols → binary search (L271–288)
4. Minimal spilling (~210–230 regs; GH200 has 256 KB/SM regs)

**Launch**:
```cuda
numBlocks = (numElements + 255) / 256;
cvfem_hex_assembly_kernel_tensor<..., 256><<<numBlocks, 256, 0, stream>>>(...);
```
Thread-per-element model. No intra-block sync.

**Measured**: 4.60 ms for 10^8 hex on GH200.

**Roofline**:
- ~3000 FLOPs/elem (12 SCS × shape evals + Jacobian + advection + diffusion)
- ~500 B/elem mem (with AoS nodes)
- 6 FLOP/byte; memory roof = 3.9 TB/s × 6 = 23 TFLOPS
- Actual = 10^8 × 3000 / 4.6 ms ≈ 65 TFLOPS → **21% of FP64 peak**
- Latency- and atomic-bound, not memory-bound

### 6.2 TensorColored — `mars_cvfem_hex_kernel_tensor_colored.hpp:25–225`

Eliminates atomics. Same as tensor but `+=` direct writes (lines 195, 209,
217). Driven by per-color element list:
```cuda
const int linearIdx = blockIdx.x * blockDim.x + threadIdx.x;
const int elemIdx = d_elemList[linearIdx];
```

Dispatch loop (`mars_cvfem_assembler.hpp:165–188`):
```cuda
for (int c = 0; c < col->numColors; ++c) {
    int cBlocks = (col->colorSize(c) + blockSize - 1) / blockSize;
    cvfem_hex_assembly_kernel_tensor_colored<256>
        <<<cBlocks, blockSize, 0, stream>>>(col->colorStart(c), ...);
}
```
**Speedup**: 10–15% on contended assembly.

### 6.3 TensorAoS — `mars_cvfem_hex_kernel_tensor_aos.hpp:26–216`

One 72-byte AoS read per node vs 9 separate loads. **3× L2 utilization**
improvement with same runtime (gather-bound case).

### 6.4 TensorPerip — `mars_cvfem_hex_kernel_tensor_perip.hpp:33–225`

Pre-compute all 64 CSR positions before SCS loop:
```cuda
int pos[64];   // L98
for (int i = 0; i < 8; ++i) {
    int row_start = matrix->rowPtr[dofs[i]];
    for (int j = 0; j < 8; ++j) {
        pos[i*8 + j] = binary_search(...);
    }
}
```
Then in 12-iter SCS loop:
```cuda
int pL = pos[nodeL*8 + n];
if (pL >= 0) atomicAdd(&matrix->values[pL], diff_coeff);
```
Saves 128 regs (lhs[64] doubles → 64 ints), occupancy +10%, ~5% slower
overall.

**__launch_bounds__(256, 2) variant** (L277–319): forces 2 CTAs/SM. Worse
on GH200 (prefers wide occupancy).

### 6.5 SmemCache — `mars_cvfem_hex_kernel_smem_cache.hpp:96–349`

Block hash-table dedupes ~400 unique nodes across 256 elements, caches in
shared memory:

**Phase 1–2** (L134–165): hash-table insert per node ID via `atomicCAS`:
```cuda
int h = smem_hash(nid, HT_MASK);
while (true) {
    int prev = atomicCAS(&s.hash_keys[h], -1, nid);
    if (prev == -1 || prev == nid) break;
    h = (h + 1) & HT_MASK;
}
```

**Phase 3–4** (L169–206): sequential slots, cooperative load of node data:
```cuda
double* base = s.node_data + sl * 9;
base[0] = d_x[nid]; base[1] = d_y[nid]; ...
s.row_start[sl] = matrix->rowPtr[dof];
```

**Phase 5–6** (L212–349): SCS loop reads from smem (~5 cyc) instead of L2
(~30 cyc).

Smem footprint: ~53 KB; with `__launch_bounds__(256, 2)`: 106 KB << 228 KB
GH200 limit.

### 6.6 WmmaTensor — `mars_cvfem_hex_kernel_tensor_wmma.hpp:64–387`

Diffusion = matrix product S^T × D:
- S[12×8]: per-SCS, nodeL ↔ +1, nodeR ↔ −1
- D[12×8]: per (SCS, node) diff_coeff
- C = S^T[8×12] × D[12×8] = lhs_diff[8×8]

Uses FP64 WMMA 8×8×4 (SM80+). 1 warp / element, 8 elements / block.

Three 8×8×4 mma_sync calls = 1536 FLOPs in ~30 cycles.

**Reality**: FP64 WMMA = 67 TFLOPS on GH200 vs 312 TFLOPS scalar peak.
Matmul too small (1536 FLOP) → latency-dominated. **2–3× slower than
scalar tensor.** Not a win at this size.

### 6.7 WgmmaTensor — Hopper-only

Async WGMMA m64n8k4. 2 elements/warp, 75% util. Incomplete; experimental.

## 7. Coloring — `mars_cvfem_coloring.hpp:48–157`

Greedy (Jones-Plassmann variant on CPU pre-pass).

1. **Build node→elements CSR** (L64–81): count, scan, fill
2. **Greedy color** (L88–122):
```cpp
std::vector<int> elemColor(numElements, -1);
std::vector<int> usedBuf(32, 0);
int tick = 1;
for (size_t e = 0; e < numElements; ++e) {
    for (int n = 0; n < 8; ++n) {
        size_t v = conn[n][e];
        for (int k = nodeElemOff[v]; k < nodeElemOff[v+1]; ++k) {
            int nb = nodeElem[k];
            if (nb != e && elemColor[nb] >= 0)
                usedBuf[elemColor[nb]] = tick;
        }
    }
    int c = 0;
    while (usedBuf[c] == tick) ++c;
    elemColor[e] = c;
    if (c >= numColors) numColors = c + 1;
    ++tick;
}
```
**Tick trick**: avoid resetting usedBuf each iteration → O(numElements)
total instead of O(numElements × maxColor).

3. **Permutation** (L124–140), **upload to device** (L142–156).

Result: typically 8–12 colors unstructured; exactly 8 for structured hex
(parity i+j+k).

## 8. Dispatch & Config — `mars_cvfem_assembler.hpp`

### Enum (L21–34)
```cpp
enum class CvfemKernelVariant {
    Original, Optimized, Shmem, Team, Tensor, TensorColored, TensorAoS,
    TensorPerip, TensorPeripLb2, SmemCache, WmmaTensor, WgmmaTensor
};
```

### Config (L41–47)
```cpp
struct Config {
    int blockSize = 256;
    CvfemKernelVariant variant = CvfemKernelVariant::Original;
    cudaStream_t stream = 0;
    const CvfemColoringData* coloring = nullptr;   // required for TensorColored
    const NodeData* nodeData = nullptr;             // required for TensorAoS
};
```

### Dispatch (L49–309): `assembleGraphLump()` parses variant, computes
grid, launches kernel.

**SmemCache special case** — needs explicit cudaFuncSetAttribute:
```cuda
cudaFuncSetAttribute(scPtr, cudaFuncAttributeMaxDynamicSharedMemorySize, smemSz);
cudaFuncSetAttribute(scPtr, cudaFuncAttributePreferredSharedMemoryCarveout, 100);
scPtr<<<numBlocks, blockSize, smemSz, stream>>>(...);
```

## 9. Bottleneck Analysis (Tensor on GH200, 10^8 hex)

Estimated breakdown of 4.60 ms:

| Phase | Time | % | Notes |
|-------|------|---|-------|
| Gather (nodes + connectivity) | 0.8 ms | 17% | scattered L2 reads |
| Jacobian + dndx | 1.2 ms | 26% | det+invert latency-bound, no ILP |
| SCS loop (adv + diff) | 1.5 ms | 33% | atomicAdd contention |
| CSR scatter | 1.1 ms | 24% | binary search + atomic |

**Optimization targets**:
1. TensorColored — drop atomics → ~10% faster
2. TensorAoS — reduce gather L2 traffic → ~5% faster
3. SmemCache — gather → smem reads → ~5% faster
4. **Async MPI overlap** — overlap halo with compute → bigger win at scale

Not worth pursuing: WmmaTensor (matmul too small), WgmmaTensor (only Hopper
async, marginal at best).

## 10. Variant Selection Quickref

| Scenario | Best Variant |
|----------|--------------|
| Standard 10^9 hex | Tensor |
| Atomic contention | TensorColored |
| Gather-bound | TensorAoS |
| Reg-pressure | TensorPerip |
| L2-bound nodes | SmemCache |
| FP64 tensor cores | not recommended (WmmaTensor too small) |
| Hopper async | WgmmaTensor (experimental) |

## 11. File Index
```
fem/
  mars_cvfem_assembler.hpp                       L17–502
  mars_cvfem_hex_kernel.hpp                      L9–481  (constants + baseline)
  mars_cvfem_hex_kernel_graph.hpp                graph variant
  mars_cvfem_hex_kernel_optimized.hpp            L56–476 (incl team L268)
  mars_cvfem_hex_kernel_shmem.hpp                L20–266
  mars_cvfem_hex_kernel_tensor.hpp               L29–294 (current best)
  mars_cvfem_hex_kernel_tensor_colored.hpp       L24–228
  mars_cvfem_hex_kernel_tensor_aos.hpp           L25–219
  mars_cvfem_hex_kernel_tensor_perip.hpp         L30–322
  mars_cvfem_hex_kernel_tensor_wmma.hpp          L64–387
  mars_cvfem_hex_kernel_smem_cache.hpp           L37–352
  mars_cvfem_coloring.hpp                        L44–157
  mars_cvfem_node_data.hpp                       L25–86
  mars_sparse_matrix.hpp
  mars_sparsity_builder.hpp                      L30–273
  mars_h1_fe_space.hpp
```
