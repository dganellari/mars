# MARS Solvers & Sparse Matrix Infrastructure — Deep Reference

`backend/distributed/unstructured/solvers/` + `fem/mars_sparse_matrix.hpp`,
`fem/mars_sparsity_builder.hpp`, `fem/mars_h1_fe_space.hpp`. For dispatcher
and ElementDomain context see [02_mars_overview.md](02_mars_overview.md).

## 1. Solver Inventory

| File | Type | Key Capability |
|------|------|----------------|
| `mars_cg_solver.hpp` | CG | cuSPARSE SpMV, cuBLAS dot, halo callback hook |
| `mars_cg_solver_with_preconditioner.hpp` | PCG | abstract `Preconditioner` interface |
| `mars_distributed_cg_solver.hpp` | distrib CG | MPI-aware wrapper over CG |
| `mars_bicgstab_solver.hpp` | BiCGSTAB | non-symmetric, 6 work vecs, breakdown detection |
| `mars_gmres_solver.hpp` | GMRES(restart) | Krylov basis, Givens, Hessenberg QR |
| `mars_hypre_pcg_solver.hpp` | Hypre PCG | IJMatrix/IJVector + BoomerAMG |
| `mars_hypre_amg_preconditioner.hpp` | precond | standalone BoomerAMG wrapper |
| `mars_hypre_pcg_eliminated_solver.hpp` | Hypre + elim | DOF elimination strategy |
| `mars_hypre_pcg_simple_solver.hpp` | Hypre baseline | comparison vs eliminated |
| `mars_distributed_hypre_pcg_solver.hpp` | distrib Hypre | ParCSR partitioning |

## 2. CG Solver — `mars_cg_solver.hpp`

### Members
```cpp
int maxIter_;                              // default 1000
RealType tolerance_;                       // default 1e-10
bool verbose_;
std::function<void(Vector&)> haloExchangeCallback_;
int ownedSize_;                            // > 0 → MPI_Allreduce on dot products
cublasHandle_t cublasHandle_;
cusparseHandle_t cusparseHandle_;
```

### Constructor (L25–34)
```cpp
ConjugateGradientSolver(int maxIter = 1000, RealType tolerance = 1e-10)
    : maxIter_(maxIter), tolerance_(tolerance), verbose_(true),
      haloExchangeCallback_(nullptr) {
    cublasCreate(&cublasHandle_);
    cusparseCreate(&cusparseHandle_);
}
```

### Setup
- No matrix factorization (Krylov)
- Diagonal extracted lazily for Jacobi (L203–222 `extractDiagonalKernel`)
- Diag safety: `diag[i] = 1.0 if |diag[i]| < 1e-14 else diag[i]` (L58–65)

### Solve loop — L129–192

```cpp
for (int iter = 0; iter < maxIter_; ++iter) {
    if (haloExchangeCallback_) haloExchangeCallback_(p);   // 1. halo
    spmv(A, p, Ap);                                         // 2. SpMV
    RealType pAp = dot(p, Ap);                              // 3. dot (Allreduce if ownedSize_>0)
    RealType alpha = rho / pAp;
    axpy(alpha,  p,  x, x);                                 // 4. AXPY
    axpy(-alpha, Ap, r, r);
    jacobiPrecondition(diag, r, z);                         // 5. precond
    RealType r_norm = std::sqrt(dot(r, r));                 // 6. converge check
    if (r_norm / b_norm < tolerance_) return true;
    RealType beta = rho_new / rho;
    axpbyPartial(1.0, z, beta, p, m);                       // 7. update p
    rho = rho_new;
}
```

### Halo exchange — L132–134, L202–205
```cpp
void setHaloExchangeCallback(std::function<void(Vector&)> cb) { haloExchangeCallback_ = cb; }
```
For rectangular matrices (m owned rows, n total cols), `p` sized to `n`.
Callback fills `p[m:n]` with neighbor values before SpMV.

### SpMV (cuSPARSE) — L215–271
```cpp
cusparseCreateCsr(&matA, A.numRows(), A.numCols(), A.nnz(),
                  A.rowOffsetsPtr(), A.colIndicesPtr(), A.valuesPtr(),
                  CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                  CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F);
cusparseCreateDnVec(&vecX, x.size(), x.data(), CUDA_R_64F);
cusparseCreateDnVec(&vecY, y.size(), y.data(), CUDA_R_64F);

size_t bufferSize = 0;
cusparseSpMV_bufferSize(handle, NON_TRANSPOSE, &alpha, matA, vecX,
                        &beta, vecY, CUDA_R_64F, ALG_DEFAULT, &bufferSize);
void* buffer; cudaMalloc(&buffer, bufferSize);
cusparseSpMV(handle, NON_TRANSPOSE, &alpha, matA, vecX, &beta, vecY,
             CUDA_R_64F, ALG_DEFAULT, buffer);
cudaFree(buffer);
```

**⚠ Bug 1**: `CUSPARSE_INDEX_32I` limits to ~2.1B nonzeros. At 10^9
elements on 1000 ranks this is fine per-rank, but watch for
single-rank matrix indexing.

**⚠ Bug 2**: alloc/free of `buffer` per call. Cache it.

### Dot product — L279–307
```cpp
size_t n = (ownedSize_ > 0) ? size_t(ownedSize_) : std::min(x.size(), y.size());
RealType local = 0;
cublasDdot(cublasHandle_, n, x.data(), 1, y.data(), 1, &local);
if (ownedSize_ > 0) {
    int initialized = 0; MPI_Initialized(&initialized);
    if (initialized) {
        RealType global = 0;
        auto type = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
        MPI_Allreduce(&local, &global, 1, type, MPI_SUM, MPI_COMM_WORLD);
        return global;
    }
}
return local;
```

### Convergence — L163–182
Relative tolerance: `||r||_2 / ||b||_2 < tol`.

## 3. BiCGSTAB — `mars_bicgstab_solver.hpp:92–208`

For non-symmetric / saddle-point / convection-dominated.

```cpp
Vector r, r0, p, v, s, t;   // 6 work vecs
for (int iter = 0; iter < maxIter_; ++iter) {
    spmv(A, p, v);
    RealType alpha = rho / dot(r0, v);
    s = r - alpha*v;
    if (||s|| / ||b|| < tol) { x += alpha*p; return true; }
    spmv(A, s, t);
    RealType omega = dot(t, s) / dot(t, t);
    x += alpha*p + omega*s;
    r = s - omega*t;
    if (|r0^T r| < 1e-30) return false;            // breakdown (L102, 137, 178)
    rho_new = dot(r0, r);
    beta = (rho_new / rho) * (alpha / omega);
    p = r + beta*(p - omega*v);
    rho = rho_new;
}
```

vs CG:
- 2 SpMV/iter
- No symmetry needed
- Breakdown possible
- Memory: 6 vecs vs 4

## 4. GMRES(restart) — `mars_gmres_solver.hpp:39–231`

For ill-conditioned, singular, or when CG stalls.

```cpp
std::vector<Vector> V(restart_ + 1);   // (m+1) × n basis

for (int outer = 0; outer < maxIter_ / restart_; ++outer) {
    r = b - A*x;
    if (usePreconditioner) r = M^{-1} r;
    beta = ||r||;
    V[0] = r / beta;
    s[0] = beta;

    for (j = 0; j < restart_; ++j) {                 // Arnoldi
        w = A * V[j];
        if (usePreconditioner) w = M^{-1} w;
        for (i = 0; i <= j; ++i) {                   // mod Gram-Schmidt
            H[i][j] = w^T V[i];
            w -= H[i][j] * V[i];
        }
        H[j+1][j] = ||w||;
        V[j+1] = w / H[j+1][j];
        // Givens rotations on H[i:i+1, j], i = 0..j-1
        // New rotation
        cs[j] = H[j][j] / h_norm; sn[j] = H[j+1][j] / h_norm;
        s[j+1] = -sn[j]*s[j]; s[j] = cs[j]*s[j];
    }
    // Back-solve H y = s; x += V*y
}
```

Memory: (m+1) × n. Typical m=30. Per-iter: 1 SpMV + O(j²) inner products.

## 5. Preconditioners

### Jacobi (diagonal) — embedded in CG files

```cpp
class JacobiPreconditioner : public Preconditioner<...> {
    void setup(const Matrix& A) override {
        diag_ = A.getDiagonal();
        thrust::transform(diag_.begin(), diag_.end(), diag_.begin(),
                         [] __device__(RealType d) {
                             return std::abs(d) > 1e-14 ? d : RealType(1.0);
                         });
    }
    void apply(const Vector& r, Vector& z) override {
        thrust::transform(r.begin(), r.end(), diag_.begin(), z.begin(),
                         thrust::divides<RealType>());
    }
};
```

Build O(nnz). Apply O(numRows). Excellent for diagonally dominant.

### AMG (Hypre BoomerAMG) — `mars_hypre_amg_preconditioner.hpp`

```cpp
void setup(const Matrix& A) override {
    HYPRE_IJMatrixCreate(comm_, ilower_, iupper_, ilower_, iupper_, &A_hypre_);
    HYPRE_IJMatrixSetObjectType(A_hypre_, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(A_hypre_);
    for (int i = 0; i < m_; ++i) {
        // ... HYPRE_IJMatrixSetValues row-by-row
    }
    HYPRE_IJMatrixAssemble(A_hypre_);
    HYPRE_IJMatrixGetObject(A_hypre_, (void**)&parcsr_A_);

    HYPRE_BoomerAMGCreate(&precond_);
    HYPRE_BoomerAMGSetCoarsenType(precond_, 6);     // Falgout
    HYPRE_BoomerAMGSetRelaxType(precond_, 6);       // SymGS / hybrid
    HYPRE_BoomerAMGSetNumSweeps(precond_, 1);
    HYPRE_BoomerAMGSetMaxLevels(precond_, 20);
    HYPRE_BoomerAMGSetup(precond_, parcsr_A_, NULL, NULL);
}
```

Build O(nnz × log L). Apply O(L × nnz_coarse), 5–15 FLOPs/nnz typical.
Robust for SPD elliptic.

## 6. Hypre Integration — `mars_hypre_pcg_solver.hpp:46–78`

```cpp
template<typename KeyType>
bool solve(const Matrix& A, const Vector& b, Vector& x,
           IndexType globalDofStart, IndexType globalDofEnd,
           IndexType globalColStart, IndexType globalColEnd,
           const std::vector<KeyType>& localToGlobalDof);
```

- Row partition: `[globalDofStart, globalDofEnd)` — owned rows
- Col partition: `[globalColStart, globalColEnd)` — ghost cols, monotonic across ranks
- For rectangular matrices: row indices in owned range, column indices in
  global DOF space, off-diagonal remote entries properly indexed

## 7. Sparse Matrix CSR — `mars_sparse_matrix.hpp`

```cpp
template<typename IndexType, typename RealType, typename AcceleratorTag>
class SparseMatrix {
    using DeviceVectorIndex = typename mars::VectorSelector<IndexType, AcceleratorTag>::type;
    using DeviceVectorReal  = typename mars::VectorSelector<RealType,  AcceleratorTag>::type;
private:
    IndexType numRows_;        // m (owned DOFs in distrib)
    IndexType numCols_;        // n (owned + ghost)
    IndexType nnz_;
    DeviceVectorIndex d_rowOffsets_;   // [m+1]
    DeviceVectorIndex d_colIndices_;   // [nnz]
    DeviceVectorReal  d_values_;       // [nnz]
};
```

`IndexType` typically `int32_t`. **Limit at ~2.1B nonzeros** (cuSPARSE
INDEX_32I). int64 path needs CUSPARSE_INDEX_64I and `IndexType=int64_t`.

### Atomic add — L239–250
```cpp
template<typename RealType, typename IndexType>
__device__ inline void atomicAddSparseEntry(
    RealType* values, const IndexType* colIndices,
    IndexType rowStart, IndexType rowEnd,
    IndexType col, RealType value)
{
    int idx = findColumnIndex(colIndices, rowStart, rowEnd, col);
    if (idx >= 0) atomicAdd(&values[idx], value);
    // silent drop if column not in pattern
}
```
Column must exist in sparsity. `buildGraphSparsity` from connectivity
guarantees all (i,j) pairs within element nodes exist.

### Distributed off-diagonal layout
```
A = [ A_own  A_ghost ]
    [m × m  m × n_ghost]
```
Halo exchange fills `x[m:n]` before SpMV; `y[0:m]` accumulates contributions
from all columns.

## 8. Sparsity Builder — `mars_sparsity_builder.hpp`

### Graph (7 NNZ/row) — `buildGraphSparsity`
```cpp
size_t numEdges = numElements * 12 * 2;
buildSparsityEdgeListKernel<<<...>>>(
    d_conn0..7, d_nodeToDof, numElements,
    d_edgeRow.data(), d_edgeCol.data());

// Diag entries
thrust::sequence(d_diagRow.begin(), d_diagRow.end());
thrust::sequence(d_diagCol.begin(), d_diagCol.end());

// Merge edges + diag
// Remove invalid (-1)
// thrust::sort by (row, col)
// thrust::unique → nnz

// Histogram + scan → rowPtr
thrust::device_vector<int> d_rowCounts(numDofs, 0);
thrust::for_each(d_allRow.begin(), d_allRow.end(),
    [counts] __device__ (int row) { atomicAdd(&counts[row], 1); });
thrust::exclusive_scan(d_rowCounts.begin(), d_rowCounts.end(), d_rowPtr.begin());
```

### Full (27 NNZ/row)
```cpp
// emit all 8×8 = 64 pairs per element
for (int i = 0; i < 8; ++i)
    for (int j = 0; j < 8; ++j) {
        edgeListRow[base + i*8 + j] = dofs[i];
        edgeListCol[base + i*8 + j] = dofs[j];
    }
```

| Metric | Graph (7 NNZ) | Full (27 NNZ) |
|--------|---------------|---------------|
| Assembly | edge flux only | full element stiffness |
| SpMV | sparser (faster) | denser |
| Use case | Poisson, advection-diffusion | elasticity, coupled |

## 9. DOF Mapping — `mars_h1_fe_space.hpp`

```cpp
H1FESpace(Domain& domain, int order = 1)
    : domain_(domain), order_(order)
{
    const auto& ownership = domain_.getNodeOwnershipMap();
    size_t numNodes = domain_.getNodeCount();
    thrust::host_vector<uint8_t> h_ownership(numNodes);
    thrust::copy(/* D2H copy */);
    numDofs_ = 0;
    for (size_t i = 0; i < numNodes; ++i) {
        if (h_ownership[i] == 1 || h_ownership[i] == 2) numDofs_++;
        // ghost (==0) NOT counted
    }
}
```

Ownership semantics:
- `1`: owned (contributes to RHS)
- `2`: shared (on partition boundary; needed for assembly, dual-owned)
- `0`: ghost (read-only copy)

## 10. Boundary Conditions — `mars_boundary_conditions.hpp`

### Dirichlet by row replacement — L57–92
```cuda
__global__ void applyDirichletKernel(
    const IndexType* boundaryDofs, IndexType numBoundaryDofs, RealType bcValue,
    const IndexType* rowOffsets, const IndexType* colIndices,
    RealType* values, RealType* rhs, IndexType numRows)
{
    IndexType idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numBoundaryDofs) return;
    IndexType row = boundaryDofs[idx];
    rhs[row] = bcValue;
    for (IndexType i = rowOffsets[row]; i < rowOffsets[row+1]; ++i) {
        values[i] = (colIndices[i] == row) ? RealType(1.0) : RealType(0.0);
    }
}
```

**⚠ Symmetry concern**: Row replacement preserves A[i, :] for interior
rows but A[j, i] for boundary col j is non-zero, so matrix becomes
non-symmetric. Two-step elimination (kernels at L95–225) restores symmetry.

### Two-step elimination — L95–225
```cuda
// Step 1: subtract boundary contributions from interior RHS
__global__ void eliminateBCStep1Kernel(...) {
    if (!isBoundaryRow(row, boundaryDofs)) {
        for (each entry in row) {
            if (isBoundaryDof(col, boundaryDofs))
                rhs[row] -= values[i] * bcValue;
        }
    }
}

// Step 2: zero boundary rows/cols
__global__ void eliminateBCStep2Kernel(...) {
    // boundary rows → identity
    // interior rows → zero columns at boundary DOFs
}
```
Result: block-diagonal with boundary DOFs decoupled.

## 11. Per-Iteration Cost (10^6 elem/rank, 7-NNZ, 4 ranks)

| Op | Cost | Notes |
|----|------|-------|
| Halo exchange | 1–2 ms | latency-bound, ~100 KB ghost @ 50 GB/s |
| SpMV (cuSPARSE) | 2–3 ms | 7M nonzeros, memory-bound, ~16 GB/s utilization |
| Dot product | 0.1 ms local + 0.5 ms Allreduce | latency-bound |
| Precond apply (Jacobi) | 0.2 ms | thrust::transform |
| AXPY ×3 | 0.3 ms | cuBLAS small-vec, latency-bound |
| **Total** | **4–7 ms** | dominated by SpMV + halo |

Convergence: CG SPD elliptic ~√κ iters. κ=100–1000 → 10–30 iters typical.
Total solve: 40–210 ms.

Weak scaling p=4: SpMV flat, halo +1 ms, Allreduce +0.3 ms (log p),
~30% MPI overhead per iter.

## 12. Identified Bugs

1. **CUSPARSE_INDEX_32I** (L227 `mars_cg_solver.hpp`) limits to 2.1B nnz.
   Fix: parameterize, use INDEX_64I when `IndexType=int64_t`.
2. **SpMV buffer alloc/free per call** (L233–241). Fix: cache descriptor +
   buffer per solver instance.
3. **BC row replacement breaks symmetry**. Fix: use elimination Step1+Step2
   when symmetry critical.
4. **Halo column off-diagonal**: cstone halo exchange must run before SpMV
   (driven by `setHaloExchangeCallback`). If not set, off-rank columns
   read stale data. Sparsity builder uses `numTotalDofs` (owned + ghost)
   to ensure cols exist locally.

## 13. Solver Selection Quickref

| Problem | Best |
|---------|------|
| SPD elliptic, well-conditioned | CG |
| Non-symmetric, steady | BiCGSTAB |
| Ill-conditioned, singular | GMRES(30) |
| Large distributed elliptic | Hypre PCG + BoomerAMG |

## 14. Canonical Solve — `mars_cvfem_graph.cu`

### Domain setup (L119–132)
```cpp
ElementDomain<ElemTag, RealType, KeyType, cstone::GpuTag> domain(
    meshFile, rank, numRanks, true, bucketSize);
size_t nodeCount = domain.getNodeCount();
size_t elementCount = domain.getElementCount();
```

### DOF map (L138–144)
```cpp
cstone::DeviceVector<int> d_nodeToDof(nodeCount);
int numDofs = buildDofMappingGpu<KeyType>(
    d_nodeOwnership.data(), d_nodeToDof.data(), nodeCount);
int numTotalDofs = static_cast<int>(nodeCount);
```

### Sparsity (L149–189)
```cpp
const auto& d_conn = domain.getElementToNodeConnectivity();
cstone::DeviceVector<int> d_rowPtr(numTotalDofs + 1);
cstone::DeviceVector<int> d_colInd, d_diagPtr(numTotalDofs);
int nnz = CvfemSparsityBuilder<KeyType>::buildGraphSparsity(
    std::get<0>(d_conn).data(), ..., std::get<7>(d_conn).data(),
    elementCount, d_nodeToDof.data(), numTotalDofs,
    d_rowPtr.data(), nullptr, nullptr, 0);
d_colInd.resize(nnz);
CvfemSparsityBuilder<KeyType>::buildGraphSparsity(
    ..., d_colInd.data(), d_diagPtr.data(), 0);
```

### Assembly (L199–240+)
```cpp
domain.cacheNodeCoordinates();
const auto& d_x = domain.getNodeX();
// ... allocate d_rhs, d_values, build CSR matrix wrapper ...
// dispatch: CvfemHexAssembler::assemble(d_matrix, d_x, ..., kernelVariant);
```

### BC (from `mars_ex1_poisson.cu`, L223–250)
```cpp
std::vector<bool> isBoundaryDOF(totalLocalDofs, false);
for (size_t i = 0; i < h_x.size(); ++i)
    if (std::abs(h_x[i] - 0.0) < 1e-6 || std::abs(h_x[i] - 8.0) < 1e-6)
        isBoundaryDOF[i] = true;
std::vector<KeyType> boundaryDofs;
for (size_t i = 0; i < isBoundaryDOF.size(); ++i)
    if (isBoundaryDOF[i]) boundaryDofs.push_back(i);
BoundaryConditionHandler<...> bc;
bc.applyDirichlet(fes, K, b, boundaryDofs, 0.0);
```

### Solve (from `mars_ex1_poisson.cu`, L280–310)
```cpp
ConjugateGradientSolver<RealType, IndexType, AcceleratorTag> solver(maxIter, tolerance);
solver.setVerbose(true);
solver.setOwnedSize(numDofs);
// solver.setHaloExchangeCallback([&](Vector& p) { dofHandler.updateGhostDofValues(p); });
cstone::DeviceVector<RealType> x(numTotalDofs, 0.0);
bool converged = solver.solve(K, b, x);
```

## 15. File Index
```
solvers/
  mars_cg_solver.hpp                       L25–307 (single + distrib base)
  mars_cg_solver_with_preconditioner.hpp   PCG abstraction
  mars_distributed_cg_solver.hpp           MPI wrapper
  mars_bicgstab_solver.hpp                 L92–208
  mars_gmres_solver.hpp                    L39–231
  mars_hypre_pcg_solver.hpp                L46–78
  mars_hypre_amg_preconditioner.hpp        BoomerAMG wrapper
  mars_hypre_pcg_eliminated_solver.hpp     DOF elim variant
  mars_hypre_pcg_simple_solver.hpp         baseline
  mars_distributed_hypre_pcg_solver.hpp    distrib Hypre

fem/
  mars_sparse_matrix.hpp                   CSRMatrix struct, atomicAdd helpers
  mars_sparsity_builder.hpp                L30–273 graph + full
  mars_h1_fe_space.hpp                     H1FESpace, DOF mapping

utils/
  mars_boundary_conditions.hpp             L57–225 row replacement + elimination
```
