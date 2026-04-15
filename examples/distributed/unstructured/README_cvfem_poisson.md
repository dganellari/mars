# MARS CVFEM Poisson Solver

Complete PDE solver using optimized CVFEM assembly for the Poisson equation.

## Problem

Solves `-Δu = f` with `u = 0` on boundary (Dirichlet BC).

This is analogous to MFEM examples ex0/ex1 but uses **CVFEM** (Control Volume Finite Element Method) instead of traditional FEM.

## Features

- **Optimized CVFEM assembly** with tensor kernel (4.60ms on GH200, 43% faster than baseline)
- **Geometric boundary detection** (identifies boundary nodes automatically)
- **Preconditioned Conjugate Gradient solver** with Jacobi preconditioning
- **GPU-accelerated** throughout (assembly + solve)
- **MPI support** for distributed meshes
- **Flexible kernel selection** (tensor, shmem, optimized, original)

## Usage

```bash
# Build (assuming CUDA environment)
cd gpu && make mars_cvfem_poisson

# Run on a hex mesh
mpirun -np 1 bin/mars_cvfem_poisson --mesh=block_128.exo

# Use different kernel variant
mpirun -np 1 bin/mars_cvfem_poisson --mesh=block_128.exo --kernel=tensor

# Adjust solver parameters
mpirun -np 1 bin/mars_cvfem_poisson --mesh=block_128.exo --source=2.0 --tol=1e-12

# Multi-GPU (if available)
mpirun -np 4 bin/mars_cvfem_poisson --mesh=block_128.exo --kernel=tensor
```

## Command-Line Options

### Required
- `--mesh=FILE` - Mesh file (.mesh or .exo format)

### Optional
- `--kernel=VARIANT` - Kernel variant: `tensor` (default, fastest), `shmem`, `optimized`, `original`
- `--source=VALUE` - Source term f in `-Δu = f` (default: 1.0)
- `--max-iter=N` - CG maximum iterations (default: 1000)
- `--tol=VALUE` - CG convergence tolerance (default: 1e-10)
- `--block-size=N` - CUDA block size (default: 256)

## Output

The solver reports:
1. **Assembly performance** - Time to build matrix and RHS
2. **Solver convergence** - CG iterations and residual norms
3. **Solution statistics** - Min, max, and L2 norm of solution

Example output:
```
========================================
MARS CVFEM Poisson Solver
========================================
Problem: -Δu = 1.0, u = 0 on boundary
Mesh: block_128.exo
Kernel: tensor
MPI ranks: 1
========================================

Mesh loaded:
  Nodes:    2146689
  Elements: 2097152

DOF mapping:
  Owned DOFs: 2146689

Sparsity pattern:
  NNZ: 57960075
  Avg NNZ/row: 27

Assembling system...
Assembly completed in 4.60 ms

Boundary conditions:
  Boundary DOFs: 132096

Solving with CG...
Iteration : 0,   ||r||_B = 1.234e+00,   ||r||_2 = 1.234e+00
Iteration : 1,   ||r||_B = 5.432e-01,   ||r||_2 = 5.432e-01
...
CG converged in 47 iterations
Average reduction factor = 0.891

========================================
Solver Results
========================================
Converged: YES
Solve time: 123.45 ms
========================================

Solution statistics:
  Min:  0.000000e+00
  Max:  1.234567e-02
  L2 norm: 5.678901e+00
========================================
```

## Implementation Details

### CVFEM Formulation

The Poisson equation `-Δu = f` is solved in CVFEM form:
```
∇·(γ∇φ) = f
```

where:
- `φ = u` (solution)
- `γ = 1` (diffusion coefficient)
- `β = 0` (no advection term)
- `f` = source term

### Matrix Assembly

Uses `CvfemHexAssembler::assembleFull()` which:
1. Loads element data (vectorized for performance)
2. Computes SCS area vectors (12 per hex8 element)
3. Assembles full 8×8 element matrices
4. Scatters to global CSR matrix with atomics

### Boundary Conditions

Geometric detection finds nodes on domain boundary:
```cpp
onBoundary = (x == xmin || x == xmax ||
              y == ymin || y == ymax ||
              z == zmin || z == zmax)
```

Applied by:
- Zeroing matrix rows for boundary DOFs
- Setting diagonal to 1
- Setting RHS to 0 (for `u = 0` BC)

### Linear Solver

Preconditioned Conjugate Gradient with:
- **Jacobi preconditioner**: `M = diag(A)`
- **cuSPARSE** for SpMV
- **cuBLAS** for vector operations
- **MFEM-style convergence reporting**

## Differences from MFEM ex1

| Aspect | MFEM ex1 | MARS CVFEM Poisson |
|--------|----------|---------------------|
| Method | Traditional FEM (H1 space) | CVFEM (control volumes) |
| Elements | Tet/Hex with shape functions | Hex8 with SCS integration |
| Assembly | Element-wise, scatter-add | Vectorized kernel, atomics |
| Sparsity | ~7-19 NNZ/row (variable) | 27 NNZ/row (full hex8) |
| Solver | Hypre BoomerAMG + PCG | GPU-native Jacobi PCG |
| Performance | CPU-based | Full GPU (assembly + solve) |

## Expected Results

For a unit cube `[0,1]³` with `f = 1`:
- **Converges** in ~30-60 iterations (depends on mesh size)
- **Solution range**: `[0, umax]` where `umax ≈ O(h²)` for mesh size `h`
- **Residual**: Should reach tolerance `<1e-10`

## Performance Notes

- **Tensor kernel is fastest** (~4.60ms assembly for 2M elements on GH200)
- **Assembly dominates** for small meshes; solver dominates for large meshes
- **MPI scaling** depends on element distribution and ghost exchange

## Building

Add to CMakeLists.txt or Makefile:
```cmake
add_executable(mars_cvfem_poisson mars_cvfem_poisson.cu)
target_link_libraries(mars_cvfem_poisson mars cusparse cublas)
```

## See Also

- [CVFEM_QUICKSTART.md](../../CVFEM_QUICKSTART.md) - Quick start guide for CVFEM assembly
- [CVFEM_OPTIMIZATION_SUMMARY.md](../../CVFEM_OPTIMIZATION_SUMMARY.md) - Kernel optimization details
- [mars_cvfem_graph.cu](mars_cvfem_graph.cu) - Graph-based assembly benchmark
- [mars_cvfem_full.cu](mars_cvfem_full.cu) - Full assembly benchmark
- [mars_ex1_poisson.cu](mars_ex1_poisson.cu) - Traditional FEM Poisson solver
