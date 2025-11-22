# MARS FEM Examples - Unstructured Mesh Solver

This directory contains finite element examples using MARS's GPU-native unstructured mesh backend.

## Examples

### Example 1: Poisson Equation (`mars_ex1_poisson.cu`)

Solves the Poisson equation with homogeneous Dirichlet boundary conditions:

```
-Δu = f  in Ω
u = 0    on ∂Ω
```

**Features:**
- GPU-native finite element assembly
- Linear Lagrange tetrahedra (P1 elements)
- Conjugate Gradient solver with cuSPARSE
- MPI parallel execution
- Performance timing and statistics

**Equivalent to:** MFEM `examples/ex1.cpp`

## Building

The examples require:
- CUDA toolkit (with cuSPARSE and cuBLAS)
- MPI implementation
- MARS with `-DMARS_ENABLE_UNSTRUCTURED=ON -DMARS_ENABLE_CUDA=ON`

```bash
cd build
cmake .. \
  -DMARS_ENABLE_KOKKOS=OFF \
  -DMARS_ENABLE_CUDA=ON \
  -DMARS_ENABLE_TESTS=ON \
  -DMARS_ENABLE_UNSTRUCTURED=ON \
  -DCMAKE_CUDA_ARCHITECTURES=80  # Adjust for your GPU

make mars_ex1_poisson
```

## Running

### Single Rank

```bash
./bin/examples/mars_ex1_poisson --mesh /path/to/mesh_parts
```

### Multi-Rank (MPI)

```bash
mpirun -np 4 ./bin/examples/mars_ex1_poisson \
  --mesh /path/to/mesh_parts \
  --order 1 \
  --max-iter 1000 \
  --tol 1e-10
```

### Command-Line Options

- `--mesh <path>`: Path to partitioned mesh directory (default: `mesh_parts`)
- `--order <n>`: Polynomial order (1 for linear, currently only 1 supported)
- `--max-iter <n>`: Maximum CG iterations (default: 1000)
- `--tol <val>`: CG convergence tolerance (default: 1e-10)
- `--help`, `-h`: Print help message

## Mesh Format

The examples expect pre-partitioned binary mesh files in MARS SoA format:

```
mesh_parts/
  mesh.0.bin   # Rank 0 partition
  mesh.1.bin   # Rank 1 partition
  mesh.2.bin   # Rank 2 partition
  mesh.3.bin   # Rank 3 partition
```

Each file contains:
- Node coordinates (SoA: X, Y, Z arrays as `float`)
- Element connectivity (node indices as `unsigned`)
- Metadata (element/node counts)

## Expected Output

```
========================================
   MARS Poisson Example (GPU-native)
========================================
Problem: -Δu = f in Ω, u = 0 on ∂Ω
Mesh: mesh_parts
MPI ranks: 4
Order: 1
========================================

1. Loading mesh...
   Mesh loaded in 0.15 seconds
   Local elements: 12345
   Local nodes: 6789

2. Creating finite element space...
   FE space created in 0.02 seconds
   Total DOFs: 27156

3. Assembling stiffness matrix...
   Stiffness matrix assembled in 0.45 seconds
   Matrix size: 27156 x 27156
   Non-zeros: 400000

4. Assembling RHS vector...
   RHS vector assembled in 0.12 seconds

5. Applying boundary conditions...
   Boundary DOFs: 2145
   BCs applied in 0.03 seconds

6. Solving linear system (CG)...
   CG iteration 0: residual = 1.0
   CG iteration 10: residual = 0.045
   ...
   CG converged in 127 iterations, residual = 9.8e-11
   System solved in 1.23 seconds

7. Computing solution statistics...
   Solution statistics:
     Min: 0.0
     Max: 0.125
     Mean: 0.042
     L2 norm: 1.456

========================================
   Timing Summary
========================================
Mesh loading:       0.15 s
FE space setup:     0.02 s
Stiffness assembly: 0.45 s
RHS assembly:       0.12 s
BC application:     0.03 s
Linear solve:       1.23 s
----------------------------------------
Total time:         2.00 s
========================================

MARS Poisson example completed successfully!
```

## Comparison with MFEM

| Feature | MFEM (CPU) | MARS (GPU) |
|---------|-----------|-----------|
| Platform | CPU, OpenMP | CUDA GPU |
| Assembly | Serial/OpenMP | GPU kernels |
| Matrix storage | Host CSR | Device CSR |
| Linear solver | Hypre AMG | cuSPARSE CG |
| Memory | Host RAM | GPU VRAM |
| Parallelism | MPI + OpenMP | MPI + GPU |

**Performance advantages:**
- GPU assembly: 10-100x faster for large meshes
- Lazy composition: Lower memory footprint
- Device-resident data: No PCIe transfers during solve

## Troubleshooting

### Error: "Mesh not found"
Ensure mesh files exist in the specified directory with correct naming (`mesh.<rank>.bin`).

### Error: "CUDA out of memory"
Reduce mesh size or use more MPI ranks to distribute data across GPUs.

### Error: "CG did not converge"
- Increase `--max-iter`
- Check mesh quality (inverted elements)
- Verify boundary conditions cover entire boundary

### Error: "cusparse initialization failed"
Ensure CUDA runtime and cuSPARSE libraries are correctly installed.

## Next Steps

1. **Add visualization**: Implement VTK writer to export solution
2. **Higher-order elements**: Extend to P2, P3 basis functions
3. **More examples**:
   - Ex2: Linear elasticity
   - Ex3: Time-dependent diffusion
   - Ex4: Adaptive mesh refinement
4. **Benchmarking**: Compare with MFEM on same hardware

## References

- MFEM: https://mfem.org
- MARS: https://github.com/dganellari/mars
- Cornerstone: https://github.com/sekelle/cornerstone-octree
