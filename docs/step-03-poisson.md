# Step 3: Poisson's equation

Step 3 solves a partial differential equation on the GPU: Poisson's equation

    −Δu = f  in Ω = (0,1)³,    u = 0  on the boundary,

with linear finite elements on tetrahedra and the conjugate gradient method (CG). The right-hand side
is chosen so that the exact solution is known,

    u = sin(πx) sin(πy) sin(πz),    f = 3π² u,

so the program can print its own error. Refining the mesh then shows how fast the error goes down.

The program is `examples/distributed/unstructured/mars_step03_poisson.cu`. It follows the same
pattern as the larger `mars_ex1_poisson` example: mesh, degrees of freedom, matrix, right-hand side,
boundary conditions, solve.

## Run it

```bash
for n in 8 16 32; do
  python3 scripts/generate_tet_cube.py --nx $n --ny $n --nz $n --output tet$n
  mpirun -np 4 ./build/examples/distributed/unstructured/mars_step03_poisson --mesh=tet$n
done
```

| Mesh | h | unknowns | CG iterations | max error | error ratio |
|---|---|---:|---:|---:|---:|
| `tet8` | 1/8 | 729 | 17 | 2.520e-2 | |
| `tet16` | 1/16 | 4,913 | 39 | 6.394e-3 | 3.94 |
| `tet32` | 1/32 | 35,937 | 76 | 1.604e-3 | 3.99 |

Halving the mesh size divides the error by 4: the method is second-order accurate, as expected for
linear elements. The errors are the same on one rank and on four, and on NVIDIA and AMD GPUs.

## The code

### Degrees of freedom

```cpp
cstone::DeviceVector<int> nodeToDof(nodeCount);
const int numOwned = buildDofMappingGpu<uint64_t>(domain.getNodeOwnershipMap().data(), nodeToDof.data(),
                                                  nodeCount);
```

With linear elements there is one unknown per mesh node. `buildDofMappingGpu` numbers the nodes this
rank owns from 0 to `numOwned - 1` and its ghost nodes ([step 1](step-01-mesh.md)) after them, so
every vector below holds the owned values first and the ghost values at the end.

### The matrix and the right-hand side

```cpp
SparseMatrix<int, double, cstone::execution::Gpu> A;
TetStiffnessAssembler<double, uint64_t>().assemble(domain, A, nodeToDof.data(), numOwned);

cstone::DeviceVector<double> b(numOwned);
TetMassAssembler<double, uint64_t>().assembleRHS(domain, b, Source{}, nodeToDof.data(), numOwned);
```

The stiffness matrix is A_ij = ∫ ∇φ_i · ∇φ_j and the load vector b_i = ∫ f φ_i, where φ_i is the
hat function of node i. Each rank assembles the rows
of the nodes it owns; a row's columns can be ghost nodes too. The matrix is stored on the GPU in CSR
format (row offsets, column indices, values).

`Source` is a small struct whose `operator()(x, y, z)` returns f. The assembler calls it on the GPU
at its quadrature points, so it is marked `__host__ __device__`.

### Boundary conditions

```cpp
if (!onBoundary) return;
for (int j = rowPtr[i]; j < rowPtr[i + 1]; ++j)
    val[j] = col[j] == i ? 1.0 : 0.0;
rhs[i] = 0.0;
```

A node on one of the six faces of the cube gets the equation u_i = 0: its row becomes a row of the
identity and its right-hand side zero. The loop runs over all local nodes on the GPU and skips ghosts,
because the owner of a node sets its row.

### Conjugate gradients

```cpp
ConjugateGradientSolver<double, int, cstone::execution::Gpu> cg(2000, 1e-10);
cg.setOwnedSize(numOwned);
if (numRanks > 1)
    cg.setHaloExchangeCallback(
        [&domain, dofMap](cstone::DeviceVector<double>& p) { domain.exchangeNodeHalo(p, dofMap); });
converged = cg.solve(A, b, u);
```

CG stops when the residual has dropped by a factor of 10¹⁰ (it starts from u = 0), or after 2000
iterations. On one rank
that is all. On several ranks, each matrix-vector product A·p reads the ghost entries of p, whose
values belong to other ranks. The callback runs before every product and brings them from their
owners through the domain's node halo exchange. This is the one line that makes the solve parallel.

### The error

```cpp
return i < numOwned ? fabs(uh[i] - exactSolution(x[n], y[n], z[n])) : 0.0;
```

For every owned node the program compares the computed value with the exact solution, takes the
largest difference on the GPU (`thrust::transform_reduce` with `thrust::maximum`) and the largest over
all ranks (`MPI_Allreduce` with `MPI_MAX`). The release checks `marsReleaseStep03Poisson8_npN` and
`marsReleaseStep03Poisson16_npN` require the errors of the table within 5%.

## Try this

- Change the exact solution to sin(2πx) sin(2πy) sin(2πz), and f to 12π² u: the
  error is larger on the same mesh, and still drops by 4 per refinement.
- Solve with `1e-6` instead of `1e-10`: the error stops improving once the solver error is larger
  than the discretization error.

Next: [step 4](step-04-cvfem.md) solves the same problem with control volumes, the method of the
flow solver.
