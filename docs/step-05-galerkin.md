# Step 5: your own operator

Steps 3 and 4 used operators that ship with MARS. Step 5 writes its own: the element matrix and
vector of the Poisson problem are about fifteen lines of the program, and MARS assembles them on the
GPU, on hexahedra or on tetrahedra. The same pattern works for any operator whose element matrix you
can write down: a reaction term, a variable coefficient, a mass matrix for time stepping.

The problem is the one of steps 3 and 4, −Δu = f on the unit cube with u = sin(πx) sin(πy) sin(πz).
The program is `examples/distributed/unstructured/mars_step05_galerkin.cu`.

## Run it

```bash
for n in 8 16 32; do
  python3 scripts/generate_hex_cube.py -n $n --output hex$n
  mpirun -np 4 ./build/examples/distributed/unstructured/mars_step05_galerkin --mesh=hex$n
done
# tet16 from step 3
mpirun -np 4 ./build/examples/distributed/unstructured/mars_step05_galerkin --mesh=tet16 --elements=tet
```

| Mesh | elements | h | unknowns | CG iterations | max error | error ratio |
|---|---|---|---:|---:|---:|---:|
| `hex8` | Q1 | 1/8 | 729 | 7 | 2.605e-2 | |
| `hex16` | Q1 | 1/16 | 4,913 | 15 | 6.447e-3 | 4.04 |
| `hex32` | Q1 | 1/32 | 35,937 | 31 | 1.608e-3 | 4.01 |
| `tet16` | P1 | 1/16 | 4,913 | 39 | 6.394e-3 | |

Trilinear hexahedra are second-order accurate, like linear tetrahedra. On tetrahedra the result is
that of step 3 to all printed digits: the same method, written by hand instead of called.

## The code

### The element kernel

```cpp
template<typename ElementTag>
struct PoissonElement
{
    static constexpr int N = ElementTag::NodesPerElement;

    __host__ __device__ void operator()(const double (&X)[N][3], double (&K)[N][N], double (&F)[N]) const
    {
        ElementValues<ElementTag, double> v;
        for (int q = 0; q < v.numQuadPoints; ++q)
        {
            v.reinit(X, q);
            const double f = 3.0 * M_PI * M_PI * exactSolution(v.x[0], v.x[1], v.x[2]);
            for (int i = 0; i < N; ++i)
            {
                F[i] += f * v.phi[i] * v.JxW;
                for (int j = 0; j < N; ++j)
                    K[i][j] += (v.grad[i][0] * v.grad[j][0] + v.grad[i][1] * v.grad[j][1]
                              + v.grad[i][2] * v.grad[j][2]) * v.JxW;
            }
        }
    }
};
```

The kernel gets the coordinates of the element's N corners (8 for a hexahedron, 4 for a tetrahedron)
and adds the element matrix K_ij = ∫ ∇φ_i · ∇φ_j and vector F_i = ∫ f φ_i into K and F, which start
at zero. The integrals are sums over quadrature points.

`ElementValues` does the geometry, as `FEValues` does in deal.II. After `reinit(X, q)` it holds, at
quadrature point q:

| member | meaning |
|---|---|
| `phi[i]` | the value of basis function i |
| `grad[i][k]` | its derivative in direction k, in physical coordinates |
| `x[k]` | the point itself |
| `JxW` | the quadrature weight times the volume factor \|det J\| |

The kernel runs on the GPU, one thread per element, so it is marked `__host__ __device__` and uses
only fixed-size arrays.

### Assembly

```cpp
SparseMatrix<int, double, cstone::execution::Gpu> A;
cstone::DeviceVector<double> b;
assembleElements(domain, A, b, nodeToDof.data(), numOwned, PoissonElement<ElementTag>{});
```

`assembleElements` (in `fem/mars_element_assembly.hpp`) builds the sparsity, calls the kernel for
every element the rank holds, and adds each K and F into the rows the rank owns. Rows of ghost nodes
are skipped: their owners add them. After the call A and b are ready for the boundary conditions and
the solver, which are the code of step 3.

### Hexahedra or tetrahedra

```cpp
elements == "tet" ? solve<TetTag>(meshDir, rank, numRanks) : solve<HexTag>(meshDir, rank, numRanks);
```

The element type is a template parameter of the domain, of `ElementValues` and of the kernel, so the
same source serves both.

## Try this

- Add a reaction term, −Δu + u = f: add `v.phi[i] * v.phi[j] * v.JxW` to K and
  `(3π² + 1) u_exact` as f. The error behaves as before.
- Multiply the gradient term by a coefficient k(x) evaluated at `v.x`: a variable-coefficient
  problem in one more line.
