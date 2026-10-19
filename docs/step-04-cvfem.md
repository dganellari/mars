# Step 4: control volumes

Step 4 solves the same problem as [step 3](step-03-poisson.md),

    −Δu = f  in Ω = (0,1)³,    u = 0  on the boundary,    u = sin(πx) sin(πy) sin(πz),

with a different method: the control-volume finite element method (CVFEM) on hexahedra. This is the
discretization of MARS's Navier–Stokes solver, so this step is the bridge from the textbook problem to
the flow examples later on.

The program is `examples/distributed/unstructured/mars_step04_cvfem.cu`. Compared with step 3 only
the matrix and the right-hand side change; the degrees of freedom, the boundary conditions, the solver
and the error are the same code.

## Run it

```bash
for n in 8 16 32; do
  python3 scripts/generate_hex_cube.py -n $n --output hex$n
  mpirun -np 4 ./build/examples/distributed/unstructured/mars_step04_cvfem --mesh=hex$n
done
```

| Mesh | h | unknowns | CG iterations | max error | error ratio |
|---|---|---:|---:|---:|---:|
| `hex8` | 1/8 | 729 | 8 | 1.295e-2 | |
| `hex16` | 1/16 | 4,913 | 23 | 3.219e-3 | 4.02 |
| `hex32` | 1/32 | 35,937 | 52 | 8.036e-4 | 4.01 |

The error drops by 4 per refinement: second order, like the finite elements of step 3. The errors are
the same on one rank and on four. The release checks `marsReleaseStep04Cvfem16_npN` and
`marsReleaseStep04Cvfem32_npN` require the errors of the table within 5%.

## The idea

The finite element method of step 3 asks that the equation hold on average against every hat function
φ_i. CVFEM asks instead that it hold on average over a small box around each node, its *control
volume*. Each hexahedron is cut into eight sub-boxes, one per corner, and node i's control volume is
made of the sub-boxes at node i from all elements that touch it.

Integrating −Δu = f over the control volume V_i and using the divergence theorem gives

    −∮ ∇u · n dS = ∫ f dV    over the boundary and the inside of V_i.

The left side is the flux of u out of the box. Inside each element, u is the trilinear interpolant of
its corner values, so the flux through each of the 12 inner faces of the element is a linear
combination of the 8 corner values. Adding those up gives row i of the matrix. The right side is
approximated by f(x_i) V_i.

Each balance is exact for its box: what leaves one control volume enters its neighbour. That is why
flow solvers like the method: mass is conserved to round-off, on every mesh.

## The code

### The matrix

```cpp
SparseMatrix<int, double, cstone::execution::Gpu> A;
assembleCvfemDiffusion(domain, A, nodeToDof.data(), numOwned);
```

`assembleCvfemDiffusion` (in `fem/mars_cvfem_diffusion.hpp`) builds the sparsity (each node couples
to the 27 nodes of the elements around it), computes the area vectors of the inner faces, and runs the
same GPU assembly kernel as the Navier–Stokes solver, with diffusion only. A fifth argument sets a
constant diffusion coefficient γ, for −∇·(γ∇u) = f.

### The right-hand side

```cpp
cstone::DeviceVector<double> volume;
cvfemControlVolumes(domain, volume);
...
if (i < numOwned) rhs[i] = 3.0 * M_PI * M_PI * exactSolution(x[n], y[n], z[n]) * V[n];
```

`cvfemControlVolumes` gives every node the volume of its control volume: each hexahedron adds one
eighth of its volume to each corner. The right-hand side is then one GPU loop over the nodes.

## Try this

- Pass `2.0` as γ to `assembleCvfemDiffusion`: the computed u halves, because the exact solution of
  −∇·(2∇u) = f is u/2.
- Run step 3 and step 4 with the same n and compare the errors and how fast they drop.

Next: [step 5](step-05-galerkin.md) writes its own element operator.
