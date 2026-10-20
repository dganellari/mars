# Step 6: the heat equation

The steps so far solved one linear system. Step 6 solves a time-dependent problem, the heat equation

    ∂u/∂t = Δu  in Ω = (0,1)³,    u = 0  on the boundary,    u(0) = u0 = sin(πx) sin(πy) sin(πz).

Its exact solution keeps the shape of u0 and decays: u(t) = e^(−3π²t) u0. The program prints the
error at the end time, so you can see how the error depends on the time step.

The program is `examples/distributed/unstructured/mars_step06_heat.cu`. It builds on
[step 5](step-05-galerkin.md): the matrices come from element kernels through `assembleElements`.

## Run it

```bash
python3 scripts/generate_hex_cube.py -n 32 --output hex32
for dt in 0.02 0.01 0.005; do
  mpirun -np 4 ./build/examples/distributed/unstructured/mars_step06_heat --mesh=hex32 --dt=$dt
done
```

| dt | time steps | CG iterations (all steps) | max error at t = 0.1 | error ratio |
|---|---:|---:|---:|---:|
| 0.02 | 5 | 141 | 4.581e-2 | |
| 0.01 | 10 | 262 | 2.285e-2 | 2.01 |
| 0.005 | 20 | 440 | 1.131e-2 | 2.02 |

Halving Δt halves the error: at this mesh size the error comes from the time steps, and implicit
Euler is first order. On the coarser `hex16` mesh the errors are within 4% of these (4.538e-2,
2.244e-2, 1.092e-2): the time-step error dominates, the mesh adds little. The errors are the same on one
rank and on four, and on NVIDIA and AMD GPUs. The release checks `marsReleaseStep06Heat10_npN` and
`marsReleaseStep06Heat20_npN` require the `hex16` errors at Δt = 0.01 and 0.005 within 5%.

## Implicit Euler

Galerkin finite elements in space turn the equation into M du/dt = −K u, with the mass matrix
M_ij = ∫ φ_i φ_j and the stiffness matrix K_ij = ∫ ∇φ_i · ∇φ_j. Implicit Euler replaces du/dt by
(u^(n+1) − u^n)/Δt and evaluates the right side at the new time:

    (M + Δt K) u^(n+1) = M u^n.

Every step is one linear system with the same matrix. Implicit Euler is stable for any Δt, and
first-order accurate: halving Δt halves the time error.

## The code

### Two element kernels

```cpp
K[i][j] += (v.phi[i] * v.phi[j] + dt * (v.grad[i][0] * v.grad[j][0] + ...)) * v.JxW;   // M + Δt K
K[i][j] += v.phi[i] * v.phi[j] * v.JxW;                                               // M
```

`HeatElement` gives the element matrix of M + Δt K, `MassElement` that of M. Both are assembled once,
before the time loop; the element vectors are not needed.

### The time loop

```cpp
SpmvPlan<double, Matrix, Vector> massTimes(handle, M, u, b);
for (int step = 0; step < steps; ++step)
{
    if (numRanks > 1) domain.exchangeNodeHalo(u, dof);
    massTimes.apply(u, b);
    ...                       // b = 0 on the boundary
    cg.solve(A, b, u);
}
```

The right-hand side M u^n needs u at the ghost nodes of the rank, whose values belong to other
ranks: `exchangeNodeHalo` brings them from their owners. `SpmvPlan` sets up the sparse product
y = M x once (the cuSPARSE or hipSPARSE descriptors and work buffer) and applies it at every step.
The number of steps is the smallest that keeps each step at most `--dt` and ends exactly at
`--time`; the program prints the step it used.

## Try this

- Make Δt large, say 0.05: the solution still decays, because implicit Euler is stable, but the
  error grows.
- Crank–Nicolson, (M + Δt/2 K) u^(n+1) = (M − Δt/2 K) u^n, is second order: change the two element
  kernels and watch the error drop by 4 when Δt halves.

Next: [step 7](step-07-cavity.md) solves the Navier–Stokes equations.
