# Steady conduction around a heated wing

This example uses the high-order scalar FEM operator and cell-wise p-multigrid.
It solves a physical boundary-value problem, rather than recovering a prescribed
interior field. The analytic slab and wing case have passed native NVIDIA GPU/MPI validation.

The public `large_hex_only.e` mesh describes the volume **outside** the wing.
The wing is a temperature-controlled surface; the surrounding volume is an
idealized stationary conducting medium. This is not conduction inside a metal
wing and not a flow or conjugate heat-transfer simulation.

## Model

Solve `-div(k grad T) = 0` with no internal heat generation:

| Exodus side set | Boundary condition |
|---|---|
| `blade` | T = 353.15 K (80 degrees Celsius) |
| `bottom`, `in`, `out`, `side`, `top` | T = 293.15 K (20 degrees Celsius) |
| `sym` | Insulated: `-k grad T . n = 0` |

The chosen conductivity is 1 W/(m K), and one mesh coordinate unit is interpreted
as one metre. These define an idealized model; they are not measured wing or air
properties. Advection, buoyancy and radiation are excluded. The inlet/outlet names
are geometric labels here, not flow boundary conditions.

The temperature difference is 60 K. With `theta = (T-293.15)/60`, the hot and cold
values are one and zero. Constant conductivity cancels from the temperature solve
and multiplies the reported heat rates.

Q3 conforming elements use five-point Gauss quadrature in each direction. The weak
form is `integral(k grad v . grad theta) = 0` for tests vanishing on the hot and cold
boundaries. The insulated plane contributes no boundary load. Dirichlet boundary
sets are nonempty, so the reduced diffusion operator is positive definite.

The implementation forms a continuous nodal lift `g` with the imposed boundary
values and zero interior values, then solves `K w = -K g` on the free nodes, with
`w=0` on Dirichlet nodes. The computed solution is `theta = g+w`. No analytic
interior temperature is supplied. Boundary marks are summed across element copies
and MPI ranks before updating the homogeneous masks used by every polynomial level.

`--heat-case wing` is specific to these seven named boundary sets. The prepared kit
pins the original mesh hash and independently verifies the complete boundary-face
classification. The node-set reconstruction is not a general Hex8 face-label API.

## Validation and output

The independent CPU reference assembles Q3 matrices on warped, rotated hexes. A unit
slab with hot x=0, cold x=1 and insulated remaining faces has the exact solution
`theta=1-x` and total heat rate 60 W for the stated units and conductivity. This
checks the lifting sign, mixed boundary conditions, assembly and heat accounting.
A wrong-sign control must fail. The GPU gate runs this slab on one and eight ranks
before the wing. Both slab runs and the eight-rank wing run passed.

The solver reports the true assembled free-node residual, boundary error, temperature
range, heat balance and the energy identity. Summed Dirichlet reactions give heat
entering at the wing and heat leaving the cold enclosure. These are discrete weak
heat rates, not a separately integrated surface-gradient error estimate. The energy
check uses `theta^T K theta = g^T K theta` up to the free residual. The kit requires
relative heat/energy defects <=1e-7, true residual <=1e-10 and boundary error <=1e-12.
Temperature extrema and the fraction of unique global GLL nodes outside the imposed
range by more than 1e-6 of the temperature difference are reported as diagnostics.
They do not fail the wing solve: high-order FEM does not guarantee a discrete maximum
principle. The fraction is a nodal count, not a physical volume fraction; exported
slice statistics separately count element-local copies. No values are clipped, and
the plot includes the measured range. A converged solve with overshoot still requires
physical interpretation before publication. The analytic slab retains its strict
profile and heat-rate gates.

The optional slice export copies only owned elements intersecting a chosen plane
and rectangle. It writes actual computed theta, geometry and the boundary lift after
all solve checks and timers. Rendering converts theta to degrees Celsius and uses
linear display subcells between GLL values; the solve remains Q3.

The example reuses ElementDomain, the distributed Exodus reader, generated BP3
operators, node-line smoothing, p=3 to p=1 transfers and BoomerAMG. No Navier–Stokes
code is needed. Full-application speedups and mesh-converged engineering heat loads
are not established by this example's implementation or CPU checks.

## Accepted GPU result

On eight GH200 GPUs, the 15,114,920-element wing case converged in 173 iterations
with a true assembled relative residual of 1.707e-14. The discrete heater power
and enclosure outflow were both 600.663 W, with relative imbalance 3.85044e-13,
energy defect 3.77769e-13 and zero prescribed-boundary error. The analytic slab
passed on one and eight GPUs, with maximum normalized temperature error
6.71685e-14. The original manufactured-solution mode also passed.

The checks used the same executable (SHA256
`9083968d0c2e27f9d6596c6c60a5c50b2747e0945f32ec9d6e9f0cf2f67f3278`)
and the public mesh SHA256
`d4b2db1b089ab0f726ab84c8b82b3521f817a72e6812989d821a3fe51f8d361c`.
The exported slice contains 28,533 uniquely owned elements; geometry, finite
values and boundary lifting were checked independently. These results establish
discrete consistency and conservation, not mesh-converged continuum accuracy.

Run the independent host reference with NumPy installed:

```bash
python3 tests/reference/heat_conduction/reference.py
```

With a CUDA/MPI/Hypre MARSIR build and the public mesh available, the wing driver
options are:

```bash
mars_marsir_bp --ptx-dir <generated-ptx-directory> --spec bp3 --p 3 \
  --mesh <large_hex_only.e> --domain --heat-case wing --mg --amg \
  --h-levels 0 --smoother lines --lines-ptx --amg-cycles 3 --amg-consolidate -1 \
  --tol 1e-13 --true-tol 1e-10 --maxit 2000 --reps 1 --solves 1
```

Use the site's MPI launcher and GPU binding for eight ranks. Optional
`--export-slice <existing-directory>/slice` writes one binary field and JSON
metadata file per rank after the numerical checks pass. Native HIP execution of
this example has not been validated.
