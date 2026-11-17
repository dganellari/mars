# Mixed flow interface methods

The mixed flow driver solves steady, constant-density incompressible flow with
nodal Hex8/Tet4 CVFEM and SIMPLE pressure correction. Velocity is in m/s,
pressure is physical pressure in Pa, density is in kg/m³ and dynamic viscosity
is in Pa·s when mesh coordinates are metres. The momentum stress is
`-p I + mu (grad(u) + grad(u)^T)`; mass flux is positive out of the first side.
The stabilized mass flux, rather than raw velocity divergence, enters continuity.
With the default `--profile mars`, an open outlet retains signed quadrature fluxes; the backflow rule closes a whole
face. Clipping individual contributions would make stored mass disagree with the
pressure equation. Convergence requires normalized continuity L1 as well as both
nonlinear residuals; conservative mode also requires physical boundary balance.

Select the interface method explicitly when recording a comparison:

```sh
mars_mixed_nonconformal_flow --partitioned-mesh INPUT --interface compatibility --advection upwind --output RESULT.csv
mars_mixed_nonconformal_flow --partitioned-mesh INPUT --interface conservative --advection upwind --output RESULT.csv
```

`compatibility` is the default. It retains OpenAccel-style independent-side
quadrature and interface flux terms, using Tri3/Quad4 face-parametric donor
scores. Equal scores use persistent face IDs within a band anchored to the
minimum score. The original OpenAccel executable retains the first equal-score
candidate encountered; matching its incidental search ordering is not promised.
Captured PI donor tables are diagnostic evidence, not production input.
Compatibility also uses OpenAccel's reduced implicit stencil: edge neighbours
in the volume and complete parent coupling at interfaces. Coefficients outside
that stencil are added to the same row's diagonal block, while the full flux
residual is retained. This changes the momentum diagonal used by Rhie–Chow.
The full CSR allocation remains so a genuinely missing column still fails;
this is coefficient folding, not a matrix-storage optimization.

`legacy` retains the earlier MARS containment-based donor score with deterministic
ties and the full implicit stencil. It is a regression control, not another
OpenAccel reference. Conservative mode also retains the full stencil.

`conservative` intersects the two faces' nodal median-dual patches. It evaluates
one flux at each common overlap centroid and scatters the same value and Jacobian
with opposite signs to the two owner rows. The unique patch owner computes the
flux; ghosts receive that value through the halo. This enforces interface
cancellation independently of which rank owns each side. The caller must provide
complete face-pair candidates, stable face IDs and both complete parent elements.

The conservative geometry currently supports static, coincident, planar convex
Tri3/Quad4 faces, including Tet4/Tet4, Hex8/Tet4 and Hex8/Hex8 parents. Warped,
folded, nonopposing or uncovered surfaces are rejected. It does not model curved,
moving, sliding or finite-gap interfaces. Centroid quadrature integrates physical
affine traces exactly; nonlinear flux accuracy still requires refinement checks.

`mixed_overlap.cpp` checks dual-area coverage, affine interpolation and gradients,
orientation, nonmatching refinement and local flux cancellation.
`mixed.cpp` checks flux derivatives and donor ordering. Integration acceptance
additionally requires complete surface coverage, physical boundary mass balance,
converged nonlinear residuals and one-/multi-rank field agreement on CPU and GPU.
OpenAccel field parity is a separate compatibility check, not the acceptance
target for the changed conservative discretization.

## Corrected OpenAccel profile

Build with `MARS_ENABLE_HYPRE=ON` and `MARS_ENABLE_PETSC=ON`, supplying shared
real-FP64 PETSc through `MARS_PETSC_ROOT` (and `MARS_PETSC_ARCH` for a source
build). The standalone `tests/nonconformal` CMake project also builds
`mars_mixed_nonconformal_flow_host` with those options. Set `HYPRE_DIR` there
to the directory containing `HYPREConfig.cmake`.

```sh
python3 scripts/prepare_mixed_nonconformal_mesh.py hybrid_hex_tet.exo --ranks 1 --output mesh-r1
mars_mixed_nonconformal_flow --partitioned-mesh mesh-r1 --profile openaccel --interface compatibility --advection high-resolution --iterations 5000 --output result
```

Use a separately prepared partition directory for each MPI size. The driver
retains the PI channel's prescribed inlet profile, material constants and
boundary labels; it does not read arbitrary OpenAccel YAML decks. This is a
stationary planar Hex8/Tet4 channel workflow, not a general pump/FSI driver.

`openaccel` uses the corrected component-major limiter, PETSc FGMRES with
block-Jacobi/ILU(0), and Hypre GMRES/Falgout. Momentum targets are relative
1e-6, absolute 1e-12 and 20 iterations; pressure targets are 1e-5, 1e-10 and
100. Both solves check the true residual of the original assembled operator.
Selected interface parents define the reference sparse allocation, including
stored zeros; any omitted nonzero coefficient fails. Matrix layout matters to
these preconditioners even when the underlying operator is unchanged.

These reference preconditioners execute on the CPU. GPU builds retain device
geometry, assembly, field updates and halos, but explicitly stage each assembled
matrix/RHS to persistent host buffers and return the correction. This is an
accuracy comparison profile, not a GPU performance profile. A CPU-only Hypre
library can run it; the default GPU profile requires GPU-enabled Hypre.

The profile also reproduces the reference's pointwise outlet clipping and its
momentum/pressure stopping rule. Its success message says **reference equation
convergence**. Reported continuity and boundary imbalance remain visible and
are not conservation acceptance. Combining it with `--interface conservative`
or `legacy` is rejected. Defaults retain signed outlet fluxes and continuity
acceptance.

Donors are constructed from geometry, never read from a captured map. The mesh
preparer assigns persistent face identity from the original element and side
ordinal, independent of side-set listing order. Selection uses the smallest
face score followed by the smallest identity within an absolute 1e-14 band.
The historical OpenAccel search resolves ties by traversal order, which can
change across platforms. Replay-free field parity therefore requires a corrected
reference with the same deterministic tie rule; it is not a promise to reproduce
every historical executable's incidental donor choices. Archive both reference
source changes and the input identity with any comparison.

With `MARS_MIXED_PARTITIONS` pointing to a directory containing `mesh-r1`,
`mesh-r2` and `mesh-r4`, the optional host CTests run the complete reference
profile at all three rank counts. CPU tests do not replace a native GPU run.

The PI case passed native CUDA validation on 1/2/4 ranks. CPU and GPU exported
fields were byte-identical at each matching rank count. One-rank MARS and the
corrected reference converged at iteration 1077, with maximum differences
3.26e-10 m/s in velocity and 6.62e-9 Pa in pressure across all 2662 nodes.
Use raw reference output (`corrected_boundary_values: false`) and allow enough
outer iterations; the original 1000-iteration ceiling stops before convergence.
No boundary masking, pressure offset or captured donor map enters this check.

The reference profile is not rank-invariant: block-Jacobi partitions change the
iteration trajectory, and outlet closure depends on that history. In the PI
host checks, one/two ranks match the one-rank reference, while four ranks converge
to a different field. Compare CPU and GPU at matching rank counts; do not infer
cross-rank field agreement from the convergence tests or widen the reference
field tolerance to hide this difference. Conservative/default rank-agreement
checks remain separate.
