# Mixed flow interface methods

The mixed flow driver solves steady, constant-density incompressible flow with
nodal Hex8/Tet4 CVFEM and SIMPLE pressure correction. Velocity is in m/s,
pressure is physical pressure in Pa, density is in kg/m³ and dynamic viscosity
is in Pa·s when mesh coordinates are metres. The momentum stress is
`-p I + mu (grad(u) + grad(u)^T)`; mass flux is positive out of the first side.
The stabilized mass flux, rather than raw velocity divergence, enters continuity.
An open outlet retains signed quadrature fluxes; the backflow rule closes a whole
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
