# Free-slip wall checks

`NSStepper::slipNodes[a]` lists the local nodes on a free-slip (symmetry) plane normal to axis `a`.
Their velocity component `a` is 0, and the other components stay free (zero shear, the natural
condition).

**Precedence:** a node with a full Dirichlet condition (wall, inlet, extra node) keeps it. A node on
two slip planes (an edge) has both normal components fixed.

**Where it acts** (`setupSlipWalls`, `applySlipWalls`): per component, the rows and columns of
that axis's slip DOFs become identity rows of a copy of the step's velocity matrix, and their RHS is
0. The copy keeps a symmetric matrix symmetric for CG, and a zero value needs no lift. After the
corrector, the component is set to 0 at those nodes again.

**What needs nothing new:** pressure. Slip walls are closed like walls (no flux through them), and
the divergence on the RHS sees zero normal velocity.

**Restrictions:** the pump BCs only (`BCKind::Pump`); anything else throws `std::invalid_argument`
on every rank.

`initFaceFluxFromVelocity(s)` fills the face-flux projection's stored flux from the current
velocity. Call it after setting an initial velocity; setup leaves the flux at 0.

```sh
cmake -S tests/slip_walls -B build-slip-walls -DCMAKE_PREFIX_PATH=/path/to/candidate/install
cmake --build build-slip-walls
ctest --test-dir build-slip-walls --output-on-failure
```

The test case is the 2D Taylor-Green vortex, an exact Navier-Stokes solution in [0, π]² with
free-slip walls on its four sides:

    u = sin x cos y e^(-2νt),  v = -cos x sin y e^(-2νt),  p = (cos 2x + cos 2y)/4 e^(-4νt)

It runs in a slab two cells thick with free slip on every face, ν = 0.1, from the exact field at
t = 0, for 25 steps of 0.002. The pressure is pinned at one interior node (`fluxNeumann`).

| mode | checks |
|---|---|
| `{tet,hex}_ff` | the face-flux projection and its convection (implicit advection) |
| `reject` | slip nodes with `BCKind::Channel` |

The flow modes run on 16² and 32² cells. Every step checks:
- the pressure solve succeeds;
- the fields are finite;
- the normal component is exactly 0 on every slip plane.

At the end:
- the 16² relative L2 velocity error is below 5e-2;
- the error falls at least as h^1.5;
- the kinetic energy matches the exact `e^(-4νt)` to 2%.

**Not covered: the nodal (DDT) projection in a closed box.** A box with no openings needs the
interior pin (`fluxNeumann`). On the nodal projection, its pressure solve does not converge on this
case:
- with Jacobi, GMRES stalls at its 2000-iteration cap;
- with BoomerAMG, it returns a solution of size 2.5e9 for a near-zero RHS, which MARS flags as
  contaminated by the null space.

The pin evidently does not reach the operator that is solved; the face-flux operator applies the pin
row itself. This is independent of the slip walls, and a closed box is not the pump's configuration.
