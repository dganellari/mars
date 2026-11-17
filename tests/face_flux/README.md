# Face-flux projection checks

`useFaceFlux` keeps a flux per (element, sub-control surface). The corrector updates it with exactly
K's face form, and the pressure outlets keep their continuity equation: Robin rows closed by the
boundary face, `u · A_b`.

**Supported**, checked collectively by `validateFaceFlux` before a step changes anything:
- Hypre's assembled pressure operator, either DDT assembled into `Apre`, or K with the divT gradient
  (`useLegacyGradient = false`). The outlet rows assume the corrector's difference-form gradient; the
  legacy gradient carries its own boundary closure and would double the outlet correction. The outlet
  rows are row-only, so the solve uses GMRES; the `MARS_HYPRE_KRYLOV` override does not apply.
- **A pressure gauge:** eligible outlets (owned pressure-outlet nodes with free velocity and an
  outward area), whose rows fix the level, and no pin; or no eligible outlet and a pin, which is then
  a real row of the operator.

**Rejected:** matrix-free DDT CG, `MARS_HYPRE_USE_DDT`, the Rhie-Chow/VMS options, the FEM
projection, periodic boundaries, an outlet beta. Each raises `std::invalid_argument` on every rank.

**Each opening flux enters once.** The eligible outlets' `u · A_b` comes from the face flux. The
inlet and ineligible outlet nodes take theirs from the SCS term (`useOpeningFluxSource = false`) or
from the prescribed source (`true`). The final control-volume balance adds the same terms as the
pressure RHS.

**The self-check** runs once, before the first face-flux solve, and aborts on a mismatch. It compares
K's face form with the assembled operator on every row but the pin. On the outlet rows it includes
`-β g · A_b`, with `g` the corrector's own gradient.

```sh
cmake -S tests/face_flux -B build-face-flux -DCMAKE_PREFIX_PATH=/path/to/candidate/install
cmake --build build-face-flux
ctest --test-dir build-face-flux --output-on-failure
```

Use a CUDA or HIP installation of this candidate with Hypre enabled. Configure the `MPIEXEC_*`
CMake variables for the system's launcher.

The fifteen tests run five modes at 1, 2 and 4 ranks. The boxes are 4³ cubes: five Tet4 each, or
one Hex8 with warped interior nodes. Each has inflow at x0, a do-nothing outlet at x1 and no-slip
walls, with the openings' areas given per node as the rectangular dual of the face grid.

| mode | checks |
|---|---|
| `tet_scs`, `hex_scs` | `useOpeningFluxSource = false`, where the outlet flux used to enter twice |
| `tet_source`, `hex_source` | `useOpeningFluxSource = true`, as ldp_fsi runs |
| `reject` | K with the legacy gradient, CG, `MARS_HYPRE_USE_DDT`, each set on one rank only; and no outlet without a pin |

In the four flow modes, every step checks that the pressure solve converges, every control volume
balances to 1e-7 of the inflow, the inlet carries Uinf times its area, and the openings' fluxes sum
to zero to 1e-7. The first step also passes the operator self-check, including the outlet rows.
