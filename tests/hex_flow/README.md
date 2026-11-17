# Hex flow checks

The pump stepper runs on HEX8 meshes (`NSStepper<..., HexTag>`):
- exact sub-control volumes in the lumped mass;
- 2×2×2 Gauss convection and SUPG;
- the penalised projection.

Setup refuses an inverted or degenerate hex (`validateHexGeometry`). An element whose Jacobian
determinant at its centre or at a Gauss point, or one of whose sub-control volumes, is not finite and
positive raises `std::invalid_argument` on every rank. The kernels no longer leave such an element out
silently. Pressure stabilization: PSPG is tet-only; Bochev–Dohrmann is refused with the penalised
projection.

```sh
cmake -S tests/hex_flow -B build-hex-flow -DCMAKE_PREFIX_PATH=/path/to/candidate/install
cmake --build build-hex-flow
ctest --test-dir build-hex-flow --output-on-failure
```

Use a CUDA or HIP installation of this candidate with Hypre enabled. Configure the `MPIEXEC_*`
CMake variables for the system's launcher.

The mesh is a box of 4³ hexes with warped interior nodes: inflow at x0, a do-nothing outlet at x1,
no-slip walls. It runs the path ldp_fsi uses: Hypre, assembled DDT, divT gradient, implicit
convection, SUPG, and the penalised projection with β = 1/2.

| test | checks |
|---|---|
| `hex_flow_{1,2,4}` | six steps; on each, the pressure solve converges below the iteration cap and the fields stay finite; after the first, the projection removes at least 90 % of the predictor's divergence |
| `hex_invalid_{1,2,4}` | one inverted hex, on one rank; setup refuses it on every rank |
| `hex_flow_ranks_{2,4}` | the kinetic energy after six steps equals the 1-rank value to 1e-6 (`tests/release/check_rank_value.py`) |

The shared NS solver's hex regressions are master's release tests (`ctest -L release`). They cover
the lumped mass, which now uses the exact hex sub-volumes.
