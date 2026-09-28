# Rectangular-duct SIMPLE validation

This benchmark checks the production Tet4 SIMPLE solver against the analytical fully developed
laminar flow in a rectangular duct. It compares the velocity profile and the pressure gradient,
treats the entrance development behind the uniform inlet explicitly, and runs a refinement
study with 1/2/4-rank comparisons. Everything is synthetic and public: `duct_mesh.py` generates
the meshes, and no input data or case files are read.

Only this directory is new. No production solver, halo, limiter or global CMake code changes.
The CUDA runs use the unchanged `mars_segregated_simple` executable; the host tests use the
unchanged production `simple_partition` and `DistributedSimpleRunner` (host build).

## Mathematical contract

### Problem, as the production code poses it

- **Domain:** Ω = (0, L) × (−W/2, W/2) × (−H/2, H/2), with L = 7, W = 2, H = 1 by default.
- **Inlet** (side set `inlet`, x = 0): uniform velocity U along the inward normal, so the mass
  flux is exactly ρUWH.
- **Walls** (side set `walls`): no-slip.
- **Outlet** (side set `outlet`, x = L): the area-mean pressure is held at `--outlet-pressure`,
  and the tangential viscous traction is kept.
- **Equations:** steady, incompressible, laminar, with the production defaults ρ = 1, μ = 0.1,
  U = 0.1, so Re_Dh = ρUD_h/μ = 1.333 with D_h = 2WH/(W+H) = 4/3.

### Fully developed solution: the Boussinesq series, not the planar parabola

Far from both ends, u = (u(y,z), 0, 0) and p = p₀ − Gx, where

    μ (u_yy + u_zz) = −G,   u = 0 on the four walls,   ∫u dA = U W H.

Let a be the shorter half-side and b the longer one. Let s be the coordinate across the short
side and t the one along the long side, and let λᵢ = iπ/(2a) for odd i:

    u = G/(2μ) (a² − s²) − 16a²G/(μπ³) Σᵢ (−1)^((i−1)/2) i⁻³ cos(λᵢ s) cosh(λᵢ t)/cosh(λᵢ b)
    U = a² G K/(3μ),   K = 1 − 192a/(π⁵ b) Σᵢ tanh(iπb/(2a))/i⁵,   G = 3μU/(a² K)

For W:H = 2:1 this gives:

| Quantity | Value |
|---|---|
| K | 0.686045031 |
| G | 0.174915632 Pa/m |
| u_c/U | 1.991796 |
| Fanning fRe | 15.548056 (Shah & London tabulate 15.54806) |

The first term of the series is the planar Poiseuille parabola. For a finite duct the parabola
is wrong everywhere:

- it gives G = 12μU/H² = 0.12, 31% low;
- it gives 1.5U on the side walls, where the duct solution is zero.

It appears only as a negative control that the comparator must reject.

**Numerics** (`duct_analytic.py`; `duct_analytic.hpp` is an independent C++ implementation):

- The series terms decay like exp(−λᵢ(b − |t|)). The equivalent expansion with the sides
  exchanged decays like exp(−iπ(a − |s|)/(2b)), and each point uses the faster of the two.
- Cosh ratios use a form that cannot overflow.
- Summation stops below 1e-18 of the leading scale.

**Checks** (`test_duct_analytic.py`, `analytic_check.cpp`), none of which use the series
itself as a reference:

- a 4th-order finite-difference residual of μΔu + G below 1e-6·G;
- no-slip in the limit at 1e-9·a from each wall;
- exact symmetry, and invariance under rotating the duct;
- the two expansions agree to 1e-13;
- Gauss quadrature of u reproduces the closed-form flow rate to 1e-9 (measured 1e-14);
- the literature values: fRe for 1:1, 1:2 and 1:4, and u_max/u_mean for the square;
- the parallel-plate limit;
- the Python and C++ implementations agree to 1e-13.

### Why a discrete fully developed state exists

- **Mesh:** hexes of hx = 2h and hy = hz = h, h = H/cells, each split into six positively
  oriented Kuhn tets around the hex diagonal. The mesh is invariant under x → x + hx, so the
  discrete equations have an exactly x-invariant solution. Its flow rate is exactly ρUWH,
  because the inlet flux is exact and the scheme conserves mass.
- **Consequence:** away from the ends, the discrete solution converges to this state. The
  window's profile measures the discretization of the fully developed problem, independent of
  how the inlet corners are treated.

### Entrance development with the uniform inlet

The plug inflow differs from the developed profile by O(U). At low Re that difference decays
like Stokes eigenmodes.

- **Decay rate:** for the plane channel of half-width a, the slowest symmetric mode decays like
  exp(−κx/a), with κ = Re(z₁)/2 = 2.1062 and z₁ = 4.2124 + 2.2507i the first root of
  sin z + z = 0. The mode is complex, so the decay of the max-norm is not monotone (the host
  runs show this).
- **Durst estimate:** the Durst et al. correlation with D = D_h gives L_e = 0.843 at
  Re_Dh = 1.33. Its advective share is only about 0.1.
- **A-priori window:**
  - start ≥ max(a·ln(10⁴)/κ, 2·L_Durst) = max(2.186, 1.686);
  - end ≤ L − a·ln(10⁴)/κ, because the outlet also disturbs the flow upstream;
  - the bounds are rounded inwards to multiples of H/2, so every level shares the window
    [2.5, 4.5];
  - the reference section is x = 3.5.

  A window inside this region fails. Higher-Re runs must lengthen the duct (`--length`) or pass
  `--window`, and the Durst term then moves the start.
- **A-posteriori check** (measured on the discrete solution, not assumed): every window
  section must equal the reference section to within 10% of the run's own max profile error,
  with a floor of 1e-6·u_c. The section-mean pressure must be linear to within 10% of the
  run's own |G error|, so development contaminates the measured error by less than 10%.
- **Reported, not gated:** the development length (profile within 1e-2 or 1e-3·u_c,
  centerline within 1%) and the outlet influence length.

### Expected discrete behaviour (production code as it is)

- **Upwind advection** is first order. On the Kuhn diagonals it adds cross-stream diffusion of
  relative size ρUh/(2μ) = h/2 at the defaults; the first-order advective terms cancel by
  translation invariance.
- **The wall is weak.** The no-slip wall adds a tangential shear 4μ|A_s|²/V·u_node, where A_s
  is the sample area and V the tet volume (`mars_segregated_simple.hpp:75`, active after the
  first outer iteration). The wall node therefore carries roughly the velocity found about h/4
  inside the fluid. This is O(h) slip, and the effective section is larger, so G_h < G.
- **Expectation:** first order for the profile and for G, with G underestimated. The orders are
  pre-asymptotic (below 1) on coarse levels.

### Error measures (reference section x = 3.5 unless stated)

- **E₂:** ‖u_h − u‖/‖u‖ on the section's Kuhn triangulation, with u_h in P1 and the Dunavant
  degree-4 rule.
- **Max errors:** E∞ = max over nodes of |u_h − u|/u_c, and E∞,int over nodes off the walls.
- **Other velocity measures:** wall slip = max |u_h|/u_c on wall nodes; transverse
  = max √(v² + w²)/U; flow ratio = the P1 section mean of u_h divided by U (reported only,
  because it mixes interpolation and slip terms of opposite sign).
- **Pressure gradient:** G_h = −(least-squares slope of the section-mean pressure over the
  window planes), and e_G = (G_h − G)/G.

### Acceptance

**Per run** (`duct_compare.py run`, and every run of a study):

1. Fields cover every lattice node once, with the lattice coordinates to 1e-12 and finite
   values. The mesh description matches the Exodus SHA-256.
2. The run used the compared ρ, μ and U (read from the log's control line), and the inlet mass
   flux is ρUWH to 1e-9.
3. The run converged. The log has `CONVERGED`. The final momentum and continuity residuals are
   at most `--residual-tol`, the mass balance at most `--mass-tol`, and du, dp and dflux at most
   `--change-tol` (all 1e-8 in the study). Cancellation is at most 1e-10, and no outlet face is
   closed or changed.
4. The window satisfies the a-priori bounds and holds at least 3 node planes.
5. The window passes the a-posteriori test above.

**Rank parity:** at every level, the 2- and 4-rank fields equal the 1-rank fields node by node:
max |Δu|/U ≤ 1e-6 and max |Δp|/(ρU²) ≤ 1e-6, with absolute pressure and no mean removed.

**Refinement** (at least three levels, each doubling `cells`; fewest-rank runs):

1. **Monotone:** E₂, E∞,int, |e_G|, wall slip and transverse velocity all decrease at every
   step.
2. **Order:** the observed order of E₂ and of |e_G| on the finest pair is at least 0.7 (formal
   order 1, with some pre-asymptotic margin).
3. **GCI** (Roache), which uses only the solutions:
   - **Order:** the three finest levels give an observed order p in [0.5, 3] and a band
     1.25·|f − m|/(2^p − 1).
   - **G:** the analytic G must lie inside the band around the finest G_h.
   - **Profile:** the same holds for the section profile, where p comes from
     RMS(u_m − u_c)/RMS(u_f − u_m) on the coarse nodes and the band is compared with the finest
     RMS error on the medium nodes.

   This rejects convergence to a wrong limit, such as a biased G or a flow-rate floor, once it
   exceeds the discretization uncertainty. Monotonicity rejects non-convergence, such as the
   planar parabola.

`test_duct_compare.py` exercises every criterion with synthetic fields (analytic plus c·h^p
terms, entrance and outlet transients, written in the production CSV formats). The clean first-
and second-order families pass. Each of the following defects must fail, and each fails for its
stated reason:

- the planar parabola;
- convergence to 1.04·G at orders 0.87 and 0.78;
- G oscillating around the analytic value;
- a 2% flow-rate floor;
- order 0.3;
- stagnant transverse velocity;
- an entrance transient inside the window, and an outlet transient inside the window;
- a window inside the a-priori region;
- a 1e-5 velocity or pressure mismatch between rank counts;
- a missing, duplicate or NaN node, or another mesh;
- unconverged metrics, or a `NOT CONVERGED` log;
- outlet reversal;
- another μ, or another inflow;
- two levels only;
- a changed Exodus file.

Rank differences of 5e-8 pass.

## Files

| File | Purpose |
|---|---|
| `duct_analytic.py`, `duct_analytic.hpp` | Series solution, K, G, fRe, planar control, entrance estimates (Python drives the comparator; C++ cross-check) |
| `duct_mesh.py` | Lattice, Kuhn tets, `inlet`/`outlet`/`walls` side sets. Writes Exodus II (netCDF 64-bit offset, standard library only) and a JSON description with SHA-256 |
| `duct_mesh.hpp` | C++ mirror of the lattice, plus a test slab partition in ElementDomain's shape |
| `duct_host_run.cpp` | Host CPU/MPI run: slab partition, production `simple_partition`, production `DistributedSimpleRunner`. Uses the production option parser, and writes the same CSV files and `CONVERGED` line as `mars_segregated_simple` |
| `duct_host_solver.hpp` | Test-only gathered GMRES(60) with ILU(0), relative true residual 1e-12 (the gates' dense LU does not scale to these sizes) |
| `duct_compare.py` | Comparator: `run` (one result) and `study` (all runs, parity, refinement, GCI) |
| `run_host_study.py` | Meshes, host runs and study in one command (used by ctest) |
| `duct_exodus_check.cpp`, `test_duct_exodus.py` | Production native reader (`read_simple_mesh`) on generated files, which must return the C++ lattice (netCDF only) |
| `test_duct_analytic.py`, `analytic_check.cpp`, `test_duct_mesh.py`, `test_duct_compare.py` | Tests described above; Python 3.6 compatible (vermin: minimum 3.3) |

## Local results (executed on this branch: host CPU, OpenMPI 4.1, 4 cores)

Test suite (`ctest -LE long`, 6 tests, about 3 minutes):

| Test | Result |
|---|---|
| analytic (C++, Python, cross-check) | pass |
| mesh (topology, orientation, tags, netCDF bytes, determinism, identical C++ mirror) | pass |
| comparator (27 cases) | pass |
| production Exodus reader on 1 and 2 ranks | pass |
| coarse duct on 1/2/4 ranks | pass |

The production reader test is built because netCDF 4.9.2 was installed locally for it.

Host study with the production runner, upwind, default controls and tolerances 1e-8. The
refinement to 16 cells is the `long` ctest.

| cells | ranks | iterations | profile L2 | max interior | G error | wall slip | transverse/U | contamination |
|---|---|---|---|---|---|---|---|---|
| 4 | 1, 2, 4 | 3338 | 1.2546e-01 | 8.3026e-02 | −2.2024e-01 | 1.9604e-01 | 4.002e-02 | 2.91e-03 |
| 8 | 1, 2, 4 | 3139 | 7.4461e-02 | 6.2211e-02 | −1.3985e-01 | 1.0857e-01 | 1.108e-02 | 1.85e-04 |

HOST_REFINEMENT_PLACEHOLDER

- **Rank parity** (2 and 4 ranks against 1): max |Δu|/U ≤ 2.3e-15 and max |Δp|/(ρU²) ≤ 8e-13 at
  both levels, with identical iteration counts.
- **Entrance decay** (nz = 8): the section deviation falls by about 5 per 0.5 H; plane-channel
  Stokes theory gives 8.2 per 0.5 H, which the finer mesh approaches.
  - development to 1e-3·u_c is complete at x = 2.0, and the 99% centerline at x = 1.25;
  - the outlet influences the flow up to 1.5 upstream at the 1e-3 level;
  - the window contamination is 1.9e-4·u_c, which is 0.2% of the measured max error.
- **High-resolution advection** (nz = 4, 1 rank): converges after 3502 iterations, with profile
  L2 1.23e-1 and G error −2.27e-1. The weak wall treatment dominates both schemes at this
  resolution.

## Daint (not executed here)

From the configured MARS CUDA/Hypre build directory (for example `mars-v010-check/build-hypre`),
after checking out `cstone-simple-duct`. Only the existing production target is built; no CMake
injection is needed. The meshes are generated on the login node (nz = 32 takes 2 s and 28 MB).

```bash
git fetch origin cstone-simple-duct && git checkout cstone-simple-duct &&
cmake --build . --target mars_segregated_simple --parallel 4
```

Refinement study, levels 8/16/32 on 1/2/4 GPUs of one node (9 runs). `advection=upwind` is
the production default and the gated case. Rerun the identical block with
`advection=high-resolution` for the limited scheme, which lands in its own directory:

```bash
(
set -uo pipefail
unset MARS_OWNERSHIP MARS_HALO_FACTOR MARS_NODEHALO_ALLOW_INCONSISTENT
advection=upwind
duct=../tests/reference/openaccel/simple_duct
duct_run=$(mktemp -d /capstor/scratch/cscs/gandanie/simple-duct-$advection-XXXXXX)
printf 'Results: %s\n' "$duct_run"
git rev-parse HEAD > "$duct_run/mars-revision.txt"
sha256sum ./examples/distributed/unstructured/mars_segregated_simple > "$duct_run/sha256.txt"
for cells in 8 16 32; do
  python3 "$duct/duct_mesh.py" --cells "$cells" --output "$duct_run/duct-$cells.exo" || exit 1
  for np in 1 2 4; do
    srun --account=csstaff --time=01:00:00 --nodes=1 --ntasks-per-node="$np" \
      --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
      ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_segregated_simple \
      --mesh "$duct_run/duct-$cells.exo" --inlet-ss inlet --outlet-ss outlet --wall-ss walls \
      --advection "$advection" --rho 1 --mu 0.1 --inlet-velocity 0.1 --outlet-pressure 0 \
      --iterations 30000 --report-every 500 \
      --residual-tol 1e-8 --mass-tol 1e-8 --change-tol 1e-8 \
      --output-prefix "$duct_run/duct-$cells-$np" \
      2>&1 | tee "$duct_run/duct-$cells-$np.log" || echo "duct-$cells-$np exit $?" >> "$duct_run/run-failures.txt"
  done
done
python3 "$duct/duct_compare.py" study "$duct_run" --levels 8,16,32 --ranks 1,2,4 \
  --report "$duct_run/study.md"
printf 'Results: %s\n' "$duct_run"
)
```

A pass requires `**PASS**` at the end of `study.md`, which the study also prints. A
`run-failures.txt` file means a run did not converge or aborted, and its log says why. The study
only reads, so it can be rerun on the same directory, for example with `--window` or a subset of
`--levels`/`--ranks`. The meshes are regenerated per directory and checked against their
SHA-256.

Optional finest level, to check the asymptotic range (mesh 225 MB, 1.9M nodes, 11M Tet4). Pass
the directory that the block above printed:

```bash
(
set -uo pipefail
unset MARS_OWNERSHIP MARS_HALO_FACTOR MARS_NODEHALO_ALLOW_INCONSISTENT
duct_run=/capstor/scratch/cscs/gandanie/simple-duct-upwind-REPLACE   # "Results:" line above
advection=upwind                                                     # as in that run
duct=../tests/reference/openaccel/simple_duct
cells=64
python3 "$duct/duct_mesh.py" --cells "$cells" --output "$duct_run/duct-$cells.exo" || exit 1
for np in 1 2 4; do
  srun --account=csstaff --time=03:00:00 --nodes=1 --ntasks-per-node="$np" \
    --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
    ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_segregated_simple \
    --mesh "$duct_run/duct-$cells.exo" --inlet-ss inlet --outlet-ss outlet --wall-ss walls \
    --advection "$advection" --rho 1 --mu 0.1 --inlet-velocity 0.1 --outlet-pressure 0 \
    --iterations 30000 --report-every 500 \
    --residual-tol 1e-8 --mass-tol 1e-8 --change-tol 1e-8 \
    --output-prefix "$duct_run/duct-$cells-$np" \
    2>&1 | tee "$duct_run/duct-$cells-$np.log" || echo "duct-$cells-$np exit $?" >> "$duct_run/run-failures.txt"
done
python3 "$duct/duct_compare.py" study "$duct_run" --levels 16,32,64 --ranks 1,2,4 --report "$duct_run/study-64.md"
)
```

## Scope and limits

- **Case:** only the production default controls (Re_Dh = 1.33) and the 2:1 section have been
  run. Other aspect ratios (`--width/--height`), lengths, μ, ρ and U are supported end to end.
  The comparator reads the controls it is given and checks them against the run log.
- **Host solver:** the host linear solver is a test oracle. The CUDA path solves with Hypre
  (rtol 1e-12) and checks the true residual (1e-13 + 1e-10·|b|) as usual.
- **Evidence:** rank parity on the host is at round-off. On the GPU, Hypre's preconditioner
  depends on the partition, so differences at the linear-solve tolerance are expected, and
  1e-6 leaves room for them. Nothing here has run on a GPU.
