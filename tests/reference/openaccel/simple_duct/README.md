# Rectangular-duct SIMPLE validation

This benchmark checks the production Tet4 SIMPLE solver against the analytical fully developed
laminar flow in a rectangular duct. It compares the velocity profile and the pressure gradient,
treats the entrance development behind the uniform inlet explicitly, and runs a refinement
study with 1/2/4-rank comparisons. Everything is synthetic and public: `duct_mesh.py` generates
the meshes, and no external case data is required.

This suite is isolated to this directory. It does not change the production solver, halo,
limiter or global CMake configuration.
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

### The discrete fully developed state: an assumption, checked per run

- **Mesh:** hexes of hx = 2h and hy = hz = h, h = H/cells, each split into six positively
  oriented Kuhn tets around the hex diagonal. The mesh is invariant under x → x + hx.
- **What the symmetry gives:** a discrete state that is invariant under that shift, with a
  linear pressure drop, is *compatible* with the discrete equations. If one is reached, its
  stabilized mass flux is ρUWH, because the inlet flux is exact and the scheme conserves mass.
  This does not imply that the P1 velocity section integral equals UWH.
- **What it does not prove:**
  - that such a discrete state exists for this nonlinear SIMPLE discretization;
  - that it is unique;
  - that the computed solution approaches it away from the ends.

  None of this is proven here. The benchmark assumes it. The window indicators test it on
  each run: they are consistent with it, but they cannot prove it.
- **Consequence, under that assumption:** the window's profile measures the discretization of
  the fully developed problem, independent of how the inlet corners are treated.

### Entrance development with the uniform inlet

The plug inflow differs from the developed profile by O(U). At low Re that difference decays
like Stokes eigenmodes.

- **Heuristic decay rate:** use the plane-channel estimate exp(−κx/a), where a is the
  half-width, κ = Re(z₁)/2 = 2.1062 and z₁ = 4.2124 + 2.2507i solves sin z + z = 0.
  This selects a window margin; it is not a decay bound for the rectangular duct.
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
- **A-posteriori indicators** (measured on the discrete solution): every window section must
  equal the reference section to within 10% of the run's own max profile error, with a floor
  of 1e-6·u_c. The section-mean pressure must be linear to within 10% of the run's own |G
  error|. These are indicators, not a bound on contamination: they detect entrance or outlet
  transients that vary across the window, but a transient that is nearly uniform over the
  window would pass them. That case is what the a-priori margin addresses, and the margin
  rests on the plane-channel Stokes estimate, not on a proof for the duct.
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

1. The evidence is complete. The Exodus file named by the mesh description exists next to it
   and matches its SHA-256. The fields cover every lattice node once, with the lattice
   coordinates to 1e-12 and finite values; `-fields.csv` or the distributed parts listed by
   `-fields.json` are both accepted. The metrics and the log exist. A recorded exit status of 0
   exists: `PREFIX.exit` from the GPU recipe, or the manifest of `run_host_study.py`. A
   missing file or a nonzero exit fails; nothing is skipped.
2. The run is the one it is labelled as. The log has exactly one `SIMPLE Tet4, R ranks ...,
   <scheme>, laminar` header and one final line. Both give the labelled rank count R. The
   scheme is the one compared (`--advection`, default upwind). The final line's iteration
   equals the last row of the metrics.
3. The log states ρ, μ and U on its control line, and they are the compared values. A missing
   control line or value fails. The inlet mass flux is ρUWH to 1e-9.
4. The run converged. The log has `CONVERGED`. The final momentum and continuity residuals are
   at most `--residual-tol`, the mass balance at most `--mass-tol`, and du, dp and dflux at most
   `--change-tol` (all 1e-8 in the study). Cancellation is at most 1e-10, and no outlet face is
   closed or changed.
5. The window satisfies the a-priori bounds and holds at least 3 node planes.
6. The window passes both a-posteriori indicators above.

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
- a changed Exodus file;
- missing evidence: no Exodus file, no Exodus name in the description, no log, or a log
  without the control line or without μ;
- a wrong run identity:
  - files labelled 2 ranks whose log is a 1-rank run ending at another iteration;
  - a log ending at another iteration than the metrics;
  - the other advection scheme;
  - two concatenated logs;
  - a log without its header;
- no exit record, a nonzero exit (137) despite complete outputs, or a malformed exit record;
- a missing distributed field part.

Rank differences of 5e-8 pass, distributed field parts pass, and high-resolution runs pass
when compared as high-resolution.

## Files

| File | Purpose |
|---|---|
| `duct_analytic.py`, `duct_analytic.hpp` | Series solution, K, G, fRe, planar control, entrance estimates (Python drives the comparator; C++ cross-check) |
| `duct_mesh.py` | Lattice, Kuhn tets, `inlet`/`outlet`/`walls` side sets. Writes Exodus II (netCDF 64-bit offset, standard library only) and a JSON description with SHA-256 |
| `duct_mesh.hpp` | C++ mirror of the lattice, plus a test slab partition in ElementDomain's shape. Coordinates are `(j*W)/ny - W/2` in both languages: no multiply feeds an add, so an FMA-contracting compiler (GCC's default on aarch64) rounds exactly like Python. Tests still allow 4·eps·max(L,W,H) on coordinates only; topology and side sets must match exactly |
| `duct_host_run.cpp` | Host CPU/MPI run: slab partition, production `simple_partition`, production `DistributedSimpleRunner`. Uses the production option parser, and writes the same CSV files and `CONVERGED` line as `mars_segregated_simple` |
| `duct_host_solver.hpp` | Test-only gathered GMRES(60) with ILU(0), relative true residual 1e-12 (the gates' dense LU does not scale to these sizes) |
| `duct_compare.py` | Comparator: `run` (one result, `--ranks` required) and `study` (all runs, parity, refinement, GCI). Checks each run's identity (ranks, advection, final iteration) and recorded exit status |
| `run_host_study.py`, `test_run_host_study.py` | Meshes, host runs and study in one command (used by ctest). Each run records a manifest: the complete launch command (resolved launcher, rank flag and count, launcher arguments, executable, driver arguments), the SHA-256 of launcher, executable and mesh, the exit status and the output SHA-256. Results are reused only when it matches. A missing launcher or executable, or a nonzero exit, fails |
| `duct_exodus_check.cpp`, `test_duct_exodus.py` | Production native reader (`read_simple_mesh`) on generated files, which must return the C++ lattice (netCDF only) |
| `test_duct_analytic.py`, `analytic_check.cpp`, `test_duct_mesh.py`, `test_duct_compare.py`, `test_run_host_study.py` | Tests described above; Python 3.6 compatible (checked with vermin) |

## Reported host results (source branch: CPU, OpenMPI 4.1, 4 cores)

Test suite (`ctest -LE long`, 7 tests, about 3 minutes). It was built with `-mfma
-ffp-contract=fast` to reproduce the FMA contraction that aarch64 GCC applies by default:

| Test | Result |
|---|---|
| analytic (C++, Python, cross-check) | pass |
| mesh (topology, orientation, tags, netCDF bytes, determinism, C++ mirror) | pass |
| comparator (44 cases: missing evidence, run identity, exit records, distributed parts) | pass |
| production Exodus reader on 1 and 2 ranks | pass |
| study script (13 cases: manifests, launcher and executable identity, stale or edited caches, nonzero exit) | pass |
| coarse duct on 1/2/4 ranks, fresh | pass |

The production reader test is built because netCDF 4.9.2 was installed locally for it. The same
7 tests also pass with this directory copied onto `cstone` 092297d and built against its current
production headers; there, the fresh host runs converge at iteration 3338 on 1, 2 and 4 ranks.

**Coordinate portability, reproduced:** with the same contracting flags, the 1d2037c formula
`-W/2 + j*hy` compiled to an FMA, and its cells = 6 coordinates differed from Python. The
current formula `(j*W)/ny - W/2` compiles without FMA and matches Python bitwise. For the
power-of-two levels (4 to 64) both formulas are exact, so those meshes are unchanged.

Host study with the production runner, upwind, default controls and tolerances 1e-8. The
three-level refinement uses the 4-rank runs, because 16 cells was run on 4 ranks only; at 4
and 8 cells, 1, 2 and 4 ranks agree to round-off. These runs predate exit records. Their exit
codes, recorded at run time by the local driver script (all 0), were transcribed into `.exit`
files, and the result below passes the identity and exit checks.

| cells | ranks | iterations | profile L2 | max interior | G error | wall slip | transverse/U | window indicator |
|---|---|---|---|---|---|---|---|---|
| 4 | 1, 2, 4 | 3338 | 1.2546e-01 | 8.3026e-02 | −2.2024e-01 | 1.9604e-01 | 4.002e-02 | 2.91e-03 |
| 8 | 1, 2, 4 | 3139 | 7.4461e-02 | 6.2211e-02 | −1.3985e-01 | 1.0857e-01 | 1.108e-02 | 1.85e-04 |
| 16 | 4 | 4484 | 4.1965e-02 | 3.8364e-02 | −8.0783e-02 | 5.8314e-02 | 2.971e-03 | 8.15e-05 |

**Host refinement verdict: FAIL**, on one criterion with the thresholds unchanged. The three-level
G order is 0.44, below the required 0.5.

- **G does not pass:** G_h goes 0.13639 → 0.15045 → 0.16079, and the successive differences
  shrink only by 1.36. Richardson with p = 0.44 overshoots to 0.1894 (+8.3%). The error-based
  orders rise, from 0.66 to 0.79, which looks pre-asymptotic rather than like a wrong limit. That
  is an interpretation, not a demonstration: convergence alone does not close refinement.
- **Everything else passes:**
  - all errors decrease monotonically;
  - the finest-pair orders are 0.83 (profile L2) and 0.79 (G), both above 0.7;
  - the profile GCI passes, with order 0.85 and finest RMS error 2.99e-2 inside the band
    3.91e-2 (the extrapolated RMS error is 6.1e-3);
  - each run passes its per-run checks.

The host levels (h = H/4 to H/16) are coarse. Whether G reaches the asymptotic range is the open
question for the Daint levels 8/16/32 and 16/32/64. The iteration counts also grow with
refinement: at 16 cells the momentum residual falls about 100× per 1000 iterations, against
about 900× at 8 cells.

- **Rank parity** (2 and 4 ranks against 1): max |Δu|/U ≤ 2.3e-15 and max |Δp|/(ρU²) ≤ 8e-13 at
  both levels, with identical iteration counts.
- **Entrance decay** (nz = 8): the section deviation falls by about 5 per 0.5 H; plane-channel
  Stokes theory gives 8.2 per 0.5 H, which the finer mesh approaches.
  - development to 1e-3·u_c is complete at x = 2.0, and the 99% centerline at x = 1.25;
  - the outlet influences the flow up to 1.5 upstream at the 1e-3 level;
  - the section indicator in the window is 1.9e-4·u_c, which is 0.2% of the measured max error.
- **High-resolution advection** (nz = 4, 1 rank): converges after 3502 iterations, with profile
  L2 1.23e-1 and G error −2.27e-1. The weak wall treatment dominates both schemes at this
  resolution.

## Integration checks

The seven short tests passed on macOS against `cstone` 092297d, with C++20,
`-Wall -Wextra -Werror -ffp-contract=fast`, and the native netCDF reader enabled.
Fresh coarse host runs converged at iteration 3338 on 1/2/4 ranks and passed field parity.
All Python files also parse with Python 3.6 grammar. The longer 4/8/16 refinement study was
not repeated during integration; the reported failure above remains open. GPU
results reported after integration are recorded below.

## Reported Daint progress (2026-09-28)

The user-provided log for `simple-duct-upwind-Ng67FX/duct-16-4` passes the per-run
checks after 4484 iterations on four GPUs: profile L2 error 4.1965e-2 and
G = 0.16078541 Pa/m, 8.0783% below the analytic value. These agree with the
reported host values to the shown precision. This is not a refinement verdict.

The 32-cell, one-GPU run then failed in its first pressure correction
(job 4877248, nid005534). The diagnostic rerun (job 4883615, nid005690) returned
after 27/2000 GMRES iterations with no Hypre error, but both the reported and
recomputed relative residual were 1.26351e-12, above the wrapper's 1e-12 limit.
The absolute residual was 5.386251e-15. It already satisfied SIMPLE's independent
linear acceptance criterion, `||b-Ax|| <= 1e-13 + 1e-10*||b||`.

The wrapper now uses that same mixed acceptance criterion for SIMPLE, while
Hypre still targets 1e-12. This changes wrapper acceptance, not the final
linear, nonlinear, field-parity or refinement thresholds. Other callers keep
their prior acceptance policy. The real CPU Hypre regression passes with
ASan/UBSan; the 32-cell CUDA rerun and full 8/16/32 study remain pending.
Preserve the completed runs; no successful fine-grid refinement is claimed yet.

## Daint commands

Run on `cstone`. From the configured MARS CUDA/Hypre build directory (for example `mars-v010-check/build-hypre`), build only the existing production
target; no CMake injection is needed:

```bash
git switch cstone && git pull --ff-only &&
test -f ../tests/reference/openaccel/simple_duct/duct_compare.py &&
cmake --build . --target mars_segregated_simple --parallel 4
```

The meshes are generated on the login node; nz = 32 takes 2 s and 28 MB. Each srun's exit
status goes to `duct-<cells>-<np>.exit`, and the comparator rejects a run whose status is
missing or nonzero, even when its fields and log look complete.

Refinement study, levels 8/16/32 on 1/2/4 GPUs of one node (9 runs). `advection=upwind` is the
production default and the gated case. Rerun the identical block with
`advection=high-resolution` for the limited scheme; it lands in its own directory, and the
comparator then requires high-resolution logs:

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
      2>&1 | tee "$duct_run/duct-$cells-$np.log"
    printf '%s\n' "${PIPESTATUS[0]}" > "$duct_run/duct-$cells-$np.exit"
  done
done
python3 "$duct/duct_compare.py" study "$duct_run" --levels 8,16,32 --ranks 1,2,4 \
  --advection "$advection" --report "$duct_run/study.md"
study_status=$?
printf 'Results: %s\n' "$duct_run"
exit "$study_status"
)
```

A pass requires `**PASS**` at the end of `study.md`, which the study also prints. The study only
reads, so it can be rerun on the same directory, for example with `--window` or a subset of
`--levels`/`--ranks`. Each directory holds its own meshes, and they are checked against their
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
    2>&1 | tee "$duct_run/duct-$cells-$np.log"
  printf '%s\n' "${PIPESTATUS[0]}" > "$duct_run/duct-$cells-$np.exit"
done
python3 "$duct/duct_compare.py" study "$duct_run" --levels 16,32,64 --ranks 1,2,4 \
  --advection "$advection" --report "$duct_run/study-64.md"
)
```

`--field-output distributed` (per-rank CSV parts plus `-fields.json`) is also accepted by the
comparator, if gathered output becomes too large.

## Scope and limits

- **Case:** only the production default controls (Re_Dh = 1.33) and the 2:1 section have been
  run. Other aspect ratios (`--width/--height`), lengths, μ, ρ and U are supported end to end.
  The comparator reads the controls it is given and checks them against the run log.
- **Host solver:** the host linear solver is a test oracle. The CUDA path solves with Hypre
  (rtol 1e-12) and checks the true residual (1e-13 + 1e-10·|b|) as usual.
- **Evidence:** rank parity on the host is at round-off. On the GPU, Hypre's preconditioner
  depends on the partition, so differences at the linear-solve tolerance are expected, and
  1e-6 leaves room for them. GPU progress is recorded above; the complete GPU
  refinement and rank-parity study has not passed.
