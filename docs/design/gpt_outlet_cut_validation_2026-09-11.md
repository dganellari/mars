# Signed cut flux and empty-opening coverage

Author: GPT/Codex. Date: 2026-09-11.
Status: diagnostic fix and new public gate implemented; host checks pass.
The revised CUDA kernels and corner-opening fixture still need Daint execution.

## Existing Daint evidence

The user supplied completion markers for eight-step public-channel runs on one,
two, and four ranks. GPT parsed all eight two-rank and four-rank records with the
production log checker. All compared physical moments/fluxes agree between those
two runs within 1.49e-14. The supplied one-rank excerpt contains step 8 and its PASS,
not the full history; its final moments agree with the other partitions to roundoff.

The four-rank run's maximum continuity RMS is 1.553e-13 /s; maximum absolute
boundary imbalance is 2.024e-14 volume/s. All 24 four-rank Jacobian checks pass,
with worst relative error 2.471e-10. BDF startup, history, and halo checks pass.
All supplied configurations report zero empty-opening ranks, so the separate
`--require-empty` coverage condition is not yet satisfied.

These results validate the tested discrete correction and partition invariants.
They do not establish continuum accuracy, long-time stability, or SIMPLE parity.
The true-J implementation is described in
[the Krylov handoff](gpt_outlet_true_j_krylov_2026-09-10.md).

## Confirmed cut-reporter defect

The PDE, nodal unknowns, Tet4 quadrature, frozen coefficients, and Chorin/BDF2
updates are unchanged. This change affects diagnostic flux sums only.

For an interior SCS flux q, production continuity adds +q to row L and -q to row R.
For the node subset `chi_i = 1[x_i <= cut]`, summing its rows gives

```
F_cut = sum_edges (chi_L - chi_R) * q_LR
sum_low_rows R_i = F_cut + sum_boundary_samples_in_low_rows q_boundary
```

The old reporter used `abs(chi_L-chi_R)`: it selected crossing edges but never
reversed the sign when R was the low-coordinate endpoint. Its claim to give an
oriented flux was wrong. This is a graph/control-volume cut, not a geometric
plane-intersection integral.

An independent calculation on the generated public mesh, with constant axial
velocity 0.1 and zero stabilization, reproduces old flux 0.0833333333333333 at
x=1,2,3. The signed sum is 0.1 at every cut: exactly the observed 5/6 factor.
This is a host calculation, not a rerun of the solved GPU fields.

The production reporter now uses `outlet_cut_weight(left,right,cut)` and respects
`VmsFluxCtx::keepSmooth`. Previously it always included the reconstructed gradient,
even when assembly excluded it. No assembly or solve kernel was changed. Reporter
scratch is persistent, and the existing global reduction remains unconditional.

## Tests and new coverage fixture

- 66 host checks execute the production cut-weight helper, including reversed
  edges, endpoint ties, same-side cancellation, and the Kuhn-cube 5/6 regression.
  Clang C++20, strict warnings, ASan and UBSan: PASS.
- Existing 1,435 outlet evaluator/algebra host checks: PASS with the same checks.
- Four log-parser test groups, including wrong cut signs/values, missing records,
  incorrect opening area, and missing empty-rank coverage: PASS.
- CUDA gate adds 180 checks comparing the actual cut kernel against actual
  interior-continuity row sums on one/two-tet fixtures: raw/scalar/nodal flux,
  smooth/compact-only modes, three axes, six cuts. Now 573 checks per rank;
  execution of the new checks is pending.
- The full public channel gate checks cuts x=1,2,3 against both inlet flux and
  independently reduced low-side continuity rows every step. It emits
  `[outlet-cut]` records, with identity tolerance 1e-10 and flux tolerance 1e-8.

`corner_openings.exo` is generated from the same public box and tetrahedra. Each
opening occupies `0 <= y,z <= 0.25` on its end face; the rest is wall. Each opening
has area 0.0625. This concentrates opening facets into corner regions so a four-rank
partition can exercise ranks with no opening facets. Coverage is measured and
required at runtime, not presumed from the mesh or number of ranks.

The generator read-back verifies positive tetrahedra, volume 4, complete exterior
coverage, outward normals, opening areas 0.0625 each, and wall area 17.875. The
original full-opening fixture regenerates byte-identically. The new fixture SHA256:
`d6a89b019562c291e5b13c877e4c285080bcf9a2cfff2a6b411ae8682054c750`.

This corner-opening case is a short algebra/communication test, with the driver's
existing inlet/wall corner precedence. It is not an analytic duct-flow accuracy
test. It retains the public volume and pressure conventions but has a different
source flux from the original full-opening fixture.

## Next Daint runs

Pull `cstone`; rebuild `mars_outlet_kernel_gate` and `mars_pump`. From `daint-gpu/`,
run only the updated four-rank kernel gate and the new two-step coverage case.
There is no need to repeat the original eight-step 1/2/4-rank baselines.

```bash
srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=4 --export=ALL --kill-on-bad-exit=1 \
  ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_outlet_kernel_gate
```

Require PASS (573 checks per rank), then:

```bash
set -o pipefail
MARS_SOLVE_TRACE=1 srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=4 --export=ALL --kill-on-bad-exit=1 \
  ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_pump \
  --mesh=../tests/data/public_outlet_channel/corner_openings.exo --inlet-ss=inlet --outlet-ss=outlet \
  --solver=hypre --vms-stab --rc-implicit --outlet=do-nothing --outlet-beta=0.05 --outlet-channel-check \
  --outlet-channel-opening-width=0.25 --outlet-channel-require-empty \
  --rho=1 --nu=0.1 --inlet-velocity=0.1 --dt=0.01 --num-steps=2 --source-ramp-steps=4 \
  --relax-u=0.3 --relax-mass=1 --tol=1e-11 --max-iter=2000 \
  --outlet-max-corrections=200 --outlet-rtol=1e-9 --outlet-div-tol=1e-8 --outlet-flux-tol=1e-10 \
  2>&1 | tee outlet-channel-empty-4.log

python3 ../scripts/outlet_channel_check.py --require-empty \
  --coverage-log outlet-channel-empty-4.log \
  outlet-channel-1.log outlet-channel-2.log outlet-channel-4.log
```

The parser compares the original three histories and separately validates the new
coverage case, including its smaller source area and cut identities. It does not
compare the different physical fixtures' field moments. If the new partition still
has no empty-opening rank, the driver stops explicitly; do not remove the coverage
requirement to turn that result into a pass.
