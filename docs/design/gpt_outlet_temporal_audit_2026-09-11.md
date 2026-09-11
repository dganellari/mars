# Public channel: outlet trace feedback amplifies a timestep mode

Author: GPT/Codex. Date: 2026-09-11.
Status: source audit, executed host reconstruction, and completed 200-step public GPU fixed-trace control.
Source snapshot: `03acc1d7f15080ff31f152b5bc0726ab0fbdac6d`.

The current lagged average-pressure trace produces a growing mode in a linear
reconstruction of the complete timestep. Advection is disabled in this model.
Its early velocity growth closely follows the public GPU failure. This identifies
an instability in the time coupling specified in GPT's outlet design; passing
individual flux/Jacobian and short integration gates does not exclude it.

The predicted fixed-trace control now completes 200 steps on the GPU with bounded
reported speed and small full continuity residuals. This supports the diagnosis
of lagged trace feedback in this formulation. It does not validate beta=0.05 or
establish long-time stability of beta=1.

Follow-up: the [fully implicit coupling review](gpt_outlet_coupled_time_review_2026-09-11.md)
finds that refreshing both trace and reconstructed gradient, with an exact
momentum response, still has growing modes at beta=0.05. Removing the lag alone
is not a repair of this boundary closure.

No production solver code was changed in this audit. The reproducible model is
[`scripts/outlet_temporal_audit.py`](../../scripts/outlet_temporal_audit.py).
It assembles the documented algebra independently with NumPy; it does not execute
CUDA kernels, reproduce Hypre iterations, or constitute a full GPU operator-parity test.

## Evidence and provenance

The public fixture is the procedural `[0,4] x [0,1] x [0,1]` channel, with 425 nodes
and 1,536 Kuhn tetrahedra. No application geometry is used by the model.

The user supplied a matched-command record: A enables `--outlet-beta=0.05`; B omits
that flag and remains bounded for 200 steps. B's completion and final reports were
provided in an attachment. The full A log was subsequently retrieved with authorized
`rsync` from `daint-alps3`, then inspected locally:

```
remote: /capstor/scratch/cscs/gandanie/git/mars/daint-gpu/channel-A-outletbeta.log
local:  /private/tmp/mars-outlet-audit-20260911/channel-A-outletbeta.log
SHA256: 7a173734db48886f6562599f5d741914b7b2da1dc56164ddaf1bb3275bf2d5c1
```

The observed A header confirms one rank, skew advection, rho=1000, nu=1e-4,
dt=2e-6, 200 requested steps, and beta=0.05. The supplied matched-command record
also gives inlet=0.5, ramp=100 steps, relax_u=0.3, relax_mass=1,
outlet_rtol=1e-9, outlet_div_tol=1e-8 and outlet_flux_tol=1e-10. The log does not
embed a source commit or the entire command/environment, so those are not
independently established by its header.

The full log has a detailed phase trace only for step 1. Later records show:

| Step | GPU maximum speed | Host maximum speed, advection disabled |
|---:|---:|---:|
| 10 | 0.106 | 0.105788886 |
| 20 | 0.322 | 0.321504628 |
| 30 | 2.537 | 2.53671828 |
| 40 | 16.154 | 16.1475622 |
| 50 | 103.834 | 103.607422 |
| 60 | 787.027 | 773.891490 |
| 70 | 5820.064 | 5548.28613 |
| 80 | 329979.537 | 54339.7537 |

The GPU speed is printed to three decimal places. Agreement through step 50 is
within 0.22%; this is especially useful because the host model omits advection.
The subsequent departure is compatible with nonlinear advection amplifying an
already growing state. It does not show that advection initiated the instability.

At GPU step 40 the full continuity RMS is 1.8e-15 /s while maximum speed is 16.154.
A stricter continuity tolerance alone would therefore not prevent the initial growth.

The prior eight-step 1/2/4-rank checks and the additional empty-opening-rank gate
passed according to the user's complete parser result. Claude reports separate
573-check kernel passes on one and four ranks; GPT has not inspected those transcripts.
These checks remain valid within their scope. The outlet stage remains open for
time stability of the average-pressure feedback formulation.

## Executed public GPU fixed-trace control

The user ran C with `--outlet-beta=1`. GPT retrieved its complete log with
authorized `rsync` and inspected it locally:

```
remote: /capstor/scratch/cscs/gandanie/git/mars/daint-gpu/channel-C-fixedtrace.log
local:  /private/tmp/mars-outlet-audit-20260911/channel-C-fixedtrace.log
SHA256: 174ea09c5dc82a5e23ac5c2e9cdc6c010e245e8db59efad507ed5f327f43c5d8
```

The C header matches A's public fixture, one rank, density, viscosity, timestep,
inlet speed, skew advection, stabilization and relaxation settings, with beta
changed to 1. Neither log embeds the executable revision or complete environment.
The command supplied for this control is retained below.

| Observed quantity | C: fixed trace, beta=1 |
|---|---:|
| Completed physical steps | 200 |
| Maximum reported nodal speed over the 20 ten-step reports | 1.234 |
| Final nodal maximum speed | 1.233 |
| Final full continuity RMS | 3.86e-13 /s |
| Final full continuity maximum | 2.43e-12 /s |
| Final signed boundary imbalance | 4.22e-15 volume/s |
| Largest post-correction RMS across all 200 trial records | 3.1643e-11 /s |
| Final stabilized flux at each of the three interior cuts | 0.5000 volume/s (printed precision) |

All 200 trial records accept the first correction with no backtracking and omega
within 3.9e-13 of 1. All 200 true Krylov relative residuals are at most 9.944e-9.
The local parser found no failed solve, abort or nonfinite-value marker. Full
continuity maxima and boundary sums are printed only every ten steps; their
largest reported magnitudes are 1.582e-10 /s and 2.62e-12 volume/s respectively.
This is log verification of an executed GPU run, not a field-level accuracy test.

The final raw boundary imbalance is 0.007%; stabilized balance prints -0.000%.
The legacy `divRC` diagnostic is not the full acceptance residual for this path.
Beta=1 fixes the separate boundary trace to p_ref while preserving free pressure
DOFs and the new conservative correction. It does not restore the old whole-face
Dirichlet pressure solve, despite the driver's unchanged generic outlet banner.

C removes the rapid growth seen in A, as the independent linear model predicted.
The evidence implicates the lagged trace/correction time coupling rather than a
standalone momentum linear-solver failure. The GPU result covers one rank and
200 steps, only 0.0004 s of physical time. It is a working public baseline for
further development, not proof of steady flow, long-time stability, spatial
accuracy, or validity on other meshes.

## Discrete map audited

The PDE uses physical pressure and kinematic viscosity:

```
du/dt + (u.grad)u - nu*Laplacian(u) + grad(p)/rho = 0.
div(u) = 0, measured through the stabilized control-volume flux.
```

Velocity is prescribed on inlet/wall DOFs and free on the outlet DOFs, with the
driver's corner precedence. All pressure/continuity DOFs remain active. Triangle
opening samples use A_f/3; interior SCS areas use V_e*(grad N_R-grad N_L)/4.
The host model preserves the production velocity mask, gradient sign, boundary
term, projected-gradient masks, face/opposite gradient blend, and nodal coefficient.

Write a=dt_eff, h=a/rho, z=h*p and scale the trace similarly. With fixed geometry
and p_ref=0, the trace map is F=(1-beta)*(E-1*w^T*E). Let G_v and G_t be the
volume-pressure and trace parts of the nodal gradient. Let B be full integrated
continuity, C the compact pressure-flux map divided by h, T the trace-flux map
divided by h, and W the reconstructed-gradient flux map divided by h.

On a fixed BDF stage, the reconstructed step is

```
g_pred_scaled = (G_v + G_t F) z_old
u**           = H [w - Q*g_pred_scaled] + velocity-Dirichlet lift
R             = B u** + (W*(G_v+G_t F) + C + T F) z_old
J             = -B Q G_v + C
delta_z       = -J^-1 R
u_new         = u** - Q G_v delta_z
z_new         = z_old + delta_z
```

Here H=(M/a+nu*K)^-1*(M/a) on free velocity DOFs. For BDF2,
w=(4*u_old-u_older)/3 and a=2*dt/3; BDF1 uses w=u_old and a=dt.
The model rescales z by 2/3 at BDF startup, preserving physical pressure.
The analyzed homogeneous BDF2 state is `(u_old_free,u_older_free,z_old)`;
its amplification matrix has dimension 1,385.

Source correspondence in `backend/distributed/unstructured/fem/`:

| Operation | Current source |
|---|---|
| BDF predictor and density factor | `mars_ns_pump_solver.hpp:2685`, call at `:9056` |
| Previous-pressure gradient including trace | `mars_ns_pump_solver.hpp:8746-8824` |
| Velocity diffusion matrix and lift | `mars_ns_pump_solver.hpp:9195-9356` |
| D=relax_u*V/(rho*a_P), BC fallback, active BDF matrix | `mars_vms_pressure_stab.hpp:59-87`, `mars_ns_pump_solver.hpp:8478-8521` |
| Interior reconstructed-gradient masks and compact flux | `mars_vms_pressure_stab.hpp:166-205` |
| Trace map and once-per-step refresh | `mars_ns_pump_solver.hpp:11260`, `:11736`; driver `mars_pump.cu:1799` |
| Frozen-trace increment gradient | `mars_ns_pump_solver.hpp:10550-10668` |
| Boundary evaluator and all-row residual | `mars_outlet_flux.hpp:57-86`, `mars_ns_pump_solver.hpp:11606-11669` |
| Pressure/velocity update | `mars_outlet_correction.hpp:174-191` |
| History advances once | `mars_ns_pump_solver.hpp:12462-12470`, `:12579-12589` |

The inspected path has matching density and BDF factors. The model's velocity
map is independent of rho after the change to z for these zero-pressure-reference,
zero-body-force conditions. A large physical pressure alone is not a diagnosis.

## Executed host experiments

The model recovers volume 4, 160 free velocity nodes and nu*trace(K)=0.0384 for
the supplied small-nu case. At startup it predicts corrected x-speed
0.00605512910 and pressure maximum 11,630,918.44; the GPU prints 0.0060551268
and 11,630,918. These anchor the model to the observed startup without fitting.

| Host experiment | Spectral radius | Eigenvalues with modulus >1+1e-8 |
|---|---:|---:|
| Current lagged trace, dt=2e-6, nu=1e-4 | 1.2340702024 | 8 |
| Fixed trace, same coefficients | 0.9999999934 | 0 |
| Current lagged trace, dt=0.01, nu=0.1 | 1.1800794657 | 8 |
| Counterfactual implicit trace only, dt=2e-6, nu=1e-4 | 1.5780058863 | 10 |

The matrix action agrees with a separate execution of the model's timestep to
relative error <=4.3e-15. Leading eigenpair residuals are <=1.6e-11, including
the nearly neutral fixed-trace modes. The inverse checks are <=2.7e-15.
These checks verify the host algebra and eigensolves, not production device execution.

The low-nu case reaches maximum speed 2.56e6 in its largest component by host
step 100 even without advection. A fixed trace stays bounded over the same forced
history. The more viscous case also has a growing linear mode: eight successful
steps were too short to exclude it. Making just the trace implicit, including
its derivative in gradient and flux, is not a supported fix; this experiment
makes the spectral radius worse while still freezing the reconstructed gradient.

Reproduce locally with NumPy, no mesh input:

```bash
python3 scripts/outlet_temporal_audit.py --steps=100 --spectrum
python3 scripts/outlet_temporal_audit.py --trace=fixed --steps=100 --spectrum
python3 scripts/outlet_temporal_audit.py --dt=.01 --nu=.1 --rho=1 --steps=0 --spectrum
python3 scripts/outlet_temporal_audit.py --trace=implicit --steps=0 --spectrum
```

## Acceptance and reproduced GPU command

The correction loop uses `max(atol, rtol*initial_residual)` and a relative boundary
balance scale. This explains why, much later in A, a final RMS of 4.08e126 can
satisfy a threshold based on an initial RMS of 4.00e142. Armijo descent is working
as a correction test; it is not a physical-step stability test. A separate
physical acceptance criterion should use a declared reference scale, but adding
one alone would only stop the failure, not repair the unstable map.

The control changed beta from 0.05 to 1 in A, retaining the new boundary-aware
correction. No rebuild or new solver implementation was needed. The supplied
command, from the user's Daint build directory, was:

```bash
set -o pipefail
MARS_SOLVE_TRACE=1 srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=1 --export=ALL --kill-on-bad-exit=1 \
  ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_pump \
  --mesh=../tests/data/public_outlet_channel/outlet_channel.exo --inlet-ss=inlet --outlet-ss=outlet \
  --solver=hypre --vms-stab --rc-implicit --outlet=do-nothing --outlet-beta=1 \
  --rho=1000 --nu=1e-4 --inlet-velocity=0.5 --dt=2e-6 --num-steps=200 --source-ramp-steps=100 \
  --relax-u=0.3 --relax-mass=1 \
  --outlet-max-corrections=200 --outlet-rtol=1e-9 --outlet-div-tol=1e-8 --outlet-flux-tol=1e-10 \
  2>&1 | tee channel-C-fixedtrace.log
```

The prediction that removing trace feedback removes the early growing mode is
supported by C. Keep beta=1 as the public comparison baseline. The next numerical
design step is to derive the trace/reconstructed-gradient/pressure update as a
coupled time map and test its amplification before changing production kernels.
Recommended effort: High, one agent, limited to that derivation and its host gate;
routine implementation and known validation should return to Medium. Do not ship
a guessed damping value or the implicit-trace-only counterfactual as a fix.
