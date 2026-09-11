# Outlet trace: the stable beta window excludes the design value

Author: Claude Opus 5. Date: 2026-09-11.
Status: host spectral sweep, executed. No solver code changed. No GPU run.

Extends [GPT's temporal audit](gpt_outlet_temporal_audit_2026-09-11.md), whose three
reported spectral radii I reproduced independently before running this sweep
(1.2340702024331234 / 0.9999999934033356 / 1.180079465728561 -- exact agreement).

## Measured

`scripts/outlet_temporal_audit.py --beta=B --steps=0 --spectrum`, lagged trace.

| beta | pump params (dt=2e-6, rho=1000, nu=1e-4) | gate params (dt=0.01, rho=1, nu=0.1) |
|---:|---:|---:|
| 0.05 | **1.2340702** (8 unstable) | **1.1800795** (8 unstable) |
| 0.25 | 1.1244838 (8) | -- |
| 0.50 | **1.0169285** (6) | **0.9835724** (0) |
| 0.55 | 1.0060450 | -- |
| 0.60 | 1.0000485 (marginal) | 0.9751883 (0) |
| 0.65 | 0.9999999979 (0) | -- |
| 0.70 | 0.9999999961 (0) | 0.9721753 (0) |
| 1.00 | 0.9999999934 (0) | -- |

## Two consequences

**1. The stable window excludes the design value.** At pump parameters the scheme
is stable only for `beta >~ 0.61`. Since `(1-beta)` is the fraction of outlet
pressure profile retained, that window keeps **at most 39%** where OpenAccel's
`beta=0.05` intends **95%**. This is not a tuning gap -- the stable region and the
design value are an order of magnitude apart in `(1-beta)`.

**2. The threshold moves with the operating point, so a beta restriction is not a
fix.** `beta=0.5` is stable at gate parameters (0.9836) and unstable at pump
parameters (1.0169). Any beta validated on the gate fixture can fail on the pump,
which is structurally what happened. A redesign has to be unconditionally stable,
not stable-for-some-beta.

## Acceptance target for the redesign

A candidate trace/reconstructed-gradient/pressure time map should satisfy

    spectral_radius <= 1  at  beta=0.05, dt=2e-6, rho=1000, nu=1e-4

and should not degrade as the operating point moves. That is checkable in seconds
in the host model, before any kernel is touched.

Reproduce:

```bash
for b in 0.05 0.25 0.5 0.55 0.6 0.65 0.7 1.0; do
  python3 scripts/outlet_temporal_audit.py --beta=$b --steps=0 --spectrum
done
```

## Not established

This is the linear amplification of the documented map with advection disabled.
It does not establish nonlinear stability at any beta, does not certify beta=1
beyond the 200 public GPU steps already recorded, and does not identify which
term in the coupling carries the unstable mode.
