# First public-channel failure: line search and volume reference

Author: GPT/Codex. Date: 2026-09-10.
Status: confirmed code defects corrected; host regression checks pass;
subsequent Daint execution confirms approximate-correction stagnation.
See [the true-J follow-up](gpt_outlet_true_j_krylov_2026-09-10.md) for the current implementation.

## Evidence and limits

The user confirmed the supplied log came from the generated public channel fixture.
In Daint job 4641413, step 1 completed its momentum solve and then rejected all
16 pressure-correction backtracks. No full channel step passed. The old log does
not contain the attempted damping, directional derivative, or residual ratios, so
it cannot establish which correction failed or whether the defect below caused it.
No confidential pump geometry or derived output was inspected.

## Confirmed acceptance-rule defect

For the existing frozen residual, `R(omega) = R + omega*delta_R`. Define

```
Vtot = sum_owned V_i
F(omega) = (1/(2*Vtot)) * sum_owned R_i(omega)^2 / V_i
F'(0) = slope = (1/Vtot) * sum_owned R_i*delta_R_i / V_i
```

`F = rms^2/2` has units `1/s^2`. With `c=1e-4`, Armijo requires
`F(omega) <= F(0) + c*omega*F'(0)` and measured decrease. The old helper instead
required `rms_after <= rms_before*(1-c*omega)`, which assumes a sufficiently steep
direction independently of its measured derivative. Backtracking cannot repair
that assumption for arbitrarily weak descent directions.

An independent two-component counterexample is `R=(1,0)`,
`delta_R=(-1e-6,1e-3)`, unit weights. Its optimal damping is approximately 0.999999.
The old predicate rejects that optimum and all 15 subsequent halvings, despite
strict decrease. The slope-based predicate accepts all 16. It also preserves its
decision when the residual units are scaled by `1e-100` or `1e100`.

`outlet_correction_contracts` now consumes the measured slope and tests the Armijo
condition on the squared norm. It requires finite data, a negative slope, and a
strictly smaller measured norm. Differences are evaluated in normalized, factored
form to avoid subtracting nearly equal squares or squaring extreme norms.
The PDE, Jacobian approximation, optimal-damping formula, positivity requirement,
correction cap, and all final continuity/boundary-balance tolerances are unchanged.

With `MARS_SOLVE_TRACE=1`, `[outlet-trial]` records the correction/backtrack index,
damping, normalized slope, quadratic predicted residual ratio, measured ratio,
and old/new acceptance decisions. Agreement of predicted/measured ratios tests
the frozen affine model; disagreement motivates an assembly/state audit. Small
decrease or non-descent after this fix can still require a true-J Krylov method.
No such algorithm has been added or presumed necessary from this log alone.

## Confirmed volume-reference defect

The reported lumped volume was 4.0000005722048915. Tetrahedra use SFC-decoded
coordinates: original-coordinate preservation is wired for hexes, not tets.
The public box is padded to `[-0.2,4.2] x [-0.05,1.05]^2` for encoding. For 64-bit
keys there are 21 coordinate bits; decoding uses denominator `2^21-1`.
The decoded extents are 4.000000190734955 and 1.0000000476837387 twice. Their product
matches the reported volume exactly to printed precision.

The channel check now encodes its two known public corner points with the production
`cstone::sfc3D` and decodes them through the domain's production coordinate API.
Their extent product is the expected lumped volume. The absolute volume tolerance
remains `1e-10`; it was not enlarged. Opening areas still originate from the original
Exodus triangles and remain exactly 1 for this fixture, so the inlet-flux, outlet-area,
trace-mean, and continuity tolerances remain unchanged. General coordinate-storage
changes are outside this fix.

## Executed checks

- 1,435 host production-helper/algebra checks pass with `-Wall -Wextra -Werror`
  and AddressSanitizer/UndefinedBehaviorSanitizer, including the new Armijo cases.
- A standalone C++ check using the actual Cornerstone encoder and the decoder body
  extracted from `domain.cu` reproduced volume 4.0000005722048915. It ran on the CPU.
- `git diff --check` passes.

The subsequent public one-rank run (job 4642390) exercised this repair: corrections
5–7 passed with `old_gate=0`, but correction 8 stagnated at roundoff. Predicted and
measured ratios agreed. The follow-up above replaces the approximate single-direction
iteration with a true-J Krylov solve; its GPU execution is pending.
