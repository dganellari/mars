# Pressure gauge checks

The pump stepper updates the pressure incrementally, p^{n+1} = p^n + φ, with φ the potential its
Poisson problem solves (the standard incremental pressure correction).

**The mean removal.** Since af88fd3d it removes φ's mean after the solve when the opening-flux
source is on, as it does for periodic boxes. On the matrix-free DDT path a single pin controls the
constant mode only weakly, and the opening-flux source excites it. The removal shifts p by a
constant every step, so where the pressure already has a gauge, it drifts off it. A gauge is a
Dirichlet pressure row: a prescribed opening pressure via `outletPressureLift`, or a masked outlet. The velocity is unaffected, since it sees ∇φ only.

`NSStepper::pressureIncrementMean` chooses:

| value | removes φ's mean |
|---|---|
| `Auto` (default) | on the matrix-free DDT path as before; on the assembled paths (Hypre DDT, K) only if no rank has a Dirichlet pressure row (`pressureHasGauge`). A pin alone fixes an arbitrary level, so pin-only runs keep the removal. |
| `Always` | whenever the periodic or opening-flux condition holds (the behaviour before this option) |
| `Never` | only if the pressure has neither Dirichlet rows nor a pin |

The RHS's mean removal on the matrix-free DDT path is unchanged. That is where it stabilises the
solve.

```sh
cmake -S tests/pressure_gauge -B build-pressure-gauge -DCMAKE_PREFIX_PATH=/path/to/candidate/install
cmake --build build-pressure-gauge
ctest --test-dir build-pressure-gauge --output-on-failure
```

The case is the public Tet4 box of `tests/outlet_pressure` with 6³ cubes, inflow and a prescribed
outlet pressure T = 10, with the opening-flux source on. It uses Hypre and the assembled DDT, takes 5
steps, and runs at ranks 1, 2 and 4.

| mode | checks |
|---|---|
| `auto` | the outlet keeps p = T to 1e-9 relative (`pressureHasGauge` is true) |
| `always` | the outlet's pressure moves off T by more than 1% (the old behaviour) |
| `velocity` | Auto and Always give the same velocity to 1e-8 relative; their pressures differ by a constant |
