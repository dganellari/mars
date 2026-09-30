# Independent geometry and outlet probes

The straight public channel establishes agreement for its geometry and controls.
The bent-inlet host/MPI gates compare MARS implementations with each other; they
are not an independent OpenAccel check. These four short references separate the
geometry and outlet choices:

| Case | Mesh | Outlet |
| --- | --- | --- |
| straight-average | Original pinned public channel | Average pressure, beta 0.05 |
| straight-static | Same channel | Constant static pressure, beta 1 |
| warped-average | Warped, rotated public channel | Average pressure, beta 0.05 |
| warped-static | Same warped channel | Constant static pressure, beta 1 |

All use high-resolution advection, `linear_linear` velocity interpolation,
rho=1 kg/m³, mu=0.1 Pa s, inward speed 0.1 m/s, outlet pressure 0 Pa, and zero
initial velocity/pressure. Relaxation and pseudo-time remain 0.3/0.3/0.75 and
0.01 s. Each runs exactly 20 outer iterations and saves every state. This is a
matched-iteration test, not a convergence test.

For the public coordinates (x,y,z), the warp is
`(x + 0.2*y*y + 0.1*z*z, y + 0.05*x*x, z)`, followed by the proper rotation
`[[.8,-.6,0],[.48,.64,-.6],[.36,.48,.8]]`. Coefficients multiplying squared
coordinates have units of inverse metres. Connectivity, node IDs and side sets
are unchanged. The runner rejects any inverted or collapsed tetrahedron.

`scripts/openaccel_simple_geometry_probe.py` accepts only the checksum-pinned
public channel and deck from the completed public capture, and requires that
capture's exact executable. It derives the four cases itself; it accepts no
user-supplied replacement geometry or physical parameters. Legacy instrumentation
is disabled, so the changed cases cannot be mistaken for frozen captures of the
original fixture. Each `probe.json` records input, executable and output hashes,
the command, revision, iteration coverage and exit status.

## Daint reference command

Use the existing **OpenAccel environment** that ran the previous public
reference. Keep the MARS build environment in its separate terminal. No build is
needed. The library check must succeed before launching the single-rank reference
job. The Python environment needs numpy and netCDF4.

```bash
(
set -euo pipefail
root=/capstor/scratch/cscs/gandanie/git
git -C "$root/mars-v010-check" pull --ff-only
python3 -c 'import numpy, netCDF4'
exe="$root/OpenAccel-reference-updates-IxgJIp/source/build/openaccel-3D.exe"
links=$(ldd "$exe")
if printf '%s\n' "$links" | grep -q 'not found'; then
  printf '%s\n' "$links"
  echo 'Restore the OpenAccel runtime environment before launching.'
  exit 1
fi
run=$(mktemp -d /capstor/scratch/cscs/gandanie/simple-geometry-reference-XXXXXX)
mkdir "$run/tmp"
export TMPDIR="$run/tmp"
printf 'Public results: %s\n' "$run"
git -C "$root/mars-v010-check" rev-parse HEAD > "$run/mars-revision.txt"
printf '%s\n' "$links" > "$run/reference-libraries.txt"
srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=1 \
  --cpus-per-task=1 --cpu-bind=cores --export=ALL --kill-on-bad-exit=1 \
  env OMP_NUM_THREADS=1 OMP_PROC_BIND=close OMP_PLACES=cores \
  python3 "$root/mars-v010-check/scripts/openaccel_simple_geometry_probe.py" \
  --capture "$root/OpenAccel-reference-updates-IxgJIp/run" \
  --executable "$exe" --output "$run/reference" \
  2>&1 | tee "$run/launch.log"
printf 'Public results: %s\n' "$run"
)
```

The four result directories contain public synthetic data only. Comparing every
saved state against the same MARS kernels locates the earliest discrepancy. A
failure restricted to the static cases points toward outlet semantics; one
restricted to warped cases points toward geometry-sensitive operations. These
are localization clues, not proofs of cause. Passing all four leaves other
controls, longer histories and reference-run provenance to investigate.

## Current evidence

Ten local harness checks pass, including unknown input rejection, pinned binary
identity, independent outlet/geometry controls, launch failure and missing states.
Generation against the actual pinned public mesh preserves every non-coordinate
variable. The minimum positive element determinant is 0.015625 m³ before the warp
and 0.013505859375 m³ afterward. Neither determinant is the element volume, which
is one sixth of it.

All four cases also complete 20 local steps through the production arithmetic in
the host `SimpleRunner` with its test-only dense linear solver. The control case
matches all first 20 states of the existing independent public OpenAccel run:
maximum velocity difference/U is 6.46e-13, and maximum absolute pressure
difference/(rho U²) is 1.04e-10. This checks the local comparison path; it does
not establish agreement for the other three cases.

Independent OpenAccel execution and subsequent field comparisons are pending.
No production kernel, solver tolerance or private input changed. These probes do
not establish GPU execution, multi-rank scaling, pump parity or convergence.
