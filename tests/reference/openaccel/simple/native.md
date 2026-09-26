# Native ElementDomain input for the public SIMPLE channel

This gate replaces the prepared text topology with the same public Exodus mesh,
loaded through the MARS Exodus reader and partitioned/ordered by `ElementDomain`.
It still runs one rank with the existing fixed channel controls. General meshes,
distributed SIMPLE and changed flow regimes are separate milestones.

The adapter validates a single Tet4 block and the inlet/outlet/walls exterior
coverage, maps original nodes to native SFC indices on-device, verifies bijective
node/element mappings, and resolves boundary face ordinals against native
connectivity. Coordinates retain their original doubles: quantized SFC coordinates
must not change the discrete geometry. The runtime constructs its existing device
CSR from native connectivity. No reference solution enters the solve.

For a node permutation P, the required relation is A_native = P A_prepared P^T
(with the three-component extension for momentum), b_native = P b_prepared, and
x_native = P x_prepared. Element ordering may change summation roundoff. Face
identities and orientation must survive the mapping because they index flux and
reversal history. The final CSV's `node` is the original Exodus storage row; the
comparison applies `node_num_map`, or the Exodus default IDs if absent.

File input and its initial validation are host operations. Native mapping and
connectivity conversion run on-device. Only small setup error checks and the
explicit final output cross back to the host. The reader is restricted to the
public profile; this is not a new general-purpose Exodus reader.

## Build and direct launch on Daint

From the existing MARS CUDA/Hypre build directory, preserve its current cache:

```bash
git pull --ff-only && cmake -S .. -B . && cmake --build . --parallel 4
./examples/distributed/unstructured/mars_segregated_simple_algebra_check
```

Use the existing converged OpenAccel reference; no OpenAccel build or run is needed.
Set `reference_run` to its actual location. Do not assume `$SCRATCH` identifies the
filesystem containing an existing checkout.

```bash
reference_run=/path/to/OpenAccel-simple-converged-20260925-151406
native_run=$(mktemp -d "$PWD/simple-native-XXXXXX")
cp "$reference_run/channel.exo" "$native_run/channel.exo"
printf 'Native results: %s\n' "$native_run"
```

With the working MARS CUDA environment active:

```bash
set -o pipefail
srun --account=csstaff --time=00:20:00 --nodes=1 --ntasks-per-node=1 \
  --export=ALL --kill-on-bad-exit=1 \
  ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_segregated_simple \
  --mesh "$native_run/channel.exo" --mesh-format exodus \
  --output-prefix "$native_run/channel" --iterations 2000 --report-every 20 \
  --residual-tol 1e-6 --mass-tol 1e-6 --change-tol 1e-6 \
  2>&1 | tee "$native_run/run.log"
```

The run must pass the native bijection/exterior-face checks and all nonlinear
criteria. Then use a Python environment with numpy/netCDF4 on the login node:

```bash
python3 ../scripts/openaccel_simple_convergence.py compare \
  --reference "$reference_run" --mars "$native_run" \
  --native-mesh "$native_run/channel.exo" --output "$native_run/comparison.json"
git rev-parse HEAD > "$native_run/mars-revision.txt"
sha256sum ./examples/distributed/unstructured/mars_segregated_simple \
  "$native_run/channel.exo" > "$native_run/sha256.txt"
```

The comparator requires the pinned public mesh hash, the native input path in the
run log, complete node coverage, original convergence gates and the same final
velocity/absolute-pressure tolerances as the prepared path. It does not subtract a
pressure mean. CPU mapping/comparator tests do not establish native CUDA execution;
the native run and final comparison above are the acceptance gate.
