# MFEM Beam-Tet Validation Example

Solves `-Δu = 1` with `u = 0` on boundary using MFEM's beam-tet.mesh.

## Usage

```bash
# Build
cd cpu && make

# Run (auto-converts MFEM mesh on first run)
mpirun -np 1 bin/mars_ex_beam_tet --mesh beam-tet.mesh
```

That's it! The example automatically:
1. Detects MFEM format
2. Calls `mfem_to_mars_binary` (built by CMake)
3. Converts with refinement
4. Caches result for next time

## Manual Conversion (Optional)

```bash
bin/mfem_to_mars_binary input.mesh output_dir refinement_levels
mpirun -np 1 bin/mars_ex_beam_tet --mesh output_dir
```

## Expected Results

- Converges in ~20-50 iterations
- Residual < 1e-10
- Matches MFEM ex1p exactly
