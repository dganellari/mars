# Public outlet channel fixture

Author: GPT/Codex. Date: 2026-09-10. Status: generated and read-back validated.

`outlet_channel.exo` is generated solely from the rectangular domain
`[0,4] x [0,1] x [0,1]`, with 16 x 4 x 4 Cartesian cells split into six positive
Kuhn tetrahedra each. It contains 425 nodes and 1,536 tetrahedra. It is public
synthetic geometry, independent of all pump meshes and case data.

The named Exodus side sets are `inlet` (x=0), `outlet` (x=4), and `walls` (the four
remaining sides). The generator verifies the serialized file: positive Jacobians,
total volume 4, exact exterior-face coverage, outward normals, opening areas 1 each,
wall area 16, and zero total boundary area vector.

Regenerate from the repository root with numpy and netCDF4 available:

```bash
python3 scripts/generate_outlet_channel.py --output tests/data/public_outlet_channel/outlet_channel.exo
```

The committed Exodus file is sufficient for Daint runs; regeneration and Python
packages are not required there. See the shared
[integration test instructions](../../../docs/design/gpt_outlet_channel_integration_2026-09-10.md).
