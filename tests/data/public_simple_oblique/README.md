# Public oblique SIMPLE channel

`channel.exo` is generated from an 8 x 2 x 2 rectangular grid, split into
192 Tet4 elements with 81 nodes. It has no external geometry source.
The box before rotation is [0,4] x [0,1] x [0,1]. The rigid transform is

```
x' = sqrt(1/2) * (x-y) + 2
y' = (x+y)/2 - sqrt(1/2)*z - 1
z' = (x+y)/2 + sqrt(1/2)*z + 3
```

Side sets are `feed` (inlet), `exit` (outlet), and two no-slip wall sets:
`casing` (original y walls) and `cover` (original z walls). Inward inlet
velocity is U times `(sqrt(1/2), 1/2, 1/2)`. This checks normal-based inlet
prescription independently of axis-aligned boundary names or coordinates.

Regenerate to a new filename using the host test target
`mars_simple_write_oblique_fixture`, built by
`tests/reference/openaccel/distributed_simple/CMakeLists.txt` when netCDF
is available. Its only argument is the output filename. Production runs read
the committed Exodus file directly; no generator or Python is needed on Daint.
