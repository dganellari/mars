# SFEM writer checks

The writer keeps the existing SFEM raw layout: float32 coordinates, int32
connectivity, scalar or three-component nodal fields in float32 or float64,
and one time entry per frame. It numbers nodes by ascending SFC key and writes
each owner's value once. Per-frame fields and one-time coordinates are packed
on the device before host transfer; mesh keys and connectivity are staged once.

Calls are collective on `MPI_COMM_WORLD`. Every rank supplies the same directory,
step, time, precision and ordered field schema. Each component spans all local
nodes. Names must be unique and must not collide after component suffixes are
added. The domain layout must remain fixed; the writer checks domain identity
and node count but cannot detect every in-place repartition. `appendTimes`
requires the same mesh and fresh step numbers; this is visualization output,
not an atomic checkpoint or a restart compatibility check.

This remains a bounded root gather, requiring total node and element counts at
most `INT_MAX` and sufficient root memory. It does not support coincident nodes
that require distinct identities. Count/descriptor/identity failures, MPI errors,
device transfer errors and filesystem errors abort the communicator. A successful
return means root closed all frame files. An aborted frame may leave partial
files and must be discarded.

```sh
cmake -S tests/sfem_writer -B build-sfem-writer \
    -DCMAKE_PREFIX_PATH=/path/to/candidate/install
cmake --build build-sfem-writer
ctest --test-dir build-sfem-writer --output-on-failure
```

Point `CMAKE_PREFIX_PATH` at a CUDA or HIP installation of this candidate.
The 12 checks use the real MARS domain and device storage with public procedural
Tet4/Hex8 boxes, poisoned ghosts, scalar/vector fields, both precisions and two
frames at 1/2/4 ranks, plus bad-shape and filesystem failures. Configure the
`MPIEXEC_*` CMake variables for the allocated system. The runner requires an
explicit success or expected-failure marker; launcher failures and timeouts fail.
