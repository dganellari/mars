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
cmake -S tests/sfem_writer -B build-sfem-writer
cmake --build build-sfem-writer
ctest --test-dir build-sfem-writer --output-on-failure
```

The 27 host/MPI tests compile an unchanged copy of the actual writer header with
small test-only device/domain shims. They cover Tet4/Hex8, both field precisions,
scalar/vector values, numbering, multiple frames and appended times at 1/2/4
ranks; four ranks includes empty-owner ranks. Poisoned ghosts and transfer byte
counts check owned-only field staging. Failure cases cover malformed fields,
rank disagreements, filename collisions, invalid connectivity, duplicate/missing
owners, and injected transfer/write/close errors. Count bounds need no huge
allocation. An unexpected exit, launcher failure or timeout cannot pass a
negative test. These checks exercise MPI and host I/O, not GPU execution.

Enable `MARS_WRITER_GPU_CHECKS=ON` with `CMAKE_PREFIX_PATH` pointing at a GPU
installation of this candidate to add 12 native checks. They use only public
procedural Tet4/Hex8 boxes, poisoned ghosts, scalar/vector fields, both precisions
and two frames at 1/2/4 ranks, plus bad-shape and filesystem failures. Configure
the `MPIEXEC_*` CMake variables for the allocated system. Native CUDA/HIP and
multi-node validation remain separate from the host shim suite.
