# Backend layout

The v0.2 source layout separates the two Kokkos users from native Cornerstone AMR.
The public CMake targets and build-option names stay the same.

| Directory | Responsibility | Enable with |
| --- | --- | --- |
| `backend/kokkos/common/` | Shared Kokkos views, mesh interfaces and utilities | `MARS_ENABLE_KOKKOS` |
| `backend/kokkos/bisection/` | Older Kokkos mesh storage and bisection | `MARS_ENABLE_KOKKOS` and `MARS_ENABLE_AMR_BACKEND` |
| `backend/distributed/structured/` | Structured mesh, DOF and SFC implementation | `MARS_ENABLE_KOKKOS` and `MARS_ENABLE_DISTRIBUTED_BACKEND` |
| `backend/distributed/communication/` | Shared distributed communication | `MARS_ENABLE_DISTRIBUTED_BACKEND` |
| `backend/distributed/unstructured/` | Native CUDA/HIP mesh management using Cornerstone | `MARS_ENABLE_UNSTRUCTURED` |
| `backend/distributed/unstructured/amr/` | Native conforming Tet4 and nonconforming Tet4/Hex8 refinement | Part of the unstructured backend |

`MARS_ENABLE_AMR_BACKEND` retains its historical name for compatibility. It enables
only the older Kokkos bisection module; it does not control native Cornerstone AMR.
The serial backend is unchanged. The move does not change support status: see
[Structured Kokkos](Structured-Kokkos.md), [AMR](AMR-Module.md) and
[Known limitations](../KNOWN_LIMITATIONS.md).

## Include compatibility

Link consumers to `Mars::mars` or `Mars::mars_unstructured` as before. Flat header
names such as `mars.hpp` and `mars_distributed_mesh_kokkos.hpp` still work through
the exported include directories. New source-relative includes use the canonical
paths above.

The following old source-relative prefixes remain supported for v0.2:

| Previous prefix | Canonical prefix |
| --- | --- |
| `kokkos/` | `backend/kokkos/common/` |
| `backend/amr/` | `backend/kokkos/bisection/` |
| `backend/distributed/mesh/` | `backend/distributed/structured/mesh/` |
| `backend/distributed/mesh_manager/` | `backend/distributed/structured/mesh_manager/` |
| `backend/distributed/sfc/` | `backend/distributed/structured/sfc/` |

CMake generates forwarding headers for enabled modules in the build tree and
installs them at the previous paths. Each forwarder includes the canonical header;
there is one implementation. Both path forms can be included in the same program.
The forwarding include directory is exported by the build target, so source-tree
consumers should also use the CMake target rather than hand-written include paths.

The structured export includes Kokkos common code and structured/shared distributed
code, while continuing to exclude the older bisection and unstructured backends.
Existing sealed validation snapshots keep their original layout and source hashes.

## Validation

The installed structured consumer includes both old and new path forms, creates a
mesh, and checks owned/ghost field values. Run it from a relocated installation in
both direct and host-staged MPI modes on each supported execution backend. The
optional Kokkos bisection consumer checks its existing terminal-edge witness.
Unstructured installation and AMR checks remain separate from these Kokkos tests.

Directory relocation does not change numerical implementations. Native results
from an older layout remain numerical evidence for those bytes, but the final
package still requires build/install/consumer checks on CPU, CUDA and HIP.
