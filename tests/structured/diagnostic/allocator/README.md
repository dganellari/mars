# CUDA allocation and MPI registration probe

This standalone executable contains no MARS headers, mesh code or communication
plans. It compares host memory, `cudaMalloc`, `cudaMallocAsync`, and a
`Kokkos::View<unsigned long*, Kokkos::CudaSpace>`. Raw CUDA modes do not initialize
Kokkos. The two View builds use the original and corrected Kokkos installations.

Every mode uses the same graph communicator, `MPI_UNSIGNED_LONG`, receive-first
posting, tag zero, peer segments and nine checked exchanges. Producers finish
before MPI, and buffers remain alive through `MPI_Waitall`. Padding must remain
unchanged. CUDA pointer and memory-pool attributes are recorded, rather than
inferring allocation type from the numerical pointer value.

The runner performs 60 independent launches: two/four ranks on one/two nodes,
three allocation/message shapes and five allocator/build combinations. Shapes
include 32 KiB, 48 KiB with an interior offset, and the observed failing 64 KiB
allocation with 2,048-word messages. These straddle the pinned Kokkos 4.4.01
40 kB allocator threshold. Raw CUDA pointers have the same 128-byte displacement
as that version's managed View header. Each failed launch is retained; later
allocators still run. Steps are bounded to two minutes.

Ordinary CUDA passing while pooled CUDA and pooled Views fail, together with
the corrected View passing, supports an allocator/transport incompatibility.
Failure of ordinary CUDA, absence of the expected pool attributes, or all modes
passing leaves that diagnosis unresolved. A successful probe does not validate
MARS; run the full structured correctness, install and performance gates with
both reference and candidate rebuilt against the same corrected Kokkos.

The existing profiled setup timings are diagnostic only. The full performance
gate remains uninstrumented with its original workloads and limits. Historical
failed jobs are retained even when a later run passes.
