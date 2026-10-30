# AMD graph assembly (HIP `tensor`)

Status: the HIP `tensor` graph assembly since 053620d9. It runs `tensor_faces` with 16 threads per
element and 32 elements per block, and adds up each block's contributions to shared rows in shared
memory before the global atomics: 1.35e8 atomics instead of 2.34e8 on hex192. CUDA keeps its
one-thread-per-element `tensor` kernel, because on H100 the aggregation doubles the time.

## Measured

Public hex192, one rank, `mars_cvfem_graph --iterations=20`, the `Assembler` line, run time,
5 alternating rounds (2026-10-30). MI300A on Beverin (ROCm 7.2, SPX mode), H100 on Daint (CUDA 12.9).

| GPU | Kernel | ms |
| --- | --- | ---: |
| MI300A | `tensor` (= `tensor_faces_agg512`) | 8.13-8.28 |
| H100 | `tensor` | 6.60-7.37 |

On gfx942 the kernel uses 94 VGPRs, no spills and 28.9 KB of shared memory per block.

## Where the MI300A time goes

A timing-only ablation (not kept in the tree) replaced the global updates of the same kernel. Only
the atomic update assembles the right matrix.

| Global update | MI300A ms |
| --- | ---: |
| atomic (shipped) | 8.13-8.19 |
| atomic, workgroup scope | 8.13-8.16 |
| plain `+=` (races) | 7.30-7.37 |
| plain store | 6.74-6.78 |
| none (all work, no writes) | 6.37-6.42 |

- Every FP64 atomic goes past the L2: `TCC_EA0_ATOMIC` = `TCC_ATOMIC` = 1.35e8 per assembly.
  Workgroup scope compiles to the same instruction (`global_atomic_add_f64`, no cache bits) and
  changes neither the time nor the counters.
- Without its atomics the kernel runs within the H100 range. The atomic operation itself costs
  about 1.4 ms.
- On H100 the L2 serves 96.6% of the atomic sectors (ncu `lts`). There, the same kernel's
  non-atomic variants were slower than the atomic one, which is not explained, so the ablation does
  not split the H100 time.

## Tried, not adopted

- Block order: the blocks in 6 contiguous groups, one per XCD if workgroups are dispatched
  round-robin. L2 hits rise from 1.7e7 to 4.2e7 and memory reads fall from 5.9e7 to 3.4e7, but the
  assembly takes 8.24-8.32 ms against 8.15-8.18 ms. With 7 groups (control): 8.46-8.51 ms.
- 64 elements per block: 1.19e8 atomics instead of 1.35e8, but 9.39-9.46 ms, and 7.18-7.22 ms even
  without writes.

Open: an atomic-free scatter. Its saving is bounded by the gap between the atomic and the plain `+=`
update (about 0.85 ms), before its own cost.
