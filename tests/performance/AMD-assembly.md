# Graph assembly on MI300A and H100 (`tensor`)

Status: shipped in v0.1.1. `tensor` takes a row map, built once on the GPU
(`fem/mars_cvfem_row_map.hpp`), and runs `tensor_rows` on AMD and `tensor_map` on NVIDIA. Without a
map (`--no-row-map`, or a caller that leaves `Config::rowMap` unset), AMD runs the aggregated
`tensor_faces` kernel (since 053620d9) and NVIDIA the one-thread `tensor` kernel. On both GPUs the
release tests passed, and so did an entry-by-entry check of every kernel against the original one on
1 and 4 ranks, the `tensor` kernels with and without the map.

## Row map

For each block of 32 elements: the CSR rows it touches, sorted, and the place of each element node
among them (a byte, so at most 255 rows per block). For each element row: the CSR offset that each of
the element's 8 columns adds to, 4 bits each. A column outside the pattern goes to the diagonal, as
in the lumped assembly. 40 bytes per element and 4 per block row: 352.7 MiB on hex192 (7.1e6 elements,
2.1e7 block rows).

- `tensor_rows`: 16 threads per element compute the element matrix as `tensor_faces` does. The block
  adds up its contributions in shared memory, 9 doubles per row (up to 8 entries and the right-hand
  side), then adds each row to the matrix with atomics, consecutive threads on consecutive entries.
- `tensor_map`: the one-thread `tensor` kernel with the map's offsets instead of a search in each CSR
  row, and without the ownership and diagonal loads.

Building the map stops with a message when a CSR row has more than 8 entries (a node and its edge
neighbours; 7 on a structured hex mesh) or a block touches more than 255 rows. `--no-row-map` runs the
kernels without it.

## Measured

Public hex192, one rank, `mars_cvfem_graph --iterations=20`, the `Assembler` line, run time,
5 alternating rounds (2026-10-31). MI300A on Beverin (ROCm 7.2, SPX mode), H100 on Daint (CUDA 12.9).

| Kernel | MI300A ms | H100 ms |
| --- | ---: | ---: |
| `tensor` | 4.56-4.59 | 4.87-4.92 |
| `tensor --no-row-map` | 8.25-8.29 | 6.97-7.61 |
| `tensor_rows` | 4.56-4.59 | 10.82-10.88 |
| `tensor_map` | 9.38-9.46 | 4.86-4.91 |

Building the map takes 31.1 ms on MI300A and 6.9 ms on H100 (run time), once per mesh and pattern.

`mars_cvfem_scale --ncells=192 --iters=20`, ms per assembly on the slowest rank, run time, 3 alternating
rounds (2026-10-31):

| Ranks | MI300A | MI300A, `--no-row-map` | H100 | H100, `--no-row-map` |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 4.61-4.63 | 8.10-8.12 | 4.98-5.00 | 8.01-8.16 |
| 4 | 1.26-1.28 | 2.15-2.16 | 1.35-1.35 | 2.11-2.17 |

- On MI300A `tensor_rows` makes 3.3e7 global atomic requests per assembly, against 1.35e8 for the
  aggregated `tensor_faces` kernel and 2.26e8 for `tensor_map` (`TCC_ATOMIC`). All of them go past
  the L2 (`TCC_EA0_ATOMIC` is the same).
- On H100 `tensor_map` executes 1.30e9 instructions and 3.6e7 global load requests, against 1.57e9 and
  8.0e7 for the kernel that searches the CSR rows (ncu). `tensor_rows` executes 3.4e9. It has 4.6 times
  fewer L2 reduction sectors than `tensor_map` (5.7e7 against 2.6e8), but the L2 serves 98% of those
  of `tensor_map`.
- On AMD `tensor_map` declares its 256 threads per block (`amdgpu_flat_work_group_size`). Without that
  the compiler assumes 1024, keeps the kernel to 128 VGPRs and spills 782: 23.8-24.1 ms.

## Why the aggregated kernel is slow on MI300A

The aggregated `tensor_faces` kernel (2026-10-30), before the row map. A timing-only ablation (not kept
in the tree) replaced its global updates. Only the atomic update assembles the right matrix.

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

- `tensor_rows` storing the rows that only one block touches, instead of adding them atomically:
  4.80-4.82 ms against 4.48-4.50 ms. The atomic requests fall only 1% (3.24e7).
- `tensor_rows` without atomics: each block writes its partial rows, and a second kernel adds them up.
  5.59-5.63 ms for both kernels, and a 1899 MiB map built in 146 ms.
- Block order, aggregated kernel: the blocks in 6 contiguous groups, one per XCD if workgroups are
  dispatched round-robin. L2 hits rise from 1.7e7 to 4.2e7 and memory reads fall from 5.9e7 to
  3.4e7, but the assembly takes 8.24-8.32 ms against 8.15-8.18 ms. With 7 groups (control):
  8.46-8.51 ms.
- 64 elements per block, aggregated kernel: 1.19e8 atomics instead of 1.35e8, but 9.39-9.46 ms, and
  7.18-7.22 ms even without writes.

Not yet: only `mars_cvfem_graph` and `mars_cvfem_scale` build the map. The other callers of `Tensor`
(`assembleCvfemDiffusion` and the drivers) keep the kernels without it.

## Local check

```bash
python3 scripts/test_cvfem_tensor_rows.py
```

It runs the map build and the kernels through a host emulator (one thread per GPU thread, a barrier
per block) under address and undefined-behaviour sanitizers. On hex meshes with scrambled element
order, DOF numbers and ownership, ghost rows, a node without a DOF and two corners sharing a DOF,
`tensor_rows`, `tensor_map` and `tensor` must match `tensor_faces` entry by entry, and both limits of
the map must stop the build. It does not replace a GPU run.
