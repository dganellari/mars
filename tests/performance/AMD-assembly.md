# CVFEM assembly on MI300A and H100 (`tensor`)

Status: shipped in v0.1.1. `tensor` takes a row map, built once on the GPU
(`fem/mars_cvfem_row_map.hpp`), for the graph and the full assembly, and runs `tensor_rows` on AMD and
`tensor_map` on NVIDIA. Without a map (`--no-row-map`, or a caller that leaves `Config::rowMap` unset)
the graph assembly runs the aggregated `tensor_faces` kernel on AMD (since 053620d9) and the one-thread
`tensor` kernel on NVIDIA, and the full assembly runs `tensor_faces` on both. On both GPUs the release
tests passed, and so did an entry-by-entry check of every kernel against the original one on 1 and 4
ranks, the `tensor` kernels with and without the map.

## Row map

For each block of 32 elements: the CSR rows it touches, sorted, and the place of each element node
among them (a byte). For each element row: the CSR offset that each of the element's 8 columns adds
to. A column outside the pattern goes to the diagonal in the graph (lumped) assembly and is dropped in
the full one.

- Columns: 4 bits each when no row has more than 15 entries and every block fits in shared memory
  (a compact map: every graph map here), else a byte each. 40 or 72 bytes per element, and 6 per block
  and per block row: 394 MiB for the graph map of hex192 (7.1e6 elements, 2.1e7 block rows), 610 MiB
  for its full map.
- `tensor_rows`: 16 threads per element compute the element matrix as `tensor_faces` does. The block
  adds up its contributions in 32 KB of shared memory, each row at its own length (the right-hand
  side, then the CSR row), then adds each row to the matrix with atomics, consecutive threads on
  consecutive entries. MI300A runs two of these blocks per CU (4 waves per SIMD, set by the registers),
  so the 32 KB keep that occupancy (compiler remarks).
- `tensor_map`: the one-thread `tensor` kernel with the map's offsets instead of a search in each CSR
  row, and without the ownership and diagonal loads.

Every mesh takes a map. A block that touches more than 254 rows, or whose rows need more than the
shared memory (3656 doubles), adds them straight to the matrix; a row longer than 254 entries is
found by a search in its CSR row. With the 28.9 KB the kernel had before (the element storage alone),
35719 of the 221184 blocks of the full map of the unstructured mesh below went straight to the matrix,
and the full assembly took 7.27-7.37 ms on MI300A; with 32 KB, 1904 blocks and 6.53-6.57 ms.

## Measured

One rank, `--iterations=20`, the `Assembler` line, run time, alternating rounds (2026-11-01,
a93365de). MI300A on Beverin (ROCm 7.2, SPX mode), H100 on Daint (CUDA 12.9). hex192 is the public
structured mesh; tri192 is unstructured, from `python3 scripts/generate_hex_trisplit.py -n 192 --nz 32`:
7.1e6 hexes, with graph rows of up to 9 entries and full rows of up to 39.

| Assembly | Mesh | MI300A ms | H100 ms |
| --- | --- | ---: | ---: |
| graph, `tensor` | hex192 | 4.81-4.83 | 4.86-4.90 |
| graph, `tensor` | tri192 | 4.93-4.96 | 4.74-4.78 |
| graph, `tensor --no-row-map` | hex192 | 8.21-8.24 | 6.60-7.38 |
| graph, `tensor --no-row-map` | tri192 | 11.63-11.72 | 7.22-7.78 |
| full, `tensor` | hex192 | 6.31-6.34 | 6.18-6.21 |
| full, `tensor` | tri192 | 6.53-6.57 | 6.40-6.41 |
| full, `tensor --no-row-map` (`tensor_faces`) | hex192 | 19.34-19.42 | 9.60-9.65 |
| full, `tensor --no-row-map` (`tensor_faces`) | tri192 | 21.18-21.31 | 10.05-10.10 |

The graph `--no-row-map` rows are from an earlier job the same day, before the 4-bit columns and the
32 KB, which that path does not use. Building the map takes
35-36 ms (graph) and 63-64 ms (full) on MI300A, and 8-9 and 17-19 ms on H100, once per mesh and
pattern. The full assembly without a map took 25.8 ms on MI300A (`tensor_early`) and 12.1 ms on H100
(the one-thread kernel) on hex192 before it ran `tensor_faces`.

The kernel of 0438d67d, which took rows of up to 8 entries only and stopped on tri192, ran the graph
assembly of hex192 in 4.55-4.58 ms on MI300A and 4.86-4.89 ms on H100 in the same rounds. On H100
`nsys` gives the same kernel time for it and for 628ef311, within 0.1% (3 runs each). The 5.6% on
MI300A is not explained: the occupancy is the same (4 waves per SIMD).

Measured at 0438d67d (2026-10-31), `mars_cvfem_scale --ncells=192 --iters=20`, ms per assembly on the
slowest rank, run time, 3 alternating rounds:

| Ranks | MI300A | MI300A, `--no-row-map` | H100 | H100, `--no-row-map` |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 4.61-4.63 | 8.10-8.12 | 4.98-5.00 | 8.01-8.16 |
| 4 | 1.26-1.28 | 2.15-2.16 | 1.35-1.35 | 2.11-2.17 |

The counters below are also from 0438d67d, hex192, graph assembly. There, `tensor_rows` took
10.82-10.88 ms on H100 and `tensor_map` 9.38-9.46 ms on MI300A: each GPU runs the other kernel.

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
- `tensor_rows` with each row padded to an odd number of doubles, as the fixed 9-double rows of
  0438d67d were to spread them over the shared-memory banks: graph hex192 6.3% behind 0438d67d instead
  of 5.9%, and more blocks too big for shared memory, so the full assembly took 6.40-6.43 ms (hex192)
  and 7.68-7.77 ms (tri192) against 6.29-6.32 and 7.27-7.37 ms, with 28.9 KB of shared memory.
- Byte columns for every map: on H100 the graph assembly of hex192 took 4.90-5.91 ms over 5 rounds,
  against 4.89-4.95 ms with 4-bit columns. On MI300A it was 3.4% behind 0438d67d in its run, and 5.9%
  with 4-bit columns in another; the two were not compared in the same run.

## Local check

```bash
python3 scripts/test_cvfem_tensor_rows.py
```

It runs the map build and the kernels through a host emulator (one thread per GPU thread, a barrier
per block) under address and undefined-behaviour sanitizers. On hex meshes with scrambled element
order, DOF numbers and ownership, ghost rows, a node without a DOF and two corners sharing a DOF, and
on fans of quads swept into hexes (rows of up to 1014 entries), graph and full,
`tensor_rows`, `tensor_map` and `tensor` must match `tensor_faces` entry by entry. The checks must take
every path of the map: blocks too big for shared memory, blocks of too many rows, searched rows, rows
longer than 8 doubles, 4-bit and byte columns. It does not replace a GPU run.
