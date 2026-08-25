# Reproducing the matrix-free talk figures

One row per slide claim: which driver produces it, how to run it, what to read
off the output. Build targets are yours to run; only launch lines are listed.

## Building

The GPU drivers below are gated on three flags -- see
[examples/distributed/unstructured/CMakeLists.txt](../examples/distributed/unstructured/CMakeLists.txt):

```
-DMARS_ENABLE_CUDA=ON -DMARS_ENABLE_UNSTRUCTURED=ON -DMARS_ENABLE_FEM_EXAMPLES=ON
```

A full known-good configure (GH200; `CMAKE_CUDA_ARCHITECTURES` is arch-specific
-- 90 for GH200/H100, 80 for A100):

```bash
cmake .. \
  -DMARS_ENABLE_KOKKOS=OFF \
  -DMARS_ENABLE_CUDA=ON \
  -DMARS_ENABLE_TESTS=ON \
  -DMARS_ENABLE_UNSTRUCTURED=ON \
  -DMARS_ENABLE_FEM_EXAMPLES=ON \
  -DMARS_ENABLE_HYPRE=ON \
  -DCMAKE_CUDA_ARCHITECTURES=90
make -j
```

`MARS_ENABLE_HYPRE` is OPTIONAL for everything in this document -- none of the
matrix-free drivers link it. It is there for the AMG-preconditioned solvers and
the AMR/NS drivers; leave it on if you build those too, drop it otherwise.
Likewise VTK and ADIOS2 are not needed here.

Two of the five drivers are host-only and are not CMake targets at all (see
sections 1 and 3); they compile with a single `g++` line.

Run the `srun` lines from your build directory; binary paths are relative to
it. Figures are redrawn by [`make_matfree_figs.py`](figures/make_matfree_figs.py)
(measured numbers are hardcoded there -- update them when you re-measure).

---

## 1. "Cheap to store" -- ~340x leaner at p=7  (fig_memory, fig_storage_per_dof)

**Source:** [mars_ho_memory_sweep.cpp](../examples/distributed/unstructured/mars_ho_memory_sweep.cpp)

Assembled CSR grows ~(2p+1)^3 per DOF (324 -> 40,500 B/DOF, p=1 -> p=7);
matrix-free shrinks (421 -> 121 B/DOF). Crossover at p=1.

Host-only, no GPU, no MPI. Not a CMake target -- compile directly:

```bash
g++ -O2 -I. examples/distributed/unstructured/mars_ho_memory_sweep.cpp -o mars_ho_memory_sweep
./mars_ho_memory_sweep
```

Reads: bytes/DOF for assembled vs matrix-free, per order. nnz counted exactly
from the DOF map, not estimated.

---

## 2. "Free to run" -- throughput flat p=1..8  (fig_throughput)

**Source:** [mars_cvfem_ho_matfree_test.cu](../examples/distributed/unstructured/mars_cvfem_ho_matfree_test.cu) &middot; [mars_cvfem_ho_compare.cu](../examples/distributed/unstructured/mars_cvfem_ho_compare.cu) &middot; kernels: [mars_cvfem_ho_matfree.hpp](../backend/distributed/unstructured/fem/mars_cvfem_ho_matfree.hpp), [mars_cvfem_ho_matfree_shfl.hpp](../backend/distributed/unstructured/fem/mars_cvfem_ho_matfree_shfl.hpp)

PA (store d_G) 6-12 GDOF/s, peak at p=3; MF (recompute) ~3x slower.
Sum-factorization makes the apply O(p^4), not O(p^6).

IMPORTANT: the published numbers are the REGISTER + WARP-SHUFFLE kernel
([`mars_cvfem_ho_matfree_shfl.hpp`](../backend/distributed/unstructured/fem/mars_cvfem_ho_matfree_shfl.hpp)), not the base kernel. That kernel keeps the
face data in registers and exchanges it with `__shfl`, which drops the
shared-memory pipe from LSU ~99% to ~50% at p=4 and roughly doubles throughput.
The base [`mars_cvfem_ho_matfree.hpp`](../backend/distributed/unstructured/fem/mars_cvfem_ho_matfree.hpp) path is the readable reference; the shfl
path is what the slide measures.

```bash
srun -A csstaff -t 00:10:00 -N1 -n1 \
  ./examples/distributed/unstructured/mars_cvfem_ho_matfree_test
```

Reads: the per-order sweep table (GDOF/s, GFLOP/s, B/DOF). Also runs the four
correctness gates (metric vs host, single-element bit-exact, A*1 = 0,
A*linear = 0) plus a sheared-mesh gate.

Against the assembled operator (crossover at p=1, where assembled SpMV wins):

```bash
srun -A csstaff -t 00:10:00 -N1 -n1 \
  ./examples/distributed/unstructured/mars_cvfem_ho_compare --ncells=256 --iters=50
```

---

## 3. "~570x accuracy per unknown"  (fig_convergence)

**Source:** [mars_cvfem_ho_convergence.cpp](../examples/distributed/unstructured/mars_cvfem_ho_convergence.cpp) &middot; [mars_cvfem_ho_mass.hpp](../backend/distributed/unstructured/fem/mars_cvfem_ho_mass.hpp) &middot; [mars_cvfem_ho_basis.hpp](../backend/distributed/unstructured/fem/mars_cvfem_ho_basis.hpp)

Manufactured solution u = sin(pi x) sin(pi y) sin(pi z), homogeneous Dirichlet,
consistent load b = M f built with the CVFEM mass operator
([`mars_cvfem_ho_mass.hpp`](../backend/distributed/unstructured/fem/mars_cvfem_ho_mass.hpp), Knaus Alg 1 -- this is the ONLY thing that header is
used for). Dense LU solve, true global L2 error by Gauss quadrature.
Expect ~p+1 rate; at a fixed ~2,200 DOF budget p=4 is ~570x more accurate
than p=1.

Host-only, small dense solves. Not a CMake target:

```bash
g++ -O2 -I. examples/distributed/unstructured/mars_cvfem_ho_convergence.cpp -o mars_cvfem_ho_convergence
./mars_cvfem_ho_convergence
```

Reads: L2 error and observed rate per order, and the accuracy-per-DOF table.

---

## 4. "It distributes for free" -- weak scaling + comm  (fig_weakscale, fig_ho_scaling, fig_comm_pergpu)

**Source:** [mars_cvfem_ho_weakscale.cu](../examples/distributed/unstructured/mars_cvfem_ho_weakscale.cu)

~98% apply weak-scaling; comm fraction falls as V^(-1/3) with per-GPU size
(51.9% -> 19.2% measured as per-GPU DOF grows 2M -> 64M). `--overlap` hides the
forward halo behind the interior apply (bit-exact vs blocking).

```bash
MARS_NODEHALO_V2=1 srun -A csstaff -t 00:20:00 -N8 -n32 \
  ~/affinity/bind_numa.sh \
  ./examples/distributed/unstructured/mars_cvfem_ho_weakscale --p=1 --overlap
```

Hold per-rank size fixed and scale ranks (1/8/32 nodes) for the weak-scaling
curve. `--irregular` runs a deliberately jittered partition (~30-40% more
cross-rank sharing) to show ownership is not cube-specific.

Reads: per-rung efficiency, comm %, aggregate GDOF/s. Gates A*1 = 0 and
A*linear = 0 fire at every rung -- if they pass, the cross-rank assembly is
correct.

---

## 5. "Toward a trillion"  (fig_trillion)

**Source:** [mars_ho_dist_apply_test.cu](../examples/distributed/unstructured/mars_ho_dist_apply_test.cu)

The correctness invariant at any scale is A*1 = 0: u = 1 on owned DOF, 0 on
ghost -> forward -> apply owned elements -> reverseAdd -> max|y| over owned DOF
must be machine zero. One mis-paired send leaves an O(1) residual at the seam.
Measured 1e-18.

Small multi-rank gate first:

```bash
MARS_NODEHALO_V2=1 srun -A csstaff -t 00:10:00 -N1 -n4 \
  ~/affinity/bind_numa.sh \
  ./examples/distributed/unstructured/mars_ho_dist_apply_test --ncells=16 --p=2
```

Scale up by raising --ncells and the node count. GPU DOF numbering
(`buildDistributedGpu`) is the DEFAULT and is ~10x faster than host
(8 s vs 80 s/rank at 625M DOF/GPU); `--host-numbering` opts out, and
`--self-check` compares the two numbering paths by permutation-invariant
quantities on the same config.

Two trillion-DOF matvecs were measured on 2048 GH200 (1.005e12 DOF,
491M DOF/GPU, 512 nodes), both bit-exact:
  - PA / store-d_G  : 111.9 B/DOF -- fast, but d_G fills HBM -> capped near 1T
  - MF / recompute  :  26.8 B/DOF -- ~3x slower per matvec, scales to ~6T
At very large global counts set `MARS_GLOBAL_BUCKETSIZE` (or rely on the
auto-scaled default) -- cstone's replicated global-tree count Allreduce
overflows a 32-bit byte count above ~2 GiB otherwise.

---

## Kernel-variant switches (mars_ho_dist_apply_test)

Which kernel actually runs is chosen by environment variable, not a CLI flag.
This matters: the PUBLISHED throughput is the shuffle kernel, which is opt-in.
Precedence: RECOMPUTE > FP32_METRIC > SHFL; NEK is excluded by the first two.

| env | selects |
|-----|---------|
| (none)                | base PA kernel, stored d_G (fp64) |
| `MARS_HO_SHFL=1`      | register + warp-shuffle PA kernel -- the published throughput |
| `MARS_HO_RECOMPUTE=1` | MF path: no d_G, metric rebuilt inline (~3x slower, ~5x higher ceiling) |
| `MARS_HO_FP32_METRIC=1` | d_G stored as float (halves the dominant DRAM stream) |
| `MARS_HO_NEK=1`       | Nek-style variant (l-invariant normal column hoisted to registers) |
| `MARS_HO_EBLOCK=<n>`  | elements per block (occupancy tuning) |

Numbering: GPU is the DEFAULT. Opt out with `--host-numbering` or
`MARS_HO_HOST_NUMBERING=1`. `--gpu-numbering` is accepted but is a no-op, and
there is NO `MARS_HO_GPU_NUMBERING` variable -- it never existed in the code.

So the two trillion rows are the same binary, different env:
  PA / store-d_G  ->  `MARS_HO_SHFL=1`
  MF / recompute  ->  `MARS_HO_RECOMPUTE=1`

---

## Not part of this story

`mars_ho_dmma_bench` benchmarks a DIFFERENT operator -- the Galerkin/spectral
Laplacian ([`mars_ho_laplacian_dmma.hpp`](../backend/distributed/unstructured/fem/mars_ho_laplacian_dmma.hpp)), not CVFEM. It shares no code with
[`mars_cvfem_ho_matfree.hpp`](../backend/distributed/unstructured/fem/mars_cvfem_ho_matfree.hpp). It is the parked FP64 tensor-core track, kept as a
benchmark. Do not quote its numbers on a CVFEM slide.

---

## Reading the operator itself

| File | Role |
|------|------|
| [mars_cvfem_ho_apply.hpp](../backend/distributed/unstructured/fem/mars_cvfem_ho_apply.hpp) | host reference -- **start here**, 165 lines |
| [mars_cvfem_ho_basis.hpp](../backend/distributed/unstructured/fem/mars_cvfem_ho_basis.hpp) | builds Btil / Dtil / D / W |
| [mars_cvfem_ho_matfree.hpp](../backend/distributed/unstructured/fem/mars_cvfem_ho_matfree.hpp) | GPU PA + MF kernels |
| [mars_cvfem_ho_matfree_shfl.hpp](../backend/distributed/unstructured/fem/mars_cvfem_ho_matfree_shfl.hpp) | register + warp-shuffle kernel (published throughput) |

Docs: [TUTORIAL_matfree.md](TUTORIAL_matfree.md) &middot; [MATFREE_SCALING_ACHIEVEMENT.md](MATFREE_SCALING_ACHIEVEMENT.md)
