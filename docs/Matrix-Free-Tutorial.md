# High-Order Matrix-Free Operator

> **Experimental in v0.1.** The drivers below check the operator apply on one GPU and on several
> ranks, but MARS v0.1 has no linear solver that uses it. See
> [Known limitations](https://github.com/dganellari/mars/blob/master/KNOWN_LIMITATIONS.md).

This page explains what the high-order (HO) matrix-free operator computes, how to build and
run the two drivers that test it, and how to read their output.

## What the operator is

The operator is the control-volume finite element (CVFEM) form of the scalar diffusion
operator `-div(grad u)` on hexahedral meshes, for polynomial orders `p = 1` to `8`. It follows
the high-order CVFEM method of Knaus (2022), Algorithm 2.

*Matrix-free* means that MARS computes `y = A x` without storing the matrix `A`. One call is one
operator application, a matrix-vector product (matvec). It is **not** a linear solve: v0.1 has
no Krylov method or preconditioner for this operator. The operator is not symmetric for
`p >= 2`, so plain CG does not apply to it; the host accuracy study below uses a dense LU solve
for this reason.

The operator is diffusion only. It has no advection term, and it is not the Navier–Stokes solver
of the [Poiseuille](poiseuille_tutorial.md) and [Taylor–Green](periodic_tgv_tutorial.md)
tutorials, which uses linear (`p = 1`) elements and assembled matrices.

## How one apply works

An element of order `p` has `(p+1)^3` nodes on Gauss–Lobatto–Legendre points. Each node is a
degree of freedom (DOF); nodes on element faces, edges and corners are shared with neighbours.

- **Reference operators.** Four small 1D matrices (`Btil`, `Dtil`, `D`, `W`) depend only on `p`.
  `buildHoCvfemOperators` (`fem/mars_cvfem_ho_basis.hpp`) builds them once; the GPU keeps them
  in constant memory.
- **Geometry.** Each element enters through a metric `detJ J^-1 J^-T` at the points of its
  sub-control surfaces: `9 p (p+1)^2` doubles per element (`d_G`). It is computed once from the
  8 corner coordinates and reused by every apply. A straight-sided hex of any shape works; the
  cross terms of the metric carry the shear.
- **Sum factorization.** The element apply is a sequence of 1D contractions along x, y and z, so
  its cost per element grows like `p^4`, not like `p^6` as for a dense
  `(p+1)^3 x (p+1)^3` element matrix.
- **Gather, apply, scatter.** For each element the kernel reads its values of `x`, applies the
  operator and adds the result into `y` with `atomicAdd`. `y` must be zero before the call.

| File (under `backend/distributed/unstructured/`) | Content |
|---|---|
| `fem/mars_cvfem_ho_basis.hpp` | the 1D reference operators |
| `fem/mars_cvfem_ho_apply.hpp` | host reference: `computeElementMetric`, `applyHoCvfemElement` |
| `fem/mars_cvfem_ho_matfree.hpp` | GPU metric and apply: `ho_cvfem_metric_perpoint_launch`, `ho_cvfem_apply_launch` |
| `fem/mars_cvfem_ho_matfree_shfl.hpp` | register and warp-shuffle variant of the apply |
| `fem/mars_ho_dof_handler.hpp`, `fem/mars_ho_dof_handler_gpu.hpp` | DOF numbering (host reference, device) |
| `fem/mars_ho_halo.hpp` | ownership of shared DOFs and the DOF halo `HoHalo` |

On several ranks three more pieces are needed:

1. **Numbering.** `buildDistributedGpu` numbers the DOFs on the device. A DOF on an edge or face
   is named by the sorted global ids of the entity's corners plus its position inside the entity
   (`DofKey`). Every rank that holds the DOF computes the same name, without communication.
2. **Ownership.** A corner DOF belongs to the rank that owns the mesh node. An edge or face DOF
   belongs to the lowest rank among the ranks that hold it. An element-interior DOF belongs to
   the rank of its element. The numbering registers corners only through the elements a rank
   owns, so it needs every node owner to own an element at its node. The older vote node
   ownership guarantees that; the default SFC ownership does not. The multi-rank drivers select
   vote ownership themselves (`MARS_OWNERSHIP=vote`, set before the domain is built); code of
   your own needs the same setting, and the halo stops with a message if it is missing.
3. **Halo.** One distributed matvec is `forward` (owners copy their values to the ghost copies),
   the apply over the elements this rank owns, then `reverseAdd` (ghost contributions are added
   to the owners).

## Build

The drivers are part of the FEM examples and need CUDA:

```bash
cmake -B build -DMARS_ENABLE_CUDA=ON -DMARS_ENABLE_FEM_EXAMPLES=ON -DCMAKE_CUDA_ARCHITECTURES=90
cmake --build build --target mars_cvfem_ho_matfree_test mars_ho_dist_apply_test -j
```

`CMAKE_CUDA_ARCHITECTURES` is 90 for H100/GH200 and 80 for A100. The commands below run from the
repository root.

## One GPU: `mars_cvfem_ho_matfree_test`

```bash
./build/examples/distributed/unstructured/mars_cvfem_ho_matfree_test
```

It takes no mesh. For each order `p = 1 ... 7` it builds a structured unit cube (from `48^3`
elements at `p = 1` down to `16^3` at `p = 7`), numbers the DOFs on the device and runs these
gates:

| Gate | Check | Pass if |
|---|---|---|
| A | GPU metric of one element against `computeElementMetric`; GPU apply of one element against `applyHoCvfemElement` | max error < 1e-12 (relative for the apply) |
| B | `A·1`: a constant field has no flux | max \|y\| < 1e-9 at every DOF |
| C | `A·x` for the linear field `x` | max \|y\| < 1e-9 at interior DOFs |
| shear | `p = 1, 2, 4` on a cube where every element has its own shear: metric and full apply against the host | metric error < 1e-12, apply < 1e-11 relative |

It then times 100 applies per order. One line per order, for example:

```
p=1 E=48 nDof=117649 nEl=110592 | metricErr=... elemRel=... | A*1=... A*lin(int)=... | PASS
    perf: ... ms/apply | ... MDOF/s | ... GB/s (useful traffic est)
```

After the shear lines it prints an order-sweep table and the verdict:

```
HO-CVFEM matrix-free GPU gate: PASS
```

The exit code is 0 if every gate passes and 1 otherwise. In the sweep table, the column
`assembled B/DOF` is an estimate with `(2p+1)^3` nonzeros on every row. Only rows of mesh
vertices have that many; see [Memory](#memory) for exact counts.

`--dof-self-check` (or `MARS_HO_DOF_SELF_CHECK=1`) first compares the device numbering with the
host `HODofHandler::build()` for `p = 1 ... 7`. The two numberings differ by a permutation, so the
check compares the DOF counts, the set of `DofKey`s and which element slots share a DOF.

## Several GPUs: `mars_ho_dist_apply_test`

```bash
mpirun -np 4 ./build/examples/distributed/unstructured/mars_ho_dist_apply_test --ncells=16 --p=3
```

Each rank generates its part of a unit cube with `ncells^3` hexes, so no mesh file is needed.
Each rank uses GPU `rank % deviceCount`. The driver builds the distributed numbering, the
ownership and the halo, and then applies the operator once to `x = 1`. The result must vanish on
every owned DOF. A wrong pairing across ranks leaves an error of order 1 at the rank
boundaries. The gate runs twice, with the host halo and with the device halo, and then the
driver times 50 matvecs:

```
p=3  owned DOF=...  max|A.1| over owned (host halo, gpu-numbering) = ...   [PASS]
p=3  device A.1 = ... [PASS]   full matvec ... MDOF/s | apply-only ... MDOF/s | comm ...%  [path=STORE-d_G]
       THROUGHPUT(p=3): apply-only ... GDOF/s | ... FLOP/DOF | ... GFLOP/s (FP64) | ... DOF/GPU  [path=STORE-d_G]
       measured: ... ms/matvec (slowest of 4 ranks, 50 iters) | sustained ... TDOF/s | ... global DOF  [path=STORE-d_G]
       STORAGE (p=3):  store-d_G = ... B/DOF  |  recompute(corners) = ... B/DOF  |  ...
```

- A gate passes when the maximum is below 1e-8.
- `full matvec` includes the halo exchange and `apply-only` does not;
  `comm = (full - apply) / full`. All times are those of the slowest rank.
- The driver exits with 0 even when a gate fails, so read the `[PASS]`/`[FAIL]` tags. If a
  rank runs out of GPU memory, the driver prints `DEVICE OOM` and stops all ranks.

| Option | Meaning |
|---|---|
| `--ncells=N` | cube with `N^3` elements (default 16) |
| `--p=P` | order, 1 to 8 (default 2) |
| `--iters=N` | number of timed matvecs (default 50) |
| `--irregular` | move interior nodes by up to a quarter of the cell size: distorted hexes and an irregular partition |
| `--sweep` | run `p = 1 ... 8` in one job, with about the same DOF count per GPU (`ncells / p` per order), and print one table |
| `--MF` or `--recompute` | do not store the metric; rebuild it from the corners in every apply |
| `--host-numbering` | number the DOFs on the host instead of the device |
| `--self-check` | number on the host and on the device, compare them, then run the gate with the device numbering |

| Environment variable | Effect |
|---|---|
| `MARS_HO_SHFL=1` | use the register and warp-shuffle kernel instead of the default kernel |
| `MARS_HO_FP32_METRIC=1` | store the metric in single precision; the arithmetic stays double |
| `MARS_HO_SKIP_HOST_GATE=1` | skip the host-halo gate; it is skipped anyway above 50M DOF per rank |
| `MARS_GLOBAL_BUCKETSIZE=N` | bucket size of the cornerstone global tree (see below) |

`--MF` turns off `MARS_HO_FP32_METRIC`, and `MARS_HO_FP32_METRIC` turns off `MARS_HO_SHFL`.
With Slurm, make sure the variables reach the ranks (`srun --export=ALL`).

The global octree of cornerstone is replicated on every rank. For large cubes the driver raises
its bucket size so that the reduction of the tree counts stays below 2 GiB, and prints
`[build] global bucketSize=...`; `MARS_GLOBAL_BUCKETSIZE` overrides the choice.

## Results you can reproduce

Both programs below run on the host (any CPU, no GPU, no MPI) in a few seconds.

### Accuracy

```bash
g++ -std=c++20 -O2 -I. examples/distributed/unstructured/mars_cvfem_ho_convergence.cpp -o mars_cvfem_ho_convergence
./mars_cvfem_ho_convergence
```

It solves `-div(grad u) = 3 pi^2 u_exact` on the unit cube with
`u_exact = sin(pi x) sin(pi y) sin(pi z)` and `u = 0` on the boundary, with the HO CVFEM operator
and a dense LU solve, and prints the relative L2 error against `u_exact`. Output (Apple clang on
macOS; another compiler can change the last digits):

```
 p=1 :  #DOFs   |   L2 error   | rate
            125 | 1.6808e-01 |   -
            343 | 7.8594e-02 | 2.26
            729 | 4.5010e-02 | 2.22
           2197 | 2.0264e-02 | 2.17

 p=2 :  #DOFs   |   L2 error   | rate
            125 | 4.3029e-02 |   -
            343 | 1.2109e-02 | 3.77
            729 | 4.9921e-03 | 3.53
           2197 | 1.4524e-03 | 3.36

 p=3 :  #DOFs   |   L2 error   | rate
            343 | 3.5038e-03 |   -
           1000 | 6.9237e-04 | 4.55
           2197 | 2.1899e-04 | 4.39

 p=4 :  #DOFs   |   L2 error   | rate
            125 | 2.0099e-03 |   -
            729 | 2.6177e-04 | 3.47
           2197 | 3.4703e-05 | 5.49
```

On the finest meshes the rate is close to `p + 1`. With the same 2197 unknowns, the `p = 4`
error is about 580 times smaller than the `p = 1` error. This holds for a smooth solution; with
corners, boundary layers or curved walls the gain is smaller.

### Memory

```bash
g++ -std=c++20 -O2 -I. examples/distributed/unstructured/mars_ho_memory_sweep.cpp -o mars_ho_memory_sweep
./mars_ho_memory_sweep
```

It counts the nonzeros of the assembled matrix exactly from the DOF map (two DOFs couple when
they share an element), at 12 bytes per nonzero (8-byte value, 4-byte column index). Its own
matrix-free column assumes one constant metric per element; the operator on this page stores the
metric at every sub-control-surface point instead. So the last column below is computed for the
operator on this page, on the same cubes: the metric (`9 p (p+1)^2` doubles per element) plus the
element-to-DOF map (`(p+1)^3` 4-byte ids per element), divided by the DOF count.

| p | cube | DOFs | nonzeros per row | assembled matrix, B/DOF | matrix-free, B/DOF |
|---|---|---|---|---|---|
| 1 | 40^3 | 68,921 | 26 | 308.5 | 297 |
| 2 | 24^3 | 117,649 | 61 | 733.3 | 165 |
| 4 | 12^3 | 117,649 | 205 | 2,462.0 | 113 |
| 7 | 6^3 | 79,507 | 685 | 8,216.6 | 93 |

Neither column counts the vectors `x` and `y`, which both need. At `p = 1` the two are about the
same size, and the 7-nonzero graph matrix of [FEM Assembly](FEM-Assembly.md) (a lumped operator,
at most 84 B/DOF for values and column indices) is smaller still: at `p = 1` use the assembled
kernels. From `p = 2` on, the assembled matrix grows quickly with `p` while the matrix-free
storage per DOF falls. For a run of `mars_ho_dist_apply_test`, the `STORAGE` line prints the
metric bytes per DOF of that run.

### GPU throughput and scaling

This page quotes no GPU throughput or scaling numbers. The two drivers print them for every
run (MDOF/s, GDOF/s, ms per matvec, the communication share), and they depend on the GPU, the
kernel variant, the order and the number of DOFs per GPU. Measure them on your system.

## Other drivers

| Target | What it does |
|---|---|
| `mars_cvfem_ho_compare` | one GPU, `p = 1`: compares the matrix-free apply with the assembled 27-nonzero matrix (parity) and with the 7-nonzero graph matrix (throughput and memory). `--ncells=N` or `--mesh=DIR`, `--iters=N` |
| `mars_cvfem_ho_weakscale` | `p = 1` distributed apply for weak-scaling runs. `--ncells=N` (global cube) or `--mesh=DIR`, `--overlap`, `--irregular` |
| `mars_ho_dist_dof_test` | distributed HO numbering and halo: global DOF count and a constant forwarded to every ghost |
| `mars_cvfem_ho_tet_test`, `mars_ho_dist_apply_tet_test` | tetrahedral HO operators (collapsed sum factorization), also experimental |

## Release tests

`ctest -L release` (see the [README](https://github.com/dganellari/mars#checking-your-build))
includes two tests of this operator:

- `marsReleaseHoMatfree` runs `mars_cvfem_ho_matfree_test` on one GPU.
- `marsReleaseHoDistApply_npN` runs `mars_ho_dist_apply_test --ncells=16 --p=3` on N ranks
  (`MARS_RELEASE_TEST_RANKS`, default 4) and requires both gates to pass.

## Reference

R. Knaus, "A fast matrix-free approach to the high-order control volume finite element method
with application to low-Mach flow", *Computers & Fluids* 239 (2022) 105408.
