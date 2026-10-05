# Quickstart: Build and Run Your First Assembly

This page takes you from a fresh clone to a GPU assembly on a cube mesh that you generate
yourself (no external data). You build MARS, generate a mesh, run the CVFEM graph-assembly
example and read its output. It is the "hello world" for the [FEM Assembly](FEM-Assembly.md)
pipeline.

---

## 1. Prerequisites

- An NVIDIA GPU with the CUDA toolkit. The example drivers need CUDA. (A HIP build for AMD GPUs
  is possible, `-DMARS_ENABLE_HIP=ON -DHIP_GPU_ARCHITECTURES=gfx942`, but it was not re-verified
  for v0.1 and it does not build the examples.)
- CMake 3.22 or newer and a C++20 compiler.
- MPI. The build needs it by default; a single-rank run can start without `mpirun`.
- Python 3 with NumPy, to generate the mesh.
- A network connection at configure time: CMake fetches cornerstone-octree, and googletest when
  tests are on and no system googletest is found.

Spack environments under `spack-envs/` give a reproducible toolchain, but a system CUDA, MPI and
CMake are enough.

---

## 2. Get the code

```bash
git clone https://github.com/dganellari/mars.git
cd mars
```

---

## 3. Build

A GPU build with the unstructured backend and the FEM examples:

```bash
cmake -B build -DMARS_ENABLE_CUDA=ON -DMARS_ENABLE_FEM_EXAMPLES=ON -DCMAKE_CUDA_ARCHITECTURES=90
cmake --build build --target mars_cvfem_graph -j
```

`CMAKE_CUDA_ARCHITECTURES` is 90 for H100/GH200 and 80 for A100. The unstructured backend is on
by default in a CUDA build. Without `--target`, `cmake --build build -j` builds the library and
all examples.

On a machine without a GPU, a plain `cmake -B build` builds only the core library; the
unstructured backend needs CUDA or HIP.

---

## 4. Generate a mesh

MARS reads a simple structure-of-arrays binary mesh: a directory with the coordinate arrays
`x.float32`, `y.float32`, `z.float32` and one connectivity column per element corner
(`i0.int64` ... `i7.int64` for hexahedra). The repository ships a generator for a structured
cube of hexahedra:

```bash
python3 scripts/generate_hex_cube.py --nx 16 --ny 16 --nz 16 --output cube16
```

This writes the directory `cube16/`: a 16 × 16 × 16 grid of the unit cube, with 4096 hexahedra
and 4913 nodes. `scripts/generate_tet_cube.py` makes a tetrahedral cube for the tet assembler
(`mars_cvfem_graph_tet`). Increase `--nx/--ny/--nz` for a bigger mesh once the small one runs.

---

## 5. Run

One GPU:

```bash
./build/examples/distributed/unstructured/mars_cvfem_graph --mesh=cube16
```

Four GPUs (MPI). It is the same binary; the space-filling curve splits the mesh across the ranks,
and each rank uses GPU `rank % deviceCount`:

```bash
mpirun -np 4 ./build/examples/distributed/unstructured/mars_cvfem_graph --mesh=cube16
```

Options (run it without `--mesh` to print them):

| Option | Meaning |
|--------|---------|
| `--mesh=DIR` | the mesh directory (required) |
| `--kernel=VARIANT` | assembly kernel, see the table below (default `original`) |
| `--iterations=N` | repeat the assembly N times, for timing (default 10) |
| `--block-size=N` | CUDA block size (default 256) |
| `--bucket-size=N` | cornerstone octree bucket size (default 64; try 32 or 16 for 16+ ranks) |
| `--bucket-size-focus=N` | focus-tree bucket size (default 8; lower = finer halo, more memory) |
| `--quiet` | do not print the timing breakdown |

`--kernel` variants (all compute the same hex assembly; see [CVFEM Kernels](CVFEM-Kernels.md)):

| Variant | Notes |
|---------|-------|
| `original` | default; reference per-element kernel |
| `optimized` | optimized element kernel |
| `shmem` | low-register, direct-scatter kernel |
| `team` | warp per element |
| `tensor` | full local matrix per thread; the recommended default on current GPUs |
| `tensor_colored` | `tensor` with graph coloring, no atomics |
| `tensor_aos` | `tensor` with packed node data |
| `tensor_perip` | `tensor` with pre-resolved CSR positions |
| `tensor_perip_lb2` | `tensor_perip` with two blocks per SM |
| `smem_cache` | shared-memory node cache |
| `wmma_tensor` | FP64 tensor cores, WMMA (experimental) |
| `wgmma_tensor` | FP64 tensor cores, WGMMA, Hopper only (experimental) |

---

## 6. Read the output

The example prints one line per pipeline phase, so you can follow the
[assembly stages](FEM-Assembly.md) in order. On one rank:

```
PHASE: domain constructed
PHASE: domain synced
PHASE: nodeCount=4913 elementCount=4096
PHASE: DOF mapping done, numDofs=4913 numTotalDofs=4913
PHASE: starting sparsity build
PHASE: sparsity done, nnz=...
```

- **domain constructed / synced**: the mesh was read and split across the GPUs (cornerstone
  octree). On several ranks each rank now holds a contiguous space-filling-curve range of
  elements plus a halo.
- **nodeCount / elementCount**: the mesh size that rank 0 sees, owned plus halo.
- **numDofs**: the equation rows rank 0 owns (`buildDofMappingGpu`). On one rank it equals the
  node count. On several ranks every node is owned by exactly one rank, so the owned counts of
  all ranks add up to 4913.
- **nnz**: the number of nonzeros of the CSR structure. The graph stencil has at most 7 nonzeros
  per row: the node itself and its 6 edge neighbours. Boundary rows have fewer.

After `--iterations` assemblies the example prints one summary line:

```
Assembler: MARS CVFEM Hex (Graph+Lump original)....: ... milliseconds @ ... GB/s (average of 10 samples) [matrix norm: ..., rhs norm: ...]
```

The time is the average assembly time; the matrix and right-hand-side norms are reduced over all
ranks. The release test `marsReleaseHexAssembly` checks that these norms, and the assembled
rows themselves, are the same on 1 and on N ranks.

At this point you hold an assembled sparse system in GPU memory: the place where a solver takes
over.

---

## 7. What just happened

You ran the whole [FEM Assembly](FEM-Assembly.md) pipeline:

```
mesh (cube16)  ->  partition (SFC)  ->  node-to-DOF map  ->  CSR sparsity  ->  assemble  ->  [your solver]
```

Everything after reading the mesh ran on the GPU. Only the mesh and the rank count change for a
bigger problem.

---

## 8. Where to go next

- **Check your build:** `ctest -L release` (see the
  [README](https://github.com/dganellari/mars#checking-your-build)).
- **The pipeline you just ran:** [FEM Assembly](FEM-Assembly.md).
- **The assembly kernels:** [CVFEM Kernels](CVFEM-Kernels.md).
- **A full CFD application:** the [Poiseuille channel-flow tutorial](poiseuille_tutorial.md)
  (incompressible Navier–Stokes), then the [Taylor–Green vortex](periodic_tgv_tutorial.md) for
  periodic boundaries.
- **High-order operators:** the [matrix-free tutorial](Matrix-Free-Tutorial.md).
