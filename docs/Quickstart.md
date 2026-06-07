# Quickstart — Build and Run Your First Assembly

This walks you from a fresh checkout to a running GPU assembly in a few minutes,
using a cube mesh you generate yourself (no external data needed). By the end you
will have built MARS, generated a mesh, run the CVFEM graph-assembly example, and
read its output. It is the "hello world" for the
[FEM Assembly](FEM-Assembly.md) pipeline.

---

## 1. Prerequisites

- A CUDA (NVIDIA) or HIP (AMD) GPU, with the matching toolkit. MARS is GPU-native;
  the assembly runs on the device.
- CMake ≥ 3.20 and a C++20 compiler.
- MPI (for multi-rank runs; single-rank works without `mpirun`).
- Python 3 with NumPy (only to generate the example mesh).

Dependencies are managed with Spack environments under `spack-envs/` if you want a
reproducible toolchain, but a system CUDA + CMake is enough to follow along.

---

## 2. Get the code

```bash
git clone https://github.com/dganellari/mars.git
cd mars
```

---

## 3. Build

A minimal GPU build with the unstructured backend and the FEM examples:

```bash
mkdir build && cd build
cmake .. \
  -DMARS_ENABLE_CUDA=ON \
  -DMARS_ENABLE_UNSTRUCTURED=ON \
  -DMARS_ENABLE_FEM_EXAMPLES=ON \
  -DCMAKE_CUDA_ARCHITECTURES=90        # 80 for A100, 90 for H100/GH200
make mars_cvfem_graph -j
cd ..
```

For an AMD (HIP) GPU, swap the CUDA flags:

```bash
cmake .. -DMARS_ENABLE_HIP=ON -DHIP_GPU_ARCHITECTURES=gfx942 \
         -DMARS_ENABLE_UNSTRUCTURED=ON -DMARS_ENABLE_FEM_EXAMPLES=ON
```

On a CPU-only laptop (for browsing/compiling, not GPU runs) use the Kokkos backend:
`cmake .. -DMARS_ENABLE_KOKKOS=ON`.

---

## 4. Generate a mesh

MARS reads a simple structure-of-arrays binary mesh (a directory of `x.float32`,
`y.float32`, `z.float32` coordinate arrays plus `i0.int64 …` connectivity columns).
The repo ships a generator that makes a structured cube of hexahedra:

```bash
python3 scripts/generate_hex_cube.py --nx 16 --ny 16 --nz 16 --output cube16
```

This writes a `cube16/` mesh directory (a 16×16×16 grid: 4096 hex elements, 4913
nodes). There is a tet version too — `scripts/generate_tet_cube.py` — if you want to
exercise the tet assembler instead. Bump `--nx/--ny/--nz` for a bigger mesh once the
small one runs.

---

## 5. Run

Single GPU:

```bash
./build/examples/mars_cvfem_graph --mesh=cube16
```

Four GPUs (MPI) — the same binary, the mesh is partitioned automatically by the
space-filling curve:

```bash
mpirun -np 4 ./build/examples/mars_cvfem_graph --mesh=cube16
```

Useful options (`--help` lists them all):

| Option | Meaning |
|--------|---------|
| `--mesh=DIR` | the mesh directory (required) |
| `--kernel=VARIANT` | assembly kernel: `tensor`, `shmem`, `optimized`, `team`, `original` (see [CVFEM Kernels](CVFEM-Kernels.md)) |
| `--iterations=N` | repeat the assembly N times (for timing) |
| `--block-size=N` | CUDA block size (default 256) |
| `--quiet` | suppress the per-phase prints |

---

## 6. Read the output

The example prints a line per pipeline phase, so you can watch the
[assembly stages](FEM-Assembly.md) happen in order:

```
PHASE: domain constructed
PHASE: domain synced
PHASE: nodeCount=4913 elementCount=4096
PHASE: DOF mapping done, numDofs=4913 numTotalDofs=4913
PHASE: starting sparsity build
PHASE: sparsity done, nnz=...
```

Reading it against the pipeline:

- **domain constructed / synced** — the mesh was read and partitioned across the
  GPUs (cornerstone octree). On multiple ranks each rank now owns a contiguous
  space-filling-curve chunk plus a halo.
- **nodeCount / elementCount** — the mesh size this rank sees (owned + halo).
- **DOF mapping done, numDofs** — the node→DOF numbering (`buildDofMappingGpu`); on a
  single rank `numDofs` equals the node count. On 4 ranks, the **sum of `numDofs`
  across ranks** equals the global unique node count (each boundary node owned by
  exactly one rank).
- **sparsity done, nnz** — the CSR structure is built; `nnz` is the number of
  nonzeros. For the 7-point graph stencil this is ≈ 7 × numDofs.

After this the example assembles the CSR values. At that point you hold an assembled
sparse system in GPU memory — the handoff point where you would attach a solver.

---

## 7. What just happened (the mental model)

You ran the whole [FEM Assembly](FEM-Assembly.md) pipeline:

```
mesh (cube16)  →  partition (SFC)  →  node→DOF map  →  CSR sparsity  →  assemble  →  [your solver]
```

Everything after "mesh" ran on the GPU. The same pipeline scales from this 4096-element
cube to billions of elements across a GPU cluster — only the mesh and the rank count
change.

---

## 8. Where to go next

- **Understand the pipeline you just ran:** [FEM Assembly](FEM-Assembly.md).
- **Tune the assembly kernel:** [CVFEM Kernels (GPU)](CVFEM-Kernels.md).
- **A full CFD application built on this:** the
  [Pump Flow tutorial](pump_cfd_tutorial.md) (incompressible Navier–Stokes).
- **How the mesh is read and split:**
  [Mesh Reading and Partitioning](Mesh-Reading-and-Partitioning.md),
  [SFC Mapping](SFC-Mapping.md), [Multi-Rank Support](Multi-Rank-Support.md).
