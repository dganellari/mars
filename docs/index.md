# MARS Documentation

MARS (Mesh Adaptive Refinement for Supercomputing) is an open-source, GPU-native mesh and assembly library for finite-element codes on tetrahedral
and hexahedral meshes, running on NVIDIA and AMD GPUs. It provides the distributed mesh and the
building blocks to assemble and solve any finite-element discretization: DOF numbering, sparsity,
owned and ghost rows, halo exchange, and Krylov solvers. The space-filling-curve partitioning and the
element halo come from the cornerstone-octree library. MARS ships ready-made operators for Galerkin
finite elements (linear tetrahedra; high-order matrix-free, experimental) and for the control-volume
finite element method (CVFEM), which its Navier–Stokes solver uses.

The mesh is read, partitioned and stored on the device; the DOF numbering, the sparsity pattern and
the assembly run there too, and MPI connects the GPUs.

This documentation covers **v0.1.1**. The major version is 0, so the API can change between minor
releases.

## Start here

- **[Quickstart](Quickstart.md)**: build MARS, generate a cube mesh and run a GPU assembly.
- **Step by step**, one idea per short program, for newcomers:
  [1. the mesh](step-01-mesh.md), [2. a first GPU loop](step-02-volume.md),
  [3. Poisson's equation](step-03-poisson.md), [4. control volumes](step-04-cvfem.md),
  [5. your own operator](step-05-galerkin.md), [6. the heat equation](step-06-heat.md),
  [7. the lid-driven cavity](step-07-cavity.md), [8. Poiseuille flow](poiseuille_tutorial.md),
  [9. the Taylor–Green vortex](periodic_tgv_tutorial.md),
  [10. high-order matrix-free](Matrix-Free-Tutorial.md) and [adaptive refinement](AMR-Module.md);
  appendix: [your own CMake project](step-appendix-cmake.md).

## Tutorials

- **[Poiseuille channel flow](poiseuille_tutorial.md)**: incompressible Navier–Stokes from the
  mesh to the validated result. Start here if you are new to CFD.
- **[Taylor–Green vortex](periodic_tgv_tutorial.md)**: the same solver on a periodic box, and how
  periodic points stay one unknown on any number of GPUs.
- **[High-order matrix-free operator](Matrix-Free-Tutorial.md)**: the experimental high-order
  CVFEM operator apply, on one GPU and on many.

## Reference

- **[FEM Assembly](FEM-Assembly.md)**: mesh → DOF map → sparsity → assembled CSR, on the GPU.
- **[CVFEM Kernels](CVFEM-Kernels.md)**: the hex assembly kernel variants and how to choose one.

## Status

- **[Known limitations](https://github.com/dganellari/mars/blob/master/KNOWN_LIMITATIONS.md)**:
  what is stable, what is experimental, and what is not supported yet.
- **[Changelog](https://github.com/dganellari/mars/blob/master/CHANGELOG.md)**: what changed in
  each release.

Stable in v0.1: the GPU mesh and assembly pipeline on one and on several ranks, and the
incompressible Navier–Stokes solver on hexahedral meshes (`mars_poiseuille_flow`, `mars_tgv`,
`mars_lid_driven_cavity`). Experimental, among others: the high-order matrix-free operators,
adaptive mesh refinement, and the segregated SIMPLE solver.

## The mesh in code

```cpp
#include "backend/distributed/unstructured/domain.hpp"

// Hex8 mesh, double precision, 64-bit SFC keys, GPU
using Domain = mars::ElementDomain<mars::HexTag, double, uint64_t, cstone::execution::Gpu>;

// Read the mesh, partition it along the SFC and build the cornerstone domain
Domain domain(meshDir, rank, numRanks);

// Built on first access, on the device
const auto& d_owner = domain.getNodeOwnershipMap();          // per local node: 1 owned, 0 ghost
const auto& d_conn  = domain.getElementToNodeConnectivity();  // local node ids, one column per corner

domain.cacheNodeCoordinates();
const auto& d_x = domain.getNodeX();                          // also getNodeY(), getNodeZ()
```

`examples/distributed/unstructured/mars_cvfem_graph.cu` continues from here to an assembled matrix;
the [Quickstart](Quickstart.md) runs it.
