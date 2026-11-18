# Imported Hex8 BP3 comparison

Build the existing MFEM 4.8 benchmark with `-DBP3_MESH_TOOLS=ON`, pointing
`BP3_NETCDF_INCLUDE_DIR` and `BP3_NETCDF_LIBRARY` at NetCDF-C when needed.
The library itself does not need MFEM's optional NetCDF support.

```
mfem_partition wing.e wing8. 8
mpiexec -n 8 mfem_bp3 --mesh-prefix wing8. --deform 0 --p 3 \
  --pc mg --cheb-order 2 --amg-cycles 3 --amg-settings mars \
  --rhs discrete --tol 1e-13 --true-tol 1e-10 --maxit 2000 --solves 3
```

The offline tool reads one complete mesh on one host, reconstructs exterior
faces and assigns contiguous intervals of MFEM’s Hilbert ordering to ranks. Rank element counts differ by at most one. It writes MFEM parallel
mesh files, including shared topology, and a small partition-count receipt.
Solver ranks then read only their own part. Preprocessing holds the global mesh
and is not a scalable parallel Exodus reader; report its wall time and peak host
memory separately. Existing output files are rejected.

The reader preserves file corner order and double coordinates. It accepts
multiple active Hex8 blocks, split or packed coordinates, and integer node IDs.
It rejects curved elements, inverted corner orientation and malformed IDs.
All physical exterior faces receive homogeneous Dirichlet conditions; side-set
names are intentionally irrelevant for this benchmark. Geometric refinement,
geometry deformation and continuum forcing are disabled on imported meshes.

The scalar problem is unit-coefficient diffusion with continuous Q3 unknowns,
Q1 geometry and five-point Gauss quadrature in each direction. Its discrete
right-hand side is A times the sampled sine field, zeroed on all exterior nodes.
The hierarchy is p=3 to p=1, with order-two Chebyshev/Jacobi smoothing and three
BoomerAMG cycles using the existing benchmark settings. The MARS comparison
uses Jacobi, MFEM's Chebyshev interval convention and no coarse consolidation.
MARS uses flexible CG and MFEM uses CG. They construct their Hilbert partitions independently;
whole-solve differences must not be attributed solely to the compiler.

`--audit-prefix` exports coordinates, the essential mask, sampled field, RHS,
probe and operator action, and solution. It is limited to 4096 imported elements.
CPU checks compare these fields against the procedural warped cube on 1/2/4/8
ranks. GPU execution and full-wing capacity need separate native checks.

Run the reusable CPU check with a CPU/MPI MFEM build:

```
python3 test_import.py --build /path/to/driver-build
```

It needs NumPy and netCDF4 only to create the synthetic fixture and compare
canonical fields. The production partitioner and cluster kit need neither.
