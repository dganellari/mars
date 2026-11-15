# Poisson peak with both AMR policies

`mars_amr_galerkin` compares local refinement with uniform refinement on a public
unit-cube mesh. A narrow solution peak makes extra resolution useful in one part
of the domain. Conforming LEPP and nonconforming red refinement share the same
Poisson assembly, boundary treatment and solver.

## Problem

On the dimensionless cube $\Omega=[0,1]^3$, solve the steady equation

$$
-\Delta u=f,\qquad u|_{\partial\Omega}=u_{\rm exact},\qquad
u_{\rm exact}(x)=\exp(-\alpha\lVert x-c\rVert^2),
$$

where $c=(0.3,0.4,0.45)$, $\alpha=100$, and

$$
f(x)=(6\alpha-4\alpha^2\lVert x-c\rVert^2)u_{\rm exact}(x).
$$

All runs use P1 tetrahedra, four-point tetrahedral quadrature, and PCG with
BoomerAMG at relative tolerance `1e-10`. The printed `l2` uses that same quadrature;
it is a numerical estimate of the error norm. Hex8 with Q1 and eight-point Gauss
quadrature is also supported by red refinement. There is no time integration.

With red refinement, hanging-node interpolation defines `P` and the solver uses
`P^T A P`. LEPP closes shared edge stars, so all nodes are independent geometric
nodes; `DofSpace` still handles distributed ownership and copies. Both paths
rebuild and redistribute through Cornerstone after refinement.

The marker selects element centroids within radius `0.25` of the known peak.
This demonstrates local resolution; it is not a solution-based error estimator.
`--uniform` instead marks every owned element. For LEPP, marking all tets still
means conforming edge bisection, not the eight-child split of red refinement.

## Run the comparison

Build MARS with its CUDA or HIP unstructured backend, `MARS_ENABLE_HYPRE=ON` and
`MARS_ENABLE_FEM_EXAMPLES=ON`. Use the public `cube16_tet.mesh` binary fixture,
which has unit-cube coordinates and one tet block. The example identifies its
Dirichlet boundary from the box; arbitrary geometries need their own boundary
selection. The mesh is a directory of binary arrays, not an Exodus filename.

From a build directory, use your MPI launcher with one rank per GPU. For example:

```bash
mpirun -n 4 ./examples/distributed/unstructured/mars_amr_galerkin \
  --mesh=/path/to/cube16_tet.mesh --element=tet --problem=peak \
  --refinement=red --levels=2 --audit-algebra

mpirun -n 4 ./examples/distributed/unstructured/mars_amr_galerkin \
  --mesh=/path/to/cube16_tet.mesh --element=tet --problem=peak \
  --refinement=lepp --levels=4 --audit-algebra

mpirun -n 4 ./examples/distributed/unstructured/mars_amr_galerkin \
  --mesh=/path/to/cube16_tet.mesh --element=tet --problem=peak \
  --refinement=red --uniform --levels=2 --audit-algebra

mpirun -n 4 ./examples/distributed/unstructured/mars_amr_galerkin \
  --mesh=/path/to/cube16_tet.mesh --element=tet --problem=peak \
  --refinement=lepp --uniform --levels=4 --audit-algebra
```

On Slurm, use the site's `srun --export=ALL` environment and placement options
instead of nesting `mpirun`. Installed executables live in `bin/examples/`.
The provided cluster gate builds a consumer against a relocated installation and
runs these comparisons on NVIDIA and AMD at one and eight ranks.

Keep `--alpha`, `--radius`, the initial mesh and solver tolerance fixed. Each
`[amr-galerkin] level=...` line reports global elements, independent DOFs, hanging
nodes, solver residual and solution errors. Compare `l2` **against `dofs`**, not
against refinement rounds: red refinement and LEPP do different work per round.
A useful result reaches a chosen error with fewer DOFs than its uniform control.
A mesh picture or lower residual alone does not establish that benefit.

`--problem=linear` supplies the existing affine correctness control. It should
pass on both policies before interpreting peak errors. `--audit-algebra` checks
prolongation, the assembled/composed operator and the true solve residual.
`--audit-mesh` is a red-tet diagnostic and cannot be combined with LEPP.
`--uniform` uses the default maximum-relative marking policy; it is not combined
with `--marking=octree`, whose leaf aggregation would change the control.

## Current scope

The underlying conforming and nonconforming implementations have separate native
correctness evidence in [AMR](AMR-Module.md) and
[conforming Tet4 AMR](Conforming-Tet-AMR.md). The combined peak comparison is new:
CUDA compatibility syntax checks pass locally, while native accuracy-versus-DOF
results are pending. No measured efficiency or scaling claim is made here.
