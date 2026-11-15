# Poisson peak with both AMR policies

`mars_amr_galerkin` compares local refinement with uniform refinement on a public
unit-cube mesh. A narrow solution peak makes extra resolution useful in one part
of the domain. Conforming LEPP and nonconforming red refinement share the same
Poisson assembly, boundary treatment and solver.

## Reproduce the v0.2.0 example

This walkthrough uses the exact public Exodus input from the accepted release
runs, shipped as `tests/amr/data/mfem/escher-refined.e`. It has 2688 tetrahedra,
665 nodes and labelled exterior faces. Its provenance and checksum are in the
[fixture README](https://github.com/dganellari/mars/blob/v0.2.0/tests/amr/data/mfem/README.md).
Run these commands from the repository root in your compiler/MPI environment.

You need C++20, MPI, CUDA or HIP, Hypre built for that GPU backend, and netCDF-C
for Exodus input. `HYPRE_DIR` and `NETCDF_ROOT` identify those installations
when they are not on the default search path. For NVIDIA Hopper:

```bash
cmake -S . -B build-amr -DCMAKE_BUILD_TYPE=Release \
  -DMARS_ENABLE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=90 \
  -DMARS_ENABLE_UNSTRUCTURED=ON -DMARS_ENABLE_FEM_EXAMPLES=ON \
  -DMARS_ENABLE_HYPRE=ON -DMARS_ENABLE_TESTS=ON
cmake --build build-amr --target mars_amr_galerkin --parallel 2
```

For AMD MI300, replace the CUDA options with `-DMARS_ENABLE_HIP=ON` and
`-DCMAKE_HIP_ARCHITECTURES=gfx942`, using the ROCm toolchain and HIP-enabled Hypre.
Select the architecture for your hardware. Check that CMake found netCDF before
running an Exodus case; the main [build guide](https://github.com/dganellari/mars#downloading-mars-and-its-dependencies)
covers dependency configuration.

First run an affine correctness control on one rank per GPU:

```bash
mpirun -n 4 ./build-amr/examples/distributed/unstructured/mars_amr_galerkin \
  --mesh=tests/amr/data/mfem/escher-refined.e --element=tet --boundary=side-sets \
  --refinement=lepp --bucket-size=8 --problem=linear --uniform \
  --levels=4 --tol=1e-10 --audit-algebra
```

Then solve the same quadratic manufactured Poisson problem with each policy:

```bash
mpirun -n 4 ./build-amr/examples/distributed/unstructured/mars_amr_galerkin \
  --mesh=tests/amr/data/mfem/escher-refined.e --element=tet --boundary=side-sets \
  --refinement=lepp --bucket-size=8 --problem=quadratic --error-indicator \
  --levels=4 --tol=1e-10 --audit-algebra

mpirun -n 4 ./build-amr/examples/distributed/unstructured/mars_amr_galerkin \
  --mesh=tests/amr/data/mfem/escher-refined.e --element=tet --boundary=side-sets \
  --refinement=red --bucket-size=8 --problem=quadratic --error-indicator \
  --levels=6 --tol=1e-10 --audit-algebra
```

Use `-n 1` for a single GPU; the accepted campaign also covers eight ranks.
On a cluster, run inside an allocation and use the site's MPI/GPU placement
recipe. With Slurm, use `srun --export=ALL` and the required environment and
binding flags instead of nesting `mpirun`. GPU-aware MPI is required by the
underlying distributed solver; `--staged` only stages the AMR sparse messages.

Each level reports the solve before the next refinement, so `--levels=4`
performs the initial solve plus four adapt-and-solve rounds. Watch `elements`,
`dofs`, `hanging`, `l2`, the algebra audit and the volume/boundary-area checks.
The accepted adaptive red run ends with 27328 elements, 4637 independent DOFs,
1541 hanging nodes and L2 error about 0.0380, down from about 0.1275 initially.
LEPP keeps zero hanging nodes. Compare error against independent DOFs; the
methods use different numbers and kinds of refinement steps, so this is not
an efficiency comparison. Exact-solution marking is a validation tool.

For the registered AMR regression suite, build all test targets before CTest:

```bash
cmake --build build-amr --parallel 2
ctest --test-dir build-amr -L amr --output-on-failure
ctest --test-dir build-amr -L release --output-on-failure
```

The release label includes AMR alongside the existing checks. Optional checks
that need additional fixtures retain their documented inputs. The complete
Escher PDE campaign is a separate package gate, not implied by the CTest label.
See [release validation](Release-0.2.0.md) for the accepted jobs and limits.

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
the labelled Escher Tet4 comparison below passes native NVIDIA/AMD checks,
including error versus independent DOFs. This does not validate every peak-driver
configuration. No measured efficiency or scaling claim is made here.


## Labelled non-box tetrahedra

For an initially conforming, single-block, nonperiodic Tet4 Exodus mesh whose
side sets cover the entire exterior surface, use `--boundary=side-sets` with
either `--refinement=lepp` or `--refinement=red`. Every exterior face must be labelled
and no interior face may be labelled. The initial surface is checked by exact
face incidence. Both methods inherit labels by face ancestry through refinement
and Cornerstone redistribution. Red refinement replaces each marked boundary
triangle with its four children, using the refiner's actual midpoint indices.
The example unions the side sets into one Dirichlet surface and routes membership
to the independent DOF owners through the constraint restriction.

The driver checks domain volume and each labelled surface area after every adapt.
LEPP also repeats the full face-incidence check. That conforming incidence check
is not applied after red refinement: an interior coarse face can meet four fine
faces. Red surface coverage instead follows the tested parent-face inheritance;
area conservation alone is not a proof of complete coverage.

The equation is `-laplace(u) = f`, with unit diffusion and exact Dirichlet data on
that surface. `--problem=linear` remains the affine consistency check.
`--problem=quadratic` uses `u=x*x+y*y+z*z` and `f=-6`; its squared P1 error is
integrated with a degree-four-exact tetrahedral rule. `--error-indicator` marks
from each owned element's integrated squared error against this manufactured
solution, using the configured maximum-relative threshold. It is an oracle for
validation, not an estimator available for an unknown physical solution.
`--uniform` supplies the comparison with all elements marked.

```bash
./mars_amr_galerkin --mesh=escher.e --element=tet --boundary=side-sets \
  --refinement=lepp --bucket-size=8 --problem=linear --uniform --levels=4 --audit-algebra
./mars_amr_galerkin --mesh=escher.e --element=tet --boundary=side-sets \
  --refinement=lepp --bucket-size=8 --problem=quadratic --error-indicator --levels=4 --audit-algebra
./mars_amr_galerkin --mesh=escher.e --element=tet --boundary=side-sets \
  --refinement=red --bucket-size=8 --problem=quadratic --error-indicator --levels=6 --audit-algebra
```

Compare the same manufactured problem and geometry, with uniform controls for
both policies. Use error versus independent DOFs, not error versus adaptation
rounds. `[amr-pde-wall]` reports cumulative wall time from initialization through
the current solve, including all previous solves/adaptations and enabled checks.
It takes the maximum over ranks and identifies audited runs. These times are
useful end-to-end costs, but audited validation runs do not establish a speedup.

The public MFEM Escher fixture can be converted with
`tests/amr/write_public_tet_exodus.cpp` (link netCDF): provide its input `.mesh`
and output `.e` paths. The converter normalizes tetrahedral orientation and keeps
the supplied boundary labels. No private geometry is needed.

Local checks validate exterior-face coverage on 1/4/8 MPI ranks, rejection of
missing/interior labels, quartic error integration and an independent dense CPU
P1 solve on four host LEPP rounds. The old box-boundary selection fails the affine
control on this geometry. Red boundary inheritance and redistribution also pass
on 1/4/8 local ranks, including a missing-midpoint rejection. An independent
constrained CPU solve on two 1-irregular host fixtures preserves affine fields
and reduces quadratic L2 error from 1.18 to 0.359, with nonzero hanging nodes.
Native NVIDIA and AMD validation now passes all 36 configurations on the same
2688-tet Escher fixture: two refinement methods, affine/uniform/adaptive solves,
1/4/8 ranks and direct/staged AMR messages. LEPP uses four rounds; red uses three
for affine/uniform controls and six for adaptive refinement. Six adaptive red
rounds reduce L2 error from 0.1274804 to 0.0380450, with 4637 independent DOFs and
1541 hanging nodes. Counts agree across ranks/platforms; quadratic errors agree
within relative 1e-5. Boundary area and volume agree with the input file to
floating-point precision. The unchanged gate requires at least 25% adaptive
error reduction; the earlier three-round red run reduced it only 18.68% and
remains a recorded failed control.

The native fixture is the public 42-tet Escher mesh uniformly red-refined twice
before conversion to Exodus. Use sufficient initial elements for the selected
rank count: a rebuilt partition with no owned cells is rejected. The commands
above accept a converted input mesh, but their accuracy is not guaranteed for
every input or number of rounds. Full job provenance and scope are in the
[validation record](https://github.com/dganellari/mars/blob/v0.2.0/tests/amr/README.md#physical-geometry-and-non-box-pde-validation).
