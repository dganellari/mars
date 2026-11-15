# MARS 0.2.0 validation

Version 0.2.0 combines the optional structured Kokkos backend with conforming
Tet4 LEPP and nonconforming Tet4/Hex8 adaptive refinement on Cornerstone.
The release checks below passed; the table distinguishes package validation
from the main-build CTest suites.

| Check | NVIDIA | AMD | Scope |
| --- | --- | --- | --- |
| Relocated structured package | 5358530: 28 checks | 743609: 84 checks across CPU and HIP configurations | Package consumers and supported structured capabilities |
| AMR/PDE package | Build 5358530, validation 5358903: 95 commands | 743609: 95 commands | Refinement, migration, boundaries and repeated Escher Poisson solves; 1, 4 and 8 ranks |
| Main-build AMR/mesh release suite | 5360032: 37/37 | 743925: 37/37 | AMR release registration, red-tet orientation and node-halo checks |
| Main-build core release suite | Builds 5360839/5360900, tests 5360935/5360941: 48/48 | Builds 743973/743974, tests 743975/743976: 48/48 | Assembly, numbering, Poisson, flow, tutorials and high-order regression checks; 1 and 4 ranks |

The standalone CPU/MPI AMR release suite passed 33/33 checks; the main-build
CPU AMR selection passed 5/5. The core GPU suite uses two CTest selections,
with 22 and 32 executions: six mesh fixtures run in both, giving 48 distinct
tests. No failed or skipped check is counted as a pass.

The package campaign used revision `5932a47a3`. Main-build AMR registration
used `e5a94667a` on NVIDIA and `3ab940f43` on AMD; the latter fixes the HIP test
architecture option. Both core campaigns used `3ab940f43`. The merge between
these package and main-build revisions did not change the AMR numerical sources.
Final release changes after these checks update documentation and citation
metadata only. Source manifests, command exit codes, raw log hashes and CTest
records were checked. Executable hashes were verified by the running validation
harness; binaries were not independently downloaded for a second hash check.

The initial NVIDIA package allocation completed its structured checks and
native build but timed out before running AMR. The continuation reused that
verified build. Earlier core-kit setup failures occurred before compilation;
the successful jobs above include both build stages and both test stages.

## Supported scope

AMR takes conforming, single-block, nonperiodic input. It has no coarsening or
refinement-history restart. Nonconforming boundary-label propagation is supported
for tetrahedra. Coordinates travel as node data; SFC keys continue to define
topology, partitioning, tree construction, communication and node identity.

The complete adaptive PDE example is a manufactured Poisson problem on the
public MFEM Escher tet mesh. Its indicator uses the known exact solution.
It validates refinement and repeated solves, not an error estimator for unknown
physics. Public cube, Escher, graded and thin fixtures establish the documented
scope; they do not establish arbitrary-geometry accuracy or curved CAD recovery.
Cross-platform and rank comparisons require identical mesh/DOF/hanging counts
and bounded solution-error agreement, not identical solver iterations.

No AMR scaling or performance advantage is claimed. Structured timing against
the older snapshot is reported information rather than a release gate.
The serial backend is unchanged. See [known limitations](https://github.com/dganellari/mars/blob/v0.2.0/KNOWN_LIMITATIONS.md),
[AMR](AMR-Module.md), [AMR checks](https://github.com/dganellari/mars/blob/v0.2.0/tests/amr/README.md) and
[structured Kokkos](Structured-Kokkos.md) for details.
