# Frozen Tet4 interior replay

Author: GPT/Codex. Date: 2026-09-21. Status: host validation; CUDA execution pending.

The header `fem/segregated/mars_segregated_tet_interior.hpp` implements a new
device-callable interior algebra slice. It is not wired into the existing MARS
Chorin solvers and does not implement the full SIMPLE algorithm.

The governing equations are steady incompressible mass conservation and momentum
balance with physical pressure p [Pa], density rho [kg/m^3], effective dynamic
viscosity mu [Pa s], and symmetric viscous stress mu(grad u + grad u^T).
This slice has no temporal or boundary terms, NSO, compressibility, body force,
moving mesh or rotating frame. The preparer rejects unsupported captured physics.
No relaxation factor is introduced here: supplied componentwise influence
coefficients already carry the reference's selected scaling.

For each oriented median-dual interior sample L->R, S is an area vector [m^2].
Tet4 nodal interpolation N and physical derivatives grad N are supplied for each
of the six samples. The frozen geometry and reconstruction are reference inputs:
this gate does not validate MARS geometry generation or gradient reconstruction.

Pressure volume flux is
`q = sum_j (u_ip,j - D_rhs,ip,j*(sum_n p_n*dN_n/dx_j - (g_L,j+g_R,j)/2))*S_j`.
Density is upwind reconstructed using the sign of q; mass flux is `m=rho_hr*q`.
The residual scatter is `rhs_L -= m`, `rhs_R += m`.
Pressure entries are `A_L,n += -rho_hr*sum_j D_lhs,ip,j*dN_n/dx_j*S_j`,
with the opposite sign at R. Distinct LHS/RHS D values accommodate the reference's
SIMPLEC coefficient choice without claiming that their construction is tested.

Momentum consumes the stored mass flux, uses sign-selected upwind velocity plus
its frozen deferred reconstruction, and linearizes only the upwind velocity.
The implicit stress block at row (L,i), column (n,j) is
`-mu_ip*(delta_ij*grad N_n dot S + dN_n/dx_i*S_j)`; R has the opposite sign.
Its contribution to RHS is minus that block times the current velocity.
The full 12x12 block retains transpose-stress component coupling.

These equations correspond to the pinned OpenAccel 0d69041 routines
`src/assemble/flow/segregatedFlow/{navierStokes,pressureCorrection}/*AssemblerElemTerms.cpp`.
All new computations occur in `mars::segregated::tet_interior`; the CUDA launcher
runs one element per thread. Expected outputs never enter the device input
buffer. Test-only file ingestion and the final result copy use the host;
production mesh assembly, CSR scatter, halos and performance remain future work.

`prepare.py` validates the existing schema2 captures by global node/sample IDs
and serializes them for a dependency-free C++/CUDA executable. Pressure mass flux
is recomputed. Momentum flux comparison is input-preservation only, since that
stage deliberately consumes stored flux. Matrices and RHS are recomputed for
both stages. Matrix/RHS/flux are scaled separately per block, with the existing
comparator tolerance `1e-12*max(1,max_abs_reference_field)`.

Executed locally on the actual 6144-block capture: all 577536 scalar comparisons
pass in the host build, with zero observed difference. Another 234 checks cover
affine-pressure cancellation, frozen Jacobians, SIMPLEC LHS/RHS separation,
local conservation, both upwind signs and an explicit off-component stress
entry. These host checks do not certify CUDA execution.

Build and run through the normal MARS checkout. From the existing CUDA build
directory on Daint, after pulling the commit:

```bash
cmake -S .. -B .
cmake --build . --target mars_segregated_replay mars_segregated_algebra_check -j4
bash ../tests/reference/openaccel/replay/run.sh
```

The targets use the CUDA/compiler/architecture settings of the existing MARS
CMake build and do not link the MARS mesh or linear solvers. They require
MARS_ENABLE_CUDA, MARS_ENABLE_UNSTRUCTURED and MARS_ENABLE_FEM_EXAMPLES.
The helper submits one five-minute single-GPU step, reads the named public
capture, saves a fresh directory under the selected build directory, and returns
failure on mismatch. Optional arguments are public capture exports directory
and build directory (default: current directory). The default capture path is
the completed OpenAccel public-inputs run of 2026-09-20, not a private mesh.
No OpenAccel rebuild or flow simulation is required. The former standalone
bundle and its hard-coded compiler paths are superseded by these CMake targets.
