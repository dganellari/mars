# Handoff: Tet collapsed-coordinate (Duffy) sum-factorization — what MARS actually has

## Context
This is for the MARSIR session adding **tetrahedral collapsed-coordinate sum-factorization** to the
MARS-IR codegen. The premise "MARS has a tet sum-factorization in the CUDA implementation" is **false** —
this doc is the self-contained truth of what exists, so you generate against our validated math instead
of chasing production code that isn't there.

---

## 1. The headline answer

| Question | Answer |
|---|---|
| Sum-fac Duffy tet in MARS **production/CUDA**? | **NO.** Tet CUDA = **p=1 linear** CVFEM only (4 nodes, 4×4, 1-pt Gauss, classical FEM stiffness, assembled). |
| Any collapsed/PKD/Duffy/Jacobi-poly in production? | **NO.** All `*_ho_*` sum-fac code is **HEX** (tensor product). `collapse` hits = *periodic-DOF* collapse. `jacobi` hits = Jacobi *preconditioners*. |
| Sum-fac Duffy tet **anywhere**? | **YES — only a host reference** `scratchpad/tet_ho_ref` (rung 9). Host C++, validated bit-exact, **uncommitted, not CUDA, not in the build.** |

**Trap to avoid:** p=1 affine tet is *not* "degenerate collapsed." Affine = constant linear-tet Jacobian;
that is unrelated to the Duffy cube→tet map, and there is **no sum-factorization at p=1** (dense 4×4 is
trivial). Duffy collapse is a **p>1 spectral** construct.

Production evidence: `backend/distributed/unstructured/fem/mars_cvfem_tet_kernel_full.hpp` ("linear
tetrahedra, full 4x4, 1-point Gauss rule, Linear tet stiffness from classical FEM"). The hex sum-fac
template to mimic structurally is `mars_cvfem_ho_apply.hpp` (Knaus Alg-2) + `mars_cvfem_ho_basis.hpp`.

---

## 2. The reference (the oracle)

Location: `scratchpad/tet_ho_ref/` — build any test with `clang++ -std=c++17 -O2 <t>.cpp -o x && ./x`.
It is the **Karniadakis–Sherwin collapsed (Duffy/PKD) construction** — identical to the Nektar++ tet
backward-transform (`ψᵃ/ψᵇ/ψᶜ`, ragged `p+q+r≤D`, nested sweeps). Status of each piece:

| Piece | File | Use as oracle? |
|---|---|---|
| Jacobi P_n^{(a,b)} + deriv | `jacobi.hpp` | ✅ validated |
| Gauss–Jacobi quad (+ deflation) | `gauss_jacobi.hpp` | ✅ validated |
| GLL nodes | `gll.hpp` | ✅ validated |
| Duffy tet quadrature | `tet_quad.hpp` | ✅ validated (monomial moments) |
| PKD collapsed modal basis | `pkd.hpp` | ✅ orthogonal to 1e-12 |
| sum-fac eval / gradient | `test_sumfac.cpp`, `test_sumfac_grad.cpp` | ✅ ==dense to 1e-13 |
| **full sum-fac GALERKIN Laplacian** | **`test_sumfac_op.cpp` (rung 9)** | ✅ **==dense, A·1=0, symmetric (p=2..5). THE oracle.** |
| 1D CVFEM flux balance | `test_cvfem_1d.cpp` | ✅ A·1=0, ==Galerkin@p=1, non-sym high-order |
| **3D collapsed CVFEM** | **`test_cvfem_3d.cpp` (rung 10b)** | ⚠️ fixed (Nanson `det(J_full)` measure): A·1=0 exact, but only APPROXIMATELY consistent (A·x_int ~1e-4→1e-5, shrinks with p) — Gauss can't integrate the rational collapse metric exactly. Research artifact, not the production path. |

**UPDATE — production tet HO code now exists in MARS** (Galerkin sum-fac, the recommended path):
`backend/distributed/unstructured/fem/mars_ho_tet_basis.hpp` (GJ quadrature + PKD tables + nodal
Vandermonde layer), `mars_ho_laplacian_tet.hpp` (host apply == rung 9 bit-exact + CUDA kernels),
`mars_ho_dof_handler_tet.hpp` (barycentric-key C0 DOF matching, no orientation tables). Host-validated
to machine zero (element parity p≤6; C0 assembly gates p≤5). These are the natural codegen targets.

**Generate against rung 9** (the Galerkin sum-fac — what the Nektar++ slide shows, the shared engine).
The CVFEM flux layer on top is unfinished research; do not treat `test_cvfem_3d.cpp` as a reference.

---

## 3. The exact math to match (self-contained)

**Collapse map** (unit tet, a,b,c ∈ [0,1]):  `r=a(1-b)(1-c),  s=b(1-c),  t=c`.

**PKD 1D factors** (argument `x∈[-1,1]`, coord = (1+x)/2):
- `A_p(a) = P_p^{(0,0)}(2a-1)`
- `B_pq(b) = (1-b)^p · P_q^{(2p+1,0)}(2b-1)`      → in code `w=½(1-x)`, `pow(w,p)·jacobi(q,2p+1,0,x)`
- `C_pqr(c)= (1-c)^{p+q} · P_r^{(2p+2q+2,0)}(2c-1)` → `pow(w,p+q)·jacobi(r,2p+2q+2,0,x)`
- Derivatives (chain rule, `d/dcoord = 2·d/dx`): `Ad=2·djacobi(p,0,0,x)`;
  `Bd = -p·w^{p-1}·J + w^p·2·dJ` (p>0 else 0); `Cd` same with exponent `e=p+q`.
- `djacobi(n,a,b,x) = ½(n+a+b+1)·jacobi(n-1,a+1,b+1,x)`.

**Quadrature (Duffy):** `ga=GJ(n,0,0)  gb=GJ(n,1,0)  gc=GJ(n,2,0)`; map [-1,1]→[0,1] with weight
factors `wa=½·ga.w, wb=¼·gb.w, wc=⅛·gc.w` (the `/2,/4,/8` for α=0,1,2). Point weight `wa·wb·wc`.

**Gauss–Jacobi stability gotcha (must replicate):** Newton on `jacobi(n,a,b)` collapses roots at α=2.
Fix = **Newton–Maehly deflation** `dx = f / (fp − f·Σ_{j<i} 1/(xi−x_j))`. Weights
`w_i = 1/((1-x_i²)·djacobi²)`, normalized to `μ0 = 2^{a+b+1}·Γ(a+1)Γ(b+1)/Γ(a+b+2)`.

**Collapse Jacobian** `J_c = ∂(a,b,c)/∂(r,s,t)`, `D1=1-s-t, D2=1-t`:
```
[ 1/D1 , r/D1² , r/D1² ]
[  0   , 1/D2  , s/D2² ]
[  0   ,   0   ,   1   ]
```
Finite at interior GJ nodes (never evaluated at the apex). **Metric** per quad point:
`G_c = W · J_c (B_el⁻¹ B_el⁻ᵀ) J_cᵀ`, `W = wa·wb·wc·detB`.

**Operator = Bᵀ G_c B** via 3 ragged sweeps (contract r→q→p forward; reverse for transpose), one factor
differentiated per gradient component: `(Ad,B,C)→∂/∂a`, `(A,Bd,C)→∂/∂b`, `(A,B,Cd)→∂/∂c`. Metric maps
the 3-vector collapsed gradient to the collapsed flux; transpose sweeps integrate back to modes.

**Validation numbers to reproduce (rung 9):** `max|A_sf−A_dense|` ~2e-15…6e-14, `A·1=0` exactly,
symmetry ~1e-14, for p=2..5 (modes 10/20/35/56).

---

## 4. Findings that shape the codegen

- **Same as Nektar++/Karniadakis–Sherwin.** Their tet backward-transform slide is exactly rung 9's `B`
  leg. We have `B`, `Bᵀ`, and `G_c` — the whole operator.
- **Collapse over-resolves `P_p`.** The N³ collapsed GLL grid gives **more distinct tet nodes than
  dim P_p** (p=2: 15 vs 10; grows to ×3 by p=6). Modal (PKD) is clean; a **nodal/CVFEM** scheme on the
  collapsed grid pays a 2–3× DOF penalty. Relevant if the tet path stays nodal-CVFEM like the hex one.
- **High-order CVFEM is non-symmetric** (Petrov–Galerkin) → needs GMRES, not CG. Galerkin (rung 9) is
  symmetric.
- **Structure to reuse:** collapsed tet = **degenerate hex** — keep the 3-direction sweep skeleton of
  the existing hex emitter, but with **ragged `p+q+r≤D` bounds + the warped `ψᵃ/ψᵇ/ψᶜ` basis + the
  collapse metric**. Far closer to the current hex code than a "face-first, 4-face" rewrite.

---

## 5. Verification for the codegen
Generated tet host code must match **rung 9**: `A·1=0` exactly, symmetry ~1e-14, and (bit-for-bit if you
compile our `test_sumfac_op.cpp`) `A_sf==A_dense`. The Gauss–Jacobi table **must** include the α=2
deflation or high-p tets go wrong.
