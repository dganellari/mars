# marsir-compiler — MARSIR Stage 1 (operator spec → source)

A build-time compiler that turns a **weak-form operator spec** into specialized
GPU (and host) kernels for the MARS high-order matrix-free CVFEM path. You author
the *pointwise* flux once; MARSIR synthesizes the gather → contract → integrate →
scatter Knaus Algorithm-2 sweep around it and emits CUDA/HIP source. This is the
libCEED split — you write the QFunction (our `flux`), the compiler owns the E/B/D
structure. (Naming: the project is MARSIR — formerly mirage, formerly mir; the
Stage-2 MLIR dialect in `marsir-mlir/` is still named `mir`.)

**Separate by design.** `marsir` is its own tool with zero MARS build coupling: it
emits source that MARS compiles. MARS's build closure never pulls in the
compiler (or, in a later stage, LLVM). Only the *generated* `.cuh` enter MARS,
under `backend/distributed/unstructured/marsir/generated/`.

## Pipeline

```
specs/*.op  --frontend-->  Operator  --synthesize-->  ElementApply (IR)
                                                          |
        +----------------+----------------+---------------+----------------+
        v                v                v               v                v
   ir (dump)      cuda (.cuh)      host_cpp (.hpp)   tet_galerkin      mlir_ir (mir IR,
                                                     (collapsed tet)   the Stage-2 bridge)
```

- `marsir/expr.py` — arithmetic expression AST + parser + C renderer +
  `diff()`/`simplify()` (the whole "language").
- `marsir/frontend.py` — parse a `.op` spec into an `Operator`.
- `marsir/ir.py` — synthesize the `ElementApply` schedule; record which tangential
  derivatives / metric components the authored flux actually uses (so unused work
  is pruned).
- `marsir/backends/cuda.py` — emit the device apply, structurally a twin of
  `mars::fem::ho_cvfem_apply_kernel` (PerPoint, FP64) with the flux inlined,
  plus the JVP kernel and the HIP macro veneer (one header, nvcc + hipcc).
- `marsir/backends/host_cpp.py` — emit a host-C++ twin for local numeric validation.
- `marsir/backends/tet_galerkin.py` — the collapsed-coordinate (Duffy/PKD) tet
  Galerkin apply (see `TET_SUMFAC_HANDOFF.md`).
- `marsir/backends/mlir_ir.py` — emit the operator as textual `mir` dialect IR
  for `marsir-mlir/`: `emit` (the PA chain, `--backend mlir`), `emit_full` (the
  whole Knaus apply, `--backend mlir-full`), `emit_full_batched` (the fused
  single-kernel form, `--backend mlir-full-batched`). One `.op` feeds both stages.

## Usage

```bash
# IR dump (Stage-1 view: what the compiler will emit)
python3 marsir-emit.py specs/laplacian.op --backend ir

# CUDA header (the deliverable MARS compiles)
python3 marsir-emit.py specs/laplacian.op --backend cuda \
    -o ../backend/distributed/unstructured/marsir/generated/marsir_laplacian_apply.cuh

# Host twin (used by the local validation)
python3 marsir-emit.py specs/laplacian.op --backend host -o generated/laplacian_apply_host.hpp

# mir dialect IR (Stage 2; pipe into marsir-mlir/build/tools/mir-opt/mir-opt)
python3 marsir-emit.py specs/laplacian.op --backend mlir               # PA chain
python3 marsir-emit.py specs/laplacian.op --backend mlir-full          # full Knaus apply
python3 marsir-emit.py specs/laplacian.op --backend mlir-full-batched  # fused kernel form
```

## Operator spec

```
operator laplacian {
    element = hex
    form    = diffusion               # selects the Knaus Alg-2 skeleton
    flux    = g2*deriv + g0*dt2 + g1*dt1
}
```

At each subcontrol-surface quadrature slot the flux may reference: `deriv`
(normal derivative Dtil·u), `dt1`/`dt2` (tangential derivatives D·interp), and
`g0`/`g1`/`g2` (metric: tang2 r-axis, tang1 s-axis, normal). Change the flux and
the whole kernel regenerates — e.g. `specs/laplacian_orthocube.op` uses
`flux = g2*deriv`, and the compiler drops both tangential contractions and the
interp buffer automatically.

## Local validation (no GPU)

```bash
python3 tests/validate_host.py
```

Emits the host twin, compiles it against the hand-written
`mars::fem::applyHoCvfemElement` oracle, and diffs the two on a distorted hex for
p=1..8. Current result: **bit-identical** (`parity_rel = 0`) at every order, with
the constant-nullspace invariant at machine precision. This proves the IR
captured the operator; the CUDA emission shares the same IR and flux rendering.

## Symbolic form differentiation → Jacobian / adjoint operators

The flux AST is closed (`Num/Var/Unary/Binary`), so it is closed under
differentiation. `Expr.diff()` + `simplify()` give **symbolic form
differentiation** (the same idea as UFL's `derivative()` in Firedrake/FEniCS,
not runtime autodiff): from one authored residual flux, the compiler derives the
**Jacobian-action (JVP) operator** — `J(u)·du = Σ_x (∂flux/∂x)·x_dot` over the
u-linear inputs — which rides the identical gather→contract→integrate→scatter
skeleton. This is the matrix-free tangent operator a Newton/FSI solve needs,
generated, not hand-written.

```bash
python3 marsir-emit.py specs/laplacian.op --backend ir        # shows the symbolic ∂flux/∂x + J·du
python3 marsir-emit.py specs/nldemo.op   --backend host-jvp   # emits the JVP host apply
python3 tests/validate_jvp.py                                 # numeric proof, p=1..8
```

`validate_jvp.py` checks two things: (A) for the **linear** Laplacian, `J(u)·du`
equals the primal `A·du` to reduction order (a linear operator is its own
Jacobian); (B) for a **nonlinear** demo flux, `J(u)·du` matches a central
finite difference `(F(u+εδu)−F(u−εδu))/2ε` — the independent proof that
`∂flux/∂x` is correct. Both PASS. Stays pure-Python/LLVM-free; real autodiff
(Enzyme) is reserved only for non-smooth fluxes (limiters, local sub-solves)
that symbolic diff can't reach.

## MARS integration + CUDA parity gate (Alps)

The committed generated header is consumed by the parity example
`examples/distributed/unstructured/mars_marsir_ho_apply.cu`, built with
`-DMARS_ENABLE_MARSIR=ON` (requires UNSTRUCTURED + CUDA). It runs the generated
kernel and the hand kernel on the same mesh / elemDof / metric `d_G` and diffs
them on device. **Validated on Alps GH200** (single rank, 16.7M elements):
generated primal == hand kernel at `gen_vs_hand_rel = 3.357e-16`; the JVP gate
`mars_marsir_jvp_test` gives `J_vs_FD_rel = 1.633e-10` vs central finite
differences on the nonlinear demo. Host sweeps p=1..8 are bit-identical
(`parity_rel = 0`).

## Status / roadmap

- **Done:** `.op` → IR → CUDA/HIP/host source; host parity bit-identical p=1..8;
  symbolic JVP (linear J==A ~2e-16, nonlinear vs FD ~1e-10); GH200 CUDA parity;
  collapsed-tet Galerkin backend (== dense oracle p=2..5); the mlir backends
  feeding Stage 2 (all six Stage-2 execution gates pass on their output).
- **Next:** a CUDA/HIP tet kernel from the tet backend; a proper `element = tet`
  spec with authored flux + tet symbolic diff; persist the tet reference basis
  library (`validate_tet.py` still needs `MARSIR_TET_ORACLE`); wire the generated
  JVP as the Newton matvec; the launch-parameter autotuner.
- **Stage 2 (real, alive):** the `mir` MLIR dialect in `marsir-mlir/` — lowering
  passes, FP64 DMMA tensor-core schedule, warp distribution. See
  `marsir-mlir/README.md` and `internal-notes/marsir_full_tutorial.md`.
