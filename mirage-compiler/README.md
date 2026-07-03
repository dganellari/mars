# mirage — MARS FEM operator IR / codegen

A build-time compiler that turns a **weak-form operator spec** into specialized
GPU (and host) kernels for the MARS high-order matrix-free CVFEM path. You author
the *pointwise* flux once; `mirage` synthesizes the gather → contract → integrate →
scatter Knaus Algorithm-2 sweep around it and emits CUDA (and, later, HIP)
source. This is the libCEED split — you write the QFunction (our `flux`), the
compiler owns the E/B/D structure.

**Separate by design.** `mirage` is its own tool with zero MARS build coupling: it
emits source that MARS compiles. MARS's build closure never pulls in the
compiler (or, in a later stage, LLVM). Only the *generated* `.cuh` enter MARS,
under `backend/distributed/unstructured/mirage/generated/`.

## Pipeline

```
specs/*.op  --frontend-->  Operator  --synthesize-->  ElementApply (IR)
                                                          |
                              +---------------------------+---------------------------+
                              v                           v                           v
                        backends.ir (dump)        backends.cuda (.cuh)        backends.host_cpp (.hpp)
```

- `mirage/expr.py` — arithmetic expression AST + parser + C renderer (the whole "language").
- `mirage/frontend.py` — parse a `.op` spec into an `Operator`.
- `mirage/ir.py` — synthesize the `ElementApply` schedule; record which tangential
  derivatives / metric components the authored flux actually uses (so unused work
  is pruned).
- `mirage/backends/cuda.py` — emit the device apply, structurally a twin of
  `mars::fem::ho_cvfem_apply_kernel` (PerPoint, FP64) with the flux inlined.
- `mirage/backends/host_cpp.py` — emit a host-C++ twin for local numeric validation.

## Usage

```bash
# IR dump (Stage-1 view: what the compiler will emit)
python3 mirage-emit.py specs/laplacian.op --backend ir

# CUDA header (the deliverable MARS compiles)
python3 mirage-emit.py specs/laplacian.op --backend cuda \
    -o ../backend/distributed/unstructured/mirage/generated/mirage_laplacian_apply.cuh

# Host twin (used by the local validation)
python3 mirage-emit.py specs/laplacian.op --backend host -o generated/laplacian_apply_host.hpp
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
python3 mirage-emit.py specs/laplacian.op --backend ir        # shows the symbolic ∂flux/∂x + J·du
python3 mirage-emit.py specs/nldemo.op   --backend host-jvp   # emits the JVP host apply
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
`examples/distributed/unstructured/mars_mirage_ho_apply.cu`, built with
`-DMARS_ENABLE_MIRAGE=ON` (requires UNSTRUCTURED + CUDA). It runs the generated
kernel and the hand kernel on the same mesh / elemDof / metric `d_G` and diffs
them on device — the go/no-go that the emitted CUDA compiles and runs
bit-identically to the validated kernel.

## Roadmap

- **Stage 1–2 (done):** operator spec → IR → CUDA/host source; local host parity;
  symbolic form differentiation → generated Jacobian-action (JVP) operator.
- **Next:** CUDA parity gate on Alps; order/metric sweep (p=1..8, PerPoint/Affine);
  HIP backend emit + gfx942 compile smoke; perf non-regression vs the hand kernel.
- **Stage 2 backend (deferred):** an MLIR/LLVM lowering path (EmitC first, then
  nvvm/rocdl) added behind the same front-end + IR — no changes to `.op` specs.
