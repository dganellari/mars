# Distributed High-Order Matrix-Free CVFEM in MARS — A Tutorial

## Who this is for

You are an HPC or scientific-computing engineer who knows some C++ and a little about
finite elements, and you want to understand — from the ground up — how MARS computes the
single most important operation in any PDE solver, `y = A·x`, *without ever storing the
matrix*, at high polynomial order, on thousands of GPUs. We start from "what is a matvec"
and end at the honest, in-progress trillion-degree-of-freedom frontier. No prior knowledge
of high-order methods, space-filling curves, or GPU halos is assumed; we build each idea
before we use it.

## Learning path

The tutorial climbs in difficulty. Each chapter assumes the one before it.

1. **Why** — the problem (PDEs → a giant `A·x`) and the payoff (high order + matrix-free).
   Pure intuition and the headline numbers. *Beginner.*
2. **The matrix-free operator on a single GPU** — what a matvec really is, sum-factorization,
   the GPU kernel, and how we prove it correct. *Beginner → intermediate.*
3. **Distributed DOF numbering** — naming a degree of freedom the same way on every rank,
   and deciding who owns it. The canonical key and a teachable ownership bug. *Intermediate.*
4. **The distributed matvec and the halo** — the `forward → apply → reverse-add` sandwich,
   receiver-driven exchange, GPUDirect, and the `A·1 = 0` oracle. *Intermediate → advanced.*
5. **Scaling to billions** — weak scaling, why communication self-resolves, the real memory
   wall, and a GPU-native rewrite of the host bottleneck. *Advanced.*
6. **The trillion frontier (honest, in-progress)** — the arithmetic that fits, the crash we
   hit, the expert lesson in *not* widening an integer, and what is genuinely not done yet.
   *Expert.*

A glossary and a "try it yourself" section with real run commands follow the chapters.

---

## Chapter 1 — Why: the problem and the payoff

### The one-paragraph picture of FEM/CVFEM

Almost every simulation in engineering — how heat spreads through a turbine blade, how air
flows over a wing, how a pump pushes fluid — is a *partial differential equation* (PDE) we
cannot solve with pen and paper. The **finite element method (FEM)** turns that continuous
problem into a giant system of linear equations by chopping the domain into small pieces
("elements" — tetrahedra, hexahedra) and approximating the unknown field (temperature,
velocity, pressure) by simple polynomials inside each piece. **CVFEM** (control-volume FEM)
is a close cousin tuned for fluid flow: it enforces conservation (mass, momentum) on small
control volumes built around each node, which is exactly what you want when you must not
lose or invent mass. The end product is always the same shape: a huge matrix `A` and a
vector, and the bulk of the compute is doing `y = A·x` over and over inside a linear solver.

So the whole game reduces to one question: **how do we apply `A` to a vector, fast, for an
`A` so big it spans thousands of GPUs?**

### What "high order" (`p`) means

Inside each element, the polynomial has a degree `p`. At `p=1` the field is linear
(flat-faced, like a coarse mosaic); at `p=4` it is a degree-4 polynomial that can curve and
bend within a single element. Higher `p` means each element carries far more degrees of
freedom (DOF) — the numbers we actually solve for. A `p=4` hex holds `(p+1)^3 = 125` DOF
instead of `8`.

Why pay for that? Because accuracy improves much faster than cost. High order gives roughly
**`p+1` order of convergence** — each refinement cuts the error by a steep power. The payoff
is measured directly in MARS: at equal DOF count, going from `p=1` to `p=4` is about **570×
more accurate**. That is the convergence story.

**Code jump — the convergence data behind `fig_convergence`:**
`docs/figures/make_matfree_figs.py`, the `conv` dictionary:
```python
conv = {
    1: ([125, 343, 729, 2197], [1.68e-1, 7.86e-2, 4.50e-2, 2.03e-2], ...),
    4: ([125, 729, 2197],      [2.01e-3, 2.62e-4, 3.47e-5],          ...),
}
```
Notice: read the *last column* of each row — at ~2197 DOF, `p=1` error is `2.0e-2` but `p=4`
is `3.5e-5`. Same number of unknowns, ~570× lower error. High order buys accuracy you cannot
get from just adding more low-order elements.

<img src="figures/fig_convergence.png" width="460" alt="L2 error vs DOF for p=1..4">

*L2 error vs DOF, one line per order. The slopes steepen with `p` (≈ `p+1` order). At ~2197
DOF the dotted line shows `p=1` and `p=4` differ by ~570× in error for the same unknown count.*

One thing to be clear about: "equal DOF count" means *different meshes on the same unit cube*. A
`p=1` hex contributes ~1 unknown, a `p=4` hex contributes ~64 (`(p+1)^3` nodes). So matching the
unknown count puts `p=1` on a fine mesh and `p=4` on a coarse one — at the ~2197-DOF point that is
a `12^3 = 1728`-element mesh (`p=1`) versus a `3^3 = 27`-element mesh (`p=4`). The fair question is
"for a fixed number of unknowns, which order resolves the solution better," and high order wins
because the error decays like ~`h^(p+1)`. One caveat: this `h^(p+1)` rate assumes a *smooth*
solution and accurate geometry — on a real mesh with corner singularities, boundary layers,
turbulence, or curved walls approximated by straight-sided hexes, the gain is far smaller than 570×.

### What "matrix-free" means, and why it is the whole point at high order

The obvious way to do `y = A·x` is to *build and store* the matrix `A` first (in CSR/sparse
format), then multiply. That is **assembled**. The problem: at high order, `A` becomes
enormous *per DOF*. A `p=4` element couples all 125 of its DOF to each other, so the stored
matrix balloons.

**Matrix-free** flips this: never store `A`. Instead, recompute `A·x` on the fly, element by
element, straight from the basis polynomials and the element geometry. You trade a little
extra arithmetic for an enormous saving in memory — and memory (HBM bandwidth and capacity)
is the real wall on a GPU.

**Code jump — the memory crossover behind `fig_memory`:**
`docs/figures/make_matfree_figs.py`:
```python
asm = [324, 1500, 4116, 8748, 15972, 26364, 40500]   # assembled bytes/DOF, p=1..7
mf  = [421,  221,  169,  147,   134,   126,   121]    # matrix-free bytes/DOF, p=1..7
...
ax.annotate("336$\\times$ at $p{=}7$", ...)
```
Notice two things. (1) The assembled cost *grows* with `p` (324 → 40500 bytes/DOF) because
the coupling explodes; the matrix-free cost *shrinks* (421 → 121) because you amortize
geometry over more DOF. (2) At `p=1` matrix-free is actually slightly *worse* (421 vs 324 —
that's the 0.8× crossover), but by `p=7` it is **336× leaner**. High order is exactly where
matrix-free pays off.

<img src="figures/fig_memory.png" width="460" alt="Operator memory footprint vs p">

*Bytes/DOF vs `p` (log scale). Assembled storage climbs 324 → 40500 as coupling explodes;
matrix-free falls 421 → 121. The lines cross just above `p=1`, and by `p=7` matrix-free is
336× leaner.*

There is a second, less obvious win. Done naively, recomputing `A·x` for a `p=4` element
would be a dense `125×125` multiply. **Sum-factorization** restructures it into a sequence of
small 1D operations along each axis, collapsing the cost so that throughput stays flat across
order.

**Code jump — the sum-factorized element apply:**
`backend/distributed/unstructured/fem/mars_cvfem_ho_apply.hpp`, inside `applyHoCvfemElement`:
```cpp
for (int q=0;q<n;++q) { bi+=op.Btil[l*n+q]*u[idx(dir,q,s,r)];
                        di+=op.Dtil[l*n+q]*u[idx(dir,q,s,r)]; }
```
Notice: `Btil` (interpolation) and `Dtil` (derivative) are **small 1D operators** (`n × n`,
where `n = p+1`). The kernel sweeps them along one direction at a time (`dir`), not as one big
dense matrix. That 1D-at-a-time structure is *why* the cost doesn't blow up with `p`.

**Code jump — the flat-throughput result behind `fig_throughput`:**
`docs/figures/make_matfree_figs.py`:
```python
mdofs = [4794, 7984, 7895, 5309, 4258, 4157, 4530]   # MDOF/s, p=1..7
ax.text(..., "flat $\\approx$4--8 GDOF/s  ($p{=}7\\approx p{=}1$)", ...)
```
Notice: throughput at `p=7` (4530 MDOF/s) is essentially the same as at `p=1` (4794). You get
the 570×-better accuracy of high order *for free* in throughput terms. (One honest caveat:
the kernel is memory- and occupancy-bound, so FP64 tensor cores do **not** help here — that's
a measured negative result, not a missed opportunity. Chapter 2 owns this.)

<img src="figures/fig_throughput.png" width="460" alt="Apply throughput vs p on H100">

*Apply throughput per order on a single H100. Every bar sits in the ~4–8 GDOF/s band;
`p=7` (4.5 GDOF/s) matches `p=1` (4.8 GDOF/s), so high order costs nothing in throughput.*

### Where the geometry actually lives

Matrix-free still needs *some* stored data per element: a small per-point metric that encodes
how the reference cube maps to the real, possibly-curved element. This is the one array that
sets the memory ceiling.

**Code jump — the per-element metric:**
`backend/distributed/unstructured/fem/mars_cvfem_ho_apply.hpp`:
```cpp
computeElementMetric(const HoCvfemOperators& op, const double corners[8][3])
```
Notice: the metric is computed once from the 8 corner coordinates and reused for every `A·x`.
For a *general curved mesh* it costs ~7.2 KB/element of HBM — and that single array, not the
arithmetic, is what caps a GPU at **647M DOF/GPU** at `p=4`.

### The journey ahead

This tutorial walks the same path the MARS operator took: a bit-exact single-GPU operator
(Ch. 2), made distributed without losing a digit (Ch. 3–4), scaled until communication stops
mattering (Ch. 5), and aimed at a trillion DOF on a fraction of the machine (Ch. 6). For
scale context: **MFEM won Gordon Bell 2025** with 55.5 trillion DOF on 43,520 MI300A GPUs
(~1.27B DOF/GPU); on Alps that class is ~9.3T DOF on ~9,200 GH200s. The MARS HO operator sits
in the same **~1B DOF/GPU** league. The trillion frontier itself is **in progress, not done**
— we will be precise about that in Chapter 6.

---

## Chapter 2 — The Matrix-Free Operator (Single GPU)

In Chapter 1 we framed the problem. Now we *do* the operation. Almost every iterative solver
— Conjugate Gradient, GMRES, multigrid — needs exactly one thing from the discretized PDE,
over and over again:

> Given a vector `x`, compute `y = A x`.

`A` is the stiffness matrix. For a diffusion (Laplacian) problem, `A x` answers "if the
solution were `x`, what would the residual fluxes be?" The whole game of a fast solver is
making this one operation — the **matvec** — as cheap as possible.

This chapter builds intuition for why a matrix-free matvec is not just possible but *faster*
at high order, walks through the actual GPU kernel, and ends with how we *prove* it is correct
— plus one honest negative result about tensor cores.

### 2.1 What is a matvec, really?

If you have `A` stored as a sparse matrix (CSR: values + column indices), the matvec is a
loop: for each row `i`, `y[i] = sum_j A[i][j] * x[j]`. You read every nonzero of `A` from
memory once per matvec. That is the "assembled" approach.

But where did `A` come from? It came from summing little contributions, one per mesh element.
Element `e` couples its own DOF to each other; the global `A` is just all those element
matrices `K_e` scattered into the right rows and columns and added up.

The matrix-free insight: instead of *storing* the assembled `K_e` and reading them back, we
**recompute the action of `K_e` on the fly**, element by element. The pattern for every
element is:

1. **Gather** — pull this element's slice of `x` out of the global vector (using the
   element-to-DOF map).
2. **Apply** — compute `y_e = K_e x_e` for this one element, *without ever forming `K_e`*.
3. **Scatter** — add `y_e` back into the global `y` at this element's DOF rows.

Do that for all elements and you have computed `y = A x`. Shared DOF on element boundaries get
contributions from several elements — the scatter-add handles that automatically.

Why bother? Two reasons. First, memory: storing `A` at high polynomial order is enormous (we
saw the crossover in Chapter 1). Second, and counterintuitively, the on-the-fly apply can be
made *cheaper* than reading a stored matrix — if we are clever about step 2. That cleverness
is sum-factorization.

### 2.2 The element: a tensor-product box of GLL nodes

We use **high-order** hex elements. Order `p` means each element carries `(p+1)` solution
nodes per direction, arranged on a 3D tensor-product grid: `(p+1)^3` DOF per element. The
nodes are not evenly spaced — they sit at **Gauss-Lobatto-Legendre (GLL)** points, which keeps
high-order interpolation stable.

The key structural fact: because the element is a tensor product, *every operation factors
into 1D operations along each axis*. We never need a dense `(p+1)^3 × (p+1)^3` element matrix.
We only need small 1D operators, built once on the host.

**Code jump — the 1D reference operators.**
`backend/distributed/unstructured/fem/mars_cvfem_ho_basis.hpp`, `struct HoCvfemOperators`:
```cpp
struct HoCvfemOperators {
    int p = 0;
    std::vector<double> zeta;       // (p+1) GLL nodes
    std::vector<double> xi;         // p Gauss points
    std::vector<double> Btil;       // p x (p+1)
    std::vector<double> Dtil;       // p x (p+1)
    std::vector<double> D;          // (p+1)^2 GLL derivative
    std::vector<double> W;          // (p+1)^2 integration
    std::vector<double> Deltatil;   // (p+1) x p
};
```
What to notice:
- These are *tiny* — at most `(p+1)^2`-sized. The whole set is a few KB even at `p=8`. That is
  the entire "operator." There is no element matrix.
- `zeta` are the `(p+1)` GLL solution nodes (the DOF). `xi` are the `p` Gauss quadrature
  points — the subcontrol-surface (SCS) faces where CVFEM evaluates fluxes.
- `Btil` interpolates from GLL nodes to the flux points; `Dtil` takes the derivative there;
  `D` differentiates GLL-to-GLL; `W` integrates a flux over a subcontrol interval.

This follows Knaus' high-order CVFEM formulation ([Knaus 2022](#references), SAND2022-3366J). At `p=1` it collapses to
ordinary linear CVFEM: `zeta={-1,1}`, `xi={0}`, `Btil=[1/2,1/2]`, `Dtil=[-1/2,1/2]`.

### 2.3 Sum-factorization: why O(p⁴) beats O(p⁶)

Here is the heart of the chapter. Naively, applying a dense element matrix is
`(p+1)^3 × (p+1)^3` work — that is **O(p⁶)** per element. For `p=4` that is `125 × 125 ≈
15,600` multiply-adds per element. It gets brutal fast.

Sum-factorization exploits the tensor structure. A 3D contraction over a tensor-product grid
can be done as **three sequential 1D sweeps**, one per axis. Each sweep is a small 1D matrix
`(p+1)×(p+1)` applied along one direction, across the `(p+1)^2` lines in the other two
directions. That is `(p+1)^2 × (p+1)^2` work per sweep — **O(p⁴)**, done a constant number of
times.

The payoff is dramatic and measurable: because the cost per DOF stops growing like `p²`,
**throughput stays roughly flat as the order rises** — in MARS, `p=7` runs at essentially the
same GDOF/s as `p=1` (~4–8 GDOF/s/GPU). You get much higher accuracy per DOF at no throughput
penalty. That is *the* reason to go high-order matrix-free.

**Code jump — the host reference apply (the ground truth).**
`backend/distributed/unstructured/fem/mars_cvfem_ho_apply.hpp`, `applyHoCvfemElement`. This is
Knaus Alg 2 ([Knaus 2022](#references)): loop over each direction `dir` and each SCS face `l`, and
do the three 1D contractions:
```cpp
for (int dir = 0; dir < 3; ++dir)
    for (int l=0;l<p;++l) {
        // step 1: normal interp(Btil) + deriv(Dtil) along `dir`
        for (int s=0;s<n;++s) for (int r=0;r<n;++r) {
            double bi=0, di=0;
            for (int q=0;q<n;++q) { bi+=op.Btil[l*n+q]*u[idx(dir,q,s,r)];
                                    di+=op.Dtil[l*n+q]*u[idx(dir,q,s,r)]; }
            interp[s*n+r]=bi; deriv[s*n+r]=di;
        }
        // step 2: tangential D-derivatives + flux = metric * grad
        // step 3+4: W-integrate the flux over the face (both tangential axes)
        for (int s=0;s<n;++s) for (int r=0;r<n;++r) {
            y[idx(dir,l,  s,r)] -= intf[s*n+r];   // distribute -/+ to the
            y[idx(dir,l+1,s,r)] += intf[s*n+r];   // two nodes bounding face l
        }
    }
```
What to notice:
- Every inner `q`-loop is length `n=(p+1)` — a 1D contraction. No `O(p⁶)` tensor ever appears.
- `idx(dir,nrm,t1,t2)` re-strides the *same* flat array so "the normal axis" is whichever of
  x/y/z `dir` points at. One array, read three ways.
- The last two lines are the CVFEM finite-volume step: a flux on subcontrol face `l` leaves
  node `l` and enters node `l+1`. That `-/+` pair is the discrete divergence — and it is
  exactly why `A·1 = 0` later (a constant has zero flux, so every face contributes `-c+c=0`).

This host version is deliberately simple and slow. It exists to be the bit-exact reference the
GPU kernel is checked against.

### 2.4 The per-element metric: where the geometry lives

The 1D operators above are defined on the *reference* cube `[-1,1]³`. Real elements are
stretched, sheared, possibly curved. The geometry enters through one object: the **metric**
`G`, computed per element from its 8 corner coordinates via the Jacobian. It is
`detJ · J⁻¹J⁻ᵀ` at each flux point — it tells the apply how a reference-space gradient maps to
a physical-space flux, including cross-axis coupling when the element is not axis-aligned.

**Code jump — the metric.**
`mars_cvfem_ho_apply.hpp`, `computeElementMetric`:
```cpp
// Gvec = metric * e_dir = detJ * (J^{-1}J^{-T})[:,dir]
double gvec[3];
for (int a=0;a<3;++a) { double v=0; for(int k=0;k<3;++k) v+=Ji[a][k]*Ji[dir][k]; gvec[a]=det*v; }
auto& g = G[((size_t)(dir*p + l)*n + s)*n + r];
g[2]=gvec[dir]; g[1]=gvec[t1axis[dir]]; g[0]=gvec[t2axis[dir]];
```
What to notice:
- `g[2]` is the **normal** coefficient, `g[0]`/`g[1]` are the two **tangential** (cross-term)
  coefficients. On an orthogonal cube the cross terms are zero and only `g[2]` survives — the
  "Affine" fast path.
- This metric is computed **once** and reused across every matvec. The apply kernel never
  touches geometry again — it just reads `G`. The apply is on the solver's hot path; the
  metric build is not.
- The flip side, and the reason for our per-GPU memory ceiling: `G` is `3·p·(p+1)²` vec3s per
  element of HBM (~7.2 KB/elem at `p=4`). That array is what caps us at ~647M DOF/GPU at
  `p=4`.

**This pattern has a name: partial assembly.** Store the geometric metric per quadrature point,
then apply the operator via sum-factorization at solve time — never assemble a global sparse
matrix and never recompute geometry on the hot path. It is exactly the method MFEM uses (store
only the per-point data `D`, evaluate the basis/gather operators on the fly under the tensor-
product structure — [MFEM partial assembly](#references)), and the method the 2025 Gordon Bell
winner used for its 3D wave forward operator. The lineage is direct: our `d_G` is byte-for-byte
Knaus's per-element metric `G` (`N_el × 3p(p+1)² × 3`), and Knaus's high-order CVFEM scheme
([Knaus 2022](#references), SAND2022-3366J) is itself partial
assembly — it stores `vol`, `m_dot`, `A`, `G` per element, `O(p³)` total, and calls it a
"memory-efficient scheme." So this is **the same method as MFEM and Knaus, only a different
discretization**: CVFEM subcontrol-volume fluxes instead of variational Galerkin. (Worth being
precise: "matrix-free" in all three means *no assembled global matrix* — both Knaus and MARS
keep the metric **persistently in DRAM**. Knaus's "registers" is only the in-kernel working
slice of one element for `p≤4`, not a separate storage tier. The truly recompute-everything
matrix-free variant stores nothing and pays the FLOPs back on every apply; in MFEM's own
benchmarks that variant *lost* on wall-clock to partial assembly. MARS's contribution over Knaus
is carrying this identical PA operator to *distributed*, all-order (`p≤7`), trillion-DOF scale.)

### 2.5 The GPU kernel: gather / apply / scatter, one thread per face slot

Now the device twin. The design constraint is *not* compute — sum-factorization made the
arithmetic cheap. The constraint is **shared-memory bandwidth and occupancy**. The kernel in
`mars_cvfem_ho_matfree.hpp` is built around three decisions:

1. The tiny reference operators (`Btil`, `Dtil`, `D`, `W`) live in **`__constant__` memory**,
   not shared. They are read-only, identical for every element, and the constant cache
   broadcasts one value to a whole warp in a single transaction — at zero shared-memory cost.
   This is the single biggest occupancy win.
2. **One thread per tangential `(s,r)` face slot.** That thread owns its entire normal column,
   so the normal-direction contraction and the `-/+` scatter to nodes `l, l+1` need *no*
   cross-thread communication and *no* `__syncthreads`. Only the two tangential contractions
   exchange data, through small shared face buffers.
3. **Pack `E` elements per block** at low order, so a `p=1` element (only 4 face slots) does
   not starve an SM.

**Code jump — the gather (and what `dof < 0` means).**
`mars_cvfem_ho_matfree.hpp`, `ho_cvfem_apply_kernel`:
```cpp
for (int l = laneInElem; l < N3; l += threadsPerElem) {
    my_u[l] = (valid && edof[l] >= 0) ? d_u[edof[l]] : RealType(0);
    my_y[l] = RealType(0);
}
__syncthreads();
```
What to notice:
- `edof` is this element's DOF map; the gather pulls `x` into shared memory `my_u`. A negative
  DOF is a constrained/absent node and contributes 0 — same convention as the host.
- `my_y` (the per-element accumulator) is zeroed here and *cannot* alias `my_u`: we keep
  reading `u` in directions 1 and 2 while `y` is already accumulating from direction 0.

**Code jump — the apply inner step (sum-factorization on-device).**
Same kernel, step 1 of the `(dir, l)` sweep:
```cpp
for (int slot = laneInElem; slot < NN; slot += threadsPerElem, ++it) {
    const int s = slot / N, r = slot % N;
    RealType bi = 0, di = 0;
    #pragma unroll
    for (int q = 0; q < N; ++q) {
        RealType uq = my_u[ho_idx<N, NN>(dir, q, s, r)];
        bi += c_Btil[l * N + q] * uq;   // interp, GLL -> flux point
        di += c_Dtil[l * N + q] * uq;   // normal derivative
    }
    faceA[s * NP + r] = bi;     // interp -> shared face buffer (tangential D needs it)
    deriv_cache[it]   = di;     // deriv stays in a register (only this thread needs it)
}
__syncthreads();               // interp now visible to the tangential contractions
```
What to notice:
- This is line-for-line the host `step 1`, but `c_Btil`/`c_Dtil` come from constant memory and
  the loop is over face *slots*, one per thread.
- `interp` must go to shared memory because the *tangential* D-contraction reads other
  threads' values. `deriv` never crosses threads, so it stays in a register — cheaper. This
  split is why the kernel needs only two small shared face buffers instead of a full `n³`
  tensor in smem.

**Code jump — the scatter (hazard-free, atomics only at the global edge).**
Step 4 of the same sweep, then the final scatter:
```cpp
my_y[ho_idx<N, NN>(dir, l,     s, r)] -= intf;
my_y[ho_idx<N, NN>(dir, l + 1, s, r)] += intf;
...
// at the end, additive scatter to global:
if (e < numElements) {
    for (int l = laneInElem; l < N3; l += threadsPerElem) {
        int dof = edof[l];
        if (dof >= 0) atomicAdd(&d_y[dof], my_y[l]);
    }
}
```
What to notice:
- *Within* the element, thread `(s,r)` owns the whole normal column, so it is the only writer
  of `my_y[idx(dir,l,s,r)]` — the `-/+` across the serial `l`-loop never races. No
  `__syncthreads` for the write itself.
- The only atomics are the *global* scatter-add at the very end, because shared DOF on element
  faces receive contributions from neighboring elements. That is the irreducible cost of
  assembling a global `y` from element pieces.

### 2.6 How we test it: A·1 = 0 and A·linear = 0

A matrix-free operator is dangerous: a bug produces *some* number, not a crash. So we lean on
two physics-grounded patch tests that any correct diffusion operator must pass. Both live in
`mars_cvfem_ho_matfree_test.cu`.

**Test B — A·1 = 0 (the constant nullspace).** Diffusion of a constant field produces no flux.
So if we feed `x = 1` everywhere, the output must be machine zero at *every* DOF — interior
and boundary alike. This catches sign errors, mis-wired scatter, wrong operator indexing —
anything that breaks the discrete divergence.
```cpp
std::vector<double> ones(nDof, 1.0);
... gpuApply<P>(d_u, d_y, d_elemDof, d_G, nEl, nDof) ...
double nullMax=0; for (long i=0;i<nDof;++i) nullMax = std::max(nullMax, std::abs(y1[i]));
...
bool okB = (nullMax < 1e-9);
```
This works structurally because of the `-/+` flux distribution we saw: a constant gives equal
flux into and out of every subcontrol face, so each face's contribution cancels exactly. In
the distributed setting this same check becomes the headline gate (Chapter 4).

**Test C — A·linear = 0 at interior DOF (consistency).** A linear field `x = x_coordinate` has
constant gradient, so its Laplacian is zero — meaning `A x` must vanish at every *interior*
DOF (boundaries carry the flux that balances the linear ramp, so we exclude them):
```cpp
double interMax=0;
for (long i=0;i<nDof;++i) if(!bdry[i]) interMax=std::max(interMax,std::abs(yx[i]));
...
bool okC = (interMax < 1e-9);
```
Test B checks the nullspace; Test C checks that the operator reproduces *linear* fields
exactly — together they pin down a consistent Laplacian. The test also runs **Gate A**: the
GPU metric kernel must reproduce the host `computeElementMetric` bit-for-bit
(`metricErr < 1e-12`), and a single-element GPU apply must match the host
`applyHoCvfemElement` to `< 1e-12` relative — the only allowed difference is floating-point
reduction order.

### 2.7 An honest negative result: FP64 tensor cores do not help

The natural next thought on a GH200 is: the inner contractions are small matrix multiplies —
feed them to the FP64 tensor cores (DMMA). We built and measured that path. **It does not
help.**

The reason is exactly the design constraint from §2.5: this operator is **shared-memory-
bandwidth and occupancy bound, not compute bound.** Sum-factorization already shrank the
arithmetic to `O(p⁴)`; the bottleneck is moving small face buffers through shared memory and
keeping enough blocks resident. Tensor cores accelerate the part that was already cheap, while
doing nothing for the part that is actually limiting — and they add register/shared pressure
that *hurts* occupancy. The sibling DMMA Galerkin kernel measured 90 KB smem/block, capping it
to 2 blocks/SM (~12.5% occupancy); that is the wall, and matrix-multiply throughput is not
what you hit it with. `ncu` confirms the diagnosis directly: ~15–19% compute utilization against
71–82% L1/shared throughput, and the DMMA path gave only ~1.06×. The apply itself already runs
at ~60% of peak HBM — it is bandwidth-efficient, just not compute-bound.

One honest caveat so this is not mistaken for a hardware verdict: it is a **formulation** choice,
not a limit of the tensor cores. MFEM *does* profit from FP64 DMMA — by **batching many elements
into one large dense GEMM**, which is compute-bound at `p≥4`. MARS runs **one warp per element**,
so each contraction is a tiny GEMM that the tensor cores cannot saturate. The lever that would
turn DMMA on is the batched-GEMM layout (plus the affine-metric path of Ch. 5/6), not the
hardware. We chose warp-per-element for its occupancy and locality; the consequence is that
tensor cores stay idle.

This is worth stating plainly because it is easy to oversell. The win in this operator is
**sum-factorization plus an occupancy-first kernel layout** (constant-memory operators, one
thread per face slot, minimal smem). Tensor cores are a measured no-op here. Knowing *which*
resource binds you — and resisting the urge to throw the shiny hardware feature at the wrong
one — is the engineering lesson of this chapter.

### Recap

- A matvec `y = A x` is gather → apply → scatter, per element, with **no stored matrix**.
- High-order elements are tensor products of GLL nodes, so the apply factors into small **1D**
  contractions: **sum-factorization** turns `O(p⁶)` into `O(p⁴)`, keeping throughput flat.
- Geometry enters once through the per-element **metric** `G`, reused across every matvec; it
  is also what sets the per-GPU memory ceiling.
- The GPU kernel is built for **occupancy**: constant-memory operators, one thread per face
  slot, sync-free normal column, atomics only at the final global scatter.
- Correctness is proven by **A·1 = 0** (nullspace) and **A·linear = 0** (consistency),
  bit-checked against a host reference.
- **FP64 tensor cores are a measured negative result** — the operator is bandwidth/occupancy
  bound, not compute bound.

---

## Chapter 3 — Distributed DOF Numbering

### 3.1 Why high order makes numbering the hard part

At first order (`p=1`) life is simple: every degree of freedom sits on a mesh **node**, and
the mesh already gives every node a global id. The numbering is done before you start.

High order breaks that comfort. For a tensor-product hex of order `p`, the `(p+1)^3` DOF live
on four *kinds* of geometric entities:

- **corners** — 8 of them, shared by up to 8 elements
- **edges** — `p-1` interior DOF each, shared by several elements
- **faces** — `(p-1)^2` interior DOF each, shared by exactly 2 elements
- **interior** — `(p-1)^3` DOF, owned by one element alone

**Code jump** — `mars_ho_dof_handler.hpp`, class header comment:
```
// DOF layout (global):
//   [0, nCorner)                              corner DOFs (= corner node ids)
//   [nCorner, +nEdge*(p-1))                   edge-interior DOFs
//   [+, +nFace*(p-1)^2)                       face-interior DOFs
//   [+, +nElem*(p-1)^3)                       element-interior DOFs
```
Notice the global DOF id space is a concatenation of four blocks. Corner DOF *reuse* the
existing node ids, so corner continuity is free. Everything after `nCorner` is new and must be
invented consistently.

The central requirement of a distributed solver is **continuity**: two elements that share an
edge or face must agree, DOF-for-DOF, on the global id of every shared DOF — even when those
elements live on *different ranks* and were numbered by *different processes that never talked
first*. If they disagree, the matrix-free apply scatters a contribution to the wrong slot, and
the operator is silently wrong.

So the problem splits into two questions:

1. **Identity** — how do two elements (possibly on two ranks) name the *same* shared DOF
   identically? → the canonical key.
2. **Ownership** — among all the ranks that hold a shared DOF, which *one* owns it? →
   min-rank-among-holders.

### 3.2 The canonical key: name a DOF by its corners, never by coordinates

The tempting shortcut is to identify a high-order DOF by its physical `(x,y,z)` location and
hash that. This **fails** at scale: the SFC (space-filling curve) quantizes coordinates, and
at `p=7` interior nodes get packed so densely that two genuinely-different DOF can quantize to
the same key. The header calls this "the SFC quantization wall at p=7".

The fix is topological, not geometric. A shared entity is *defined* by the global corner ids
it connects, and those corner ids are already globally consistent (they came from the P1
mesh). So:

- an **edge** is named by its 2 endpoint global corner ids, *sorted*
- a **face** is named by its 4 corner global ids, *sorted*
- the position *within* the entity is named by a canonical index

Sorting is what makes the key rank-independent: element A may traverse an edge low→high and
element B high→low, but `{min, max}` is the same array on both sides.

**Code jump** — `mars_ho_dof_handler.hpp`, `struct DofKey`:
```cpp
// kind: 0 corner, 1 edge, 2 face, 3 interior. g* = sorted global defining-corner
// ids (unused = -1); pos = index within the edge/face/interior block. Two ranks
// that share a DOF compute the SAME key, so the halo matches by it.
struct DofKey { int kind; long g0, g1, g2, g3; int pos; };
```
Notice the key carries no coordinates at all — just `kind`, up to four sorted corner gids, and
a slot `pos`. This is the whole identity contract: *equal key ⇒ same DOF*, computed
independently on any rank.

Here is the edge key being built. Watch the sort and the matching `pos` flip:

**Code jump** — `mars_ho_dof_handler.hpp`, `buildDistributed`, the `cnt == 2` (edge) branch:
```cpp
long gA = cornerGid[cA], gB = cornerGid[cB];
long klo = (gA<=gB)?gA:gB, khi = (gA<=gB)?gB:gA;
int  pos = (gA<=gB) ? (t-1) : (pm1-t);   // canonical low->high
...
key = DofKey{1, klo, khi, -1, -1, pos};
```
Two things to notice. First, `{klo, khi}` is the sorted endpoint pair — orientation-free.
Second, `pos` is *reversed* (`pm1-t`) exactly when this element walks the edge high→low, so the
t-th interior node still lands in the same canonical slot the neighbour assigns it. The key is
identical on both sides.

Faces are the same idea with 4 sorted ids and a 2D canonical frame (`hex_face_canonical_pos`)
so the `(p-1)^2` interior slots agree under any of the 8 ways two hexes can relabel a shared
quad. Interior DOF get `kind=3` keyed on the element itself — they are never shared, so their
"key" only needs to be locally unique.

### 3.3 Ownership, and the bug that orphans DOF

Identity tells us *which DOF are the same*. Ownership picks *the single rank responsible for
it* — the rank that holds the authoritative value, sums contributions into it during
reverse-add, and broadcasts it during forward. Exactly one owner, agreed by everyone.

For two of the four kinds this is already solved:

- **corners** inherit the P1 cstone node owner. cstone guarantees the owner contains the node,
  so it is consistent and correct.
- **interior** DOF are element-local; owner = the element's owner. Never shared.

Edges and faces are the hard case, and here is the teachable bug.

**The wrong heuristic: "lowest global corner owns it."** It sounds reasonable — pick the
smallest corner gid of the edge/face and let *that corner's owner* own the whole entity. It is
wrong, and wrong in a way that passes single-rank tests and only breaks distributed.

The failure: the owner of the lowest corner may **not contain the edge or face at all**. A
corner can be shared by 8 elements spread across many ranks; the rank that happens to own that
corner node might hold none of the elements touching this particular edge. You have just
assigned ownership to a rank that has no slot for the DOF — it is **orphaned**. No one's
reverse-add lands on it; the value is garbage.

**Code jump** — `mars_ho_dof_handler.hpp`, `buildDistributed` doc comment:
```
//   edge/face-> PROVISIONAL myRank. The lowest-global-corner heuristic is WRONG
//               here: that corner's owner may have the corner but not the edge/face,
//               orphaning the DOF. resolveHoDofOwnership() (a peer key exchange)
//               sets the real owner = min rank among the ranks that actually hold
//               it. dofShared[d]=1 flags edge/face DOF whose defining corners are
//               all P1-shared -> candidates for that resolution.
```
Notice the design decision: `buildDistributed` does **not** try to guess the final owner. It
sets `owner = myRank` provisionally and flags the DOF in `dofShared`, deferring the real choice
to a step that has actual cross-rank information.

The flag itself, in the edge branch:

**Code jump** — `mars_ho_dof_handler.hpp`, `buildDistributed`, edge branch (continued):
```cpp
owner  = myRank;
shared = (sharedCorner[cA] && sharedCorner[cB]) ? 1 : 0;
```
A DOF is a *candidate* for resolution only when **all** its defining corners are P1-shared
(`sharedCorner`). If any corner is interior to this rank, no other rank can hold the entity, so
the DOF is purely local and `shared` stays 0 — we never pay communication for it.

### 3.4 The fix: min-rank-among-holders by peer key exchange

The correct owner must be a rank that **actually holds the DOF**. The set of holders of an
edge/face is exactly the set of ranks whose local elements touch it — and those ranks are
mutual mesh neighbours, so they are already P1 halo peers. So the rule is:

> **owner = the minimum rank among all ranks that hold this DOF.**

`min` is a deterministic tie-break that every holder computes to the same answer, and because
it is chosen *from the holders*, it can never orphan. The implementation is a small all-to-peers
exchange of the canonical keys.

**Code jump** — `mars_ho_halo.hpp`, `resolveHoDofOwnership` doc comment:
```
// Resolve ownership of shared HO edge/face DOF by MIN-RANK-AMONG-HOLDERS.
// ... owner = min(current, peer rank). The owner is thus
// always a rank that CONTAINS the DOF -> no orphans, and all holders agree (the
// holders of an edge/face are mutual neighbours).
```

Each rank packs its shared keys, sends them to every peer, and for every peer key that matches
one of its own shared DOF, lowers that DOF's owner toward the peer's rank if the peer's rank is
smaller:

**Code jump** — `mars_ho_halo.hpp`, `resolveHoDofOwnership`, the resolution loop:
```cpp
for (int i=0;i<np;++i) {
    int pr = peers[i];
    for (auto& k : peerKeys[i]) {
        auto it = myShared.find(k);
        if (it != myShared.end() && pr < dofOwner[it->second]) dofOwner[it->second] = pr;
    }
}
```
Notice the logic is symmetric: a key only updates ownership when *both* this rank and peer `pr`
have it (`it != end`), i.e. both are genuine holders. The `pr < dofOwner` comparison drives
every holder of a given key to the same minimum, with no central coordinator and no global
numbering pass. Provisional `myRank` from `buildDistributed` is simply the starting value of
that running minimum.

Because the matching is keyed on `DofKey` — the corner-sorted identity from §3.2 — a peer's key
finds the local DOF if and only if it is genuinely the same DOF. Identity and ownership reuse
the exact same key. That is the whole point of making the key canonical.

### 3.5 `dofShared` vs `dofBoundary`: two flags, two jobs

The numbering produces two boolean masks that look similar but serve different purposes, and
conflating them costs either correctness or memory.

**Code jump** — `mars_ho_dof_handler.hpp`, member declarations:
```cpp
std::vector<uint8_t> dofShared;   // 1 if edge/face DOF on a rank boundary
                                  // -> candidate for ownership RESOLUTION
std::vector<uint8_t> dofBoundary; // 1 if ANY shared DOF (corner/edge/face on
                                  // a rank boundary). Superset of dofShared
                                  // ... the halo keys/maps over dofBoundary,
                                  // NOT all numDof (else a numDof-sized
                                  // std::map OOMs the host at scale).
```
- **`dofShared`** is the *narrow* set: only edge/face DOF whose ownership is still undecided.
  `resolveHoDofOwnership` iterates exactly this set. Corners are excluded because their owner is
  already correct from cstone.
- **`dofBoundary`** is the *wide* set: every DOF that ever crosses a rank boundary, including
  shared corners. This is the superset, and it is what the **halo** keys and maps over.

**Code jump** — `mars_ho_dof_handler.hpp`, `buildDistributed`, write-back:
```cpp
dofShared[dof]   = (uint8_t)shared;            // edge/face -> ownership resolution
dofBoundary[dof] = (uint8_t)(boundary | shared); // any shared DOF -> halo keys/maps
```
Notice `boundary` carries the shared-corner case (a corner can be on a rank boundary yet need
no ownership resolution), while `shared` carries the edge/face case. The OR gives the halo its
full exchange set.

Why two sets and not one? Scale. The halo builds a `std::map` from key → local DOF. If it keyed
*all* `numDof`, that map is roughly 78 GB at 650M DOF/GPU and OOMs the host instantly. Keying
only `dofBoundary` makes it **O(surface area)** instead of O(volume) — and the surface-to-volume
ratio shrinks as the per-GPU problem grows, which is exactly why communication self-resolves
from ~22% down to ~4% as DOF/GPU climbs (Chapter 5). The memory discipline here is not a
micro-optimization; it is what lets the run reach billions of DOF per GPU at all.

### 3.6 Putting it together

The numbering pipeline, in order:

1. `build()` — number the local dense `elemDof` (corners reuse node ids; edges/faces deduped by
   sorted-corner key in `std::map`; interior is element-local).
2. `buildDistributed()` — tag every DOF with a canonical `DofKey`, a provisional owner (correct
   for corners and interiors, `myRank` for edges/faces), and the `dofShared` / `dofBoundary`
   masks.
3. `resolveHoDofOwnership()` — one peer key exchange turns the provisional edge/face owners into
   the true **min-rank-among-holders**, with no orphans and unanimous agreement.

After this, every shared DOF has exactly one owner, every rank names it identically, and the
halo (Chapter 4) can build a receiver-driven, truncation-free exchange straight off the same
`DofKey`.

The single most transferable idea: **identity before ownership, and both from one canonical
key.** Name the shared thing the same way on every rank first; only then decide who owns it —
and decide it from the set that actually holds it, never from a proxy like "lowest corner."

---

## Chapter 4 — The Distributed Matvec and the Halo

In the single-GPU chapters we treated `y = A·x` as a closed loop: every DOF the kernel touches
lives in one device array, the kernel reads it, the kernel writes it. The moment we split the
mesh across GPUs that assumption breaks. This chapter is about the *one* thing that breaks and
the *one* communication pattern that fixes it — and how we prove the fix is correct with a
single number.

### 4.1 Why a distributed matvec needs a halo at all

Picture two ranks that share a mesh face. The DOF sitting *on* that shared face are physically
one DOF, but each rank stores its own copy in its own local arrays. When rank A loops over its
owned elements and applies the element operator, the elements touching the seam need the value
of `x` at those shared DOF — including the contribution that *rank B's* elements will make to
them.

So a distributed matvec is not one operation, it is a sandwich:

1. **forward-halo(x)** — pull neighbour-owned DOF values into my ghost slots, so my `x` is
   complete over every element I own.
2. **apply OWNED elements** — run the exact same single-GPU kernel, but only over the elements
   this rank owns.
3. **reverse-add(y)** — my elements deposited partial sums into ghost slots that I don't own;
   send those back and *add* them onto the true owner. After this, every owned `y` slot holds
   the full sum of all element contributions, from every rank.
4. **dot products → MPI_Allreduce** — Krylov inner products are global sums over *owned* DOF
   only (ghosts would double-count).

The key mental model: a DOF on a partition boundary is *owned* by exactly one rank and
*ghosted* by the others. `forward` is "owner broadcasts the input"; `reverseAdd` is "ghosts
return their output contribution." Interior DOF (strictly inside one element) are never shared,
so they never enter the halo at all — only corners, edges, and faces on the seam do.

**Code jump — `examples/.../mars_ho_dist_apply_test.cu`, `runDistApply`, the A·1 matvec:**
```cpp
std::vector<double> u(nDof, 0.0);
for (long d = 0; d < nDof; ++d) if (dof.dofOwner[d] == rank) u[d] = 1.0;
halo.forward(u);                                  // ghosts <- owner's 1
cudaMemcpy(d_u, u.data(), ...);
cudaMemset(d_y, 0, ...);
ho_cvfem_apply_launch<double, P>(d_u, d_y, d_elemDof, d_G, nEl);   // OWNED elems
...
halo.reverseAdd(y);                               // ghost y summed into owners
```
Notice `u` starts as 1 *only on owned DOF*; `forward` is what makes the ghost slots also 1. The
apply launch is byte-for-byte the single-GPU kernel — distribution lives entirely in
`forward`/`reverseAdd`, not in the operator.

### 4.2 Receiver-driven construction: making `A.send[B] == B.recv[A]` true by design

The dangerous part of any halo is not the exchange — it is *building the send/recv lists so
they agree*. If rank A thinks it should send 5 DOF to B but B only expects 3, MPI either
truncates (silent data loss) or aborts (`MPI_ERR_TRUNCATE`). Symmetry of the two lists is the
whole ballgame.

The trap is to build the lists from *ownership*: "I own this shared DOF, so I'll send it to
everyone who might want it." That is **sender-driven**, and it over-claims — the sender guesses
the receiver's interest and the guess is wrong at corners where three or four partitions meet.

We do the opposite. **Receiver-driven**: the rank that *needs* a ghost value is the one that
initiates. Each rank looks at its ghost DOF (`dofOwner != me`), packs the canonical `DofKey` of
each, and *requests* it from the owner. The owner doesn't decide what to send — it answers
requests. Because every send entry exists *only* because a matching request arrived,
`A.send[B] == B.recv[A]` is true by construction, not by hope.

**Code jump — `mars_ho_halo.hpp`, `HoHalo::build`, the ghost-request pass:**
```cpp
for (int d = 0; d < numDof; ++d) {
    if (!dofBoundary[d]) continue;
    int o = dofOwner[d];
    if (o == myRank || o < 0) continue;          // owned or interior
    auto it = peerIdx.find(o);
    if (it == peerIdx.end()) continue;           // owner not a neighbour
    recvLocal[it->second].push_back(d);          // I will RECEIVE this DOF
    reqKeys[it->second].push_back(packKey(dofKey[d]));  // ...and I ask for it by key
}
```
What to notice: the recv list is built *first*, from my own ghosts. The send list does not
exist yet — it will be whatever requests arrive.

The requests then go out, and the owner turns each received key back into one of its local DOF
— *in the order the keys arrived*. That received order is exactly the requester's recv-slot
order, so forward and reverse stay aligned slot-for-slot without any extra sorting.

**Code jump — `mars_ho_halo.hpp`, `HoHalo::build`, building the send list from requests:**
```cpp
for (int i = 0; i < np; ++i) {
    sendLocal[i].reserve(gotKeys[i].size());
    for (auto& k : gotKeys[i]) {
        auto it = keyToLocal.find(k);
        if (it != keyToLocal.end()) sendLocal[i].push_back(it->second);
    }
}
```
The send list is literally "the keys someone asked me for, mapped to my local DOF." There is no
independent ownership-based guess to disagree with the receiver. That is why there is no
truncation and no over-claim.

One scaling subtlety worth pausing on: we key *only boundary DOF* into the lookup map. Keying
all `numDof` would build a `std::map` with one entry per DOF — about 78 GB at 650M DOF/GPU —
and OOM the host. Since only surface DOF are ever exchanged, the `dofBoundary[d]` gate keeps
construction O(surface), not O(volume):

**Code jump — `mars_ho_halo.hpp`, `HoHalo::build`, the boundary-only key map:**
```cpp
std::map<std::array<long,6>, int> keyToLocal;
for (int d = 0; d < numDof; ++d) if (dofBoundary[d]) keyToLocal.emplace(packKey(dofKey[d]), d);
```

### 4.3 Forward and reverse are mirror images

Once the CSR send/recv lists exist, the two exchange directions are nearly identical code —
they just swap the roles of the send and recv offsets, and swap *overwrite* for *add*.

**Code jump — `mars_ho_halo.hpp`, `HoHalo::forward` and `reverseAdd`:**
```cpp
void forward(std::vector<RealType>& vec) const {        // owner -> ghost: OVERWRITE
    ...
    for (int i = 0; i < ns; ++i) sbuf[i] = vec[sendDof_[i]];
    exchangeVals(sbuf, rbuf, sendOffsets_, recvOffsets_, 0x484b);
    for (int i = 0; i < nr; ++i) vec[recvDof_[i]] = rbuf[i];   // =
}
void reverseAdd(std::vector<RealType>& vec) const {     // ghost -> owner: ADD
    ...
    for (int i = 0; i < ns; ++i) sbuf[i] = vec[recvDof_[i]];
    exchangeVals(sbuf, rbuf, recvOffsets_, sendOffsets_, 0x484c);  // offsets swapped
    for (int i = 0; i < nr; ++i) vec[sendDof_[i]] += rbuf[i];  // +=
}
```
Two things to notice. First, the offset arguments to `exchangeVals` are swapped between the two
— reverse flows along the same edges, opposite direction. Second, the *only* operator
difference is `=` versus `+=`. Forward overwrites because each ghost has exactly one owner (one
value is authoritative). Reverse accumulates because one owned DOF may collect partial sums from
several ghosting peers.

The exchange itself is the standard non-blocking pattern — all `Irecv` posted, then all
`Isend`, then one `Waitall`. Because the peer lists are symmetric, there is no deadlock and no
unexpected message.

### 4.4 The device path: full-GPU per matvec, GPUDirect exchange

The host path above is the *correctness* path — simple, easy to reason about, used by the
gates. But on the scaling runs we cannot afford to copy `x` and `y` to the host every
iteration. The device path keeps everything in HBM: a CUDA kernel gathers the boundary values
straight out of the device solution vector into a contiguous send buffer, MPI sends *device
pointers* (GPUDirect / CUDA-aware MPI), and a second kernel scatters the received buffer back.

**Code jump — `mars_ho_halo.hpp`, gather/scatter kernels:**
```cpp
template<typename RealType>
__global__ void hoHaloGatherKernel(const RealType* vec, const int* idx, RealType* buf, int n)
{ int i = blockIdx.x*blockDim.x + threadIdx.x; if (i < n) buf[i] = vec[idx[i]]; }

template<typename RealType>
__global__ void hoHaloScatterAddKernel(RealType* vec, const int* idx, const RealType* buf, int n)
{ int i = blockIdx.x*blockDim.x + threadIdx.x; if (i < n) atomicAdd(&vec[idx[i]], buf[i]); }
```
The reverse direction uses `atomicAdd`, not a plain store — same reason as the host `+=`: a
single owned DOF can receive contributions from multiple peers landing in overlapping slots, so
the accumulation must be atomic.

**Code jump — `mars_ho_halo.hpp`, `HoHalo::forwardDevice`:**
```cpp
if (nSend_) { ... hoHaloGatherKernel<<<g,b,0,stream>>>(d_vec, d_sendDof_, d_sendBuf_, nSend_); }
cudaStreamSynchronize(stream);
... MPI_Irecv(d_recvBuf_+recvOffsets_[p], c, mpiT, peers_[p], 0x4860, ...);   // device ptr
... MPI_Isend(d_sendBuf_+sendOffsets_[p], c, mpiT, peers_[p], 0x4860, ...);   // device ptr
if(!rq.empty()) MPI_Waitall(...);
if (nRecv_) { ... hoHaloScatterKernel<<<g,b,0,stream>>>(d_vec, d_recvDof_, d_recvBuf_, nRecv_); }
```
The MPI buffers are `d_sendBuf_` / `d_recvBuf_` — raw device pointers handed directly to MPI.
The only host involvement is launching three kernels and posting the requests; the solution
vector itself never leaves HBM. The single `cudaStreamSynchronize` before MPI is the necessary
fence — the send buffer must be fully gathered before it goes on the wire.

The device path and host path are *the same exchange* expressed twice; the driver runs both and
checks they produce bit-identical results, so you can debug on the simple one and ship the fast
one.

### 4.5 The correctness oracle: A·1 = 0 across ranks

How do you know a cross-rank assembly is right, without a reference solution to compare against?
Use a property the operator must satisfy *by construction*. The unconstrained high-order
Laplacian sums every row to zero — physically, a constant field has zero diffusion. So if we
set `x = 1` everywhere and apply, the result must be `0` at every DOF. A single wrong send/recv
pairing leaves an O(1) residual right at the partition interface, because the seam DOF didn't
collect their neighbour's contribution. **A·1 = 0 is a sharp, automatic detector of exactly the
bug distribution introduces.**

**Code jump — `examples/.../mars_ho_dist_apply_test.cu`, the gate reduction:**
```cpp
double locMax = 0; long locOwned = 0;
for (long d = 0; d < nDof; ++d)
    if (dof.dofOwner[d] == rank) { locMax = std::max(locMax, std::abs(y[d])); ++locOwned; }
double gMax = 0; MPI_Allreduce(&locMax, &gMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
...
printf("max|A.1| over owned ... = %.3e   [%s]", gMax, gMax < 1e-8 ? "PASS" : "FAIL");
```
The max is taken over *owned* DOF only — a ghost slot is allowed to hold a partial sum; only the
owner holds the final, must-be-zero value. The `MPI_Allreduce(MPI_MAX)` turns it into one global
number.

In practice this gate passes at **1.7e-18 / 2.5e-17 / 2.8e-17 for p = 2 / 3 / 4 on 4 ranks** —
machine zero, i.e. the cross-rank assembly reproduces the single-rank result to the last bit.
The same A·1 check is run a second time through the *device* halo (`forwardDevice` /
`reverseAddDevice`), and the device result is bit-identical to the host result, which is how we
trust the GPUDirect path on the big runs.

One last connection back to Chapter 3, because it is where the subtle bug hid. The send/recv
lists are only correct if `dofOwner` is correct, and for shared edge/face DOF "owner = the
lowest global corner's owner" is **wrong** — it orphans DOF whose corner-owner doesn't actually
contain the edge. The fix is *min-rank-among-holders*: every holder broadcasts its shared-DOF
keys, and ownership goes to the smallest rank that genuinely contains the DOF.

**The takeaway for this chapter:** a distributed matrix-free matvec is the single-GPU kernel
wrapped in `forward → apply-owned → reverseAdd`, with global dot-products reduced over owned DOF.
The halo is correct because it is *receiver-driven* (`A.send[B] == B.recv[A]` by construction,
so MPI can never truncate or over-claim), it is *cheap* because only boundary DOF are keyed and
exchanged, it is *fast* because the device path gathers/scatters in HBM and sends device pointers
over GPUDirect, and it is *trusted* because `A·1 = 0` across ranks is a sharp, reference-free
oracle.

---

## Chapter 5 — Scaling to Billions

We have a high-order matrix-free operator that is bit-exact and runs at a flat ~4–8 GDOF/s per
GPU regardless of polynomial order. That is the single-GPU story. This chapter is about what
happens when you put 8, 64, or a few thousand of those GPUs together and push the global problem
toward a trillion degrees of freedom.

The headline to carry through: **at scale, the things that hurt small problems stop mattering,
and the thing that limits you becomes memory, not communication.** We earn each claim with a
measurement and the code that produced it.

### 5.1 The shape of the problem: weak scaling

There are two ways to scale. *Strong* scaling holds the global problem fixed and adds GPUs —
each GPU does less work, and eventually communication dominates. *Weak* scaling holds the
**work per GPU** fixed and adds GPUs to grow the global problem. For a trillion-DOF run you
cannot fit the problem on one GPU, so weak scaling is the only honest frame.

Two efficiencies matter, and they are very different curves:

- **Apply efficiency** — just the operator (`ho_cvfem_apply_launch`), no halo. Pure local
  compute. Perfect weak-scaling keeps this at 100%.
- **Full-matvec efficiency** — forward-halo → apply → reverse-add. Includes every MPI message.
  This is what an actual Krylov iteration pays.

**Code jump** — `docs/figures/make_matfree_figs.py`, figure 7 (`fig_ho_scaling`):
```python
ho_gpus      = [1, 8]
ho_apply_eff = [100.0, 98.1]   # apply-only, ~4.2M DOF/GPU held (5965 -> 46811 MDOF/s)
ho_full_eff  = [100.0, 78.4]   # full matvec, blocking halo (5790 -> 36319 MDOF/s)
```
Notice: apply weak-scales at **98%** — the operator itself barely notices it is distributed. The
full matvec sits at **78%** at the *small* per-GPU size (~4.2M DOF). That ~20-point gap is the
halo, and it **self-resolves as per-GPU size grows** (4M→11M DOF/GPU pushes full-matvec
78%→87%). We see why next.

<img src="figures/fig_ho_scaling.png" width="460" alt="HO matrix-free distributed weak scaling p=4">

*p=4, ~4.2M DOF/GPU held fixed (device==host halo bit-exact at 1e-18). Operator apply holds
98% from 1→8 GPU; the full matvec sits at 78% (blocking halo). Growing per-GPU size to
4→11M DOF closes that gap: comm 22→15%, full matvec 78→87%.*

At the largest run we did — **40 billion DOF on 64 GPUs at p=4** — apply held ~6 GDOF/s/GPU, the
*same* rate as a single GPU (≈100% apply weak-scaling), and the full matvec was ~90% efficient.
So adding 64× the GPUs cost ~10% on the operation that matters.

The same picture holds on the larger 1/8/32-GPU GH200 ladder (2M elements/GPU held, 2M → 64M):

<img src="figures/fig_weakscale.png" width="460" alt="Weak scaling on GH200 across 1/8/32 GPUs">

*Apply weak-scales at 95% to 64M DOF / 32 GPU. The full matvec degrades faster (50% blocking,
56% with comm overlap) — overlap buys ~1.22× at this small per-GPU size, where comm is still
a large fraction.*

### 5.2 Why communication self-resolves

The intuition is geometry. A rank owns a 3D block of the mesh. The **compute** it does scales
with the block's *volume*, V. The **communication** it does — the halo — scales with the block's
*surface area*, which grows like V^(2/3). So the ratio of comm to compute scales like:

```
comm / compute  ~  V^(2/3) / V  =  V^(-1/3)
```

As you make each GPU's block bigger (more DOF/GPU), this ratio *shrinks*. Communication does not
vanish — it becomes a smaller and smaller fraction of a growing compute bill. And there is a
second gift from 3D: the number of neighbours a block can touch **saturates at 26** (6 faces + 12
edges + 8 corners). Past that, growing the problem adds volume but no new peers.

**Code jump** — `docs/figures/make_matfree_figs.py`, figure 5 (`fig_comm_pergpu`):
```python
pergpu_dof = [2.05e6, 8.06e6, 16.10e6, 32.16e6, 64.39e6]   # measured per-GPU DOF
comm_meas  = [51.9, 41.9, 33.5, 28.6, 19.2]                # measured comm %
kfit, pfit = 200.0, -0.36   # fit comm/apply = k * V^p ; p ~ -1/3
```
Notice: the fitted exponent is **−0.36**, essentially the V^(−1/3) the geometry predicts. Comm
goes 52% → 19% purely by making each GPU's job bigger. The fit extrapolates to the
trillion-relevant point: at ~300M DOF/GPU, comm is ~17%.

<img src="figures/fig_comm_pergpu.png" width="460" alt="Comm fraction vs DOF per GPU with V^-1/3 model">

*Measured comm fraction on 8 GPUs (52% → 19% as per-GPU DOF grows 2M → 64M) against the
V^(−1/3) surface-to-volume model. The dashed line marks the trillion-relevant ~300M DOF/GPU,
where comm drops to ~17%.*

**Code jump** — `examples/distributed/unstructured/mars_ho_dist_apply_test.cu`, `runDistApply`:
```cpp
double mdofsFull  = (double)gOwned / maxFull  / 1e6;
double mdofsApply = (double)gOwned / maxApply / 1e6;
double commFrac   = (maxFull - maxApply) / maxFull * 100.0;
```
Notice: comm fraction is *measured*, not modelled — the wall-clock gap between the full matvec
(`forwardDevice → apply → reverseAddDevice`) and apply-only, both timed over 50 iterations after
5 warmups. Across the weak-scaling ladder this number falls 22% → 15% → 4% as per-GPU DOF grows
4M → 11M → 625M. The lesson for a trillion-DOF run: **do not fight the halo, feed the GPU.**
Bigger per-GPU blocks make the comm problem disappear on their own.

This is also why the halo is, deliberately, *blocking* in this code. Overlapping comm with
compute buys you the most when comm is a large fraction — and that is precisely the regime
(small per-GPU) we are scaling *away* from. At the sizes a trillion run actually uses, there is
little left to overlap.

### 5.3 The real wall: per-GPU memory

If communication self-resolves, what actually stops you? **HBM.** And specifically, one array.

The matrix-free operator does not store a matrix, but it does store a per-quadrature-point
geometric metric tensor, `d_G`. For high order this dominates the footprint.

**Code jump** — `examples/distributed/unstructured/mars_ho_dist_apply_test.cu`, `runDistApply`:
```cpp
const size_t gLen = nEl * (size_t)(3 * P * n * n) * 3;   // n = P+1
...
if (me==cudaSuccess) me = cudaMalloc(&d_G, sizeof(double) * gLen);
```
Notice: `d_G` is `nElem × (3·P·(P+1)²) × 3` doubles. The `3·P·(P+1)²` is the number of
sub-control-volume faces per element in Knaus's scheme; the trailing `×3` is the metric vector
per face. At p=4 this is ~7.2 KB of HBM **per element** — far more than the DOF vectors
`d_u`/`d_y`. This array, not the messages, sets the ceiling.

How do we know the ceiling is ~647M DOF/GPU at p=4? We did not estimate it — we **probed it**.
The allocation is checked, and an out-of-memory is reported as a measurement, not a crash:

**Code jump** — `examples/distributed/unstructured/mars_ho_dist_apply_test.cu`, `runDistApply`:
```cpp
if (me != cudaSuccess) {
    double needGB = (sizeof(double)*gLen + sizeof(double)*2.0*nDof
                     + sizeof(int)*(double)nEl*N3 + sizeof(double)*(double)nEl*24) / 1e9;
    if (rank == 0)
        printf("p=%d  nDof=%ld nEl=%zu  DEVICE OOM (needs ~%.1f GB, d_G=%.1f GB): %s"
               "  -- per-GPU ceiling exceeded\n",
               P, nDof, nEl, needGB, sizeof(double)*gLen/1e9, cudaGetErrorString(me));
    cudaGetLastError();   // clear the sticky error so the next p can still run
    ...
    return;   // clean exit -- not a null-deref segfault
}
```
Notice three deliberate choices. First, the OOM message **breaks out `d_G`'s share** of the
total — that is how we attribute the ceiling to the metric. Second, `cudaGetLastError()` clears
the sticky error so a sweep over p can keep going. Third, it `return`s cleanly instead of
dereferencing a null pointer — so an OOM probe at scale gives you a number, not a corrupted run.
Sweeping `nEl` upward until this fires is exactly how 647M DOF/GPU was measured.

This ceiling is the good news for the trillion target. At 647M DOF/GPU, a 10¹² problem needs only
~1,550 GPUs (~390 nodes), versus ~9,300 GPUs for a p=1 *element*-based trillion:

**Code jump** — `docs/figures/make_matfree_figs.py`, figure 6 (`fig_trillion`):
```python
ceil_M = [70, 77, 108, 647]   # per-GPU ceilings: assembly CSR(p1), MF apply(p1), domain, HO MF apply(p4)
need_M = 1e12 / 10752 / 1e6    # ~93M DOF/GPU to fit 1e12 on full Alps
```
Notice the ceiling climbs from 70M (stored CSR) to 647M (high-order matrix-free) — a ~9× density
gain, and the direct reason high order is the right tool for the trillion frontier. We are well
above the 93M/GPU that a full-Alps trillion would require.

<img src="figures/fig_trillion.png" width="460" alt="Per-GPU capacity by mode, path to a trillion">

*Validated per-GPU ceilings: assembled CSR 70M, p=1 matrix-free 77M, domain decomposition
108M, and p=4 HO matrix-free 647M DOF/GPU. The dashed line is the ~93M/GPU a 10¹² problem
needs on full Alps — HO matrix-free clears it ~7×, putting a trillion DOF in ~1,550 GPUs.*

(For context, MFEM won Gordon Bell 2025 at 55.5 trillion DOF on 43,520 GPUs — ~1.27B DOF/GPU.
Our operator sits in the ~1B DOF/GPU class. The trillion run itself is **in progress, not done**;
Chapter 6 is honest about why.)

The lever to raise the ceiling further is in *what* `d_G` stores. Today we store the **general
per-point metric** — a full `detJ·J⁻¹J⁻ᵀ` at every flux point — because the Jacobian varies
inside the element. If the element geometry is **affine** (a parallelepiped: cube, sheared box,
uniformly stretched), the Jacobian is constant, so the metric is the *same* at every point and
collapses to ~9 doubles per element instead of `3·p·(p+1)²` vec3s — roughly **100× leaner at
p=4**. We do not exploit this. The 647M figure is therefore the *general-geometry* number, with
headroom left on the table for any structured or affine subregion. The affine-metric path is the
same lever that would also unlock batched-GEMM tensor cores (§2.7): both want the metric to be
one small constant per element, not a per-point array.

### 5.4 What about unstructured meshes and tets?

Two questions come up immediately: does the 647M DOF/GPU number assume a cube, and does any of
this work for tetrahedra?

**Distorted hexes: yes, for free.** The metric is built from the 8 corner coordinates via a
**trilinear Jacobian per element** (`mars_cvfem_ho_apply.hpp`), so it handles *any* warped,
sheared, or stretched hex — the cross-term coefficients `g[0]`/`g[1]` from §2.4 are exactly the
non-orthogonality. We do **not** assume or exploit cube geometry anywhere on the storage path.
So the scaling numbers in this chapter are the **fully-unstructured-hex** numbers: 647M DOF/GPU
transfers to a real warped hex mesh at zero extra cost (and is conservative — the affine lever
above would only help structured regions). The cube in the test driver is a convenience for
generating elements, not a shortcut the kernel relies on.

**Tets: no — they need a different kernel.** The entire speed argument rests on the **tensor-
product** structure: GLL nodes on a `(p+1)³` lattice are what let sum-factorization split a 3D
contraction into three 1D sweeps (§2.3). A tetrahedron has no such product structure, so this
kernel does not apply to tets at all. High-order tets need a separate scheme (collapsed-
coordinate / Bernstein bases), which is its own MARS track. This is not a MARS-specific
limitation — the same hex/tensor-product constraint binds MFEM's high-order partial assembly.

### 5.5 The straggler: when setup, not solve, kills the run

There is a subtler scaling wall that has nothing to do with the GPU. The DOF numbering —
assigning a global identity to every corner, edge, face, and interior node — originally ran **on
the host**, using `std::map` to deduplicate shared edges and faces. At ~10M elements/rank that is
~180M edge+face entries to insert, and it cost **~79 seconds per rank**.

In an MPI job, the slowest rank is the job. One unlucky rank with a heavier partition becomes a
*straggler* and stalls a 256-rank launch in a barrier before the solve even starts. The test
prints exactly this so the bottleneck is visible:

**Code jump** — `examples/distributed/unstructured/mars_ho_dist_apply_test.cu`, `runDistApply`:
```cpp
double tb0 = MPI_Wtime();
dof.buildDistributed(D.elemCorners, (long)D.nodeCount, P, D.cornerGid,
                     D.cornerOwner, D.elemOwner, rank, D.sharedCorner);
double tb1 = MPI_Wtime();
if (rank == 0)
    printf("[build] buildDistributed (host) %.3fs (numDof=%ld)\n", tb1-tb0, dof.numDof);
```
Notice: the numbering is wrapped in its own timer, separate from the apply. At scale this
`[build]` line, not the matvec, was the wall-clock the user stared at.

The fix is to recognize what the `std::map` dedup *is*: a sort + unique + binary-search lookup.
That is a textbook thrust pipeline. So we wrote a GPU-native numbering that does the identical
job on the device.

### 5.6 GPU-native numbering: sort, unique, scatter

The core move is to turn each element's 12 edges and 6 faces into a packed, sorted key, then let
thrust find the unique entities:

**Code jump** — `backend/distributed/unstructured/fem/mars_ho_dof_handler_gpu.hpp`,
`struct HoEntityKey`:
```cpp
// edge: hi = (c0<<32)|c1 (c0<=c1), lo = 0
// face: hi = (c0<<32)|c1, lo = (c2<<32)|c3 (c0<=c1<=c2<=c3)
struct HoEntityKey {
    uint64_t hi, lo;
    __host__ __device__ bool operator<(const HoEntityKey& o) const {
        return hi < o.hi || (hi == o.hi && lo < o.lo);
    }
};
```
Notice: an edge or face is identified by its **sorted** local corner ids, packed into two 64-bit
words. Sorting the corners makes the key orientation-independent — the same physical edge gets
the same key no matter which element or direction it is seen from. That is the GPU equivalent of
the host's `std::map<std::array<int,2>>` keyed on sorted corners.

Deduplication is then three thrust calls — emit, sort, unique:

**Code jump** — `backend/distributed/unstructured/fem/mars_ho_dof_handler_gpu.hpp`,
`buildDistributedGpu` (STAGE 1):
```cpp
// ... for_each emits 12 edge keys per element into d_edgeAll ...
thrust::sort(d_edgeAll.begin(), d_edgeAll.end());
d_uEdge.resize(d_edgeAll.size());
auto uend = thrust::unique_copy(d_edgeAll.begin(), d_edgeAll.end(), d_uEdge.begin());
nEdge = (long)(uend - d_uEdge.begin());
```
Notice: this is the whole dedup — sort brings identical keys adjacent, `unique_copy` collapses
them, and the count of survivors *is* the number of distinct edges. The same three lines repeat
for faces. What took ~79s in `std::map` is now a sort the GPU eats in milliseconds.

The final stage walks every (element, i, j, k) node, classifies it (corner / edge / face /
interior by how many indices are on a boundary), finds its entity id by binary search over the
unique keys, and **scatters** the per-DOF identity:

**Code jump** — `backend/distributed/unstructured/fem/mars_ho_dof_handler_gpu.hpp`,
`buildDistributedGpu` (STAGE 3):
```cpp
p_dofOwner[dofId] = owner;
p_dofKind[dofId]  = kind;
p_dofG0[dofId]    = gg0;  // ... g1,g2,g3,pos ...
if (shared)              atomicOr(&p_dofShared[dofId], 1);
if (boundary || shared)  atomicOr(&p_dofBoundary[dofId], 1);
```
Notice the race reasoning, which is the crux of why this is correct without locks: **many
elements write the same `dofId`, but they all write the same value** — owner, kind, key, and pos
are pure functions of *global* ids, so a plain store is race-safe (every writer agrees). Only the
boundary/shared *flags* are an OR across writers, so only those use `atomicOr`. This is the kind
of "what is actually shared between threads?" analysis that lets GPU code drop locks safely.

### 5.7 Trust but verify: the A/B self-check

A faster numbering is worthless if it is a *different* numbering. There is one subtlety: the host
assigns local edge/face ids in `std::map` *insertion* order, while the GPU assigns them in
*sorted-key* order. So the local DOF ids legitimately differ — by a permutation. The cross-rank
identity is the `DofKey`, never the local id. So the correct equivalence is not element-wise; it
is: **same counts, and same multiset of `DofKey`s**.

**Code jump** — `examples/distributed/unstructured/mars_ho_dist_apply_test.cu`,
`selfCheckNumbering`:
```cpp
std::sort(hk.begin(), hk.end());   // hk = host DofKeys packed to array<long,6>
std::sort(gk.begin(), gk.end());   // gk = gpu  DofKeys
long mismatch = 0;
for (size_t i = 0; i < hk.size(); ++i)
    if (hk[i] != gk[i]) { ++mismatch; ... }
```
Notice: it sorts *both* key lists and compares them as multisets — permutation-invariant by
construction. It also checks `numDof`/`nEdge`/`nFace`, the shared/boundary *totals* (not
per-id), and the owner histogram. All are quantities that survive a relabeling of local ids.

The result: **zero DofKey mismatches** at p=2/3/4, and the end-to-end A·1 = 0 apply gate passes
with the GPU numbering — proving it is not just structurally equal but *operationally*
identical. Critically, the GPU path is **opt-in**; the validated host numbering stays the
default:

**Code jump** — `examples/distributed/unstructured/mars_ho_dist_apply_test.cu`, `main`:
```cpp
const char* envGpu = std::getenv("MARS_HO_GPU_NUMBERING");
bool useGpu = cliGpu || (envGpu && std::atoi(envGpu) != 0);
Numbering mode = cliSelfCheck ? Numbering::SelfCheck
               : useGpu       ? Numbering::Gpu
                              : Numbering::Host;   // default unchanged
```
Notice: a new, faster path lands *alongside* the proven one, flag-gated
(`MARS_HO_GPU_NUMBERING=1` or `--gpu-numbering`), and ships with its own A/B harness. That is how
you remove a bottleneck at scale without betting the run on unverified code.

### 5.8 Takeaways

- **Weak-scale, and feed the GPU.** Apply weak-scales at ~98–100% to 40B DOF/64 GPUs; full
  matvec ~90%. The gap is the halo.
- **Communication self-resolves.** Surface/volume gives comm/compute ∼ V^(−1/3) (measured
  exponent −0.36); peer count caps at 26 in 3D. Comm falls 22%→4% just by growing per-GPU size.
- **Memory, not messages, is the wall.** The per-point metric `d_G` (~7.2 KB/elem at p=4) sets a
  ~647M DOF/GPU ceiling, measured by a checked-malloc OOM probe. High order raises this ceiling
  ~9× over stored CSR, putting a trillion DOF within ~1,550 GPUs.
- **Setup can be the straggler.** Host `std::map` numbering cost ~79s/rank and stalled a 256-rank
  launch. The GPU rewrite (sort + unique + scatter) is bit-equivalent (zero DofKey mismatches,
  A·1 passes), flag-gated, and leaves the host default untouched.
- **Trillion is in progress, not done.** The >4.3B-element domain-build crash turned out to be a
  2 GiB Allreduce message in cstone's global tree (not a count-value overflow); the fix is a
  coarser global `bucketSize` (`MARS_GLOBAL_BUCKETSIZE`), free for the apply. The remaining gap to
  a demonstrated trillion is the distributed *solve*, not the operator (Chapter 6).

---

## Chapter 6 — The Trillion Frontier (honest, in-progress)

We have an operator that is bit-exact, sum-factorized, and weak-scales to 40 billion DOF on 64
GPUs. The natural question for a PASC talk and a Gordon Bell submission is: *what stands between
us and a trillion?* This chapter is deliberately honest. Part of it is arithmetic that already
works. Part of it is a crash we hit at scale, and — more importantly — the **expert lesson** in
why the obvious fix was the wrong fix. We close with what is genuinely *not* done: the
distributed solve.

### 6.1 The arithmetic: why high order is the trillion path that *fits*

Start with the counting, because it reframes the whole problem.

A single hex element of order `p` carries `(p+1)^3` GLL nodes, but those nodes are shared with
neighbours. The *unique* DOF contributed per element, in the interior of a large mesh, is `p^3`
(one element "owns" its `p^3` lattice and shares the boundary lattice points). So:

| order p | DOF/element (p³) | elements for 10¹² DOF |
|---|---|---|
| 1 | 1 | ~10¹² (a trillion **elements**) |
| 2 | 8 | ~1.25 × 10¹¹ |
| 4 | 64 | **~1.56 × 10¹⁰** (15.6B elements) |
| 7 | 343 | ~2.9 × 10⁹ |

This is the central intuition of the whole project. **The element count — not the DOF count — is
what stresses the mesh infrastructure** (the octree, the SFC partition, the halo, the
per-element metric arrays). Raising `p` lets us hit a trillion DOF with *two orders of magnitude
fewer elements*. At p=4, a trillion DOF is only 15.6B elements. That is why the matrix-free HO
operator is the trillion path that actually fits Alps: at the measured 647M DOF/GPU ceiling, 10¹²
DOF needs ~1,550 GPUs (~390 nodes) — under half of Alps — whereas a p=1 *element*-trillion would
need ~9,300 GPUs.

The flat-throughput result from Chapter 2 is what makes this free: p=7 runs at essentially the
same GDOF/s as p=1, so packing more DOF into fewer elements costs the operator nothing.

### 6.2 Where the elements live: the cstone domain

Every MARS element is identified by a Hilbert SFC key, and the global decomposition is owned by a
`cstone::Domain`. The element-count → DOF-count translation above only matters because cstone is
what scales the *elements*.

**Code jump** — `backend/distributed/unstructured/domain.hpp`, the domain type alias:
```cpp
using DomainType = cstone::Domain<KeyType, RealType, AcceleratorTag>;
...
int bucketSize = 64, unsigned bucketSizeFocus = 8,
```
Notice two things. `KeyType` is the SFC key type (64-bit Hilbert keys — that part is *not* the
bottleneck). And `bucketSize = 64`: cstone subdivides any octree node holding more than 64
particles. Hold that number; it is the key to the whole crash story below.

cstone is a cosmology code. In its world "particle" = our "element" (one SFC key per element).
The leaf node counts it computes are how it decides where to refine the global tree and how to
balance work across ranks.

### 6.3 The crash: a >4.3B-global-element run dies inside the cstone domain build

When we pushed past **4.3 billion global elements**, the run crashed *inside* the cstone domain
build — not in our operator, not in our halo. 4.3B ≈ 2³², and that number is a giant flashing
arrow pointing at 32-bit integers. The tempting diagnosis writes itself: *"cstone counts
elements in `unsigned`; a global element count over 2³² overflows; migrate cstone to `uint64`
and move on."*

That diagnosis is wrong, and understanding *why* is the most transferable lesson in this entire
tutorial.

### 6.4 The expert lesson: read the library's design before patching it

cstone's `unsigned` (uint32) node counts are not an oversight. They are a **deliberate,
documented design contract**. The counts are a *refinement signal*, never an index.

**Code jump** — `cornerstone-octree/include/cstone/tree/csarray.hpp`, the `updateOctree` doc
comment:
```
 *    It is sensible to assume that the bucket size of the tree is much smaller than 2^32,
 *    and thus it is ok to use 32-bit integers for the node counts, because if the node count
 *    happens to be bigger than 2^32 for a node, this node will anyway be divided until the
 *    node count is smaller than the bucket size. We just have to make sure to prevent overflow,
 *    in MPI_Allreduce, therefore, maxCount should be set to 2^32/numRanks - 1 ...
```
Read that carefully. A node count is *only* ever compared against `bucketSize` (=64) to decide
*split or don't split*. If a node truly held more than 2³² elements, it would have been
subdivided long before — by construction no surviving leaf ever needs a count near 2³². The
count is a thermostat reading, not an address. A uint32 thermostat is fine.

**Code jump** — `cornerstone-octree/include/cstone/tree/csarray.hpp`, `calculateNodeCount`.
Notice the count is computed in `size_t` and then *deliberately saturated* down to the cap:
```cpp
size_t count = rangeEnd - rangeStart;
return stl::min(count, maxCount);   // saturate, do not overflow
```
And the cap travels through the MPI reduction with a custom saturating operator, so even the
*sum across ranks* cannot wrap:

**Code jump** — `cornerstone-octree/include/cstone/tree/update_mpi_gpu.cuh`, `sumCapped`:
```cpp
inline void sumCapped(void* inP, void* inoutP, int* len, MPI_Datatype*) {
    ...
    inout[i] = a + b >= std::max(a, b) ? a + b : std::numeric_limits<unsigned>::max();
}
```
That ternary *is* the contract: if `a + b` would wrap below `max(a,b)`, clamp to `2^32-1`
instead. The `expectOverflows` flag in `updateOctreeGlobalGpu` chooses this op precisely
*because* the design anticipates capped counts at scale (`maxCount = 2^32/numRanks - 1`). Three
independent places in cstone cooperate to make uint32 counts correct: per-node saturation, the
`maxCount` contract, and the saturating allreduce.

So the conclusion: **a blanket `uint32 → uint64` migration of cstone would not be a fix — it
would be a fork.** It would fight a design that is internally consistent and correct, balloon the
node-count memory across the whole tree, and break compatibility with upstream. We rejected the
sledgehammer.

The discipline here generalizes far past cstone. When a number near a power-of-two coincides with
a crash, the lazy move is to widen the type. The expert move is to **read the library's design
first**: is this quantity an *index* (must not saturate) or a *signal* (designed to saturate)? In
cstone it is a signal. The 2³² coincidence is a red herring; the count *values* are fine. The
real cause is one level up, and §6.5 has it.

### 6.5 The real blocker, and the one-line fix: a 2 GiB Allreduce message

Pinning the actual failure with `CUDA_LAUNCH_BLOCKING=1` put it on the **global-tree count
Allreduce** — the same `sumCapped` collective from §6.4, but the problem is the *message size*,
not the count values. The replicated global octree at a cube-2500 run has ~805M leaf nodes, each
a `uint32` count. That is `805M × 4 = 3.22 GB` in one `MPI_Allreduce` — past `2³¹` bytes. Cray
MPICH's chunked recursive-doubling implementation carries the per-chunk byte count in a 32-bit
field, so it truncates by one uint32 and corrupts the reduction. Note the distinction from §6.4:
each individual count is a fine uint32; it is the **number of nodes times 4 bytes** that overflows
a 32-bit *byte* counter inside MPICH. This is a **scale threshold on the global element count**
(it sets the global tree size), not a rank-count limit and not a `LocalIndex`/uint32-count
overflow on our side.

**Code jump** — `examples/distributed/unstructured/mars_ho_dist_apply_test.cu`, global bucketSize:
```cpp
// global bucketSize auto-scales with ncells^3 so the global tree stays under
// the ~2 GiB Allreduce limit; MARS_GLOBAL_BUCKETSIZE overrides.
```
The fix needs no cstone patch and no apply-side change: **make the global tree coarser** by
raising cstone's global `bucketSize`. A bigger bucket means fewer leaf nodes, so the count
message shrinks. We auto-scale it with `ncells³` (cube-2500 → bucketSize 128, message ~1.5 GB,
safely under 2 GiB), with `MARS_GLOBAL_BUCKETSIZE` as an env override. Crucially this is
**free for the operator**: the apply is entirely local, so a coarser *global* partition tree
does not touch matvec throughput at all — it only changes how the elements are bucketed for the
domain decomposition. The lesson for distributed scale: a collective that worked at a billion can
silently die at a few billion purely on **message bytes**, and the cheapest fix is often to make
the *replicated* structure coarser, not to widen a type.

### 6.6 Where MARS sits: the Gordon Bell yardstick

For calibration: MFEM won **Gordon Bell 2025** with **55.5 trillion DOF on 43,520 MI300A GPUs**
(~1.27B DOF/GPU); on Alps GH200 the comparable figure is ~9.3T DOF on ~9,200 GPUs. Our HO
matrix-free operator, at a measured **647M DOF/GPU** clean (p=4, the per-point metric array `d_G`
~7.2 KB/elem is the wall), sits squarely in the **~1B DOF/GPU class** — the same league as the
GB-winning code's per-GPU density. The affine-cube metric (512× less `d_G`) would push us to
~1–2B DOF/GPU. We are operator-competitive per GPU; the gap to a *demonstrated* trillion is the
solve below — not operator throughput, and no longer the domain build (§6.5).

A second calibration is on *method*, and it lands on the same choice MARS made. The **2025
Gordon Bell winner** (a Cascadia-tsunami digital twin, arXiv:2504.16344) discretizes its 3D
acoustic-gravity forward operator with **partial assembly** — verbatim, PA "stores an
asymptotically optimal amount of data: O(1) per degree of freedom," and the authors explicitly
chose it over the fully matrix-free option MFEM also supports ("matrix-free assembly where no
data is stored and all computations are done on the fly"), "for faster time-to-solution." That is
exactly the trade MARS makes in §2.4: store the metric once, apply via sum-factorization, never
recompute geometry on the hot path. So the operator class behind the last two Gordon Bell-scale
results — MFEM-PA and this one — is the same class as the MARS HO operator; the differentiators
are our CVFEM discretization and the distributed all-order reach, not the assembly strategy.

### 6.7 What is NOT done: the distributed solve

This is the honest boundary of the work. We have a distributed **matvec**. We do not yet have a
distributed **solve**.

A matvec is one application of the operator. A linear solve is hundreds of them wrapped in a
Krylov method with a preconditioner, and *each* piece must be made distributed and
halo-correct. What exists today versus what remains:

- **Done:** distributed HO DOF numbering (`HODofHandler::buildDistributed`, with the
  min-rank-among-holders ownership from Chapter 3), the receiver-driven halo (`HoHalo`), and the
  full distributed matvec `forward-halo(x) → apply OWNED elements → reverse-add(y)`, validated by
  `A·1 = 0` at ~1e-17 on 4 ranks. The GPU-native numbering path (`buildDistributedGpu`) removes
  the host straggler.
- **Not done — the roadmap:**
  1. **FlexGMRES** over the distributed operator. The matvec is ready; what's missing is the
     distributed inner products (`MPI_Allreduce` over *owned* DOF only, using the same ownership
     map that made `A·1=0` correct) and the Arnoldi/restart machinery on top.
  2. **A low-order preconditioner.** High-order matrix-free operators are notoriously
     ill-conditioned; FlexGMRES alone will stall. The standard remedy is a **low-order-refined
     (LOR) AMG** preconditioner: build a p=1 element on the *same* GLL nodes, assemble that sparse
     operator (cheap and AMG-friendly), and use AMG on it to precondition the high-order
     matrix-free solve. None of this is built yet.
  3. **Robust domain build past 4.3B elements** — the 2 GiB Allreduce blocker is fixed with the
     auto-scaling global `bucketSize` (§6.5); larger-scale validation runs are the remaining work.

The intellectually honest framing for the talk: MARS has a **trillion-class operator** —
bit-exact, sum-factorized, weak-scaling, in the ~1B DOF/GPU density band of the Gordon Bell
winner. The trillion-DOF *solve* is the next milestone, and its critical path is FlexGMRES +
LOR-AMG over the operator we already trust, plus hardening the cstone domain build at the
4-billion-element frontier. **We claim the operator. We do not yet claim the trillion.**

---

## Glossary

- **DOF (degree of freedom)** — one scalar unknown the solver computes. A `p=4` hex carries
  `(p+1)^3 = 125` DOF; many are shared with neighbours.
- **p (polynomial order)** — the degree of the polynomial inside each element. `p=1` is linear;
  higher `p` curves within an element and converges at roughly order `p+1`.
- **GLL (Gauss-Lobatto-Legendre) nodes** — the `(p+1)` non-uniform points per direction where
  high-order solution DOF live; chosen for interpolation stability.
- **matvec** — the operation `y = A·x`, the inner loop of every Krylov solver. The whole tutorial
  is about doing it fast.
- **matrix-free** — computing `A·x` by recomputing element contributions on the fly, never
  storing `A`. Trades a little arithmetic for a large memory saving; the win grows with `p`.
- **sum-factorization** — exploiting the tensor-product structure of an element to do a 3D
  contraction as three 1D sweeps, turning `O(p⁶)` per element into `O(p⁴)`. Keeps throughput flat
  across order.
- **halo** — the boundary DOF a rank does not own but needs (ghosts), plus the exchange that
  fills them (`forward`) and returns their contributions (`reverseAdd`).
- **weak-scaling** — adding GPUs while holding work-per-GPU fixed, so the global problem grows.
  The honest frame for problems too big for one GPU. (Contrast: *strong-scaling* fixes the global
  problem and adds GPUs.)

---

## Try it yourself

The distributed apply test (`mars_ho_dist_apply_test`) is the live driver behind every gate in
this tutorial. Build it with the FEM examples (see `CLAUDE.md`), then run on Alps with `srun`.
The binary path is relative from the build directory.

**1. The headline gate — A·1 = 0 across ranks (host numbering, the default):**
```bash
srun -N1 -n4 ./examples/distributed/unstructured/mars_ho_dist_apply_test --ncells=16 --p=2
```
Expect a `[build] buildDistributed (host) ...` line and a `max|A.1| over owned ... PASS` line at
~1e-17. Try `--p=3` and `--p=4` to reproduce the `1.7e-18 / 2.5e-17 / 2.8e-17` ladder.

**2. The GPU-native numbering path (opt-in, leaves host default untouched):**
```bash
# via environment flag
MARS_HO_GPU_NUMBERING=1 srun -N1 -n4 ./examples/distributed/unstructured/mars_ho_dist_apply_test --ncells=16 --p=2
# or via CLI flag
srun -N1 -n4 ./examples/distributed/unstructured/mars_ho_dist_apply_test --gpu-numbering --ncells=16 --p=2
```

**3. The A/B self-check — host vs GPU numbering bit-equivalence:**
```bash
srun -N1 -n4 ./examples/distributed/unstructured/mars_ho_dist_apply_test --self-check --ncells=16 --p=2
```
Expect `[self-check] DofKey multiset: 0 / N mismatched ... MATCH`, matching `numDof/nEdge/nFace`
and owner histogram, a `build time: host ... gpu ... speedup ...x` line, and then the A·1 gate
re-run with the GPU numbering.

**4. Comm self-resolution — grow per-GPU size and watch comm fraction fall.** Increase `--ncells`
(and ranks/nodes) so per-GPU DOF climbs; the final line reports measured `commFrac`. Across the
ladder it falls 22% → 15% → 4% as per-GPU DOF grows 4M → 11M → 625M.

**5. Probe the per-GPU memory ceiling.** Push `--ncells`/`--p` until you hit the checked-malloc
OOM. Instead of a segfault you get
`p=... DEVICE OOM (needs ~X GB, d_G=Y GB) -- per-GPU ceiling exceeded`, attributing the wall to
`d_G`. Sweeping upward until this fires is exactly how the 647M DOF/GPU figure was measured.

> Note (Alps run convention): on multi-GPU nodes each rank binds one device
> (`cudaSetDevice(rank % devCount)`), and `srun` needs `--export=ALL` for the environment flags to
> reach the ranks.

---

## References

The MARS HO operator is the tensor-product / sum-factorized partial-assembly method, with a CVFEM
discretization. The references below are the sources the tutorial cites directly.

- **[Knaus 2022]** R. Knaus, "A fast matrix-free approach to the high-order control volume finite
  element method with application to low-Mach flow," *Computers & Fluids*, vol. 239, art. 105408,
  2022. DOI: [10.1016/j.compfluid.2022.105408](https://doi.org/10.1016/j.compfluid.2022.105408).
  Sandia report SAND2022-3366J (OSTI: <https://www.osti.gov/biblio/1870437>). The high-order CVFEM
  matrix-free method this operator is based on: the tensor-product (sum-factorized) element
  residual evaluation and the stored per-element metric `G`.

- **[MFEM partial assembly]** the partial-assembly / sum-factorization approach for tensor-product
  high-order elements:
  - MFEM performance & partial assembly page: <https://mfem.org/performance/>.
  - J. Andrej, N. Atallah, J.-P. Bäcker, J. Camier, D. Copeland, V. Dobrev, Y. Dudouit, T. Duswald,
    B. Keith, D. Kim, T. Kolev, B. Lazarov, K. Mittal, W. Pazner, S. Petrides, S. Shiraiwa,
    M. Stowell, V. Tomov, "High-performance finite elements with MFEM," arXiv:2402.15940, 2024.
    <https://arxiv.org/abs/2402.15940>.
  - R. Anderson et al., "MFEM: A Modular Finite Element Methods Library," *Computers & Mathematics
    with Applications*, vol. 81, pp. 42–74, 2021. DOI:
    [10.1016/j.camwa.2020.06.009](https://doi.org/10.1016/j.camwa.2020.06.009).
  - (background) J. Brown et al., "libCEED: Fast algebra for high-order element-based
    discretizations," *Journal of Open Source Software*, vol. 6, no. 63, art. 2945, 2021. DOI:
    [10.21105/joss.02945](https://doi.org/10.21105/joss.02945).

- **[GB 2025]** S. Henneking, S. Venkat, V. Dobrev, J. Camier, T. Kolev, M. Fernando,
  A.-A. Gabriel, O. Ghattas, "Real-time Bayesian inference at extreme scale: A digital twin for
  tsunami early warning applied to the Cascadia subduction zone," arXiv:2504.16344, 2025
  (2025 ACM Gordon Bell Prize). <https://arxiv.org/abs/2504.16344>. The MFEM-PA forward operator
  cited for the per-GPU-density and method comparison in Chapters 1 and 6.

- **[FP64 tensor cores]** J. Tu, I. Karlin, J. Camier, V. Dobrev, T. Kolev, S. Henneking,
  O. Ghattas, "Accelerating High-Order Finite Element Simulations at Extreme Scale with FP64 Tensor
  Cores," arXiv:2603.09038, 2026. <https://arxiv.org/abs/2603.09038>. Context for the FP64 DMMA
  discussion in §2.7.
