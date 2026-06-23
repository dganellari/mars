# Talk-Prep Brief — Distributed High-Order Matrix-Free CVFEM (PASC / Gordon Bell)

This is the presenter's cheat sheet for the matrix-free scaling talk. It gives the slide arc, the
exact numbers to say out loud, the one honest framing of trillion status, the negative result to
own up front, and a prepared Q&A. Companion deep-dive: `docs/TUTORIAL_matfree.md`. Figures:
`docs/figures/make_matfree_figs.py` (fig_memory / fig_throughput / fig_convergence /
fig_comm_pergpu / fig_trillion / fig_ho_scaling).

---

## The 60-second elevator version (say this if you have one minute)

"Every PDE solver does one thing in its inner loop: `y = A·x`. We built that operation
matrix-free — never storing the matrix — at high polynomial order, on GPUs. Sum-factorization
makes it run at the same throughput for p=7 as for p=1, while high order gives ~570× more accuracy
per unknown. We made it distributed without losing a bit, weak-scaled it to 40 billion DOF on 64
GH200s, and showed that communication self-resolves as you grow each GPU's job. The operator sits
in the ~1-billion-DOF-per-GPU class — the same league as the 2025 Gordon Bell winner. The
trillion-DOF *solve* is the next milestone; we are not claiming a trillion yet."

---

## Slide narrative arc (5 slides + title)

The arc is a logic chain: **high order is worth it → matrix-free is how you afford it → it stays
fast → it distributes for free → it scales toward a trillion.** Each slide answers the objection
the previous one raises.

### Slide 1 — Memory crossover: *why matrix-free at all* (fig_memory)
- The objection it kills: "high order is too expensive to store."
- Show assembled bytes/DOF *growing* with p (324 → 40,500) vs matrix-free *shrinking*
  (421 → 121).
- **Say:** "At p=7, matrix-free is 336× leaner per DOF. The assembled matrix gets worse with
  order; matrix-free gets better."
- One honest aside: at p=1 matrix-free is slightly *worse* (the 0.8× crossover). Don't hide it —
  it makes the p=7 win credible.

### Slide 2 — Flat throughput: *high order is free to run* (fig_throughput)
- The objection it kills: "high order must be slow."
- Show MDOF/s essentially flat across p=1..7 (~4794 at p=1, ~4530 at p=7).
- **Say:** "p=7 runs at the same rate as p=1. That's sum-factorization — O(p⁴), not O(p⁶)."
- This is the slide where you name **sum-factorization** as the hero.

### Slide 3 — Convergence: *high order is worth it* (fig_convergence)
- The payoff for slides 1–2.
- Show: at ~2197 DOF, p=1 error 2.0e-2 vs p=4 error 3.5e-5.
- **Say:** "Same number of unknowns, ~570× lower error. That accuracy is what you're buying with
  the order — and slides 1 and 2 just showed it's cheap to store and free to run."

### Slide 4 — Weak scaling + comm self-resolves: *it distributes for free* (fig_ho_scaling +
fig_comm_pergpu)
- Show apply weak-scaling ~98–100%, full matvec ~78% at small per-GPU rising to ~90% at 40B/64
  GPUs.
- Show comm fraction falling 52% → 19% (measured) as per-GPU DOF grows, fitted exponent −0.36 ≈
  the V^(−1/3) surface/volume law.
- **Say:** "The operator doesn't notice it's distributed — apply weak-scales at ~98%. The halo is
  the only gap, and it *self-resolves*: bigger per-GPU blocks have less relative communication,
  because surface-to-volume shrinks like V^(−1/3) and the neighbour count saturates at 26 in 3D.
  Don't fight the halo — feed the GPU."

### Slide 5 — Trillion path + per-GPU ceiling: *where this is going* (fig_trillion)
- Show per-GPU ceilings: CSR p1 70M → MF p1 77M → HO MF p4 **647M**.
- Show the consequence: 10¹² DOF in ~1,550 GPUs (~390 nodes, <half of Alps) vs ~9,300 for p1
  elements.
- **Say:** "The wall isn't communication — it's HBM, specifically the per-point metric array. At
  647M DOF/GPU, a trillion DOF needs under half of Alps. High order raised that ceiling ~9× over
  stored CSR." Then deliver the honest status (next section) immediately.

---

## Exact headline numbers (say these, don't round away the precision)

- **Accuracy:** ~570× more accurate per DOF, p=1 → p=4 (at ~2197 DOF: 2.0e-2 vs 3.5e-5). ~p+1
  convergence order.
- **Memory crossover:** 0.8× at p=1 (slightly worse), **336× leaner at p=7** (matrix-free vs
  assembled bytes/DOF).
- **Throughput:** flat **~4–8 GDOF/s/GPU** across p=1..7 (p=7 ≈ p=1). Sum-factorization =
  O(p⁴) not O(p⁶).
- **Correctness gate:** A·1 = 0 at **1.7e-18 / 2.5e-17 / 2.8e-17** for p=2/3/4 on 4 ranks; device
  halo bit-identical to host.
- **Scaling:** **40 billion DOF on 64 GH200s** (p=4). Apply ~6 GDOF/s/GPU = the single-GPU rate
  (~100% apply weak-scaling). Full matvec ~90% efficient.
- **Comm self-resolution:** 22% → 15% → 4% as per-GPU DOF grows 4M → 11M → 625M (measured
  exponent −0.36).
- **Per-GPU ceiling:** **647M DOF/GPU at p=4** (set by the metric array `d_G`, ~7.2 KB/elem).
- **Trillion economics:** 10¹² DOF (p=4) = 15.6B elements ≈ **~1,550 GPUs (~390 nodes)** vs
  ~9,300 GPUs for p=1 elements.
- **GB yardstick:** MFEM won Gordon Bell 2025 at **55.5T DOF on 43,520 MI300A GPUs (~1.27B
  DOF/GPU)**; on Alps that class ≈ 9.3T DOF on ~9,200 GH200s. MARS HO operator = **~1B DOF/GPU
  class**.
- **Host straggler fix:** host `std::map` numbering ~79s/rank → GPU-native (sort+unique+scatter),
  **zero DofKey mismatches**, flag-gated, host default unchanged.

---

## The ONE honest framing of trillion status

Use this exact framing whenever trillion comes up. Do not let a slide or a question drift it.

> **"We have a trillion-CLASS operator, not a trillion-DOF result. The operator is bit-exact,
> sum-factorized, and weak-scales to 40 billion DOF on 64 GPUs at the ~1-billion-DOF-per-GPU
> density of the Gordon Bell winner. The trillion-DOF run is in progress: a beyond-4.3-billion-
> element configuration currently crashes inside the cstone domain build, and we're still pinning
> the cause. We claim the operator and the scaling. We do not claim a trillion."**

What this prevents: someone tweeting "MARS did a trillion DOF." You did not, and the figure
(fig_trillion) is a *projection* from a measured ceiling, not a measured run. Always say
"projected ~1,550 GPUs," never "we ran it on 1,550 GPUs."

---

## The negative result to own (don't let a reviewer find it for you)

**FP64 tensor cores (DMMA) do not help this operator.** Put it on a slide or say it in the comm
slide. Framing:

> "We built the FP64 tensor-core path and measured it. It doesn't help — and that's not a failure,
> it's a diagnosis. This operator is shared-memory-bandwidth and occupancy bound, not compute
> bound. Sum-factorization already made the arithmetic cheap; tensor cores accelerate the part
> that was already cheap and add register/shared pressure that *hurts* occupancy. The DMMA
> Galerkin variant hit 90 KB smem/block → 2 blocks/SM → ~12.5% occupancy. Knowing which resource
> binds you is the engineering win; throwing the shiny feature at the wrong bottleneck is the
> trap."

Owning this *increases* credibility — it signals you measured rather than assumed, and it
pre-empts the inevitable "did you try tensor cores?" question.

---

## Anticipated Q&A

**Q: How do you compare to MFEM / why should I care given they won Gordon Bell?**
A: MFEM is the yardstick, and we cite it as such. They demonstrated 55.5T DOF at ~1.27B DOF/GPU
on MI300A. Our *operator* sits in the same ~1B DOF/GPU density class on GH200 — same league per
GPU. The difference is honest scope: they have a full demonstrated solve at trillion+; we have a
demonstrated operator + matvec to 40B and a trillion *projection*. We're not claiming to have beaten
them; we're showing a CVFEM matrix-free operator is competitive per-GPU and laying out the path to
the solve.

**Q: Why CVFEM and not standard Galerkin FEM?**
A: CVFEM enforces conservation (mass/momentum) on control volumes around each node — it does not
lose or invent mass, which is exactly what you need for the fluid problems this targets (pumps,
channels, the coupled NS work). We follow Knaus's high-order CVFEM formulation. The sum-
factorization and matrix-free machinery is the same idea you'd use in spectral-element Galerkin;
the conservation property is why we chose the control-volume form. (Honest aside if pressed: high-
order CVFEM is less standard than HO Galerkin, so the operator-design + correctness gates were
non-trivial — hence the A·1=0 / A·linear=0 patch tests.)

**Q: Is the solve done? Can you actually solve a system, not just apply the operator?**
A: No — and we're explicit about that. We have a distributed *matvec*, validated by A·1=0 to
machine zero on 4 ranks. A *solve* is hundreds of matvecs in a Krylov method with a
preconditioner. The roadmap is FlexGMRES (distributed inner products over owned DOF, reusing the
ownership map that made A·1=0 correct) plus a low-order-refined (LOR) AMG preconditioner, because
HO matrix-free operators are ill-conditioned and GMRES alone will stall. None of that is built
yet. The operator is the trustworthy foundation it sits on.

**Q: What's the per-GPU limit, and what sets it?**
A: ~647M DOF/GPU at p=4 on GH200, and it's memory, not compute or communication. The matrix-free
apply stores one geometric array — the per-point metric `d_G`, ~7.2 KB/element at p=4. That, not
the messages, is the wall. We measured it with a checked-malloc OOM probe that attributes the
ceiling to `d_G` specifically. Good news for trillion: high order raised this ceiling ~9× over
stored CSR. (If pressed on headroom: an affine-cube metric stores 512× less and would push us to
~1–2B DOF/GPU, but that only applies to axis-aligned elements.)

**Q: 4.3 billion ≈ 2³². Why not just make cstone use uint64 and be done?**
A: Because cstone's uint32 counts are correct *by design*, not by accident. The node count is a
refinement *signal* — it's only ever compared to the bucket size (64) to decide split-or-not. If a
node held more than 2³² particles it would have been subdivided long before; no surviving leaf ever
needs a count near 2³². cstone *deliberately* saturates: per-node `min(count, maxCount)`, a
documented `maxCount = 2^32/numRanks − 1` contract, and a custom saturating MPI allreduce
(`sumCapped`). A blanket uint64 migration wouldn't be a fix — it'd be a fork that fights a
consistent design and balloons tree memory. The 2³² coincidence is a red herring; the real crash
is something on *our* side of the boundary (likely a genuine 32-bit index/allocation in the domain
build), which we're pinning with CUDA_LAUNCH_BLOCKING. The lesson: read whether a quantity is an
*index* or a *signal* before widening its type.

**Q (if asked): What changed to make this run at scale at all?**
A: Two memory disciplines. (1) The halo keys only *boundary* DOF, not all DOF — that's O(surface)
not O(volume); keying all DOF would build a ~78 GB host map at 650M DOF/GPU. (2) The DOF numbering
moved from host `std::map` (~79s/rank, a straggler that stalled 256-rank launches) to a GPU-native
sort+unique+scatter, validated bit-equivalent (zero DofKey mismatches) and flag-gated so the proven
host path stays default.

---

## Landmines to avoid on stage

- Do **not** say "we did a trillion." Say "trillion-class operator" / "projected to ~1,550 GPUs."
- Do **not** oversell tensor cores — own the negative result instead.
- Do **not** present fig_trillion as a measured run — it's a projection from the 647M ceiling.
- Do **not** claim a solve — only the operator + matvec are validated.
- If the A·1 numbers come up, they are machine-zero (1e-17), i.e. distribution reproduces the
  single-rank result *to the last bit* — that's the strong claim you *can* make confidently.
