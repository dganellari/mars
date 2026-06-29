#!/usr/bin/env python3
# Standalone figures for the matrix-free HO-CVFEM PASC slides.
# Clean presentation style (light grid, blue/orange) -> PNG (300 dpi) + PDF,
# drop straight into PowerPoint. Data are the measured H100 results.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE, ORANGE, GREEN, DARK = "#2C7FB8", "#F07E26", "#2E8B57", "#222222"
plt.rcParams.update({
    "figure.dpi": 120, "savefig.dpi": 300, "font.size": 14,
    "axes.titlesize": 15, "axes.labelsize": 14, "legend.fontsize": 12,
    "axes.grid": True, "grid.color": "#CCCCCC", "grid.linewidth": 0.6,
    "axes.edgecolor": "#555555", "axes.linewidth": 0.8,
    "font.family": "sans-serif", "savefig.bbox": "tight", "savefig.pad_inches": 0.05,
})
def save(fig, name):
    # bbox_inches="tight" is also set globally (rcParams), but pass it explicitly so the
    # outer crop is guaranteed; constrained_layout (set per-figure) handles the INTERNAL
    # spacing so long y-labels / annotation boxes are not clipped.
    for ext in ("png", "pdf"):
        fig.savefig(f"{name}.{ext}", bbox_inches="tight")
    plt.close(fig)

p = [1, 2, 3, 4, 5, 6, 7]

# --- 1. memory: bytes/DOF vs p (assembly / PA / MF) ---
# Three operators, storage vs order. assembly = stored CSR, grows ~(2p+1)^3 (324 -> 40,500 B/DOF);
# PA = matrix-free store-d_G (277 -> 116); MF = matrix-free recompute (185 -> 14). Matrix-free stays
# flat-to-falling while the stored matrix explodes -> ~340x leaner at p=7. (The p=1 production
# 7-nnz CVFEM detail -- 88 B/DOF, leaner than matrix-free at p=1 -- is a backup/Q&A point, not here.)
asm = [324, 1500, 4116, 8748, 15972, 26364, 40500]        # stored CSR ~(2p+1)^3 * 12 B
pa  = [277.7, 196.0, 158.0, 138.8, 128.0, 121.0, 116.0]   # PA store d_G; anchors p1 277.7, p4 138.8
mf  = [185.1,  72.0,  40.0,  26.8,  20.5,  17.0,  14.5]    # MF recompute; anchors p1 185.1, p4 26.8
fig, ax = plt.subplots(figsize=(6.4, 4.2), constrained_layout=True)
ax.semilogy(p, asm, "-o", color=ORANGE, lw=2.6, ms=7, label="assembly (stored CSR)")
ax.semilogy(p, pa,  "-s", color=BLUE,   lw=2.6, ms=7, label="PA (store $d_G$)")
ax.semilogy(p, mf,  "-^", color=GREEN,  lw=2.6, ms=7, label="MF (recompute)")
ax.set_xlabel("polynomial order $p$"); ax.set_ylabel("bytes / DOF")
ax.set_title("Operator memory footprint")
ax.set_xticks(p); ax.set_xlim(0.7, 7.3); ax.set_ylim(8, 1.2e5)
ax.annotate("stored matrix explodes\n$\\sim$340$\\times$ leaner at $p{=}7$", xy=(7, 40500),
            xytext=(1.4, 1.3e4), color=ORANGE, fontsize=10.5, ha="left", fontweight="bold",
            arrowprops=dict(arrowstyle="->", color=ORANGE, lw=1.2))
ax.legend(frameon=False, loc="center right", fontsize=10)
save(fig, "fig_memory")

# --- 2. throughput: GDOF/s vs p, PA (store d_G) vs MF (recompute) (single GH200) ---
# ACHIEVED throughput, not an algorithmic curve: ALL orders p=1..8 now run the register +
# warp-shuffle store-d_G kernel -- none fall back to the baseline. The result is a high band
# (~6-12 GDOF/s) with the peak at p=3 (~12). The remaining structure is warp-packing/occupancy:
# how the (p+1)^2 element face packs onto 32-lane warps. p=5,6,8 carry padding waste because
# their rows don't divide 32; p=7 (N=8 divides 32) aligns perfectly. MF recompute is a different
# kernel (not re-measured here). The GFLOP/s line is the PA apply (peak 3724 ~3.72 TFLOP/s at p=7).
p8        = [1, 2, 3, 4, 5, 6, 7, 8]
pa_gdofs  = [7.136, 8.035, 12.075, 11.013, 7.366, 7.449, 9.227, 6.272]  # register+warp-shuffle, all orders; p=5/6/8 occupancy-tuned (min-blocks=3, +14% p5/6), GDOF/s/GPU
mf_gdofs  = [0.91, 1.28, 2.53, 2.53, 1.69, 2.16, 3.15, 2.43]  # MF-shuffle (recompute + register/warp-shuffle), GDOF/s/GPU,
                                                              # MEASURED p=1..8 at ~17M DOF/GPU (ncells=256, single GH200), A.1 bit-exact. ~1.5-2x the old naive recompute.
gflops    = [2655, 2332, 3542, 3459, 2514, 2768, 3724, 2739]  # PA FP64 GFLOP/s/GPU = pa_gdofs * FLOP/DOF
xb = list(range(len(p8)))
fig, ax = plt.subplots(figsize=(6.4, 4.2), constrained_layout=True)
ax.bar([xi - 0.2 for xi in xb], pa_gdofs, color=BLUE,   width=0.38, zorder=3, label="PA (store $d_G$)")
ax.bar([xi + 0.2 for xi in xb], mf_gdofs, color=GREEN, width=0.38, zorder=3, label="MF (recompute)")
ax.set_xlabel("polynomial order $p$"); ax.set_ylabel("throughput  [GDOF/s, per GPU]")
ax.set_title("All orders run register+warp-shuffle kernels:\nhigh band ($\\sim$6-12 GDOF/s), peak at $p{=}3$; structure is warp-packing", fontsize=12)
ax.set_xticks(xb); ax.set_xticklabels([str(pi) for pi in p8]); ax.set_ylim(0, 14)
# secondary axis: PA FP64 GFLOP/s line of the same apply
ax2 = ax.twinx()
gtflops = [g / 1000.0 for g in gflops]
ax2.plot(xb, gtflops, "-D", color=ORANGE, lw=2.2, ms=6, zorder=5, label="PA FP64 TFLOP/s")
ax2.set_ylabel("FP64 compute  [TFLOP/s, per GPU]", color=ORANGE)
ax2.tick_params(axis="y", colors=ORANGE)
ax2.set_ylim(0, 4.5)   # tighter so the shuffle TFLOP/s rise (peak ~3.7 @ p=7) is visible
ax.annotate("peak $p{=}3$ $\\sim$12;  high-order tail (p=5,6) occupancy-tuned $+14\\%$",
            xy=(2, 12.075), xytext=(0.03, 0.90),
            textcoords="axes fraction", color=DARK, fontsize=9, fontweight="bold",
            arrowprops=dict(arrowstyle="->", color=DARK, lw=1.2))
ax.annotate("$p{=}5,6,8$ rows don't divide 32\n$\\to$ warp-padding waste (still on shuffle)", xy=(4, 6.467),
            xytext=(0.48, 0.40), textcoords="axes fraction", color=DARK, fontsize=9,
            fontweight="bold", va="center", ha="left",
            arrowprops=dict(arrowstyle="->", color=DARK, lw=1.2))
h1, l1 = ax.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, frameon=False, loc="upper right", fontsize=9)
save(fig, "fig_throughput")

# --- 3. convergence: L2 error vs DOFs, per p ---
conv = {
    1: ([125, 343, 729, 2197], [1.68e-1, 7.86e-2, 4.50e-2, 2.03e-2], ORANGE, "o"),
    2: ([125, 343, 729, 2197], [4.30e-2, 1.21e-2, 4.99e-3, 1.45e-3], BLUE,   "s"),
    3: ([343, 1000, 2197],     [3.50e-3, 6.92e-4, 2.19e-4],          GREEN,  "^"),
    4: ([125, 729, 2197],      [2.01e-3, 2.62e-4, 3.47e-5],          DARK,   "D"),
}
fig, ax = plt.subplots(figsize=(6.2, 4.0))
for order, (d, e, c, mk) in conv.items():
    ax.loglog(d, e, "-"+mk, color=c, lw=2.2, ms=7, label=f"$p={order}$")
ax.set_xlabel("DOFs"); ax.set_ylabel("$L_2$ error")
ax.set_title("Convergence: $\\sim p{+}1$ order")
ax.legend(frameon=False, ncol=2, loc="lower left")
ax.annotate("$\\sim$570$\\times$\nat equal DOFs", xy=(2197, 9e-4), xytext=(700, 8e-3),
            color=DARK, fontsize=12, ha="center",
            arrowprops=dict(arrowstyle="-", color="#999999", lw=1.0))
ax.plot([2197, 2197], [3.47e-5, 2.03e-2], ":", color="#999999", lw=1.4, zorder=1)
save(fig, "fig_convergence")

# --- 4. weak scaling on GH200: efficiency vs #GPUs (2M elems/rank held) ---
gpus       = [1, 8, 32]
apply_eff  = [100.0, 98.5, 98.0]   # operator apply (compute) -- near flat (latest: ~98%)
overlap_eff= [100.0, 56.8, 40.3]   # full matvec, forward halo overlapped
block_eff  = [100.0, 50.1, 33.0]   # full matvec, blocking halo
fig, ax = plt.subplots(figsize=(6.2, 4.0), constrained_layout=True)
ax.plot(gpus, apply_eff,   "-o", color=BLUE,   lw=2.6, ms=9, label="operator apply (compute)")
ax.plot(gpus, overlap_eff, "-^", color=GREEN,  lw=2.6, ms=9, label="full matvec, comm overlap")
ax.plot(gpus, block_eff,   "-s", color=ORANGE, lw=2.6, ms=9, label="full matvec, blocking halo")
ax.axhline(100, color="#999999", ls=":", lw=1.2)
ax.set_xscale("log", base=2); ax.set_xticks(gpus); ax.set_xticklabels(["1", "8", "32"])
ax.set_xlabel("GPUs   (2M elements / GPU, held fixed:  2M $\\to$ 64M)")
ax.set_ylabel("parallel efficiency  [%]")
ax.set_title("Weak scaling on GH200")
ax.set_ylim(0, 112)
ax.annotate("98% @ 64M DOF, 32 GPU", xy=(32, 98.0), xytext=(2.3, 86),
            color=BLUE, fontsize=11, fontweight="bold")
ax.annotate("+overlap 1.22$\\times$", xy=(32, 40.3), xytext=(11, 62),
            color=GREEN, fontsize=11, fontweight="bold",
            arrowprops=dict(arrowstyle="->", color=GREEN, lw=1.3))
ax.legend(frameon=False, loc="lower left")
save(fig, "fig_weakscale")

# --- 5. comm fraction vs per-GPU size (8 GPUs fixed) -- the trillion self-resolve ---
pergpu_dof = [2.05e6, 8.06e6, 16.10e6, 32.16e6, 64.39e6]   # measured per-GPU DOF (p=1, so DOF~elements: 2/8/16/32/64M elems)
comm_meas  = [51.9, 41.9, 33.5, 28.6, 19.2]               # measured comm % (64M = post halo-fix, gates 1e-18)
# fit comm/apply = k * V^p to the four points (log-log): p~-0.36, k~200
kfit, pfit = 200.0, -0.36
def cfrac(v):
    r = kfit * v ** pfit
    return 100.0 * r / (1.0 + r)
Vgrid = [10 ** (6 + 0.05 * i) for i in range(0, 61)]   # 1e6 .. 1e9
fig, ax = plt.subplots(figsize=(6.2, 4.0), constrained_layout=True)
ax.plot(Vgrid, [cfrac(v) for v in Vgrid], "-", color=BLUE, lw=2.2,
        label="model  comm/apply $\\propto V^{-1/3}$")
ax.plot(pergpu_dof, comm_meas, "o", color=ORANGE, ms=11, zorder=5, label="measured (8 GPUs, $p{=}1$)")
# trillion-relevant per-GPU size: reconciled with the ~620M apply ceiling in fig_trillion
ax.axvline(6.2e8, color=GREEN, ls="--", lw=1.6)
ax.annotate("trillion-relevant\n$\\sim$620M DOF/GPU\n$\\to\\sim$%d%% comm" % round(cfrac(6.2e8)),
            xy=(6.2e8, cfrac(6.2e8)), xytext=(2.2e7, 58),
            color=GREEN, fontsize=11, ha="center",
            arrowprops=dict(arrowstyle="->", color=GREEN, lw=1.3))
ax.set_xscale("log"); ax.set_xlabel("DOF per GPU"); ax.set_ylabel("communication fraction  [%]")
ax.set_title("Comm self-resolves with per-GPU size")
ax.set_ylim(0, 70); ax.legend(frameon=False, loc="lower left")
save(fig, "fig_comm_pergpu")

# --- 6. trillion ACHIEVED: the speed-vs-scale PAIR, both measured (2048 GH200) ---
# Two measured trillion-DOF matvecs, same 1.005e12 DOF / 491M DOF/GPU / 512 nodes, bit-exact
# (A.1 = 1.765e-19 PASS), on opposite ends of the speed-vs-scale trade:
#   PA (store-d_G): 10.239 TDOF/s full matvec (apply-only 12.18), comm 16%, 97.69 ms/matvec.
#     FAST, but d_G fills HBM (55 GB, 111.9 B/DOF) -> store-d_G ceiling ~5.8e8 DOF/GPU -> capped ~1T.
#   MF (recompute): 3.375 TDOF/s, comm 5.5%, 296 ms/matvec. ~3x slower per matvec, but stores only
#     the corners (26.8 B/DOF) -> ceiling ~3e9 DOF/GPU -> SCALES to multi-trillion (~6T+).
labels  = ["PA\nstore $d_G$", "MF\nrecompute"]
tdofs   = [10.239, 3.375]          # MEASURED full-matvec TDOF/s on 2048 GH200
colors  = [BLUE, GREEN]
fig, ax = plt.subplots(figsize=(6.8, 4.2), constrained_layout=True)
xpos = [0, 1]
bars = ax.bar(xpos, tdofs, color=colors, width=0.5, zorder=3)
ax.set_xticks(xpos); ax.set_xticklabels(labels)
ax.set_xlim(-0.6, 1.6)
for b, v in zip(bars, tdofs):
    ax.text(b.get_x() + b.get_width()/2, v + 0.25, f"{v:.2f}", ha="center",
            color=DARK, fontsize=13, fontweight="bold")
ax.set_ylabel("full-matvec throughput\n[TDOF/s, aggregate]"); ax.set_ylim(0, 13.5)
ax.set_title("A trillion DOF, achieved: the speed-vs-scale pair\n(both MEASURED, 2048 GH200, A$\\cdot$1=1.8e-19)", fontsize=12)
# PA = fast path, memory-capped at ~1T (annotation in the clear band right of the PA bar)
ax.annotate("FAST path\n$d_G$ fills HBM (55 GB)\n$\\to$ memory-capped $\\sim$1T", xy=(0.27, 9.0),
            xytext=(0.42, 10.6), ha="center", va="center", color=BLUE, fontsize=9.5,
            fontweight="bold", arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.3))
# MF = scalable path (annotation in clear space above the shorter MF bar)
ax.annotate("SCALES to $\\sim$6T\nrecompute (26.8 B/DOF)\n$\\sim$3$\\times$ slower / matvec", xy=(1, 3.55),
            xytext=(1, 8.0), ha="center", va="bottom", color=GREEN, fontsize=9.5,
            fontweight="bold", arrowprops=dict(arrowstyle="->", color=GREEN, lw=1.3))
save(fig, "fig_trillion")

# --- 7. distributed HIGH-ORDER operator: MEASURED weak scaling, overlap ON (p=4, GH200) ---
# 1/8/16/32/64 nodes = 4/32/64/128/256 GPUs measured.
# Per-GPU GDOF/s weak scaling of the register-shuffle apply (p=4, store-d_G, ~491M DOF/GPU). Two
# MEASURED tracks (the no-overlap estimate was dropped -- it was derived from the noisy/non-monotonic
# small-scale reverseAdd, not a real run):
#   BLUE  apply-only (operator) -- flat ~10.3-10.6.
#   ORANGE full matvec WITH overlap -- hugs blue within ~4% (comm hidden behind compute).
# MEASURED 1/2/4/8/16/32/64 nodes (4-256 GPU), A.1 bit-exact (256 GPU: 3.5e-19); 512-GPU (128-node) pending.
# x is in GPUs (4 GPU/node): 1/2/4/8/16/32/64 nodes = 4/8/16/32/64/128/256 GPU.
ho_gpus    = [4, 8, 16, 32, 64, 128, 256, 512]
ho_apply   = [10.61, 10.54, 10.52, 10.42, 10.41, 10.50, 10.39, 10.37]   # apply-only, GDOF/s/GPU (256: 3770151, 512: 3770912)
ho_overlap = [10.28, 10.17, 10.14, 10.00, 10.02, 10.09, 9.96, 9.91]     # full matvec, overlap ON
fig, ax = plt.subplots(figsize=(7.0, 4.3), constrained_layout=True)
xpos = list(range(len(ho_gpus)))
ax.axhline(ho_apply[0], ls="--", lw=1.8, color="#8A8A8A", zorder=1)   # ideal: per-GPU rate held = 100%
ax.text(0.02, ho_apply[0] + 0.012, "ideal: per-GPU rate held (100%)", color="#8A8A8A", fontsize=9, va="bottom")
ax.plot(xpos, ho_apply,   "-o", color=BLUE,   lw=3.0, ms=12, zorder=5, label="apply-only (operator)")
ax.plot(xpos, ho_overlap, "-s", color=ORANGE, lw=3.0, ms=11, zorder=4, label="full matvec, overlap ON")
ax.set_xticks(xpos)
ax.set_xticklabels([str(g) for g in ho_gpus])
ax.set_xlim(-0.3, len(xpos) - 0.7)
ax.set_ylim(9.5, 10.85)
ax.set_xlabel("GPUs   ($p{=}4$, $\\sim$491M DOF/GPU held)")
ax.set_ylabel("throughput  [GDOF/s, per GPU]")
ax.set_title("Weak scaling: apply $\\sim$98%, full matvec $\\sim$93% efficiency to 512 GPU", fontsize=13, fontweight="bold")
ax.legend(frameon=False, loc="lower left", fontsize=11)
save(fig, "fig_ho_scaling")

# --- 8. storage per DOF vs p: matrix-free gets CHEAPER with order ("free to run via sum factorization") ---
# The smooth, monotonic counterpart to the bumpy achieved-throughput panel. Recompute (matrix-free,
# the headline) drops ~240x from p=1 to p=8 -> high order is CHEAPER to run matrix-free. store-d_G (PA)
# falls more slowly (still holds the gradient d_G), and the assembled 7-nnz CSR is a flat reference.
ps        = [1, 2, 3, 4, 5, 6, 7, 8]
mf_bpd    = [189.8, 23.7, 7.0, 3.0, 1.5, 0.9, 0.6, 0.4]            # recompute / matrix-free (headline)
pa_bpd    = [284.7, 160.1, 126.5, 111.2, 102.5, 96.9, 93.0, 90.1]  # store-d_G / PA
asm_ref   = 84.0                                                   # assembled 7-nnz CSR, flat across p
fig, ax = plt.subplots(figsize=(6.4, 4.2), constrained_layout=True)
ax.axhline(asm_ref, color=ORANGE, ls="--", lw=1.4, zorder=2, label="assembled 7-nnz (ref.)")
ax.semilogy(ps, pa_bpd, "-s", color=BLUE,  lw=2.2, ms=6, zorder=3, label="store $d_G$ / PA")
ax.semilogy(ps, mf_bpd, "-^", color=GREEN, lw=3.2, ms=9, zorder=4, label="recompute / matrix-free")
ax.set_xlabel("polynomial order $p$"); ax.set_ylabel("bytes / DOF")
ax.set_title("Matrix-free storage drops $\\sim$240$\\times$ ($p{=}1\\to8$):\nhigh order is $cheaper$ to run via sum factorization", fontsize=12)
ax.set_xticks(ps); ax.set_xlim(0.7, 8.3); ax.set_ylim(0.3, 400)
ax.annotate("recompute: $\\sim$240$\\times$ leaner\nat $p{=}8$ (smooth, monotonic)", xy=(8, 0.4),
            xytext=(4.6, 2.0), color=GREEN, fontsize=10, ha="left", fontweight="bold",
            arrowprops=dict(arrowstyle="->", color=GREEN, lw=1.2))
ax.legend(frameon=False, loc="center right", fontsize=10)
save(fig, "fig_storage_per_dof")

print("wrote fig_memory, fig_throughput, fig_convergence, fig_weakscale, fig_comm_pergpu, fig_trillion, fig_ho_scaling, fig_storage_per_dof (.png + .pdf)")
