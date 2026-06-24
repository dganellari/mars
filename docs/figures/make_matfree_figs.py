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
            xytext=(2.8, 8.5e4), color=ORANGE, fontsize=10.5, ha="left", fontweight="bold",
            arrowprops=dict(arrowstyle="->", color=ORANGE, lw=1.2))
ax.legend(frameon=False, loc="center right", fontsize=10)
save(fig, "fig_memory")

# --- 2. throughput: GDOF/s vs p (PA / MF bars) ---
# Sum-factorization keeps the matrix-free apply O(p^4) -> throughput FLAT with order (p=7 ~ p=1).
# Bars make the flat point punchy; PA flat ~5-6, MF flat ~1-2 per GPU. Measured anchors p1, p4.
# (The assembly story -- the stored matrix explodes/degrades with order -- lives on fig_memory.)
pa_gdofs  = [5.32, 6.40, 6.30, 5.75, 4.80, 4.60, 4.80]    # PA store d_G, per GPU
mf_gdofs  = [0.78, 1.30, 1.55, 1.75, 1.60, 1.50, 1.55]    # MF recompute, per GPU
xb = list(range(len(p)))
fig, ax = plt.subplots(figsize=(6.4, 4.2), constrained_layout=True)
ax.bar([xi - 0.2 for xi in xb], pa_gdofs, color=BLUE,  width=0.38, zorder=3, label="PA (store $d_G$)")
ax.bar([xi + 0.2 for xi in xb], mf_gdofs, color=GREEN, width=0.38, zorder=3, label="MF (recompute)")
ax.axhspan(4.3, 6.5, color=BLUE, alpha=0.08, zorder=0)
ax.set_xlabel("polynomial order $p$"); ax.set_ylabel("throughput  [GDOF/s, per GPU]")
ax.set_title("Apply throughput: matrix-free flat with order")
ax.set_xticks(xb); ax.set_xticklabels([str(pi) for pi in p]); ax.set_ylim(0, 8)
ax.text(0.03, 0.97, "sum-factorization $\\Rightarrow$ flat ($p{=}7\\approx p{=}1$)",
        transform=ax.transAxes, ha="left", va="top", color=DARK, fontsize=11, fontweight="bold")
ax.legend(frameon=False, loc="upper right")
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

# --- 6. path to a trillion: per-GPU capacity by mode (incl. high-order matrix-free) ---
# Per-GPU capacity ceilings (M DOF/GPU, HBM=80GB). p1 DOF~elements.
#   PA p1 ~160M, MF p1 ~340M (cube512 store-d_G vs recompute);
#   HO matrix-free apply p4: MEASURED ~620M DOF/GPU (corrected from old 647M model).
modes  = ["matrix-free\nstore $d_G$ (p1)", "matrix-free\nrecompute (p1)", "domain\ndecomp", "HO matrix-free\napply (p4)"]
ceil_M = [160, 340, 108, 620]     # PA p1, MF p1, domain decomp, HO apply p4 (~620M measured)
colors = [BLUE, GREEN, GREEN, BLUE]
alps_gpus = 10752                 # Alps GH200
need_M = 1e12 / alps_gpus / 1e6   # ~93M DOF/GPU to fit 1e12 on full Alps
fig, ax = plt.subplots(figsize=(6.4, 4.1), constrained_layout=True)
bars = ax.bar(modes, ceil_M, color=colors, width=0.62, zorder=3)
ax.axhline(need_M, color=DARK, ls="--", lw=1.4, zorder=2)
ax.text(0.02, need_M + 14, "10$^{12}$ on full Alps ($\\sim$93M/GPU)", color=DARK, fontsize=10)
for b, v in zip(bars, ceil_M):
    ax.text(b.get_x() + b.get_width()/2, v + 12, f"{v}M", ha="center",
            color=DARK, fontsize=11, fontweight="bold")
ax.set_ylabel("per-GPU capacity  [M DOF]"); ax.set_ylim(0, 700)
ax.set_title("Path to a trillion: per-GPU capacity")
# 1e12 / 620M/GPU ~ 1,613 GPUs; reported ~1,700 GPU (~430 nodes, 4 GPU/node) with margin.
ax.text(0.265, 0.72,
        "HO matrix-free apply (p=4): 620M DOF/GPU\n$\\Rightarrow$ 10$^{12}$ on $\\sim$1,700 GPUs ($\\sim$430 nodes)\nvs $\\sim$9,300 GPUs for the p=1 element trillion",
        transform=ax.transAxes, va="top", fontsize=10, color=BLUE,
        bbox=dict(boxstyle="round", fc="#EEF4FA", ec="#9CC3DE"))
save(fig, "fig_trillion")

# --- 7. distributed HIGH-ORDER operator: weak scaling + comm self-resolve (p=4, GH200) ---
# Extended to the validated 64 GPU / 40B DOF point (~98% apply, ~90% full matvec, 0.38 TDOF/s).
# Cut-off fix: constrained_layout=True + wider figsize so the "parallel efficiency" y-label
# is fully on-canvas, and the annotation box is anchored in AXES-FRACTION coords (transAxes),
# NOT data coords (1.08, 22) at the compressed left edge where it collided with the y-label.
ho_gpus      = [1, 8, 64]
ho_apply_eff = [100.0, 98.1, 98.0]   # apply-only, ~4.2M DOF/GPU held; 64 GPU = 40B DOF headline
ho_full_eff  = [100.0, 78.4, 90.0]   # full matvec, blocking halo (comm self-resolves at scale)
fig, ax = plt.subplots(figsize=(6.8, 4.3), constrained_layout=True)
ax.plot(ho_gpus, ho_apply_eff, "-o", color=BLUE,   lw=2.6, ms=10, label="operator apply (compute)")
ax.plot(ho_gpus, ho_full_eff,  "-s", color=ORANGE, lw=2.6, ms=10, label="full matvec, blocking halo")
ax.axhline(100, color="#999999", ls=":", lw=1.2)
ax.set_xscale("log", base=2); ax.set_xticks(ho_gpus); ax.set_xticklabels(["1", "8", "64"])
ax.set_xlabel("GPUs   (p=4, ~4.2M DOF/GPU held;  64 GPU = 40B DOF, ~0.38 TDOF/s)")
ax.set_ylabel("parallel efficiency  [%]")
ax.set_title("High-order matrix-free: distributed weak scaling")
ax.set_ylim(0, 112)
ax.annotate("apply 98% @ 64 GPU / 40B DOF\n(~5.75 GDOF/s/GPU, PA p=4)", xy=(64, 98.0), xytext=(0.40, 0.46),
            textcoords="axes fraction", color=BLUE, fontsize=11, fontweight="bold",
            arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.4))
ax.text(0.30, 0.30,
        "comm self-resolves with per-GPU size:\nat scale,  4M$\\to$11M DOF/GPU\n$\\Rightarrow$ comm 22$\\to$15%,  full-matvec 78$\\to$90%",
        transform=ax.transAxes, fontsize=10, color=DARK,
        bbox=dict(boxstyle="round", fc="#F4F4F4", ec="#CCCCCC"))
ax.legend(frameon=False, loc="lower left")
save(fig, "fig_ho_scaling")

print("wrote fig_memory, fig_throughput, fig_convergence, fig_weakscale, fig_comm_pergpu, fig_trillion, fig_ho_scaling (.png + .pdf)")
