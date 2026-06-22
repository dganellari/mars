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
    for ext in ("png", "pdf"):
        fig.savefig(f"{name}.{ext}")
    plt.close(fig)

p = [1, 2, 3, 4, 5, 6, 7]

# --- 1. memory crossover: bytes/DOF vs p ---
asm = [324, 1500, 4116, 8748, 15972, 26364, 40500]
mf  = [421, 221, 169, 147, 134, 126, 121]
fig, ax = plt.subplots(figsize=(6.2, 4.0))
ax.semilogy(p, asm, "-o", color=ORANGE, lw=2.4, ms=7, label="assembled (stored matrix)")
ax.semilogy(p, mf,  "-s", color=BLUE,   lw=2.4, ms=7, label="matrix-free")
ax.set_xlabel("polynomial order $p$"); ax.set_ylabel("bytes / DOF")
ax.set_title("Operator memory footprint")
ax.set_xticks(p); ax.set_xlim(0.7, 7.3)
ax.annotate("336$\\times$ at $p{=}7$", xy=(7, 40500), xytext=(4.6, 30000),
            color=ORANGE, fontsize=13, fontweight="bold",
            arrowprops=dict(arrowstyle="->", color=ORANGE, lw=1.4))
ax.legend(frameon=False, loc="center right")
save(fig, "fig_memory")

# --- 2. throughput vs p (flat) ---
mdofs = [4794, 7984, 7895, 5309, 4258, 4157, 4530]
fig, ax = plt.subplots(figsize=(6.2, 4.0))
ax.bar(p, [m/1000 for m in mdofs], color=BLUE, width=0.6, zorder=3)
ax.axhspan(4.3, 8.0, color=BLUE, alpha=0.08, zorder=0)
ax.set_xlabel("polynomial order $p$"); ax.set_ylabel("throughput  [GDOF/s]")
ax.set_title("Matrix-free apply throughput (H100)")
ax.set_xticks(p); ax.set_ylim(0, 9)
ax.text(4.0, 8.4, "flat $\\approx$4--8 GDOF/s  ($p{=}7\\approx p{=}1$)",
        ha="center", color=DARK, fontsize=12)
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
apply_eff  = [100.0, 96.9, 95.0]   # operator apply (compute) -- near flat
overlap_eff= [100.0, 56.8, 40.3]   # full matvec, forward halo overlapped
block_eff  = [100.0, 50.1, 33.0]   # full matvec, blocking halo
fig, ax = plt.subplots(figsize=(6.2, 4.0))
ax.plot(gpus, apply_eff,   "-o", color=BLUE,   lw=2.6, ms=9, label="operator apply (compute)")
ax.plot(gpus, overlap_eff, "-^", color=GREEN,  lw=2.6, ms=9, label="full matvec, comm overlap")
ax.plot(gpus, block_eff,   "-s", color=ORANGE, lw=2.6, ms=9, label="full matvec, blocking halo")
ax.axhline(100, color="#999999", ls=":", lw=1.2)
ax.set_xscale("log", base=2); ax.set_xticks(gpus); ax.set_xticklabels(["1", "8", "32"])
ax.set_xlabel("GPUs   (2M elements / GPU, held fixed:  2M $\\to$ 64M)")
ax.set_ylabel("parallel efficiency  [%]")
ax.set_title("Weak scaling on GH200")
ax.set_ylim(0, 112)
ax.annotate("95% @ 64M DOF, 32 GPU", xy=(32, 95.0), xytext=(2.3, 85),
            color=BLUE, fontsize=11, fontweight="bold")
ax.annotate("+overlap 1.22$\\times$", xy=(32, 40.3), xytext=(11, 62),
            color=GREEN, fontsize=11, fontweight="bold",
            arrowprops=dict(arrowstyle="->", color=GREEN, lw=1.3))
ax.legend(frameon=False, loc="lower left")
save(fig, "fig_weakscale")

# --- 5. comm fraction vs per-GPU size (8 GPUs fixed) -- the trillion self-resolve ---
pergpu_dof = [2.05e6, 8.06e6, 16.10e6, 32.16e6, 64.39e6]   # measured per-GPU DOF (2/8/16/32/64M elems)
comm_meas  = [51.9, 41.9, 33.5, 28.6, 19.2]               # measured comm % (64M = post halo-fix, gates 1e-18)
# fit comm/apply = k * V^p to the four points (log-log): p~-0.36, k~200
kfit, pfit = 200.0, -0.36
def cfrac(v):
    r = kfit * v ** pfit
    return 100.0 * r / (1.0 + r)
Vgrid = [10 ** (6 + 0.05 * i) for i in range(0, 61)]   # 1e6 .. 1e9
fig, ax = plt.subplots(figsize=(6.2, 4.0))
ax.plot(Vgrid, [cfrac(v) for v in Vgrid], "-", color=BLUE, lw=2.2,
        label="model  comm/apply $\\propto V^{-1/3}$")
ax.plot(pergpu_dof, comm_meas, "o", color=ORANGE, ms=11, zorder=5, label="measured (8 GPUs)")
ax.axvline(3.0e8, color=GREEN, ls="--", lw=1.6)
ax.annotate("trillion-relevant\n$\\sim$300M DOF/GPU\n$\\to\\sim$17%% comm" % (),
            xy=(3.0e8, cfrac(3.0e8)), xytext=(2.2e7, 58),
            color=GREEN, fontsize=11, ha="center",
            arrowprops=dict(arrowstyle="->", color=GREEN, lw=1.3))
ax.set_xscale("log"); ax.set_xlabel("DOF per GPU"); ax.set_ylabel("communication fraction  [%]")
ax.set_title("Comm self-resolves with per-GPU size")
ax.set_ylim(0, 70); ax.legend(frameon=False, loc="lower left")
save(fig, "fig_comm_pergpu")

# --- 6. path to a trillion: per-GPU capacity by mode vs the 1e12-on-Alps line ---
modes  = ["assembly\n(CSR)", "matrix-free\napply", "domain\ndecomp (DD)"]
ceil_M = [70, 77, 108]            # validated per-GPU ceilings (M elements/GPU)
alps_gpus = 10752                 # Alps GH200
need_M = 1e12 / alps_gpus / 1e6   # per-GPU needed for 1e12 on full Alps (~93M)
fig, ax = plt.subplots(figsize=(6.2, 4.0))
bars = ax.bar(modes, ceil_M, color=[ORANGE, ORANGE, GREEN], width=0.62, zorder=3)
ax.axhline(need_M, color=DARK, ls="--", lw=1.6, zorder=2)
ax.text(2.45, need_M + 2, "10$^{12}$ on full Alps ($\\sim$%.0fM/GPU)" % need_M,
        ha="right", va="bottom", color=DARK, fontsize=11)
for b, v in zip(bars, ceil_M):
    ax.text(b.get_x() + b.get_width()/2, v + 1.5, f"{v}M", ha="center",
            color=DARK, fontsize=12, fontweight="bold")
ax.set_ylabel("per-GPU capacity  [M elements]"); ax.set_ylim(0, 125)
ax.set_title("Path to a trillion (Alps, $\\sim$10.7k GH200)")
ax.text(0.03, 0.97, "DD reaches 10$^{12}$ at $\\sim$9,300 GPUs\n(validated 108.7M/GPU clean, gates 1e-18)",
        transform=ax.transAxes, va="top", fontsize=10, color=GREEN)
save(fig, "fig_trillion")

# --- 7. distributed HIGH-ORDER operator: weak scaling + comm self-resolve (p=4, GH200) ---
ho_gpus      = [1, 8]
ho_apply_eff = [100.0, 98.1]   # apply-only, ~4.2M DOF/GPU held (5965 -> 46811 MDOF/s)
ho_full_eff  = [100.0, 78.4]   # full matvec, blocking halo (5790 -> 36319 MDOF/s)
fig, ax = plt.subplots(figsize=(6.2, 4.0))
ax.plot(ho_gpus, ho_apply_eff, "-o", color=BLUE,   lw=2.6, ms=10, label="operator apply (compute)")
ax.plot(ho_gpus, ho_full_eff,  "-s", color=ORANGE, lw=2.6, ms=10, label="full matvec, blocking halo")
ax.axhline(100, color="#999999", ls=":", lw=1.2)
ax.set_xscale("log", base=2); ax.set_xticks(ho_gpus); ax.set_xticklabels(["1", "8"])
ax.set_xlabel("GPUs   (p=4, ~4.2M DOF/GPU held;  device==host halo bit-exact 1e-18)")
ax.set_ylabel("parallel efficiency  [%]")
ax.set_title("High-order matrix-free: distributed weak scaling")
ax.set_ylim(0, 112)
ax.annotate("apply 98% @ 8 GPU\n(5.9 GDOF/s/GPU = single-GPU rate)", xy=(8, 98.1), xytext=(1.12, 58),
            color=BLUE, fontsize=11, fontweight="bold")
ax.text(1.08, 22,
        "comm self-resolves with per-GPU size:\nat 8 GPU,  4M$\\to$11M DOF/GPU\n$\\Rightarrow$ comm 22$\\to$15%,  full-matvec 78$\\to$87%",
        fontsize=10, color=DARK,
        bbox=dict(boxstyle="round", fc="#F4F4F4", ec="#CCCCCC"))
ax.legend(frameon=False, loc="lower left")
save(fig, "fig_ho_scaling")

print("wrote fig_memory, fig_throughput, fig_convergence, fig_weakscale, fig_comm_pergpu, fig_trillion, fig_ho_scaling (.png + .pdf)")
