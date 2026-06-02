#!/usr/bin/env python3
"""Plot per-step convergence diagnostics from a MARS pump NS run.

Usage:
    python plot_convergence.py run.log [out.png]

Reads the stdout captured from a `mars_pump` run (argv[1]) and parses the
per-step monitor lines, e.g.:

    Step    50  ft=0.12  u_rms=1.2e-01  u_max=3.4  uMax/U=0.80  d(u_rms)=5.0e-04 \\
        div*L/U=2.1e-03  cg_p=42  pres_r0=1.3e+02  pres_res=8.7e-07  cg_uvw=12

It writes a 4-panel PNG (argv[2], default convergence.png) so a CFD engineer
can see at a glance whether the run is marching toward steady state.

WHY iterations-trending-down means convergence: the pressure CG is warm-started
each step with phi from the previous step. As the flow nears steady state the
pressure field changes less from one step to the next, so that warm start is
already close to the new solution and the CG needs fewer iterations to drive the
residual under tolerance. A sustained downward trend in cg_p (the pressure CG
iteration count) is therefore a direct, cheap proxy for "the solution has
stopped changing" -- i.e. convergence to steady state. The companion panels
(non-dim divergence and the steady-state velocity residual d(u_rms)) should
trend down for the same reason.
"""

import re
import sys

import matplotlib

# Agg backend: render to file without a display so this works headless over
# SSH / on a login node where no X server is attached.
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# One tolerant regex per field. We search each token independently instead of
# matching the whole line so that optional fields (divRC*L/U, dt) and any
# early line that is missing a field do not break the parse.
_NUM = r"([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)"
_FIELDS = {
    "step":    re.compile(r"\bStep\s+([0-9]+)"),
    "ft":      re.compile(r"\bft=" + _NUM),
    "cg_p":    re.compile(r"\bcg_p=([0-9]+)"),
    "div":     re.compile(r"\bdiv\*L/U=" + _NUM),
    "d_urms":  re.compile(r"\bd\(u_rms\)=" + _NUM),
    "pres_r0": re.compile(r"\bpres_r0=" + _NUM),
    "pres_res": re.compile(r"\bpres_res=" + _NUM),
}


def parse_log(path):
    """Return a dict of field-name -> list of (step, value), skipping lines or
    fields that are absent. Only lines that carry a step index are considered."""
    series = {k: [] for k in _FIELDS if k != "step"}
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if "Step" not in line:
                continue
            m_step = _FIELDS["step"].search(line)
            if not m_step:
                continue
            step = int(m_step.group(1))
            for name, rx in _FIELDS.items():
                if name == "step":
                    continue
                m = rx.search(line)
                if not m:
                    continue  # robust to missing fields on early lines
                val = int(m.group(1)) if name == "cg_p" else float(m.group(1))
                series[name].append((step, val))
    return series


def _plot_panel(ax, data, title, ylabel):
    """Draw one log-y panel. Empty / all-non-positive series get a placeholder
    so the figure still renders cleanly."""
    ax.set_title(title)
    ax.set_xlabel("step")
    ax.set_ylabel(ylabel)
    ax.grid(True, which="both", ls=":", alpha=0.5)
    # log scale needs strictly positive y; drop non-positive samples.
    pts = [(s, v) for (s, v) in data if v > 0]
    if not pts:
        ax.text(0.5, 0.5, "no data", ha="center", va="center",
                transform=ax.transAxes)
        return
    steps = [s for s, _ in pts]
    vals = [v for _, v in pts]
    ax.set_yscale("log")
    ax.plot(steps, vals, marker=".", ms=3, lw=1)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    log_path = sys.argv[1]
    out_path = sys.argv[2] if len(sys.argv) > 2 else "convergence.png"

    series = parse_log(log_path)

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle(f"MARS pump convergence -- {log_path}", fontsize=13)

    _plot_panel(axes[0, 0], series["cg_p"],
                "Pressure CG iterations/step (down = converging)",
                "cg_p (iters)")
    _plot_panel(axes[0, 1], series["div"],
                "Non-dim divergence (max)",
                "div*L/U")
    _plot_panel(axes[1, 0], series["d_urms"],
                "Steady-state residual d(u_rms) (down = converging)",
                "d(u_rms)")
    _plot_panel(axes[1, 1], series["pres_r0"],
                "Pressure RHS magnitude |r0|/step",
                "pres_r0 (|r0|)")

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_path, dpi=120)
    n = len(series["cg_p"])
    print(f"parsed {n} step lines from {log_path} -> wrote {out_path}")


if __name__ == "__main__":
    main()
