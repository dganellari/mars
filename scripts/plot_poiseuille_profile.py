#!/usr/bin/env python3
"""Plot solved Poiseuille points on the exact plane-Poiseuille parabola.

The exact curve is the Wikipedia "Plane Poiseuille flow" solution
(https://en.wikipedia.org/wiki/Hagen-Poiseuille_equation):

    u(y) = G/(2 mu) * y * (h - y),   G = -dp/dx = 8 mu U_max / h^2

drawn from the closed-form formula (NOT fitted to the data). The dots are the
solved velocity at one fixed x station, read from the CSV written by
mars_poiseuille_flow (<prefix>_profile.csv: y, u_solved, u_analytic).

Usage:
    python3 scripts/plot_poiseuille_profile.py poiseuille_profile.csv [out.png]
        [--umax 1.5] [--mu 0.01]

Defaults match the validation case: rho=1, mu=0.01, U_mean=1 -> U_max=1.5,
h taken from the wall extents in the CSV. The script cross-checks its formula
against the driver's analytic column and warns if they disagree.
"""

import argparse
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv_path")
    ap.add_argument("out_path", nargs="?", default=None)
    ap.add_argument("--umax", type=float, default=1.5,
                    help="exact centerline speed U_max = 1.5*U_mean (default 1.5)")
    ap.add_argument("--mu", type=float, default=0.01,
                    help="dynamic viscosity, only for the G annotation (default 0.01)")
    args = ap.parse_args()
    out_path = args.out_path or args.csv_path.replace(".csv", ".png")

    y, u_solved, u_driver = [], [], []
    with open(args.csv_path) as f:
        for row in csv.DictReader(f):
            y.append(float(row["y"]))
            u_solved.append(float(row["u_solved"]))
            u_driver.append(float(row["u_analytic"]))
    y = np.array(y)
    u_solved = np.array(u_solved)
    u_driver = np.array(u_driver)

    # Walls from the data extents (the plane includes the wall nodes, u=0).
    y0, y1 = y.min(), y.max()
    h = y1 - y0
    G = 8.0 * args.mu * args.umax / h**2

    # Exact Wikipedia form, in wall coordinates yw = y - y0 in [0, h]:
    # u = G/(2 mu) * yw * (h - yw). With G as above this equals
    # U_max * 4 yw (h - yw) / h^2.
    def exact(yv):
        yw = yv - y0
        return G / (2.0 * args.mu) * yw * (h - yw)

    # Consistency: the driver's analytic column must be the same curve.
    # Tolerance 1e-4: the driver's wall extents come from a slab around the
    # probe plane and differ from this plane's by a few 1e-6 -- harmless.
    drift = float(np.max(np.abs(exact(y) - u_driver)))
    if drift > 1e-4 * args.umax:
        print(f"WARNING: driver analytic column differs from the closed form "
              f"by up to {drift:.3e} -- check U_max / wall extents")

    rms = float(np.sqrt(np.mean((u_solved - exact(y)) ** 2)))

    yy = np.linspace(y0, y1, 400)
    fig, ax = plt.subplots(figsize=(5, 4.2))
    ax.plot(exact(yy), yy, "-", color="tab:blue", lw=2,
            label=f"exact: u=G/(2μ)·y(h−y), G={G:.3g}")
    ax.plot(u_solved, y, "o", color="tab:red", ms=3.5, mfc="none",
            label="MARS (solved)")
    ax.set_xlabel("u (streamwise velocity)")
    ax.set_ylabel("y (wall-normal)")
    ax.set_title(f"Plane Poiseuille profile, RMS = {rms:.2e}")
    ax.legend(loc="center left", fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"wrote {out_path}")
    print(f"  h={h:.4g}  U_max={args.umax}  mu={args.mu}  G=8*mu*U_max/h^2={G:.4g}")
    print(f"  RMS vs exact: {rms:.3e}")


if __name__ == "__main__":
    main()
