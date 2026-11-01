#!/usr/bin/env python3
"""Host study on blocks of a real Hex8 Exodus mesh (the NASA wing tip): the smoothers of
scripts/bp3_smoother_study.py on the elements within --layers vertex-neighbour layers of a seed
element, Dirichlet on the block's boundary, b = A u with u = sin(pi x) sin(pi y) sin(pi z) on the
block's bounding box mapped to the unit cube. The coarse solve is exact, so this shows whether the
mesh's own elements make the point smoother need many iterations, and what the other smoothers do.
Run: python3 scripts/bp3_wing_block_study.py --mesh large_hex_only.e --seeds e0,e1 [--layers 6]
     [--smoothers jacobi,element,nodeline]   (seeds: zero-based element indices in file order)
"""
import argparse
import time

import netCDF4 as nc
import numpy as np

import bp3_smoother_study as S

FRAMES = np.array([[0, 1, 3, 4], [1, 2, 0, 5], [2, 3, 1, 6], [3, 0, 2, 7],
                   [4, 7, 5, 0], [5, 4, 6, 1], [6, 5, 7, 2], [7, 6, 4, 3]])
EDGES = np.array([[[0, 1], [3, 2], [4, 5], [7, 6]], [[0, 3], [1, 2], [4, 7], [5, 6]],
                  [[0, 4], [1, 5], [2, 6], [3, 7]]])


def shapes(p):   # edge aspect ratio and minimum scaled Jacobian, as scripts/hex_mesh_quality.py
    vec = p[:, EDGES[:, :, 1]] - p[:, EDGES[:, :, 0]]
    lengths = np.linalg.norm(vec, axis=3).mean(axis=2)
    a, b, c = (p[:, FRAMES[:, i]] - p[:, FRAMES[:, 0]] for i in (1, 2, 3))
    sj = np.einsum("nki,nki->nk", np.cross(a, b), c) / (
        np.linalg.norm(a, axis=2) * np.linalg.norm(b, axis=2) * np.linalg.norm(c, axis=2))
    return lengths.max(1) / lengths.min(1), sj.min(1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", required=True)
    ap.add_argument("--seeds", required=True)
    ap.add_argument("--layers", type=int, default=6)
    ap.add_argument("--p", type=int, default=3)
    ap.add_argument("--smoothers", default="jacobi,element,nodeline")
    ap.add_argument("--pre", type=int, default=2)
    ap.add_argument("--post", type=int, default=2)
    ap.add_argument("--rtol", type=float, default=1e-10)
    ap.add_argument("--maxit", type=int, default=3000)
    args = ap.parse_args()

    t0 = time.time()
    f = nc.Dataset(args.mesh)
    X = np.stack([f.variables[c][:].data for c in ("coordx", "coordy", "coordz")], 1)
    conn = np.concatenate([f.variables[f"connect{b + 1}"][:].data for b in range(len(f.dimensions["num_el_blk"]))]) - 1
    flat = conn.ravel(); order = np.argsort(flat, kind="stable")
    start = np.searchsorted(flat[order], np.arange(len(X) + 1))
    print(f"{args.mesh}: {len(conn)} Hex8, {len(X)} nodes, read and indexed in {time.time() - t0:.0f} s; "
          f"p = {args.p} -> 1, Chebyshev ({args.pre},{args.post}), exact p = 1 solve, FCG to true residual "
          f"{args.rtol:.0e}", flush=True)
    for seed in (int(s) for s in args.seeds.split(",")):
        elems = {seed}; frontier = [seed]
        for _ in range(args.layers):
            nodes = np.unique(conn[frontier])
            near = set((np.concatenate([order[start[v]:start[v + 1]] for v in nodes]) // 8).tolist())
            frontier = list(near - elems); elems |= near
        B = np.array(sorted(elems))
        corners, vid = X[conn[B]], conn[B]
        ar, sj = shapes(corners)
        fl, cl = S.levels(corners, vid, args.p)
        lo, hi = fl.xn.reshape(-1, 3).min(0), fl.xn.reshape(-1, 3).max(0)
        u = fl.mask * np.prod(np.sin(np.pi * (fl.xn - lo) / (hi - lo)), -1)
        ar0, sj0 = shapes(X[conn[[seed]]])
        print(f"  seed {seed} (aspect ratio {ar0[0]:.0f}, scaled Jacobian {sj0[0]:.3f}): {len(B)} elements, "
              f"{len(fl.inner)} free nodes; aspect ratio median {np.median(ar):.0f} max {ar.max():.0f}, scaled "
              f"Jacobian min {sj.min():.3f}", flush=True)
        for name in args.smoothers.split(","):
            t1 = time.time()
            prec, label = S.smoother(name, fl, corners, vid)
            it, err, lmax = S.solve(fl, cl, prec, u, args.pre, args.post, args.rtol, args.maxit)
            print(f"    {label}: {it} it (err {err:.1e}, lmax {lmax:.2f}, {time.time() - t1:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
