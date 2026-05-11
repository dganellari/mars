#!/usr/bin/env python3
"""Generate a tet cube mesh in MARS binary format (i0..i3.int64 + x,y,z.float32).

Each axis-aligned unit hex is split into 6 tets via Kuhn's decomposition
(all 6 share the main diagonal v0-v6). All tets are right-handed (positive
Jacobian); node ordering is fixed up per-tet by swapping nodes 2 and 3 when
the signed volume would be negative.

The vertex grid is identical to scripts/generate_hex_cube.py — same coord
files, just different connectivity files (4 columns instead of 8, 6x as
many elements).
"""

import argparse
import numpy as np
import os
import time


def generate_tet_cube(nx, ny, nz, output_dir, use_int64=True):
    num_nodes = (nx + 1) * (ny + 1) * (nz + 1)
    num_hex   = nx * ny * nz
    num_elems = 6 * num_hex

    print(f"Generating {nx}x{ny}x{nz} tet cube (Kuhn 6-tet split): "
          f"{num_elems:,} tets from {num_hex:,} hexes, {num_nodes:,} nodes")
    os.makedirs(output_dir, exist_ok=True)

    dtype = np.int64 if use_int64 else np.int32
    ext   = "int64" if use_int64 else "int32"

    # -- coords --
    t0 = time.time()
    print("Writing coordinates...")
    inv_x = np.float32(1.0 / nx)
    inv_y = np.float32(1.0 / ny)
    inv_z = np.float32(1.0 / nz)
    ix = np.arange(nx + 1, dtype=np.float32) * inv_x
    iy = np.arange(ny + 1, dtype=np.float32) * inv_y
    iz = np.arange(nz + 1, dtype=np.float32) * inv_z

    xs = np.broadcast_to(ix, (nz + 1, ny + 1, nx + 1)).reshape(-1)
    ys = np.broadcast_to(iy[:, None], (ny + 1, nx + 1)).reshape(-1)
    ys = np.tile(ys, nz + 1)
    zs = np.repeat(iz, (ny + 1) * (nx + 1))

    xs.astype(np.float32, copy=False).tofile(os.path.join(output_dir, "x.float32"))
    ys.astype(np.float32, copy=False).tofile(os.path.join(output_dir, "y.float32"))
    zs.astype(np.float32, copy=False).tofile(os.path.join(output_dir, "z.float32"))
    del xs, ys, zs
    print(f"  coords: {time.time() - t0:.2f}s")

    # -- connectivity --
    print("Writing connectivity (4 columns x 6 tets per hex)...")
    t0 = time.time()
    _write_connectivity(nx, ny, nz, output_dir, dtype, ext)
    print(f"  connectivity: {time.time() - t0:.2f}s")

    print(f"Done. Output in {output_dir}")
    print(f"  Elements: {num_elems:,}")
    print(f"  Nodes:    {num_nodes:,}")
    print(f"  Format:   {ext} + float32")


# Kuhn 6-tet decomposition of a unit hex with corner numbering:
#   v0=(0,0,0) v1=(1,0,0) v2=(1,1,0) v3=(0,1,0)
#   v4=(0,0,1) v5=(1,0,1) v6=(1,1,1) v7=(0,1,1)
# Each tet has positive signed volume with this ordering (verified below).
KUHN_TETS = np.array([
    [0, 1, 2, 6],
    [0, 2, 3, 6],
    [0, 3, 7, 6],
    [0, 7, 4, 6],
    [0, 4, 5, 6],
    [0, 5, 1, 6],
], dtype=np.int64)


def _verify_kuhn_orientation():
    """Sanity check: every Kuhn tet has positive signed volume."""
    corners = np.array([
        [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
        [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],
    ], dtype=np.float64)
    for t in KUHN_TETS:
        p = corners[t]
        v1 = p[1] - p[0]; v2 = p[2] - p[0]; v3 = p[3] - p[0]
        det = np.dot(np.cross(v1, v2), v3)
        if det <= 0:
            raise RuntimeError(f"Kuhn tet {t} has non-positive det: {det}")


def _write_connectivity(nx, ny, nz, output_dir, dtype, ext):
    _verify_kuhn_orientation()

    nx1 = nx + 1
    ny1 = ny + 1
    plane = ny1 * nx1

    # Hex corner offsets in the global node grid (i,j,k) -> linear index
    # consistent with generate_hex_cube.py:
    #   n0 = k*plane + j*nx1 + i, then standard hex ordering
    hex_off = np.array([
        0,                    # v0 = (i,   j,   k)
        1,                    # v1 = (i+1, j,   k)
        nx1 + 1,              # v2 = (i+1, j+1, k)
        nx1,                  # v3 = (i,   j+1, k)
        plane,                # v4 = (i,   j,   k+1)
        plane + 1,            # v5 = (i+1, j,   k+1)
        plane + nx1 + 1,      # v6 = (i+1, j+1, k+1)
        plane + nx1,          # v7 = (i,   j+1, k+1)
    ], dtype=np.int64)

    # Hex base-node indices
    i = np.arange(nx, dtype=np.int64)
    j = np.arange(ny, dtype=np.int64)
    k = np.arange(nz, dtype=np.int64)
    K, J, I = np.meshgrid(k, j, i, indexing='ij')
    n0 = (K * plane + J * nx1 + I).reshape(-1)   # shape (num_hex,)
    del K, J, I

    num_hex = n0.size
    # For each of 6 tets, build its 4 connectivity columns
    # Resulting layout: tets ordered such that the 6 tets of hex h are
    # consecutive elements [6*h, 6*h+1, ..., 6*h+5]
    #
    # Memory: 6 * num_hex * 4 * 8 bytes; for 1024^3 = ~48 GB — chunked
    # version may be needed there but for cube16/64/128 it's fine.
    cols = [np.empty(6 * num_hex, dtype=np.int64) for _ in range(4)]
    for ti in range(6):
        tet = KUHN_TETS[ti]
        for c in range(4):
            cols[c][ti::6] = n0 + hex_off[tet[c]]

    from concurrent.futures import ThreadPoolExecutor
    def _write_one(c, col):
        col.astype(dtype, copy=False).tofile(os.path.join(output_dir, f"i{c}.{ext}"))
    with ThreadPoolExecutor(max_workers=4) as ex:
        list(ex.map(_write_one, range(4), cols))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Generate tet cube mesh in MARS binary format")
    p.add_argument("--nx", type=int, default=16)
    p.add_argument("--ny", type=int, default=16)
    p.add_argument("--nz", type=int, default=16)
    p.add_argument("--output", type=str, default="tet_cube.mesh")
    p.add_argument("--int32", action="store_true", help="Use int32 instead of int64")
    args = p.parse_args()

    generate_tet_cube(args.nx, args.ny, args.nz, args.output, use_int64=not args.int32)
