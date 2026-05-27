#!/usr/bin/env python3
"""Generate a hex cube mesh in MARS binary format (i*.int64 + x.float32).

Vectorized numpy version: all coordinate and connectivity arrays are built
with numpy ops (no Python loops over elements/nodes). For a 1024^3 cube
(~1.07B elements) this avoids ~hours of interpreter overhead and is
limited mainly by disk write bandwidth.

Memory note: at 1024^3 each int64 connectivity column is 8 GB, so 8 columns
peak at 64 GB resident. If that exceeds RAM, use --chunked which streams
connectivity in z-slabs and uses only a few GB at a time.
"""

import argparse
import numpy as np
import os
import time


def generate_hex_cube(nx, ny, nz, output_dir, use_int64=True, chunked=False):
    num_nodes = (nx + 1) * (ny + 1) * (nz + 1)
    num_elems = nx * ny * nz

    print(f"Generating {nx}x{ny}x{nz} hex cube: {num_elems:,} elements, {num_nodes:,} nodes")
    os.makedirs(output_dir, exist_ok=True)

    dtype = np.int64 if use_int64 else np.int32
    ext   = "int64" if use_int64 else "int32"

    # ----- Coordinates: build with broadcasted arange + tofile -----
    t0 = time.time()
    print("Writing coordinates...")
    # i runs fastest; ix has shape (nz+1, ny+1, nx+1) with i values
    inv_x = np.float32(1.0 / nx)
    inv_y = np.float32(1.0 / ny)
    inv_z = np.float32(1.0 / nz)

    # Use broadcasting; full arrays of shape (nz+1, ny+1, nx+1) flattened to 1-D
    ix = np.arange(nx + 1, dtype=np.float32) * inv_x
    iy = np.arange(ny + 1, dtype=np.float32) * inv_y
    iz = np.arange(nz + 1, dtype=np.float32) * inv_z

    # x varies fastest, then y, then z. Use np.broadcast_to + tile patterns.
    # x: tile ix by (nz+1)*(ny+1)
    xs = np.broadcast_to(ix, (nz + 1, ny + 1, nx + 1)).reshape(-1)
    # y: each y repeated (nx+1) times, then sequence repeated (nz+1) times
    ys = np.broadcast_to(iy[:, None], (ny + 1, nx + 1)).reshape(-1)
    ys = np.tile(ys, nz + 1)
    # z: each z repeated (ny+1)*(nx+1) times
    zs = np.repeat(iz, (ny + 1) * (nx + 1))

    xs.astype(np.float32, copy=False).tofile(os.path.join(output_dir, "x.float32"))
    ys.astype(np.float32, copy=False).tofile(os.path.join(output_dir, "y.float32"))
    zs.astype(np.float32, copy=False).tofile(os.path.join(output_dir, "z.float32"))
    del xs, ys, zs
    t1 = time.time()
    print(f"  coords: {t1 - t0:.2f}s")

    # ----- Connectivity: 8 columns, fully vectorized -----
    print("Writing connectivity...")
    t0 = time.time()

    if chunked:
        _write_connectivity_chunked(nx, ny, nz, output_dir, dtype, ext)
    else:
        _write_connectivity_dense(nx, ny, nz, output_dir, dtype, ext)

    t1 = time.time()
    print(f"  connectivity: {t1 - t0:.2f}s")

    print(f"Done. Output in {output_dir}")
    print(f"  Elements: {num_elems:,}")
    print(f"  Nodes:    {num_nodes:,}")
    print(f"  Format:   {ext} + float32")


def _write_connectivity_dense(nx, ny, nz, output_dir, dtype, ext):
    """Build all 8 connectivity columns in memory and write each.

    Memory: 8 * num_elems * sizeof(dtype) bytes peak (8 * 8 * num_elems for int64).
    For 1024^3 that's 64 GB — only viable on big-RAM nodes.
    """
    nx1 = nx + 1
    ny1 = ny + 1

    # Element index space: (nz, ny, nx). Build i, j, k once.
    # n0 = k*(ny+1)*(nx+1) + j*(nx+1) + i
    # offsets fit in int64 always; cast at the end if needed
    i = np.arange(nx, dtype=np.int64)
    j = np.arange(ny, dtype=np.int64)
    k = np.arange(nz, dtype=np.int64)

    # Broadcast to (nz, ny, nx). Indexing 'ij' keeps k fastest-outermost.
    K, J, I = np.meshgrid(k, j, i, indexing='ij')
    n0 = K * (ny1 * nx1) + J * nx1 + I  # shape (nz, ny, nx)
    n0 = n0.reshape(-1)
    del K, J, I

    plane = ny1 * nx1  # nodes per z-plane

    # Compute the 8 corners as offsets from n0.
    # n0 = (k,j,i)
    # n1 = (k,j,i+1)        -> n0 + 1
    # n2 = (k,j+1,i+1)      -> n0 + nx1 + 1
    # n3 = (k,j+1,i)        -> n0 + nx1
    # n4..n7: same as n0..n3 but with k+1 -> + plane
    cols = [
        n0,
        n0 + 1,
        n0 + nx1 + 1,
        n0 + nx1,
        n0 + plane,
        n0 + plane + 1,
        n0 + plane + nx1 + 1,
        n0 + plane + nx1,
    ]

    # Write 8 columns in parallel via threads. Each tofile releases the GIL
    # while doing the actual write, so threading is fine (no need for processes).
    from concurrent.futures import ThreadPoolExecutor
    def _write_one(c, col):
        col.astype(dtype, copy=False).tofile(os.path.join(output_dir, f"i{c}.{ext}"))
    with ThreadPoolExecutor(max_workers=8) as ex:
        list(ex.map(_write_one, range(8), cols))


def _write_connectivity_chunked(nx, ny, nz, output_dir, dtype, ext):
    """Stream connectivity in z-slabs to keep peak memory bounded.

    Per slab uses ~ 8 * ny * nx * sizeof(dtype) bytes. For nx=ny=1024 that's
    ~64 MB per slab. Each output file is opened in append-binary mode and
    grown one slab at a time.
    """
    nx1 = nx + 1
    ny1 = ny + 1
    plane = ny1 * nx1

    # Open all 8 files for binary write
    files = [open(os.path.join(output_dir, f"i{c}.{ext}"), "wb") for c in range(8)]

    try:
        # Pre-build the (j,i) grid (constant across z-slabs)
        i = np.arange(nx, dtype=np.int64)
        j = np.arange(ny, dtype=np.int64)
        J, I = np.meshgrid(j, i, indexing='ij')
        n0_layer = J * nx1 + I  # shape (ny, nx)
        n0_layer = n0_layer.reshape(-1)

        for k in range(nz):
            n0 = (k * plane + n0_layer)  # shape (ny*nx,)
            cols = [
                n0,
                n0 + 1,
                n0 + nx1 + 1,
                n0 + nx1,
                n0 + plane,
                n0 + plane + 1,
                n0 + plane + nx1 + 1,
                n0 + plane + nx1,
            ]
            for c, col in enumerate(cols):
                col.astype(dtype, copy=False).tofile(files[c])
    finally:
        for f in files:
            f.close()


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Generate hex cube mesh in MARS binary format")
    p.add_argument("--nx", type=int, default=100)
    p.add_argument("--ny", type=int, default=100)
    p.add_argument("--nz", type=int, default=100)
    p.add_argument("--output", type=str, default="hex_cube_mesh")
    p.add_argument("--int32", action="store_true", help="Use int32 instead of int64")
    p.add_argument("--chunked", action="store_true",
                   help="Stream connectivity in z-slabs (lower peak memory)")
    args = p.parse_args()

    generate_hex_cube(args.nx, args.ny, args.nz, args.output,
                      use_int64=not args.int32, chunked=args.chunked)
