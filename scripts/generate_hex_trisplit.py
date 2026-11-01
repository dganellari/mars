#!/usr/bin/env python3
"""Generate an unstructured all-hex mesh in MARS binary format (i*.int64 + x.float32).

The unit square is cut into n x n squares, each square into 2 triangles, and each triangle into 3 quads
(corner, edge midpoint, centroid, edge midpoint). The quads are swept through nz layers of the unit cube.
The nodes of the triangle grid then have 6 quads around them, so 8 edge neighbours in 3D: graph rows of 9
entries and full rows of 39, against 7 and 27 on a structured hex mesh. 6 n^2 nz hexes: n = 192 and
nz = 32 give 7.1 million, as many as hex192.
"""
import argparse
import os

import numpy as np


def generate(n, nz, output_dir):
    v = lambda i, j: i + (n + 1) * j                     # triangle grid nodes
    nv = (n + 1) ** 2
    hx = lambda i, j: nv + i + n * j                     # midpoints of the edges along x, i < n
    nh = n * (n + 1)
    vy = lambda i, j: nv + nh + i + (n + 1) * j          # along y, j < n
    dg = lambda i, j: nv + 2 * nh + i + n * j            # of the diagonals (i, j)-(i+1, j+1)
    ct = lambda i, j, t: nv + 2 * nh + n * n + 2 * (i + n * j) + t  # triangle centroids
    per_layer = nv + 2 * nh + 3 * n * n

    i, j = np.meshgrid(np.arange(n + 1), np.arange(n + 1), indexing='xy')
    x2 = np.zeros(per_layer)
    y2 = np.zeros(per_layer)
    x2[v(i, j)], y2[v(i, j)] = i / n, j / n
    i, j = np.meshgrid(np.arange(n), np.arange(n + 1), indexing='xy')
    x2[hx(i, j)], y2[hx(i, j)] = (i + .5) / n, j / n
    i, j = np.meshgrid(np.arange(n + 1), np.arange(n), indexing='xy')
    x2[vy(i, j)], y2[vy(i, j)] = i / n, (j + .5) / n
    i, j = np.meshgrid(np.arange(n), np.arange(n), indexing='xy')
    x2[dg(i, j)], y2[dg(i, j)] = (i + .5) / n, (j + .5) / n
    # triangle 0: (i,j), (i+1,j), (i+1,j+1); triangle 1: (i,j), (i+1,j+1), (i,j+1), both counterclockwise
    x2[ct(i, j, 0)], y2[ct(i, j, 0)] = (i + 2 / 3) / n, (j + 1 / 3) / n
    x2[ct(i, j, 1)], y2[ct(i, j, 1)] = (i + 1 / 3) / n, (j + 2 / 3) / n

    # the 3 counterclockwise quads of each triangle: corner, midpoint, centroid, midpoint
    i, j = i.ravel(), j.ravel()
    quads = []
    for t, (a, b, c, ab, bc, ca) in enumerate([
            (v(i, j), v(i + 1, j), v(i + 1, j + 1), hx(i, j), vy(i + 1, j), dg(i, j)),
            (v(i, j), v(i + 1, j + 1), v(i, j + 1), dg(i, j), hx(i, j + 1), vy(i, j))]):
        g = ct(i, j, t)
        quads += [np.stack([a, ab, g, ca]), np.stack([b, bc, g, ab]), np.stack([c, ca, g, bc])]
    quads = np.concatenate(quads, axis=1)                # 4 x 6 n^2

    os.makedirs(output_dir, exist_ok=True)
    layers = np.arange(nz + 1)
    np.tile(x2, nz + 1).astype(np.float32).tofile(os.path.join(output_dir, 'x.float32'))
    np.tile(y2, nz + 1).astype(np.float32).tofile(os.path.join(output_dir, 'y.float32'))
    np.repeat(layers / nz, per_layer).astype(np.float32).tofile(os.path.join(output_dir, 'z.float32'))
    k = np.arange(nz)
    for c in range(8):
        col = quads[c % 4][None, :] + per_layer * (k[:, None] + c // 4)
        col.astype(np.int64).ravel().tofile(os.path.join(output_dir, f'i{c}.int64'))
    print(f'{6 * n * n * nz} hexes, {per_layer * (nz + 1)} nodes in {output_dir}')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('-n', type=int, default=192, help='squares per side of the triangle grid')
    p.add_argument('--nz', type=int, default=32, help='layers')
    p.add_argument('--output', default='hex_trisplit')
    a = p.parse_args()
    generate(a.n, a.nz, a.output)
