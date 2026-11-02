#!/usr/bin/env python3
# Element lines of a Hex8 Exodus mesh, as the node-line smoother builds them: each element's thin
# axis is the local axis with the shortest mean edge; an element and its neighbour are chained
# when their shared face lies across the thin axis of both. Reports the line lengths, and how
# many chain links the driver's partition (Morton order of the centroids, equal chunks) cuts.
# Usage: hex_mesh_lines.py mesh.e [ranks,...]
import sys
import numpy as np
import netCDF4 as nc

corner_ref = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                       [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]])
EDGES = np.array([[[0, 1], [3, 2], [4, 5], [7, 6]], [[0, 3], [1, 2], [4, 7], [5, 6]],
                  [[0, 4], [1, 5], [2, 6], [3, 7]]])


def spread(v):   # 21 bits, two zeros between (mars_marsir_bp.cu morton_spread)
    v = v.astype(np.uint64) & np.uint64(0x1fffff)
    for s, m in ((32, 0x1f00000000ffff), (16, 0x1f0000ff0000ff), (8, 0x100f00f00f00f00f),
                 (4, 0x10c30c30c30c30c3), (2, 0x1249249249249249)):
        v = (v | (v << np.uint64(s))) & np.uint64(m)
    return v


def components(nb):   # labels of the paths and cycles of a graph of degree at most 2 (nb: -1 = none)
    lab = np.arange(len(nb)); idx = np.arange(len(nb))
    n0 = np.where(nb[:, 0] >= 0, nb[:, 0], idx); n1 = np.where(nb[:, 1] >= 0, nb[:, 1], idx)
    rounds = 0
    while True:
        new = np.minimum(lab, np.minimum(lab[n0], lab[n1]))
        new = new[new]; new = new[new]
        rounds += 1
        if np.array_equal(new, lab):
            return lab, rounds
        lab = new


def report(name, lab):
    sizes = np.bincount(lab)[lab]   # each element's line length
    bins = [(1, 1), (2, 9), (10, 99), (100, 999), (1000, 10 ** 12)]
    parts = ", ".join(f"{lo}-{hi if hi < 10 ** 12 else ''}: {100.0 * np.mean((sizes >= lo) & (sizes <= hi)):.1f}%"
                      for lo, hi in bins)
    print(f"{name}: {len(np.unique(lab))} lines, longest {sizes.max()}, elements by their line's length: {parts}; "
          f"element-weighted mean length {sizes.mean():.0f}")
    return sizes


path = sys.argv[1]
ranks = [int(r) for r in sys.argv[2].split(",")] if len(sys.argv) > 2 else [64, 512, 4096]
f = nc.Dataset(path)
X = np.stack([f.variables[c][:].data for c in ("coordx", "coordy", "coordz")], 1)
conn = np.concatenate([f.variables[f"connect{b + 1}"][:].data for b in range(len(f.dimensions["num_el_blk"]))]) - 1
E = len(conn)
axis = np.empty(E, np.int8); ar = np.empty(E); cent = np.empty((E, 3))
for s in range(0, E, 1 << 20):
    p = X[conn[s:s + (1 << 20)]]
    L = np.linalg.norm(p[:, EDGES[:, :, 1]] - p[:, EDGES[:, :, 0]], axis=3).mean(2)
    axis[s:s + len(p)] = np.argmin(L, 1); ar[s:s + len(p)] = L.max(1) / L.min(1); cent[s:s + len(p)] = p.mean(1)

# The two faces across each element's thin axis; a face listed by two elements chains them.
sel = np.stack([corner_ref[:, k] == t for k in range(3) for t in (-1, 1)])   # (6, 8) corner masks
faces = np.empty((E, 2, 4), np.int64)
for side, t in enumerate((-1, 1)):
    for k in range(3):
        m = axis == k
        faces[m, side] = conn[m][:, sel[2 * k + (t > 0)]]
faces = np.sort(faces.reshape(-1, 4), 1)
key_hi = (faces[:, 0] << 24) | faces[:, 1]; key_lo = (faces[:, 2] << 24) | faces[:, 3]
order = np.lexsort((key_lo, key_hi))
same = (key_hi[order[1:]] == key_hi[order[:-1]]) & (key_lo[order[1:]] == key_lo[order[:-1]])
a, b = order[:-1][same], order[1:][same]   # (element * 2 + side) of both ends of each link
nb = np.full((E, 2), -1, np.int64)
nb[a // 2, a % 2] = b // 2; nb[b // 2, b % 2] = a // 2
links = len(a)
print(f"{path}: {E} Hex8; {links} chain links ({100.0 * np.mean((nb >= 0).any(1)):.1f}% of elements have one)")
lab, rounds = components(nb)
sizes = report("element lines", lab)
for lo, hi in ((1, 3), (3, 10), (10, 100), (100, 1e9)):
    m = (ar >= lo) & (ar < hi)
    print(f"  aspect ratio {lo:g}-{hi:g}: {100.0 * m.mean():.1f}% of elements, in lines of length >= 2: "
          f"{100.0 * np.mean(sizes[m] >= 2):.1f}%, median line length {np.median(sizes[m]):.0f}")

lo, hi = X.min(0), X.max(0)
fr = np.clip((cent - lo) / (hi - lo) * 2097152.0, 0, 2097151).astype(np.uint64)
key = spread(fr[:, 0]) << np.uint64(2) | spread(fr[:, 1]) << np.uint64(1) | spread(fr[:, 2])
rank_order = np.argsort(key, kind="stable")
for P in ranks:
    owner = np.empty(E, np.int64)
    first = [E // P * r + min(r, E % P) for r in range(P + 1)]
    for r in range(P):
        owner[rank_order[first[r]:first[r + 1]]] = r
    cut = nb.copy()
    m = cut >= 0
    cut[m & (owner[np.maximum(cut, 0)] != owner[:, None])] = -1
    ncut = (links * 2 - (cut >= 0).sum()) // 2
    print(f"{P} ranks: {ncut} of {links} links cut ({100.0 * ncut / links:.1f}%); ", end="")
    report("pieces", components(cut)[0])
