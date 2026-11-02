#!/usr/bin/env python3
"""The deformed cube of mars_marsir_bp as an Exodus Hex8 mesh, scrambled, for the --mesh gate.

The vertices are those of the built-in mesh (x += a sin(pi x) sin(pi y) sin(pi z) (1,1,1)),
so mars_marsir_bp --mesh on this file solves the same discrete system as the built-in cube:
the same iterations and, to rounding, the same history. What the reader and the tables must
undo is scrambled: node and element numbering are random permutations, every element's corners
are listed in one of the 24 proper rotations of its frame, and --mirrored of them in a
left-handed order (the driver mirrors those back).

Run: python3 scripts/write_bp_cube_exodus.py --ne 16 [--deform 0.05] [--mirrored 0.1] out.e
"""
import argparse

import numpy as np
from netCDF4 import Dataset

ap = argparse.ArgumentParser()
ap.add_argument("--ne", type=int, default=16)
ap.add_argument("--deform", type=float, default=0.05)
ap.add_argument("--mirrored", type=float, default=0.1, help="fraction of elements listed left-handed")
ap.add_argument("--seed", type=int, default=1)
ap.add_argument("out")
args = ap.parse_args()
rng = np.random.default_rng(args.seed)
ne = args.ne

# Vertices of the lattice, deformed as lattice_vertex in mars_marsir_bp.cu.
g = np.arange(ne + 1) / ne
X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
bump = args.deform * np.sin(np.pi * X) * np.sin(np.pi * Y) * np.sin(np.pi * Z)
coords = np.stack([X + bump, Y + bump, Z + bump], -1).reshape(-1, 3)
lattice_id = lambda i, j, k: (i * (ne + 1) + j) * (ne + 1) + k

# VTK = Exodus Hex8 corner order (cellwise::corner_bits).
bits = np.array([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)])
corner_of_bits = {tuple(b): c for c, b in enumerate(bits)}


def frames(proper):
    """Corner relabelings: corner c of the new frame is corner out[c] of the element."""
    out = []
    for perm in [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]:
        inversions = sum(perm[a] > perm[b] for a in range(3) for b in range(a + 1, 3))
        for flips in range(8):
            f = [(flips >> d) & 1 for d in range(3)]
            if (inversions + sum(f)) % 2 != (0 if proper else 1):
                continue
            out.append([corner_of_bits[tuple(bits[c][perm[d]] ^ f[d] for d in range(3))] for c in range(8)])
    return out


proper, improper = frames(True), frames(False)
assert len(proper) == 24 and len(improper) == 24

conn = []
for i in range(ne):
    for j in range(ne):
        for k in range(ne):
            corners = [lattice_id(i + b[0], j + b[1], k + b[2]) for b in bits]
            frame = (improper if rng.random() < args.mirrored else proper)[rng.integers(24)]
            conn.append([corners[frame[c]] for c in range(8)])
conn = np.array(conn, dtype=np.int64)

node_perm = rng.permutation(len(coords))   # new id of old node
new_coords = np.empty_like(coords)
new_coords[node_perm] = coords
conn = node_perm[conn][rng.permutation(len(conn))]

with Dataset(args.out, "w", format="NETCDF3_64BIT_OFFSET") as f:
    f.createDimension("num_dim", 3)
    f.createDimension("num_nodes", len(new_coords))
    f.createDimension("num_elem", len(conn))
    f.createDimension("num_el_blk", 1)
    f.createDimension("num_el_in_blk1", len(conn))
    f.createDimension("num_nod_per_el1", 8)
    c = f.createVariable("connect1", "i4", ("num_el_in_blk1", "num_nod_per_el1"))
    c.elem_type = "HEX8"
    c[:] = conn + 1
    for a, name in enumerate(["coordx", "coordy", "coordz"]):
        f.createVariable(name, "f8", ("num_nodes",))[:] = new_coords[:, a]
print(f"{args.out}: {len(conn)} elements, {len(new_coords)} nodes, ~{args.mirrored:.0%} listed left-handed")
