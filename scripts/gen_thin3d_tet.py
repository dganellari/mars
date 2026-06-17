#!/usr/bin/env python3
"""Generate thin-3D tetrahedral meshes for the Phase-0 coupled-solver benchmarks.

Two 2D benchmarks (lid-driven skewed cavity, backward-facing step) meshed as
one-cell-thick 3D tet slabs that MARS's AmrManager<TetTag> reads. Boundary
conditions are tagged GEOMETRICALLY in the driver (meshio cannot write Exodus
side-sets), so this only emits node coordinates + tet connectivity.

Each structured hex cell is Kuhn-split into 6 tets; node dedup at block
interfaces is by rounded coordinate. All tet volumes are forced positive
(MARS skips V<=0 tets, which would tear the mesh).

  python3 gen_thin3d_tet.py skew --beta 45 --n 32 --out skew_b45_n32.exo
  python3 gen_thin3d_tet.py bfs  --n 40 --out bfs_n40.exo
"""
import argparse
import math
import numpy as np

try:
    import meshio
except ImportError:
    raise SystemExit("meshio required: pip3 install meshio")
try:
    import netCDF4
except ImportError:
    raise SystemExit("netCDF4 required for Exodus output: pip3 install netCDF4")

# VTK hex corner order for a lattice cell (i,j,k):
_HEX = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
        (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
# 6-tet Kuhn/Freudenthal split sharing the 0-6 diagonal:
_KUHN = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6),
         (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]


def _tet_vol(p0, p1, p2, p3):
    return np.dot(np.cross(p1 - p0, p2 - p0), p3 - p0) / 6.0


def hexes_to_tet_mesh(hex_corner_coords):
    """hex_corner_coords: (Nhex, 8, 3) float64 -> (points, tets) with dedup."""
    flat = np.asarray(hex_corner_coords, dtype=np.float64).reshape(-1, 3)
    # dedup by rounding to 1e-9 (lattice coords are exact to that)
    keys = np.round(flat, 9)
    uniq, inv = np.unique(keys, axis=0, return_inverse=True)
    points = uniq.astype(np.float64)
    hexconn = inv.reshape(-1, 8)

    tets = []
    nbad = 0
    for h in hexconn:
        for t in _KUHN:
            n = [int(h[t[0]]), int(h[t[1]]), int(h[t[2]]), int(h[t[3]])]
            v = _tet_vol(points[n[0]], points[n[1]], points[n[2]], points[n[3]])
            if v < 0:                       # flip to positive orientation
                n[1], n[2] = n[2], n[1]
                nbad += 1
            tets.append(n)
    tets = np.asarray(tets, dtype=np.int64)
    return points, tets, nbad


def gen_skew(beta_deg, n, dz=None):
    """Unit parallelogram cavity: (x0,y0) in [0,1]^2 -> (x0+y0*cos b, y0*sin b)."""
    b = math.radians(beta_deg)
    cb, sb = math.cos(b), math.sin(b)
    if dz is None:
        dz = 1.0 / n
    hexes = []
    for j in range(n):
        for i in range(n):
            corners = []
            for (di, dj, dk) in _HEX:
                x0 = (i + di) / n
                y0 = (j + dj) / n
                corners.append([x0 + y0 * cb, y0 * sb, dk * dz])
            hexes.append(corners)
    return np.asarray(hexes, dtype=np.float64)


def gen_bfs(n, h=1.0, ratio=1.94, Li=5.0, Lo=20.0, dz=None):
    """Backward-facing step: inlet channel (height h, x in [-Li,0], y in [S,H])
    over the downstream channel (height H=ratio*h, x in [0,Lo], y in [0,H])."""
    H = ratio * h
    S = H - h                      # step height
    dy = h / n                     # uniform spacing from the inlet channel
    if dz is None:
        dz = dy
    nyH = int(round(H / dy))       # downstream rows
    nyh = n                        # inlet rows
    nxi = max(1, int(round(Li / dy)))
    nxo = max(1, int(round(Lo / dy)))

    hexes = []

    def block(x0, x1, y0, y1, nx, ny):
        out = []
        xs = np.linspace(x0, x1, nx + 1)
        ys = np.linspace(y0, y1, ny + 1)
        for j in range(ny):
            for i in range(nx):
                corners = []
                for (di, dj, dk) in _HEX:
                    corners.append([xs[i + di], ys[j + dj], dk * dz])
                out.append(corners)
        return out

    hexes += block(-Li, 0.0, S, H, nxi, nyh)     # inlet channel
    hexes += block(0.0, Lo, 0.0, H, nxo, nyH)     # downstream full height
    return np.asarray(hexes, dtype=np.float64)


def write_exodus(points, tets, out):
    """Write the minimal Exodus that MARS's reader consumes: separate coordx/y/z
    and a 1-based connect1 block (mars_read_exodus_mesh.hpp reads coordx/coordy/
    coordz at :114-116 and subtracts 1 from connect1 at :172). meshio's writer
    emits a combined `coord` instead, which MARS rejects with 'Variable not found'."""
    nN, nE, npe = len(points), len(tets), tets.shape[1]
    ds = netCDF4.Dataset(out, "w", format="NETCDF3_64BIT_OFFSET")
    ds.createDimension("num_nodes", nN)
    ds.createDimension("num_elem", nE)
    ds.createDimension("num_dim", 3)
    ds.createDimension("num_el_blk", 1)
    ds.createDimension("num_el_in_blk1", nE)
    ds.createDimension("num_nod_per_el1", npe)
    cx = ds.createVariable("coordx", "f8", ("num_nodes",))
    cy = ds.createVariable("coordy", "f8", ("num_nodes",))
    cz = ds.createVariable("coordz", "f8", ("num_nodes",))
    cx[:] = points[:, 0]
    cy[:] = points[:, 1]
    cz[:] = points[:, 2]
    conn = ds.createVariable("connect1", "i4", ("num_el_in_blk1", "num_nod_per_el1"))
    conn[:, :] = (tets + 1).astype(np.int32)        # Exodus connectivity is 1-based
    ds.close()


def write_mesh(points, tets, out):
    if out.endswith(".exo") or out.endswith(".e"):
        write_exodus(points, tets, out)
    else:                                            # .vtu etc. for local viz
        meshio.write(out, meshio.Mesh(points=points, cells=[("tetra", tets)]))


def report(points, tets, nbad, out):
    vols = np.array([_tet_vol(points[t[0]], points[t[1]], points[t[2]], points[t[3]])
                     for t in tets])
    print(f"[gen] {out}: nodes={len(points)} tets={len(tets)} "
          f"flipped={nbad} minV={vols.min():.3e} maxV={vols.max():.3e}")
    if vols.min() <= 0:
        raise SystemExit("ERROR: non-positive tet volume after orientation fix")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sk = sub.add_parser("skew")
    sk.add_argument("--beta", type=float, default=45.0)
    sk.add_argument("--n", type=int, default=32)
    sk.add_argument("--out", required=True)
    bf = sub.add_parser("bfs")
    bf.add_argument("--n", type=int, default=40)
    bf.add_argument("--out", required=True)
    a = ap.parse_args()

    if a.cmd == "skew":
        hexes = gen_skew(a.beta, a.n)
    else:
        hexes = gen_bfs(a.n)
    points, tets, nbad = hexes_to_tet_mesh(hexes)
    report(points, tets, nbad, a.out)
    write_mesh(points, tets, a.out)


if __name__ == "__main__":
    main()
