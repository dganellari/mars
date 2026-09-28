#!/usr/bin/env python3
"""Write a synthetic 4 x 1 x 1 Tet4 duct for performance runs; never reads a case mesh."""
import argparse
from pathlib import Path
import sys

import numpy as np
from netCDF4 import Dataset

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / "scripts"))
from generate_tet_cube import KUHN_TETS


def write(path, nx, ny, nz):
    if min(nx, ny, nz) < 1 or 6 * nx * ny * nz > np.iinfo(np.int32).max // 16:
        raise ValueError("positive dimensions within the native reader's index capacity required")
    path = Path(path)
    if path.exists():
        raise ValueError("output exists")
    # Vertex IDs and element order match the small public reference generator.
    k, j, i = np.meshgrid(np.arange(nz + 1), np.arange(ny + 1), np.arange(nx + 1), indexing="ij")
    xyz = np.column_stack((4. * i.ravel() / nx, j.ravel() / ny, k.ravel() / nz))
    plane = (nx + 1) * (ny + 1)
    k, j, i = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx), indexing="ij")
    bases = (k * plane + j * (nx + 1) + i).ravel()
    offsets = np.array([0, 1, nx + 2, nx + 1, plane, plane + 1, plane + nx + 2, plane + nx + 1])
    tets = (bases[:, None, None] + offsets[KUHN_TETS][None, :, :]).reshape(-1, 4)
    faces = ((0, 1, 3), (1, 2, 3), (0, 3, 2), (0, 2, 1))
    patches = [[], [], []]
    # Planar tests select only exterior faces, without an O(elements) Python dictionary.
    for side, face in enumerate(faces, 1):
        points = xyz[tets[:, face]]
        inlet = np.all(points[:, :, 0] == 0, axis=1)
        outlet = np.all(points[:, :, 0] == 4, axis=1)
        walls = np.any(np.all(points[:, :, 1:] == 0, axis=1) |
                       np.all(points[:, :, 1:] == 1, axis=1), axis=1)
        for patch, mask in zip(patches, (inlet, outlet, walls)):
            elements = np.flatnonzero(mask) + 1
            patch.append(np.column_stack((elements, np.full(len(elements), side))))
    with Dataset(path, "w", format="NETCDF3_64BIT_OFFSET") as mesh:
        mesh.title = "PUBLIC_SYNTHETIC_SIMPLE_PERFORMANCE_CHANNEL"
        mesh.api_version = np.float32(7.22)
        mesh.version = np.float32(7.22)
        mesh.floating_point_word_size = np.int32(8)
        mesh.file_size = np.int32(1)
        for name, count in dict(num_dim=3, num_nodes=len(xyz), num_elem=len(tets), num_el_blk=1,
                                num_el_in_blk1=len(tets), num_nod_per_el1=4, num_side_sets=3, len_name=33).items():
            mesh.createDimension(name, count)
        for axis, name in enumerate("xyz"):
            mesh.createVariable("coord" + name, "f8", ("num_nodes",))[:] = xyz[:, axis]
        connect = mesh.createVariable("connect1", "i4", ("num_el_in_blk1", "num_nod_per_el1"))
        connect.elem_type = "TETRA4"
        connect[:] = tets + 1
        for prefix, dimension, ids in (("eb", "num_el_blk", [1]), ("ss", "num_side_sets", [1, 2, 3])):
            prop = mesh.createVariable(prefix + "_prop1", "i4", (dimension,))
            prop.setncattr("name", "ID")
            prop[:] = ids
            mesh.createVariable(prefix + "_status", "i4", (dimension,))[:] = 1
        names = np.zeros((3, 33), dtype="S1")
        for number, (name, pieces) in enumerate(zip(("inlet", "outlet", "walls"), patches), 1):
            entries = np.concatenate(pieces)
            dimension = "num_side_ss" + str(number)
            mesh.createDimension(dimension, len(entries))
            for column, kind in enumerate(("elem", "side")):
                mesh.createVariable(kind + "_ss" + str(number), "i4", (dimension,))[:] = entries[:, column]
            names[number - 1, :len(name)] = np.frombuffer(name.encode("ascii"), dtype="S1")
        mesh.createVariable("ss_names", "S1", ("num_side_sets", "len_name"))[:] = names


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    for name, default in (("nx", 64), ("ny", 16), ("nz", 16)):
        parser.add_argument("--" + name, type=int, default=default)
    args = parser.parse_args()
    write(args.output, args.nx, args.ny, args.nz)
    print("Wrote synthetic channel: " + str(args.output))
