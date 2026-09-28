#!/usr/bin/env python3
"""Synthetic Tet4 rectangular-duct mesh for native SIMPLE (Exodus II, netCDF 64-bit offset).

Duct x in [0, L] (inlet side set "inlet" at x=0, "outlet" at x=L), y in [-W/2, W/2],
z in [-H/2, H/2], side set "walls" on the four sides. nz = cells across the height, ny = W/H * nz,
nx = L/(stretch*H) * nz: hexes of hx = stretch*h, hy = hz = h, each split into six positively
oriented Kuhn tets around the hex diagonal. The mesh is invariant under x -> x + hx, so the
discrete equations admit an exactly x-invariant (fully developed) state.

Numbering (0-based; Exodus stores +1): node i + (nx+1)(j + (ny+1)k); element t + 6 (i + nx (j + ny k)),
t the Kuhn path. The fields CSV of mars_segregated_simple reports this node number, so the
comparator recovers (i, j, k) without coordinate matching (and checks coordinates anyway).

Writes FILE.exo and FILE.json (lattice parameters, counts, SHA-256). Standard library only,
Python 3.6 compatible; `--dump FILE` writes the canonical text compared against the C++ mirror.
"""
import argparse
import array
import hashlib
import json
import os
import struct
import sys

FORMAT = "mars-simple-duct-v1"
# MARS tet_face_node (= Exodus TETRA4 side numbering minus one).
FACE_NODES = ((0, 1, 3), (1, 2, 3), (0, 3, 2), (0, 2, 1))
KUHN = ((0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0))
SIDE_SETS = ("inlet", "outlet", "walls")


class Lattice(object):
    def __init__(self, cells, length=7.0, width=2.0, height=1.0, stretch=2.0):
        self.length, self.width, self.height, self.stretch = float(length), float(width), float(height), float(stretch)
        ny, nx = cells * width / height, cells * length / (stretch * height)
        self.nz, self.ny, self.nx = int(cells), int(round(ny)), int(round(nx))
        if (cells < 2 or cells % 2 or abs(ny - self.ny) > 1e-9 or self.ny % 2 or abs(nx - self.nx) > 1e-9
                or self.nx < 1 or min(length, width, height, stretch) <= 0):
            raise ValueError("duct lattice: need even cells, W/H*cells even, L*cells/(stretch*H) integer")
        self.hx, self.hy, self.hz = self.length / self.nx, self.width / self.ny, self.height / self.nz

    @property
    def nodes(self):
        return (self.nx + 1) * (self.ny + 1) * (self.nz + 1)

    @property
    def elements(self):
        return 6 * self.nx * self.ny * self.nz

    def node(self, i, j, k):
        return i + (self.nx + 1) * (j + (self.ny + 1) * k)

    def ijk(self, g):
        i = g % (self.nx + 1)
        g //= self.nx + 1
        return i, g % (self.ny + 1), g // (self.ny + 1)

    def x(self, i):
        return i * self.hx

    def y(self, j):
        return -self.width / 2 + j * self.hy

    def z(self, k):
        return -self.height / 2 + k * self.hz

    def boundary_faces(self):
        return {"inlet": 2 * self.ny * self.nz, "outlet": 2 * self.ny * self.nz,
                "walls": 4 * self.nx * (self.ny + self.nz)}

    def kuhn(self):
        """Corner offsets of the six tets, positively oriented (same rule as duct_mesh.hpp)."""
        tets = []
        for axes in KUHN:
            c = [0, 0, 0]
            corners = [tuple(c)]
            for s in range(2):
                c[axes[s]] += 1
                corners.append(tuple(c))
            corners.append((1, 1, 1))
            h = (self.hx, self.hy, self.hz)
            m = [[(corners[r + 1][q] - corners[0][q]) * h[q] for q in range(3)] for r in range(3)]
            det = (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                   + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))
            if det < 0:
                corners[2], corners[3] = corners[3], corners[2]
            tets.append(corners)
        return tets

    def coordinates(self):
        xs, ys, zs = array.array("d"), array.array("d"), array.array("d")
        for k in range(self.nz + 1):
            for j in range(self.ny + 1):
                for i in range(self.nx + 1):
                    xs.append(self.x(i))
                    ys.append(self.y(j))
                    zs.append(self.z(k))
        return xs, ys, zs

    def connectivity(self):
        """Element-major 0-based node ids, 4 per element."""
        sj, sk = self.nx + 1, (self.nx + 1) * (self.ny + 1)
        offsets = [c[0] + sj * c[1] + sk * c[2] for tet in self.kuhn() for c in tet]
        out = array.array("i")
        for k in range(self.nz):
            for j in range(self.ny):
                base = sj * j + sk * k
                for i in range(self.nx):
                    corner = base + i
                    out.extend([corner + o for o in offsets])
        return out

    def side_sets(self):
        """{name: [(element, ordinal)]} in element order; a face lies on one plane of its hex."""
        planes = {}   # (tet type, ordinal) -> (axis, 0|1)
        for t, tet in enumerate(self.kuhn()):
            for o, face in enumerate(FACE_NODES):
                pts = [tet[n] for n in face]
                for axis in range(3):
                    values = set(p[axis] for p in pts)
                    if len(values) == 1:
                        planes[(t, o)] = (axis, values.pop())
        sets = dict((name, []) for name in SIDE_SETS)
        n = (self.nx, self.ny, self.nz)
        for k in range(self.nz):
            for j in range(self.ny):
                for i in range(self.nx):
                    hex_id = i + self.nx * (j + self.ny * k)
                    index = (i, j, k)
                    for (t, o), (axis, side) in sorted(planes.items()):
                        if (side == 0 and index[axis] != 0) or (side == 1 and index[axis] != n[axis] - 1):
                            continue
                        name = ("inlet" if side == 0 else "outlet") if axis == 0 else "walls"
                        sets[name].append((6 * hex_id + t, o))
        for name in SIDE_SETS:
            sets[name].sort()
        return sets


# ---------------------------------------------------------------- netCDF 64-bit offset writer
NC_CHAR, NC_INT, NC_FLOAT, NC_DOUBLE = 2, 4, 5, 6
SIZES = {NC_CHAR: 1, NC_INT: 4, NC_FLOAT: 4, NC_DOUBLE: 8}


def _pad(n):
    return (4 - n % 4) % 4


def _name(s):
    b = s.encode("ascii")
    return struct.pack(">i", len(b)) + b + b"\0" * _pad(len(b))


def _values(kind, values):
    if kind == NC_CHAR:
        raw = values.encode("ascii") if isinstance(values, str) else bytes(values)
        count = len(raw)
    else:
        code = {NC_INT: "i", NC_FLOAT: "f", NC_DOUBLE: "d"}[kind]
        count = len(values)
        raw = struct.pack(">%d%s" % (count, code), *values)
    return struct.pack(">ii", kind, count) + raw + b"\0" * _pad(len(raw))


def _attributes(attrs):
    if not attrs:
        return struct.pack(">ii", 0, 0)
    out = struct.pack(">ii", 12, len(attrs))
    for name, kind, values in attrs:
        out += _name(name) + _values(kind, values)
    return out


def _big_endian(data):
    data = array.array(data.typecode, data)
    if sys.byteorder == "little":
        data.byteswap()
    return data.tobytes()


def write_netcdf(path, dims, gattrs, variables):
    """dims: [(name, length)]; variables: [(name, [dim names], type, attrs, data)]. CDF-2 layout."""
    index = dict((name, i) for i, (name, _) in enumerate(dims))
    lengths = dict(dims)
    payloads = []
    for name, vdims, kind, attrs, data in variables:
        count = 1
        for d in vdims:
            count *= lengths[d]
        if kind == NC_CHAR:
            raw = bytes(data)
        else:
            raw = _big_endian(data)
        if len(raw) != count * SIZES[kind]:
            raise ValueError("variable %s: %d bytes for %d values" % (name, len(raw), count))
        payloads.append(raw + b"\0" * _pad(len(raw)))

    def header(begins):
        out = b"CDF\x02" + struct.pack(">i", 0)
        out += struct.pack(">ii", 10, len(dims)) + b"".join(_name(n) + struct.pack(">i", l) for n, l in dims)
        out += _attributes(gattrs)
        out += struct.pack(">ii", 11, len(variables))
        for (name, vdims, kind, attrs, _), payload, begin in zip(variables, payloads, begins):
            out += _name(name) + struct.pack(">i", len(vdims)) + b"".join(struct.pack(">i", index[d]) for d in vdims)
            out += _attributes(attrs) + struct.pack(">i", kind)
            out += struct.pack(">i", min(len(payload), 2 ** 32 - 1)) + struct.pack(">q", begin)
        return out
    size = len(header([0] * len(variables)))   # begin offsets have fixed width
    begins, offset = [], size
    for payload in payloads:
        begins.append(offset)
        offset += len(payload)
    with open(path, "xb") as f:
        f.write(header(begins))
        for payload in payloads:
            f.write(payload)


def _names(names, width=33):
    out = bytearray()
    for n in names:
        b = n.encode("ascii")
        if len(b) >= width:
            raise ValueError("name too long: " + n)
        out += b + b"\0" * (width - len(b))
    return out


def write_exodus(lattice, path):
    xs, ys, zs = lattice.coordinates()
    conn = lattice.connectivity()
    one = array.array("i", (v + 1 for v in conn))
    sets = lattice.side_sets()
    dims = [("len_name", 33), ("num_dim", 3), ("num_nodes", lattice.nodes), ("num_elem", lattice.elements),
            ("num_el_blk", 1), ("num_side_sets", len(SIDE_SETS)), ("num_el_in_blk1", lattice.elements),
            ("num_nod_per_el1", 4)]
    dims += [("num_side_ss%d" % (s + 1), len(sets[name])) for s, name in enumerate(SIDE_SETS)]
    gattrs = [("title", NC_CHAR, "Public synthetic SIMPLE rectangular duct (%s)" % FORMAT),
              ("api_version", NC_FLOAT, [7.22]), ("version", NC_FLOAT, [7.22]),
              ("floating_point_word_size", NC_INT, [8]), ("file_size", NC_INT, [1])]
    variables = [
        ("coor_names", ["num_dim", "len_name"], NC_CHAR, [], _names(["x", "y", "z"])),
        ("coordx", ["num_nodes"], NC_DOUBLE, [], xs),
        ("coordy", ["num_nodes"], NC_DOUBLE, [], ys),
        ("coordz", ["num_nodes"], NC_DOUBLE, [], zs),
        ("eb_names", ["num_el_blk", "len_name"], NC_CHAR, [], _names(["fluid"])),
        ("eb_status", ["num_el_blk"], NC_INT, [], array.array("i", [1])),
        ("eb_prop1", ["num_el_blk"], NC_INT, [("name", NC_CHAR, "ID")], array.array("i", [1])),
        ("connect1", ["num_el_in_blk1", "num_nod_per_el1"], NC_INT, [("elem_type", NC_CHAR, "TETRA4")], one),
        ("ss_names", ["num_side_sets", "len_name"], NC_CHAR, [], _names(SIDE_SETS)),
        ("ss_status", ["num_side_sets"], NC_INT, [], array.array("i", [1] * len(SIDE_SETS))),
        ("ss_prop1", ["num_side_sets"], NC_INT, [("name", NC_CHAR, "ID")], array.array("i", range(1, len(SIDE_SETS) + 1))),
    ]
    for s, name in enumerate(SIDE_SETS):
        variables.append(("elem_ss%d" % (s + 1), ["num_side_ss%d" % (s + 1)], NC_INT, [],
                          array.array("i", (e + 1 for e, _ in sets[name]))))
        variables.append(("side_ss%d" % (s + 1), ["num_side_ss%d" % (s + 1)], NC_INT, [],
                          array.array("i", (o + 1 for _, o in sets[name]))))
    write_netcdf(path, dims, gattrs, variables)
    return sets


def canonical(lattice):
    """Byte-exact text of the mesh (shared with duct_mesh_dump in C++)."""
    xs, ys, zs = lattice.coordinates()
    conn = lattice.connectivity()
    lines = ["duct %d %d %d %.17g %.17g %.17g" % (lattice.nx, lattice.ny, lattice.nz, lattice.length, lattice.width, lattice.height)]
    lines += ["n %.17g %.17g %.17g" % (xs[g], ys[g], zs[g]) for g in range(lattice.nodes)]
    lines += ["e %d %d %d %d" % tuple(conn[4 * e:4 * e + 4]) for e in range(lattice.elements)]
    for name, faces in sorted(lattice.side_sets().items()):
        lines += ["s %s %d %d" % (name, e, o) for e, o in faces]
    return "\n".join(lines) + "\n"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cells", type=int, required=True, help="cells across the height (even)")
    p.add_argument("--length", type=float, default=7.0)
    p.add_argument("--width", type=float, default=2.0)
    p.add_argument("--height", type=float, default=1.0)
    p.add_argument("--stretch", type=float, default=2.0, help="hx / h")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--output", help="new FILE.exo; FILE.json is written next to it")
    group.add_argument("--dump", help="write the canonical text instead")
    o = p.parse_args(argv)
    try:
        lattice = Lattice(o.cells, o.length, o.width, o.height, o.stretch)
        if o.dump:
            with open(o.dump, "x") as f:
                f.write(canonical(lattice))
            return 0
        if not o.output.endswith(".exo"):
            raise ValueError("--output must end in .exo")
        meta = o.output[:-4] + ".json"
        if os.path.exists(o.output) or os.path.exists(meta):
            raise ValueError("output exists; choose a fresh name")
        sets = write_exodus(lattice, o.output)
        info = {"format": FORMAT, "cells": lattice.nz, "nx": lattice.nx, "ny": lattice.ny, "nz": lattice.nz,
                "length": lattice.length, "width": lattice.width, "height": lattice.height, "stretch": lattice.stretch,
                "hx": lattice.hx, "hy": lattice.hy, "hz": lattice.hz, "nodes": lattice.nodes, "elements": lattice.elements,
                "side_sets": dict((n, len(sets[n])) for n in SIDE_SETS), "exodus": os.path.basename(o.output),
                "sha256": sha256(o.output)}
        expected = lattice.boundary_faces()
        if any(len(sets[n]) != expected[n] for n in SIDE_SETS):
            raise RuntimeError("side-set counts differ from the lattice formula")
        with open(meta, "x") as f:
            json.dump(info, f, indent=1, sort_keys=True)
            f.write("\n")
        print("Wrote %s: %d nodes, %d Tet4, inlet/outlet/walls %d/%d/%d faces"
              % (o.output, lattice.nodes, lattice.elements, len(sets["inlet"]), len(sets["outlet"]), len(sets["walls"])))
    except (OSError, ValueError, RuntimeError) as e:
        print("ERROR: %s" % e, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
