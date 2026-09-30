#!/usr/bin/env python3
"""Compare two Poiseuille --comparison-output frames (.pvtu) node by node.

Pieces are merged by exact Float64 coordinates. A node written by several
ranks (ghost copies) must carry identical values in every piece. Velocity is
scaled by --u and pressure by rho*U^2. Exit 1 if the node sets differ or a
scaled difference exceeds --tol.
"""
import argparse
import math
import os
import sys
import xml.etree.ElementTree as ET

FIELDS = ("u", "v", "w", "p")


def read_piece(path):
    grid = ET.parse(path).getroot().find("UnstructuredGrid/Piece")
    points = [float(t) for t in grid.find("Points/DataArray").text.split()]
    arrays = {a.get("Name"): a for a in grid.find("PointData").findall("DataArray")}
    values = {f: [float(t) for t in arrays[f].text.split()] for f in FIELDS}
    n = int(grid.get("NumberOfPoints"))
    if len(points) != 3 * n or any(len(values[f]) != n for f in FIELDS):
        raise ValueError("{}: array lengths do not match NumberOfPoints".format(path))
    for i in range(n):
        yield tuple(points[3 * i:3 * i + 3]), tuple(values[f][i] for f in FIELDS)


def load(pvtu):
    root = ET.parse(pvtu).getroot()
    pieces = [p.get("Source") for p in root.find("PUnstructuredGrid").findall("Piece")]
    if not pieces:
        raise ValueError("{}: no pieces".format(pvtu))
    nodes = {}
    for piece in pieces:
        for xyz, vals in read_piece(os.path.join(os.path.dirname(pvtu), piece)):
            if not all(math.isfinite(v) for v in xyz + vals):
                raise ValueError("{}: nonfinite value at {}".format(piece, xyz))
            if nodes.setdefault(xyz, vals) != vals:
                raise ValueError("{}: copies of node {} differ between pieces".format(piece, xyz))
    return nodes, len(pieces)


def main(argv=None):
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument("reference")
    a.add_argument("candidate")
    a.add_argument("--u", type=float, default=1.0)
    a.add_argument("--rho", type=float, default=1.0)
    a.add_argument("--tol", type=float, default=1e-6)
    o = a.parse_args(argv)
    try:
        ref, ref_pieces = load(o.reference)
        got, got_pieces = load(o.candidate)
    except (OSError, ValueError, KeyError, AttributeError, ET.ParseError) as e:
        print("FAIL: {}".format(e))
        return 1
    if ref.keys() != got.keys():
        print("FAIL: node sets differ ({} vs {} nodes)".format(len(ref), len(got)))
        return 1
    vel = pres = 0.0
    for xyz, r in ref.items():
        g = got[xyz]
        vel = max(vel, math.hypot(math.hypot(r[0] - g[0], r[1] - g[1]), r[2] - g[2]) / o.u)
        pres = max(pres, abs(r[3] - g[3]) / (o.rho * o.u * o.u))
    ok = vel <= o.tol and pres <= o.tol
    print("{}: nodes={} pieces={}/{} max velocity/U={:.3e} max pressure/(rho U^2)={:.3e} tol={:g}".format(
        "PASS" if ok else "FAIL", len(ref), ref_pieces, got_pieces, vel, pres, o.tol))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
