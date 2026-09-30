#!/usr/bin/env python3
"""Compare two SIMPLE field CSVs or distributed JSON manifests (node,x,y,z,u,v,w,p) node by node: a one-rank and a P-rank run.

Velocity is scaled by U and pressure by rho*U^2 (defaults U=.1, rho=1).
No pressure mean is removed. Exit 1 if coordinates differ or a scaled error exceeds --tol.
"""
import argparse
import csv
import math
import json
from pathlib import Path
import sys

def load(path):
    path = Path(path)
    expected = None
    paths = [path]
    if path.suffix == ".json":
        manifest = json.loads(path.read_text())
        if not isinstance(manifest, dict) or manifest.get("format") != "mars-simple-fields-v1":
            raise ValueError("invalid SIMPLE field manifest")
        expected, parts = manifest.get("nodes"), manifest.get("parts")
        if type(expected) is not int or expected <= 0 or not isinstance(parts, list) or not parts:
            raise ValueError("invalid manifest node count or parts")
        if any(not isinstance(p, str) or Path(p).name != p or not p.endswith(".csv") for p in parts):
            raise ValueError("manifest parts must be CSV filenames")
        if len(set(parts)) != len(parts):
            raise ValueError("duplicate manifest part")
        paths = [path.parent / p for p in parts]
    rows = {}
    for part in paths:
        with open(part, newline="") as f:
            reader = csv.DictReader(f)
            columns = ["node", "x", "y", "z", "u", "v", "w", "p"]
            if reader.fieldnames != columns:
                raise ValueError("{}: expected columns {}".format(part, ",".join(columns)))
            for r in reader:
                node = int(r["node"])
                if node < 0 or node in rows or None in r:
                    raise ValueError("{}: invalid or duplicate node, or extra columns".format(part))
                values = {k: float(r[k]) for k in columns[1:]}
                if not all(math.isfinite(v) for v in values.values()):
                    raise ValueError("{}: nonfinite field or coordinate".format(part))
                rows[node] = values
    if not rows:
        raise ValueError("{}: no field records".format(path))
    if expected is not None and (len(rows) != expected or any(n not in rows for n in range(expected))):
        raise ValueError("manifest does not cover every source node exactly once")
    return rows

def main(argv=None):
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument("reference"); a.add_argument("candidate")
    a.add_argument("--tol", type=float, default=1e-6)
    a.add_argument("--rho", type=float, default=1.)
    a.add_argument("--inlet-velocity", type=float, default=.1)
    o = a.parse_args(argv)
    try:
        if not math.isfinite(o.tol) or o.tol <= 0:
            raise ValueError("tolerance must be finite and positive")
        if not all(math.isfinite(v) and v > 0 for v in (o.rho, o.inlet_velocity)):
            raise ValueError("density and velocity scale must be finite and positive")
        pressure_scale = o.rho * o.inlet_velocity * o.inlet_velocity
        if not math.isfinite(pressure_scale) or pressure_scale <= 0:
            raise ValueError("invalid pressure scale")
        ref, got = load(o.reference), load(o.candidate)
    except (OSError, ValueError, TypeError, KeyError, csv.Error) as e:
        print("FAIL: {}".format(e)); return 1
    if sorted(ref) != sorted(got):
        print("FAIL: node sets differ"); return 1
    coord = vel = pres = 0.0
    for n, r in ref.items():
        g = got[n]
        coord = max(coord, *(abs(r[k] - g[k]) for k in "xyz"))
        # Nested hypot avoids overflow from squaring a large but finite difference.
        vel = max(vel, math.hypot(math.hypot(r["u"] - g["u"], r["v"] - g["v"]), r["w"] - g["w"]) / o.inlet_velocity)
        pres = max(pres, abs(r["p"] - g["p"]) / pressure_scale)
    ok = coord <= 1e-12 and vel <= o.tol and pres <= o.tol
    print(f"{'PASS' if ok else 'FAIL'}: nodes={len(ref)} max|dx|={coord:.3e} "
          f"max velocity/U={vel:.3e} max pressure/(rho U^2)={pres:.3e} tol={o.tol:g}")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
