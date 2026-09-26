#!/usr/bin/env python3
"""Compare two SIMPLE field CSVs (node,x,y,z,u,v,w,p) node by node: a one-rank and a P-rank run.

Velocity is scaled by U=0.1 m/s and pressure by rho*U^2=0.01 Pa, as in the OpenAccel comparison.
No pressure mean is removed. Exit 1 if coordinates differ or a scaled error exceeds --tol.
"""
import argparse
import csv
import math
import sys

def load(path):
    rows = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        columns = ["node", "x", "y", "z", "u", "v", "w", "p"]
        if reader.fieldnames != columns:
            raise ValueError("{}: expected columns {}".format(path, ",".join(columns)))
        for r in reader:
            node = int(r["node"])
            if node < 0 or node in rows or None in r:
                raise ValueError("{}: invalid or duplicate node, or extra columns".format(path))
            values = {k: float(r[k]) for k in columns[1:]}
            if not all(math.isfinite(v) for v in values.values()):
                raise ValueError("{}: nonfinite field or coordinate".format(path))
            rows[node] = values
    if not rows:
        raise ValueError("{}: no field records".format(path))
    return rows

def main(argv=None):
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument("reference"); a.add_argument("candidate")
    a.add_argument("--tol", type=float, default=1e-6)
    o = a.parse_args(argv)
    try:
        if not math.isfinite(o.tol) or o.tol <= 0:
            raise ValueError("tolerance must be finite and positive")
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
        vel = max(vel, math.hypot(math.hypot(r["u"] - g["u"], r["v"] - g["v"]), r["w"] - g["w"]) / 0.1)
        pres = max(pres, abs(r["p"] - g["p"]) / 0.01)
    ok = coord <= 1e-12 and vel <= o.tol and pres <= o.tol
    print(f"{'PASS' if ok else 'FAIL'}: nodes={len(ref)} max|dx|={coord:.3e} "
          f"max velocity/U={vel:.3e} max pressure/(rho U^2)={pres:.3e} tol={o.tol:g}")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
