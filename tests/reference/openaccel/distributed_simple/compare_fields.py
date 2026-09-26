#!/usr/bin/env python3
"""Compare two SIMPLE field CSVs (node,x,y,z,u,v,w,p) node by node: a one-rank and a P-rank run.

Velocity is scaled by U=0.1 m/s and pressure by rho*U^2=0.01 Pa, as in the OpenAccel comparison.
No pressure mean is removed. Exit 1 if coordinates differ or a scaled error exceeds --tol.
"""
import argparse, csv, sys

def load(path):
    with open(path, newline="") as f:
        rows = {int(r["node"]): r for r in csv.DictReader(f)}
    return rows

def main():
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument("reference"); a.add_argument("candidate")
    a.add_argument("--tol", type=float, default=1e-6)
    o = a.parse_args()
    ref, got = load(o.reference), load(o.candidate)
    if sorted(ref) != sorted(got):
        print("FAIL: node sets differ"); return 1
    coord = vel = pres = 0.0
    for n, r in ref.items():
        g = got[n]
        coord = max(coord, *(abs(float(r[k]) - float(g[k])) for k in "xyz"))
        vel = max(vel, sum((float(r[k]) - float(g[k])) ** 2 for k in "uvw") ** 0.5 / 0.1)
        pres = max(pres, abs(float(r["p"]) - float(g["p"])) / 0.01)
    ok = coord <= 1e-12 and vel <= o.tol and pres <= o.tol
    print(f"{'PASS' if ok else 'FAIL'}: nodes={len(ref)} max|dx|={coord:.3e} "
          f"max velocity/U={vel:.3e} max pressure/(rho U^2)={pres:.3e} tol={o.tol:g}")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
