#!/usr/bin/env python3
"""Run a CVFEM graph-assembly driver on 1 rank and on N ranks and compare the norms it prints.

The graph kernels assemble owned rows only, so the global matrix and RHS norms must not depend
on the rank count. A missing halo contribution or a node owned by two ranks changes them.

    check_rank_invariance.py --np 4 --mpiexec mpiexec --numproc-flag -n -- <driver> <args...>
"""
import argparse
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_rows import compare as compare_rows  # noqa: E402

NORMS = re.compile(r"matrix norm:\s*([-+0-9.eE]+),\s*rhs norm:\s*([-+0-9.eE]+)")


def run(cmd, row_prefix=None):
    print("+ " + " ".join(cmd), flush=True)
    env = dict(os.environ)
    if row_prefix:
        env["MARS_ROW_DUMP"] = row_prefix
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env)
    sys.stdout.write(res.stdout)
    if res.returncode != 0:
        sys.exit(f"FAIL: exit code {res.returncode}")
    m = NORMS.search(res.stdout)
    if not m:
        sys.exit("FAIL: no 'matrix norm: ..., rhs norm: ...' line in the output")
    return float(m.group(1)), float(m.group(2))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--np", type=int, required=True)
    p.add_argument("--mpiexec", required=True)
    p.add_argument("--numproc-flag", default="-n")
    p.add_argument("--preflags", default="", help="extra mpiexec flags, space separated")
    # The drivers print norms with 7 significant digits.
    p.add_argument("--rtol", type=float, default=1e-5)
    p.add_argument("--rows", metavar="DIR",
                   help="also dump owned rows (MARS_ROW_DUMP) into DIR and compare them node by node")
    p.add_argument("driver", nargs=argparse.REMAINDER)
    a = p.parse_args()
    driver = a.driver[1:] if a.driver and a.driver[0] == "--" else a.driver
    if not driver:
        sys.exit("usage: missing driver command after --")

    def launch(n):
        return [a.mpiexec, a.numproc_flag, str(n)] + a.preflags.split() + driver

    rows1 = rowsn = None
    if a.rows:
        os.makedirs(a.rows, exist_ok=True)
        rows1, rowsn = os.path.join(a.rows, "np1"), os.path.join(a.rows, f"np{a.np}")
        for pre in (rows1, rowsn):
            for fn in [f for f in os.listdir(a.rows) if f.startswith(os.path.basename(pre) + ".rank")]:
                os.remove(os.path.join(a.rows, fn))
    m1, r1 = run(launch(1), rows1)
    mn, rn = run(launch(a.np), rowsn)
    print(f"1 rank : matrix norm {m1:.6e}  rhs norm {r1:.6e}")
    print(f"{a.np} ranks: matrix norm {mn:.6e}  rhs norm {rn:.6e}")

    if not m1 > 0.0:
        sys.exit("FAIL: zero matrix norm on 1 rank (nothing was assembled)")
    ok = True
    for name, ref, val in (("matrix", m1, mn), ("rhs", r1, rn)):
        denom = max(abs(ref), abs(val))
        rel = abs(val - ref) / denom if denom > 0.0 else 0.0
        if rel > a.rtol:
            print(f"FAIL: {name} norm differs by {rel:.2e} (rtol {a.rtol:.0e})")
            ok = False
    if a.rows:
        problems = compare_rows(rows1, rowsn, 1e-10)
        if problems:
            print("FAIL: " + "\n".join(problems))
            ok = False
    if not ok:
        sys.exit(1)
    print("PASS: rank-count invariant" + (" (norms and every owned row)" if a.rows else " (norms)"))


if __name__ == "__main__":
    main()
