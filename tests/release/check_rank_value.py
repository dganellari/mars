#!/usr/bin/env python3
"""Run a solver example on 1 rank and on N ranks and require the same printed value.

A periodic point, a ghost node or a halo element handled wrongly on N ranks changes the
solution, so a value such as the kinetic energy must not depend on the rank count. The last
match of --regex in each run is compared; both runs must exit 0 and print a finite value.

    check_rank_value.py --np 4 --regex 'KE=([-+0-9.eE]+)' --mpiexec srun --numproc-flag -n -- <driver> <args...>
"""
import argparse
import math
import re
import subprocess
import sys


def run(cmd, pattern):
    print("+ " + " ".join(cmd), flush=True)
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    sys.stdout.write(res.stdout)
    if res.returncode != 0:
        sys.exit(f"FAIL: exit code {res.returncode}")
    found = pattern.findall(res.stdout)
    if not found:
        sys.exit(f"FAIL: pattern {pattern.pattern!r} not found in the output")
    value = float(found[-1])
    if not math.isfinite(value):
        sys.exit(f"FAIL: non-finite value {found[-1]}")
    return value


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--np", type=int, required=True)
    p.add_argument("--regex", required=True)
    p.add_argument("--mpiexec", required=True)
    p.add_argument("--numproc-flag", default="-n")
    p.add_argument("--preflags", default="", help="extra mpiexec flags, space separated")
    # The examples print 11 significant digits; the runs agree to all of them in practice.
    p.add_argument("--rtol", type=float, default=1e-8)
    p.add_argument("driver", nargs=argparse.REMAINDER)
    a = p.parse_args()
    driver = a.driver[1:] if a.driver and a.driver[0] == "--" else a.driver
    if not driver:
        sys.exit("usage: missing driver command after --")
    pattern = re.compile(a.regex)

    def launch(n):
        return [a.mpiexec, a.numproc_flag, str(n)] + a.preflags.split() + driver

    v1 = run(launch(1), pattern)
    vn = run(launch(a.np), pattern)
    denom = max(abs(v1), abs(vn))
    rel = abs(vn - v1) / denom if denom > 0.0 else 0.0
    print(f"1 rank : {v1:.10e}")
    print(f"{a.np} ranks: {vn:.10e}  (relative difference {rel:.2e}, rtol {a.rtol:.0e})")
    if rel > a.rtol:
        sys.exit("FAIL: the value depends on the rank count")
    print("PASS")


if __name__ == "__main__":
    main()
