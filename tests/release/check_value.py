#!/usr/bin/env python3
"""Run a command, require exit code 0, and check that a number it prints lies in [lo, hi].

    check_value.py --regex 'Max:\\s*([-+0-9.eE]+)' --lo 0.0534 --hi 0.0590 -- <command...>
"""
import argparse
import re
import subprocess
import sys


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--regex", required=True, help="pattern whose first group is the number")
    p.add_argument("--lo", type=float, required=True)
    p.add_argument("--hi", type=float, required=True)
    p.add_argument("cmd", nargs=argparse.REMAINDER)
    a = p.parse_args()
    cmd = a.cmd[1:] if a.cmd and a.cmd[0] == "--" else a.cmd
    if not cmd:
        sys.exit("usage: missing command after --")

    print("+ " + " ".join(cmd), flush=True)
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    sys.stdout.write(res.stdout)
    if res.returncode != 0:
        sys.exit(f"FAIL: exit code {res.returncode}")
    m = re.search(a.regex, res.stdout)
    if not m:
        sys.exit(f"FAIL: pattern {a.regex!r} not found in the output")
    val = float(m.group(1))
    if not (a.lo <= val <= a.hi):
        sys.exit(f"FAIL: {val:.6e} outside [{a.lo:.6e}, {a.hi:.6e}]")
    print(f"PASS: {val:.6e} in [{a.lo:.6e}, {a.hi:.6e}]")


if __name__ == "__main__":
    main()
