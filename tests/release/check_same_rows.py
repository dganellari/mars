#!/usr/bin/env python3
"""Assemble two meshes that differ only in element corner numbering and require identical rows.

    check_same_rows.py --rows DIR --mesh-a A --mesh-b B -- <launcher...> <driver> --mesh={mesh} ...

"{mesh}" in the command is replaced by each mesh in turn. Correct assembly depends on geometry
only, so every owned row (matched by node SFC key) must agree; a wrong reference-to-physical
gradient transform shows up as soon as the element axes are not aligned with x, y, z.
"""
import argparse
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_rows import compare as compare_rows  # noqa: E402


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--rows", required=True, metavar="DIR")
    p.add_argument("--mesh-a", required=True)
    p.add_argument("--mesh-b", required=True)
    p.add_argument("--rtol", type=float, default=1e-9)
    p.add_argument("cmd", nargs=argparse.REMAINDER)
    a = p.parse_args()
    cmd = a.cmd[1:] if a.cmd and a.cmd[0] == "--" else a.cmd
    if not cmd or not any("{mesh}" in c for c in cmd):
        sys.exit("usage: command after -- must contain {mesh}")

    os.makedirs(a.rows, exist_ok=True)
    prefixes = []
    for tag, mesh in (("a", a.mesh_a), ("b", a.mesh_b)):
        prefix = os.path.join(a.rows, tag)
        for fn in [f for f in os.listdir(a.rows) if f.startswith(tag + ".rank")]:
            os.remove(os.path.join(a.rows, fn))
        run = [c.replace("{mesh}", mesh) for c in cmd]
        print("+ " + " ".join(run), flush=True)
        res = subprocess.run(run, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                             env=dict(os.environ, MARS_ROW_DUMP=prefix))
        sys.stdout.write(res.stdout)
        if res.returncode != 0:
            sys.exit(f"FAIL: exit code {res.returncode} on {mesh}")
        prefixes.append(prefix)

    problems = compare_rows(prefixes[0], prefixes[1], a.rtol)
    if problems:
        print("FAIL: " + "\n".join(problems))
        sys.exit(1)
    print("PASS: identical rows for both corner numberings")


if __name__ == "__main__":
    main()
