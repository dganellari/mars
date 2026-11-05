#!/usr/bin/env python3
"""Place one-rank and distributed tests inside explicit scheduler steps."""
import argparse
import os
import sys


def command(args, layout):
    if args.ranks not in (1, 2, 4) or layout not in ("single", "multi"):
        raise ValueError("Expected 1/2/4 ranks and single/multi layout")
    nodes = min(2, args.ranks) if layout == "multi" else 1
    return [args.launcher, "-n", str(args.ranks), "--nodes=" + str(nodes),
            "--ntasks-per-node=" + str(args.ranks // nodes), "--export=ALL"] + args.launcher_arg + args.binary


def main():
    if sys.argv[1:] == ["--version"]:
        print("MARS structured scheduler launcher 1")
        return
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-n", dest="ranks", type=int, required=True)
    parser.add_argument("--launcher", required=True)
    parser.add_argument("--launcher-arg", action="append", default=[])
    parser.add_argument("binary", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if not args.binary:
        parser.error("Missing executable")
    argv = command(args, os.environ.get("MARS_GATE_LAYOUT", "single"))
    os.execv(argv[0], argv)


if __name__ == "__main__":
    main()
