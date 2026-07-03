#!/usr/bin/env python3
"""mirage-emit: compile an operator spec to an IR dump or GPU/host source.

    python3 mirage-emit.py specs/laplacian.op --backend ir
    python3 mirage-emit.py specs/laplacian.op --backend cuda -o generated/laplacian_apply.cuh
    python3 mirage-emit.py specs/laplacian.op --backend host -o generated/laplacian_apply_host.hpp
"""

import argparse
import sys

from mirage import parse_spec_file, synthesize, dump
from mirage.backends import cuda, host_cpp, mlir_ir


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("spec", help="path to the .op operator spec")
    ap.add_argument("--backend", choices=["ir", "cuda", "host", "host-jvp", "mlir"],
                    default="ir",
                    help="ir = text dump (default); cuda/host = primal source; "
                         "host-jvp = Jacobian-action host source (symbolic diff); "
                         "mlir = mir-dialect partial-assembly MLIR")
    ap.add_argument("-o", "--out", help="output file (default: stdout)")
    args = ap.parse_args(argv)

    op = parse_spec_file(args.spec)
    ea = synthesize(op)

    if args.backend == "ir":
        text = dump(ea) + "\n"
    elif args.backend == "cuda":
        text = cuda.emit(ea)
    elif args.backend == "host-jvp":
        text = host_cpp.emit_jvp(ea)
    elif args.backend == "mlir":
        text = mlir_ir.emit(ea)
    else:
        text = host_cpp.emit(ea)

    if args.out:
        with open(args.out, "w") as f:
            f.write(text)
        print(f"wrote {args.out}", file=sys.stderr)
    else:
        sys.stdout.write(text)


if __name__ == "__main__":
    main()
