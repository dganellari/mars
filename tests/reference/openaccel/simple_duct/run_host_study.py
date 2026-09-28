#!/usr/bin/env python3
"""Host CPU/MPI duct study: meshes, SIMPLE runs on 1/2/4 ranks, comparator.

  run_host_study.py OUTDIR --host-run BUILD/duct_host_run --levels 4,8,16 --ranks 1,2,4
                    [--mpiexec mpiexec --numproc-flag -n --mpi-args ...] [--parity-only]

For every level it writes OUTDIR/duct-<cells>.exo/.json with duct_mesh.py (the Exodus file is
the one the CUDA runs read; the host driver builds the identical lattice itself, which the
mesh tests check), runs duct_host_run on each rank count with the study tolerances, and then
duct_compare.py study. Existing results in OUTDIR are reused, never overwritten. Exit 0 only
when the study passes.
"""
import argparse
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import duct_compare as dc  # noqa: E402
import duct_mesh as dm  # noqa: E402

TOLERANCES = ["--residual-tol", "1e-8", "--mass-tol", "1e-8", "--change-tol", "1e-8"]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("outdir")
    p.add_argument("--host-run", required=True)
    p.add_argument("--levels", required=True)
    p.add_argument("--ranks", default="1,2,4")
    p.add_argument("--mpiexec", default="mpiexec")
    p.add_argument("--numproc-flag", default="-n")
    p.add_argument("--mpi-args", default="", help="extra launcher arguments, space separated")
    p.add_argument("--iterations", default="30000")
    p.add_argument("--parity-only", action="store_true", help="per-run and rank checks without refinement")
    o = p.parse_args(argv)
    levels, ranks = [int(x) for x in o.levels.split(",")], [int(x) for x in o.ranks.split(",")]
    os.makedirs(o.outdir, exist_ok=True)
    for c in levels:
        exo = os.path.join(o.outdir, "duct-%d.exo" % c)
        if not os.path.exists(exo) and dm.main(["--cells", str(c), "--output", exo]) != 0:
            return 1
        for r in ranks:
            prefix = os.path.join(o.outdir, "duct-%d-%d" % (c, r))
            if os.path.exists(prefix + "-fields.csv"):
                continue
            command = [o.mpiexec, o.numproc_flag, str(r)] + o.mpi_args.replace(";", " ").split() + [
                o.host_run, "--cells", str(c), "--output-prefix", prefix, "--iterations", o.iterations,
                "--report-every", "500"] + TOLERANCES
            print(" ".join(command), flush=True)
            with open(prefix + ".log", "w") as log:
                status = subprocess.call(command, stdout=log, stderr=subprocess.STDOUT)
            with open(prefix + ".log") as log:
                print("  exit %d; %s" % (status, log.read().strip().split("\n")[-1]), flush=True)
    report = os.path.join(o.outdir, "study.md")
    if os.path.exists(report):
        os.remove(report)
        os.remove(os.path.splitext(report)[0] + ".json")
    args = ["study", o.outdir, "--levels", o.levels, "--ranks", o.ranks, "--report", report]
    if o.parity_only:
        args.append("--parity-only")
    return dc.main(args)


if __name__ == "__main__":
    sys.exit(main())
