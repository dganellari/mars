#!/usr/bin/env python3
"""Host CPU/MPI duct study: meshes, SIMPLE runs on 1/2/4 ranks, comparator.

  run_host_study.py OUTDIR --host-run BUILD/duct_host_run --levels 4,8,16 --ranks 1,2,4
                    [--mpiexec mpiexec --numproc-flag -n --mpi-args ...] [--parity-only]

For every level it writes OUTDIR/duct-<cells>.exo/.json with duct_mesh.py (the Exodus file is
the one the CUDA runs read; the host driver builds the identical lattice itself, which the mesh
tests check), runs duct_host_run on each rank count with the study tolerances, and then runs
duct_compare.py study.

Every run leaves PREFIX.run.json: executable path and SHA-256, driver arguments, rank count,
mesh SHA-256, exit status and the SHA-256 of its fields, metrics and log. A result in OUTDIR is
reused only when that manifest matches the requested run and the files are unchanged. Anything
else is an error and is neither rerun nor overwritten: no manifest, another executable or
argument list, or edited outputs. The exception is --rerun-stale (used by ctest after rebuilds),
which deletes and reruns results that carry this tool's manifest; results without one remain an
error. A missing executable, a nonzero exit (recorded, so it persists on reuse) or a failed
launch fails the study. Exit 0 only when every run exited 0 and the study passes.
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import duct_compare as dc  # noqa: E402
import duct_mesh as dm  # noqa: E402

TOLERANCES = ["--residual-tol", "1e-8", "--mass-tol", "1e-8", "--change-tol", "1e-8"]
OUTPUTS = ("-fields.csv", "-metrics.csv", ".log")
MANIFEST = "mars-simple-duct-run-v1"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def identity(executable, digest, arguments, ranks, mesh_digest):
    return {"format": MANIFEST, "executable": executable, "executable_sha256": digest, "arguments": arguments,
            "ranks": ranks, "mesh_sha256": mesh_digest}


def cached(prefix, wanted, rerun_stale=False):
    """None when there is nothing to reuse; else the recorded exit status of a verified matching
    result. Raises ValueError for unverified, mismatched or edited results, except that with
    rerun_stale a result carrying this tool's manifest is deleted (its four files only) and
    reported as None, so it is run again."""
    manifest = prefix + ".run.json"
    present = [prefix + s for s in OUTPUTS if os.path.exists(prefix + s)]
    if not os.path.exists(manifest):
        if present:
            raise ValueError("%s: results without a manifest (%s); use a fresh OUTDIR" % (prefix, ", ".join(present)))
        return None
    with open(manifest) as f:
        record = json.load(f)
    if record.get("format") != MANIFEST:
        raise ValueError("%s: not a run_host_study manifest" % manifest)
    problem = None
    differs = [k for k in wanted if record.get(k) != wanted[k]]
    if differs:
        problem = "cached result differs in %s from the requested run" % ", ".join(differs)
    for suffix, digest in sorted(record.get("outputs", {}).items()):
        if problem is None and (not os.path.exists(prefix + suffix) or sha256(prefix + suffix) != digest):
            problem = "%s changed or vanished after the run" % suffix
    if problem is None:
        return record["exit"]
    if not rerun_stale:
        raise ValueError("%s: %s; use a fresh OUTDIR or --rerun-stale" % (prefix, problem))
    print("stale %s (%s): rerunning" % (prefix, problem), flush=True)
    for path in present + [manifest]:
        os.remove(path)
    return None


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("outdir")
    p.add_argument("--host-run", required=True)
    p.add_argument("--levels", required=True)
    p.add_argument("--ranks", default="1,2,4")
    p.add_argument("--mpiexec", default="mpiexec")
    p.add_argument("--numproc-flag", default="-n")
    p.add_argument("--mpi-args", default="", help="extra launcher arguments, space or semicolon separated")
    p.add_argument("--iterations", default="30000")
    p.add_argument("--parity-only", action="store_true", help="per-run and rank checks without refinement")
    p.add_argument("--rerun-stale", action="store_true",
                   help="delete and rerun results whose manifest no longer matches (never results without one)")
    o = p.parse_args(argv)
    levels, ranks = [int(x) for x in o.levels.split(",")], [int(x) for x in o.ranks.split(",")]
    executable = os.path.abspath(o.host_run)
    if not (os.path.isfile(executable) and os.access(executable, os.X_OK)):
        print("FAIL: host executable %s does not exist or is not executable" % executable)
        return 1
    digest = sha256(executable)
    os.makedirs(o.outdir, exist_ok=True)
    failed = []
    for c in levels:
        exo = os.path.join(o.outdir, "duct-%d.exo" % c)
        try:
            if not os.path.exists(exo) and dm.main(["--cells", str(c), "--output", exo]) != 0:
                raise ValueError("mesh generation failed")
            dc.load_mesh(exo[:-4] + ".json")   # the description must match the file
        except (OSError, ValueError, KeyError) as e:
            print("FAIL: mesh %d: %s" % (c, e))
            return 1
        mesh_digest = sha256(exo)
        for r in ranks:
            prefix = os.path.join(o.outdir, "duct-%d-%d" % (c, r))
            arguments = ["--cells", str(c), "--output-prefix", prefix, "--iterations", o.iterations,
                         "--report-every", "500"] + TOLERANCES
            wanted = identity(executable, digest, arguments, r, mesh_digest)
            try:
                status = cached(prefix, wanted, o.rerun_stale)
            except (OSError, ValueError, KeyError) as e:
                print("FAIL: %s" % e)
                return 1
            if status is None:
                command = [o.mpiexec, o.numproc_flag, str(r)] + o.mpi_args.replace(";", " ").split() + [executable] + arguments
                print(" ".join(command), flush=True)
                with open(prefix + ".log", "w") as log:
                    try:
                        status = subprocess.call(command, stdout=log, stderr=subprocess.STDOUT)
                    except OSError as e:
                        log.write("ERROR: launch failed: %s\n" % e)
                        status = 127
                record = dict(wanted, command=command, exit=status,
                              outputs=dict((s, sha256(prefix + s)) for s in OUTPUTS if os.path.exists(prefix + s)))
                with open(prefix + ".run.json", "w") as f:
                    json.dump(record, f, indent=1, sort_keys=True)
                with open(prefix + ".log") as log:
                    lines = log.read().strip().split("\n")
                print("  exit %d; %s" % (status, lines[-1] if lines else ""), flush=True)
            else:
                print("reused %s (exit %d, manifest verified)" % (prefix, status), flush=True)
            if status != 0:
                failed.append("duct-%d-%d exit %d" % (c, r, status))
    report = os.path.join(o.outdir, "study.md")
    for stale in (report, os.path.splitext(report)[0] + ".json"):
        if os.path.exists(stale):
            os.remove(stale)
    args = ["study", o.outdir, "--levels", o.levels, "--ranks", o.ranks, "--report", report]
    if o.parity_only:
        args.append("--parity-only")
    verdict = dc.main(args)
    if failed:
        print("FAIL: runs exited nonzero: %s" % ", ".join(failed))
        return 1
    return verdict


if __name__ == "__main__":
    sys.exit(main())
