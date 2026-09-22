#!/usr/bin/env python3
"""run only the bundled public smoke case; preserve logs and exit status."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run(bundle, executable, output, capture_interior=False):
    bundle, executable, output = (p.resolve() for p in (bundle, executable, output))
    manifest = json.loads((bundle / "manifest.json").read_text())
    if manifest["fixture"] != "public_channel" or manifest["iterations"] != 2:
        raise ValueError("this runner accepts only the two-iteration public channel smoke case")
    names = ["input.i", "channel.exo", "run_openaccel_public.py"]
    if capture_interior:
        names += ["openaccel_reference_check.py", "contract_v1.json", "provenance.json"]
        if manifest.get("require_node_capture"):
            names += ["openaccel_node_check.py"]
    for name in names:
        if sha256(bundle / name) != manifest["sha256"][name]:
            raise ValueError("bundle checksum mismatch: " + name)
    # One wrapper owns the log/output directory. Multi-rank runs need a separate launcher.
    for key in ("SLURM_NTASKS", "OMPI_COMM_WORLD_SIZE", "PMI_SIZE"):
        if int(os.environ.get(key, "1")) != 1:
            raise ValueError("smoke wrapper requires one MPI rank: " + key)
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise ValueError("executable is missing or not executable: " + str(executable))
    output.mkdir(parents=True, exist_ok=False)
    for name in ("input.i", "channel.exo", "manifest.json"):
        shutil.copyfile(bundle / name, output / name)
    environment = os.environ.copy()
    # A smoke run must not accidentally inherit a previous export destination.
    environment.pop("MARS_OPENACCEL_EXPORT_DIR", None)
    environment.pop("MARS_OPENACCEL_PUBLIC_FIXTURE", None)
    if capture_interior:
        (output / "exports").mkdir()
        environment["MARS_OPENACCEL_EXPORT_DIR"] = str(output / "exports")
        environment["MARS_OPENACCEL_PUBLIC_FIXTURE"] = "public_channel"
        shutil.copyfile(bundle / "provenance.json", output / "provenance.json")
    environment["OMP_NUM_THREADS"] = "1"
    record = {"fixture": "public_channel",
              "binary": str(executable), "binary_sha256": sha256(executable),
              "bundle": manifest, "python_version": sys.version,
              "numerical_parity_verified": False,
              "convergence_verified": False, "status": "started"}
    result_path = output / "run.json"
    result_path.write_text(json.dumps(record, indent=2) + "\n")
    with (output / "run.log").open("w") as log:
        # Daint's python3 may be 3.6; the `text` alias requires 3.7.
        process = subprocess.Popen([str(executable), "-i", "input.i"], cwd=str(output),
                                   env=environment, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, universal_newlines=True,
                                   errors="replace")
        with process.stdout:
            for line in process.stdout:
                log.write(line)
                log.flush()
                print(line, end="", flush=True)
        returncode = process.wait()
    log_text = (output / "run.log").read_text()
    iterations = [int(x) for x in re.findall(r"^Iter = (\d+)\s*$", log_text, re.MULTILINE)]
    completed = returncode == 0 and iterations == [1, 2]
    record.update(returncode=returncode, iteration_headers=iterations,
                  status="runtime_smoke_completed" if completed else "runtime_smoke_failed")
    if capture_interior:
        record["coverage"] = "local-interior-only"
        record["full_contract_passed"] = False
        if completed:
            try:
                from openaccel_reference_check import load_dump
                capture = load_dump(output / "exports",
                                    require_inputs=bool(manifest.get("require_frozen_inputs", False)))
                expected = {(stage, call) for stage in ("momentum.interior", "pressure.interior")
                            for call in (1, 2)}
                keys = capture["blocks"]
                if ({key[:2] for key in keys} != expected
                        or any(sum(key[:2] == pair for key in keys) != 1536 for pair in expected)
                        or len(capture["hashes"]) != 4
                        or capture["producer"] != "openaccel"
                        or capture["signature"][0] != "public_channel"):
                    raise ValueError("expected two calls per stage, one rank, 1536 public elements each")
                record["export_sha256"] = capture["hashes"]
                record["frozen_input_blocks"] = len(capture["inputs"])
                record["status"] = "interior_capture_completed"
                if manifest.get("require_node_capture"):
                    from openaccel_node_check import load_nodes
                    nodes, node_hashes = load_nodes(output / "exports" / "nodes")
                    if {call for stage, call, node in nodes} != {1, 2} or any(
                            sum(s == stage and c == call for s, c, n in nodes) != 425
                            for stage in range(3) for call in (1, 2)):
                        raise ValueError("expected two node calls with 425 public nodes each")
                    if not any(stage == 3 for stage, call, node in nodes):
                        raise ValueError("public boundary relaxation was not captured")
                    record.update(status="node_capture_completed", coverage="local-interior-and-steady-nodes",
                                  node_records=len(nodes), node_sha256=node_hashes)
            except (ValueError, KeyError, TypeError, OSError, OverflowError) as error:
                completed = False
                record.update(status="interior_capture_failed", capture_error=str(error))
                print("Capture check failed: " + str(error), file=sys.stderr)
    result_path.write_text(json.dumps(record, indent=2) + "\n")
    print("Public smoke: " + record["status"] + "; log: " + str(output / "run.log"))
    print("Two iterations do not establish convergence or numerical parity.")
    return 0 if completed else (returncode if 0 < returncode < 256 else 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--capture-interior", action="store_true",
                        help="require actual interior exports from the instrumented public reference")
    args = parser.parse_args()
    try:
        return run(args.bundle, args.executable, args.output, args.capture_interior)
    except (OSError, ValueError, KeyError) as error:
        print("ERROR: " + str(error), file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
