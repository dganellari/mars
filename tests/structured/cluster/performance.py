#!/usr/bin/env python3
"""Compare identical checked workloads against the pinned pre-follow-up source."""
import argparse
import json
import math
from pathlib import Path
import statistics
import subprocess
import sys

sys.dont_write_bytecode = True
import run

PREFIX = "MARS_STRUCTURED_PERF "
CASES = {"morton-quad-rect": (128, 96, 0), "hilbert-quad-rect": (128, 96, 0),
         "hilbert-quad-full": (128, 128, 0), "morton-hex-rect": (32, 24, 16),
         "hilbert-hex-rect": (32, 24, 16), "hilbert-hex-full": (32, 32, 32)}
METRICS = {"mesh_seconds": 50.e-6, "dof_seconds": 50.e-6, "halo_seconds": 5.e-6}


def audit(text, ranks, backend, transport, nodes, threads):
    receipts = run.audit_receipts(text, backend, ranks, transport, nodes, threads)
    rows = [json.loads(line.split(PREFIX, 1)[1]) for line in text.splitlines() if PREFIX in line]
    if len(rows) != len(CASES) or {row.get("case") for row in rows} != set(CASES):
        raise ValueError("Missing or duplicate performance workloads")
    for row in rows:
        nx, ny, nz = CASES[row["case"]]
        if (tuple(row.get(k) for k in ("nx", "ny", "nz")) != (nx, ny, nz) or
                row.get("ranks") != ranks or row.get("iterations") != 20 or row.get("errors") != 0 or
                row.get("cells") != nx * ny * (nz or 1) or
                row.get("owned_dofs") != (nx + 1) * (ny + 1) * (nz + 1 if nz else 1) or
                any(not isinstance(row.get(k), (int, float)) or not math.isfinite(row[k]) or row[k] <= 0 for k in METRICS)):
            raise ValueError("Performance workload failed correctness, shape or timing checks")
    return {"receipts": receipts, "cases": {row["case"]: row for row in rows}}


def compare(samples):
    comparisons = []
    for layout in sorted({s["layout"] for s in samples}):
        for ranks in (1, 2, 4):
            if layout == "multi" and ranks == 1:
                continue
            for case in CASES:
                for metric, allowance in METRICS.items():
                    values = {}
                    for variant in ("reference", "candidate"):
                        group = [s for s in samples if s["layout"] == layout and s["ranks"] == ranks and s["variant"] == variant]
                        if len(group) != 5 or sorted(s["round"] for s in group) != list(range(5)):
                            raise ValueError("Missing or duplicate performance rounds")
                        values[variant] = [s["audit"]["cases"][case][metric] for s in group]
                    reference, candidate = (statistics.median(values[v]) for v in ("reference", "candidate"))
                    comparisons.append({"layout": layout, "ranks": ranks, "case": case, "metric": metric,
                        "reference_seconds": reference, "candidate_seconds": candidate,
                        "ratio": candidate / reference, "absolute_allowance_seconds": allowance,
                        "pass": candidate <= 1.10 * reference + allowance})
    if not comparisons:
        raise ValueError("No performance comparisons")
    return comparisons


def execute(command, log, env):
    with log.open("w") as stream:
        result = run.run_child(command, stdout=stream, stderr=subprocess.STDOUT, env=env)
    if result.returncode:
        raise ValueError("Command failed: " + str(log))


def validate(args):
    if args.work.exists():
        raise ValueError("Use a new performance result directory")
    args.work.mkdir(parents=True)
    result = {"status": "FAILED", "samples": [], "commands": [], "driver_files": run.driver_files(), "policy":
              "five alternating reference/candidate rounds; median <= 1.10*reference + 50us setup/DOFs or 5us halo"}
    try:
        candidate = json.loads((args.candidate / "result.json").read_text())
        if candidate.get("status") != "PASS":
            raise ValueError("Candidate capability/install validation must pass first")
        for path, checksum in candidate["artifacts"].items():
            if run.digest(path) != checksum:
                raise ValueError("Validated candidate artifact changed: " + path)
        result["source_manifest_sha256"] = run.verify_snapshot(args.source)
        if candidate["source_manifest_sha256"] != result["source_manifest_sha256"]:
            raise ValueError("Candidate source differs from capability run")
        result["reference_manifest_sha256"] = run.verify_snapshot(args.reference)
        result["candidate_result_sha256"] = run.digest(args.candidate / "result.json")
        result["backend"], result["transport"] = candidate["backend"], candidate["transport"]
        env = run.execution_environment(candidate["environment"])
        env["OMP_NUM_THREADS"] = str(args.threads)
        env["KOKKOS_NUM_THREADS"] = str(args.threads)
        env["OMP_PROC_BIND"] = "false"
        result["environment"] = run.environment_receipt(env)
        steps = dict(json.loads((args.candidate / "commands.json").read_text()))
        if args.reuse:
            prepared = json.loads((args.reuse / "result.json").read_text())
            if prepared.get("status") != "BUILT_NOT_RUN":
                raise ValueError("Reuse requires an unrun prepared performance build")
            for key in ("source_manifest_sha256", "reference_manifest_sha256", "candidate_result_sha256"):
                if prepared[key] != result[key]:
                    raise ValueError("Prepared provenance mismatch: " + key)
            binaries = {key: Path(value) for key, value in prepared["binary_paths"].items()}
            for path, checksum in prepared["binaries"].items():
                if run.digest(path) != checksum:
                    raise ValueError("Prepared performance binary changed")
            result["prepared_result_sha256"] = run.digest(args.reuse / "result.json")
        else:
            reference_work = args.work / "reference"
            for phase in ("configure", "build", "install"):
                command = [value.replace(str(args.source), str(args.reference)).replace(str(args.candidate), str(reference_work)) for value in steps[phase]]
                if phase == "build":
                    command = command[:command.index("--target") + 1] + ["mars"]
                result["commands"].append(command)
                execute(command, args.work / ("reference-" + phase + ".log"), env)
            binaries = {}
            for variant, install in (("reference", reference_work / "install"), ("candidate", args.candidate / "install")):
                build = args.work / (variant + "-performance")
                command = steps["consumer-configure"][:]
                command[command.index("-S") + 1] = str(args.source / "tests/structured/performance")
                command[command.index("-B") + 1] = str(build)
                command = ["-DCMAKE_PREFIX_PATH=" + str(install) if value.startswith("-DCMAKE_PREFIX_PATH=") else value for value in command]
                result["commands"].append(command)
                execute(command, args.work / (variant + "-configure.log"), env)
                command = [steps["build"][0], "--build", str(build), "--parallel", "8"]
                result["commands"].append(command)
                execute(command, args.work / (variant + "-build.log"), env)
                binaries[variant] = build / "mars_structured_performance"
        result["binary_paths"] = {key: str(value) for key, value in binaries.items()}
        result["binaries"] = {str(p): run.digest(p) for p in binaries.values()}
        if args.prepare_only:
            result["status"] = "BUILT_NOT_RUN"
            return result
        listing = json.loads((args.candidate / "consumer-list.log").read_text())
        templates = {int(test["name"].split("mpi")[-1]) if "mpi" in test["name"] else 1: test["command"] for test in listing["tests"]}
        layouts = ("single", "multi") if args.multi_node else ("single",)
        for layout in layouts:
            env["MARS_GATE_LAYOUT"] = layout
            for ranks in (1, 2, 4):
                if layout == "multi" and ranks == 1:
                    continue
                for repeat in range(5):
                    order = ("reference", "candidate") if repeat % 2 == 0 else ("candidate", "reference")
                    for variant in order:
                        command = [str(binaries[variant]) if Path(value).name == "mars_structured_consumer" else value for value in templates[ranks]]
                        log = args.work / ("{}-np{}-{}-{}.log".format(layout, ranks, repeat, variant))
                        result["commands"].append(command)
                        execute(command, log, env)
                        checked = audit(log.read_text(), ranks, result["backend"], result["transport"],
                                        2 if layout == "multi" else 1, args.threads)
                        result["samples"].append({"layout": layout, "ranks": ranks, "round": repeat,
                                                  "variant": variant, "audit": checked})
        result["comparisons"] = compare(result["samples"])
        if run.verify_snapshot(args.source) != result["source_manifest_sha256"] or run.verify_snapshot(args.reference) != result["reference_manifest_sha256"]:
            raise ValueError("Source changed during performance run")
        if run.driver_files() != result["driver_files"]:
            raise ValueError("Performance drivers changed during run")
        result["binaries"] = {str(p): run.digest(p) for p in binaries.values()}
        result["status"] = "PASS" if all(x["pass"] for x in result["comparisons"]) else "TIMING_REJECTED"
    except Exception as error:
        result["error"] = str(error)
        raise
    finally:
        result["artifacts"] = {str(p): run.digest(p) for p in args.work.rglob("*") if p.is_file() and (p.suffix == ".log" or p.name == "CMakeCache.txt")}
        (args.work / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("candidate", "reference", "source", "work"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--multi-node", action="store_true")
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--reuse", type=Path)
    args = parser.parse_args()
    for name in ("candidate", "reference", "source", "work"):
        setattr(args, name, getattr(args, name).resolve())
    if args.reuse:
        args.reuse = args.reuse.resolve()
    if args.threads < 1 or (args.prepare_only and args.reuse):
        parser.error("Invalid threads or prepare/reuse combination")
    try:
        result = validate(args)
    except Exception as error:
        print("FAILED: " + str(error), file=sys.stderr)
        return 1
    print(result["status"] + ": " + str(args.work / "result.json"))
    return 0 if result["status"] in ("PASS", "BUILT_NOT_RUN") else 1


if __name__ == "__main__":
    sys.exit(main())
