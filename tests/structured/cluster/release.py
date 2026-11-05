#!/usr/bin/env python3
"""Prepare once, then run the structured CPU/CUDA/HIP distributed release gates."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.dont_write_bytecode = True
import run


def configurations(backend):
    return [(space, mode) for space in (("serial", "openmp") if backend == "cpu" else (backend,))
            for mode in ("direct", "staged")]


def execute(argv, log, env=None):
    print("Running " + log.stem, flush=True)
    with log.open("w") as output:
        process = run.run_child(argv, stdout=output, stderr=subprocess.STDOUT, env=env)
    if process.returncode:
        raise ValueError("Failed " + str(log))


def prepare(args, result):
    kit = args.source / "tests/structured/cluster"
    common = ["--cxx", str(args.cxx), "--arch", args.arch, "--jobs", str(args.jobs)]
    if args.gpu_compiler:
        common += ["--gpu-compiler", str(args.gpu_compiler)]
    for backend in ("serial", "openmp") if args.backend == "cpu" else (args.backend,):
        deps = args.work / (backend + "-kokkos")
        execute([sys.executable, str(kit / "build_deps.py"), "--source", str(args.source.parent / "dependencies/kokkos"),
                 "--work", str(deps), "--backend", backend] + common, args.work / (backend + "-kokkos.log"))
        for transport in ("direct", "staged"):
            name = backend + "-" + transport
            target = args.work / name
            threads = args.threads if backend == "openmp" else 1
            argv = [sys.executable, str(kit / "run.py"), "--source", str(args.source), "--work", str(target),
                    "--backend", backend, "--transport", transport, "--threads", str(threads),
                    "--kokkos-dir", str(deps / "install/lib/cmake/Kokkos"),
                    "--gtest-source", str(args.source.parent / "dependencies/googletest"),
                    "--launcher", str(kit / "mpi_launch.py"), "--launcher-arg=--launcher=" + str(args.launcher)]
            argv += ["--launcher-arg=--launcher-arg=" + value for value in args.launcher_arg] + common
            execute(argv, args.work / (name + ".log"))
            execute([sys.executable, str(kit / "performance.py"), "--candidate", str(target),
                     "--source", str(args.source), "--reference", str(args.reference),
                     "--work", str(args.work / (name + "-performance-build")), "--threads", str(threads),
                     "--prepare-only"], args.work / (name + "-performance-build.log"))
            result["configurations"].append(name)
    if args.backend == "cpu":
        build = args.work / "kokkos-off"
        argv = ["cmake", "-S", str(args.source), "-B", str(build), "-G", "Ninja", "-DCMAKE_BUILD_TYPE=Release",
                "-DCMAKE_CXX_COMPILER=" + str(args.cxx), "-DMARS_ENABLE_KOKKOS=OFF", "-DMARS_ENABLE_MPI=ON",
                "-DMARS_ENABLE_CUDA=OFF", "-DMARS_ENABLE_HIP=OFF", "-DMARS_ENABLE_UNSTRUCTURED=OFF",
                "-DMARS_ENABLE_SERIAL_BACKEND=OFF", "-DMARS_ENABLE_AMR_BACKEND=OFF",
                "-DMARS_ENABLE_EXAMPLES=OFF", "-DMARS_ENABLE_TESTS=OFF"]
        execute(argv, args.work / "kokkos-off-configure.log")
        execute(["cmake", "--build", str(build), "--target", "mars", "--parallel", str(args.jobs)],
                args.work / "kokkos-off-build.log")
        result["kokkos_off_core_build"] = "PASS"
    result["status"] = "PREPARED_NOT_RELEASE_VALIDATED"


def validate(args, result):
    # Updated orchestration can reuse an immutable, already validated source build.
    kit = Path(__file__).resolve().parent
    prepared = json.loads((args.prepared / "result.json").read_text())
    if (prepared.get("status") != "PREPARED_NOT_RELEASE_VALIDATED" or
            prepared.get("source_manifest_sha256") != result["source_manifest_sha256"] or
            prepared.get("reference_manifest_sha256") != result["reference_manifest_sha256"] or
            prepared.get("backend") != args.backend):
        raise ValueError("Preparation/source/backend mismatch")
    result["prepared_result_sha256"] = run.digest(args.prepared / "result.json")
    for backend, transport in configurations(args.backend):
        name = backend + "-" + transport
        target = args.prepared / name
        checked = json.loads((target / "result.json").read_text())
        if checked.get("status") != "PASS" or checked.get("source_manifest_sha256") != result["source_manifest_sha256"]:
            raise ValueError("Preparation did not validate " + name)
        for path, checksum in dict(checked["artifacts"], **checked["dependencies"]["kokkos_files"]).items():
            if run.digest(path) != checksum:
                raise ValueError("Prepared artifact changed: " + path)
        env = run.execution_environment(checked["environment"])
        threads = args.threads if backend == "openmp" else 1
        env.update(MARS_GATE_LAYOUT="multi", OMP_NUM_THREADS=str(threads), KOKKOS_NUM_THREADS=str(threads), OMP_PROC_BIND="false")
        commands = dict(json.loads((target / "commands.json").read_text()))
        record = {"name": name, "multinode": {}, "environment": run.environment_receipt(env)}
        for phase, expected, gtest in (("tests", run.SUITES, True), ("consumer-tests", run.CONSUMER, False)):
            argv = commands[phase][:]
            xml = args.work / (name + "-" + phase + ".xml")
            argv[argv.index("--output-junit") + 1] = str(xml)
            execute(argv, args.work / (name + "-" + phase + ".log"), env)
            record["multinode"][phase] = run.audit_junit(xml, expected, backend, transport, gtest, 2, threads)
        perf = args.work / (name + "-performance")
        argv = [sys.executable, str(kit / "performance.py"), "--candidate", str(target), "--source", str(args.source),
                "--reference", str(args.reference), "--work", str(perf), "--threads", str(threads), "--multi-node",
                "--reuse", str(args.prepared / (name + "-performance-build"))]
        # A timing rejection must not hide the other transport's correctness results.
        try:
            execute(argv, args.work / (name + "-performance.log"), env)
        except ValueError:
            if not (perf / "result.json").is_file():
                raise
        performance = json.loads((perf / "result.json").read_text())
        record["performance"] = performance["status"]
        record["performance_result_sha256"] = run.digest(perf / "result.json")
        result["configurations"].append(record)
    result["status"] = "PASS" if all(x["performance"] == "PASS" for x in result["configurations"]) else "FAILED_OR_TIMING_REJECTED"
    result["scope"] = "Structured gates only; complete release also requires the unstructured GPU regression jobs"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("source", "reference", "work"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--prepared", type=Path)
    parser.add_argument("--phase", choices=("prepare", "validate"), required=True)
    parser.add_argument("--backend", choices=("cpu", "cuda", "hip"), required=True)
    parser.add_argument("--cxx", type=Path)
    parser.add_argument("--gpu-compiler", type=Path)
    parser.add_argument("--arch", default="native")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--launcher", type=Path)
    parser.add_argument("--launcher-arg", action="append", default=[])
    args = parser.parse_args()
    for name in ("source", "reference", "work", "prepared"):
        if getattr(args, name):
            setattr(args, name, getattr(args, name).resolve())
    if args.work.exists() or args.threads < 2 or args.jobs < 1:
        parser.error("Fresh work directory, at least two OpenMP threads and positive build jobs required")
    if args.phase == "prepare" and (not args.cxx or not args.launcher or (args.backend != "cpu" and not args.gpu_compiler)):
        parser.error("Preparation needs compiler and launcher paths")
    if args.phase == "validate" and not args.prepared:
        parser.error("Validation needs --prepared")
    args.work.mkdir(parents=True)
    result = {"status": "FAILED", "backend": args.backend, "configurations": [], "threads": args.threads,
              "driver_files": run.driver_files()}
    try:
        result["source_manifest_sha256"] = run.verify_snapshot(args.source)
        result["reference_manifest_sha256"] = run.verify_snapshot(args.reference)
        (prepare if args.phase == "prepare" else validate)(args, result)
        if run.verify_snapshot(args.source) != result["source_manifest_sha256"] or run.verify_snapshot(args.reference) != result["reference_manifest_sha256"]:
            raise ValueError("Source changed during run")
        if run.driver_files() != result["driver_files"]:
            raise ValueError("Validation drivers changed during run")
    except Exception as error:
        result.update(status="FAILED", error=str(error))
        print("FAILED: " + str(error), file=sys.stderr)
    result["artifacts"] = {str(p): run.digest(p) for p in args.work.glob("*") if p.is_file() and p.suffix in (".xml", ".log")}
    (args.work / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(result["status"] + ": " + str(args.work / "result.json"))
    return 0 if result["status"] in ("PASS", "PREPARED_NOT_RELEASE_VALIDATED") else 1


if __name__ == "__main__":
    sys.exit(main())
