#!/usr/bin/env python3
"""Build and validate one native structured Kokkos configuration offline."""

import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import signal
import sys
import xml.etree.ElementTree as ET


GTEST_SHA = "15460959cbbfa20e66ef0b5ab497367e47fc0a04"
SUITES = {"morton_test": 1, "hilbert_test": 1, "mesh_manager_test": 1,
          "dof_handler_test": 1, "mesh_ghost_layer_test": 1}
for _name in ("mesh_manager_test", "dof_handler_test", "mesh_ghost_layer_test"):
    for _ranks in (2, 4):
        SUITES["{}_mpi{}".format(_name, _ranks)] = _ranks
CONSUMER = {"structured_install": 1, "structured_install_mpi2": 2, "structured_install_mpi4": 4}
RECEIPT = "MARS_STRUCTURED_RECEIPT "
TEST_RECEIPT = "MARS_STRUCTURED_TEST_RECEIPT "
GTEST_CASES = {"morton_test": 4, "hilbert_test": 5, "mesh_manager_test": 14,
               "dof_handler_test": 17, "mesh_ghost_layer_test": 23}
TEST_ENVIRONMENT = ("OMP_NUM_THREADS", "MPICH_GPU_SUPPORT_ENABLED", "MPICH_GPU_IPC_ENABLED",
                    "MARS_STRUCTURED_EXPECT_BACKEND")
RUNTIME_ENVIRONMENT = TEST_ENVIRONMENT + ("SLURM_JOB_ID", "SLURM_JOB_NODELIST", "SLURM_NTASKS",
    "SLURM_STEP_ID", "CUDA_VISIBLE_DEVICES", "HIP_VISIBLE_DEVICES", "ROCR_VISIBLE_DEVICES", "MARS_GATE_LAYOUT")


def execution_environment(recorded):
    # Scheduler identity and device visibility belong to the current allocation.
    env = dict(os.environ)
    for key in list(env):
        if key.startswith(("MARS_STRUCTURED_", "GTEST_")):
            del env[key]
    env.update({key: recorded[key] for key in TEST_ENVIRONMENT if key in recorded})
    return env


def environment_receipt(env):
    return {key: env[key] for key in RUNTIME_ENVIRONMENT if key in env}


def driver_files():
    directory = Path(__file__).resolve().parent
    return {str(directory / name): digest(directory / name)
            for name in ("run.py", "release.py", "performance.py", "mpi_launch.py")}


def digest(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def verify_snapshot(source):
    manifest_path = source.parent / "source-manifest.json"
    manifest = json.loads(manifest_path.read_text())
    expected = manifest["files"]
    actual = {str(path.relative_to(source)) for path in source.rglob("*") if path.is_file()}
    if actual != set(expected):
        raise ValueError("Source snapshot has missing or extra files")
    for name, checksum in expected.items():
        path = source / name
        if Path(name).is_absolute() or ".." in Path(name).parts or path.is_symlink() or digest(path) != checksum:
            raise ValueError("Source snapshot checksum mismatch: " + name)
    return digest(manifest_path)


def verify_dependency(source, revision):
    manifest_path = source / "dependency-manifest.json"
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("revision") != revision:
        raise ValueError("Wrong dependency revision: " + str(source))
    files = manifest["files"]
    actual = {str(path.relative_to(source)) for path in source.rglob("*") if path.is_file() and path != manifest_path}
    if actual != set(files):
        raise ValueError("Dependency snapshot has missing or extra files: " + str(source))
    for name, checksum in files.items():
        path = source / name
        if Path(name).is_absolute() or ".." in Path(name).parts or path.is_symlink() or digest(path) != checksum:
            raise ValueError("Dependency source checksum mismatch: " + name)
    return digest(manifest_path)


def cmake_value(text, variable):
    match = re.search(r"\bSET\(\s*" + re.escape(variable) + r"\s+([^)]*)\)", text, re.I)
    return match.group(1).strip().strip('"') if match else ""


def unblock_child_signal():
    # Some batch environments block SIGCHLD, which prevents CMake from reaping its tools.
    signal.pthread_sigmask(signal.SIG_UNBLOCK, {signal.SIGCHLD})


def child_signal_policy():
    if not hasattr(signal, "pthread_sigmask"):
        return {"available": False}
    mask = signal.pthread_sigmask(signal.SIG_BLOCK, set())
    return {"available": True, "parent_sigchld_blocked": signal.SIGCHLD in mask,
            "child_sigchld_unblocked": True}


def run_child(command, **kwargs):
    # The kit is single threaded; keep the parent mask and all other child signals unchanged.
    if hasattr(signal, "pthread_sigmask"):
        kwargs["preexec_fn"] = unblock_child_signal
    return subprocess.run(command, **kwargs)


def tool_version(executable, log_path, required=True):
    resolved = shutil.which(str(executable))
    if not resolved:
        if required:
            raise ValueError("Executable not found: " + str(executable))
        return {"path": str(executable), "available": False}
    path = Path(resolved).absolute()
    with log_path.open("w") as log:
        process = run_child([str(path), "--version"], stdout=log, stderr=subprocess.STDOUT)
    banner = log_path.read_text()
    if required and (process.returncode or not banner.strip()):
        raise ValueError("Could not identify compiler/tool version: " + str(path))
    return {"path": str(path), "resolved_binary": str(path.resolve()), "sha256": digest(path),
            "version_exit": process.returncode, "version_banner": banner}


def is_gpu(backend):
    return backend in ("cuda", "hip")


def verify_dependencies(args):
    version = (args.kokkos_dir / "KokkosConfigVersion.cmake").read_text()
    common = (args.kokkos_dir / "KokkosConfigCommon.cmake").read_text()
    config = (args.kokkos_dir / "KokkosConfig.cmake").read_text()
    if cmake_value(version, "PACKAGE_VERSION") != "4.4.1":
        raise ValueError("The gate requires Kokkos 4.4.01 (package version 4.4.1)")
    if args.backend.upper() not in cmake_value(common, "Kokkos_DEVICES").split(";"):
        raise ValueError("Kokkos does not export the requested execution space")
    if is_gpu(args.backend):
        if "COMPILE_AS_CMAKE_LANGUAGE" not in cmake_value(common, "Kokkos_OPTIONS").split(";"):
            raise ValueError("Rebuild Kokkos with Kokkos_ENABLE_COMPILE_AS_CMAKE_LANGUAGE=ON")
        if cmake_value(config, "Kokkos_COMPILE_LANGUAGE") != args.backend.upper():
            raise ValueError("Kokkos does not export the matching native GPU language")
        architectures = cmake_value(common, "Kokkos_{}_ARCHITECTURES".format(args.backend.upper())).split(";")
        if args.arch not in architectures:
            raise ValueError("Kokkos target architectures do not contain --arch=" + args.arch)
    gtest_manifest = verify_dependency(args.gtest_source, GTEST_SHA)
    dependency_files = {}
    # Header, package and library hashes identify the installed dependency actually used.
    prefix = args.kokkos_dir.parent.parent.parent
    if prefix in (Path("/usr"), Path("/usr/local"), Path("/opt/homebrew"), Path("/opt/rocm")):
        raise ValueError("Use a dedicated Kokkos install prefix; shared system prefixes are not accepted")
    for directory in (prefix / "include", prefix / "lib", prefix / "lib64"):
        if directory.is_dir():
            for path in sorted(directory.rglob("*")):
                if path.is_file():
                    dependency_files[str(path)] = digest(path)
    if not dependency_files:
        raise ValueError("No installed Kokkos dependency files found")
    return {"kokkos_version": "4.4.01", "kokkos_files": dependency_files,
            "googletest_sha": GTEST_SHA, "googletest_manifest_sha256": gtest_manifest}


def audit_receipts(output, backend, ranks, transport, min_nodes=1, threads=1):
    receipts = []
    for line in output.splitlines():
        if RECEIPT in line:
            receipts.append(json.loads(line.split(RECEIPT, 1)[1]))
    if len(receipts) != ranks or sorted(receipt.get("rank") for receipt in receipts) != list(range(ranks)):
        raise ValueError("Missing or duplicate per-rank runtime receipts")
    hosts = {}
    gpu = is_gpu(backend)
    for receipt in receipts:
        if (str(receipt.get("backend", "")).upper() != backend.upper() or receipt.get("ranks") != ranks or
                receipt.get("kokkos_version") != 40401 or
                receipt.get("transport") != transport or
                receipt.get("mpi_device_buffer") is not (gpu and transport == "direct") or
                not receipt.get("hostname") or not receipt.get("mpi_library") or
                not isinstance(receipt.get("concurrency"), int) or receipt["concurrency"] < threads):
            raise ValueError("Runtime receipt does not match execution space, threads or transport")
        if gpu:
            if (not isinstance(receipt.get("device"), int) or receipt["device"] < 0 or
                    not receipt.get("device_id") or receipt["device_id"] == "cpu" or
                    not receipt.get("device_name") or receipt["device_name"] == "cpu" or
                    not isinstance(receipt.get("runtime_version"), int) or receipt["runtime_version"] <= 0 or
                    not isinstance(receipt.get("driver_version"), int) or receipt["driver_version"] <= 0):
                raise ValueError("Missing native GPU runtime/device evidence")
        elif (receipt.get("device") != -1 or receipt.get("device_id") != "cpu" or
              receipt.get("device_name") != "cpu" or receipt.get("runtime_version") != 0):
            raise ValueError("CPU run used an unexpected device")
        hosts.setdefault(receipt["hostname"], []).append(receipt)
    if len(hosts) < min(min_nodes, ranks):
        raise ValueError("MPI test did not span the required number of hosts")
    for records in hosts.values():
        if (any(record.get("local_ranks") != len(records) for record in records) or
                sorted(record.get("local_rank") for record in records) != list(range(len(records))) or
                (gpu and len({record["device_id"] for record in records}) != len(records))):
            raise ValueError("Local-rank or physical-device assignment is inconsistent")
    return receipts


def audit_junit(path, expected, backend, transport, gtest, min_nodes=1, threads=1):
    root = ET.parse(path).getroot()
    cases = list(root.iter("testcase"))
    names = [case.attrib.get("name") for case in cases]
    if len(names) != len(expected) or set(names) != set(expected):
        raise ValueError("CTest results do not contain the exact required test matrix")
    results = {}
    for case in cases:
        name = case.attrib["name"]
        if case.find("failure") is not None or case.find("error") is not None or case.find("skipped") is not None:
            raise ValueError("CTest failed or skipped: " + name)
        if case.attrib.get("status", "run") not in ("run", "passed"):
            raise ValueError("CTest did not run: " + name)
        output = "\n".join(element.text or "" for element in case.findall("system-out"))
        if gtest:
            records = [json.loads(line.split(TEST_RECEIPT, 1)[1])
                       for line in output.splitlines() if TEST_RECEIPT in line]
            cases = GTEST_CASES[re.sub(r"_mpi[24]$", "", name)]
            if (len(records) != expected[name] or
                    sorted(record.get("rank") for record in records) != list(range(expected[name])) or
                    any(record.get("ranks") != expected[name] or record.get("cases") != cases or
                        record.get("passed") != cases or record.get("skipped") != 0 for record in records)):
                raise ValueError("GoogleTest did not pass the complete suite on every rank: " + name)
        results[name] = audit_receipts(output, backend, expected[name], transport, min_nodes, threads)
    return results


def test_binary(test):
    basename = "mars_structured_consumer" if test["name"] == "structured_install" else test["name"]
    executables = [Path(argument) for argument in test["command"] if Path(argument).name == basename]
    if len(executables) != 1:
        raise ValueError("CTest command does not identify exactly one expected binary: " + test["name"])
    return executables[0]


def plan(args):
    gpu = args.backend.upper()
    build, install, consumer = args.work / "build", args.work / "install", args.work / "consumer"
    launcher = ["-DMPIEXEC_EXECUTABLE=" + str(args.launcher),
                "-DMPIEXEC_NUMPROC_FLAG=" + args.rank_flag,
                "-DMPIEXEC_PREFLAGS=" + ";".join(args.launcher_arg)]
    common = ["-DCMAKE_BUILD_TYPE=Release", "-DCMAKE_CXX_STANDARD=20", "-DCMAKE_CXX_COMPILER=" + str(args.cxx),
              "-DKokkos_DIR=" + str(args.kokkos_dir), "-DMARS_STRUCTURED_TEST_RANKS=2;4",
              "-DMARS_STRUCTURED_FORCE_MPI=ON"] + launcher
    if is_gpu(args.backend):
        common += ["-DCMAKE_{}_COMPILER={}".format(gpu, args.gpu_compiler),
                   "-DCMAKE_{}_STANDARD=20".format(gpu), "-DCMAKE_{}_STANDARD_REQUIRED=ON".format(gpu),
                   "-DCMAKE_{}_ARCHITECTURES={}".format(gpu, args.arch)]
    if args.backend == "cuda":
        common.append("-DCMAKE_CUDA_HOST_COMPILER=" + str(args.cxx))
    configure = [args.cmake, "-S", str(args.source), "-B", str(build), "-G", "Ninja"] + common + [
        "-DCMAKE_INSTALL_PREFIX=" + str(install), "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
        "-DMARS_ENABLE_KOKKOS=ON", "-DMARS_ENABLE_MPI=ON", "-DMARS_ENABLE_DISTRIBUTED_BACKEND=ON",
        "-DMARS_ENABLE_CUDA=" + ("ON" if args.backend == "cuda" else "OFF"),
        "-DMARS_ENABLE_HIP=" + ("ON" if args.backend == "hip" else "OFF"),
        "-DMARS_ENABLE_CUDAUVM=OFF",
        "-DMARS_ENABLE_RDMA=" + ("ON" if args.transport == "direct" else "OFF"),
        "-DMARS_ENABLE_UNSTRUCTURED=OFF", "-DMARS_ENABLE_SERIAL_BACKEND=OFF", "-DMARS_ENABLE_AMR_BACKEND=OFF",
        "-DMARS_ENABLE_EXAMPLES=OFF", "-DMARS_ENABLE_FEM_EXAMPLES=OFF", "-DMARS_ENABLE_KOKKOS_KERNELS=OFF",
        "-DMARS_ENABLE_ADIOS2=OFF", "-DMARS_ENABLE_VTK=OFF", "-DMARS_ENABLE_TESTS=ON",
        "-DCMAKE_DISABLE_FIND_PACKAGE_GTest=TRUE", "-DFETCHCONTENT_FULLY_DISCONNECTED=ON",
        "-DFETCHCONTENT_SOURCE_DIR_GOOGLETEST=" + str(args.gtest_source)]
    return [
        ("configure", configure),
        ("build", [args.cmake, "--build", str(build), "--parallel", str(args.jobs), "--target", "mars"] + list(SUITES)[:5]),
        ("test-list", [args.ctest, "--test-dir", str(build), "-L", "structured", "--show-only=json-v1"]),
        ("tests", [args.ctest, "--test-dir", str(build), "-L", "structured", "--parallel", "1", "--no-tests=error", "--output-on-failure",
                   "--test-output-size-passed", "1048576", "--test-output-size-failed", "1048576",
                   "--verbose", "--output-junit", str(args.work / "tests.xml")]),
        ("install", [args.cmake, "--install", str(build)]),
        ("consumer-configure", [args.cmake, "-S", str(args.source / "tests/structured/install-consumer"),
                                "-B", str(consumer), "-G", "Ninja", "-DCMAKE_PREFIX_PATH=" + str(install)] + common),
        ("consumer-build", [args.cmake, "--build", str(consumer), "--parallel", str(args.jobs)]),
        ("consumer-list", [args.ctest, "--test-dir", str(consumer), "-L", "structured-install", "--show-only=json-v1"]),
        ("consumer-tests", [args.ctest, "--test-dir", str(consumer), "-L", "structured-install", "--parallel", "1", "--no-tests=error",
                            "--test-output-size-passed", "1048576", "--test-output-size-failed", "1048576",
                            "--output-on-failure", "--verbose", "--output-junit", str(args.work / "consumer.xml")]),
    ]


def run(args):
    if args.work.exists():
        raise ValueError("Use a new work directory for each backend and transport")
    args.work.mkdir(parents=True)
    status = {"status": "FAILED", "backend": args.backend, "transport": args.transport,
              "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
              "host": platform.platform(), "child_signal_policy": child_signal_policy()}
    env = execution_environment({})
    env["MARS_STRUCTURED_EXPECT_BACKEND"] = {"serial": "Serial", "openmp": "OpenMP"}.get(args.backend, args.backend.upper())
    env["MARS_GATE_LAYOUT"] = "single"
    env["OMP_NUM_THREADS"] = str(getattr(args, "threads", 1))
    env["OMP_PROC_BIND"] = "false"
    env["KOKKOS_NUM_THREADS"] = str(getattr(args, "threads", 1))
    status["environment"] = environment_receipt(env)
    try:
        status["source_manifest_sha256"] = verify_snapshot(args.source)
        status["dependencies"] = verify_dependencies(args)
        commands = plan(args)
        (args.work / "commands.json").write_text(json.dumps(commands, indent=2) + "\n")
        if args.plan_only:
            status["status"] = "PLANNED_NOT_VALIDATED"
            return status
        status["tools"] = {name: tool_version(executable, args.work / (name + "-version.log"))
                           for name, executable in (("cxx", args.cxx), ("gpu-compiler", args.gpu_compiler),
                                                    ("cmake", args.cmake), ("ctest", args.ctest), ("launcher", args.launcher)) if executable}
        lists = {}
        for name, command in commands:
            print("Running " + name, flush=True)
            with (args.work / (name + ".log")).open("w") as log:
                result = run_child(command, stdout=log, stderr=subprocess.STDOUT, env=env)
            if result.returncode:
                raise ValueError("{} failed (exit {}); inspect {}.log".format(name, result.returncode, name))
            if name == "configure":
                cache = (args.work / "build/CMakeCache.txt").read_text()
                match = re.search(r"^MPI_CXX_COMPILER:[^=]+=(.+)$", cache, re.M)
                if not match:
                    raise ValueError("Configured build did not record MPI_CXX_COMPILER")
                status["mpi_compiler"] = tool_version(match.group(1), args.work / "mpi-cxx-version.log", required=False)
            if name in ("test-list", "consumer-list"):
                listing = json.loads((args.work / (name + ".log")).read_text())
                expected = SUITES if name == "test-list" else CONSUMER
                if {test["name"] for test in listing["tests"]} != set(expected) or len(listing["tests"]) != len(expected):
                    raise ValueError("CTest registration is missing or duplicates required tests")
                lists[name] = listing
            if name == "tests":
                status["tests"] = audit_junit(args.work / "tests.xml", SUITES, args.backend, args.transport, True, threads=getattr(args, "threads", 1))
            if name == "consumer-tests":
                status["consumer"] = audit_junit(args.work / "consumer.xml", CONSUMER, args.backend, args.transport, False, threads=getattr(args, "threads", 1))
        if getattr(args, "multi_node", False):
            env["MARS_GATE_LAYOUT"] = "multi"
            for phase, expected, gtest in (("tests", SUITES, True), ("consumer-tests", CONSUMER, False)):
                command = list(dict(commands)[phase])
                xml = args.work / (phase + "-multinode.xml")
                command[command.index("--output-junit") + 1] = str(xml)
                with (args.work / (phase + "-multinode.log")).open("w") as log:
                    process = run_child(command, stdout=log, stderr=subprocess.STDOUT, env=env)
                if process.returncode:
                    raise ValueError("Multi-node " + phase + " failed")
                status[phase + "-multinode"] = audit_junit(xml, expected, args.backend, args.transport,
                                                           gtest, 2, getattr(args, "threads", 1))
        status["artifacts"] = {}
        for listing in lists.values():
            for test in listing["tests"]:
                if test["name"] in ("morton_test", "hilbert_test", "mesh_manager_test", "dof_handler_test", "mesh_ghost_layer_test", "structured_install"):
                    path = test_binary(test)
                    status["artifacts"][str(path)] = digest(path)
        for directory in (args.work / "build", args.work / "consumer"):
            for name in ("CMakeCache.txt", "compile_commands.json"):
                path = directory / name
                if path.is_file():
                    status["artifacts"][str(path)] = digest(path)
        for path in sorted((args.work / "install").rglob("*")):
            if path.is_file():
                status["artifacts"][str(path)] = digest(path)
        if verify_snapshot(args.source) != status["source_manifest_sha256"]:
            raise ValueError("Source snapshot changed during validation")
        if verify_dependencies(args) != status["dependencies"]:
            raise ValueError("Dependency files changed during validation")
        status["status"] = "PASS"
        return status
    except Exception as error:
        status["error"] = str(error)
        raise
    finally:
        artifacts = status.setdefault("artifacts", {})
        for pattern in ("*.log", "*.xml", "commands.json", "build/CMakeCache.txt", "consumer/CMakeCache.txt"):
            for path in args.work.glob(pattern):
                if path.is_file():
                    artifacts[str(path)] = digest(path)
        status["finished_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        (args.work / "result.json").write_text(json.dumps(status, indent=2, sort_keys=True) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("source", "work", "cxx", "kokkos-dir", "gtest-source", "launcher"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--backend", choices=("cuda", "hip", "serial", "openmp"), required=True)
    parser.add_argument("--transport", choices=("direct", "staged"), default="direct")
    parser.add_argument("--gpu-compiler", type=Path)
    parser.add_argument("--arch", default="native")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--multi-node", action="store_true")
    parser.add_argument("--rank-flag", default="-n")
    parser.add_argument("--launcher-arg", action="append", default=[])
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--cmake", default="cmake")
    parser.add_argument("--ctest", default="ctest")
    parser.add_argument("--plan-only", action="store_true", help="Check inputs and write commands without building or running")
    args = parser.parse_args()
    for name in ("source", "work", "cxx", "gpu_compiler", "kokkos_dir", "gtest_source", "launcher"):
        if getattr(args, name) is None:
            continue
        path = getattr(args, name).expanduser()
        setattr(args, name, path.absolute() if name in ("cxx", "gpu_compiler", "launcher") else path.resolve())
    if args.threads <= 0 or (args.backend == "openmp" and args.threads < 2):
        parser.error("OpenMP needs at least two threads; all thread counts must be positive")
    if is_gpu(args.backend) and (args.gpu_compiler is None or args.arch == "native"):
        parser.error("GPU builds need --gpu-compiler and --arch")
    if args.jobs <= 0 or any(character in args.arch for character in ";\n\r"):
        parser.error("Use a positive job count and one explicit GPU architecture")
    for name in ("cxx", "gpu_compiler", "launcher"):
        if getattr(args, name) is not None and not getattr(args, name).is_file():
            parser.error("--{} must name an existing executable".format(name.replace("_", "-")))
    if args.launcher.name == "srun" and "--export=ALL" not in args.launcher_arg:
        args.launcher_arg.insert(0, "--export=ALL")
    try:
        result = run(args)
    except Exception as error:
        print("FAILED: " + str(error), file=sys.stderr)
        return 1
    print(result["status"] + ": " + str(args.work / "result.json"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
