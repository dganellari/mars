#!/usr/bin/env python3
"""Build the pinned Kokkos snapshot with native CMake GPU language support."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.dont_write_bytecode = True
from run import digest, tool_version, verify_dependency, run_child, child_signal_policy

KOKKOS_SHA = "15dc143e5f39949eece972a798e175c4b463d4b8"
ARCHITECTURES = {("cuda", "80"): "AMPERE80", ("cuda", "86"): "AMPERE86",
                 ("cuda", "90"): "HOPPER90", ("hip", "gfx90a"): "AMD_GFX90A",
                 ("hip", "gfx942"): "AMD_GFX942"}


def commands(args):
    build, install = args.work / "build", args.work / "install"
    configure = [args.cmake, "-S", str(args.source), "-B", str(build), "-G", "Ninja",
                 "-DCMAKE_BUILD_TYPE=Release", "-DCMAKE_CXX_STANDARD=20",
                 "-DCMAKE_CXX_COMPILER=" + str(args.cxx), "-DCMAKE_INSTALL_PREFIX=" + str(install),
                 "-DCMAKE_INSTALL_LIBDIR=lib", "-DKokkos_ENABLE_SERIAL=ON",
                 "-DKokkos_ENABLE_OPENMP=" + ("ON" if args.backend == "openmp" else "OFF"),
                 "-DKokkos_ENABLE_CUDA=" + ("ON" if args.backend == "cuda" else "OFF"),
                 "-DKokkos_ENABLE_HIP=" + ("ON" if args.backend == "hip" else "OFF"),
                 "-DKokkos_ENABLE_TESTS=OFF", "-DKokkos_ENABLE_EXAMPLES=OFF", "-DKokkos_ENABLE_BENCHMARKS=OFF"]
    if args.backend in ("cuda", "hip"):
        if (args.backend, args.arch) not in ARCHITECTURES:
            raise ValueError("Unchecked GPU architecture: " + args.arch)
        gpu = args.backend.upper()
        configure += ["-DCMAKE_{}_COMPILER={}".format(gpu, args.gpu_compiler),
                      "-DCMAKE_{}_STANDARD=20".format(gpu), "-DCMAKE_{}_STANDARD_REQUIRED=ON".format(gpu),
                      "-DCMAKE_{}_ARCHITECTURES={}".format(gpu, args.arch),
                      "-DKokkos_ENABLE_COMPILE_AS_CMAKE_LANGUAGE=ON",
                      "-DKokkos_ARCH_{}=ON".format(ARCHITECTURES[(args.backend, args.arch)])]
    if args.backend == "cuda":
        configure.append("-DCMAKE_CUDA_HOST_COMPILER=" + str(args.cxx))
    steps = [("configure", configure),
             ("build", [args.cmake, "--build", str(build), "--parallel", str(args.jobs)]),
             ("install", [args.cmake, "--install", str(build)])]
    return steps


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("source", "work", "cxx"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--backend", choices=("cuda", "hip", "serial", "openmp"), required=True)
    parser.add_argument("--arch", default="native")
    parser.add_argument("--gpu-compiler", type=Path)
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--cmake", default="cmake")
    parser.add_argument("--plan-only", action="store_true")
    args = parser.parse_args()
    for name in ("source", "work", "cxx", "gpu_compiler"):
        if getattr(args, name) is None:
            continue
        path = getattr(args, name).expanduser()
        setattr(args, name, path.absolute() if name in ("cxx", "gpu_compiler") else path.resolve())
    if args.work.exists() or args.jobs <= 0:
        parser.error("Use a fresh work directory and a positive job count")
    args.work.mkdir(parents=True)
    result = {"status": "FAILED", "backend": args.backend, "arch": args.arch, "revision": KOKKOS_SHA, "child_signal_policy": child_signal_policy()}
    try:
        result["source_manifest_sha256"] = verify_dependency(args.source, KOKKOS_SHA)
        steps = commands(args)
        (args.work / "commands.json").write_text(json.dumps(steps, indent=2) + "\n")
        if args.plan_only:
            result["status"] = "PLANNED_NOT_VALIDATED"
        else:
            result["tools"] = {name: tool_version(executable, args.work / (name + "-version.log"))
                               for name, executable in (("cxx", args.cxx), ("gpu-compiler", args.gpu_compiler),
                                                        ("cmake", args.cmake)) if executable}
            for name, command in steps:
                with (args.work / (name + ".log")).open("w") as log:
                    process = run_child(command, stdout=log, stderr=subprocess.STDOUT)
                if process.returncode:
                    raise ValueError("Kokkos {} failed; inspect {}.log".format(name, name))
            if verify_dependency(args.source, KOKKOS_SHA) != result["source_manifest_sha256"]:
                raise ValueError("Kokkos source changed during build")
            result["cache_sha256"] = digest(args.work / "build/CMakeCache.txt")
            result["status"] = "PASS"
    except Exception as error:
        result["error"] = str(error)
    (args.work / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(result["status"] + ": " + str(args.work / "result.json"))
    return 1 if result["status"] == "FAILED" else 0


if __name__ == "__main__":
    sys.exit(main())
