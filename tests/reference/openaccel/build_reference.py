#!/usr/bin/env python3
"""build the pinned capture checkout with the working build's dependencies."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def git(source, *args):
    return subprocess.check_output(["git", "-C", str(source), *args],
                                   universal_newlines=True).strip()


def read_cache(path):
    values = {}
    for line in path.read_text().splitlines():
        if not line or line.startswith(("#", "//")) or "=" not in line or ":" not in line:
            continue
        key, value = line.split("=", 1)
        values[key.split(":", 1)[0]] = value
    return values


def configure_arguments(cache):
    required = ("CMAKE_C_COMPILER", "CMAKE_CXX_COMPILER", "MPI_C_COMPILER", "MPI_CXX_COMPILER",
                "Trilinos_DIR", "HDF5_DIR")
    for key in required:
        if not cache.get(key) or cache[key].endswith("-NOTFOUND"):
            raise ValueError("working build cache is missing " + key)
    yaml_keys = ("YAML-CPP_DIR",) if cache.get("YAML-CPP_DIR") and not cache["YAML-CPP_DIR"].endswith("-NOTFOUND") else (
        "YAML_DIR", "YAML_CPP_LIBRARIES", "YAML_CPP_INCLUDE_DIR")
    for key in yaml_keys:
        if not cache.get(key) or cache[key].endswith("-NOTFOUND"):
            raise ValueError("working build cache is missing " + key)
    keys = set(required) | set(yaml_keys) | {"CMAKE_Fortran_COMPILER", "CMAKE_BUILD_TYPE", "CMAKE_C_FLAGS", "CMAKE_CXX_FLAGS",
                           "CMAKE_C_FLAGS_RELEASE", "CMAKE_CXX_FLAGS_RELEASE",
                           "CMAKE_EXE_LINKER_FLAGS", "MPI_CXX_SKIP_MPICXX"}
    keys.update(key for key in cache if key.startswith("WITH_") and not key.endswith("-ADVANCED"))
    keys.update(key for key in ("YAML_DIR", "YAML_CPP_LIBRARIES", "YAML_CPP_INCLUDE_DIR")
                if cache.get(key) and not cache[key].endswith("-NOTFOUND"))
    arguments = ["-D" + key + "=" + cache[key] for key in sorted(keys) if key in cache]
    # This public capture deck uses Trilinos only. Keep optional solvers out of discovery.
    arguments += ["-DSPATIAL_DIM=3", "-DWITH_TRILINOS_SOLVER=ON",
                  "-DWITH_CANONICAL_SUFFIX=OFF",
                  "-DCMAKE_DISABLE_FIND_PACKAGE_PkgConfig=ON", "-DPKG_CONFIG_FOUND=FALSE",
                  "-DCMAKE_DISABLE_FIND_PACKAGE_HYPRE=ON", "-DHYPRE_FOUND=FALSE"]
    return arguments


def build(args):
    bundle = Path(__file__).resolve().parent
    source, baseline, destination = (p.resolve() for p in (args.source, args.baseline_build, args.destination))
    provenance = json.loads((bundle / "provenance.json").read_text())
    patch = bundle / "instrumentation.patch"
    if hashlib.sha256(patch.read_bytes()).hexdigest() != provenance["patch_sha256"]:
        raise ValueError("instrumentation patch checksum mismatch")
    if git(source, "rev-parse", "HEAD") != provenance["reference_revision"]:
        raise ValueError("reference source is not at the recorded pin")
    if git(source, "status", "--porcelain", "--untracked-files=no", "--ignore-submodules=untracked"):
        raise ValueError("tracked reference source is dirty; preserve it and use the clean pin")
    if git(source / "src/solver", "rev-parse", "HEAD") != provenance["solver_revision"]:
        raise ValueError("solver submodule is not at the recorded pin")
    if destination.exists():
        raise ValueError("destination exists; resume its build explicitly, or choose a fresh destination")
    cache = read_cache(baseline / "CMakeCache.txt")
    if Path(cache["CMAKE_HOME_DIRECTORY"]).resolve() != source:
        raise ValueError("baseline CMake cache belongs to a different source checkout")
    configure = configure_arguments(cache)
    for key in ("CMAKE_C_COMPILER", "CMAKE_CXX_COMPILER", "MPI_C_COMPILER", "MPI_CXX_COMPILER",
                "Trilinos_DIR", "HDF5_DIR"):
        if not Path(cache[key]).exists():
            raise ValueError("working dependency is unavailable; restore its uenv: " + cache[key])
    for key in ("CMAKE_Fortran_COMPILER", "YAML-CPP_DIR", "YAML_DIR", "YAML_CPP_LIBRARIES", "YAML_CPP_INCLUDE_DIR"):
        if cache.get(key) and not cache[key].endswith("-NOTFOUND") and not Path(cache[key]).exists():
            raise ValueError("working dependency is unavailable; restore its uenv: " + cache[key])
    subprocess.check_call(["git", "-C", str(source), "apply", "--check", str(patch)])
    subprocess.check_call(["git", "-C", str(source), "worktree", "add", "--detach",
                           str(destination), provenance["reference_revision"]])
    subprocess.check_call(["git", "-C", str(destination), "submodule", "update", "--init", "--recursive"])
    if git(destination / "src/solver", "rev-parse", "HEAD") != provenance["solver_revision"]:
        raise ValueError("new solver checkout differs from pin")
    subprocess.check_call(["git", "-C", str(destination), "apply", str(patch)])
    build_dir = destination / "build"
    command = ["cmake", "-S", str(destination), "-B", str(build_dir),
               "-G", "Unix Makefiles"] + configure
    (destination / "reference-build.json").write_text(json.dumps({
        "baseline_cache_sha256": hashlib.sha256(
            (baseline / "CMakeCache.txt").read_bytes()).hexdigest(),
        "provenance": provenance, "configure_command": command}, indent=2)+"\n")
    subprocess.check_call(command)
    subprocess.check_call(["cmake", "--build", str(build_dir), "--parallel", str(args.jobs)])
    print("Built instrumented reference: " + str(build_dir / "openaccel-3D.exe"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--baseline-build", type=Path, required=True)
    parser.add_argument("--destination", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("jobs must be positive")
    try:
        build(args)
    except (ValueError, KeyError, OSError, subprocess.CalledProcessError) as error:
        parser.exit(1, "ERROR: " + str(error) + "\n")
