#!/usr/bin/env python3
"""Build the Kokkos-OFF GPU configuration and audit the complete existing release suite."""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
import run

CSTONE = "2a95e63b0ca5454712ab83f27505055c078887dc"
TARGETS = ("mars_cvfem_graph", "mars_cvfem_graph_tet", "mars_cvfem_ho_matfree_test", "mars_cvfem_poisson",
           "mars_ex1_poisson", "mars_ho_dist_apply_test", "mars_lid_driven_cavity", "mars_poiseuille_flow",
           "mars_step01_mesh", "mars_step02_volume", "mars_step03_poisson", "mars_step04_cvfem",
           "mars_step05_galerkin", "mars_step06_heat", "mars_tgv")
NAMES = {"marsReleaseMeshHex", "marsReleaseMeshTet", "marsReleaseMeshHexScaled", "marsReleaseMeshHexCycled",
         "marsReleaseMeshHex32", "marsReleaseMeshTet16", "marsReleaseHexAssembly", "marsReleaseHexAssemblyScaled",
         "marsReleaseHexNumbering", "marsReleaseTetAssembly", "marsReleasePoisson", "marsReleasePoissonCycled",
         "marsReleasePoisson_np4", "marsReleaseCavity", "marsReleaseChannel", "marsReleaseTgv", "marsReleaseHoMatfree",
         "marsReleaseHoDistApply_np4"}
for name in ("Ex1Poisson", "TgvDecay", "Poiseuille", "Step01", "Step02Volume", "Step02Moment", "Step03Poisson8",
             "Step03Poisson16", "Step04Cvfem16", "Step04Cvfem32", "Step05Hex16", "Step05Hex32", "Step05Tet8",
             "Step06Heat10", "Step06Heat20"):
    NAMES.update("marsRelease" + name + "_np" + str(rank) for rank in (1, 4))


def audit_junit(path):
    cases = list(ET.parse(path).getroot().iter("testcase"))
    if len(cases) != 48 or {case.get("name") for case in cases} != NAMES:
        raise ValueError("Release suite must contain exactly the existing 48 tests, including Hypre and netCDF cases")
    for case in cases:
        if (case.get("status", "run") not in ("run", "passed") or
                any(case.find(tag) is not None for tag in ("failure", "error", "skipped"))):
            raise ValueError("Required release test failed/skipped: " + case.get("name", "unknown"))
    return sorted(NAMES)


def execute(command, log, env):
    print("Running " + log.stem, flush=True)
    with log.open("w") as stream:
        result = run.run_child(command, stdout=stream, stderr=subprocess.STDOUT, env=env)
    if result.returncode:
        raise ValueError("Failed " + str(log))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("source", "work", "cxx", "gpu-compiler", "launcher", "hypre", "netcdf"):
        p.add_argument("--" + name, type=Path, required=True)
    p.add_argument("--backend", choices=("cuda", "hip"), required=True)
    p.add_argument("--arch", required=True)
    p.add_argument("--launcher-arg", action="append", default=[])
    p.add_argument("--jobs", type=int, default=8)
    p.add_argument("--multi-node", action="store_true")
    args = p.parse_args()
    args.source, args.work = args.source.resolve(), args.work.resolve()
    if args.work.exists():
        p.error("Use a fresh regression directory")
    args.work.mkdir(parents=True)
    result = {"status": "FAILED", "backend": args.backend, "commands": []}
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", OMP_NUM_THREADS="1", MPICH_GPU_SUPPORT_ENABLED="1", MPICH_GPU_IPC_ENABLED="0")
    env["MARS_GATE_LAYOUT"] = "single"
    try:
        result["source_manifest_sha256"] = run.verify_snapshot(args.source)
        deps = args.source.parent / "dependencies"
        result["cornerstone_manifest_sha256"] = run.verify_dependency(deps / "cornerstone", CSTONE)
        result["gtest_manifest_sha256"] = run.verify_dependency(deps / "googletest", run.GTEST_SHA)
        libraries = {}
        for name, prefix, pattern in (("hypre", args.hypre, "libHYPRE.so"), ("netcdf", args.netcdf, "libnetcdf.so")):
            found = [prefix / folder / pattern for folder in ("lib", "lib64") if (prefix / folder / pattern).is_file()]
            if not found:
                raise ValueError("Missing " + pattern + " under " + str(prefix))
            libraries[name] = found[0].resolve()
        config = (args.hypre / "include/HYPRE_config.h").read_text()
        if not re.search(r'^#define HYPRE_USING_' + args.backend.upper() + r'\b', config, re.M):
            raise ValueError("Hypre does not enable the requested GPU backend")
        result["external_dependencies"] = {str(path): run.digest(path) for path in list(libraries.values()) +
            [args.hypre / "include/HYPRE_config.h", args.netcdf / "include/netcdf.h"]}
        execute([sys.executable, "-c", "import numpy; print(numpy.__version__)"], args.work / "numpy.log", env)
        result["tools"] = {name: run.tool_version(exe, args.work / (name + "-version.log"))
                           for name, exe in (("cxx", args.cxx), ("gpu", args.gpu_compiler), ("cmake", "cmake"))}
        build = args.work / "build"
        gpu = args.backend.upper()
        flags = ["--launcher=" + str(args.launcher)] + ["--launcher-arg=" + flag for flag in args.launcher_arg]
        command = ["cmake", "-S", str(args.source), "-B", str(build), "-G", "Ninja", "-DCMAKE_BUILD_TYPE=Release",
            "-DCMAKE_CXX_COMPILER=" + str(args.cxx), "-DCMAKE_{}_COMPILER={}".format(gpu, args.gpu_compiler),
            "-DCMAKE_{}_ARCHITECTURES={}".format(gpu, args.arch), "-DHIP_GPU_ARCHITECTURES=" + args.arch,
            "-DMARS_ENABLE_CUDA=" + ("ON" if gpu == "CUDA" else "OFF"), "-DMARS_ENABLE_HIP=" + ("ON" if gpu == "HIP" else "OFF"),
            "-DMARS_ENABLE_KOKKOS=OFF", "-DMARS_ENABLE_MPI=ON", "-DMARS_ENABLE_UNSTRUCTURED=ON", "-DMARS_ENABLE_HYPRE=ON",
            "-DMARS_ENABLE_FEM_EXAMPLES=ON", "-DMARS_ENABLE_TESTS=ON", "-DMARS_ENABLE_EXAMPLES=OFF",
            "-DMARS_ENABLE_MARSIR=OFF", "-DMARS_ENABLE_DEV_EXAMPLES=OFF", "-DMARS_ENABLE_SEGREGATED=OFF",
            "-DMARS_ENABLE_SERIAL_BACKEND=OFF", "-DMARS_ENABLE_AMR_BACKEND=OFF", "-DMARS_ENABLE_ADIOS2=OFF",
            "-DFETCHCONTENT_FULLY_DISCONNECTED=ON", "-DFETCHCONTENT_SOURCE_DIR_CORNERSTONE_FETCH=" + str(deps / "cornerstone"),
            "-DCMAKE_DISABLE_FIND_PACKAGE_GTest=TRUE", "-DFETCHCONTENT_SOURCE_DIR_GOOGLETEST=" + str(deps / "googletest"),
            "-DCMAKE_DISABLE_FIND_PACKAGE_PkgConfig=ON", "-DHYPRE_INCLUDE_DIR=" + str(args.hypre / "include"),
            "-DHYPRE_LIBRARY=" + str(libraries["hypre"]), "-DNETCDF_INCLUDE_DIRS=" + str(args.netcdf / "include"),
            "-DNETCDF_LIBRARY=" + str(libraries["netcdf"]), "-DPython3_EXECUTABLE=" + sys.executable,
            "-DMPIEXEC_EXECUTABLE=" + str(args.source / "tests/structured/cluster/mpi_launch.py"), "-DMPIEXEC_NUMPROC_FLAG=-n",
            "-DMPIEXEC_PREFLAGS=" + ";".join(flags), "-DMARS_RELEASE_TEST_RANKS=4"]
        if gpu == "CUDA":
            command.append("-DCMAKE_CUDA_HOST_COMPILER=" + str(args.cxx))
        result["commands"].append(command)
        execute(command, args.work / "configure.log", env)
        command = ["cmake", "--build", str(build), "--parallel", str(args.jobs), "--target"] + list(TARGETS)
        result["commands"].append(command)
        execute(command, args.work / "build.log", env)
        execute(["ctest", "--test-dir", str(build), "-L", "release", "--show-only=json-v1"], args.work / "inventory.json", env)
        inventory = json.loads((args.work / "inventory.json").read_text())
        if len(inventory["tests"]) != 48 or {t["name"] for t in inventory["tests"]} != NAMES:
            raise ValueError("Required release tests not registered")
        result["tests"] = {}
        for layout in (("single", "multi") if args.multi_node else ("single",)):
            env["MARS_GATE_LAYOUT"] = layout
            xml = args.work / ("release-" + layout + ".xml")
            command = ["ctest", "--test-dir", str(build), "-L", "release", "--parallel", "1", "--no-tests=error",
                       "--output-on-failure", "--verbose", "--output-junit", str(xml)]
            result["commands"].append(command)
            execute(command, args.work / ("release-" + layout + ".log"), env)
            result["tests"][layout] = audit_junit(xml)
        binaries = [path for path in build.rglob("*") if path.is_file() and path.name in TARGETS]
        if {path.name for path in binaries} != set(TARGETS):
            raise ValueError("Missing release driver binaries")
        result["binaries"] = {str(path): run.digest(path) for path in binaries}
        if run.verify_snapshot(args.source) != result["source_manifest_sha256"]:
            raise ValueError("Release sources changed")
        if run.verify_dependency(deps / "cornerstone", CSTONE) != result["cornerstone_manifest_sha256"]:
            raise ValueError("Cornerstone sources changed")
        for path, checksum in result["external_dependencies"].items():
            if run.digest(path) != checksum:
                raise ValueError("External dependency changed")
        result["status"] = "PASS"
    except Exception as error:
        result["error"] = str(error)
        print("FAILED: " + str(error), file=sys.stderr)
    result["artifacts"] = {str(path): run.digest(path) for path in args.work.rglob("*") if path.is_file() and
                           (path.suffix in (".xml", ".log") or path.name in ("CMakeCache.txt", "inventory.json"))}
    (args.work / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(result["status"] + ": " + str(args.work / "result.json"))
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
