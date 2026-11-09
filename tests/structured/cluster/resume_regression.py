#!/usr/bin/env python3
"""Run a missing two-node release suite using checksum-pinned existing executables."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

sys.dont_write_bytecode = True
import regression
import run


def verify_files(files):
    for name, checksum in files.items():
        if run.digest(Path(name)) != checksum:
            raise ValueError("Prepared file changed: " + name)


def quote(value):
    value = str(value)
    level = "="
    while "]" + level + "]" in value:
        level += "="
    return "[" + level + "[" + value + "]" + level + "]"


def mirror_inventory(inventory, old_build, destination, binaries):
    """Keep executables in place and put every test output in the new run directory."""
    old_build, destination = str(old_build), str(destination)

    def relocate(value):
        if isinstance(value, list):
            return [relocate(item) for item in value]
        if value == old_build:
            return destination
        if isinstance(value, str) and value not in binaries:
            return value.replace(old_build + "/", destination + "/")
        return value

    tests = []
    text = []
    for test in inventory["tests"]:
        updated = {"name": test["name"], "command": relocate(test["command"]),
                   "properties": [{"name": prop["name"], "value": relocate(prop["value"])}
                                  for prop in test["properties"]]}
        text.append("add_test(" + " ".join(quote(arg) for arg in [updated["name"]] + updated["command"]) + ")")
        for prop in updated["properties"]:
            value = prop["value"]
            if prop["name"] == "WORKING_DIRECTORY":
                Path(value).mkdir(parents=True, exist_ok=True)
            if isinstance(value, list):
                value = ";".join(str(item) for item in value)
            text.append("set_tests_properties(" + quote(updated["name"]) + " PROPERTIES " +
                        quote(prop["name"]) + " " + quote(value) + ")")
        tests.append(updated)
    # The original CMake configure created these output directories before registering tests.
    for rank in (1, 4):
        (Path(destination) / ("release_poiseuille_np" + str(rank))).mkdir(parents=True, exist_ok=True)
    (Path(destination) / "CTestTestfile.cmake").write_text("\n".join(text) + "\n")
    return tests


def canonical_tests(tests):
    return {test["name"]: {"command": test["command"],
                           "properties": {prop["name"]: prop["value"] for prop in test["properties"]}}
            for test in tests}


def resume(expected, work):
    work = Path(work).absolute()
    if work.exists():
        raise ValueError("Use a fresh continuation directory")
    work.mkdir(parents=True)
    result = {"status": "FAILED", "backend": expected["backend"], "prepared": expected["prepared"]}
    try:
        prepared = Path(expected["prepared"])
        source = Path(expected["source"])
        dependencies = source.parent / "dependencies"
        if run.verify_snapshot(source) != expected["source_manifest_sha256"]:
            raise ValueError("Prepared source manifest changed")
        for name, record in expected["dependencies"].items():
            if run.verify_dependency(dependencies / name, record["revision"]) != record["manifest_sha256"]:
                raise ValueError("Prepared dependency changed: " + name)
        verify_files(expected["files"])
        prior = json.loads((prepared / "result.json").read_text())
        if prior.get("backend") != expected["backend"] or prior.get("source_manifest_sha256") != expected["source_manifest_sha256"]:
            raise ValueError("Prepared backend or source does not match the continuation")
        if prior.get("error") != "Single-rank halo regression did not execute its test":
            raise ValueError("This continuation requires the recorded truncated-output audit failure")
        result["single_node"] = regression.audit_junit(prepared / "release-single.xml", prepared / "release-single.log")
        inventory = json.loads((prepared / "inventory.json").read_text())
        if {test["name"] for test in inventory["tests"]} != regression.NAMES or len(inventory["tests"]) != 50:
            raise ValueError("Prepared release inventory changed")
        original_build = prepared / "build"
        testdir = work / "ctest"
        testdir.mkdir()
        tests = mirror_inventory(inventory, original_build, testdir, expected["binaries"])
        env = dict(os.environ, MARS_GATE_LAYOUT="multi", OMP_NUM_THREADS="1",
                   MPICH_GPU_SUPPORT_ENABLED="1", MPICH_GPU_IPC_ENABLED="0", PYTHONDONTWRITEBYTECODE="1")
        listing = run.run_child(["ctest", "--test-dir", str(testdir), "-L", "release", "--show-only=json-v1"],
                                env=env, capture_output=True, text=True, check=True)
        (work / "inventory.json").write_text(listing.stdout)
        observed = json.loads(listing.stdout)["tests"]
        if canonical_tests(observed) != canonical_tests(tests):
            raise ValueError("Continuation changed test commands or properties")
        xml = work / "release-multi.xml"
        log = work / "release-multi.log"
        command = ["ctest", "--test-dir", str(testdir), "-L", "release", "--parallel", "1", "--no-tests=error",
                   "--test-output-size-passed", "1048576", "--test-output-size-failed", "1048576",
                   "--output-on-failure", "--verbose", "--output-junit", str(xml)]
        result["command"] = command
        regression.execute(command, log, env)
        result["multi_node"] = regression.audit_junit(xml, log)
        verify_files(expected["files"])
        if run.verify_snapshot(source) != expected["source_manifest_sha256"]:
            raise ValueError("Prepared sources changed during validation")
        for name, record in expected["dependencies"].items():
            if run.verify_dependency(dependencies / name, record["revision"]) != record["manifest_sha256"]:
                raise ValueError("Prepared dependency changed during validation: " + name)
        result.update(status="PASS", binaries=expected["binaries"], source_manifest_sha256=expected["source_manifest_sha256"])
    except Exception as error:
        result["error"] = str(error)
    result["artifacts"] = {str(path): run.digest(path) for path in work.glob("*") if path.is_file()}
    (work / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(result["status"] + ": " + str(work / "result.json"))
    if result.get("error"):
        print(result["error"], file=sys.stderr)
    return 0 if result["status"] == "PASS" else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expected", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    expected = json.loads(args.expected.read_text())
    return resume(expected, args.work)


if __name__ == "__main__":
    sys.exit(main())
