#!/usr/bin/env python3
"""Export only the source needed by the optional structured backend."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess


ROOT_FILES = {"CMakeLists.txt", "LICENSE", "LICENSE.md", "LICENSE.txt", "examples/CMakeLists.txt"}
ROOTS = ("cmake/", "base/", "core/", "io/", "kokkos/", "backend/distributed/", "tests/structured/")
SUFFIXES = {".cmake", ".in", ".hpp", ".h", ".hxx", ".cpp", ".cc", ".c", ".cu", ".cuh", ".py", ".sh"}
DEPENDENCIES = {"kokkos": "15dc143e5f39949eece972a798e175c4b463d4b8",
                "googletest": "15460959cbbfa20e66ef0b5ab497367e47fc0a04"}


def eligible(name):
    path = Path(name)
    if path.is_absolute() or ".." in path.parts:
        return False
    if any(part.startswith(".") for part in path.parts):
        return False
    if "pump" in name.lower() or "unstructured" in path.parts:
        return False
    if name in ROOT_FILES:
        return True
    if not name.startswith(ROOTS):
        return False
    return (path.name == "CMakeLists.txt" or path.suffix in SUFFIXES or
            (name.startswith("tests/structured/cluster/") and path.name == "README.md"))


def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def export(source, destination):
    source = source.resolve()
    destination = destination.resolve()
    if destination.exists():
        raise ValueError("Export destination must not exist")
    names = subprocess.check_output(
        ["git", "-C", str(source), "ls-files", "-z", "--cached", "--others", "--exclude-standard"],
    ).decode().split("\0")
    selected = sorted({name for name in names if name and eligible(name)})
    required = {"CMakeLists.txt", "backend/distributed/testsuite/StructuredTests.cmake",
                "tests/structured/install-consumer/main.cpp", "tests/structured/mars_structured_runtime.hpp"}
    if not required.issubset(selected):
        raise ValueError("Missing required structured sources: " + ", ".join(sorted(required - set(selected))))
    # Check names and links before opening any file. Build trees and notes are never copied.
    for name in selected:
        path = source / name
        if path.is_symlink() or not path.is_file():
            raise ValueError("Expected a regular source file: " + name)
        if any(parent.is_symlink() for parent in path.parents if parent != source and source in parent.parents):
            raise ValueError("Source parent is a symbolic link: " + name)
    destination.mkdir(parents=True)
    files = {}
    for name in selected:
        target = destination / "source" / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source / name, target)
        files[name] = digest(target)
    manifest = {
        "format": 1,
        "scope": "optional structured Kokkos only; no serial, AMR, unstructured, meshes or internal notes",
        "git_head": subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip(),
        "git_branch": subprocess.check_output(["git", "-C", str(source), "branch", "--show-current"], text=True).strip(),
        "files": files,
    }
    (destination / "source-manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def export_dependency(source, destination, revision):
    head = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    dirty = subprocess.check_output(["git", "-C", str(source), "status", "--porcelain"], text=True).strip()
    if head != revision or dirty:
        raise ValueError("Dependency must be a clean checkout at " + revision)
    names = subprocess.check_output(["git", "-C", str(source), "ls-files", "-z"]).decode().split("\0")
    files = {}
    destination.mkdir(parents=True)
    for name in sorted(name for name in names if name):
        path = source / name
        if path.is_symlink() or not path.is_file() or ".." in Path(name).parts:
            raise ValueError("Dependency source is not a regular file: " + name)
        target = destination / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        files[name] = digest(target)
    (destination / "dependency-manifest.json").write_text(
        json.dumps({"format": 1, "revision": revision, "files": files}, indent=2, sort_keys=True) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--destination", type=Path, required=True)
    parser.add_argument("--kokkos-source", type=Path, required=True)
    parser.add_argument("--gtest-source", type=Path, required=True)
    args = parser.parse_args()
    manifest = export(args.source, args.destination)
    export_dependency(args.kokkos_source, args.destination / "dependencies/kokkos", DEPENDENCIES["kokkos"])
    export_dependency(args.gtest_source, args.destination / "dependencies/googletest", DEPENDENCIES["googletest"])
    print("Exported {} source files to {}".format(len(manifest["files"]), args.destination))


if __name__ == "__main__":
    main()
