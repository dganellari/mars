#!/usr/bin/env python3
"""Export the committed unstructured release sources with the checked structured overlay."""
import argparse
import json
from pathlib import Path
import subprocess

import run

ROOTS = ("backend/", "base/", "cmake/", "core/", "examples/", "io/", "scripts/", "tests/")
SUFFIXES = {".cmake", ".in", ".hpp", ".h", ".hxx", ".cpp", ".cc", ".c", ".cu", ".cuh", ".py", ".sh"}


def eligible(name):
    p = Path(name)
    if name in ("CMakeLists.txt", "LICENSE.md", "tests/data/poiseuille/poiseuille_hex_14k_elem.e"):
        return True
    if not name.startswith(ROOTS) or "data" in p.parts:
        return False
    return p.name == "CMakeLists.txt" or p.suffix in SUFFIXES


def export(repo, structured, destination):
    if destination.exists():
        raise ValueError("Use a new export destination")
    overlay = json.loads((structured.parent / "source-manifest.json").read_text())
    run.verify_snapshot(structured)
    revision = overlay["git_head"]
    entries = subprocess.check_output(["git", "-C", str(repo), "ls-tree", "-rz", revision]).split(b"\0")
    selected = []
    for entry in entries:
        if not entry:
            continue
        metadata, raw_name = entry.split(b"\t", 1)
        name = raw_name.decode()
        mode, kind, oid = metadata.decode().split()
        if eligible(name):
            if kind != "blob" or mode not in ("100644", "100755"):
                raise ValueError("Nonregular release source: " + name)
            selected.append((name, oid, mode))
    destination.mkdir(parents=True)
    files = {}
    for name, oid, mode in selected:
        path = destination / "source" / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(subprocess.check_output(["git", "-C", str(repo), "cat-file", "blob", oid]))
        path.chmod(0o755 if mode == "100755" else 0o644)
        files[name] = run.digest(path)
    for name, checksum in overlay["files"].items():
        source = structured / name
        target = destination / "source" / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(source.read_bytes())
        target.chmod(source.stat().st_mode & 0o777)
        files[name] = checksum
    manifest = {"format": 1, "git_head": revision, "files": files,
                "structured_manifest_sha256": run.digest(structured.parent / "source-manifest.json"),
                "scope": "Committed release code with structured overlay; excludes working AMR edits, private documents and geometry; only shipped Poiseuille mesh included"}
    (destination / "source-manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("repo", "structured", "destination"):
        p.add_argument("--" + name, type=Path, required=True)
    args = p.parse_args()
    manifest = export(args.repo, args.structured, args.destination)
    print("Exported {} files".format(len(manifest["files"])))


if __name__ == "__main__":
    main()
