#!/usr/bin/env python3
"""Regenerate the committed kernel files and check they are byte-identical to git HEAD.

The pipelines write into generated/, so git is the reference: a change to a pass or a
pipeline that alters any committed PTX shows up as a diff. Run from marsir-mlir/ after
building build/tools/mir-opt:

  python3 test/check_generated.py                    all files in MANIFEST, against git HEAD
  python3 test/check_generated.py bp3 dss            only entries whose output starts with these
  python3 test/check_generated.py --snapshot=DIR     also copy every regenerated file to DIR
  python3 test/check_generated.py --against=DIR      compare with DIR instead of git HEAD

--snapshot/--against guard a refactor: snapshot with the unchanged compiler, then check each
step against it. Some committed files are older than the compiler (not regenerated after a
schedule change), so HEAD alone cannot tell a refactor's diff from that drift.

Files under generated/ not in MANIFEST are research artifacts of the older make_ptx_*.sh
scripts; they are listed as not covered, never silently skipped. On failure the regenerated
files are left in place for inspection; `git checkout generated/` restores them."""
import filecmp
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# output file -> pipeline arguments (python3 test/<script> <args>)
MANIFEST = {}
for spec in ("bp1", "bp3", "bp1_tri", "bp3_tri"):
    for p in (1, 3, 7):
        MANIFEST[f"{spec}_p{p}_sm90.ptx"] = ["bp_gpu_pipeline.py", spec, str(p)]
MANIFEST["hl_full_p7_sm90.ptx"] = ["hl_gpu_pipeline.py", "7"]
MANIFEST["hl_full_p7_nochain_sm90.ptx"] = ["hl_gpu_pipeline.py", "7", "--no-chain"]
for p in (1, 3):
    MANIFEST[f"hl_cellwise_p{p}_sm90.ptx"] = ["hl_gpu_pipeline.py", str(p), "--cellwise"]
for p in (1, 3, 7):
    MANIFEST[f"hl_cellwise_p{p}_f32_sm90.ptx"] = ["hl_gpu_pipeline.py", str(p), "--cellwise", "--f32"]
for p in (1, 3, 7):
    MANIFEST[f"dss_p{p}_sm90.ptx"] = ["dss_gpu_pipeline.py", str(p)]
MANIFEST["precond_p7_sm90.ptx"] = ["dss_gpu_pipeline.py", "7", "--precond"]
MANIFEST["bicgstab_p7_sm90.ptx"] = ["dss_gpu_pipeline.py", "7", "--bicgstab"]
MANIFEST["pmg_sm90.ptx"] = ["dss_gpu_pipeline.py", "--pmg"]
MANIFEST["pmg_f32_sm90.ptx"] = ["dss_gpu_pipeline.py", "--pmg", "--f32"]
for p in (3, 7):
    MANIFEST[f"lines_p{p}_sm90.ptx"] = ["dss_gpu_pipeline.py", str(p), "--lines"]
for p in (1, 3, 7):
    for b in (8, 32):
        MANIFEST[f"block_dss_p{p}_b{b}_sm90.ptx"] = ["block_dss_gpu_pipeline.py", str(p), str(b)]
MANIFEST["tet_laplacian_sm90.ptx"] = ["tet_gpu_pipeline.py"]


def main():
    opt = {a.split("=", 1)[0]: a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--")}
    only = [a for a in sys.argv[1:] if not a.startswith("--")]
    snapshot, against = opt.get("--snapshot"), opt.get("--against")
    if snapshot:
        os.makedirs(snapshot, exist_ok=True)
    gen = os.path.join(ROOT, "generated")
    uncovered = sorted(f for f in os.listdir(gen) if f not in MANIFEST)
    failed = []
    for out, args in MANIFEST.items():
        if only and not any(out.startswith(o) for o in only):
            continue
        run = subprocess.run([sys.executable, os.path.join(HERE, args[0]), *args[1:]], cwd=ROOT,
                             capture_output=True, text=True)
        missing = not os.path.exists(os.path.join(gen, out))
        if against:
            same = not missing and filecmp.cmp(os.path.join(gen, out), os.path.join(against, out), shallow=False)
        else:
            same = subprocess.run(["git", "diff", "--quiet", "--", os.path.join("generated", out)], cwd=ROOT).returncode == 0
        if snapshot and not missing:
            shutil.copy2(os.path.join(gen, out), os.path.join(snapshot, out))
        if run.returncode != 0 or not same or missing:
            why = (f"pipeline exit {run.returncode}" if run.returncode else "missing" if missing
                   else f"differs from {against or 'HEAD'}")
            failed.append(out)
            print(f"FAIL {out}: {why}  (test/{' '.join(args)})")
            if run.returncode:
                print("     " + "\n     ".join((run.stderr or run.stdout).strip().splitlines()[-5:]))
        else:
            print(f"ok   {out}")
    print(f"not covered (no pipeline command): {', '.join(uncovered)}")
    print("GENERATED CHECK:", "PASS" if not failed else f"FAIL ({len(failed)})")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
