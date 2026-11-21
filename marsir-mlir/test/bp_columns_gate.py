#!/usr/bin/env python3
# BP3 in the column form (test/bp_columns.py, docs/design/marsir_bp3_columns_2026-11-21.md) in
# the CPU warp emulator on both targets, with the launch shapes of the MI300A kit: p = 3 with 32
# threads per element, p = 7 with 128. Also checks that the 3-D buffers left shared memory (the
# gfx942 build: shared bytes per block, no private-memory spill).
# Run from marsir-mlir/ after building build/tools/mir-opt:  python3 test/bp_columns_gate.py [--tsan]

import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
TMP = os.path.join(ROOT, "build")
tsan = ["--emu-tsan"] if "--tsan" in sys.argv else []
# Shared bytes per block of the gfx942 kernels (p = 3: two elements per 64-thread warp).
SHARED = {3: 1680, 7: 5328}
failed = []
for p, lanes in ((3, 32), (7, 128)):
    mir = os.path.join(TMP, f"bp3_columns_p{p}.mir")
    subprocess.run([sys.executable, os.path.join(HERE, "bp_columns.py"), str(p), mir], check=True, capture_output=True)
    opts = [f"--mir-file={mir}", "--columns", "--lane-owned", f"--lanes={lanes}", "--no-even-odd"]
    for target, cmd in (("nvidia", [os.path.join(HERE, "bp_gpu_pipeline.py")]),
                        ("amd", [os.path.join(HERE, "amd_target.py"), "bp_gpu_pipeline"])):
        r = subprocess.run([sys.executable, *cmd, "bp3", str(p), *opts, "--emulate", *tsan],
                           cwd=ROOT, capture_output=True, text=True)
        last = (r.stdout.strip().splitlines() or [""])[-1]
        ok = r.returncode == 0 and "PASS" in last and "DATA RACE" not in r.stdout + r.stderr
        err = re.search(r"max.err. = \S+", last)
        print(f"  {'ok  ' if ok else 'FAIL'} p={p} {target:6s} {err.group(0) if ok and err else (r.stderr or r.stdout)[-200:]}")
        failed += [] if ok else [f"p={p} {target}"]
    r = subprocess.run([sys.executable, os.path.join(HERE, "amd_target.py"), "bp_gpu_pipeline", "bp3", str(p), *opts],
                       cwd=ROOT, capture_output=True, text=True)
    s = os.path.join(ROOT, "generated", f"bp3_p{p}_gfx942.s")
    asm = open(s).read() if os.path.exists(s) else ""
    shared = re.search(r"\.group_segment_fixed_size: (\d+)", asm)
    private = re.search(r"\.private_segment_fixed_size: (\d+)", asm)
    ok = r.returncode == 0 and shared and int(shared.group(1)) == SHARED[p] and private and private.group(1) == "0"
    print(f"  {'ok  ' if ok else 'FAIL'} p={p} gfx942 shared {shared.group(1) if shared else '?'} B per block "
          f"(expected {SHARED[p]}), private {private.group(1) if private else '?'} B")
    failed += [] if ok else [f"p={p} resources"]
    if os.path.exists(s):
        os.remove(s)
print("BP3 COLUMNS GATE:", "PASS" if not failed else f"FAIL ({', '.join(failed)})")
sys.exit(1 if failed else 0)
