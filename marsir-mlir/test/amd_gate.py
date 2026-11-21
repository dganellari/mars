#!/usr/bin/env python3
# The AMD (gfx942) kernels in the CPU warp emulator: every family test/amd_target.py builds, run
# with the AMD target (64-lane wavefronts, MFMA, the AMD launch shapes) against the same references
# as the NVIDIA gates. amd_target.py also rejects any alloca outside a kernel's entry block (a
# dynamic stack allocation per element on AMDGPU, which the emulator cannot see), so a kernel that
# builds here has none. Not covered: execution on the GPU (Beverin kits, marsir-mlir/AMD.md).
# Run from marsir-mlir/ after building build/tools/mir-opt:  python3 test/amd_gate.py [--tsan]

import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

CASES = [(f"{spec} p={p}", ["bp_gpu_pipeline", spec, str(p)])
         for spec in ("bp1", "bp3", "bp1_tri", "bp3_tri") for p in (1, 3, 7)]
CASES += [("hl full p=7 (MFMA tensor cores)", ["7"])]

tsan = ["--emu-tsan"] if "--tsan" in sys.argv else []
failed = []
for name, args in CASES:
    r = subprocess.run([sys.executable, os.path.join(HERE, "amd_target.py"), *args, "--emulate", *tsan],
                       cwd=ROOT, capture_output=True, text=True)
    last = (r.stdout.strip().splitlines() or [""])[-1]
    ok = r.returncode == 0 and "PASS" in last and "DATA RACE" not in r.stdout
    print(f"  {'ok  ' if ok else 'FAIL'} {name:34s} {last[last.find('max|err|'):][:40] if ok else (r.stderr or r.stdout).strip()[-200:]}")
    if not ok:
        failed.append(name)
print("AMD EMULATOR GATE:", "PASS" if not failed else f"FAIL ({len(failed)})")
sys.exit(1 if failed else 0)
