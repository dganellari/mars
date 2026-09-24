#!/usr/bin/env python3
"""Per-commit regression tracker for MARS example runs.

Modelled on OpenAccel's tools/python/regression_tests/quick_test.py: run the cases, store one
JSON record per commit, and report each quantity as a delta against the previous record. PASS/FAIL
is the exit code; the numbers are reported descriptively rather than thresholded, because a change
in mfRms or div*L/U is information, not automatically a failure.

Where ours differs, deliberately: OpenAccel tracks solver iterations only. We also capture the
physics -- mfRms, Q_out/Q_in, div*L/U, u_max, u_rms. Those are the numbers this project keeps
comparing across configurations, and on 2026-09-05 a baseline had to be re-run from scratch
because it existed only in a terminal that had scrolled away.

  ./scripts/mars_regression.py --list
  ./scripts/mars_regression.py --case pump_baseline
  ./scripts/mars_regression.py --compare            # last two records, no run
"""

import argparse
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

RESULTS = Path(__file__).resolve().parent / "regression_results"

# Each case is a name plus the argv tail. The launcher, mesh paths and rank count come from the
# environment so this file carries nothing machine-specific and nothing confidential.
#   MARS_REG_LAUNCH  e.g. "srun --account=csstaff --time=60 --nodes=1 --ntasks-per-node=1 $HOME/affinity/bind_numa.sh"
#   MARS_REG_BIN     e.g. "./examples/distributed/unstructured/mars_pump"
#   MESH_BIG / INLET_BIG / OUTLET_BIG  as already used interactively
CASES = {
    # The production reference. Everything else is judged against this.
    "pump_baseline": [
        "--solver=hypre", "--pspg", "--bj", "--rho=1000", "--nu=1e-4",
        "--inlet-velocity=0.5", "--opening-flux-source", "--outlet=do-nothing",
        "--source-ramp-steps=5000", "--dt=2e-6", "--num-steps=20000", "--steady-tol=1e-5",
    ],
    # Explicit Rhie-Chow on the flux, scalar tau. Kills the checkerboard, wrecks mass.
    "pump_vms_scalar": [
        "--solver=hypre", "--vms-stab", "--bj", "--rho=1000", "--nu=1e-4",
        "--inlet-velocity=0.5", "--opening-flux-source", "--outlet=do-nothing",
        "--source-ramp-steps=5000", "--dt=2e-6", "--num-steps=20000", "--steady-tol=1e-5",
    ],
    # Same, with OpenAccel's per-node D off the momentum diagonal.
    "pump_vms_diag": [
        "--solver=hypre", "--vms-stab", "--bj", "--rho=1000", "--nu=1e-4",
        "--inlet-velocity=0.5", "--opening-flux-source", "--outlet=do-nothing",
        "--source-ramp-steps=5000", "--dt=2e-6", "--num-steps=20000", "--steady-tol=1e-5",
    ],
}
CASE_ENV = {"pump_vms_scalar": {"MARS_VMS_GLOBAL_TAU": "1"}}

# Last occurrence wins: these are printed every report interval, and the converged value is the
# one that matters.
PATTERNS = {
    "u_rms":    re.compile(r"u_rms=([0-9.eE+-]+)"),
    "u_max":    re.compile(r"u_max=([0-9.eE+-]+)"),
    "div_L_U":  re.compile(r"div\*L/U=([0-9.eE+-]+)"),
    "cg_p":     re.compile(r"cg_p=(-?\d+)"),
    "ratio":    re.compile(r"ratio=([0-9.eE+-]+)"),
    "mfRms":    re.compile(r"mean-free p RMS=([0-9.eE+-]+)"),
    "steps":    re.compile(r"Pump run complete:\s+(\d+) steps"),
    "ms_step":  re.compile(r"\(([0-9.]+) ms/step\)"),
}


def git_commit():
    return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], text=True).strip()


def parse(out):
    vals = {}
    for key, rx in PATTERNS.items():
        hits = rx.findall(out)
        if hits:
            try:
                vals[key] = float(hits[-1])
            except ValueError:
                pass
    # A run that aborted is not a result, whatever numbers it printed on the way down.
    vals["diverged"] = "[diverged]" in out
    vals["steady"] = "[steady]" in out
    return vals


def run_case(name, dry):
    launch = os.environ.get("MARS_REG_LAUNCH", "")
    binary = os.environ.get("MARS_REG_BIN", "./examples/distributed/unstructured/mars_pump")
    mesh = os.environ.get("MESH_BIG")
    if not mesh:
        sys.exit("MESH_BIG is not set (also needs INLET_BIG / OUTLET_BIG)")
    cmd = launch.split() + [binary] + CASES[name] + [
        f"--mesh={mesh}",
        f"--inlet-ss={os.environ['INLET_BIG']}",
        f"--outlet-ss={os.environ['OUTLET_BIG']}",
    ]
    env = dict(os.environ, MARS_SOLVE_TRACE="1", MARS_CHECKER="1", **CASE_ENV.get(name, {}))
    if dry:
        print(" ".join(cmd))
        return None
    t0 = time.time()
    proc = subprocess.run(cmd, env=env, capture_output=True, text=True)
    out = proc.stdout + proc.stderr
    log = RESULTS / f"{git_commit()}_{name}.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text(out)   # the log ALWAYS survives, which is the point
    rec = parse(out)
    rec.update(rc=proc.returncode, wall_s=round(time.time() - t0, 1), log=log.name)
    return rec


def load_records():
    RESULTS.mkdir(parents=True, exist_ok=True)
    recs = []
    for f in sorted(RESULTS.glob("*.json"), key=lambda p: p.stat().st_mtime):
        recs.append((f.stem, json.loads(f.read_text())))
    return recs


def store(commit, cases):
    RESULTS.mkdir(parents=True, exist_ok=True)
    f = RESULTS / f"{commit}.json"
    # Fold into the existing record: a single-case run must not drop the other cases from this
    # commit (the same trap OpenAccel's runner documents).
    merged = json.loads(f.read_text()) if f.exists() else {}
    merged.update(cases)
    f.write_text(json.dumps(merged, indent=2, sort_keys=True))
    return f


def label(cur, base, lower_is_better=True):
    if cur is None or base is None:
        return "new" if base is None else "n/a"
    if base == 0:
        return "same" if cur == 0 else f"{cur:+.3g}"
    rel = (cur - base) / abs(base)
    if abs(rel) < 1e-3:
        return "same"
    better = (rel < 0) == lower_is_better
    return f"{'better' if better else 'worse'} {rel*100:+.1f}%"


def report(cases, prev):
    cols = ["mfRms", "ratio", "div_L_U", "u_max", "steps", "ms_step"]
    # ratio's target is 1.0, not zero, so "lower is better" is meaningless for it.
    lower = {"mfRms": True, "div_L_U": True, "ms_step": True, "u_max": False, "steps": True}
    for name, rec in sorted(cases.items()):
        status = "FAIL" if rec.get("rc") else ("DIVERGED" if rec.get("diverged") else "pass")
        print(f"\n{name}  [{status}]  {rec.get('wall_s')}s  -> {rec.get('log')}")
        base = (prev or {}).get(name, {})
        for c in cols:
            cur, old = rec.get(c), base.get(c)
            tag = "" if c == "ratio" else f"   {label(cur, old, lower.get(c, True))}"
            shown = f"{cur:.4g}" if isinstance(cur, float) else str(cur)
            print(f"    {c:<9} {shown:>12}{tag}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", action="append", help="case name; repeatable. default: all")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--compare", action="store_true", help="report the last two records, run nothing")
    ap.add_argument("--dry-run", action="store_true", help="print the commands only")
    ap.add_argument("--no-store", action="store_true")
    a = ap.parse_args()

    if a.list:
        for k in CASES:
            print(k)
        return 0

    recs = load_records()
    if a.compare:
        if len(recs) < 2:
            print(f"need two records, have {len(recs)} in {RESULTS}")
            return 1
        (_, prev), (commit, cur) = recs[-2], recs[-1]
        print(f"comparing {commit} against the previous record")
        report(cur, prev)
        return 0

    names = a.case or list(CASES)
    prev = recs[-1][1] if recs else {}
    out = {}
    for n in names:
        if n not in CASES:
            sys.exit(f"unknown case {n}; --list to see them")
        r = run_case(n, a.dry_run)
        if r:
            out[n] = r
    if a.dry_run or not out:
        return 0
    report(out, prev)
    if not a.no_store:
        print(f"\nstored {store(git_commit(), out)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
