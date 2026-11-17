#!/usr/bin/env python3
"""Record each CUDA compile and bound compiler failures without changing optimization."""

import hashlib
import json
import os
from pathlib import Path
import resource
import signal
import subprocess
import sys
import time


def main():
    log_dir = Path(sys.argv[1])
    command = sys.argv[2:]
    if not command:
        raise SystemExit("Missing compiler command")
    memory_gib = int(os.environ.get("MARS_AMR_COMPILE_MEMORY_GIB", "64"))
    timeout = float(os.environ.get("MARS_AMR_COMPILE_TIMEOUT", "600"))
    if memory_gib <= 0 or timeout <= 0:
        raise SystemExit("Compiler memory and time bounds must be positive")
    source = next((arg for arg in command if arg.endswith(".cu")), "cuda-command")
    key = hashlib.sha256(json.dumps(command).encode()).hexdigest()[:12]
    log_dir.mkdir(parents=True, exist_ok=True)
    receipt = log_dir / (Path(source).name + "-" + key + ".json")
    limited = sys.platform.startswith("linux")
    _, inherited_hard = resource.getrlimit(resource.RLIMIT_AS)
    requested = memory_gib * 1024**3
    bound = requested if inherited_hard == resource.RLIM_INFINITY else min(requested, inherited_hard)

    def limit_memory():
        resource.setrlimit(resource.RLIMIT_AS, (bound, bound))
        resource.setrlimit(resource.RLIMIT_CORE, (0, 0))

    result = {"source": source, "command": command, "status": "RUNNING",
              "memory_limit_gib": memory_gib, "address_limit_applied": limited,
              "core_dumps_disabled": limited,
              "effective_address_limit_bytes": bound if limited else None,
              "timeout_seconds": timeout}
    receipt.write_text(json.dumps(result, indent=2) + "\n")
    start = time.monotonic()
    try:
        process = subprocess.Popen(command, start_new_session=True,
                                   preexec_fn=limit_memory if limited else None)
        try:
            code = process.wait(timeout=timeout)
            result["status"] = "PASS" if code == 0 else "FAILED"
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            code = 124
            result["status"] = "TIMEOUT"
    except (OSError, subprocess.SubprocessError) as error:
        code = 127
        result["status"] = "LAUNCH_FAILED"
        result["error"] = str(error)
    result["exit_code"] = code
    result["elapsed_seconds"] = time.monotonic() - start
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    result["peak_child_rss_bytes"] = int(usage.ru_maxrss * (1024 if limited else 1))
    receipt.write_text(json.dumps(result, indent=2) + "\n")
    print(f"CUDA compile {Path(source).name}: {result['status']}, "
          f"{result['elapsed_seconds']:.1f}s, "
          f"peak child RSS {result['peak_child_rss_bytes'] / 1024**3:.3f} GiB", flush=True)
    return code if code >= 0 else 128 - code


if __name__ == "__main__":
    raise SystemExit(main())
