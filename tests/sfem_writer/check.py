"""Run one isolated writer case; a timeout or launcher failure never passes."""

import pathlib
import shutil
import subprocess
import sys

mode, directory, *command = sys.argv[1:]
path = pathlib.Path(directory)
if path.exists():
    if path.is_dir():
        shutil.rmtree(path)
    else:
        path.unlink()
result = subprocess.run(
    command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=120
)
print(result.stdout, end="")
if mode == "pass":
    valid = result.returncode == 0 and any(
        marker in result.stdout
        for marker in ("PASS: SFEM host/MPI check", "PASS: SFEM GPU/MPI check")
    )
else:
    valid = (
        result.returncode != 0
        and "FAULT_INJECTION_ARMED" in result.stdout
        and "SfemWriter: collective output failed" in result.stdout
        and "TEST_FAILURE:" not in result.stdout
    )
sys.exit(0 if valid else 1)
