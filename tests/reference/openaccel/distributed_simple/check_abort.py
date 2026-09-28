"""An expected MPI abort must return promptly and report the injected fault."""
import subprocess
import sys

try:
    result = subprocess.run(sys.argv[1:], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            universal_newlines=True, timeout=20)
except subprocess.TimeoutExpired:
    sys.exit("FAIL: one-rank exchange error stranded a peer")
print(result.stdout)
expected = {"destroy-in-flight": "exchange destroyed before end",
            "double-begin": "another exchange is in flight",
            "reverse-in-flight": "another exchange is in flight",
            "metadata-in-flight": "another exchange is in flight",
            "end-without-begin": "end called without begin"}.get(
                sys.argv[-1], "invalid field pointer or buffer stride")
if result.returncode == 0 or expected not in result.stdout:
    sys.exit("FAIL: missing expected exchange abort")
print("PASS: one-rank error terminated the MPI job")
