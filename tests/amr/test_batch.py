import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


class BatchTest(unittest.TestCase):
    def run_batch(self, backend="cpu", cluster_status=0, audit_status=0, mutated=False):
        # Keep the fake spool separate from the bundle, as Slurm does.
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as temporary:
            root = Path(temporary)
            bundle = root / "bundle with spaces"
            source = bundle / "source/tests/amr"
            source.mkdir(parents=True)
            spool = root / "spool/job123"
            spool.mkdir(parents=True)
            script = spool / "slurm_script"
            text = Path(__file__).with_name("batch.sh").read_text()
            if mutated:
                text = text.replace("bundle=$2", 'bundle=$(cd "$(dirname "$0")" && pwd)')
            script.write_text(text)
            (bundle / "manifest.json").write_text('{"fixture": true}\n')
            trace = root / "trace.txt"
            (bundle / "verify.py").write_text(
                'import os\nfrom pathlib import Path\n'
                'with Path(os.environ["BATCH_TRACE"]).open("a") as f: f.write("verify\\n")\n'
            )
            (source / "cluster.sh").write_text(
                '#!/usr/bin/env bash\nset -eu\n'
                'printf "cluster %s\\n" "$1" >> "$BATCH_TRACE"\n'
                '[[ -d $3 ]]\n[[ $1 != cpu || -d $4 ]]\n'
                'mkdir -p "$2"\nexit "$BATCH_CLUSTER_STATUS"\n'
            )
            (bundle / "cornerstone").mkdir()
            (bundle / "kokkos").mkdir()
            (bundle / "audit.py").write_text(
                'import os\nfrom pathlib import Path\n'
                'with Path(os.environ["BATCH_TRACE"]).open("a") as f: f.write("audit\\n")\n'
                'raise SystemExit(int(os.environ["BATCH_AUDIT_STATUS"]))\n'
            )
            result = subprocess.run(
                ["bash", str(script), backend, str(bundle)], cwd=spool,
                env={**os.environ, "SLURM_JOB_ID": "123", "BATCH_TRACE": str(trace),
                     "BATCH_CLUSTER_STATUS": str(cluster_status), "BATCH_AUDIT_STATUS": str(audit_status)},
                capture_output=True, text=True,
            )
            if mutated:
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(trace.exists())
                return
            expected = cluster_status or audit_status
            self.assertEqual(result.returncode, expected, result.stderr)
            work = root / f"{bundle.name}-results-{backend}-123"
            self.assertEqual((work / "batch-exit.txt").read_text().strip(), str(expected))
            self.assertEqual(json.loads((work / "source-manifest.json").read_text()), {"fixture": True})
            events = ["verify", "cluster " + backend, "verify"]
            if not cluster_status:
                events.append("audit")
            self.assertEqual(trace.read_text().splitlines(), events)

    def test_spooled_cpu(self):
        self.run_batch()

    def test_spooled_cuda(self):
        self.run_batch("cuda")

    def test_spooled_hip(self):
        self.run_batch("hip")

    def test_build_failure(self):
        self.run_batch(cluster_status=7)

    def test_audit_failure(self):
        self.run_batch(audit_status=9)

    def test_script_directory_regression(self):
        self.run_batch(mutated=True)


if __name__ == "__main__":
    unittest.main()
