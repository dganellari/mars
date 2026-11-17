import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class CompilerLauncherTest(unittest.TestCase):
    def run_compiler(self, program, timeout="10"):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as directory:
            root = Path(directory)
            compiler = root / "compiler.py"
            compiler.write_text(program)
            logs = root / "logs"
            command = [sys.executable, str(compiler), "source with spaces.cu"]
            result = subprocess.run(
                [sys.executable, str(Path(__file__).with_name("compile_cuda.py")), str(logs)] + command,
                env={**os.environ, "MARS_AMR_COMPILE_TIMEOUT": timeout,
                     "MARS_AMR_COMPILE_MEMORY_GIB": "2"}, capture_output=True, text=True,
            )
            receipts = list(logs.glob("*.json"))
            self.assertEqual(len(receipts), 1, result.stderr)
            receipt = json.loads(receipts[0].read_text())
            self.assertEqual(receipt["command"], command)
            self.assertEqual(receipt["source"], "source with spaces.cu")
            self.assertEqual(receipt["address_limit_applied"], sys.platform.startswith("linux"))
            return result, receipt

    def test_success_and_limit(self):
        result, receipt = self.run_compiler(
            "import resource,sys\n"
            "if sys.platform.startswith('linux'):\n"
            " assert resource.getrlimit(resource.RLIMIT_AS) == (2*1024**3, 2*1024**3)\n"
            " assert resource.getrlimit(resource.RLIMIT_CORE) == (0, 0)\n"
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(receipt["status"], "PASS")

    def test_compiler_failure(self):
        result, receipt = self.run_compiler("raise SystemExit(7)\n")
        self.assertEqual(result.returncode, 7)
        self.assertEqual(receipt["status"], "FAILED")

    def test_timeout(self):
        result, receipt = self.run_compiler("import time\ntime.sleep(60)\n", timeout="0.1")
        self.assertEqual(result.returncode, 124)
        self.assertEqual(receipt["status"], "TIMEOUT")

    def test_signal(self):
        result, receipt = self.run_compiler("import os,signal\nos.kill(os.getpid(), signal.SIGTERM)\n")
        self.assertEqual(result.returncode, 128 + 15)
        self.assertEqual(receipt["exit_code"], -15)

    def test_missing_compiler(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as directory:
            root = Path(directory)
            result = subprocess.run(
                [sys.executable, str(Path(__file__).with_name("compile_cuda.py")),
                 str(root / "logs"), str(root / "missing-compiler"), "test.cu"],
                capture_output=True, text=True,
            )
            self.assertEqual(result.returncode, 127)
            receipt = json.loads(next((root / "logs").glob("*.json")).read_text())
            self.assertEqual(receipt["status"], "LAUNCH_FAILED")
            self.assertIn("error", receipt)


if __name__ == "__main__":
    unittest.main()
