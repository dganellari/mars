#!/usr/bin/env python3
"""run_host_study.py: cached results are reused only when verified, and failures propagate.

A fake launcher and a fake host executable replay prepared synthetic outputs (see
test_duct_compare.synthetic), so the tests take seconds and exercise only the orchestration:
a missing executable, a launch failure and a nonzero exit fail; a verified cache is reused; a
cache from another executable, other arguments, edited outputs or without a manifest is an
error and is not rerun, except that --rerun-stale replaces results carrying a manifest (never
results without one).
"""
import os
import shutil
import stat
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_host_study as rhs  # noqa: E402
import test_duct_compare as tdc  # noqa: E402

LAUNCHER = """#!%s
import os, sys
os.execv(sys.argv[3], sys.argv[3:])   # drop "-n R"
"""
FAKE = """#!%s
# fake duct_host_run (%s): replay prepared outputs for --output-prefix
import os, shutil, sys
args = sys.argv[1:]
prefix = args[args.index("--output-prefix") + 1]
source = os.path.join(os.environ["FAKE_SOURCE"], os.path.basename(prefix))
for suffix in ("-fields.csv", "-metrics.csv"):
    shutil.copyfile(source + suffix, prefix + suffix)
with open(source + ".log") as f:
    sys.stdout.write(f.read())
sys.exit(int(os.environ.get("FAKE_EXIT", "0")))
"""


def script(path, text):
    with open(path, "w") as f:
        f.write(text)
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR)
    return path


class HostStudy(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix="duct-host-study-")
        self.source = os.path.join(self.dir, "source")
        os.makedirs(self.source)
        tdc.synthetic(self.source, levels=(4,), ranks=(1, 2))
        self.launcher = script(os.path.join(self.dir, "launch"), LAUNCHER % sys.executable)
        self.fake = script(os.path.join(self.dir, "fake_run"), FAKE % (sys.executable, "A"))
        self.out = os.path.join(self.dir, "out")
        os.environ["FAKE_SOURCE"] = self.source
        os.environ.pop("FAKE_EXIT", None)

    def tearDown(self):
        os.environ.pop("FAKE_EXIT", None)
        shutil.rmtree(self.dir)

    def run_study(self, executable=None, launcher=None, extra=()):
        return rhs.main([self.out, "--host-run", executable or self.fake, "--levels", "4", "--ranks", "1,2",
                         "--mpiexec", launcher or self.launcher, "--parity-only"] + list(extra))

    def test_clean_run_passes_and_reuses_verified_results(self):
        self.assertEqual(self.run_study(), 0)
        self.assertTrue(os.path.exists(os.path.join(self.out, "duct-4-2.run.json")))
        os.environ["FAKE_SOURCE"] = os.path.join(self.dir, "nowhere")   # a rerun would fail
        self.assertEqual(self.run_study(), 0)

    def test_missing_executable_fails_without_running(self):
        self.assertEqual(self.run_study(executable=os.path.join(self.dir, "no_such_binary")), 1)
        self.assertFalse(os.path.exists(os.path.join(self.out, "duct-4-1-fields.csv")))

    def test_launch_failure_fails(self):
        self.assertEqual(self.run_study(launcher=os.path.join(self.dir, "no_such_mpiexec")), 1)

    def test_nonzero_exit_fails(self):
        os.environ["FAKE_EXIT"] = "2"   # outputs are complete and valid; the exit status is not
        self.assertEqual(self.run_study(), 1)
        os.environ.pop("FAKE_EXIT")
        self.assertEqual(self.run_study(), 1)   # the recorded failure is not forgotten on reuse

    def test_cache_from_another_executable_is_rejected(self):
        self.assertEqual(self.run_study(), 0)
        other = script(os.path.join(self.dir, "fake_run_b"), FAKE % (sys.executable, "B"))
        self.assertEqual(self.run_study(executable=other), 1)

    def test_cache_from_other_arguments_is_rejected(self):
        self.assertEqual(self.run_study(), 0)
        self.assertEqual(self.run_study(extra=["--iterations", "10"]), 1)

    def test_edited_outputs_are_rejected(self):
        self.assertEqual(self.run_study(), 0)
        with open(os.path.join(self.out, "duct-4-2-fields.csv"), "a") as f:
            f.write("\n")
        self.assertEqual(self.run_study(), 1)

    def test_results_without_manifest_are_rejected(self):
        self.assertEqual(self.run_study(), 0)
        os.remove(os.path.join(self.out, "duct-4-1.run.json"))
        self.assertEqual(self.run_study(), 1)
        self.assertEqual(self.run_study(extra=["--rerun-stale"]), 1)   # never deleted: no manifest
        self.assertTrue(os.path.exists(os.path.join(self.out, "duct-4-1-fields.csv")))

    def test_rerun_stale_replaces_a_mismatched_result(self):
        self.assertEqual(self.run_study(), 0)
        other = script(os.path.join(self.dir, "fake_run_b"), FAKE % (sys.executable, "B"))
        self.assertEqual(self.run_study(executable=other, extra=["--rerun-stale"]), 0)
        self.assertEqual(self.run_study(executable=other), 0)   # now verified for B
        self.assertEqual(self.run_study(), 1)                  # and stale for A


if __name__ == "__main__":
    unittest.main()
