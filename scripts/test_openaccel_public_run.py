"""Runner failure-path tests using an explicitly fake executable, not CFD validation."""
import contextlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from run_openaccel_public import run, sha256


class PublicRunTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.bundle = self.root / "bundle"
        self.bundle.mkdir()
        (self.bundle / "input.i").write_text("# fake runner fixture\n")
        (self.bundle / "channel.exo").write_bytes(b"fake mesh; never passed to OpenAccel")
        shutil.copyfile(Path(__file__).with_name("run_openaccel_public.py"),
                        self.bundle / "run_openaccel_public.py")
        (self.bundle / "manifest.json").write_text(json.dumps({
            "fixture": "public_channel", "iterations": 2,
            "sha256": {name: sha256(self.bundle / name)
                       for name in ("input.i", "channel.exo", "run_openaccel_public.py")}}))
        self.executable = self.root / "fake-solver"
        self.output = self.root / "result"

    def execute(self, body, capture_interior=False):
        self.executable.write_text("#!/bin/sh\n" + body)
        self.executable.chmod(0o700)
        with patch.dict(os.environ, {}, clear=True), contextlib.redirect_stdout(io.StringIO()):
            return run(self.bundle, self.executable, self.output, capture_interior)

    def enable_capture_bundle(self):
        manifest_path = self.bundle / "manifest.json"
        manifest = json.loads(manifest_path.read_text())
        for name in ("openaccel_reference_check.py", "contract_v1.json", "provenance.json"):
            (self.bundle / name).write_text("{}")
            manifest["sha256"][name] = sha256(self.bundle / name)
        manifest_path.write_text(json.dumps(manifest))

    def test_uninstrumented_binary_cannot_pass_capture(self):
        self.enable_capture_bundle()
        with contextlib.redirect_stderr(io.StringIO()):
            code = self.execute('test -d "$MARS_OPENACCEL_EXPORT_DIR" || exit 9\n'
                                'test "$MARS_OPENACCEL_PUBLIC_FIXTURE" = public_channel || exit 9\n'
                                "printf 'Iter = 1\\nIter = 2\\n'\n", True)
        self.assertEqual(code, 1)
        result = json.loads((self.output / "run.json").read_text())
        self.assertEqual(result["returncode"], 0)
        self.assertEqual(result["status"], "interior_capture_failed")
        self.assertIn("no exports", result["capture_error"])
        self.assertFalse(result["full_contract_passed"])

    def test_capture_requires_hashed_instrumentation_bundle(self):
        with self.assertRaises((KeyError, OSError)):
            self.execute("exit 0\n", True)
        self.assertFalse(self.output.exists())

    def boundary_bundle(self):
        self.enable_capture_bundle()
        path = self.bundle / 'manifest.json'
        manifest = json.loads(path.read_text())
        manifest['require_boundary_capture'] = True
        name = 'openaccel_boundary_check.py'
        (self.bundle / name).write_text('# fake hashed checker fixture\n')
        manifest['sha256'][name] = sha256(self.bundle / name)
        path.write_text(json.dumps(manifest))
        return dict(blocks={(s,c,n): {} for s in ('pressure.interior','momentum.interior')
                            for c in (1,2) for n in range(1536)},
                    hashes={str(i):'hash' for i in range(4)}, inputs={},
                    producer='openaccel',signature=['public_channel'])

    def test_required_boundary_cannot_silently_pass_without_exports(self):
        capture = self.boundary_bundle()
        with patch('openaccel_reference_check.load_dump',return_value=capture), contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(self.execute("printf 'Iter = 1\\nIter = 2\\n'\n",True),1)
        result = json.loads((self.output/'run.json').read_text())
        self.assertIn('missing boundary stages',result['capture_error'])

    def test_boundary_completion_still_does_not_claim_parity(self):
        capture = self.boundary_bundle()
        boundary = ({(s,c,21+s%3):{} for s in range(6) for c in (1,2)}, {'boundary':'hash'})
        with patch('openaccel_reference_check.load_dump',return_value=capture), \
                patch('openaccel_boundary_check.load_boundary',return_value=boundary):
            self.assertEqual(self.execute("printf 'Iter = 1\\nIter = 2\\n'\n",True),0)
        result = json.loads((self.output/'run.json').read_text())
        self.assertEqual(result['status'],'boundary_capture_completed')
        self.assertEqual(result['coverage'],'frozen-interior-boundary-blocks')
        self.assertFalse(result['numerical_parity_verified'])
        self.assertFalse(result['full_contract_passed'])

    def test_success_does_not_claim_convergence(self):
        self.assertEqual(self.execute("printf 'Iter = 1\nIter = 2\n'\n"), 0)
        result = json.loads((self.output / "run.json").read_text())
        self.assertEqual(result["status"], "runtime_smoke_completed")
        self.assertFalse(result["convergence_verified"])
        self.assertFalse(result["numerical_parity_verified"])

    def update_bundle(self):
        capture = self.boundary_bundle()
        path = self.bundle / 'manifest.json'
        manifest = json.loads(path.read_text())
        manifest.pop('require_boundary_capture')
        manifest['require_update_capture'] = True
        name = 'openaccel_update_check.py'
        (self.bundle / name).write_text('# fake hashed checker fixture\n')
        manifest['sha256'][name] = sha256(self.bundle / name)
        path.write_text(json.dumps(manifest))
        return capture

    def test_update_capture_requires_exports(self):
        capture = self.update_bundle()
        with patch('openaccel_reference_check.load_dump',return_value=capture), contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(self.execute("printf 'Iter = 1\\nIter = 2\\n'\n",True),1)
        result = json.loads((self.output/'run.json').read_text())
        self.assertEqual(result['status'],'interior_capture_failed')

    def test_update_capture_does_not_claim_numerical_parity(self):
        capture = self.update_bundle()
        updates = {(s,i,n,0): {} for s,count in enumerate((425,425,425,9216,96,96,32,96,96,1))
                   for i in (1,2) for n in range(count)}
        with patch('openaccel_reference_check.load_dump',return_value=capture), \
                patch('openaccel_update_check.load_updates',return_value=(updates,{'updates':'hash'})):
            self.assertEqual(self.execute("printf 'Iter = 1\\nIter = 2\\n'\n",True),0)
        result = json.loads((self.output/'run.json').read_text())
        self.assertEqual(result['status'],'update_capture_completed')
        self.assertFalse(result['numerical_parity_verified'])
        self.assertFalse(result['full_contract_passed'])

    def test_solver_failure_preserves_exit_and_log(self):
        self.assertEqual(self.execute("echo 'missing library' >&2\nexit 7\n"), 7)
        self.assertIn("missing library", (self.output / "run.log").read_text())

    def test_launch_uses_python36_compatible_arguments(self):
        real_popen = subprocess.Popen

        def python36_popen(args, *, cwd, env, stdout, stderr,
                           universal_newlines=False, errors=None):
            self.assertIsInstance(cwd, str)
            self.assertTrue(universal_newlines)
            return real_popen(args, cwd=cwd, env=env, stdout=stdout, stderr=stderr,
                              universal_newlines=universal_newlines, errors=errors)

        with patch("run_openaccel_public.subprocess.Popen", python36_popen):
            self.assertEqual(self.execute("printf 'Iter = 1\\nIter = 2\\n'\n"), 0)

    def test_zero_exit_without_both_iterations_fails(self):
        self.assertEqual(self.execute("echo 'Iter = 1'\n"), 1)

    def test_modified_mesh_refused_before_execution(self):
        (self.bundle / "channel.exo").write_bytes(b"changed")
        with self.assertRaisesRegex(ValueError, "checksum mismatch"):
            self.execute("exit 0\n")
        self.assertFalse(self.output.exists())

    def test_existing_output_preserved(self):
        self.output.mkdir()
        marker = self.output / "keep"
        marker.write_text("previous result")
        with self.assertRaises(FileExistsError):
            self.execute("exit 0\n")
        self.assertEqual(marker.read_text(), "previous result")

    def test_multiple_ranks_refused(self):
        with patch.dict(os.environ, {"SLURM_NTASKS": "4"}), self.assertRaisesRegex(ValueError, "one MPI rank"):
            run(self.bundle, self.executable, self.output)
        self.assertFalse(self.output.exists())


if __name__ == "__main__":
    unittest.main()
