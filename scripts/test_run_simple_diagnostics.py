import json
from pathlib import Path
import signal
import stat
import subprocess
import sys
import tempfile
import time
import unittest

import run_simple_diagnostics as capture


class CaptureTests(unittest.TestCase):
    def paths(self, root):
        return [root / name for name in ('SECRET.log', 'SECRET.exit', 'public.json')]

    def command(self, paths, source):
        return [sys.executable, str(Path(capture.__file__).resolve()),
                '--log', str(paths[0]), '--exit-file', str(paths[1]), '--output', str(paths[2]),
                '--', sys.executable, '-c', source]

    def test_status_and_private_bytes(self):
        for code, completion, expected in (
                (0, 'CONVERGED iterations=10 ranks=4 exchange_rounds=41', 'converged'),
                (2, 'NOT CONVERGED: iteration limit iterations=10 ranks=4 exchange_rounds=41', 'iteration_limit'),
                (255, 'ERROR: prepared Hypre: true residual evaluation failed', 'failed')):
            with self.subTest(code=code), tempfile.TemporaryDirectory() as work:
                paths = self.paths(Path(work))
                source = ('import sys; sys.stdout.buffer.write(b"SECRET \\xff\\n"); '
                          'sys.stdout.flush(); print({!r}, file=sys.stderr); sys.exit({})').format(completion, code)
                result = subprocess.run(self.command(paths, source), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                self.assertEqual(result.returncode, code, result.stderr)
                self.assertEqual(paths[1].read_text(), str(code) + '\n')
                self.assertIn(b'SECRET \xff\n', paths[0].read_bytes())
                summary = json.loads(paths[2].read_text())
                self.assertEqual(summary['capture_status'], 'finished')
                self.assertEqual(summary['run_status'], expected)
                self.assertNotIn('SECRET', paths[2].read_text())
                self.assertNotIn(b'SECRET', result.stdout + result.stderr)
                for path in paths:
                    self.assertEqual(stat.S_IMODE(path.stat().st_mode), 0o600)

    def test_live_snapshot_before_command_exits(self):
        with tempfile.TemporaryDirectory() as work:
            root = Path(work)
            paths = self.paths(root)
            release = root / 'release'
            source = ('import pathlib,time; print("ERROR: prepared Hypre: matrix refresh failed or changed the ParCSR object", flush=True); '
                      'p=pathlib.Path({!r})\nwhile not p.exists(): time.sleep(0.02)\n').format(str(release))
            process = subprocess.Popen(self.command(paths, source), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            try:
                deadline = time.monotonic() + 10
                while True:
                    if paths[2].exists():
                        text = paths[2].read_text()
                        if text:
                            summary = json.loads(text)
                            if summary.get('software_errors'):
                                break
                    self.assertIsNone(process.poll())
                    self.assertLess(time.monotonic(), deadline, 'live report not published')
                    time.sleep(0.02)
                self.assertEqual(summary['capture_status'], 'running')
                self.assertEqual(summary['run_status'], 'failed')
                self.assertEqual(paths[1].read_text(), '')
                release.touch()
                stdout, stderr = process.communicate(timeout=10)
                self.assertEqual(process.returncode, 1, stderr)
                self.assertEqual(paths[1].read_text(), '0\n')
                # A successful launcher cannot erase an earlier application error.
                self.assertEqual(json.loads(paths[2].read_text())['run_status'], 'failed')
            finally:
                release.touch()
                process.communicate(timeout=10)

    def test_launch_failure_is_redacted(self):
        with tempfile.TemporaryDirectory() as work:
            paths = self.paths(Path(work))
            code = capture.capture(['/SECRET/missing-command'], *paths)
            self.assertEqual(code, 127)
            summary = json.loads(paths[2].read_text())
            self.assertEqual(summary['capture_status'], 'launch_failed')
            self.assertEqual(summary['run_status'], 'failed')
            self.assertNotIn('SECRET', paths[2].read_text())

    def test_existing_outputs_prevent_launch(self):
        for occupied in range(3):
            with self.subTest(occupied=occupied), tempfile.TemporaryDirectory() as work:
                root = Path(work)
                paths = self.paths(root)
                paths[occupied].write_text('SECRET original')
                marker = root / 'launched'
                command = self.command(paths, 'import pathlib; pathlib.Path({!r}).touch()'.format(str(marker)))
                result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(marker.exists())
                self.assertEqual(paths[occupied].read_text(), 'SECRET original')
                self.assertNotIn(b'SECRET', result.stdout + result.stderr)
                self.assertNotIn(work.encode(), result.stdout + result.stderr)

    def test_signal_exit_is_not_success(self):
        with tempfile.TemporaryDirectory() as work:
            paths = self.paths(Path(work))
            source = 'import os,signal; os.kill(os.getpid(),signal.SIGTERM)'
            result = subprocess.run(self.command(paths, source), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.assertEqual(result.returncode, 128 + signal.SIGTERM)
            self.assertEqual(json.loads(paths[2].read_text())['run_status'], 'failed')


if __name__ == '__main__':
    unittest.main()
