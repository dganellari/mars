"""Synthetic preparation tests; no solver or cluster commands are executed."""

import contextlib
import io
import json
import unittest

from netCDF4 import Dataset

import prepare_simple_snapshot_probe as probe
from simple_snapshot_compare import digest
import test_simple_snapshot_compare as snapshots


class ProbeTests(unittest.TestCase):
    def setUp(self):
        self.fixture = snapshots.SnapshotTests()
        self.fixture.setUp()
        f = self.fixture
        self.root = f.root
        (self.root / 'flow-metrics.csv').write_text(f.metrics.read_text())
        f.log.write_text(f.log.read_text() + 'linear_cache=1 halo_overlap=1 field_output=distributed profile=0\n')
        self.saved = json.loads(f.case.read_text())['arguments']
        (self.root / 'args.nul').write_bytes(b'\0'.join(value.encode() for value in self.saved) + b'\0')
        self.exe = self.root / 'synthetic-executable'
        self.exe.write_text('synthetic file; never executed\n')
        self.exe.chmod(0o700)
        (self.root / 'executable.sha256').write_text(digest(self.exe) + '  /original/private/path\n')
        self.output_dir = self.root / 'probe'
        self.public = self.root / 'preparation-public.json'

    def tearDown(self):
        self.fixture.tearDown()

    def run_probe(self, expected=0, **overrides):
        values = dict(baseline=self.root, case=self.fixture.case, reference_dir=self.fixture.reference,
                      executable=self.exe, output_dir=self.output_dir, output=self.public)
        values.update(overrides)
        argv = []
        for key, value in values.items():
            argv.append('--' + key.replace('_', '-'))
            if value is not True:
                argv.append(str(value))
        with contextlib.redirect_stdout(io.StringIO()) as stdout, contextlib.redirect_stderr(io.StringIO()) as stderr:
            code = probe.main(argv)
        self.assertEqual(code, expected)
        result = json.loads(self.public.read_text())
        self.assertTrue(all(type(value) in (str, bool) for value in result.values()))
        for text in (stdout.getvalue(), stderr.getvalue(), self.public.read_text()):
            self.assertNotIn(str(self.root), text)
            self.assertNotIn('/original/private/path', text)
            self.assertNotIn('987.654321', text)
        if expected:
            self.assertFalse((self.output_dir / 'args.nul').exists())
        return result

    def test_same_binary_arguments_controls_and_earliest_positive_state(self):
        result = self.run_probe()
        self.assertEqual(result['preparation_status'], 'ready')
        plan = json.loads((self.output_dir / 'probe.json').read_text())
        self.assertEqual(plan['iteration'], 50)
        self.assertEqual(plan['ranks'], 2)
        args = (self.output_dir / 'args.nul').read_bytes().decode().split('\0')[:-1]
        self.assertEqual(args[:len(self.saved)], self.saved)
        extra = dict(zip(args[len(self.saved)::2], args[len(self.saved)+1::2]))
        self.assertEqual(extra, {'--iterations': '50', '--report-every': '1',
                                '--residual-tol': '1e-06', '--mass-tol': '1e-06', '--change-tol': '1e-06',
                                '--linear-cache': '1', '--halo-overlap': '1', '--field-output': 'distributed', '--profile': '0'})

    def test_binary_change_rejected_before_launch_files_exist(self):
        self.exe.write_text('different file')
        self.assertEqual(self.run_probe(1)['failed_check'], 'binary_changed_or_unavailable')

    def test_missing_binary_rejected(self):
        self.exe.unlink()
        self.assertEqual(self.run_probe(1)['failed_check'], 'binary_changed_or_unavailable')

    def test_explicit_binary_change_records_both_hashes(self):
        old_hash = digest(self.exe)
        self.exe.write_text('rebuilt implementation')
        result = self.run_probe(allow_executable_change=True)
        self.assertFalse(result['binary_matches_baseline'])
        self.assertTrue(result['executable_change_allowed'])
        plan = json.loads((self.output_dir / 'probe.json').read_text())
        self.assertEqual(plan['baseline_executable_sha256'], old_hash)
        self.assertEqual(plan['executable_sha256'], digest(self.exe))
        self.assertTrue(plan['executable_change_allowed'])

    def test_binary_opt_in_still_checks_reference_deck(self):
        self.exe.write_text('rebuilt implementation')
        self.fixture.deck.write_text(self.fixture.deck.read_text() + '\n# altered\n')
        self.assertEqual(self.run_probe(1, allow_executable_change=True)['failed_check'], 'saved_arguments_or_controls')

    def test_binary_opt_in_still_checks_logged_controls(self):
        self.exe.write_text('rebuilt implementation')
        f = self.fixture
        f.log.write_text(f.log.read_text().replace('alpha_p=0.3', 'alpha_p=0.1'))
        self.assertEqual(self.run_probe(1, allow_executable_change=True)['failed_check'], 'saved_arguments_or_controls')

    def test_binary_opt_in_still_requires_executable(self):
        self.exe.unlink()
        self.assertEqual(self.run_probe(1, allow_executable_change=True)['failed_check'], 'binary_changed_or_unavailable')

    def test_failed_baseline_rejected(self):
        self.fixture.exit.write_text('143')
        self.assertEqual(self.run_probe(1)['failed_check'], 'baseline_completion')

    def test_nul_arguments_must_match_case(self):
        (self.root / 'args.nul').write_bytes(b'--relax-p\x000.1\0')
        self.assertEqual(self.run_probe(1)['failed_check'], 'saved_arguments_or_controls')

    def test_reference_deck_change_rejected(self):
        self.fixture.deck.write_text(self.fixture.deck.read_text() + '\n# altered\n')
        self.assertEqual(self.run_probe(1)['failed_check'], 'saved_arguments_or_controls')

    def test_mars_logged_control_change_rejected(self):
        f = self.fixture
        f.log.write_text(f.log.read_text().replace('alpha_p=0.3', 'alpha_p=0.1'))
        self.assertEqual(self.run_probe(1)['failed_check'], 'saved_arguments_or_controls')

    def test_runtime_options_preserved(self):
        f = self.fixture
        f.log.write_text(f.log.read_text().replace('linear_cache=1 halo_overlap=1', 'linear_cache=0 halo_overlap=0'))
        self.run_probe()
        args = json.loads((self.output_dir / 'probe.json').read_text())['arguments']
        values = dict(zip(args[::2], args[1::2]))
        self.assertEqual(values['--linear-cache'], '0')
        self.assertEqual(values['--halo-overlap'], '0')

    def test_missing_runtime_record_rejected(self):
        f = self.fixture
        f.log.write_text('\n'.join(line for line in f.log.read_text().splitlines() if not line.startswith('linear_cache=')))
        self.assertEqual(self.run_probe(1)['failed_check'], 'baseline_runtime_options')

    def test_reference_piece_times_must_match(self):
        with Dataset(str(self.fixture.ref_paths[1]), 'a') as ds:
            ds['time_whole'][1] = 40
        self.assertEqual(self.run_probe(1)['failed_check'], 'early_reference_state')

    def test_missing_reference_piece_rejected(self):
        self.fixture.ref_paths[1].unlink()
        self.assertEqual(self.run_probe(1)['failed_check'], 'early_reference_state')

    def test_no_run_beyond_the_short_budget(self):
        self.assertEqual(self.run_probe(1, max_iteration=10)['failed_check'], 'early_reference_state')

    def test_no_fractional_steady_iteration(self):
        for path in self.fixture.ref_paths:
            with Dataset(str(path), 'a') as ds:
                ds['time_whole'][1] = 1.5
        self.assertEqual(self.run_probe(1)['failed_check'], 'early_reference_state')

    def test_no_only_initial_state(self):
        paths = self.fixture.ref_paths
        self.fixture.write_exodus(paths[0], [4, 0, 2, 1], times=[0])
        self.fixture.write_exodus(paths[1], [3, 1, 5, 2], times=[0])
        self.assertEqual(self.run_probe(1)['failed_check'], 'early_reference_state')

    def test_earliest_positive_state_must_precede_baseline_final(self):
        for path in self.fixture.ref_paths:
            with Dataset(str(path), 'a') as ds:
                ds['time_whole'][:] = [0, 100, 200]
        self.assertEqual(self.run_probe(1)['failed_check'], 'early_reference_state')

    def test_no_launch_from_a_profiled_baseline(self):
        f = self.fixture
        f.log.write_text(f.log.read_text().replace('profile=0', 'profile=1'))
        self.assertEqual(self.run_probe(1)['failed_check'], 'baseline_runtime_options')


if __name__ == '__main__':
    unittest.main()
