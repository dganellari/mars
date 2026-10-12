"""Early recovery must preserve the target and count all correction work."""
import contextlib
import io
import json
from pathlib import Path
import unittest
from unittest.mock import patch

import simple_pressure_replay as replay
import simple_pressure_profile as profile
import test_simple_pressure_profile as profiles
import test_simple_pressure_replay as fixtures


class BudgetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        fixtures.ReplayTests.setUpClass()
        cls.checker = fixtures.ReplayTests.checker

    @classmethod
    def tearDownClass(cls):
        fixtures.ReplayTests.tearDownClass()

    def setUp(self):
        self.fixture = fixtures.ReplayTests()
        self.fixture.setUp(); self.addCleanup(self.fixture.doCleanups)
        self.root = self.fixture.root
        self.capture = self.fixture.captured_record()
        self.exe = self.root/'SECRET-exe'; self.exe.write_text('executable')
        config = dict(family='hypre', rtol=1e-10, atol=0., max_iterations=200,
                      options=dict(solver=dict(type='gmres'), precond=dict(type='boomeramg')))
        self.values, _ = profile.resolve(config, (0., 1e-10), profiles.measured_defaults())
        self.source = self.root/'pressure.profile'; self.source.write_text(profile.text_profile(self.values))
        capture = replay.startup.read_json(self.capture/'capture.json')
        capture['libraries'] = {}; (self.capture/'capture.json').write_text(json.dumps(capture))
        self.make_arm('baseline', 0, 100, 10, 10.)
        self.make_arm('early', 1, int(self.values['kdim']), 20, 4.)

    def make_arm(self, name, cycles, initial, correction, seconds):
        path = self.root/name; path.mkdir()
        solution = [1.+i/8. for i in range(12)]
        fixtures.write_solution(path/'result', solution, 1)
        fixtures.write_solution(path/'result/initial', [0.]*12, 1)
        fixtures.write_expansion(path/'result/expansion', solution, [0.]*12, 1)
        config = dict(self.values, recovery_rounds=4, retain_update_expansion=1)
        if cycles: config['initial_restart_cycles'] = cycles
        config_path = path/'gpu-reference.settings'; config_path.write_text(profile.text_profile(config))
        report = dict(self.values, effective_levels=3, effective_relax_3=18,
            result_iterations=initial, result_solve_error=256, result_global_error=256,
            result_fatal_error=0, result_converged=0, result_reported=1e-4,
            recovery_requested_rounds=4, recovery_rounds=1, recovery_iterations=correction,
            recovery_stop=0, recovery_controls_restored=1, recovery_audit_version=1, recovery_audit_steps=1,
            expansion_requested=1, expansion_version=1, expansion_step=1,
            initial_restart_cycles=cycles, initial_maxiter=self.values['kdim'] if cycles else 200,
            initial_controls_restored=1, solve_timing_version=1, solve_seconds=seconds)
        report['effective_relax_1'] = report.pop('relax_down')
        report['effective_relax_2'] = report.pop('relax_up')
        (path/'result/rank-000000.report').write_text(profile.text_profile(report))
        (path/'result/rank-000000.recovery').write_text(
            'initial 1 1.000001 limit 1e-10\ncorrection 1 iterations {} return 0 global 0 converged 1 reported 1e-12\n'
            'candidate 1e-11 2e-11\n'.format(correction))
        (path/'result/rank-000000.correction-audit').write_text(
            'mars-pressure-correction-audit-v1\nstep 1 rhs 1 1.000001 correction 1e-11 2e-11 ideal 1e-11 2e-11 '
            'rounding 0 1e-14 candidate 1e-11 2e-11 target 1e-10 1.000001e-10\n')
        record = dict(schema=replay.SCHEMA, backend='mars', profile='gpu-reference', exit_code=0,
            gpu_profile_source=str(self.source), recovery_rounds=4, retain_update_expansion=True,
            expansion_version=1, recovery_audit_version=1, environment={},
            command=['launcher', str(self.exe), str(self.capture/'system'), str(config_path), str(path/'result')],
            capture_sha256=replay.digest(self.capture/'capture.json'),
            inputs=replay.hashes([self.exe, self.source, config_path, self.capture/'capture.json']))
        if cycles: record['initial_restart_cycles'] = cycles
        (path/'launch-start.json').write_text(json.dumps(record))
        record['files'] = replay.hashes([p for p in (path/'result').rglob('*') if p.is_file()] + [path/'launch-start.json'])
        (path/'replay.json').write_text(json.dumps(record))

    def compare(self, status=0, checker=None):
        output = self.root/'public.json'
        if output.exists(): output.unlink()
        args = ['compare-budget', '--capture-run', str(self.capture), '--baseline-run', str(self.root/'baseline'),
                '--early-run', str(self.root/'early'), '--checker', str(checker or self.checker), '--output', str(output)]
        with patch.object(replay, 'loaded_libraries_match'), contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(replay.main(args), status)
        public = replay.startup.read_json(output)
        self.assertNotIn('SECRET', json.dumps(public)); self.assertNotIn(str(self.root), json.dumps(public))
        return public

    def refresh_outputs(self, name='early'):
        path = self.root/name; record = replay.startup.read_json(path/'replay.json')
        record['files'] = replay.hashes([Path(p) for p in record['files']])
        (path/'replay.json').write_text(json.dumps(record))

    def test_total_work_and_measured_time_with_independent_target(self):
        public = self.compare()
        self.assertEqual(public['outcome'], 'less_work_and_shorter_measured_solve')
        self.assertTrue(public['both_candidates_accepted'])
        self.assertFalse(public['nonlinear_convergence_verified'])
        self.assertTrue(public['only_initial_budget_changed'])

    def test_unavailable_checker_is_not_a_capture_failure(self):
        nonexecutable = self.root/'SECRET-nonexecutable'; nonexecutable.write_text('not executable')
        for checker in (self.root/'SECRET-missing-checker', self.root, nonexecutable):
            with self.subTest(checker=checker), patch.object(replay, 'check_candidate') as check:
                public = self.compare(1, checker)
                self.assertEqual(public['failed_check'], 'checker_executable')
                check.assert_not_called()

    def test_invalid_capture_remains_a_capture_failure(self):
        record = self.capture/'capture.json'; record.write_text('{}')
        self.assertEqual(self.compare(1)['failed_check'], 'capture_identity')

    def test_failed_pair_cannot_claim_benefit(self):
        path = self.root/'early/result/expansion/rank-000000.high'
        path.write_bytes(bytes(path.stat().st_size)); self.refresh_outputs()
        public = self.compare()
        self.assertEqual(public['outcome'], 'candidate_rejected')
        self.assertFalse(public['both_candidates_accepted'])

    def test_correction_work_is_not_hidden(self):
        path = self.root/'early/result/rank-000000.report'
        report = replay.numeric_file(path)
        report.update(recovery_iterations=200)
        path.write_text(profile.text_profile(report))
        trace = self.root/'early/result/rank-000000.recovery'
        trace.write_text(trace.read_text().replace('iterations 20', 'iterations 200'))
        self.refresh_outputs()
        public = self.compare()
        self.assertEqual(public['total_krylov_work_relation'], 'higher')
        self.assertEqual(public['outcome'], 'no_clear_work_and_time_benefit')

    def test_stale_output_fails_before_residual_check(self):
        (self.root/'early/result/rank-000000.report').write_text('SECRET corrupted')
        self.assertEqual(self.compare(1)['failed_check'], 'early_evidence')

    def test_invalid_budget_and_timing_are_rejected(self):
        path = self.root/'early/result/rank-000000.report'; original = path.read_text()
        for key, value in (('initial_controls_restored', 0), ('initial_maxiter', 200),
                           ('result_iterations', 201), ('solve_seconds', 0), ('solve_seconds', float('nan'))):
            with self.subTest(key=key, value=value):
                path.write_text(original); report = replay.numeric_file(path); report[key] = value
                path.write_text(profile.text_profile(report)); self.refresh_outputs()
                self.assertEqual(self.compare(1)['failed_check'], 'early_budget_metadata')

    def test_changed_environment_is_not_controlled_pair(self):
        path = self.root/'early'
        for name in ('launch-start.json', 'replay.json'):
            record = replay.startup.read_json(path/name); record['environment'] = {'MARS_HYPRE_FLEXGMRES':'1'}
            (path/name).write_text(json.dumps(record))
        self.refresh_outputs()
        self.assertEqual(self.compare(1)['failed_check'], 'controlled_pair')

    def test_early_request_without_retention_never_launches(self):
        args = ['replay', '--capture-run', str(self.capture), '--backend', 'mars', '--profile', 'gpu-reference',
                '--gpu-profile-pair', str(self.root), '--executable', str(self.exe),
                '--initial-restart-cycles', '1', '--output-dir', str(self.root/'invalid-request')]
        with patch.object(replay.subprocess, 'run') as launch, contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(replay.main(args), 1)
            launch.assert_not_called()
        self.assertEqual(replay.startup.read_json(self.root/'invalid-request/public.json')['failed_check'], 'profile_selection')


if __name__ == '__main__':
    unittest.main()
